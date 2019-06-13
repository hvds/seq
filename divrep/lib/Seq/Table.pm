package Seq::Table;
use strict;
use warnings;

use parent 'DBIx::Class';
__PACKAGE__->load_components(qw{ BitField Core });

my %types = (
    id => {
        data_type => 'int',
        extra => { unsigned => 1 },
        is_auto_increment => 1,
    },
    uint => {
        data_type => 'int',
        extra => { unsigned => 1 },
    },
    bigint => {
        data_type => 'text',
        munge => [
            sub { "$_[0]" },
            sub { use Math::GMP; Math::GMP->new($_[0]) },
        ],
    },
    float => {
        data_type => 'float',
    },
    modlist => {
        data_type => 'text',
        munge => [
            sub { join ' ', ${ $_[0] } },
            sub { [ split /\s/, $_[0] ] },
        ],
    },
    flags => {
        data_type => 'int',
        extra => { unsigned => 1 },
        default_value => 0,
        bitfield => sub { [ split /\s+/, $_[0] ] },
        args => 'bitfield',
    },
);

sub define {
    my($class, $name, $table, $fields) = @_;
    my(%col, %inflate, @key, @order);

    for my $text (@$fields) {
        my($maybe, $key, $type, $args, $name) = $text =~ m{
            ^
            (maybe \s+)?
            (key \s+)?
            (\w+)
            (?: \s* \( ( .* ) \) )? \s+
            (\w+)
            \z
        }x or die "Can't parse field definition '$text'";
        my %def = %{ $types{$type} // die "Unknown type '$type' in '$text'" };
        $def{is_nullable} = 1 if defined $maybe;
        push @key, $name if defined $key;
        $inflate{$name} = delete $def{munge} if $def{munge};
        if (my $which = delete $def{args}) {
            defined $args or die "No args given for type '$type' in '$text'";
            $def{$which} = $def{$which}->($args);
        } else {
            !defined $args or die "Type '$type' takes no args in '$text'";
        }
        $col{$name} = \%def;
        push @order, $name;
    }

    $class->load_components(qw{ BitField Core });
    $class->source_name($name);
    $class->table($table);
    $class->add_columns(%col{@order});
    while (my($col, $munge) = each %inflate) {
        $class->inflate_column($col => {
            deflate => $munge->[0],
            inflate => $munge->[1],
        });
    }
    $class->set_primary_key(@key);
    $class->resultset_class('DBIx::Class::ResultSet::BitField');
    return;
}

sub Dump {
    my($self) = @_;
    use Data::Dumper;
    return Dumper($self->{_column_data});
}
