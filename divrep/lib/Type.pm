package Type;
use strict;
use warnings;

my %types = (
    t => 'Type::TauSeq',
    a => 'Type::AddSeq',
    s => 'Type::Semip',
);

sub new {
    my($class, $typename, $c) = @_;
    my $type = $types{$typename} // die "Unknown typename '$typename'";
    eval qq{use $type; 1} or die $@;
    my $self = bless {
        c => $c,
        n => $c->n,
        f => $c->f,
        min => $c->min,
    }, $type;
    $self->init;
    return $self;
}

sub check_exceptions {
    return;
}

sub check_fixed {
    return;
}

for my $method (qw{
    name func_value func_name func func_target
    apply_m to_test test_target
}) {
    no strict 'refs';
    *$method = sub { die "$_[0] must provide implementation of sub $method" };
}

sub c { $_[0]->{c} }

for my $method (qw{ n f min }) {
    no strict 'refs';
    *$method = sub { $_[0]->{$method} };
}

1;
