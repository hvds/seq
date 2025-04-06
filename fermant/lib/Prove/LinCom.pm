package Prove::LinCom;
use strict;
use warnings;
use Math::BigRat;

=head1 NAME

Prove::LinCom - (min, max) pairs of linear combinations of variables

=head1 DESCRIPTION

Blessed arrayref of C< [min, max] > pairs, where each of min, max is
a linear combination of C< (1, a, b, ...) > as an arrayref of integers.

=cut

sub new {
    my($class, $data) = @_;
    return bless $data, $class;
}

sub new_splice {
    my($class, $vars, $vals) = @_;
    my @n;
    @n[@$vars] = @$vals;
    $_ //= 0 for @n;
    return $class->new(\@n);
}

sub new_clean {
    my($self, $data) = @_;
    pop @$data while @$data && $data->[-1] == 0;
    return bless $data, ref($self);
}

sub copy {
    return ref($_[0])->new([ @{ $_[0] } ]);
}

sub negate {
    return ref($_[0])->new([ map -$_, @{ $_[0] } ]);
}

sub is_const {
    my($self) = @_;
    return @$self < 2 ? 1 : 0;
}

sub const {
    my($self) = @_;
    die if @$self > 1;
    return $self->[0] // 0;
}

sub popone {
    my($self) = @_;
    my @n = @$self[0 .. $#$self - 1];
    return $self->new_clean(\@n);
}

sub size {
    my($self) = @_;
    return @$self ? $#$self : 0;
}

sub top {
    my($self) = @_;
    die if @$self < 2;
    return +($#$self, $self->[-1]);
}

sub all_pair {
    my($self) = @_;
    return @$self > 1
        ? map([ $_, $self->[$_] ], grep $self->[$_], 0 .. $#$self)
        : [ 0, $self->[0] // 0 ];
}

sub addi {
    my($self, $i, $c) = @_;
    my @n = @$self;
    push @n, 0 while $#n < $i;
    $n[$i] += $c;
    return $self->new_clean(\@n);
}

sub divc {
    my($self, $c) = @_;
    my $mul = Math::BigRat->new(1) / $c;
    return $self->new_clean([ map $_ * $mul, @$self ]);
}

sub addmul {
    my($self, $other, $c) = @_;
    #die if @$other > (@$self || 1);
    my @n = @$self;
    $n[$_] += $other->[$_] * $c for 0 .. $#$other;
    return $self->new_clean(\@n);
}

sub parse {
    my($class, $e, $vi) = @_;
    return $class->new(_parse($e, $vi));
}

sub _parse {
    my($e, $vi) = @_;
    my $t = $e->type;
    my @a;
    if ($t eq 'pluslist') {
        for (@{ $e->args }) {
            my $b = _parse($_, $vi);
            $a[$_] += $b->[$_] // 0 for 0 .. $#$b;
        }
    } elsif ($t eq 'negate') {
        my $b = _parse($e->args->[0], $vi);
        $a[$_] -= $b->[$_] // 0 for 0 .. $#$b;
    } elsif ($t eq 'name') {
        my $i = $vi->{$e->name} // die "Unknown name @{[ $e->name ]}";
        $a[$i] = 1;
    } elsif ($t eq 'integer') {
        $a[0] += $e->rat->numify;
    } elsif ($t eq 'mullist') {
        my $args = $e->args;
        if (@$args == 2 && $args->[0]->is_const) {
            my $c = $args->[0]->rat->numify;
            my $b = _parse($e->args->[1], $vi);
            $a[$_] += $c * ($b->[$_] // 0) for 0 .. $#$b;
        } else {
            die "Unexpected mullist @{[ $e->bracketed ]}\n";
        }
    } else {
        die "Unexpected type '$t' in @{[ $e->bracketed ]}\n";
    }
    return \@a;
}

1;
