package Type::Semip;
use strict;
use warnings;

use parent qw{ Type };
use Math::Prime::Util qw{ factor_exp };

sub init { return }

# Count number of primes dividing n (with multiplicity), given factorization
sub nprimes {
    my($px) = @_;
    my $sum = 0;
    $sum += $_->[1] for @$px;
    return $sum;
}

sub name { 'semip' }

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $d + $k * $n;
}

sub func_name { 'nprimes' }

sub func { nprimes([ factor_exp($_[1]) ]) }

sub func_target { 2 }

# No result > m can be divisible by m if m is a semiprime
sub apply_m {
    my($self, $m, $fm) = @_;
    if (nprimes($fm) == 2) {
        my $c = $self->c;
        my $n = $self->n;
        for (0 .. $self->f - 1) {
            my $off = -$n * $_;
            $c->suppress($m, $off % $m, $m + $off);
        }
    }
    return;
}

sub to_test {
    my($self) = @_;
    return [ 0 .. $self->f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $n = $self->n;
    return [ "$k", sub { nprimes([ factor_exp($_[0] + $n * $k) ]) == 2 } ];
}

1;
