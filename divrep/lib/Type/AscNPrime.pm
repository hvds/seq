package Type::AscNPrime;
use strict;
use warnings;

use parent qw{ Type };
use Math::Prime::Util qw{ factor_exp next_prime };
use Memoize;

=head2 Type::AscNPrime

nprime(d + kn) = k + 1

=cut

sub init { return }

# Count number of primes dividing n (with multiplicity), given factorization
sub nprimes {
    my($px) = @_;
    my $sum = 0;
    $sum += $_->[1] for @$px;
    return $sum;
}

sub name { 'ascnp' }
sub dbname { 'ascnp' }

sub smallest { 1 }

memoize('gprio', NORMALIZER => sub { "$_[1]" });
sub gprio {
    # we're really only interested in n=1 here (A072875)
    my($self, $n) = @_;
    return -(
        (log($n) / log(2)) ** 10
    );
}

sub ming { 1 }

sub maxg {
    my($self, $n) = @_;
    die "cannot maxg(0)" unless $n;
    # Unknown limit, except in terms of what we can practically reach
    100;
}

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $d + $k * $n;
}

sub func_name { 'nprimes' }

sub func { nprimes([ factor_exp($_[1]) ]) }

sub func_target { $_[1] + 1 }

sub apply_m {
    my($self, $m, $fm) = @_;
    my $v = nprimes($fm);

    # if m has v primes, it cannot divide d + kn: k < v - 1, and cannot divide
    # d + (v-1)n if greater than m
    my $c = $self->c;
    my $n = $self->n;
    my $lim = $self->f - 1;
    $lim = $v - 1 if $lim > $v - 1;
    for (0 .. $lim) {
        my $off = -$n * $_;
        my $min = ($_ == $v - 1) ? $m + $off : 0;
        $c->suppress($m, $off % $m, $min);
    }
    return;
}

sub to_testf {
    my($self, $f) = @_;
    return [ 0 .. $f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $n = $self->n;
    return [ "$k", sub { nprimes([ factor_exp($_[0] + $n * $k) ]) == $k + 1 } ];
}

1;
