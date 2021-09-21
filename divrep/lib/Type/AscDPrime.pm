package Type::AscDPrime;
use strict;
use warnings;

use parent qw{ Type };
use Math::Prime::Util qw{ factor_exp gcd };
use Memoize;

=head2 Type::AscDPrime

dprime(d + kn) = k + 1

=cut

sub init { return }

# Count number of distinct primes dividing n, given factorization
sub dprimes { scalar @{ $_[0] } }

sub name { 'ascdp' }
sub dbname { 'ascdp' }

sub smallest { 1 }

memoize('gprio', NORMALIZER => sub { "$_[1]" });
sub gprio {
    # we're really only interested in n=1 here (A086560)
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

sub func_name { 'dprimes' }

sub func { dprimes([ factor_exp($_[1]) ]) }

sub func_target { $_[1] + 1 }

sub apply_m {
    my($self, $m, $fm) = @_;
    my $v = dprimes($fm);

    # if m has v distinct primes, it cannot divide d + kn: k < v - 1
    my $c = $self->c;
    my $n = $self->n;
    my $lim = $self->f - 1;
    $lim = $v - 2 if $lim > $v - 2;
    for (0 .. $lim) {
        my $off = -$n * $_;
        $c->suppress($m, $off % $m, 0);
    }

    # Let r = rad(m). If r | d + (v - 1)n then (d + (v - 1)n) / r
    # can be divisible only by the primes that divide r, so
    # gcd((d + (v - 1)n) / r, r) != 1.
    if ($v - 1 <= $self->f - 1) {
        my $r = 1;
        $r *= $_->[0] for @$fm;
        if ($m % ($r * $r) == 0) {
            my $mdr = $m / $r;
            for my $w (1 .. $r - 1) {
                next if gcd($w, $r) > 1;
                my $off = $mdr * $w - $n * ($v - 1);
                $c->suppress($m, $off % $m, $off);
            }
        }
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
    return [ "$k", sub { dprimes([ factor_exp($_[0] + $n * $k) ]) == $k + 1 } ];
}

1;
