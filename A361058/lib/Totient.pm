package Totient;

use strict;
use warnings;
use Math::Prime::Util qw{ factor_exp is_prime divisors };
use Math::GMP;

=head1 Totient::totient

Return a list of integers x such that phi(x) = n.

Results are returned as strings, so they can easily be instantiated as
bigints if needed without loss of precision.

=cut

our $DEBUG = 0;
my $ztwo = Math::GMP->new(2);

# Assume the 5 known Fermat primes are the only ones that exist in any
# range we'll be asked to deal with.
my @fermat = (
    [ 1, 3 ], [ 2, 5 ], [ 4, 17 ], [ 8, 257 ], [ 16, 65537 ],
);
sub totient {
    my($n) = @_;
    return (1, 2) if $n == 1;
    return () if $n & 1;
    return (3, 4, 6) if $n == 2;

    my @stack = ([ $n, 1, {}]);
    my %result;
    while (@stack) {
        my($n, $mult, $used) = @{ pop @stack };
        if ($n == 1) {
            $result{$mult * $_} = 1 for (1, 2);
            next;
        }
        next if $n & 1;
        if ($n == 2) {
            $result{$mult * $_} = 1 for ($used->{3} ? (4) : (3, 4, 6));
            next;
        }

        my $fn = [ factor_exp($n) ];
        # first try using any odd prime that divides n
        for my $i (reverse 1 .. $#$fn) {
            my($pi, $ei) = @{ $fn->[$i] };
            next if $used->{$pi};
            $used->{$pi} = 1;
            next if $n % ($pi - 1);
            my $xu = Math::GMP->new($pi);
            my $nu = $n / ($pi - 1);
            for my $eu (0 .. $ei) {
                push @stack, [ $nu, $xu * $mult, { %$used } ];
                $xu *= $pi;
                $nu /= $pi;
            }
        }

        # any remaining results must be squarefree other than powers of 2
        my $np1 = $n + 1;
        if (!$used->{$np1} && is_prime($np1)) {
            $result{$mult * $_} = 1 for ($np1, $np1 * 2);
        }

        # any remaining result is also composite
        # if n is 2 mod 4, no further solutions possible
        my $e2 = $fn->[0][1];
        next if $e2 == 1;

        if (@$fn > 1) {
            my $base = $ztwo * $fn->[-1][0];
            my $rest = $n / $base;
            for my $d (divisors($rest)) {
                last if $d == $rest;
                my $pi = $base * $d + 1;
                next if $used->{$pi} || !is_prime($pi);
                $used->{$pi} = 1;
                push @stack, [ $rest / $d, $pi * $mult, { %$used } ];
            }
            next;
        }

        # what's left is 2^e2; assume only known fermat primes are relevant
        for (@fermat) {
            my($pow, $prime) = @$_;
            last if $pow > $e2;
            next if $used->{$prime};
            $used->{$prime} = 1;
            push @stack, [ $n / ($prime - 1), $prime * $mult, { %$used } ];
        }
        $result{$mult * $_} = 1 for ($n * 2);
        next;
    }
    return sort { $a <=> $b } keys %result;
}

1;
