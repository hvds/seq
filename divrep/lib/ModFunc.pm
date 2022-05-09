package ModFunc;
use strict;

use Math::GMP ();
use Math::Prime::Util qw{ sqrtmod chinese factor_exp };
use List::Util qw{ reduce };

require Exporter;
our @ISA = qw/ Exporter /;
our @EXPORT_OK = qw/
    quadres is_residue quadvec
    mod_combine
    gcd
    allsqrtmod allsqrtmodfact allsqrtmodpk
/;

my $zero = Math::GMP->new(0);
sub MBI { defined($_[0]) ? Math::GMP->new(@_) : undef }

{
    my %q;

    #
    # Returns Legendre ( n / m ), which is 0 if gcd(n, m) > 1,
    # 1 if n is a quadratic residue (mod m), -1 otherwise.
    #
    sub quadres {
        my($n, $m) = @_;
        return 0 if gcd($n, $m) > 1;
        my $v = quadvec($m);
        return vec($v, ($n % $m), 1) ? 1 : -1;
    }

    #
    # Returns true if Legendre ( n / m ) is 0 or 1.
    #
    sub is_residue {
        my($n, $m) = @_;
        my $v = quadvec($m);
        return vec($v, ($n % $m), 1);
    }

    #
    # Returns an m-bit vector, such that bit i is true iff is_residue(i, m).
    #
    sub quadvec {
        my($m) = @_;
        ${ $q{$m} ||= do {
            my $v = "";
            vec($v, (($_ * $_) % $m), 1) = 1 for 0 .. $m - 1;
            \$v;
        } };
    }
}


#
# Given constraints (a, b) and (c, d) representing n=a(mod b), n=c(mod d),
# return (e, f) such that n=e(mod f) encapsulates both constraints.
#
sub mod_combine {
    my($a, $b, $c, $d) = @_;
    my $g = gcd($b, $d);
    if ($g > 1) {
        die <<DIE unless ($a % $g) == ($c % $g);
Inconsistent values: can't combine $a(mod $b) with $c(mod $d)
DIE
        while ($g > 1) {
            my $gb = gcd($b / $g, $g);
            my $gd = gcd($d / $g, $g);
            if ($gb > 1) {
                # p|gb constrained to higher power in (a, b) than in (c, d)
                $d /= $gb;
                $c %= $d;
            } elsif ($gd > 1) {
                # p|gd constrained to higher power in (c, d) than in (a, b)
                $b /= $gd;
                $a %= $b;
            } else {
                # p|g constrained to same power in both; pick one
                $b /= $g;
                $a %= $b;
            }
            $g = gcd($b, $d);
        }
    }
    # with the gcd reduced, we can safely take inverses
    my $bi = $b->bmodinv($d);
    my $di = $d->bmodinv($b);
    my $ab = (($a * $di) % $b) * $d;
    my $cd = (($c * $bi) % $d) * $b;
    my $bd = $b * $d;
    return (($ab + $cd) % $bd, $bd);
}

sub gcd {
    my $n0 = shift;
    $n0 = Math::GMP->new($n0) unless ref $n0;
    return reduce { $a->bgcd($b) } ($n0, @_);
}

=head2 allsqrtmodpk( $a, $p, $k )

Given integer C<a>, prime C<p> and positive integer C<k>, return a list of
the square roots C<x> of C<a mod p^k> having C<< 0 <= x < p^k >>. If no
square root exists, an empty list is returned.

=cut

sub allsqrtmodpk {
    my($a, $p, $k) = @_;
    my $pk = $p ** $k;
    unless ($a % $p) {
        unless ($a % $pk) {
            # if p^k divides a, we need the square roots of zero, satisfied by
            # ip^j with 0 <= i < p^{floor(k/2)}, j = p^{ceil(k/2)}
            my $low = $p ** ($k >> 1);
            my $high = ($k & 1) ? ($low * $p) : $low;
            return map $high * $_, 0 .. $low - 1;
        }
        my $a2 = $a / $p;
        return () if $a2 % $p;      # p divides a, p^2 does not
        my $pj = $pk / $p;
        return map {
            my $qp = $_ * $p;
            map $qp + $_ * $pj, 0 .. $p - 1;
        } allsqrtmodpk($a2 / $p, $p, $k - 2);
    }
    my $q = MBI(sqrtmod($a, $pk));
    return () unless defined $q;
    return +($q, $pk - $q) if $p != 2;
    return +($q) if $k == 1;
    return +($q, $pk - $q) if $k == 2;
    my $pj = $pk / $p;
    my $q2 = ($q * ($pj - 1)) % $pk;
    return +($q, $pk - $q, $q2, $pk - $q2);
}

=head2 allsqrtmodfact( $a, $n, $f )

Given integers C<a> and C<n>, and the factorization C<f> of C<n> as returned
by C<factor_exp(n)>, return a list of the square roots C<x> of C<a mod n>
having C<< 0 <= x < n >>. If no square root exists, an empty list is returned.

=cut

sub allsqrtmodfact {
    my($a, $n, $f) = @_;
    return +($zero) unless @$f;
    my($p, $k) = @{ $f->[0] };
    $p = MBI($p);
    my @q = allsqrtmodpk($a, $p, $k);
    return @q if @$f == 1;
    my $pk = $p ** $k;
    my $n2 = $n / $pk;
    return map {
        my $q2 = $_;
        map MBI(chinese([ $q2, $n2 ], [ $_, $pk ])), @q;
    } allsqrtmodfact($a, $n2, [ @$f[1 .. $#$f] ]);
}

=head2 allsqrtmod( $a, $n )

Given integers C<a> and C<n>, return a list of the square roots C<x> of
C<a mod n> having C<< 0 <= x < n >>. If no square root exists, an empty
list is returned.

=cut

sub allsqrtmod {
    my($a, $n) = @_;
    return allsqrtmodfact($a, $n, [ factor_exp($n) ]);
}

1;
