package ModFunc;
use strict;

use Math::GMP ();
use Math::Prime::Util qw{
    factor factor_exp chinese is_prime
    sqrtmod invmod divmod
};
use List::Util qw{ reduce };

# allsqrtmod() and allrootmod() and associated functions adapted from
# implementations in Math::Prime::Util::PP

require Exporter;
our @ISA = qw{ Exporter };
our @EXPORT_OK = qw{
    quadres is_residue quadvec
    mod_combine
    gcd
    allsqrtmod allsqrtmodfact allsqrtmodpk
    allrootmod allrootmod_kprime allrootmod_prime allrootmod_prime_power
      allrootmod_cprod ts_prime
};

sub _numeric { $a <=> $b }
sub MBI { defined($_[0]) ? Math::GMP->new(@_) : undef }
my $zero = MBI(0);
my $zone = MBI(1);
my $ztwo = MBI(2);

sub powmod { MBI(&Math::Prime::Util::powmod(@_)) }

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
    return sort _numeric @q if @$f == 1;
    my $pk = $p ** $k;
    my $n2 = $n / $pk;
    return sort _numeric map {
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
    return () if $n == 0;
    return allsqrtmodfact($a, $n, [ factor_exp($n) ]);
}

=head2 ts_prime($a, $k, $p, $refzeta)

"Tonelli-Shanks kth roots  alternate version"

=cut

sub ts_prime {
    my($a, $k, $p, $refzeta) = @_;
    my($e, $r) = ($zero, $p - 1);
    ++$e, $r /= $k while !($r % $k);

    my $ke = ($p - 1) / $r;
    my $x = powmod($a, MBI(invmod($k % $r, $r)), $p);
    my $B = (powmod($x, $k, $p) * MBI(invmod($a, $p))) % $p;
    my($T, $y, $t, $A) = ($ztwo, $zone);

    while ($y == 1) {
        $t = powmod($T, $r, $p);
        $y = powmod($t, $ke / $k, $p);
        ++$T;
    }

    while ($ke != $k) {
        $ke /= $k;
        $T = $t;
        $t = powmod($t, $k, $p);
        $A = powmod($B, $ke / $k, $p);
        while ($A != 1) {
            $x = ($x * $T) % $p;
            $B = ($B * $t) % $p;
            $A = ($A * $y) % $p;
        }
    }
    $$refzeta = $t;
    return $x;
}

=head2 allrootmod_cprod(\@r1, $n1, \@r2, $n2)

Given coprime integers C<n1>, C<n2> and arrayrefs C<r1>, C<r2>, such that
elements of C<r1> and C<r2> are C<k>'th roots of some C<a> mod C<n1> and
mod C<n2> respectively, returns a new list of the C<k>'th roots of C<a>
mod C<n1 * n2>.

=cut

sub allrootmod_cprod {
    my($r1, $n1, $r2, $n2) = @_;
    my $n = $n1 * $n2;
    my $inv = MBI(invmod($n1, $n2))
            // die "allrootmod_cprod() has undefined inverse";
    return map {
        my $q1 = $_;
        map +($q1 + $n1 * (($inv * ($_ - $q1)) % $n2)) % $n, @$r2;
    } @$r1;
}

=head2 allrootmod_prime($a, $k, $p)

Given integer C<a> and prime C<k> and C<p>, returns a list of the C<k>'th
roots C<x> of C<a mod p> having C<< 0 <= x < p >>. If no C<k>'th root
exists, an empty list is returned.

=cut

sub allrootmod_prime {
    my($a, $k, $p) = @_;
    $a %= $p;
    return ($a) if $p == 2 || $a == 0;

    my $pm = $p - 1;
    my $g = gcd($k, $pm);

    # If co-prime, there is exactly one root.
    return powmod($a, MBI(invmod($k % $pm, $pm)), $p)
            if $g == 1;

    # Check generalized Euler's criterion
    return () if powmod($a, $pm / $g, $p) != 1;

    # Special case for p=3 for performance
    return ($zone, $ztwo) if $p == 3;

    # Call a Tonelli-Shanks solver that also allows us to get all the roots
    my $z;
    my $r = MBI(ts_prime($a, $k, $p, \$z));
    die "allrootmod: failed to find root"
            if $z == 0 || powmod($r, $k, $p) != $a;
    $z = MBI($z);
    my @roots = ($r);
    my $r2 = ($r * $z) % $p;
    while ($r2 != $r && @roots < $k) {
        push @roots, $r2;
        $r2 = ($r2 * $z) % $p;
    }
    die "allrootmod: excess roots found"
            if $r2 != $r;
    return @roots;
}

=head2 allrootmod_prime_power($a, $k, $p, $x)

Given integers C<a>, C<x> and prime C<k>, C<p>, returns a list of the C<k>'th
roots C<r> of C<a mod p^x> having C<< 0 <= r < p^x >>. If no C<k>'th root
exists, an empty list is returned.

=cut

sub allrootmod_prime_power {
    my($a, $k, $p, $x) = @_;
    return allrootmod_prime($a, $k, $p) if $x == 1;
    my $px = $p ** $x;
    my $pk = $p ** $k;

    $a %= $px;
    if ($a == 0) {
        my $t = ($x - 1) / $k + 1;
        my $pt = $p ** $t;
        my $pr = $p ** ($x - $t);
        return map +(($_ * $pt) % $px), 0 .. $pr - 1;
    }

    if (($a % $pk) == 0) {
        my $apk = $a / $pk;
        my $pe1 = $p ** ($k - 1);
        my $pek = $p ** ($x - $k + 1);
        return map {
            my $rp = ($_ * $p) % $px;
            map +($_ * $pek + $rp) % $px, 0 .. $pe1 - 1;
        } allrootmod_prime_power($apk, $k, $p, $x - $k);
    }
    return () if ($a % $p) == 0;

    my $pxp = $px * $p;
    my $ered = ($p > 2 || $x < 5)
        ? ($x + 1) >> 1
        : ($x + 3) >> 1;
    my @r = allrootmod_prime_power($a, $k, $p, $ered);

    if ($k != $p) {
        return map {
            my $s = $_;
            my $t = powmod($s, $k - 1, $px);
            my $t1 = ($a - powmod($t, $s, $px)) % $px;
            my $t2 = ($k * $t) % $px;
            my $g = gcd($t1, $t2);
            ($s + MBI(divmod($t1 / $g, $t2 / $g, $px))) % $px;
        } @r;
    }

    my @r2;
    for my $s (@r) {
        my $t = powmod($s, $k - 1, $pxp);
        my $t1 = ($a - $t * $s) % $pxp;
        my $t2 = ($k * $t) % $pxp;
        my $g = gcd($t1, $t2);
        my $r = ($s + MBI(divmod($t1 / $g, $t2 / $g, $px))) % $px;
        push @r2, $r if powmod($r, $k, $px) == $a;
    }
    my $pxm = $px / $p;
    my %r;
    for my $r (@r2) {
        for my $j (0 .. $k - 1) {
            $r{($r * ($j * $pxm + 1)) % $px} = undef;
        }
    }
    return map MBI($_), keys %r;
}

=head2 allrootmod_kprime( $a, $k, $n, $nf )

Given integers C<a>, C<k> and C<n> and factorization C<nf> with C<k>
prime, return a list of the C<k>'th roots C<x> of C<a mod n> having
C<< 0 <= x < n >>. If no C<k>'th root exists, an empty list is returned.

=cut

sub allrootmod_kprime {
    my($a, $k, $n, $nf) = @_;
    return allsqrtmodfact($a, $n, $nf) if $k == 2;

    my $m = $zone;
    my @roots;
    for (@$nf) {
        my($p, $x) = @$_;
        $p = MBI($p);
        my @r2 = ($x == 1)
            ? allrootmod_prime($a, $k, $p)
            : allrootmod_prime_power($a, $k, $p, $x);
        return () unless @r2;
        my $px = $p ** $x;
        @roots = @roots
            ? allrootmod_cprod(\@roots, $m, \@r2, $px)
            : @r2;
        $m *= $px;
    }
    return @roots;
}

=head2 allrootmod( $a, $k, $n )

Given integers C<a>, C<k> and C<n>, return a list of the C<k>'th roots
C<x> of C<a mod n> having C<< 0 <= x < n >>. If no C<k>'th root exists,
an empty list is returned.

=cut

sub allrootmod {
    my($a, $k, $n) = @_;
    return () if $n == 0;
    $a %= $n;
    return () if $k <= 0 && $a == 0;
    if ($k < 0) {
        $a = MBI(invmod($a, $n));
        return () unless defined($a) && $a > 0;
        $k = -$k;
    }
    return ($a) if $n <= 2 || $k == 1;
    return +($a == 1) ? map(MBI($_), (0 .. $n - 1)) : () if $k == 0;

    my $nf = [ factor_exp($n) ];
    if (is_prime($k)) {
        return sort _numeric allrootmod_kprime($a, $k, $n, $nf);
    }

    my @roots = ($a);
    for my $kp (factor($k)) {
        @roots = map allrootmod_kprime($_, $kp, $n, $nf), @roots;
    }
    return sort _numeric @roots;
}

1;
