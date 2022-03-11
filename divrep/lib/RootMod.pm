package RootMod;
use strict;
use warnings;

require Exporter;
our @ISA = qw{ Exporter };
our @EXPORT_OK = qw{
    allrootmod countrootmod
};

use Math::Prime::Util qw{
    factor factor_exp is_prime
    addmod mulmod divmod powmod invmod
    chinese gcd kronecker valuation
};

*allrootmod = Math::Prime::Util->can('allrootmod') // \&_allrootmod;

BEGIN {
  use constant MPU_MAXBITS     => (~0 == 4294967295) ? 32 : 64;
}

use Math::GMP ();
my $zone = Math::GMP->new(1);

sub countrootmod {
    my($A, $k, $n) = @_;
    return scalar @{[ allrootmod($A, $k, $n) ]};
}

sub submod {
    my($a, $b, $n) = @_;
    return addmod($a, -$b, $n);
}

# Tonelli-Shanks
sub _sqrtmod_prime {
    my($a, $p) = @_;
    my($x, $q, $e, $t, $z, $r, $m, $b);
    my $Q = $p - 1;

    if (($p % 4) == 3) {
        $r = powmod($a, ($p + 1) >> 2, $p);
        return undef unless mulmod($r, $r, $p) == $a;
        return $r;
    }
    if (($p % 8) == 5) {
        $m = addmod($a, $a, $p);
        $t = powmod($m, ($p - 5) >> 3, $p);
        $z = mulmod($m, mulmod($t, $t, $p), $p);
        $r = mulmod($t, mulmod($a, submod($z, 1, $p), $p), $p);
        return undef unless mulmod($r, $r, $p) == $a;
        return $r;
    }

    # Verify Euler's criterion for odd p
    return undef if $p != 2 && powmod($a, $Q >> 1, $p) != 1;

    # Cohen Algorithm 1.5.1.  Tonelli-Shanks.
    $e = valuation($Q, 2);
    $q = $Q / ($zone << $e);
    $t = 3;
    while (kronecker($t, $p) != -1) {
        $t += 2;
        return undef if $t == 201 && !is_prime($p);
    }
    $z = powmod($t, $q, $p);
    $b = powmod($a, $q, $p);
    $r = $e;
    $q = ($q + 1) >> 1;
    $x = powmod($a, $q, $p);
    while ($b != 1) {
        $t = $b;
        for ($m = 0; $m < $r && $t != 1; $m++) {
            $t = mulmod($t, $t, $p);
        }
        $t = powmod($z, $zone << ($r - $m - 1), $p);
        $x = mulmod($x, $t, $p);
        $z = mulmod($t, $t, $p);
        $b = mulmod($b, $z, $p);
        $r = $m;
    }
    # Expected to always be true.
    return undef unless mulmod($x, $x, $p) == $a;
    return $x;
}

sub _sqrtmod_prime_power {
    my($a, $p, $e) = @_;
    my($r, $s);

    if ($e == 1) {
        $a %= $p if $a >= $p;
        return $a if $p == 2 || $a == 0;
        $r = _sqrtmod_prime($a, $p);
        return (defined $r && (mulmod($r, $r, $p) == $a) ? $r : undef);
    }

    my $n = $p ** $e;
    my $pk = $p * $p;

    return 0 if ($a % $n) == 0;

    if (($a % $pk) == 0) {
        my $apk = $a / $pk;
        $s = _sqrtmod_prime_power($apk, $p, $e-2);
        return undef unless defined $s;
        return $s * $p;
    }

    return undef if ($a % $p) == 0;

    my $ered = ($p > 2 || $e < 5) ? (($e + 1) >> 1) : (($e + 3) >> 1);
    $s = _sqrtmod_prime_power($a, $p, $ered);
    return undef unless defined $s;

    my $np = ($p == 2) ? ($n * $p) : $n;
    my $t1 = submod($a, mulmod($s, $s, $np), $np);
    my $t2 = addmod($s, $s, $np);
    my $gcd = gcd($t1, $t2);
    $r = addmod($s, divmod($t1 / $gcd, $t2 / $gcd, $n), $n);
    return ((mulmod($r, $r, $n) == ($a % $n)) ? $r : undef);
}

# helper function for allsqrtmod() - return list of all square roots of
# a (mod p^k), assuming a integer, p prime, k positive integer.
sub _allsqrtmodpk {
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
        # p divides a, p^2 does not
        my $a2 = $a / $p;
        return () if $a2 % $p;
        my $pj = $pk / $p;
        return map {
            my $qp = $_ * $p;
            map $qp + $_ * $pj, 0 .. $p - 1;
        } _allsqrtmodpk($a2 / $p, $p, $k - 2);
    }
    my $q = _sqrtmod_prime_power($a, $p, $k);
    return () unless defined $q;
    return ($q, $pk - $q) if $p != 2;
    return ($q) if $k == 1;
    return ($q, $pk - $q) if $k == 2;
    my $pj = $pk / $p;
    my $q2 = ($q * ($pj - 1)) % $pk;
    return ($q, $pk - $q, $q2, $pk - $q2);
}

# helper function for allsqrtmod() - return list of all square roots of
# a (mod p^k), assuming a integer, n positive integer > 1, f arrayref
# of [ p, k ] pairs representing factorization of n. Destroys f.
sub _allsqrtmodfact {
    my($a, $n, $f) = @_;
    my($p, $k) = @{ shift @$f };
    my @q = _allsqrtmodpk($a, $p, $k);
    return @q unless @$f;
    my $pk = $p ** $k;
    my $n2 = $n / $pk;
    return map {
        my $q2 = $_;
        map chinese([ $q2, $n2 ], [ $_, $pk ]), @q;
    } _allsqrtmodfact($a, $n2, $f);
}

sub _ts_prime {
    my($a, $k, $p, $refzeta) = @_;
    my($e, $r) = (0, $p - 1);
    while (!($r % $k)) {
        $e++;
        $r /= $k;
    }
    my $ke = ($p - 1) / $r;
    my $x = powmod($a, invmod($k % $r, $r), $p);
    my $B = mulmod(powmod($x, $k, $p), invmod($a, $p), $p);
    my($T, $y, $t, $A) = (2, 1);
    while ($y == 1) {
        $t = powmod($T, $r, $p);
        $y = powmod($t, $ke / $k, $p);
        $T++;
    }
    while ($ke != $k) {
        $ke /= $k;
        $T = $t;
        $t = powmod($t, $k, $p);
        $A = powmod($B, $ke / $k, $p);
        while ($A != 1) {
            $x = mulmod($x, $T, $p);
            $B = mulmod($B, $t, $p);
            $A = mulmod($A, $y, $p);
        }
    }
    $$refzeta = $t if defined $refzeta;
    $x;
}

sub _allrootmod_cprod {
    my($aroots1, $p1, $aroots2, $p2) = @_;
    my($t, $n, $inv);

    $n = $p1 * $p2;
    $inv = invmod($p1, $p2);
    die("CRT has undefined inverse") unless defined $inv;

    my @roots;
    for my $q1 (@$aroots1) {
        for my $q2 (@$aroots2) {
            $t = mulmod($inv, submod($q2, $q1, $p2), $p2);
            $t = addmod($q1, mulmod($p1, $t, $n), $n);
            push @roots, $t;
        }
    }
    return @roots;
}

sub _allrootmod_prime {
    my($a, $k, $p) = @_;        # prime k, prime p
    $a %= $p if $a >= $p;

    return ($a) if $p == 2 || $a == 0;

    # If co-prime, there is exactly one root.
    my $g = gcd($k, $p - 1);
    if ($g == 1) {
        my $r = powmod($a, invmod($k % ($p - 1), $p - 1), $p);
        return ($r);
    }

    # Check generalized Euler's criterion
    return () if powmod($a, ($p - 1) / $g, $p) != 1;

    # Special case for p=3 for performance
    return (1, 2) if $p == 3;

    # Call one of the general TS solvers to get all the roots.
    my $z;
    my $r = _ts_prime($a, $k, $p, \$z);
    die "allrootmod: failed to find root"
              if $z == 0 || powmod($r, $k, $p) != $a;
    my @roots = ($r);
    my $r2 = mulmod($r, $z, $p);
    while ($r2 != $r && @roots < $k) {
        push @roots, $r2;
        $r2 = mulmod($r2, $z, $p);
    }
    die "allrootmod: excess roots found" if $r2 != $r;
    return @roots;
}

sub _allrootmod_prime_power {
    my($a, $k, $p, $e) = @_;        # prime k, prime p
    return _allrootmod_prime($a, $k, $p) if $e == 1;
    my $n = $p ** $e;
    my $pk = $p ** $k;
    my @roots;

    if (($a % $n) == 0) {
        my $t = ($e - 1) / $k + 1;
        my $nt = $p ** $t;
        my $nr = $p ** ($e - $t);
        @roots = map { mulmod($_, $nt, $n) } 0 .. $nr-1;
        return @roots;
    }

    if (($a % $pk) == 0) {
        my $apk = $a / $pk;
        my $pe1 = $p ** ($k - 1);
        my $pek = $p ** ($e - $k + 1);
        my @roots2 = _allrootmod_prime_power($apk, $k, $p, $e - $k);
        for my $r (@roots2) {
            my $rp = mulmod($r, $p, $n);
            for my $j (0 .. $pe1 - 1) {
                push @roots, addmod($rp, mulmod($j, $pek, $n), $n);
            }
        }
        return @roots;
    }

    return () if ($a % $p) == 0;

    my $np = $n * $p;
    my $ered = ($p > 2 || $e < 5)
        ? (($e + 1) >> 1)
        : (($e + 3) >> 1);
    my @roots2 = _allrootmod_prime_power($a, $k, $p, $ered);

    if ($k != $p) {
        for my $s (@roots2) {
            my $t = powmod($s, $k - 1, $n);
            my $t1 = submod($a, mulmod($t, $s, $n), $n);
            my $t2 = mulmod($k, $t, $n);
            my $gcd = gcd($t1, $t2);
            my $r = addmod($s, divmod($t1 / $gcd, $t2 / $gcd, $n), $n);
            push @roots, $r;
        }
    } else {
        my @rootst;
        for my $s (@roots2) {
            my $t = powmod($s, $k - 1, $np);
            my $t1 = submod($a, mulmod($t, $s, $np), $np);
            my $t2 = mulmod($k, $t, $np);
            my $gcd = gcd($t1, $t2);
            my $r = addmod($s, divmod($t1 / $gcd, $t2 / $gcd, $n), $n);
            push @rootst, $r if powmod($r, $k, $n) == ($a % $n);
        }
        my $ndivp = $n / $p;
        my %roots;  # We want to remove duplicates
        for my $r (@rootst) {
            for my $j (0 .. $k-1) {
                $roots{ mulmod($r, addmod(mulmod($j, $ndivp, $n), 1, $n), $n) } = undef;
            }
        }
        @roots = keys(%roots);
    }
    return @roots;
}

sub _allrootmod_kprime {
    my($a, $k, $n, @nf) = @_;       # k prime, n factored into f^e,f^e,...

    return _allsqrtmodfact($a, $n, \@nf) if $k == 2;

    my $N = 1;
    my @roots;
    foreach my $F (@nf) {
        my($f,$e) = @$F;
        my $fe = $f ** $e;
        my @roots2 = ($e == 1)
            ? _allrootmod_prime($a, $k, $f)
            : _allrootmod_prime_power($a, $k, $f, $e);
        return () unless @roots2;
        if (@roots == 0) {
            @roots = @roots2;
        } else {
            @roots = _allrootmod_cprod(\@roots, $N, \@roots2, $fe);
        }
        $N *= $fe;
    }

    return @roots;
}

sub _allrootmod {
    my($A, $k, $n) = map {
        UNIVERSAL::isa($_, 'Math::GMP') ? $_ : Math::GMP->new($_)
    } @_;
    $n = -$n if $n < 0;

    return () if $n == 0;
    $A = $A % $n;

    return () if $k <= 0 && $A == 0;

    if ($k < 0) {
        $A = invmod($A, $n);
        return () unless defined $A && $A > 0;
        $k = -$k;
    }

    return ($A) if $n <= 2 || $k == 1;
    return ($A == 1) ? (0 .. $n-1) : () if $k == 0;

    my @roots;
    my @nf = factor_exp($n);

    if (is_prime($k)) {
        @roots = _allrootmod_kprime($A, $k, $n, @nf);
    } else {
        @roots = ($A);
        for my $primek (factor($k)) {
            my @rootsnew = ();
            for my $r (@roots) {
                push @rootsnew, _allrootmod_kprime($r, $primek, $n, @nf);
            }
            @roots = @rootsnew;
        }
    }

    @roots = sort { $a <=> $b } @roots;
    return @roots;
}

1;
