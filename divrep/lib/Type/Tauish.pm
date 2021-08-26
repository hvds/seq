package Type::Tauish;
use strict;
use warnings;

use parent qw{ Type };
use Math::GMP;
use Math::Prime::Util qw{ factor_exp };
use Memoize;

use ModFunc qw{ quadvec gcd };

my $zone = Math::GMP->new(1);

sub func_name { 'tau' }
sub func { tau($_[1]) }
sub func_target { $_[0]->{target} }

sub apply_m {
    my($self, $m, $fm) = @_;
    my $c = $self->c;
    my $tn = $self->func_target;

    #
    # If target is prime, the only thing we can usefully do is avoid primes
    #
    if ($tn == 2) {
        $self->series(0, $m, $m) if @$fm == 1 && $fm->[0][1] == 1;
        return;
    }

    #
    # If we know target must not divide m (maybe only beyond some maximum
    # target value), use that to deny various values (mod m).
    #
    my $max = $self->test_series($m, $tn);
    $self->series(0, $m, $max) if $max >= 0;

    #
    # target == m (mod m rad(m)) is suppressible when tau(m) ~| tau(n)
    #
    # Alternatively, when tau(n) / tau(m) is odd, we can fix any other
    # factors to be quadratic residues.
    #
    $self->test_m_rad_m($m, $tn);
    return;
}

#
# A target xp^y (mod p^{y+1}) is suppressible when y+1 ~| tau(n).
#
# More generally:
#    tau(m) ~| tau(n)  and  target is m (mod m rad(m)) => tau(target) != tau(n)
#
# However, given m = \prod p_i^a_i, if a) tau(m) > tau(n), or better:
#    a)  \not \exists k: tau(km) = tau(n)
# or:
#    b)  \exists i: tau(m / p_i^a_i) ~| tau(n)
# we will already have suppressed a superset of what this gives us.
#
# Additionally, when tau(m) | tau(n) with tau(n) / tau(m) odd, then if
# target is um (mod m rad m), gcd(u, rad m) = 1 we can suppress
# um (mod m rad m) unless u is a quadratic residue mod each prime p
# dividing m, and um (mod q m rad m) unless u is a quadratic residue mod q.
#
# We can go further when instead of q we take q^s: there are more non-residues
# mod q^s (but it means recalculating what m and rad(m) mean). For any
# q > 2, it is only useful to take q^s : gcd(pmr, q^{s+1}) = q^s since
# for an odd prime a coprime non-residue mod q^{s+1} is always a non-residue
# mod q^s. For q=2 we need a special case to cope eg with 3(mod 4).
#
sub test_m_rad_m {
    my($self, $pmr, $tn) = @_;
    my $debug = $self->debug;
    $debug > 1 && warn "m rad m: try $pmr\n";
    my($p, $m, $r) = discern_rad($pmr) or do {
        $debug > 2 && warn "m rad m: too many spares\n";
        return;
    };
    my $mf = [ factor_exp($m) ];
    my $tm = tau_factor($mf);

    my $rp = $r * $p;
    my($qm, $qr, $qrp) = ($m, $r, $rp);
    my $quad_series = sub {
        my($mod) = @_;
        my $vec = quadvec($mod);
        for my $u (1 .. $qrp - 1) {
            next if gcd($u, $qr) > 1;
            next if vec($vec, $u % $mod, 1);
            $debug > 2 && warn "m rad m: series $u * $qm, $pmr\n";
            $self->series($u * $qm, $pmr, 0);
        }
    };

    if ($p > 1) {
        $quad_series->($p) if divides_oddly($tn, $tm);
        # all remaining tests require p==1
        return;
    }

    $quad_series->($r) if divides_oddly($tn, $tm);

    if ($mf->[0][0] == 2) {
        # p=2 uniquely introduce new non-residues at higher powers; eg
        # 3 is a non-residue (mod 4).
        my $pp2 = $mf->[0][1];
        my $tm2 = $tm / ($pp2 + 1) * $pp2;
        if (divides_oddly($tn, $tm2)) {
            my $p = $mf->[0][0];
            $qm = $m / $p;
            $qrp = $rp * $p;
            my $pp = $p ** ($pp2 + 1);
            # we need to allow 2 mod 8, but not 4 mod 8
            $qr = $r;    # not $qr = $r / $p
            $quad_series->($p ** ($pp2 + 1));
        }
    }

# I think this is entirely wrong, and correcting it would make it a useless
# copy of a previous test.
#   {
#       for my $j (0 .. $#$mf) {
#           my $pp = $mf->[$j][1] + 1;
#           if (divides_oddly($tn, $tm / $pp)) {
#               my $p = $mf->[$j][0];
#               $qm = $m / $p;
#               $qr = $r / $p;
#               $qrp = $rp * $p;
#               $quad_series->($p ** $pp);
#           }
#       }
#   }

    # case a (coarse version)
    return if $tm > $tn;

    # main requirement
    return if ($tn % $tm) == 0;

    for my $i (0 .. $#$mf) {
        my $tm_i = $tm / ($mf->[$i][1] + 1);
        # case b
        return if ($tn % $tm_i) != 0;
    }

    for my $u (1 .. $r - 1) {
        next if gcd($u, $r) > 1;
        $debug && warn "m rad m: suppress $m * $u mod $pmr\n";
        $self->series($m * $u, $pmr, 0);
    }
}

sub check_fixed {
    my($self) = @_;
    my $n = $self->n;
    my $f = $self->f;
    my $tn = $self->func_target;
    #
    # Go through each k working out what factors are fixed, to find
    # opportunities for additional constraint optimizations.
    #
    my $fixpow;
    for my $k (@{ $self->to_test }) {
        my $force = $self->find_fixed($n, $tn, $k) or next;
        if ($fixpow) {
            printf "317 Ignoring secondary fix_power(%s)\n", join ', ', @$force;
            my $pell = $self->fix_pell($n, $fixpow, $force);
            if ($pell) {
                my($a, $b, $c) = @$pell;
                my $sign = ($b >= 0) ? '+' : do { $b = -$b; '-' };
                # FIXME: can we upgrade to resolving the Pell equation here?
                printf "317 satisfying Pell %s x^2 %s %s y^2 = %s\n",
                        $a, $sign, $b, $c;
            }
        } else {
            $fixpow = $force;
        }
    }
    return $fixpow;
}

#
# Given known constraints, determine all prime powers known to divide
# this target, and take advantage of that in two ways:
# - if the known factors leave the unknown factors required to be a
#   perfect square or higher power, we can apply the fix_power optimization
# - (TODO) generate a targeted tester that divides by the known factors,
#   and then uses either a fast is_prime or specific tau check
# Returns (fix_power, tester) where I<fix_power> is an arrayref of
# [ k, x, z ] suitable to pass as arguments to the fix_power method,
# and I<tester> is a subref or undef.
#
sub find_fixed {
    my($self, $n, $tn, $k) = @_;
    my $c = $self->c;
    my($float, $spare) = $self->float_spare($n, $k);
    my $fixp = gcd($spare, $float);

    # primes in $fixp will show up at a fixed power (at the level of $float)
    # and so contribute a fixed multiple to tau(target); primes in $float
    # but not in $fixp will show up at the floating power or greater.
    my($fixed_tau, $fixed_mult) = ($zone, $zone);
    my @ffloat = grep {
        my($p, $x) = @$_;
        is_fixed($p, $x, $c, $n, $fixp) ? do {
            # take the fixed power, splice out of the list of remaining
            # floating powers
            $fixed_tau *= $x + 1;
            $fixed_mult *= $p ** $x;
            0;
        } : do {
            # floating power: keep it in the list
            1;
        };
    } factor_exp($float);

    my($tneed, $sanity_check) = $tn->bdiv($fixed_tau);
    unless ($sanity_check == 0) {
        print "502 Error: fixed $fixed_tau not available in tau $tn\n";
        exit 1;
    }
    if ($tneed == 1) {
        printf <<OUT, $n, $k, $c->elapsed;
404 Error: %s + %sd must be divisible by n (%.2fs)
OUT
        exit 1;
    }
    # now handle the rest of $float
    if ($tneed & 1) {
        # we can fix a square, and maybe more - highest fixable power is
        # gcd of { p_i - 1 } over primes p_i dividing $tneed
        my $z = $tneed - 1;
        $z = gcd($z, $_->[0] - 1) for factor_exp($tneed);
        return [ $k, $fixed_mult, $z ];
    }
    return undef;
}

sub is_fixed {
    my($p, $x, $c, $n, $fixp) = @_;
    $p = Math::GMP->new($p) unless ref $p;
    if (($fixp % $p) == 0) {
        return 1;
    }
    # if p^x would float, it may still be fixed if p^{x+1} divides n
    # and 0 (mod p) is disallowed.
    if (($n % ($p ** ($x + 1))) == 0 && $c->disallowed($p, 0)) {
        return 1;
    }
    return 0;
}

#
# disallow any target being v (mod m) for target > max
#
sub series {
    my($self, $v, $m, $max) = @_;
    $max ||= 0;
    ($self->debug > 1) && warn "series ($v, $m, $max)\n";
    my $n = $self->n;
    for my $k (@{ $self->to_test }) {
        $self->suppress_k($k, $v, $m, $max);
    }
    ($self->debug > 1) && warn "series ($v, $m, $max): applied\n";
}

sub test_series {
    my($self, $m, $tn) = @_;
    my $fm = [ factor_exp($m) ];
    my $tm = tau_factor($fm);
    my $max = -1;

    #
    # If tau(m) >= tau(n) we can limit the options (mod m).
    #
    if ($tm > $tn) {
        return 0;
    } elsif ($tm == $tn) {
        $max = $m;
    }

    #
    # Let m=\prod{p_i^a_i}, m_2=\prod{p_i^b_i} where b_i is least k >= a_i
    # such that b_i + 1 divides tau(n). Then if tau(m_2) >= tau(n) we can
    # limit the options (mod m) for targets > m_2.
    #
    my $fm2 = [ map [ $_->[0], first_divisor($tn, $_->[1] + 1) - 1 ], @$fm ];
    my $m2 = unfactor_exp($fm2);
    my $tm2 = tau_factor($fm2);
# CHECKME: why do we require m2 > m?
    if ($m2 > $m) {
        if ($tm2 > $tn) {
            return 0;
        } elsif ($tm2 == $tn) {
            $max = $m2 unless $max >= 0 && $max < $m2;
        }
    }

    #
    # Let m_3=\prod{p_i^c_i}, c^i>=b_i. Valid solutions require selecting
    # k_j, m_{3j} such that tau(k_j m_{3j})=tau(n). If k_j=1 for all j,
    # we can limit the options (mod m) for targets > max(m_{3j}).
    #
    if (gcd($tm2, $tn) < $tm2) {
        my $maxm3 = 0;
        my $good = 1;
        my $iter = iter_factor($tn);
        while (my $tn_factors = $iter->()) {
            next if @$tn_factors < @$fm;
            my $fm3 = [ map [ $fm->[$_][0], $tn_factors->[$_] - 1 ], 0 .. $#$fm ];
            next if grep $fm3->[$_][1] < $fm2->[$_][1], 0 .. $#$fm;
            $good = 0, last if @$tn_factors > @$fm;
            my $m3 = unfactor_exp($fm3);
            $maxm3 = $m3 if $maxm3 < $m3;
        }
        $max = $maxm3 if $good && !($max >= 0 && $max < $maxm3);
    }
    return $max;
}

#
# Given integer s, returns (p, m, rad(m)) such that s = p m rad(m), and
# p is prime or 1. If no such split is possible (because more than one
# prime divides s with a multiplicity of one), returns nothing.
#
sub discern_rad {
    my($s) = @_;
    my $p = $zone;
    my $r = $zone;
    for (factor_exp($s)) {
        my($sp, $spp) = @$_;
        if ($spp == 1) {
            if ($p > 1) {
                return;
            } else {
                $p *= $sp;
            }
        } else {
            $r *= $sp;
        }
    }
    return +($p, $s / $p / $r, $r);
}

# 
# Given n, d, returns true if d divides n to give an odd number.
# 
sub divides_oddly {
    my($n, $d) = @_;
    my($q, $r) = $n->bdiv($d);
    return 0 if $r;
    return +($q & 1) ? 1 : 0;
}   

#
# Given factorisation of n = [ [p_1, pow_1], ... ], return rad(n) = prod(p_i)
#
sub rad_factor {
    my($f) = @_;
    my $rad = $zone;
    $rad *= $_->[0] for @$f;
    return $rad;
}

#
# Given factorisation of n = [ [p_1, pow_1], ... ], return tau(n) = prod(pow_i + 1)
#
sub tau_factor {
    my($f) = @_;
    my $tau = $zone;
    $tau *= $_->[1] + 1 for @$f;
    return $tau;
}

#
# Given n, d return the least divisor of n >= d.
# Assumed d <= n.
#
sub first_divisor {
    my($n, $d) = @_;
    ++$d while gcd($d, $n) < $d;
    return $d;
}

#
# Given [ [p_1, pow_1], ... ] return the n of which it is the factorisation.
#
sub unfactor_exp {
    my($f) = @_;
    my $n = 1;
    for (@$f) {
        $n *= $_->[0] ** $_->[1];
    }
    return $n;
}

# 
# Given n, return an iterator that yields each factorisation of n into
# parts > 1, or undef.
#
sub iter_factor {
    my $n = shift;
    my $f = [ $n ];
    return sub {
        return undef unless @$f;
        my $r = [ @$f ];
        my $spare = $zone;
        ITERFAC: while (@$f) {
            my $last = pop @$f;
            my $prod = $spare * $last;
            for my $d ($spare + 1 .. int($prod / 2)) {
                my($q, $rem) = $prod->bdiv($d);
                next if $rem;
                push @$f, $q, $d;
                last ITERFAC;
            }
            $spare = $prod;
        }
        return @$r ? $r : undef;
    };
}

sub tau {
    my $n = shift;
    my $k = $zone;
    $k *= $_->[1] + 1 for factor_exp($n);
    return $k;
}

1;
