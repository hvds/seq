package Type::AddSeq;
use strict;
use warnings;

use parent qw{ Type };
use Math::GMP;
use Math::Prime::Util qw{ factor_exp };
use ModFunc qw{ quadvec gcd };

=head1 Type::AddSeq

What is the longest possible arithmetic progression with difference
n, such that each element e has tau(e) = n? What is the minimal
starting point of such an AP for each possible length?

=cut

my $zone = Math::GMP->new(1);

sub init {
    my($self) = @_;
    $self->{target} = $self->n;
    return;
}

sub name { 'addseq' }

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $d + $k * $n;
}

sub func_name { 'tau' }

sub func { tau($_[1]) }

sub func_target { $_[0]->{target} }

sub apply_m {
    my($self, $m, $fm) = @_;
    my $c = $self->{c};
    my $tn = $self->func_target;

    #
    # If n is prime, the only thing we can usefully do is avoid primes
    #
    if ($tn == 2) {
        $self->series(0, $m, $m) if @$fm == 1 && $fm->[0][1] == 1;
        return;
    }

    #
    # If we know n+kd must not divide m (maybe only beyond some maximum
    # value for n+kd), use that to deny various values (mod m).
    #
    my $max = $self->test_series($m, $tn);
    $self->series(0, $m, $max) if $max >= 0;

    #
    # n+kd == m (mod m rad(m)) is suppressible when tau(m) ~| tau(n)
    #
    # Alternatively, when tau(n) / tau(m) is odd, we can fix any other
    # factors to be quadratic residues.
    #
    $self->test_m_rad_m($m, $tn);
    return;
}

sub to_test {
my($self) = @_;
return [ 0 .. $self->f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $n = $self->n;
    my $tau = $self->func_target;
    if ($tau == 2) {
        return [ "$k", sub { is_prime($_[0] + $n * $k) } ];
    }
    return [ "$k", sub { $tau == tau($_[0] + $n * $k) } ];
}

#
# disallow d+kn=v (mod m) for d+kn > max
#
sub series {
    my($self, $v, $m, $max) = @_;
    my $c = $self->c;
    $max ||= 0;
    ($c->debug > 1) && warn "series ($v, $m, $max)\n";
    my $n = $self->n;
    for my $k (0 .. $c->f - 1) {
        my $diff = $k * $n;
        $c->suppress($m, ($v - $diff) % $m, $max - $diff);
    }
    ($c->debug > 1) && warn "series ($v, $m, $max): applied\n";
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
    # limit the options (mod m) for n+kd > m_2.
    #
    my $fm2 = [ map [ $_->[0], first_divisor($tn, $_->[1] + 1) - 1 ], @$fm ];
    my $m2 = unfactor_exp($fm2);
    my $tm2 = tau_factor($fm2);
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
    # we can limit the options (mod m) for n+kd > max(m_{3j}).
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
# d+kn == xp^y (mod p^{y+1}) is suppressible when y+1 ~| tau(n).
#
# More generally:
#    tau(m) ~| tau(n)  and  d+kn == m (mod m rad(m)) => tau(d+kn) != tau(n)
#
# However, given m = \prod p_i^a_i, if a) tau(m) > tau(n), or better:
#    a)  \not \exists k: tau(km) = tau(n)
# or:
#    b)  \exists i: tau(m / p_i^a_i) ~| tau(n)
# we will already have suppressed a superset of what this gives us.
#
# Additionally, when tau(m) | tau(n) with tau(n) / tau(m) odd, then if
# d + kn == um (mod m rad m), gcd(u, rad m) = 1 we can suppress
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
    my $c = $self->c;
    my $debug = $c->debug;
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

sub tau {
    my $n = shift;
    my $k = $zone;
    $k *= $_->[1] + 1 for factor_exp($n);
    return $k;
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
# Return TRUE if the check c (factorizing to fm) should be skipped for n
# according to the cp option, ie if the largest factor of c not dividing
# n is greater than cp.
#
sub cp_skip {
    my($cp, $n, $fm) = @_;
    # assume the greatest factor
    my $i = $#$fm;
    # skip past those that divide n
    --$i while $i >= 0 && ($n % $fm->[$i][0]) == 0;
    return 1 if $i >= 0 && $fm->[$i][0] > $cp;
    return 0;
}
#
# Given n, d, returns true if d divides n to give an odd number.
#
sub divides_oddly {
    my($n, $d) = @_;
use Data::Dumper; use Carp; Carp::confess(Dumper(\@_)) unless ref($n);
    my($q, $r) = $n->bdiv($d);
    return 0 if $r;
    return +($q & 1) ? 1 : 0;
}

#
# Given factorisation of n = [ [p_1, pow_1], ... ], return rad(n) = prod(p_i)
#
sub rad_factor {
    my($f) = @_;
    my $rad = MBI(1);
    $rad *= $_->[0] for @$f;
    return $rad;
}

#
# Given factorisation of n = [ [p_1, pow_1], ... ], return tau(n) = prod(pow_i + 1)
#
sub tau_factor {
    my($f) = @_;
    my $tau = 1;
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

1;
