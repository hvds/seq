package Type::AddSeq;
use strict;
use warnings;

use parent qw{ Type::Tauish };
use Math::GMP;
use Math::Prime::Util qw{ factor_exp is_prime next_prime divmod };
use Memoize;

use ModFunc qw{ gcd };
use RootMod qw{ allrootmod };
*tau = \&Type::Tauish::tau;

=head1 Type::AddSeq

What is the longest possible arithmetic progression starting d and with
difference n, such that each element e=d+kn has tau(e) = n? What is the
minimal starting point of such an AP for each possible length?

=cut

my $INF = 1_000_000;

sub init {
    my($self) = @_;
    $self->{target} = $self->n;
    return;
}

sub name { 'addseq' }
sub dbname { 'addseq' }

memoize('gprio', NORMALIZER => sub { "$_[1]" });
sub gprio {
    my($self, $n) = @_;
    return -$INF if $n & 1;
    my @f = factor_exp($n);
    my $highest = $f[-1][0] // 1;
    return -(
        (log($n) / log(2)) ** 2
        + (log($highest) / log(2)) ** 3
    );
}

sub ming { 0 }

sub maxg {
    my($type, $n) = @_;
    die "can't maxg(0)" unless $n;
    return 1 if $n & 1;

    my $basep = 2;
    while (0 == ($n % $basep)) {
        $basep = next_prime($basep);
    }
    for (my $pow = 2; 1; ++$pow) {
        next if 0 == ($n % ($pow + 1));
        # if pow < n, we can have basep^pow appear once (actually as a
        # higher power), else we can't have it at all
        return +($pow < $n) ? 2 * ($basep ** $pow) : $basep ** $pow;
    }
}

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $d + $k * $n;
}

sub to_testf {
    my($self, $f) = @_;
    return [ 0 .. $f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $nk = $self->n * $k;
    my $tau = $self->func_target;
    if ($tau == 2) {
        return [ "$k", sub { is_prime($_[0] + $nk) } ];
    }
    return [ "$k", sub { $tau == tau($_[0] + $nk) } ];
}

sub float_spare {
    my($self, $n, $k) = @_;
    my $c = $self->c;
    my $kmult = $c->mult;
    my $kmod = ($c->mod_mult + $n * $k) % $kmult;
    my $float = gcd($kmult, $kmod);
    my $spare = $kmult / $float;
    return ($float, $spare);
}

#
# If d + kn = um (mod m rad m) and tau(m) divides tau(n) oddly, then
# d + (k + m)n does the same, so of u and u + n, either both are
# squares or one must share a factor with m.
#
sub more_m_rad_m {
    my($self, $m, $r, $tm) = @_;
    my $f = $self->f;
    return unless $f > $m;
    my $n = $self->n;
    return unless gcd($r, $n) == 1;

    my $mr = $m * $r;
    # We know n is even, so greatest squares with difference n must have
    # u = x^2, u + n = (x + 2)^2 => n = 4x + 4
    my $lim = ((($n + 4) >> 2) ** 2) * $m;
    for my $u (1 .. $r - 1) {
        next if gcd($u, $r) > 1 || gcd($u + $n, $r) > 1;
        for my $k (0 .. $f - $m) {
            $self->suppress_k($k, $u * $m, $mr, $lim);
            # or equivalently:
            #$self->suppress_k($k + $m, (($u + $n) % $r) * $m, $mr, $lim);
        }
    }
    return;
}

sub fix_pell {
    my($self, $n, $fix1, $fix2) = @_;
    my($k1, $x1, $z1) = @$fix1;
    my($k2, $x2, $z2) = @$fix2;
    return undef unless ($z1 & 1) == 0 && ($z2 & 1) == 0;

    my $a = $x1 ** ($z1 >> 1);
    my $b = -$x2 ** ($z2 >> 1);
    my $c = $n * ($k1 - $k2);
    my $g = gcd($a, gcd($b, $c));
    return [ map $_ / $g, ($a, $b, $c) ];
}

#
# suppress the possibility that the k'th target is v (mod m) for target > max
#
sub suppress_k {
    my($self, $k, $v, $m, $max) = @_;
    my $diff = $k * $self->n;
    $self->c->suppress($m, ($v - $diff) % $m, $max - $diff);
    return;
}

#
# Calculate floor(y) given d: floor(y) = floor(((d + kn) / x) ^ (1/z))
#
sub dtoy {
    my($self, $c, $val) = @_;
    my $base = $val + $c->n * $c->pow_k;
    return +($base / $c->pow_x)->broot($c->pow_z);
}

sub dtoceily {
    my($self, $c, $val) = @_;
    my $g = $c->pow_g;
    return $g + $self->dtoy($c, $val - $g);
}

#
# Calculate d given y: d = xy^z - kn
#
sub ytod {
    my($self, $c, $val) = @_;
    return $c->pow_x * $val ** $c->pow_z - $c->n * $c->pow_k;
}

#
# Given y == y_m (mod m) and d = xy^z - kn, return (d_s, s) as the
# value and modulus of the corresponding constraint on d, d == d_s (mod s).
# If no valid d is possible, returns s == 0.
#
sub mod_ytod {
    my($self, $c, $val, $mod) = @_;
    my($n, $k, $x, $z) = ($c->n, $c->pow_k, $c->pow_x, $c->pow_z);
    my $base = $x * $val ** $z - $n * $k;
# CHECKME: should we increase $mod if gcd($mod, $z) > 1?
    return ($base % $mod, $mod);
}

#
# Given d == d_m (mod m) and d = xy^z - kn, return an arrayref of arrayrefs
# (y_s, s) as the values and moduli of possible corresponding constraints
# on y, y == y_s (mod s),
#
sub mod_dtoy {
    my($self, $c, $val, $mod, $which) = @_;
    my($n, $k, $x, $z) = ($c->n, $c->pow_k, $c->pow_x, $c->pow_z);
    my $xyz = $val + $k * $n;
    my $g = gcd($xyz, $mod, $x);
    if ($g > 1) {
        $_ /= $g for ($xyz, $mod, $x);
    }
    my $yz = divmod($xyz, $x, $mod) // return [];
    return [ map [ $_, $mod ], allrootmod($yz, $z, $mod) ];
}

1;
