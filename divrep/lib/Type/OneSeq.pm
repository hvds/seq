package Type::OneSeq;
use strict;
use warnings;

use parent qw{ Type::Tauish };
use Math::GMP;
use Math::Prime::Util qw{ is_prime next_prime factor_exp };
use Memoize;

use ModFunc qw/ gcd /;
*tau = \&Type::Tauish::tau;

=head1 Type::OneSeq

What is the longest possible arithmetic progression starting d and with
difference 1, such that each element e=d+k has tau(e) = n? What is
the minimal starting point of such an AP for each possible length?

=cut

my $INF = 1_000_000;

sub init {
    my($self) = @_;
    my $n = $self->n;
    $self->{target} = $n if defined $n;
    return;
}

sub name { 'oneseq' }
sub dbname { 'oneseq' }

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
    for (my $pow = 2; 1; ++$pow) {
        next if 0 == ($n % ($pow + 1));
        # if pow < n, we can have basep^pow appear once (actually as a
        # higher power), else we can't have it at all
        return +($pow < $n) ? 2 * ($basep ** $pow) : $basep ** $pow;
    }
}

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $d + $k;
}

sub to_testf {
    my($self, $f) = @_;
    return [ 0 .. $f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $n = $self->n;
    my $tau = $self->func_target;
    if ($tau == 2) {
        return [ "$k", sub { is_prime($_[0] + $k) } ];
    }
    return [ "$k", sub { $tau == tau($_[0] + $k) } ];
}

sub float_spare {
    my($self, $n, $k) = @_;
    my $c = $self->c;
    my $kmult = $c->mult;
    my $kmod = ($c->mod_mult + $k) % $kmult;
    my $float = gcd($kmult, $kmod);
    my $spare = $kmult / $float;
    return ($float, $spare);
}

sub fix_pell {
    my($self, $n, $fix1, $fix2) = @_;
    my($k1, $x1, $z1) = @$fix1;
    my($k2, $x2, $z2) = @$fix2;
    return undef unless ($z1 & 1) == 0 && ($z2 & 1) == 0;

    my $a = $x1 ** ($z1 >> 1);
    my $b = -$x2 ** ($z2 >> 1);
    my $c = $k1 - $k2;
    my $g = gcd($a, gcd($b, $c));
    return [ map $_ / $g, ($a, $b, $c) ];
}

#
# suppress the possibility that the k'th target is v (mod m) for target > max
#
sub suppress_k {
    my($self, $k, $v, $m, $max) = @_;
    $self->c->suppress($m, ($v - $k) % $m, $max - $k);
    return;
}

#
# Calculate floor(y) given d: floor(y) = floor(((d + k) / x) ^ (1/z))
#
sub dtoy {
    my($self, $c, $val) = @_;
    my $base = $val + $c->pow_k;
    return +($base / $c->pow_x)->broot($c->pow_z);
}

sub dtoceily {
    my($self, $c, $val) = @_;
    my $g = $c->pow_g;
    return $g + $self->dtoy($c, $val - $g);
}

#
# Calculate d given y: d = xy^z - k
#
sub ytod {
    my($self, $c, $val) = @_;
    return $c->pow_x * $val ** $c->pow_z - $c->pow_k;
}

#
# Given y == y_m (mod m) and d = xy^z - k, return (d_s, s) as the
# value and modulus of the corresponding constraint on d, d == d_s (mod s).
# If no valid d is possible, returns s == 0.
#
sub mod_ytod {
    my($self, $c, $val, $mod) = @_;
    my($n, $k, $x, $z) = ($c->n, $c->pow_k, $c->pow_x, $c->pow_z);
    my $base = $x * $val ** $z - $k;
# CHECKME: should we increase $mod if gcd($mod, $z) > 1?
    return ($base % $mod, $mod);
}

1;
