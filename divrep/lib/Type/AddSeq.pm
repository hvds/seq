package Type::AddSeq;
use strict;
use warnings;

use parent qw{ Type::Tauish };
use Math::GMP;
use Math::Prime::Util qw{ factor_exp is_prime next_prime };
use Memoize;

use ModFunc qw{ gcd };
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
        + (log($highest) / log(2)) ** 3.5
    );
}

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

1;
