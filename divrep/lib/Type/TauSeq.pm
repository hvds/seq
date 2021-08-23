package Type::TauSeq;
use strict;
use warnings;

use parent qw{ Type::Tauish };
use Math::GMP;
use Math::Prime::Util qw{ is_prime next_prime factor_exp };
use Memoize;

use ModFunc qw/ gcd /;
*tau = \&Type::Tauish::tau;

=head1 Type::TauSeq

What is the longest possible arithmetic progression starting n and with
difference d, such that each element e=n+kd has tau(e) = tau(n)? What is
the minimal starting point of such an AP for each possible length?

=cut

sub init {
    my($self) = @_;
    my $n = $self->n;
    $self->{target} = tau($n) if defined $n;
    return;
}

sub name { 'tauseq' }
sub dbname { 'tauseq' }

memoize('gprio', NORMALIZER => sub { "$_[1]" });
sub gprio {
    my($self, $n) = @_;
    my @f = factor_exp($n);
    my $highest = $f[-1][0] // 1;
    my $prime = (@f == 1) && ($f[0][1] == 1);
    return - (
        log($n) / log(2)
        + ($prime ? 10 : 0)
        + 10 * log($highest) / log(2)
    ) / 2;
}

sub maxg {
    my($self, $n) = @_;
    return $n;
}

sub func_value {
    my($self, $n, $k, $d) = @_;
    return $n + $k * $d;
}

sub to_testf {
    my($self, $f) = @_;
    return [ 1 .. $f - 1 ];
}

sub test_target {
    my($self, $k) = @_;
    my $n = $self->n;
    my $tau = $self->func_target;

    if ($tau == 2) {
        return [ "$k", sub { is_prime($n + $_[0] * $k) } ];
    }

    my $g = gcd($k, $n);
    my $g3 = $g * $g * $g;
    if ($g > 1 && is_prime($g) && $tau == 4 && $n != $g3 && $self->min > $g3) {
        my $gn = $n / $g;
        my $gk = $k / $g;
        return [ "p$k", sub { is_prime($gn + $_[0] * $gk) } ];
    }

    return [ "$k", sub { $tau == tau($n + $_[0] * $k) } ];
}

sub float_spare {
    my($self, $n, $k) = @_;
    my $c = $self->c;
    my $kmult = $c->mult * $k;
    my $kmod = ($n + $c->mod_mult * $k) % $kmult;
    my $float = gcd($kmult, $kmod);
    my $spare = $kmult / $float;
    return ($float, $spare);
}

sub fix_pell {
    my($self, $n, $fix1, $fix2) = @_;
    my($k1, $x1, $z1) = @$fix1;
    my($k2, $x2, $z2) = @$fix2;
    return undef unless ($z1 & 1) == 0 && ($z2 & 1) == 0;

    my $a = ($k2 * $x1) ** ($z1 >> 1);
    my $b = -($k1 * $x2) ** ($z2 >> 1);
    my $c = $n * ($k2 - $k1);
    my $g = gcd($a, gcd($b, $c));
    return [ map $_ / $g, ($a, $b, $c) ];
}

#
# For certain individual cases we have proven limits that we cannot easily
# derive generically here. Avoid wasting time attempting to discover those
# unless the user asks to ignore exceptions (-ix).
#
sub check_exceptions {
    my($self) = @_;
    for (
        [ 243, 4 ], # see Constraint::Fact
        [ 2401, 3 ],
    ) {
        my($n, $f) = @$_;
        if ($self->n == $n && $self->f > $f) {
            printf <<OUT, $n, $f, $self->c->elapsed;
403 Error: f(%s) > %s known impossible by exception (%.2fs)
OUT
            exit 1;
        }
    }
}

#
# suppress the possibility that the k'th target is v (mod m) for target > max
#
sub suppress_k {
    my($self, $k, $v, $m, $max) = @_;
    my $n = $self->n;
    my $g = gcd($k, $m);
    my($gv, $gvr) = ($v - $n)->bdiv($g);
    return if $gvr != 0;    # No constraint: this modval can't happen
    my $m2 = $m / $g;
    my $inv = ($k / $g)->bmodinv($m2);
    my $v2 = ($gv * $inv) % $m2;
    ($self->debug > 2) && warn "k=$k, g=$g, gv=$gv, gvr=$gvr, m2=$m2, inv=$inv, v2=$v2\n";
# CHECKME: that "+ 1" looks unnecessarily cautious
    $self->c->suppress($m2, $v2, int(($max - $n) / $k) + 1);
    return;
}

1;
