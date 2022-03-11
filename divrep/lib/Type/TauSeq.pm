package Type::TauSeq;
use strict;
use warnings;

use parent qw{ Type::Tauish };
use Math::GMP;
use Math::Prime::Util qw{ is_prime next_prime factor_exp divmod };
use Memoize;

use ModFunc qw/ gcd /;
use RootMod qw{ allrootmod };
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

sub ming { 1 }

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
# If the order of p in n is k and p^{k+1} divides the diff, we must have
# g(n, f) = p^k g(n / p^k, f). If constraints force the diff to be so
# divisible, g(n) = p^k g(n / p^k) for the rest of the sequence, and there's
# no point recalculating them. (The -f flag forces the calculation.)
#
# FIXME: we should report this only if g(n / p^k) would fix p, else not
# all solutions for g(n / p^k) would generate solutions for g(n), and in
# particular the minimum solution may not.
#
sub check_mult {
    my($self) = @_;
    my($c, $n, $f) = ($self->c, $self->n, $self->f);
    my($mult, $mod_mult) = ($c->parent->mult, $c->parent->mod_mult);
    my $mf = [ factor_exp($mult) ];
    for (@$mf) {
        my($p, $k) = @$_;   # p^k is fixed by constraints
        next if $k < 2;
        next unless 0 == ($n % $p);
        my $pk = Math::GMP->new($p) ** $k;
        next if 0 == ($n % $pk);
        next unless 0 == ($mod_mult % $pk);
        do { $pk /= $p } while 0 != ($n % $pk);
        my $nd = $n / $pk;
        printf <<OUT, $f, $pk, $nd, $pk;
201 Dependent - %s:%s f(%s) # fix %s
OUT
        exit 0;
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

#
# Calculate floor(y) given d: floor(y) = floor(((n + kd) / x) ^ (1/z))
#
sub dtoy {
    my($self, $c, $val) = @_;
    my $base = $c->n + $c->pow_k * $val;
    return +($base / $c->pow_x)->broot($c->pow_z);
}

sub dtoceily {
    my($self, $c, $val) = @_;
    my $g = $c->pow_g;
    return $g + $self->dtoy($c, $val - $g);
}

#
# Calculate d given y: d = (xy^z - n) / k
# We expect k should always divide the expression exactly, or raise an error;
# however in rare cases that's ok, these are marked by setting $::LOOSE_CUR.
#
sub ytod {
    my($self, $c, $val) = @_;
    my $base = $c->pow_x * $val ** $c->pow_z - $c->n;
    my($div, $rem) = $base->bdiv($c->pow_k);
    if ($rem != 0) {
        use Carp;
        confess sprintf(
            "_ytod(k = %s, x = %s, z = %s => %s) not divisible by k",
            $c->pow_k, $c->pow_x, $c->pow_z, $val,
        ) unless $::LOOSE_CUR;
        # if loose is good enough, return the floor but mark it
        $div .= ':LOOSE';
    }
    return $div;
}

#
# Given y == y_m (mod m) and d = (xy^z - n) / k, return (d_s, s) as the
# value and modulus of the corresponding constraint on d, d == d_s (mod s).
# If no valid d is possible, returns s == 0.
#
sub mod_ytod {
    my($self, $c, $val, $mod) = @_;
    my($n, $k, $x, $z, $g)
            = ($c->n, $c->pow_k, $c->pow_x, $c->pow_z, $c->pow_g);
    my $base = $x * $val ** $z - $n;
    my $gbase = $base / $g;
# CHECKME: should we increase $mod if gcd($mod, $z) > 1?
    return ($gbase % $mod, $mod) if $k == $g;

    my $gk = $k / $g;
    my $g2 = gcd($gk, $mod);
    if ($g2 == 1) {
        my $inv = $gk->bmodinv($mod);
        return (($gbase * $inv) % $mod, $mod);
    } elsif (($gbase % $g2) == 0) {
        my $gmod = $mod / $g2;
        my $g2base = $gbase / $g2;
        my $g2k = $gk / $g2;
        my $inv = $g2k->bmodinv($gmod);
        return (($g2base * $inv) % $gmod, $gmod);
    } else {
        return (1, 0);
    }
}

#
# Given d == d_m (mod m) and d = (xy^z - n) / k, return an arrayref of
# arrayrefs (y_s, s) as the values and moduli of possible corresponding
# constraints on y, y == y_s (mod s),
#
sub mod_dtoy {
    my($self, $c, $val, $mod, $which) = @_;
    my($n, $k, $x, $z) = ($c->n, $c->pow_k, $c->pow_x, $c->pow_z);
    my $xyz = $val * $k + $n;
    my $g = gcd($xyz, $mod, $x);
    if ($g > 1) {
        $_ /= $g for ($xyz, $mod, $x);
    }
    my $yz = divmod($xyz, $x, $mod) // return [];
    return [ map [ $_, $mod ], allrootmod($yz, $z, $mod) ];
}

1;
