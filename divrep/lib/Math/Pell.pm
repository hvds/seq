package Math::Pell;
use strict;
use warnings;
use Math::GMP;
use Math::Prime::Util qw{
    gcd factor_exp factor divisors is_square sqrtmod is_prime chinese
};
use List::Util qw{ any };
use Math::Prime::Util::GMP;

use lib 'lib';
use ModFunc qw{ allsqrtmod };

use Exporter qw{ import };
our @EXPORT_OK = qw{
    with_signed without_signed with_limit without_limit
    full gen_pell neg_pell pell
    bxy factorx sqdiff sqsum
    cf convergents
    _pell_fund_sol _sqfree
};

sub _z { map Math::GMP->new($_), @_ }
my $zero = Math::GMP->new(0);
my $zone = Math::GMP->new(1);

{
    use Data::Dumper;
    sub _str {
        # stringify bigints in a data structure for dumpability
        !defined($_[0]) ? undef
        : !ref($_[0]) ? $_[0]
        : ref($_[0]) eq 'ARRAY' ? [ map _str($_), @{ $_[0] } ]
        : ref($_[0]) eq 'HASH' ? {
            map +($_ => _str($_[0]->{$_})), keys %{ $_[0] }
        }
        : "z_$_[0]"
    }
    sub _dump {
        local $Data::Dumper::Indent = $Data::Dumper::Sortkeys = 1;
        print Dumper(_str($_[0]));
    }
}

our $want_signed = 0;
sub with_signed {
    my($cb) = @_;
    local $want_signed = 1;
    return $cb->();
}
sub without_signed {
    my($cb) = @_;
    local $want_signed = 0;
    return $cb->();
}

our $want_limit;
sub with_limit {
    my($cb, $limit) = @_;
    local $want_limit = $limit;
    return $cb->();
}
sub without_limit {
    my($cb) = @_;
    local $want_limit = undef;
    return $cb->();
}

sub _fail {
    return sub { return +(undef, undef) };
}

sub full {
    # ax^2 + bxy +cy^2 + dx + ey + f = 0
    my($a, $b, $c, $d, $e, $f) = _z(@_);

    my $gen = [
        undef,  # abcde
        undef,  # abcd
        undef,  # abce
        undef,  # abc
        undef,  # abde
        sub { factorx($a, $b, $d, $f) },                        # abd
        undef,  # abe
        sub { factorx($a, $b, $zero, $f) },                     # ab
        undef,  # acde
        undef,  # acd
        undef,  # ace
        sub {   # ac
            return gen_pell(-$c, -$f) if $a == 1;
            my($mul) = _sqfree($a);
            my $root = ($a * $mul)->bsqrt;
            return _linear_filter(gen_pell(-$c * $mul, -$f * $mul),
                    $root, $zero, $zone, $zero);
        },
        undef,  # ade
        sub { die "Not two variable" },                         # ad
        undef,  # ae
        sub { die "Not two variable" },                         # a
        sub { _swap_filter($c, $b, $zero, $e, $d, $f) },        # bcde
        sub { _swap_filter($c, $b, $zero, $zero, $d, $f) },     # bcd
        sub { _swap_filter($c, $b, $zero, $e, $zero, $f) },     # bce
        sub { _swap_filter($c, $b, $zero, $zero, $zero, $f) },  # bc
        sub { bxy($b, $d, $e, $f) },                            # bde
        sub { bxy($b, $d, $zero, $f) },                         # bd
        sub { bxy($b, $zero, $e, $f) },                         # be
        sub { bxy($b, $zero, $zero, $f) },                      # b
        sub { _swap_filter($c, $zero, $zero, $e, $d, $f) },     # cde
        sub { _swap_filter($c, $zero, $zero, $zero, $d, $f) },  # cd
        sub { die "Not two variable" },                         # ce
        sub { die "Not two variable" },                         # c
        sub { die "Not second-order" },                         # de
        sub { die "Not second-order" },                         # d
        sub { die "Not second-order" },                         # e
        sub { die "Not second-order" },                         # none
    ]->[
        ($a == 0) * 16 + ($b == 0) * 8 + ($c == 0) * 4
                + ($d == 0) * 2 + ($e == 0)
    ];
    return $gen->() if $gen;

    # Can't handle directly, must transform to something we can.

    my $D = $b * $b - 4 * $a * $c;
    return _disc0($a, $b, $c, $d, $e, $f)
            if $D == 0;

    my $E = $b * $d - 2 * $a * $e;
    my $N = $E * $E - $D * ($d * $d - 4 * $a * $f);
    # so X^2 - DY^2 = N with X = Dy + E, Y = 2ax + by + d
    my $limit = $want_limit;
    my $iter = want_limit(sub { gen_pell($D, $N) }, undef);
    return sub {
        while (1) {
            my($X, $Y) = $iter->() // return +(undef, undef);
            my $Dx = $X - $E;
# FIXME: avoid infinite loop by tracking remainders (mod $D) and (mod $b),
# and recording when (if ever) they simultaneously hit 0 before looping.
            my $x = $Dx / $D;
            return +(undef, undef) if defined($limit) && abs($x) > $limit;
            next unless 0 == ($Dx % $D);
            my $yb = $Y - $d - 2 * $a * $x;
            next unless 0 == ($yb % $b);
            return +($x, $yb / $b);
        }
    };
}

sub _exhaust {
    my($list) = @_;
    my $i = 0;
    my $limit = $want_limit;
    return sub {
        my($x, $y) = @{ $list->[$i++] // return +(undef, undef) };
        return +(undef, undef) if defined($limit) && abs($x) > $limit;
        return +($x, $y);
    };
}

sub _swap_filter {
    my($r) = @_;
    my $limit = $want_limit;
    return sub {
        my($x, $y) = $r->();
        return +(undef, undef) unless defined($x);
        return +(undef, undef) if defined($limit) && abs($y) > $limit;
        return +($y, $x);
    };
}

sub _mul_filter {
    my($r, $f) = @_;
    return $r if $f == 1;
    my $limit = $want_limit;
    return sub {
        my($x, $y) = $r->();
        return +(undef, undef) unless defined $x;
        return +(undef, undef) if defined($limit) && abs($x * $f) > $limit;
        return +($x * $f, $y * $f);
    };
}

sub _linear_filter {
    my($r, $a, $b, $c, $d) = @_;
    # $r returns (ax + b, cy + d) pairs
    my $limit = $want_limit;
    return sub {
      retry:
        my($axb, $cyd) = $r->();
        return +(undef, undef) unless defined $axb;
        my $ax = $axb - $b;
        return +(undef, undef) if defined($limit) && abs($ax / $a) > $limit;
        goto retry if $ax % $a;
        my $cy = $cyd - $d;
        goto retry if $cy % $c;
        return +($ax / $a, $cy / $c);
    };
}

sub _2sign_filter {
    # duplicate (x, y) to ((x, y), (-x, -y))
    my($r) = @_;
    return $r unless $want_signed;
    my($x, $y);
    my $state = 0;
    return sub {
        return +(undef, undef) if !defined $state;
        while (1) {
            ($x, $y) = $r->() unless $state;
            $state = undef, return +(undef, undef) unless defined $x;
            if ($state == 0) {
                ++$state;
                return +($x, $y);
            }
            $state = 0;
            return +(-$x, -$y) if $x || $y;
        }
    };
}

sub _4sign_filter {
    # duplicate (x, y) to ((x, y), (x, -y), (-x, y), (-x, -y))
    my($r) = @_;
    return $r unless $want_signed;
    my($x, $y);
    my $state = 0;
    return sub {
        return +(undef, undef) if !defined $state;
        while (1) {
            ($x, $y) = $r->() unless $state;
            $state = undef, return +(undef, undef) unless defined $x;
            if ($state == 0) {
                ++$state;
                return +($x, $y);
            }
            if ($state == 1) {
                ++$state;
                return +($x, -$y) if $y;
            }
            if ($state == 2) {
                ++$state;
                return +(-$x, $y) if $x;
            }
            $state = 0;
            return +(-$x, -$y) if $x && $y;
        }
    };
}

sub _factor_filter {
    my($f) = @_;
    my @d = divisors(abs($f));
    my $limit = $want_limit;
    return sub {
        my $d = shift(@d) // return +(undef, undef);
        return +(undef, undef) if defined($limit) && abs($d) > $limit;
        return +($d, $f / $d);
    };
}

sub _ordered {
    my($list) = @_;
    [ sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$list ];
}

sub _cmp_rxy {
    my($left, $right) = @_;
    return +(
        abs($left->[1]) <=> abs($right->[1])
        || ($left->[1] < 0) <=> ($right->[1] < 0)
        || abs($left->[2]) <=> abs($right->[2])
        || ($left->[2] < 0) <=> ($right->[2] < 0)
    );
}

sub _push_heap {
    my($h, $o) = @_;
    my $node = @$h;
    push @$h, $o;
    while ($node) {
        my $parent = ($node - 1) >> 1;
        last if _cmp_rxy($h->[$parent], $h->[$node]) <= 0;
        @$h[$parent, $node] = @$h[$node, $parent];
        $node = $parent;
    }
    return;
}

sub _pop_heap {
    my($h) = @_;
    return undef unless @$h;
    my $value = $h->[0];
    my $switch = pop @$h;
    $h->[0] = $switch if @$h;
    my $node = 0;
    my $size = @$h;
    while (1) {
        my $child = ($node << 1) + 1;
        last if $child >= $size;
        ++$child if (
            $child + 1 < $size
            && _cmp_rxy(@$h[$child, $child + 1]) > 0
        );
        last if _cmp_rxy(@$h[$child, $node]) >= 0;
        @$h[$node, $child] = @$h[$child, $node];
        $node = $child;
    }
    return $value;
}

sub _interleave {
    my $heap = [];
    for my $r (@_) {
        my($x, $y) = $r->();
        _push_heap($heap, [ $r, $x, $y ])
                if defined $x;
    }
    return sub {
        my $o = _pop_heap($heap);
        return +(undef, undef) unless defined $o;
        my($r, $x, $y) = @$o;
        my($nx, $ny) = $r->();
        _push_heap($heap, [ $r, $nx, $ny ]) if defined $nx;
        return +($x, $y);
    };
}

sub _ceil {
    my $i = int($_[0]);
    return +($i == $_[0]) ? $i : $i + 1;
}

sub _floor {
    my $i = int($_[0]);
    return +($i == $_[0]) ? $i - 1 : $i;
}

sub _quad_ranges {
    my($a, $b, $c) = @_;
    # Inflection point k_m of the parabola occurs at k = -b/2a; zeros k_{-}
    # and k_{+} occur at (-b +/- (b^2 - 4ac))/2a so the segmentation points
    # are: -inf, k_{-}, k_m, k_{+}, +inf
    my $mid_low = (-$b) / 2 / $a;
    my $mid_high = $mid_low + 1;
    my $disc = $b * $b - 4 * $a * $c;
    if ($disc < 0) {
        return +($want_signed || $a > 0)
            ? (
                [ undef, $mid_low, -1 ],
                [ $mid_high, undef, 1 ],
            )
            : ();
    }
    my $root_low = $disc->bsqrt;
    my $left_low = int((-$b - $root_low) / 2 / $a);
    my $right_low = int((-$b + $root_low + 1) / 2 / $a);
    return $want_signed ? (
        [ undef, $left_low, -1 ],
        [ $left_low + 1, $mid_low, 1 ],
        [ $mid_high, $right_low, -1 ],
        [ $right_low + 1, undef, 1 ],
    ) : ($disc < 0) ? (
        [ $left_low + 1, $mid_low, 1 ],
        [ $mid_high, $right_low, -1 ],
    ) : (
        [ undef, $left_low, -1 ],
        [ $right_low + 1, undef, 1 ],
    );
}

sub _intersect_ranges {
    my($x, $y) = @_;
    my($xi, $yi) = (0, 0);
    my @z;
    while ($xi < @$x && $yi < @$y) {
        my($xl, $xr, $xd) = @{ $x->[$xi] };
        my($yl, $yr, $yd) = @{ $y->[$yi] };
        if (defined($xl) && defined($yr) && $xl > $yr) {
            ++$yi;
            next;
        }
        if (defined($xr) && defined($yl) && $xr < $yl) {
            ++$xi;
            next;
        }
        my $zl = defined($xl) ? defined($yl)
                ? ($xl < $yl) ? $yl : $xl : $xl : $yl;
        my $zr = defined($xr) ? defined($yr)
                ? $xr > $yr ? $yr : $xr : $xr : $yr;
        my $zd = $xd;
        push @z, [ $zl, $zr, $zd ];
        last unless defined $zr;
        ++$xi if defined($xr) && $xr == $zr;
        ++$yi if defined($yr) && $yr == $zr;
    }
    return @z;
}

sub _range_quad {
    my($k, $end, $step, $signed, $x_a, $x_b, $x_c, $y_a, $y_b, $y_c) = @_;
    ($k, $end) = ($end, $k) if $step < 0;
    my $done = 0;
    my $limit = $want_limit;
    return sub {
      retry:
        return +(undef, undef) if $done;
        my($x, $y) = (
            ($x_a * $k + $x_b) * $k + $x_c,
            ($y_a * $k + $y_b) * $k + $y_c,
        );
        return +(undef, undef) if defined($limit) && abs($x) > $limit;
        $k += $step;
        $done = 1 if defined($end) && (
            ($step > 0 && $k > $end)
            || ($step < 0 && $k < $end)
        );
        goto retry if !$signed && ($x < 0 || $y < 0);
        return +($x, $y);
    };
}

=head2 _quad_param ( $x_a, $x_b, $x_c, $y_a, $y_b, $y_c )

Return a list of generators each of which generates an ordered set of
C<(x, y)> results, so that they can be interleaved to provide a complete
ordered set of results.

The results are parametrized as integer values of C<k> for
C<(x_a k^2 + x_b k + x_c, y_a k^2 + y_b k + y_c)>. As such we want to
find related points on two parabolas, each of which has up to four
segments that need to be treated separately, so in principle we may need
to return up to 8 generators. In practice, the most we ever need is 4.

=cut

sub _quad_param {
    my($x_a, $x_b, $x_c, $y_a, $y_b, $y_c) = @_;
    die "_quad_param requires quadratics"
            unless $x_a && $y_a;
    # If a is positive, the inner two segments will give negative results,
    # or the outer two if a is negative.
    # If we do _not_ want signed results, we pick the two positive segments
    # of the x parabola, and find intersections of those with the positive
    # segments of the y parabola; if we _do_ want signed results, we just
    # take the 4 segments of the x parabola: we don't need to segment the
    # y parabola in that case.
    my @xr = _quad_ranges($x_a, $x_b, $x_c);
    if ($want_signed) {
        return map _range_quad(@$_, 1, @_), @xr;
    }
    my @yr = _quad_ranges($y_a, $y_b, $y_c);
    return map _range_quad(@$_, 0, @_), _intersect_ranges(\@xr, \@yr);
}

sub _disc0 {
    my($a, $b, $c, $d, $e, $f) = @_;
    # (2ax + by + d)^2 = y(2bd - 4ae) + (d^2 - 4af)
    my $E = 2 * $b * $d - 4 * $a * $e;
    my $F = $d * $d - 4 * $a * $f;
    my @sol = map {
        my($T) = _z($_);
        # t = 2ax + by + d = kE + T, t^2 = Ey + F
        # -> y = ((kE + T)^2 - F) / E
        my($y_a, $y_b, $y_c) = ($E, 2 * $T, ($T * $T - $F) / $E);
        # x = (t - by - d)/2a = (-by + kE + T - d)/2a
        my($x_a, $x_b, $x_c)
                = (-$b * $y_a, $E - $b * $y_b, $T - $d - $b * $y_c);
        my $x_d = 2 * $a;
        my $gx = gcd($x_a, $x_b, $x_c, $x_d);
        die "Expected x: ($x_a k^2 + $x_b k + $x_c) / $x_d to divide evenly"
                if $gx != abs($x_d);
        [ (map $_ / $x_d, ($x_a, $x_b, $x_c)), ($y_a, $y_b, $y_c) ];
    } allsqrtmod($F, abs($E));
    return _interleave(map _quad_param(@$_), @sol);
}

sub factorx {
    my($a, $b, $d, $f) = @_;
    # x(ax + by + d) = -f
    die "Not two variable" unless $b;
    my $limit = $want_limit;
    my $r = without_limit(sub {
        with_signed(sub { _2sign_filter(_factor_filter(-$f)) });
    });
    return $want_signed
        ? sub {
          retry:
            my($x, $axbyd) = $r->();
            return +(undef, undef) unless defined $x;
            return +(undef, undef) if defined($limit) && abs($x) > $limit;
            my $by = $axbyd - $a * $x - $d;
            goto retry if $by % $b;
            return +($x, $by / $b);
        }
        : sub {
          retry:
            my($x, $axbyd) = $r->();
            return +(undef, undef) unless defined $x;
            goto retry if $x < 0;
            return +(undef, undef) if defined($limit) && abs($x) > $limit;
            my $by = $axbyd - $a * $x - $d;
            goto retry if $by % $b;
            my $y = $by / $b;
            goto retry if $y < 0;
            return +($x, $y);
        };
}

sub bxy {
    my($b, $d, $e, $f) = @_;
    # bxy + dx + ey + f = 0
    if ($d == 0) {
        die "Too many solutions" if $f == 0;
        if ($e == 0) {
            # bxy = -f
            return _fail() if $f % $b;
            return _2sign_filter(_factor_filter(-$f / $b));
        }
        # (bx + e)y = -f
        my $r = without_limit(sub { _2sign_filter(_factor_filter(-$f)) });
        return _linear_filter($r, $b, $e, $zone, $zero);
    }
    if ($e == 0) {
        die "Too many solutions" if $f == 0;
        # x(by + d) = -f
        my $r = without_limit(sub { _2sign_filter(_factor_filter(-$f)) });
        return _linear_filter($r, $zone, $zero, $b, $d);
    }
    # (bx + e)(by + d) = de - bf
    my $g = $d * $e - $b * $f;
    die "Too many solutions" if $g == 0;
    my $r = without_limit(sub { _2sign_filter(_factor_filter($g)) });
    return _linear_filter($r, $b, $e, $b, $d);
}

sub gen_pell {
    # x^2 - Dy^2 = N
    my($D, $N) = @_;
    return sqdiff(1, $D->bsqrt, $N) if is_square($D);
    return sqsum(-$D, $N) if $D < 0;
    return pell($D) if $N == 1;
    return neg_pell($D) if $N == -1;

    my $neg = ($N < 0) ? 1 : 0;
    my $aN = abs($N);

    # If D is not a quadratic residue (mod |N|), there can be no solution.
    my $q = sqrtmod($D, $aN) // return _fail();
    my @match;
    for my $P_0 ($q, ($q ? $aN - $q : ())) {
        my @best;
        my($cf, $cfr) = cf($D, $P_0, $zone, $aN);
        my $conv = convergents($cf, $cfr);
        my($P_i, $Q_i, $A_i, $B_i) = ($P_0, $aN, $conv->());
        for my $cfi (@$cf, @$cfr) {
            # calculate 1 / ((P + sqrt(D)) / Q - x)
            #         = (Qx - P + sqrt(D)) / ((D - (P + Qx)^2 / Q)
            my $disc = $P_i * $P_i - $D;
            die "logic error: expect $P_i^2 == $D (mod $Q_i)"
                    if $disc % $Q_i;
            my($P_n, $Q_n) = (
                $Q_i * $cfi - $P_i,
                $cfi * ( 2 * $P_i - $Q_i * $cfi) - $disc / $Q_i,
            );
            my($A_n, $B_n) = $conv->();
            if ($Q_n == 1) {
                my $G_i = abs($aN * $A_i - $P_0 * $B_i);
                if ($G_i * $G_i - $D * $B_i * $B_i == $N) {
                    @best = ($G_i, $B_i) if !@best || $G_i < $best[0];
                }
            }
            ($P_i, $Q_i, $A_i, $B_i) = ($P_n, $Q_n, $A_n, $B_n);
        }
        # It isn't entirely clear from the paper, but I _think_ we get
        # exactly 0 or 1 fundamental solutions per $P_0.
        push @match, [ @best ] if @best;
    }
    @match = sort { $a->[0] <=> $b->[0] } @match;
    return _fail() unless @match;
    my($e, $f) = _pell_fund_sol($D);
    my $i = 0;
    my $limit = $want_limit;
    return sub {
        my $which = $match[$i];
        $i = ($i + 1) % @match;
        my($x, $y) = @$which;
        return +(undef, undef) if defined($limit) && abs($x) > $limit;
        # next x' + y' sqrt(D) = (x + y sqrt(D))(e + f sqrt(D))
        @$which = (
            $x * $e + $y * $f * $D,
            $x * $f + $y * $e,
        );
        return +($x, $y);
    };
}

sub neg_pell {
    # x^2 - Dy^2 = -1
    my($D) = @_;
    my($cf, $cfr) = cf($D, $zero, $zone, $zone);
    if (@$cfr == 0) {
        die "TODO: can't handle neg_pell($D) since $D is a square";
    }
    # negative Pell's equation has no solutions if period is even
    return _fail() if (@$cfr & 1) == 0;
    my $conv = convergents($cf, $cfr);
    my($p, $q, $ok);
    for (1 .. @$cf + @$cfr) {
        ($p, $q) = $conv->();
        $ok = 1, last if $p * $p - $D * $q * $q == -1;
    }
    return _fail() unless $ok;
    my($x, $y) = ($p, $q);
    my $limit = $want_limit;
    return sub {
        my($x1, $y1) = ($x, $y);
        return +(undef, undef) if defined($limit) && abs($x1) > $limit;
        ($x, $y) = (
            $p * $p * $x + $D * $q * $q * $x + 2 * $D * $p * $q * $y,
            $p * $p * $y + $D * $q * $q * $y + 2 * $p * $q * $x,
        );
        return +($x1, $y1);
    };
}

sub pell {
    # x^2 - Dy^2 = 1
    my($D) = @_;
    my($x0, $y0) = _pell_fund_sol($D);
    if (!defined $x0) {
        return sqdiff(1, $y0, 1);
    }
    my($x, $y) = ($zone, $zero);
    my $limit = $want_limit;
    return sub {
        my($x1, $y1) = ($x, $y);
        return +(undef, undef) if defined($limit) && abs($x1) > $limit;
        ($x, $y) = ($x * $x0 + $D * $y * $y0, $y * $x0 + $x * $y0);
        return +($x1, $y1);
    }
}

=head2 sqdiff ( $a, $b, $N )

Given positive integers C<a, b> and integer C<N>, returns an iterator
that returns pairs C< (x, y) > satisfying C< (ax)^2 - (by)^2 = N >;
if not under the effect of L<"with_signed>, restricts the results to
those with C< x >= 0, y >= 0 >.

The iterator will return C<(undef, undef)> when there are no more
solutions to return.

=cut

sub sqdiff {
    my($a, $b, $N) = @_;
    my $limit = $want_limit;
    # solutions of (ax)^2 = (by)^2 + N
    if ($N < 0) {
        # make the sign filter the outermost wrapper for consistency
        return _4sign_filter(_swap_filter(
            without_limit(sub { without_signed(sub { sqdiff($b, $a, -$N) }) })
        ));
    }
    if ($N == 0) {
        my($x, $y) = ($zero, $zero);
        return _4sign_filter(sub {
            my($xn, $yn) = ($x, $y);
            return +(undef, undef) if defined($limit) && abs($xn) > $limit;
            ($x, $y) = ($x + $b, $y + $a);
            return +($xn, $yn);
        });
    }
    if ($N & 1) {
        my @d = divisors($N);
        # keep only divisors >= sqrt(N), so we generate (x, y) >= 0 in order
        splice @d, 0, @d >> 1;
        return _4sign_filter(sub {
          retry:
            my $d = shift(@d) // return +(undef, undef);
            # (ax + by)(ax - by) = (d)(N / d)
            my($ax, $by) = (($d + $N / $d) / 2, ($d - $N / $d) / 2);
            my $x = $ax / $a;
            return +(undef, undef) if defined($limit) && abs($x) > $limit;
            goto retry if $ax % $a || $by % $b;
            return +($x, $by / $b);
        });
    }
    return _fail() if ($N % 4) == 2;
    my $qN = $N / 4;
    my @d = divisors($qN);
    # keep only divisors >= sqrt(qN), so we generate (x, y) >= 0 in order
    splice @d, 0, @d >> 1;
    return _4sign_filter(sub {
      retry:
        my $d = shift(@d) // return +(undef, undef);
        # (ax + by)(ax - by) = (2d)(N / 2d)
        my($ax, $by) = ($d + $qN / $d, $d - $qN / $d);
        my $x = $ax / $a;
        return +(undef, undef) if defined($limit) && abs($x) > $limit;
        goto retry if $ax % $a || $by % $b || $by < 0;
        return +($x, $by / $b);
    });
}

=head2 sqsum ( $D, $N )

Given positive integer C<D> and integer C<N>, returns an iterator
that returns pairs C< (x, y) > satisfying C< x^2 + Dy^2 = N >;
if not under the effect of L<"with_signed>, restricts the results to
those with C<< x >= 0, y >= 0 >>.

The iterator will return C<(undef, undef)> when there are no more
solutions to return.

=cut

sub sqsum {
    my($D, $N) = @_;
    # solutions of x^2 + Dy^2 = N
    return _fail() if $N < 0;
    return _exhaust([[ 0, 0 ]]) if $N == 0;

    ($D, my $qD) = _sqfree($D);
    if ($qD > 1) {
        return _4sign_filter(_linear_filter(
            without_limit(sub { without_signed(sub { sqsum($D, $N) }) }),
            1, 0, $qD, 0,
        ));
    }

    my @d = factor_exp($N);
    if ($D == 1) {
        my $f = $zone;
        my @sol;
        for (@d) {
            my($p, $e) = @$_;
            my $sol = _sqsum_prime($p);
            if ($sol) {
                push @sol, ($sol) x $e;
            } elsif ($e & 1) {
                return _fail();
            } else {
                push @sol, ([ $p, 0 ]) x ($e >> 1);
            }
        }

        return _4sign_filter(_mul_filter(
            without_limit(sub { _exhaust(_ordered(_combine_gauss(\@sol))) }),
            $f,
        ));
    }

    # For each z: z^2 | N, find coprime solutions to x^2 + Dy^2 = N/z^2
    # to yield (zx, zy).
    # To find coprime solutions for N, we find each T: T^2 == -D (mod N)
    # and solve for x = Ty (mod N) by setting T^2 + D = qN, x = Ty - kN
    # and then solving the resulting equation qy^2 - 2tky + Nk^2 = 1
    # by setting y/k to the convergents of t/q.

    my $maxz = $zone;
    $maxz *= $_->[0] ** ($_->[1] >> 1) for @d;

    my(@sol, %seen);
    for my $z (divisors($maxz)) {
        my $N2 = $N / ($z * $z);
        for my $T (allsqrtmod(-$D, $N2)) {
            my $q = ($T * $T + $D) / $N2;
            my $r = ($q == 1)
                ? sub { return +($zone, $zero) }
                : ($T == 0)
                    ? sub { return +($zero, $zone) }
                    : convergents(cf($zone, $T, $zero, $q));
            while (1) {
                my($y, $k) = $r->();
                return _fail() unless defined $y;
                next unless ($q * $y - 2 * $T * $k) * $y + $N2 * $k * $k == 1;
                my $x = abs($T * $y - $k * $N2);
                push @sol, [ $x * $z, $y * $z ] unless $seen{"$x:$y"}++;
                last;
            }
        }
    }
    return _4sign_filter(_exhaust(_ordered(\@sol)));
}

sub _combine_gauss {
    my($sol) = @_;
    my @result = ([ $zero, $zone ], [ $zone, $zero ]);
    for (@$sol) {
        my($x, $y) = @$_;
        my %seen;
        @result = grep !$seen{"$_->[0]:$_->[1]"}++, map {
            my($p, $q) = @$_;
            +(
                [ abs($p * $x - $q * $y), $p * $y + $q * $x ],
                [ abs($p * $x + $q * $y), abs($p * $y - $q * $x) ],
            )
        } @result;
    }
    return \@result;
}

=head2 _sqsum_prime($p)

Given a prime C<$p>, returns the arrayref C<[ x, y ]> such that
C<< x^2 + Dy^2 = p, x <= y >>, or C<undef> if there is no such pair.

L<Brillhart 1972|https://www.ams.org/journals/mcom/1972-26-120/S0025-5718-1972-0314745-6/S0025-5718-1972-0314745-6.pdf>

=cut

sub _sqsum_prime {
    my($p) = @_;
    return [ $zone, $zone ] if $p == 2;
    return undef if ($p % 4) == 3;
    my $r0 = sqrtmod($p - 1, $p);
    my $r1 = $p % $r0;
    # Brillhart has unnecessary special-case for $r1 == 1 at this stage
    ($r0, $r1) = ($r1, $r0 % $r1) while $r0 * $r0 > $p;
    return [ $r1, $r0 ];
}

sub _pell_fund_sol {
    my($D) = @_;
    my($cf, $cfr) = cf($D, $zero, $zone, $zone);
    if (@$cfr == 0) {
        # D is a square
        return +(undef, $cf->[0]);
    }
    my $conv = convergents($cf, $cfr);
    # not sure how far we have to go in the worst case, but D=58 needs
    # @$cf + @$cfr * 2 - 1
    for (1 .. @$cf + @$cfr * 2) {
        my($p, $q) = $conv->();
        return +($p, $q) if $p * $p - $D * $q * $q == 1;
    }
    die "TODO: No principle solution found for pell($D)";
}

=head2 cf ( $d, $a, $b, $c )

Calculates the periodic continued fraction for C<(a + b sqrt(d)) / c>, and
returns a list of two arrayrefs in which the first arrayref holds the
(possibly empty) initial non-repeating terms of the continued fraction,
and the second arrayref holds the repeating terms.

C<$a, $b, $c, $d> are all assumed to be integers, with C<$c> non-zero.
If C<$d> is a square or C<$b> is zero, the second arrayref will be empty
(since the expression is then rational).

=cut

sub cf {
    # a/c + b/c sqrt(d)
    my($d, $a, $b, $c) = @_;
    my @cf;
    my %seen;
    while (1) {
        return +(\@cf, []) if $c == 0 || ($a == 0 && $b == 0);
        ($a, $b, $c) = map -$_, ($a, $b, $c) if $c < 0;
        my $g = gcd($a, $b, $c);
        if ($g != 1) {
            ($a, $b, $c) = map $_ / $g, ($a, $b, $c);
        }
        my $key = "$a:$b:$c";
        last if $seen{$key};
        # Taking floor(sqrt(b^2d)) gives us the accuracy we need
        my $bsd = ($b * $b * $d)->bsqrt;
        if ($b < 0) {
            $bsd = -$bsd;
            --$bsd unless is_square($d);
        }
        my $x = int(($a + $bsd) / $c);
        push @cf, $x;
        $seen{$key} = @cf;

        # now calculate 1 / (a/c + b/c sqrt(d) - x)
        # = (ac - c^2x - bc sqrt(d)) / ((a + cx)^2 - b^2d)
        my $cx = $c * $x;
        my $amcx = $a - $cx;
        ($a, $b, $c) = (
            $c * $amcx,
            - $b * $c,
            $amcx * $amcx - $b * $b * $d,
        );
    }
    return (\@cf, [ splice @cf, $seen{"$a:$b:$c"} - 1 ]);
}

=head2 convergents ( $cf, $cfr )

Given a continued fraction as returned by C<cf()>, returns an iterator
that returns successive convergents in the form C<($a, $b)>, such that
C< $a / $b > is the value of the convergent. If the continued fraction
is finite (ie it represents a rational), then after the final pair
representing the exact rational has been returned, further calls will
yield C<(undef, undef)>.

=cut

sub convergents {
    my($cf, $cfr) = @_;
    my($p0, $p1, $q0, $q1) = ($zero, $zone, $zone, $zero);
    my $i = 0;
    return sub {
        my $next = ($i < @$cf) ? $cf->[$i]
            : @$cfr ? $cfr->[($i - @$cf) % @$cfr]
            : return +(undef, undef);
        ++$i;
        ($p0, $p1) = ($p1, $p1 * $next + $p0);
        ($q0, $q1) = ($q1, $q1 * $next + $q0);
        return +($p1, $q1);
    };
}

=head2 _sqfree ( $n )

Given a non-zero integer C<n> returns the pair C<(n/k, k)>, where C<k>
is the greatest square that divides C<n>.

=cut

sub _sqfree {
    my($n) = @_;
    my $k = $zone;
    for (factor_exp(abs($n))) {
        my($p, $e) = @$_;
        $k *= $p ** ($e >> 1) if $e > 1;
    }
    return +($n / $k / $k, $k);
}

1;
