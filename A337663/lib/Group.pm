package Group;
use strict;
use warnings;

use List::Util qw{ all min max };
use Algorithm::Loops qw{ NestedLoops };

use lib 'lib';
use Sym;

my %cache_seed;
# $seed_bits[$k] represents the possible arrangements of $k squares
# around a central square (distinct up to symmetry) representing
# the 8 possible locations as octal bits 0757.
my @seed_bits = (
    undef, # [ 0000 ], but we never need 0
    undef, # [ 0400, 0200 ], but we never need 1
    [ 0300, 0210, 0201, 0050, 0500, 0401 ],
    [ 0700, 0610, 0601, 0602, 0604, 0640, 0504, 0502, 0412, 0250 ],
    [ 0505, 0704, 0514, 0710, 0702, 0550, 0512, 0641, 0611,
            0603, 0642, 0342, 0252 ],
    [ 0057, 0147, 0156, 0155, 0153, 0117, 0253, 0255, 0345, 0307 ],
    [ 0457, 0547, 0556, 0707, 0257, 0356 ],
    [ 0776, 0775 ],
    [ 0777 ],
);

sub new {
    my($class, $x, $y, $vals, $sym, $sums) = @_;
    return bless {
        x => $x,
        y => $y,
        vals => $vals,
        sym => $sym,
        sums => $sums,
    }, $class;
}

sub x { shift->{x} }
sub y { shift->{y} }
sub vals { shift->{vals} }
sub sym { shift->{sym} }

sub str {
    my($self) = @_;
    return join ' ', map join('.', @$_), @{ $self->vals };
}

# returns a list of the distinct group objects that can be seeded by
# placing $k 1s around a single $k
sub new_seed {
    my($class, $k) = @_;
    return @{ $cache_seed{$k} //= [ map $class->new_bits($k, $_), @{
        $seed_bits[$k] // die "unexpected seed request for k = $k"
    } ] };
}

sub new_bits {
    my($class, $k, $bits) = @_;
    my @row = unpack '(A3)3', sprintf '%09b', $bits;
    my @vals = map [ map $_ + 0, split // ], @row;
    my @tvals = map { my $i = $_; [ map $vals[$_][$i], 0 .. 2 ] } 0 .. 2;
    my @col = map join('', @$_), @tvals;
    $vals[1][1] = $k;

    my $sym = 0;
    for (Sym->all) {
        $sym |= $_->bit if $_->check(\@vals);
    }

    my($x, $y) = (3, 3);
    if ($row[2] eq '000') {
        --$x;
        pop @vals;
    }
    if ($row[0] eq '000') {
        --$x;
        shift @vals;
    }
    if ($col[2] eq '000') {
        --$y;
        pop @$_ for @vals;
    }
    if ($col[0] eq '000') {
        --$y;
        shift @$_ for @vals;
    }
    return $class->new($x, $y, \@vals, $sym);
}

# returns an array indexed by sums of arrayrefs of locations [x, y]
sub sums {
    my($self) = @_;
    return $self->{sums} //= do {
        my @sums;
        my($x, $y, $sym) = ($self->x, $self->y, $self->sym);
        my $vals = $self->vals;
        my %seen;
        for my $i (-1 .. $x) {
            CELL: for my $j (-1 .. $y) {
                next if $i >= 0 && $i < $x && $j >= 0 && $j < $y
                        && $vals->[$i][$j];
                for (Sym->all_bits($sym)) {
                    my($ti, $tj) = @{ $_->transform_loc($self, [ $i, $j ]) };
                    next CELL if $seen{"$ti.$tj"};
                }
                my $sum = 0;
                for my $di (-1 .. 1) {
                    next if $i + $di < 0;
                    next if $i + $di >= $x;
                    for my $dj (-1 .. 1) {
                        next if $j + $dj < 0;
                        next if $j + $dj >= $y;
                        next if $di == 0 && $dj == 0;

                        $sum += $vals->[$i + $di][$j + $dj];
                    }
                }
                if ($sum) {
                    push @{ $sums[$sum] }, [ $i, $j ];
                    $seen{"$i.$j"} = 1 if $sym;
                }
            }
        }
        \@sums;
    };
}

# returns the symmetries retained by this group when a new unique
# value is placed at $locn
sub next_sym {
    my($self, $locn) = @_;
    my($locx, $locy) = @$locn;
    my $sym = $self->sym;
    my $next_sym = 0;
    for (Sym->all_bits($sym)) {
        my($tx, $ty) = @{ $_->transform_loc($self, $locn) };
        $next_sym |= $_->bit if "$tx.$ty" eq "$locx.$locy";
    }
    return $next_sym;
}

# return a new group object made by placing the given $k at the
# specified $locn in $self
sub place {
    my($self, $locn, $k) = @_;
    my($x, $y, $vals) = ($self->x, $self->y, $self->vals);
    my $sym = $self->next_sym($locn);
    $vals = [ map [ @$_ ], @$vals ];
    my($locx, $locy) = @$locn;
    if ($locx < 0) {
        unshift @$vals, [ (0) x $y ];
        ++$locx;
        ++$x;
    }
    if ($locx >= $x) {
        push @$vals, [ (0) x $y ];
        ++$x;
    }
    if ($locy < 0) {
        unshift @$_, 0 for @$vals;
        ++$locy;
        ++$y;
    }
    if ($locy >= $y) {
        push @$_, 0 for @$vals;
        ++$y;
    }
    $vals->[$locx][$locy] = $k;
    return Group->new($x, $y, $vals, $sym, undef);
}

# return a list of new group objects, each made by placing the given $k
# at the specified $locn in $self, as well as $use 1s arranged adjacent
# to $locn (but not adjacent to any other location in which some k > 1
# has been placed)
sub place_with {
    my($self, $locn, $k, $use) = @_;
    my($x, $y, $vals) = ($self->x, $self->y, $self->vals);
    my $sym = 0;    # FIXME: work out the symmetries for each arrangement
    $vals = [ map [ @$_ ], @$vals ];
    my($locx, $locy) = @$locn;
    if ($locx < 0) {
        unshift @$vals, [ (0) x $y ];
        ++$locx;
        ++$x;
    }
    if ($locx >= $x) {
        push @$vals, [ (0) x $y ];
        ++$x;
    }
    if ($locy < 0) {
        unshift @$_, 0 for @$vals;
        ++$locy;
        ++$y;
    }
    if ($locy >= $y) {
        push @$_, 0 for @$vals;
        ++$y;
    }
    $vals->[$locx][$locy] = $k;

    # start with the 8 squares around $locn
    my @avail = ([ 1, 1, 1 ], [ 1, 0, 1 ], [ 1, 1, 1 ]);
    for my $dx (-2 .. 2) {
        my $px = $locx + $dx;
        next if $px < 0;
        next if $px >= $x;
        for my $dy (-2 .. 2) {
            my $py = $locy + $dy;
            next if $py < 0;
            next if $py >= $y;
            next if $dx == 0 && $dy == 0;
            my $val = $vals->[$px][$py];
            if ($val == 1) {
                # just reject this point
                $avail[$dx + 1][$dy + 1] = 0
                        if $dx >= -1 && $dy >= -1
                        && $dx <= 1 && $dy <= 1;
            } elsif ($val > 1) {
                # reject any point neighbouring this
                for my $d2x ($dx - 1 .. $dx + 1) {
                    next if $d2x < -1 || $d2x > 1;
                    for my $d2y ($dy - 1 .. $dy + 1) {
                        next if $d2y < -1 || $d2y > 1;
                        $avail[$d2x + 1][$d2y + 1] = 0;
                    }
                }
            }
        }
    }
    my @actual;
    for my $px (0 .. 2) {
        for my $py (0 .. 2) {
            push @actual, [ $locx + $px - 1, $locy + $py - 1 ]
                    if $avail[$px][$py];
        }
    }

    return () if @actual < $use;

    my @result;
    NestedLoops([
        [ 0 .. $#actual ],
        (sub { [ $_ + 1 .. $#actual ] }) x ($use - 1),
    ], sub {
        my @locns = @actual[@_];
        my $next_vals = [ map [ @$_ ], @$vals ];
        my($next_x, $next_y) = ($x, $y);
        my($dx, $dy) = (0, 0);
        my $minx = min(map $_->[0], @locns);
        my $miny = min(map $_->[1], @locns);
        my $maxx = max(map $_->[0], @locns);
        my $maxy = max(map $_->[1], @locns);
        if ($maxx >= $next_x) {
            push @$next_vals, [ (0) x $next_y ];
            ++$next_x;
        }
        if ($maxy >= $next_y) {
            push @$_, 0 for @$next_vals;
            ++$next_y;
        }
        if ($minx < 0) {
            unshift @$next_vals, [ (0) x $next_y ];
            ++$next_x;
            ++$dx;
        }
        if ($miny < 0) {
            unshift @$_, 0 for @$next_vals;
            ++$next_y;
            ++$dy;
        }
        $next_vals->[ $_->[0] + $dx ][ $_->[1] + $dy ] = 1 for @locns;
        push @result, Group->new($next_x, $next_y, $next_vals, $sym, undef);
    });
    return @result;
}

sub _neighbours {
    my($self, $loc) = @_;
    my @n = map [ 0, 0, 0 ], 0 .. 2;
    my($x, $y, $vals) = ($self->x, $self->y, $self->vals);
    my($lx, $ly) = @$loc;
    for my $i (-2 .. 2) {
        next if $lx + $i < 0 || $lx + $i >= $x;
        for my $j (-2 .. 2) {
            next if $ly + $j < 0 || $ly + $j >= $y;
            my $val = $vals->[$lx + $i][$ly + $j];
            next if $val == 0;
            $n[$i + 1][$j + 1] = -1
                    if $i >= -1 && $i <= 1 && $j >= -1 && $j <= 1;
            if ($val > 1) {
                for my $dx (-1 .. 1) {
                    next if $dx + $i < -1 || $dx + $i > 1;
                    for my $dy (-1 .. 1) {
                        next if $dy + $j < -1 || $dy + $j > 1;
                        $n[$dx + $i + 1][$dy + $j + 1] ||= 1;
                    }
                }
            }
        }
    }
    return \@n;
}

# return a list of new group objects, each made by coalescing the groups
# $self and $other such that the location $sloc on $self and the location
# $oloc on $other coincide; the value $k is placed in that shared location,
# and an additional $use 1s are placed (safely) adjacent to it.
sub coalesce {
    my($self, $sloc, $other, $oloc, $k, $use) = @_;
    my($sx, $sy, $svals, $ssym)
            = ($self->x, $self->y, $self->vals, $self->sym);
    my($ox, $oy, $ovals, $osym)
            = ($other->x, $other->y, $other->vals, $other->sym);
    my($slocx, $slocy) = @$sloc;
    my($olocx, $olocy) = @$oloc;

    # first look at the 3x3 square around the common location, to see
    # if the two groups can in principle fit together around there
    # _and_ leave room for an additional $use 1s.
    # we mark -1 for spot in use, 0 for a spot available to use
    my $slocal = $self->_neighbours($sloc);
    my $olocal = $other->_neighbours($oloc);
    $_->[1][1] = 1 for ($slocal, $olocal);
    my $savail = grep $_ == 0, map @$_, @$slocal;
    my $sneed = grep $_ == -1, map @$_, @$slocal;
    my $oavail = grep $_ == 0, map @$_, @$olocal;
    my $oneed = grep $_ == -1, map @$_, @$olocal;
    return () if $savail < $oneed + $use
            || $oavail < $sneed + $use;

    # now check if any transform of $other would allow $olocal to be placed
    # over $slocal without conflict, and still leaving enough free spots
    # for us to add $use 1s.
    my @sym;
    SYM: for my $sym (Sym->all) {
        my $tlocal = $sym->transform($olocal);
        my @free;
        for (my $i = 0; $i <= 2; ++$i) {
            for (my $j = 0; $j <= 2; ++$j) {
                my $sv = $slocal->[$i][$j];
                my $tv = $tlocal->[$i][$j];
                next SYM if ($sv == -1 && $tv != 0)
                        || ($tv == -1 && $sv != 0);
                push @free, [ $i - 1, $j - 1 ] if $sv == 0 && $tv == 0;
            }
        }
        next unless @free >= $use;
        push @sym, [ $sym, \@free ];
    }

    return () unless @sym;

    # now try applying the same transform on the whole of $other, to
    # coalesce with $self. If that succeeds, we then have a solution
    # for each way we can distribute $use over $free.
    my @result;
    my($sxmin, $sxmax, $symin, $symax)
            = (-$slocx, $sx - $slocx, -$slocy, $sy - $slocy);
    SYM2: for (@sym) {
        my($sym, $free) = @$_;
        my $tvals = $sym->transform($ovals);
        my($tlocx, $tlocy) = @{ $sym->transform_loc($other, $oloc) };
        my($tx, $ty) = $sym->is_transpose ? ($oy, $ox) : ($ox, $oy);

        my($txmin, $txmax, $tymin, $tymax)
                = (-$tlocx, $tx - $tlocx, -$tlocy, $ty - $tlocy);
        my $xmin = min($sxmin, $txmin);
        my $xmax = max($sxmax, $txmax);
        my $ymin = min($symin, $tymin);
        my $ymax = max($symax, $tymax);
        my $dsx = $sxmin - $xmin;
        my $dsy = $symin - $ymin;
        my $dtx = $txmin - $xmin;
        my $dty = $tymin - $ymin;
        my $cx = $xmax - $xmin;
        my $cy = $ymax - $ymin;
        my $clocx = $slocx + $dsx;
        my $clocy = $slocy + $dsy;
        my $cvals = [ map [ (0) x $cy ], 1 .. $cx ];

        for my $i (0 .. $cx - 1) {
            for my $j (0 .. $cy - 1) {
                my $bs = ($i >= $dsx && $i <= $dsx + $sx - 1
                        && $j >= $dsy && $j <= $dsy + $sy - 1);
                my $bt = ($i >= $dtx && $i <= $dtx + $tx - 1
                        && $j >= $dty && $j <= $dty + $ty - 1);
                my $vs = $bs ? $svals->[$i - $dsx][$j - $dsy] : undef;
                my $vt = $bt ? $tvals->[$i - $dtx][$j - $dty] : undef;
                if ($bs && $bt) {
                    next SYM2 if $vs && $vt;
                    $cvals->[$i][$j] = $vs || $vt;
                } elsif (!$bs && !$bt) {
                    $cvals->[$i][$j] = 0;
                } else {
                    $cvals->[$i][$j] = $vs // $vt;
                }
            }
        }

        my $finish = sub {
            my($vals, $frees) = @_;
            my($rx, $ry) = ($cx, $cy);
            my($rlocx, $rlocy) = ($clocx, $clocy);
            for ([ $k, [ 0, 0 ] ], map [ 1, $_ ], @$frees) {
                my($val, $rel) = @$_;
                my($vx, $vy) = ($rlocx + $rel->[0], $rlocy + $rel->[1]);
                if ($vx < 0) {
                    unshift @$vals, [ (0) x $ry ];
                    ++$rx;
                    ++$vx;
                }
                if ($vx >= $rx) {
                    push @$vals, [ (0) x $ry ];
                    ++$rx;
                }
                if ($vy < 0) {
                    unshift @$_, 0 for @$vals;
                    ++$ry;
                    ++$vy;
                }
                if ($vy >= $ry) {
                    push @$_, 0 for @$vals;
                    ++$ry;
                }
                $vals->[$vx][$vy] = $val;
            }
            # FIXME: is it worth finding symmetries here?
            return Group->new($rx, $ry, $vals, 0);
        };

        my @select;
        if ($use) {
            NestedLoops([
                [ 0 .. $#$free ],
                (sub { [ $_ + 1 .. $#$free ] }) x ($use - 1),
            ], sub {
                push @select, [ @$free[@_] ];
            });
        } else {
            push @select, [];
        }

        if (@select == 1) {
            push @result, $finish->($cvals, $select[0]);
        } else {
            for my $select (@select) {
                my $rvals = [ map [ @$_ ], @$cvals ];
                push @result, $finish->($rvals, $select);
            }
        }
    }
    return @result;
}

1;
