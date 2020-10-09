package Group;
use strict;
use warnings;

use List::Util qw{ all min max };
use Algorithm::Loops qw{ NestedLoops };

our $SYM_xy; BEGIN { $SYM_xy = 0 }; # null symmetry
our $SYM_xY; BEGIN { $SYM_xY = 1 }; # x -> x, y -> -y
our $SYM_Xy; BEGIN { $SYM_Xy = 2 };
our $SYM_XY; BEGIN { $SYM_XY = 4 };
our $SYM_yx; BEGIN { $SYM_yx = 8 };
our $SYM_yX; BEGIN { $SYM_yX = 16 };
our $SYM_Yx; BEGIN { $SYM_Yx = 32 };
our $SYM_YX; BEGIN { $SYM_YX = 64 };

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

    my $sym = $SYM_xy;
    $sym |= $SYM_Xy if all { $row[$_] eq $row[2 - $_] } 0 .. 2;
    $sym |= $SYM_xY if all { $col[$_] eq $col[2 - $_] } 0 .. 2;
    $sym |= $SYM_XY if all { $row[$_] eq reverse $row[2 - $_] } 0 .. 2;
    $sym |= $SYM_yx if all { $row[$_] eq $col[$_] } 0 .. 2;
    $sym |= $SYM_yX if all { $row[$_] eq $col[2 - $_] } 0 .. 2;
    $sym |= $SYM_Yx if all { $row[$_] eq reverse($col[$_]) } 0 .. 2;
    $sym |= $SYM_YX if all { $row[$_] eq reverse($col[2 - $_]) } 0 .. 2;

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
            for my $j (-1 .. $y) {
                next if $i >= 0 && $i < $x && $j >= 0 && $j < $y
                        && $vals->[$i][$j];
                if ($sym) {
                    next if ($sym & $SYM_xY)
                        && $seen{"$i.@{[ $y - $j ]}"};
                    next if ($sym & $SYM_Xy)
                        && $seen{"@{[ $x - $i ]}.$j"};
                    next if ($sym & $SYM_XY)
                        && $seen{"@{[ $x - $i ]}.@{[ $y - $j ]}"};
                    next if ($sym & $SYM_yx)
                        && $seen{"$j.$i"};
                    next if ($sym & $SYM_yX)
                        && $seen{"$j.@{[ $x - $i ]}"};
                    next if ($sym & $SYM_Yx)
                        && $seen{"@{[ $y - $j ]}.$i"};
                    next if ($sym & $SYM_YX)
                        && $seen{"@{[ $y - $j ]}.@{[ $x - $i ]}"};
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
    my($x, $y, $sym) = ($self->x, $self->y, $self->sym);
    my($locx, $locy) = @$locn;
    my $next_sym = 0;
    while ($sym) {
        my $bit = $sym - ($sym & ($sym - 1));
        $sym -= $bit;
        my $check = {
            $SYM_Xy => sub { $locx * 2 == $x - 1 },
            $SYM_xY => sub { $locy * 2 == $y - 1 },
            $SYM_XY => sub { $locx * 2 == $x - 1 && $locy * 2 == $y - 1 },
            $SYM_yx => sub { $locx == $locy },
            $SYM_yX => sub { $locx * 2 == $x - 1 && $locy * 2 == $y - 1 },
            $SYM_Yx => sub { $locx * 2 == $x - 1 && $locy * 2 == $y - 1 },
            $SYM_YX => sub { $locx + $locy == $x - 1 },
        };
        $next_sym |= $bit if $check->{$bit}->();
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

1;
