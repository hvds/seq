package Topsu;
use strict;
use warnings;

use Algorithm::Loops qw{ NestedLoops };

sub new {
    my($class, $in) = @_;
    my $self = bless {
        in => $in,
    }, $class;
    $self->init;
    return $self;
}

sub new_fh {
    my($class, $fh) = @_;
    my %seen = ('.' => undef);
    my $next = 0;
    return $class->new([
        map {
            s/^\s+//; s/\s+\z//;
            [ map {
                exists($seen{$_}) ? $seen{$_} : ($seen{$_} = $next++)
            } split //, $_ ]
        } <$fh>
    ]);
}

sub init {
    my($self) = @_;
    my $in = $self->{in};
    my $size = @$in;
    my(@x, @y, @z, @p, @r);
    for my $x (0 .. $size - 1) {
        for my $y (0 .. $size - 1) {
            my $z = $in->[$x][$y] // next;
            ++$x[$x];
            ++$y[$y];
            ++$z[$z];
            push @p, [ $x, $y, $z ];
            if ($r[$z]) {
                push @{ $r[$z]{p} }, [ $x, $y ];
                $r[$z]{xmin} = $x if $r[$z]{xmin} > $x;
                $r[$z]{xmax} = $x if $r[$z]{xmax} < $x;
                $r[$z]{ymin} = $y if $r[$z]{ymin} > $y;
                $r[$z]{ymax} = $y if $r[$z]{ymax} < $y;
            } else {
                $r[$z] = {
                    p => [ [ $x, $y ] ],
                    xmin => $x, xmax => $x,
                    ymin => $y, ymax => $y,
                };
            }
        }
    }
    my(@b, @s);
    for (@r) {
        $b[$_->{xmin}] //= $_->{xmax};
        $b[$_->{xmin}] = $_->{xmax} if $b[$_->{xmin}] < $_->{xmax};
        $s[$_->{ymin}] //= $_->{ymax};
        $s[$_->{ymin}] = $_->{ymax} if $s[$_->{ymin}] < $_->{ymax};
    }
    my(@band, @stack);
    for (my $bi = 0; $bi < $size; ++$bi) {
        my $max = $bi;
        for (my $bj = $bi; $bj <= $max; ++$bj) {
            $max = $b[$bj] if defined($b[$bj]) && $max < $b[$bj];
        }
        push @band, [ $bi, $max ];
        $bi = $max;
    }
    for (my $si = 0; $si < $size; ++$si) {
        my $max = $si;
        for (my $sj = $si; $sj <= $max; ++$sj) {
            $max = $s[$sj] if defined($s[$sj]) && $max < $s[$sj];
        }
        push @stack, [ $si, $max ];
        $si = $max;
    }
    my(@bi, @si);
    for (0 .. $#band) {
        my $b = $band[$_];
        @bi[ $b->[0] .. $b->[1] ] = ($_) x ($b->[1] + 1 - $b->[0]);
    }
    for (0 .. $#stack) {
        my $s = $stack[$_];
        @si[ $s->[0] .. $s->[1] ] = ($_) x ($s->[1] + 1 - $s->[0]);
    }
    for (@r) {
        $_->{b} = $bi[$_->{xmin}];
        $_->{s} = $si[$_->{ymin}];
    }
    @$self{qw{ size x y z p band stack bi si region }}
            = ($size, \@x, \@y, \@z, \@p, \@band, \@stack, \@bi, \@si, \@r);
    return;
}

sub _solvable {
    my($size, $valid, $unused, $used, $p) = @_;
    if (@$used == $size - 1) {
        # solution, mark all used entries as valid
        $valid->{"$_->[0]-$_->[1]-$_->[2]"} = 1 for (@$used, $p);
        return 1;
    }
    $used = [ @$used, $p ];
    $unused = [
        grep $_->[0] != $p->[0] && $_->[1] != $p->[1] && $_->[2] != $p->[2],
                @$unused
    ];
    for (@$unused) {
        return 1 if _solvable($size, $valid, $unused, $used, $_);
    }
    return 0;
}

# True if:
# - each row, column and region has at least two points;
# - each point can form part of a complete solution;
# - the regions are connected;
# - rows in each band and columns in each stack are connected.
sub is_valid {
    my($self) = @_;
    my $size = $self->{size};
    # every row, column and region must have at least two points
    for my $i (0 .. $size - 1) {
        return 0 if grep +($_->[$i] // 0) < 2, @$self{qw{ x y z }};
    }
    return 0 unless $self->is_connected;
    # each point must form part of a complete solution
    my %valid;
    my @p = @{ $self->{p} };
    for my $p (@p) {
        next if $valid{"$p->[0]-$p->[1]-$p->[2]"};
        return 0 unless _solvable($size, \%valid, \@p, [], $p);
    }
    return 1;
}

sub is_connected {
    my($self) = @_;

    # check if regions are connected
    my(@r) = @{ $self->{region} };
    my $r0 = pop @r;
    my @pend = ([ 0, $r0 ], [ 1, $r0 ]);
    while (@pend) {
        my($dim, $rp) = @{ pop @pend };
        if ($dim == 0) {
            for my $i (reverse 0 .. $#r) {
                my $ri = $r[$i];
                next if $ri->{b} != $rp->{b};
                push @pend, [ 1, $ri ];
                splice @r, $i, 1;
            }
        } else {
            for my $i (reverse 0 .. $#r) {
                my $ri = $r[$i];
                next if $ri->{s} != $rp->{s};
                push @pend, [ 0, $ri ];
                splice @r, $i, 1;
            }
        }
    }
    return 0 if @r;

    my(@b, @s);
    for my $r (@{ $self->{region} }) {
        push @{ $b[$r->{b}] }, $r;
        push @{ $s[$r->{s}] }, $r;
    }

    # check if rows are connected within each band
    for my $bi (0 .. $#b) {
        my($bmin, $bmax) = @{ $self->{band}[$bi] };
        my $br = $b[$bi];
        next if @$br <= 1;
        my @vec = map {
            my $r = $_;
            my $v = 0;
            for my $p (@{ $r->{p} }) {
                $v |= 1 << ($p->[0] - $bmin);
            }
            $v;
        } @$br;
        my($v, $w) = (pop(@vec), 0);
        while ($v) {
            $w |= $v;
            $v = 0;
            for my $i (reverse 0 .. $#vec) {
                next if ($w & $vec[$i]) == 0;
                $v |= $vec[$i];
                splice @vec, $i, 1;
            }
        }
        return 0 if @vec;
    }

    # check if columns are connected within each stack
    for my $si (0 .. $#s) {
        my($smin, $smax) = @{ $self->{stack}[$si] };
        my $sr = $s[$si];
        next if @$sr <= 1;
        my @vec = map {
            my $r = $_;
            my $v = 0;
            for my $p (@{ $r->{p} }) {
                $v |= 1 << ($p->[1] - $smin);
            }
            $v;
        } @$sr;
        my($v, $w) = (pop @vec, 0);
        while ($v) {
            $w |= $v;
            $v = 0;
            for my $i (reverse 0 .. $#vec) {
                next if ($w & $vec[$i]) == 0;
                $v |= $vec[$i];
                splice @vec, $i, 1;
            }
        }
        return 0 if @vec;
    }
    return 1;
}

sub _canon_mkstr {
    my($self, $x, $y, $xy) = @_;
    my $size = $self->{size};
    my @set = (0 .. $size - 1);
    my $tx = [ @{ +{map +($x->[$_] => $_), @set} }{@set} ];
    my $ty = [ @{ +{map +($y->[$_] => $_), @set} }{@set} ];
    my @out;
    for (@{ $self->{p} }) {
        my($px, $py, $pz) = @$_;
        if ($xy) {
            $out[ $ty->[$py] ][ $tx->[$px] ] = $pz;
        } else {
            $out[ $tx->[$px] ][ $ty->[$py] ] = $pz;
        }
    }
    my @lookup;
    my $next = 0;
    for my $xi (0 .. $size - 1) {
        for my $yi (0 .. $size - 1) {
            my $pz = $out[$xi][$yi] // next;
            $out[$xi][$yi] = $lookup[$pz] //= $next++;
        }
    }
    return join "\xfe", map {
        join '', map {
            defined($_) ? chr($_) : "\xff"
        } @$_[0 .. $size - 1];
    } @out;
}

sub _init_canon {
    my($self) = @_;
    my $size = $self->{size};
    my $best = $self->_canon_mkstr(([ 0 .. $size - 1 ]) x 2, 0);
    my $max = 1;
    my @max;
    for my $r (@{ $self->{region} }) {
        my(@x, @y);
        for my $p (@{ $r->{p} }) {
            ++$x[ $p->[0] ];
            ++$y[ $p->[1] ];
        }
        for (0 .. $size - 1) {
            my $x = $x[$_] or next;
            if ($x > $max) {
                @max = ();
                $max = $x;
            }
            if ($x == $max) {
                push @max, [ 0, $_, $r->{b}, $r->{s} ];
            }
        }
        for (0 .. $size - 1) {
            my $y = $y[$_] or next;
            if ($y > $max) {
                @max = ();
                $max = $y;
            }
            if ($y == $max) {
                push @max, [ 1, $_, $r->{b}, $r->{s} ];
            }
        }
    }
    return {
        best => $best,
        best_at => [ ([ 0 .. $size - 1 ]) x 2, 0 ],
        max => \@max,
        size => $size,
        self => $self,
    };
}

sub _canon_x {
    my($ctx, $v, $x) = @_;
    ++$ctx->{count}[$v + $ctx->{size}];
    my %seen = map +($_ => 1), @$x;
    if (@$x) {
        my $bi = $ctx->{self}{bi}[ $x->[-1] ];
        my $band = $ctx->{self}{band}[$bi];
        my $band_unused = [ grep !$seen{$_}, $band->[0] .. $band->[1] ];
        return $band_unused if @$band_unused;
    } else {
        # require that the first region in the first row can be maximal
        my $y0 = $ctx->{y}[0];
        my $xy = $ctx->{xy};
        if ($xy == 0) {
            # pick any good row that coincides with this stack
            my $si = $ctx->{self}{si}[$y0];
            return [
                sort { $a <=> $b }
                    map $_->[1],
                        grep $_->[0] == 0 && $_->[3] == $si,
                            @{ $ctx->{max} }
            ];
        } else {
            # pick any good band that coincides with this column
            my @bands = map $_->[2],
                grep $_->[0] == 1 && $_->[1] == $y0,
                    @{ $ctx->{max} };
            return [
                sort { $a <=> $b }
                    map $_->[0] .. $_->[1],
                        @{ $ctx->{self}{band} }[ @bands ]
            ];
        }
    }
    return [ grep !$seen{$_}, 0 .. $ctx->{size} - 1 ];
}

sub _canon_y {
    my($ctx, $v, $y) = @_;
    if ($v) {
        $y = [ $ctx->{y}[0], @$y ];
    }
    ++$ctx->{count}[$v + $ctx->{size} * 2];
    my %seen = map +($_ => 1), @$y;
    if (@$y) {
        my $si = $ctx->{self}{si}[ $y->[-1] ];
        my $stack = $ctx->{self}{stack}[$si];
        my $stack_unused = [ grep !$seen{$_}, $stack->[0] .. $stack->[1] ];
        return $stack_unused if @$stack_unused;
    }
    return [ grep !$seen{$_}, 0 .. $ctx->{size} - 1 ];
}

sub _check_canon {
    my($ctx) = @_;
    my $this = $ctx->{self}->_canon_mkstr(@$ctx{qw{ x y xy }});
    if ($this lt $ctx->{best}) {
        $ctx->{best} = $this;
        $ctx->{best_at} = [ @$ctx{qw{ x y xy }} ];
    }
    return;
}

sub _final_canon {
    my($ctx) = @_;
    my($best, $best_at) = @$ctx{qw{ best best_at }};
    return {
        best => $best,
        x => $best_at->[0],
        y => $best_at->[1],
        xy => $best_at->[2],
        count => $ctx->{count},
    };
}

sub canon {
    my($self) = @_;
    return $self->{canon} //= do {
        my $size = $self->{size};
        my $ctx = _init_canon($self);
        my(@cx, @cy);
        for my $v (0 .. $size - 1) {
            push @cx, sub { _canon_x($ctx, $v, [ @_ ]) };
            push @cy, sub { _canon_y($ctx, $v, [ @_ ]) };
        }
        for my $xy (0, 1) {
            local $ctx->{xy} = $xy;
            for my $y0 (@{ _canon_y($ctx, 0, []) }) {
                local $ctx->{y} = [ $y0 ];
                my $iterx = NestedLoops(\@cx);
                while (my @x = $iterx->()) {
                    local $ctx->{x} = \@x;
                    my $itery = NestedLoops([ @cy[1 .. $#cy] ]);
                    while (my @y = $itery->()) {
                        local $ctx->{y} = [ $y0, @y ];
                        _check_canon($ctx);
                    }
                }
            }
        }
        _final_canon($ctx);
    };
}

sub canon_str {
    my($self) = @_;
    my $c = $self->canon;
    my $s = $c->{best};
    $s =~ s{(.)}{
        my $o = ord($1);
        ($o == 254) ? "\n" : ($o == 255) ? '.' : chr(ord('a') + $o)
    }ges;
    return $s;
}

sub canon_px_str {
    my($self) = @_;
    my $c = $self->canon;
    return sprintf 'x=[%s] y=[%s] xy=%s',
        join(' ', @{ $c->{x} }),
        join(' ', @{ $c->{y} }),
        $c->{xy};
}

sub canon_count {
    my($self) = @_;
    my $c = $self->canon;
    return $c->{count}[-1];
}

1;
