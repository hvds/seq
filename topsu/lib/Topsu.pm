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

sub _canon_mkstr {
    my($self, $x, $y, $z, $xy) = @_;
    my $size = $self->{size};
    my @out;
    for (@{ $self->{p} }) {
        my($px, $py, $pz) = @$_;
        if ($xy) {
            $out[ $y->[$py] ][ $x->[$px] ] = $z->[$pz];
        } else {
            $out[ $x->[$px] ][ $y->[$py] ] = $z->[$pz];
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
    my $best = $self->_canon_mkstr(([ 0 .. $size - 1 ]) x 3, 0);
    return {
        best => $best,
        size => $size,
        self => $self,
    };
}

sub _canon_z {
    my($ctx, $v, $z) = @_;
    ++$ctx->{count}[$v];
    my %seen = map +($_ => 1), @$z;
    return [ grep !$seen{$_}, 0 .. $ctx->{size} - 1 ];
}

sub _canon_x {
    my($ctx, $v, $x) = @_;
    ++$ctx->{count}[$v + $ctx->{size}];
    my %seen = map +($_ => 1), @$x;
    return [ grep !$seen{$_}, 0 .. $ctx->{size} - 1 ];
}

sub _canon_y {
    my($ctx, $v, $y) = @_;
    ++$ctx->{count}[$v + $ctx->{size} * 2];
    my %seen = map +($_ => 1), @$y;
    return [ grep !$seen{$_}, 0 .. $ctx->{size} - 1 ];
}

sub _check_canon {
    my($ctx) = @_;
    my $this = $ctx->{self}->_canon_mkstr(@$ctx{qw{ x y z xy }});
    if ($this lt $ctx->{best}) {
        $ctx->{best} = $this;
        $ctx->{best_at} = [ @$ctx{qw{ x y z xy }} ];
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
        z => $best_at->[2],
        xy => $best_at->[3],
        count => $ctx->{count},
    };
}

sub canon {
    my($self) = @_;
    return $self->{canon} //= do {
        my $size = $self->{size};
        my $ctx = _init_canon($self);
        my(@cx, @cy, @cz);
        for my $v (0 .. $size - 1) {
            push @cx, sub { _canon_x($ctx, $v, [ @_ ]) };
            push @cy, sub { _canon_y($ctx, $v, [ @_ ]) };
            push @cz, sub { _canon_z($ctx, $v, [ @_ ]) };
        }
        my $iterz = NestedLoops(\@cz);
        while (my @z = $iterz->()) {
            local $ctx->{z} = \@z;
            for my $xy (0, 1) {
                local $ctx->{xy} = $xy;
                my $iterx = NestedLoops(\@cx);
                while (my @x = $iterx->()) {
                    local $ctx->{x} = \@x;
                    my $itery = NestedLoops(\@cy);
                    while (my @y = $itery->()) {
                        local $ctx->{y} = \@y;
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
    return sprintf 'x=[%s] y=[%s] z=[%s] xy=%s',
        join(' ', @{ $c->{x} }),
        join(' ', @{ $c->{y} }),
        join(' ', @{ $c->{z} }),
        $c->{xy};
}

sub canon_count {
    my($self) = @_;
    my $c = $self->canon;
    return $c->{count}[-1];
}

1;
