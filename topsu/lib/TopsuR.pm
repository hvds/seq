package TopsuR;
use strict;
use warnings;

use Algorithm::Loops qw{ NextPermuteNum };
use List::Util qw{ max };

sub new {
    my($class, $band, $stack, $region) = @_;
    my $self = bless {
        band => $band,
        stack => $stack,
        region => $region,
    }, $class;
    $self->init;
    return $self;
}

sub init {
    my($self) = @_;
    my @in;
    for my $r (@{ $self->{region} }) {
        $in[ $r->[0] ][ $r->[1] ] = 1;
    }
    $self->{in} = \@in;
    return;
}

sub _canon_mkstr {
    my($self, $b, $s, $xy) = @_;
    my($band, $stack, $region) = @$self{qw{ band stack region }};
    my $tb = [ @{ +{map +($b->[$_] => $_), 0 .. $#$b} }{0 .. $#$b} ];
    my $ts = [ @{ +{map +($s->[$_] => $_), 0 .. $#$s} }{0 .. $#$s} ];
    my @str = $xy ? (('.' x @$b) x @$s) : (('.' x @$s) x @$b);
    for my $ri (0 .. $#$region) {
        my $r = $region->[$ri];
        my($b, $s) = ($tb->[ $r->[0] ], $ts->[ $r->[1] ]);
        if ($xy) {
            substr($str[$s], $b, 1) = chr(ord('A') - 1 + $band->[$b]);
        } else {
            substr($str[$b], $s, 1) = chr(ord('A') - 1 + $stack->[$s]);
        }
    }
    return join ":", @str;
}

sub as_str {
    my($self) = @_;
    return $self->{str} //= $self->_canon_mkstr(
        [ 0 .. $#{ $self->{band} } ], [ 0 .. $#{ $self->{stack} } ], 0
    );
}

sub canon {
    my($self) = @_;
    return $self->{canon} //= do {
        my($band, $stack, $in, $region) = @$self{qw{ band stack in region }};
        my($best, @best_at) = ('');
        my($max, @max) = 0;
        for my $r (@$region) {
            my($bi, $si) = @$r;
            my($bw, $sw) = ($band->[$bi], $stack->[$si]);
            if ($bw > $max) {
                $max = $bw;
                @max = ();
            }
            if ($sw > $max) {
                $max = $sw;
                @max = ();
            }
            push @max, [ 1, $bi, $si ] if $bw == $max;
            push @max, [ 0, $bi, $si ] if $sw == $max;
        }
        for (@max) {
            my($xy, $b0, $s0) = @$_;
            if ($xy) {
                my(@gba, @gbb);
                for my $bi (0 .. $#$band) {
                    next if $bi == $b0;
                    if ($in->[$bi][$s0]) {
                        push @{ $gba[$band->[$bi]] }, $bi;
                    } else {
                        push @{ $gbb[$band->[$bi]] }, $bi;
                    }
                }
                my $iter = _banded_permute([
                    [ $b0 ], grep defined, reverse(@gba), reverse(@gbb)
                ]);
                my $sx = [ $s0, grep $_ != $s0, 0 .. $#$stack ];
                while (my $bx = $iter->()) {
                    my $str = join ':',
                        sort { $b cmp $a }
                            split /:/,
                                $self->_canon_mkstr($bx, $sx, $xy);
                    next unless $str gt $best;
                    $best = $str;
                    @best_at = ($bx, $sx, $xy);
                }
            } else {
                my(@gsa, @gsb);
                for my $si (0 .. $#$stack) {
                    next if $si == $s0;
                    if ($in->[$b0][$si]) {
                        push @{ $gsa[$stack->[$si]] }, $si;
                    } else {
                        push @{ $gsb[$stack->[$si]] }, $si;
                    }
                }
                my $iter = _banded_permute([
                    [ $s0 ], grep defined, reverse(@gsa), reverse(@gsb)
                ]);
                my $bx = [ $b0, grep $_ != $b0, 0 .. $#$band ];
                while (my $sx = $iter->()) {
                    my $str = join ':',
                        sort { $b cmp $a }
                            split /:/,
                                $self->_canon_mkstr($bx, $sx, $xy);
                    next unless $str gt $best;
                    $best = $str;
                    @best_at = ($bx, $sx, $xy);
                }
            }
        }
        {
            best => $best,
            b => $best_at[0],
            s => $best_at[1],
            xy => $best_at[2],
        };
    };
}

sub canon_str {
    return shift->canon->{best};
}

sub _banded_permute {
    my @set = map [ sort { $a <=> $b } @$_ ], @{ $_[0] };
    my @iter = grep @$_ > 1, @set;
    my $done = 0;
    return sub {
        return undef if $done;
        my $ret = [ map @$_, @set ];
        my $ic = $#iter;
        while ($ic >= 0) {
            last if NextPermuteNum(@{ $iter[$ic] });
            --$ic;
        }
        $done = 1 if $ic < 0;
        return $ret;
    }
}

1;
