package Board;
use strict;
use warnings;

use Group;

our $count = 0;
our $feedback = 4;

sub new {
    my($class, $k, $unused, $groups) = @_;
    return bless {
        k => $k,
        unused => $unused,
        groups => $groups,
    }, $class;
}
sub k { shift->{k} }
sub unused { shift->{unused} }
sub groups { shift->{groups} }

sub str {
    my($self) = @_;
    return join '; ', map $_->str, @{ $self->groups };
}

sub recurse {
    my($self, $next_unused, $next_groups, $best, $best_board) = @_;
    my $k = $self->k;
    my $next = Board->new($k + 1, $next_unused, $next_groups);

    # record count of boards recursed through
    ++$count;

    # give some feedback of progress
    print $next->str, "\n" if $k <= $feedback;

    my($next_best, $next_best_board) = $next->try($best, $best_board);
    return +($next_best > $best ? ($next_best, $next_best_board)
        : $k > $best ? ($k, $next)
        : ($best, $best_board)
    );
}

sub try {
    my($self, $best, $best_board) = @_;
    my $k = $self->k;
    my $unused = $self->unused;
    my $groups = $self->groups;

    # try making a new group
    for (my $mult = 1; $mult * $k <= $unused; ++$mult) {
        my $next_unused = $unused - $k * $mult;
        for my $group (Group->new_seed($k)) {
            my $next_groups = [ @$groups, $group ];
            ($best, $best_board) = $self->recurse(
                $next_unused, $next_groups, $best, $best_board,
            );
        }
    }

    # try extending an existing group
    for my $i (0 .. $#$groups) {
        my $gi = $groups->[$i];
        my $sums = $gi->sums;
        my $maxsum = $#$sums;

        # try exact first
        for my $mult (1 .. int($maxsum / $k)) {
            for my $locn (@{ $sums->[$k * $mult] // [] }) {
                my $new_group = $gi->place($locn, $k);
                my $next_groups = [ @$groups ];
                $next_groups->[$i] = $new_group;
                ($best, $best_board) = $self->recurse(
                    $unused, $next_groups, $best, $best_board,
                );
            }
        }

        # try with unused
        for my $mult (1 .. int(($maxsum + $unused) / $k)) {
            for my $diff (1 .. $unused) {
                my $sum = $k * $mult - $diff;
                last if $sum < 0;
                my $locns = $sums->[$sum] // next;
                my $next_unused = $unused - $diff;
                for my $locn (@$locns) {
                    for my $new_group ($gi->place_with($locn, $k, $diff)) {
                        my $next_groups = [ @$groups ];
                        $next_groups->[$i] = $new_group;
                        ($best, $best_board) = $self->recurse(
                            $next_unused, $next_groups, $best, $best_board,
                        );
                    }
                }
            }
        }

        # try by coalesce
        for my $j ($i + 1 .. $#$groups) {
            my $gj = $groups->[$j];
            my $sumsj = $gj->sums;
            my $maxsumj = $#$sumsj;
            for my $mult (1 .. int(($maxsum + $maxsumj + $unused) / $k)) {
                my $targ = $k + $mult;
                for my $sumi (1 .. $targ - 1) {
                    my $locnsi = $sums->[$sumi] // next;
                    my $need = $targ - $sumi;
                    my $min = $need - $unused;
                    $min = 1 if $min < 1;
                    for my $sumj ($min .. $need) {
                        my $locnsj = $sumsj->[$sumj] // next;
                        my $use = $need - $sumj;
                        my $next_unused = $unused - $use;
                        for my $locni (@$locnsi) {
                            for my $locnj (@$locnsj) {
                                for my $new_group ($gi->coalesce(
                                    $locni, $gj, $locnj, $k, $use,
                                )) {
                                    my $next_groups = [ @$groups ];
                                    $next_groups->[$i] = $new_group;
                                    splice @$next_groups, $j, 1;
                                    ($best, $best_board) = $self->recurse(
                                        $next_unused, $next_groups,
                                        $best, $best_board,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return ($best, $best_board);
}

1;
