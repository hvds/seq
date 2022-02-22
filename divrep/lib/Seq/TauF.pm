package Seq::TauF;
use strict;
use warnings;
use feature qw{state};

use Math::GMP;
use List::Util qw{ first min max };

my $zero = Math::GMP->new('0');
# Assume we don't need to do anything clever to check values to this limit.
my $SIMPLE = Math::GMP->new(1 << 20);
# Assume we don't need to do much clever if expected runtime is this fast.
my $FAST = 10;
my $SLOW = 600;
my $USE_TS = 1200;
my $EXPECT_TS = 1200;

=head1 NAME

Seq::TauF

=head1 DESCRIPTION

C<f(n, k)> is defined as the smallest difference C<d> such that all of
C<< n + xd, 0 <= d < k >> have the same number of factors. Knowledge
of this function is recorded in the database table C<tauf> described
by the embedded C<Seq::TauF::Schema> class; results from this table
will generate objects of this class.

For a given C<n> and C<k>, C<f> is the smallest known difference
satisfying the criterion (but is not necessarily known to be minimal:
see the C<status> flags below). This is stored as a string representing
a L<Math::GMP> bigint.

=cut

use parent 'Seq::Table';
my($TABLE, $tauf) = ('TauF', __PACKAGE__);
$tauf->define($TABLE, 'tauf', [
    'key uint n',
    'key uint k',
    'maybe bigint f',
    'maybe bigint depend_m',
    'maybe uint depend_n',
    'maybe uint minc',
    'modlist test_order',
    'flags(complete external estimated depend impossible) status',
]);
$tauf->belongs_to(
    g => 'Seq::TauG', {
        'foreign.n' => 'self.n',
    },
);
$tauf->has_many(
    runs => 'Seq::Run', {
        'foreign.n' => 'self.n',
        'foreign.k' => 'self.k',
    },
    {
        order_by => 'runid',
        # It is useful to keep them for debugging, but may need to reconsider
        cascade_delete => 0,
    },
);

sub rprio {
    my($self, $type) = @_;
    return $type->gprio($self->n);
}

sub good {
    my($self, $db, $run, $good, $best) = @_;
    $self->f($good) if !$self->f || $self->f > $good;
    $self->complete(1);
    $self->update;
    printf "f(%s, %s) = %s\n", $self->n, $self->k, $self->f;
    if ($best > $self->k) {
        my $next = Seq::TauF->forceFor($self->g, $db, $self->k + 1);
        return $next->good($db, $run, $good, $best);
    }
    return $self->g->good($db, $best, $good);
}

sub _partial {
    my($self, $db, $good) = @_;
    if (!$self->f || $self->f > $good) {
        $self->f($good);
        $self->update;
        printf "f(%s, %s) <= %s\n", $self->n, $self->k, $self->f;
        return 1;
    }
    return 0;
}

sub partial {
    my($self, $db, $run, $good, $best) = @_;
    my $g = $self->g;
    for my $k ($db->type->ming + 1 .. $best) {
        my $this = ($k == $self->k) ? $self : Seq::TauF->forceFor($g, $db, $k);
        $this->_partial($db, $good);
    }
    return $g->partial($db, $best);
}

sub ugly {
    my($self, $db, $run) = @_;
    $self->complete(1);
    $self->impossible(1);
    $self->update;
    return $self->g->ugly($db, $self->k);
}

sub bad {
    my($self, $db, $run, $bad) = @_;
    return $self->g->bad($db, $bad);
}

sub set_minc {
    my($self, $db, $minc) = @_;
    $self->minc($minc);
    $self->update;
    return ($self);
}

sub depends {
    my($self, $db, $depend_m, $depend_n) = @_;
# FIXME: g(depend_n).max may already be less than our k, in which case now
# is probably our only opportunity to set g(n).max
    $self->depend(1);
    $self->depend_m($depend_m);
    $self->depend_n($depend_n);
    $self->complete(1);
    $self->update;
    return $self->g->depends($db, $self->k - 1, $depend_n);
}

sub update_depends {
    my($class, $db, $g) = @_;
    my $n = $g->n;
    my $dn = $g->depend_n;
    my $table = $db->resultset($TABLE);
    my @f = $g->f->all;
    my %f = map +($_->k => $_), @f;
    my $m = (first { $_->depend_m } @f)->depend_m;
    my $gd = $g->depended;
    my $ming = $g->ming;

    for my $fd ($gd->f->all) {
        next unless $fd->f;
        my $self;
        if ($self = $f{$fd->k}) {
            if ($fd->complete) {
                if (!$self->f || $self->f > $fd->f) {
                    $self->f($fd->f * $m);
                    $self->depend_m($m);
                    $self->depend_n($dn);
                    $self->depend(1);
                    $self->complete(1);
                } elsif (!$self->complete) {
                    $self->complete(1);
                }
            } else {
                next if $self->complete;
                if (!$self->f || $self->f > $fd->f) {
                    $self->f($fd->f * $m);
                    $self->depend_m($m);
                    $self->depend_n($dn);
                    $self->depend(1);
                }
            }
            $self->update;
        } else {
            $self = $table->new({
                n => $n,
                k => $fd->k,
                test_order => '',
            });
            $self->f($fd->f * $m);
            $self->depend_m($m);
            $self->depend_n($dn);
            $self->depend(1);
            $self->complete($fd->complete);
            $self->insert;
        }
        # FIXME: see TauG->partial
        $ming = $self->k if $self->complete && $ming < $self->k;
    }
    # FIXME: maybe pass max($g->checked, $gd->checked * $m)
    $g->good($db, $ming, $g->checked);
    return;
}

sub _strategy {
    my($self, $db) = @_;
    my $type = $db->type;
    my $g = $self->g;

    if ($g->checked < $SIMPLE) {
        return Seq::Run->gen(
            $self, $db, {
                optn => $g->checked || 1,
                optx => $SIMPLE,
                optc => 100,
                priority => $type->fprio($self->n, $self->k, 0),
            },
        );
    }

    # This must exist if we're not SIMPLE
    my $r = $self->lastRun($db);

    # If we've seen a fix_power once, make sure we set a min value for -c
    # so we don't lose it again.
    if ($r->fix_power && ! $self->minc) {
        return Seq::Run::BisectFP->new($db->type, $g, $self, $r->optc);
    }

    my $optn = $g->checked + 1;
    my $optx = $g->checked * 2;
    my $_n = sub { "$_[0]" + 0 };
    # last run may not have completed, in which case the actual point reached
    # will be in g->checked; but the latter may also have been advanced by
    # someone else, so take optx as the safer guess if it's the smaller
    my $prevrange = $_n->(min($g->checked, $r->optx) - $r->optn);
    my $prep = $r->preptime;
    my $run = $r->runtime * $_n->($optx + 1 - $optn) / $prevrange;
    my $expect = ($prep + $run) || 1;

    # decide what -c value to supply
    my $optc;
    if ($expect < $FAST) {
        $optc = max(100, $r->optc / 2);
    } else {
        $optc = 100 + $optx->sizeinbase_gmp(2) * 50;
        $optc /= 4 if $r->fix_power;
    }
    my $eprep = $prep * $optc / $r->optc;
    $optc = ($optc + $r->optc) / 2 if (
        ($optc > $r->optc && $prep > $run)
        || ($optc < $r->optc && $prep < $run / 10)
    );
    $optc = max($optc, $self->minc // 0);
    $optc = int($optc);
    $eprep = $prep * $optc / $r->optc;
    $expect = ($eprep + $run) || 1;

    if ($expect < $SLOW) {
        return Seq::Run->gen(
            $self, $db, {
                optn => $optn,
                optx => $optx,
                optc => $optc,
                priority => $type->fprio($self->n, $self->k, $expect),
            },
        );
    }

    my $bisect = $self->maybe_bisectg($db->type, $expect);

    # If it's slow we could try sharding, but for now just reduce the range
    my $factor = (1 + $SLOW / $expect) / 2;
    $optx = _bmulf($optx, $factor);
    # don't increase priority just because we're reducing the range
    # $expect *= $factor;

    # we want a special run to optimize test order if:
    # a) we're slow enough (USE_TS)
    # b) it would be useful (n is not prime)
    # c) we haven't already done so ($self->test_order)
    # and d) we have at least one run with this k
    if ($expect > $USE_TS
        && !$g->prime
        && !@{ $self->test_order }
        && $r->k == $self->k
    ) {
        # test order runs are slow, we need to scale our range down
        my $log = $r->logpath($type);
        my $fh;
        open($fh, '<', $log)
                or return $self->failed("Can't open $log for optimizing: $!");
        my $last;
        while (<$fh>) {
            chomp;
            my($rc) = /^(\d{3}) /
                    or return $self->failed("Can't parse log line '$_'");
            $last = $_ if $rc == 301;
        }
        close $fh;
        defined($last) or die "$log: no 301 line found";
        my($keep, $test) = $last =~ m{
            ^ 301 \s+ .*?
            keep \s+ (\d+) .*?
            seen \s+ \[ (.*?) \] \z
        }x or die "Can't parse 301 result: '$last'";
        my($sum, $i) = (0, 1);
        $sum += $_ * $i++ for split /\s+/, $test;

        # a normal run took $sum tests for $keep values; an optimizing
        # run will take $self->k * $keep tests, so scale it down,
        # with care to avoid mixing bigints and floats
        $optx = $optn + int(('' . ($optx - $optn)) * ($sum / $self->k / $keep)
                * ($EXPECT_TS / $SLOW));
        return (Seq::Run->gen(
            $self, $db, {
                optn => $optn,
                optx => $optx,
                optc => $optc,
                optimize => 1,
                priority => $type->fprio($self->n, $self->k, $expect),
            },
        ), ($bisect // ()));
    }

    return (Seq::Run->gen(
        $self, $db, {
            optn => $optn,
            optx => $optx,
            optc => $optc,
            optimize => 0,
            priority => $type->fprio($self->n, $self->k, $expect),
        },
    ), ($bisect // ()));
}

sub maybe_bisectg {
    my($self, $type, $expect) = @_;
    my $depth = 1000 * (1 + int(log($expect / $SLOW) / log(10)));
    my $g = $self->g;
    return undef if $depth <= ($g->bisected // 0);
    return Seq::Run::BisectG->new($type, $g, $self, $depth);
}

sub lastRun {
    my($self, $db) = @_;
    my @r = sort { $a->runid <=> $b->runid }
            grep !$_->optimizing, $self->runs->all;
    return $r[-1] if @r;
    my $r = Seq::Run->lastForN($db, $self->n);
    return $r if $r;
print $self->Dump;
die "No last run found";
}

sub prep {
    my($self, $db) = @_;
    my @r = grep !$_->complete, $self->runs->all;
    my @ready = grep $_->running, @r;
    return @ready if @ready;
    return @r if @r;    # not complete or running, so ready to run
    return $self->_strategy($db);
}

sub runnable { return () }

sub forceFor {
    my($class, $taug, $db, $k) = @_;
    my $n = $taug->n;
    my $self = $db->resultset($TABLE)->find({ n => $n, k => $k });
    unless ($self) {
        die "g($n) <= @{[ $taug->maxg ]}, not forcing k=$k"
                unless $k <= $taug->maxg;
        $self = $db->resultset($TABLE)->new({
            n => $n,
            k => $k,
            f => $zero,
            test_order => '',
        });
        $self->insert;
    }
    return $self;
}

sub nextFor {
    my($class, $taug, $db) = @_;
    return () if $taug->complete || $taug->depend;
    my $coll = $class->allFor($taug, $db);
    my $self = first { !$_->complete } @$coll;
    if (!$self) {
        my $last = $coll->[-1];
        my $nextk = ($last ? $last->k : $taug->ming) + 1;
        $self = $class->forceFor($taug, $db, $nextk);
        if ($last && $last->minc) {
            $self->minc($last->minc);
            $self->update;
        }
    }
    return $self;
}

sub allFor {
    my($class, $taug, $db) = @_;
    return [ $db->resultset($TABLE)->search(
        { n => $taug->n },
        { order_by => 'k' },
    )->all ];
}

# Multiply bigint first arg by NV second arg, returning bigint.
# When bmulf is not available, we risk loss of accuracy and overflow.
sub _bmulf {
    state $have_bmulf = Math::GMP->can('bmulf') ? 1 : 0;
    if ($have_bmulf) {
        return $_[0]->bmulf($_[1]);
    } else {
        my $r = int("$_[0]" * $_[1]);
        $r = sprintf '%.0f', $r if $r =~ /e/;
        return Math::GMP->new($r);
    }
}

1;
