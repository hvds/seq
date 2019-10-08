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

For a given C<n> and C<k>, C<

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

sub priority { shift->g->priority }

sub good {
    my($self, $db, $run, $good, $best) = @_;
    if ($run->optm) {
        die "TODO: handle partial result";
    }
    $self->f($good);
    $self->complete(1);
    $self->update;
    printf "f(%s, %s) = %s\n", $self->n, $self->k, $self->f;
    if ($best > $self->k) {
        my $table = $db->resultset($TABLE);
        my @extra = map $table->new({
            n => $self->n,
            k => $_,
            f => $good,
            test_order => '',
        }), $self->k + 1 .. $best;
        $_->complete(1) for @extra;
        $_->insert for @extra;
        printf "f(%s, %s) = %s\n", $_->n, $_->k, $_->f
                for @extra;
    }
    return $self->g->good($db, $best, $good);
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
        next unless $fd->complete && $fd->f;
        my $self;
        my $new = 0;
        if ($self = $f{$fd->k}) {
            next if $self->f;
        } else {
            $self = $table->new({
                n => $n,
                k => $fd->k,
                test_order => '',
            });
            $new = 1;
        }
        $self->f($fd->f * $m);
        $self->depend_m($m);
        $self->depend_n($dn);
        $self->depend(1);
        $self->complete(1);
        $new ? $self->insert : $self->update;
        $ming = $self->k if $ming < $self->k;
    }
    $g->good($db, $ming, $g->checked);
    return;
}

sub _strategy {
    my($self, $db) = @_;
    my $g = $self->g;

    if ($g->checked < $SIMPLE) {
        return Seq::Run->gen(
            $self, $db, {
                optn => $g->checked || 1,
                optx => $SIMPLE,
                optc => 100,
                priority => $g->priority,
            },
        );
    }

    # This must exist if we're not SIMPLE
    my $r = $self->lastRun($db);

    # If we've seen a fix_power once, make sure we set a min value for -c
    # so we don't lose it again.
    if ($r->fix_power && ! $self->minc) {
        return Seq::Run::BisectFP->new($g, $self, $r->optc);
    }

    my $optn = $g->checked + 1;
    my $optx = $g->checked * 2;
    my $_n = sub { "$_[0]" + 0 };
    my $prep = $r->preptime;

    # Note, must use $g->checked for the range, since the last run need
    # not have reached $r->optx.
    my $run = $r->runtime * $_n->($optx + 1 - $optn)
            / $_n->($g->checked + 1 - $r->optn);
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
    $eprep = $prep * $optc / $r->optc;
    $expect = ($eprep + $run) || 1;

    if ($expect < $SLOW) {
        return Seq::Run->gen(
            $self, $db, {
                optn => $optn,
                optx => $optx,
                optc => $optc,
                priority => $g->priority + min(0, -log($expect) / log(2)),
            },
        );
    }

    my $bisect = $self->maybe_bisectg($expect);

    # If it's slow we could try sharding, but for now just reduce the range
    my $factor = (1 + $SLOW / $expect) / 2;
    $optx = _muld($optx, $factor);
    # don't increase priority just because we're reducing the range
    # $expect *= $factor;
    return (Seq::Run->gen(
        $self, $db, {
            optn => $optn,
            optx => $optx,
            optc => $optc,
            optimize => ($expect > $USE_TS && !$g->prime),
            priority => $g->priority + min(0, -log($expect) / log(2)),
        },
    ), ($bisect // ()));
}

sub maybe_bisectg {
    my($self, $expect) = @_;
    my $depth = 1000 * (1 + int(log($expect / $SLOW) / log(10)));
    my $g = $self->g;
    return undef if $depth <= ($g->bisected // 0);
    return Seq::Run::BisectG->new($g, $self, $depth);
}

sub lastRun {
    my($self, $db) = @_;
    my @r = sort { $a->runid <=> $b->runid } $self->runs->all;
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

sub nextFor {
    my($class, $taug, $db) = @_;
    return () if $taug->complete || $taug->depend;
    my $coll = $class->allFor($taug, $db);
    my $self = first { !$_->complete } @$coll;
    if (!$self) {
        my $last = $coll->[-1];
        my $nextk = $last ? $last->k + 1 : 2;
        $self = $db->resultset($TABLE)->new({
            n => $taug->n,
            k => $nextk,
            f => $zero,
            test_order => '',
        });
        $self->minc($last->minc) if $last;
        $self->insert;
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
# When muld is not available, we risk loss of accuracy and overflow.
sub _muld {
    state $have_muld = Math::GMP->can('muld') ? 1 : 0;
    if ($have_muld) {
        return $_[0]->muld($_[1]);
    } else {
        return Math::GMP->new(
            int("$_[0]" * $_[1])
        );
    }
}

1;
