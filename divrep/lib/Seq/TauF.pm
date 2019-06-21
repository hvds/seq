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
    'flags(complete external estimated depend) status',
    'float priority',
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
);

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
            priority => 0,
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
    $self->update;
    return $self->g->ugly($db, $self->k);
}

sub bad {
    my($self, $db, $run, $bad) = @_;
    return $self->g->bad($db, $bad);
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
    my $table = $db->resultset($TABLE);
    my $self = ($g->f->all)[-1];
    my $depending = $table->search(
        {
            n => $self->depend_n,
            k => { '>=', $self->k },
        },
        { order_by => 'k' },
    )->search_bitfield({ complete => 1 });
    for my $d ($depending->all) {
        my $this = ($d->k == $self->k) ? $self : $table->new({
            n => $self->n,
            k => $d->k,
        });
        $this->f($d->f * $self->depend_m);
        $this->depend_m($self->depend_m);
        $this->depend_n($self->depend_n);
        $this->priority(0);
        $this->depend(1);
        $this->complete(1);
        ($this == $self) ? $this->update : $this->insert;
        $g->good($db, $this->k, $this->f);
    }
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
    my $r = $self->lastRun($db);
    my $optn = $g->checked + 1;
    my $optx = $g->checked * 2;
    my $_n = sub { "$_[0]" + 0 };
    my $prep = $r->preptime;
    my $run = $r->runtime * $_n->($optx + 1 - $optn)
            / $_n->($r->optx + 1 - $r->optn);
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
        my $next = @$coll ? $coll->[-1]->k + 1 : 2;
        $self = $db->resultset($TABLE)->new({
            n => $taug->n,
            k => $next,
            f => $zero,
        });
        # Hack: hardcode this until we can determine automatically for fix_power
        if ($self->n == 72) {
            $self->minc(200);
        }
        $self->priority($taug->priority);
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
