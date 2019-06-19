package Seq::Run;
use strict;
use warnings;

my $PROG = './gtauseq';
my $LOGS = './logs';

=head1 NAME

Seq::Run

=head1 DESCRIPTION

=cut

use parent 'Seq::Table';
my $TABLE = 'Run';
__PACKAGE__->define($TABLE, 'run', [
    'key id runid',
    'uint n',
    'uint k',
    'flags(complete running optimizing) status',
    'bigint optn',
    'bigint optx',
    'uint optc',
    'maybe uint optcp',
    'maybe modlist optm',
    'maybe float preptime',
    'maybe float runtime',
    'float priority',
]);
__PACKAGE__->belongs_to(
    f => 'Seq::TauF', {
        'foreign.n' => 'self.n',
        'foreign.k' => 'self.k',
    },
);

sub logpath {
    my($self) = @_;
    return sprintf '%s/%s.%s-%s',
            $LOGS, $self->n, $self->k, $self->runid;
}

sub gen {
    my($class, $tauf, $db, $args) = @_;
    my $self = $db->resultset($TABLE)->new({
        n => $tauf->n,
        k => $tauf->k,
        %$args{qw{ optn optx optc optcp optm priority }},
    });
    $self->optimizing(1) if $args->{optimizing};
    $self->insert;
    return $self;
}

sub restrategise {
    my($class, $db) = @_;
    for my $self ($db->resultset($TABLE)->search_bitfield(
        { complete => 0 },
    )->all) {
        unlink $self->logpath if $self->running;
        $self->delete;
    }
}

sub lastForN {
    my($self, $db, $n) = @_;
    return $db->resultset($TABLE)->find({
        n => $n,
    }, {
        order_by => { -desc => [qw{ k runid }] },
        rows => 1,
    });
}

sub command {
    my($self) = @_;
    return [
        $PROG,
        '-n', '' . $self->optn,
        '-x', '' . $self->optx,
# FIXME: we want them deflated
        ($self->optc ? ('-c', $self->optc) : ()),
        ($self->optcp ? ('-cp', $self->optcp) : ()),
        ($self->optm ? ('-m', $self->optm) : ()),
        $self->n,
        $self->k,
    ];
}

sub prep {
    my($self, $db) = @_;
    return () if $self->complete;
    return $self->finalize($db) if $self->running;
    return ();  # ready to run
}

sub runnable {
    my($self, $db) = @_;
    return () if $self->complete || $self->running;
    return $self;
}

sub run {
    my($self, $db) = @_;
    my $cmd = $self->command;
    my $log = $self->logpath;
    if (my $pid = fork()) {
        $self->running(1);
        $self->update;
        return $pid;
    }
    open STDOUT, '>', $log
            or die "Can't open $log for writing: $!";
    exec @$cmd
            or die "Can't exec [@$cmd]";
}

sub failed {
    my($self, $warning) = @_;
    $self->running(0);
    $self->update;
    return ();
}

sub finalize {
    my($self, $db) = @_;
    my $log = $self->logpath;
    my $fh;
    open($fh, '<', $log)
            or return $self->failed("Can't open $log for reading: $!");
    my %line;
    while (<$fh>) {
        chomp;
        my($rc) = /^(\d{3}) /
                or return $self->failed("Can't parse log line '$_'");
        push @{ $line{$rc} }, $_;
    }

    my $ren = qr{\d+(?:e\d+)?};
    my $rend = sub { $_[0] =~ s{e(\d+)$}{0 x length($1)}er };

    my($good, $bad, $ugly, $depend_m, $depend_n);
    for (@{ $line{309} // [] }) {
        /\((\d+\.\d*)s\)$/ && $self->preptime($1);
    }
    for (@{ $line{200} // [] }) {
        my($n, $k, $d, $t) = m{
            ^ 200 \s+ f\( (\d+) ,\s+ (\d+) \)
            \s+ = \s+ ($ren) \s+
            \( (\d+\.\d+) s \)
            \s* $
        }x or return $self->failed("Can't parse 200 result: '$_'");
        $n == $self->n && $k == $self->k
                or return $self->failed("(n, k) mismatch in '$_'");
        $self->runtime($t - $self->preptime);
        $good = $rend->($d);
    }
    for (@{ $line{500} // [] }) {
        my($n, $k, $d, $t) = m{
            ^ 500 \s+ f\( (\d+) ,\s+ (\d+) \)
            \s+ > \s+ ($ren) \s+
            \( (\d+\.\d+) s \)
            \s* $
        }x or return $self->failed("Can't parse 500 result: '$_'");
        $n == $self->n && $k == $self->k
                or return $self->failed("(n, k) mismatch in '$_'");
        $self->runtime($t - $self->preptime);
        $bad = $rend->($d);
    }
    for (@{ $line{402} // [] }) {
        my($t) = m{
            ^ 402 \s+ Error: \s+ all \s+ values \s+ \(.*?\) \s+ disallowed
            \s+ \( (\d+\.\d+) s \) \s* $
        }x or return $self->failed("Can't parse 402 result: '$_'");
        $self->preptime($t);
        $ugly = 1;
    }
    for (@{ $line{201} // [] }) {
        my($k, $dm, $dn) = m{
            ^ 201 \s+ Dependent \s+ - \s+
            ($ren) : (\d+) \s+ f\( (\d+) \)
        }x or return $self->failed("Can't parse 201 result: '$_'");
        $k == $self->k or return $self->failed("(n, k) mismatch in '$_'");
        ($depend_m, $depend_n) = ($rend->($dm), $dn);
    }
    my($best, $tau);
    for (@{ $line{211} // [] }) {
        my($s, $t) = m{
            ^ 211 \s+ Sequence \s+ (\d+) :
            \s+ (\d+) \s+ = \s+ tau
        }x or return $self->failed("Can't parse 211 result: '$_'");
        if ($s == 0) {
            $best = 1;
            $tau = $t;
        } else {
            ++$best if $best == $s && $tau == $t;
        }
    }
    if ($good && $best < $self->k) {
        return $self->failed("Inconsistent results: success with best $best");
    }
    unless ($good || $bad || $ugly || $depend_n) {
        return $self->failed("No valid result in $log");
    }

    $self->complete(1);
    $self->update;

    return $self->f->good($db, $self, $good, $best) if $good;
    return $self->f->bad($db, $self, $bad) if $bad;
    return $self->f->ugly($db, $self) if $ugly;
    return $self->f->depends($db, $depend_m, $depend_n) if $depend_n;
    die "panic";
}

1;
