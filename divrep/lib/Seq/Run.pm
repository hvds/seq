package Seq::Run;
use strict;
use warnings;

use Seq::Run::BisectG;
use Seq::Run::BisectFP;
use List::Util qw{ max };

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
    'uint owner',
    'flags(complete running optimizing fix_power old partial) status',
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

sub rprio {
    my($self, $type) = @_;
    return $self->priority;
}

sub logpath {
    my($self, $type) = @_;
    return sprintf '%s/%s.%s-%s',
            $type->logpath, $self->n, $self->k, $self->runid;
}

sub gen {
    my($class, $tauf, $db, $args) = @_;
    my $owner = $db->type->owner;
    my $self = $db->resultset($TABLE)->new({
        n => $tauf->n,
        k => $tauf->k,
        owner => $owner,
        %$args{qw{ optn optx optc optcp optm priority }},
    });
    $self->optimizing(1) if $args->{optimize};
    $self->insert;
    return $self;
}

sub restrategise {
    my($class, $db) = @_;
    my $owner = $db->type->owner;
    for my $self ($db->resultset($TABLE)->search_bitfield(
        { complete => 0 },
    )->search({ owner => $owner })->all) {
        unlink $self->logpath($db->type) if $self->running;
        $self->delete;
    }
}

sub lastForN {
    my($self, $db, $n) = @_;
    my $owner = $db->type->owner;
    return $db->resultset($TABLE)->find({
        n => $n,
        owner => $owner,
    }, {
        order_by => { -desc => [qw{ k runid }] },
        rows => 1,
    });
}

sub command {
    my($self) = @_;
    my $ts = join ',', @{ $self->f->test_order };
    return [
        '-n', '' . $self->optn,
        '-x', '' . $self->optx,
# FIXME: we want them deflated
        ($self->optc ? ('-c', $self->optc) : ()),
        ($self->optcp ? ('-cp', $self->optcp) : ()),
        ($self->optm ? ('-m', $self->optm) : ()),
        ($ts ? ('-ts', $ts)
            : $self->optimizing ? ('-ta') : ()),
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
    my $type = $db->type;
    my $cmd = $self->command;
    my $named = sprintf 'gt(%s,%s)', $self->n, $self->k;
    my $log = $self->logpath($type);
    if ($self->running || -e $log) {
        use Carp;
        warn Carp::longmess(sprintf"already running: %s for %s %s\n",
                $self->runid, $self->n, $self->k);
        return undef;
    }
    my $pid = $type->invoke('gtauseq', $named, $cmd, $log);
    $self->running(1);
    $self->update;
    return $pid;
}

sub failed {
    my($self, $warning) = @_;
    warn $warning;
    $self->running(0);
    $self->update;
    return ();
}

sub parse_ta {
    my($self, $type, $ta) = @_;
    my($base, @num) = split /\s+/, $ta;
    my $which = $type->to_testf($self->k);
    # higher numbers are better at rejecting candidates, so should come first
    return [ sort { $num[$b - 1] <=> $num[$a - 1] } @$which ];
}

sub finalize {
    my($self, $db) = @_;
    my $type = $db->type;
    my $funcname = $type->func_name;
    my $log = $self->logpath($type);
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
    close $fh;

    # required for the func_matches() call below, for some types
    $type->bind($self->n);

    my $ren = qr{\d+(?:e\d+)?};
    my $rend = sub { $_[0] =~ s{e(\d+)$}{0 x $1}er };

    my($good, $bad, $ugly, $depend_m, $depend_n, $fix_power, $last_fail,
            $test_order);
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
        $self->runtime($t - ($self->preptime // 0));
        $good = $rend->($d);
    }
    for (( $line{301} // [] )->[-1] // ()) {
        my($n, $k, $d, $ta) = m{
            ^ 301 \s+ After \s+ [\d\.]+s \s+ for \s+ \( (\d+) ,\s+ (\d+) \)
            \s+ reach \s+ d=(\d+) (?: :LOOSE)? \s+
            \(.*?\) \s+ seen \s+ \[ (.*?) \] \z
        }x or return $self->failed("Can't parse 301 result: '$_'");
        $n == $self->n && $k == $self->k
                or return $self->failed("(n, k) mismatch in '$_'");
        $last_fail = $d;
        if ($self->optimizing) {
            my $order = $self->parse_ta($type, $ta);
            # If no tests were found, it implies all the time is going
            # on modular tests. Set a simple test order in that case,
            # so we don't try to calculate it again.
            $test_order = @$order ? $order : $type->to_testf($self->k);
        }
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
        # 402 Error: all values ... disallowed (${time}s)
        my($t) = m{
            ^ 402 \s+ Error: .*?  \s+ \( (\d+\.\d+) s \) \s* $
        }x or return $self->failed("Can't parse 402 result: '$_'");
        $self->preptime($t);
        $ugly = 1;
    }
    for (@{ $line{403} // [] }) {
        # 403 Error: ... known impossible by exception (${time}s)
        my($t) = m{
            ^ 403 \s+ Error: .*?  \s+ \( (\d+\.\d+) s \) \s* $
        }x or return $self->failed("Can't parse 403 result: '$_'");
        $self->preptime($t);
        $ugly = 1;
    }
    for (@{ $line{404} // [] }) {
        # 404 Error: n + ${k}d must be divisible by n (${time}s)
        my($t) = m{
            ^ 404 \s+ Error: .*?  \s+ \( (\d+\.\d+) s \) \s* $
        }x or return $self->failed("Can't parse 404 result: '$_'");
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
    my $best;
    for (@{ $line{211} // [] }) {
        my($s, $t) = m{
            ^ 211 \s+ Sequence \s+ (\d+) :
            \s+ (\d+) \s+ = \s+ \Q$funcname\E
        }x or return $self->failed("Can't parse 211 result: '$_'");
        $best = $s + 1 if $type->func_matches($s, $t);
    }
    for (@{ $line{311} // [] }) {
        $fix_power = 1;
    }
    if ($good && $best < $self->k) {
        return $self->failed("Inconsistent results: success with best $best");
    }
    unless ($good || $bad || $ugly || $depend_n) {
        return $self->failed("No valid result in $log");
    }

    $self->fix_power(1) if $fix_power;
    $self->complete(1);
    $self->running(0);

    my $tauf = $self->f;
    my @result = $self->partial ? (
        $good ? $tauf->partial($db, $self, $good, $best)
        : ($bad || $ugly || $depend_n) ? ()
        : die "panic"
    ) : (
        $good ? $tauf->good($db, $self, $good, $best)
        : $bad ? do {
            # set the test order before bad() generates the next run
            if ($bad && $test_order && $self->optimizing) {
                $tauf->test_order($test_order);
                $tauf->update;
            }
            my $badm = max(map Math::GMP->new($_),
                    $bad, grep defined, $last_fail);
            $tauf->bad($db, $self, $badm);
        } : $ugly ? $tauf->ugly($db, $self)
        : $depend_n ? $tauf->depends($db, $depend_m, $depend_n)
        : die "panic"
    );

    eval { $self->update; 1 } or do {
        my $e = $@;
        print $self->Dump;
        die $e;
    };
    return @result;
}

1;
