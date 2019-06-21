package Seq::Run::BisectG;
use strict;
use warnings;

my $BISECTG = './bisect-g';
my $LOGS = './logs';

=head1 NAME

Seq::Run::BisectG

=head1 DESCRIPTION

=cut

sub new {
    my($class, $g, $f, $c) = @_;
    return bless {
        g => $g,
        f => $f,
        c => $c,
    }, $class;
}

sub g { shift->{g} }
sub f { shift->{f} }
sub c { shift->{c} }

sub logpath {
    my($self) = @_;
    return sprintf '%s/rbg%s-%s',
            $LOGS, $self->g->n, $self->c;
}

sub priority {
    my($self) = @_;
    return $self->f->priority + 1;
}

sub prep { () }

sub runnable {
    my($self, $db) = @_;
    return () if $self->g->complete
            || $self->g->prime
            || $self->g->depend;
    return $self;
}

sub command {
    my($self) = @_;
    my $g = $self->g;
    return [
        $BISECTG, $g->n, $g->ming, $g->maxg, $g->checked, $self->c,
    ];
}

sub run {
    my($self, $db) = @_;
    my $cmd = $self->command;
    my $log = $self->logpath;
    if (my $pid = fork()) {
        return $pid;
    }
    open STDOUT, '>', $log
            or die "Can't open $log for writing: $!";
    exec @$cmd
            or die "Can't exec [@$cmd]";
}

sub failed {
    my($self, $warning) = @_;
    warn $warning;
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
    my($maxg, $btime);
    for (@{ $line{200} // [] }) {
        (my($n), $maxg, $btime) = m{
            ^ 200 \s+ g\( (\d+) \) \s+ <= \s+ (\d+)
            \s+ \( ([\d\.]+) s\) \z
        }x or return $self->failed("Can't parse 200 line '$_'");
        $n == $self->g->n or return $self->failed("n mismatch in '$_'");
    }
    if ($maxg) {
        $self->g->bisect($db, $maxg, $self->c, $btime);
        return ();
    }
    my $fail = join "\n", map @{ $line{$_} }, grep /^5/, keys %line;
    return $self->failed("bisect failed: $_", $fail // 'unknown cause');
}

1;
