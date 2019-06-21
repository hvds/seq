package Seq::Run::BisectFP;
use strict;
use warnings;

my $BISECTFP = './bisect-fp';
my $LOGS = './logs';

=head1 NAME

Seq::Run::BisectFP

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
    my $f = $self->f;
    return sprintf '%s/rbfp%s.%s-%s',
            $LOGS, $f->n, $f->k, $self->c;
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
            || $self->g->depend
            || $self->f->minc;
    return $self;
}

sub command {
    my($self) = @_;
    my $g = $self->g;
    my $f = $self->f;
    return [
        $BISECTFP, $f->n, $f->k, $g->checked, $self->c,
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
    my $f = $self->f;
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
    my $minc;
    for (@{ $line{200} // [] }) {
        (my($n, $k), $minc) = m{
            ^ 200 \s+ f\( (\d+) ,\s+ (\d+) \) \s+ needs \s+ -c \s+ (\d+)
            \s+ to \s+ find \s+ fix_power
        }x or return $self->failed("Can't parse 200 line '$_'");
        $n == $f->n && $k == $f->k
                or return $self->failed("(n, k) mismatch in '$_'");
    }
    if ($minc) {
        return $self->f->set_minc($db, $minc);
    }
    my $fail = join "\n", map @{ $line{$_} }, grep /^5/, keys %line;
    return $self->failed("bisect failed: $_", $fail // 'unknown cause');
}

1;
