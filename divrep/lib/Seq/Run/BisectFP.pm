package Seq::Run::BisectFP;
use strict;
use warnings;

=head1 NAME

Seq::Run::BisectFP

=head1 DESCRIPTION

=cut

sub new {
    my($class, $type, $g, $f, $c) = @_;
    return bless {
        type => $type,
        g => $g,
        f => $f,
        c => $c,
    }, $class;
}

sub type { shift->{type} }
sub g { shift->{g} }
sub n { shift->g->n }
sub f { shift->{f} }
sub c { shift->{c} }

sub logpath {
    my($self) = @_;
    my $f = $self->f;
    return sprintf '%s/rbfp%s.%s-%s',
            $self->type->logpath, $f->n, $f->k, $self->c;
}

sub rprio {
    my($self, $type) = @_;
    return $type->gprio($self->n) + 1;
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
    return [ $f->n, $f->k, $g->checked, $self->c ];
}

sub run {
    my($self, $db) = @_;
    my $cmd = $self->command;
    my $named = sprintf 'bfp(%s, %s)', $self->n, $self->f->k;
    my $log = $self->logpath;
    return $self->type->invoke('bisect_fp', $named, $cmd, $log);
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
    for (@{ $line{402} // [] }) {
        (my($n, $k), $minc) = m{
            ^ 402 \s+ f\( (\d+) ,\s+ (\d+) \) \s+ all \s+ values \s+
            disallowed \s+ with \s+ -c \s+ (\d+) \d+
        }x or return $self->failed("Can't parse 402 line '$_'");
        $n == $f->n && $k == $f->k
                or return $self->failed("(n, k) mismatch in '$_'");
        # We'll let the harness rediscover this for itself, to get
        # a better trail of proof; for now we'll just let this force
        # minc for the subsequent attempt.
        # FALLTHROUGH
    }

    if ($minc) {
        return $self->f->set_minc($db, $minc);
    }
    my $fail = join "\n", map @{ $line{$_} }, grep /^5/, keys %line;
    return $self->failed("bisect failed: %s", $fail // 'unknown cause');
}

1;
