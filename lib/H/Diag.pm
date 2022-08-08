package H::Diag;
use strict;
use warnings;

use Exporter qw{import};
our @EXPORT_OK = qw{ diag keep_diag };

=head1 NAME

H::Diag - Simple progress diagnostics

=head1 DESCRIPTION

Provides a mechanism to show progress diagnostics by updating a single
line.

Note that on loading, this module attempts to install a C< $SIG{TSTP} >
handler such that when the process is paused it will C<keep_diag()>.
This ensures that when resumed, further diagnostics appear as they should.
If the system has no C<TSTP> signal, this will die on loading; in that
case it should be safe to catch the death and use C<diag()> and
C<keep_diag()> normally, but now without the benefit of the handler.

=head1 FUNCTIONS

=over 4

=item diag ( string )

Replace the current progress diagnostic with the supplied string.

=item keep_diag ( )

Keep the current progress diagnostic (if any); future progress will
be shown on a new line.

=back

=cut

my $s = '';

sub diag {
    # FIXME: can we do this in one go with something like "^A^K"?
    print "\x08 \x08" x length($s);
    ($s) = @_;
    local $| = 1;
    print $s;
    return;
}

sub keep_diag {
    if (length $s) {
        $s = '';
        local $| = 1;
        print "\n";
    }
    return;
}

sub _init {
    # Adapted from IPC::Signal
    my $tstp_num = do {
        use Config;
        my @num = split ' ', $Config::Config{'sig_num'};
        my @name = split ' ', $Config::Config{'sig_name'};
        my %map;
        @map{@name} = @num;
        $map{'TSTP'};
    } // die 'No signum found for TSTP';

    # Adapted from IPC::Signal::Force
    use POSIX qw(SIG_SETMASK SIG_UNBLOCK sigprocmask);
    $SIG{'TSTP'} = sub {
        keep_diag();
        my $sigset = new POSIX::SigSet $tstp_num;
        my $mask = new POSIX::SigSet;
        local $SIG{'TSTP'};
        kill 'TSTP', 0;
        sigprocmask(SIG_UNBLOCK, $sigset, $mask);
        if($mask->ismember($tstp_num)) {
            sigprocmask(SIG_SETMASK, $mask, $mask);
        }
    };
}

_init();
1;
