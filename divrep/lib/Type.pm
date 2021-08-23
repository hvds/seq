package Type;
use strict;
use warnings;
use List::Util ();

my %types = (
    t => 'Type::TauSeq',
    a => 'Type::AddSeq',
    s => 'Type::Semip',
);
my %typename = reverse %types;

my $LOGS = './logs';
my %PROG = (
    gtauseq => './gtauseq',
    bisect_fp => './bisect-fp',
    bisect_g => './bisect-g',
    find_shard => './find-shard',
);

sub new {
    my($class, $typename, $c) = @_;
    my $type = $types{$typename} // die "Unknown typename '$typename'";
    eval qq{use $type; 1} or die $@;
    my $self = bless {
        c => $c,
        n => $c ? $c->n : undef,
        f => $c ? $c->f : undef,
        min => $c ? $c->min : undef,
    }, $type;
    $self->init;
    return $self;
}

sub typename {
    my($self) = @_;
    my $class = ref($self) // $self;
    return $typename{$class} // die "Unknown type '$class'";
}

sub invoke {
    my($self, $which, $named, $args, $log) = @_;
    my $prog = $PROG{$which} // die "Unknown prog to invoke '$which'";
    my $typename = $self->typename;
    my $pid = fork();
    unless ($pid) {
        open STDOUT, '>', $log
                or die "503 Can't open $log for writing: $!";
        exec($prog $named, "-y$typename", @$args)
                or die "505 Could not exec $prog";
    }
    return $pid;
}

sub invoked {
    my($self, $which, $named, $args, $log) = @_;
    my $pid = $self->invoke($which, $named, $args, $log);
    waitpid($pid, 0);
    open my $f, '<', $log
            or die "504 Can't open $log for reading: $!";
    my %line;
    while (<$f>) {
        chomp;
        my($rc) = /^(\d{3}) /
                or die "502 Can't parse log line '$_'";
        push @{ $line{$rc} }, $_;
    }
    return \%line;
}

sub logpath {
    my($self) = @_;
    return sprintf "%s/%s", $LOGS, $self->typename;
}

sub check_exceptions {
    return;
}

sub check_fixed {
    return;
}

for my $method (qw{
    name dbname func_value func_name func func_target
    apply_m to_testf test_target
    gprio maxg
}) {
    no strict 'refs';
    *$method = sub { die "$_[0] must provide implementation of sub $method" };
}

sub dbuser { 'hv' }
sub dbpass { 'hv' }

sub c { $_[0]->{c} }
sub debug { $_[0]->{debug} }

sub to_test {
    my($self) = @_;
    return $self->to_testf($self->f);
}

for my $method (qw{ n f min }) {
    no strict 'refs';
    *$method = sub { $_[0]->{$method} };
}

sub fprio {
    my($self, $n, $k, $expect) = @_;
    my $prio = $self->gprio($n);
    if ($expect) {
        $prio += List::Util::min(0, -log($expect) / log(2));
    }
    return $prio;
}

1;
