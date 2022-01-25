package Type;
use strict;
use warnings;
use List::Util ();
use Scalar::Util qw{ weaken };

my %types = (
    t => 'Type::TauSeq',
    a => 'Type::AddSeq',
    s => 'Type::Semip',
    o => 'Type::OneSeq',
    d => 'Type::AscDPrime',
    n => 'Type::AscNPrime',
);
my %typename = reverse %types;

my $LOGS = './logs';
my %PROG = (
    gtauseq => './gtauseq',
    bisect_fp => './bisect-fp',
    bisect_g => './bisect-g',
    find_shard => './find-shard',
);

my %owner = (
    harness => 0,
    upperlim => 1,
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
        debug => $c ? $c->debug : 0,
    }, $type;
    weaken($self->{c}) if $c;
    $self->init;
    return $self;
}

sub bind {
    my($self, $n) = @_;
    $self->{n} = $n;
    $self->init;
    return;
}

sub typename {
    my($self) = @_;
    my $class = ref($self) // $self;
    return $typename{$class} // die "Unknown type '$class'";
}

sub _owner {
    my($proto, $name) = @_;
    return $owner{$name} // die "Unknown owner '$name'";
}

sub bind_owner {
    my($self, $name) = @_;
    $self->{owner_id} = $self->_owner($name);
    return;
}

sub owner {
    my($self) = @_;
    return $self->{owner_id} // die "No owner bound";
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

sub check_mult {
    return;
}

for my $method (qw{
    name dbname func_value func_name func func_target
    apply_m to_testf test_target
    gprio ming maxg smallest
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
        if ($expect < 0) {
            use Carp; Carp::confess("$expect");
        }
        $prio += List::Util::min(0, -log($expect) / log(2));
    }
    return $prio;
}

# Most types have constant func_target, but not all
sub func_matches {
    my($self, $k, $result) = @_;
    my $target = $self->func_target($k) // die "No target set in $self";
    return $target == $result;
}

1;
