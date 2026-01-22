package Shape;
use strict;
use warnings;

use Fcntl qw{ SEEK_SET };

sub new {
    my($class, $data) = @_;
    return bless $data, $class;
}

sub new_empty {
    my($class) = @_;
    my $self = $class->new({
        xmin => 0,
        xmax => -1,
        ymin => 0,
        ymax => -1,
        point => [],
        disallowed => [],
        neighbour => [('') x 2],
        count => 0,
        off => 0,
    });
    # set a seed neighbour
    $self->setv($self->neighbour, 1, 1);
    return $self;
}

sub parse_recover {
    my($class, $line, $i, $fetcher) = @_;
    chomp $line;
    my($off, $time) = $line =~ m{^A (\d+) ([\d.]+)$}
            or die "parse error '$line'";
warn "recover $line\n";
    my $struct = {};
    my $parent_line = $fetcher->($i, $off, $struct);
    my $parent = $class->parse($parent_line, $off);
warn "recover parent @{[ $parent->disp ]}\n";
    @$struct{qw{ shape time }} = ($parent, $time);
    return $struct;
}

sub parse {
    my($class, $line, $off) = @_;
    chomp $line;
    my($orig, $disallowed) = split ';', $line, 2;
    die "parse error ($line)" unless $orig =~ s{^(\d+):(\d+):}{};
    my($xspan, $yspan) = ($1, $2);
    my @pt_rows = split m{/}, $orig;
    my @dis_rows = split m{/}, $disallowed;
    # build pre-expanded
    my $self = $class->new({
        xmin => 2,
        xmax => $xspan + 1,
        ymin => 2,
        ymax => $yspan + 1,
        count => 0,
        point => [],
        disallowed => [],
        neighbour => [],
        off => $off,
    });
    my($pt, $dis, $nb) = @$self{qw{ point disallowed neighbour }};
    for my $i (0 .. $#pt_rows) {
        my $row = $pt_rows[$i];
        for my $j (0 .. length($row)) {
            next unless substr($row, $j, 1) eq '*';
            my($x, $y) = ($i + 2, $j + 2);
            $self->setv($pt, $x, $y);
            $self->unsetv($nb, $x, $y);
            for ([-1, 0], [1, 0], [0, -1], [0, 1]) {
                my($xn, $yn) = ($x + $_->[0], $y + $_->[1]);
                $self->setv($nb, $xn, $yn)
                        unless $self->hasv($pt, $xn, $yn);
            }
            ++$self->{count};
        }
    }
    for my $i (0 .. $#dis_rows) {
        my $row = $dis_rows[$i];
        for my $j (0 .. length($row)) {
            next unless substr($row, $j, 1) eq '*';
            my($x, $y) = ($i + 1, $j + 1);
            $self->setv($dis, $x, $y);
            $self->unsetv($nb, $x, $y);
        }
    }
    return $self;
}

sub count { shift->{count} }

sub point { shift->{point} }

sub neighbour { shift->{neighbour} }

sub disallowed { shift->{disallowed} }

sub clone {
    my($self) = @_;
    return ref($self)->new({
        %$self,
        point => [ @{ $self->point } ],
        disallowed => [ @{ $self->disallowed } ],
        neighbour => [ @{ $self->neighbour } ],
    });
}

sub neighbours {
    my($self) = @_;
    return $self->allv($self->neighbour);
}

sub disallow {
    my($self, $x, $y) = @_;
    $self->setv($self->disallowed, $x, $y);
    return;
}

sub writable {
    my($self) = @_;
    return sprintf "%s;%s\n", $self->disp, $self->rdd;
}

sub aligner {
    my($self, $t0) = @_;
    return sprintf "A %s %s\n", $self->{off}, sprintf '%.2f', (times)[0] - $t0;
}

sub end_writable {
    my($class, $t0) = @_;
    return sprintf "E %s\n", sprintf '%.2f', (times)[0] - $t0;
}

sub canonical {
    my($self) = @_;
    return '' if $self->count == 0;
    my($xmin, $xmax, $ymin, $ymax) = @$self{qw{ xmin xmax ymin ymax }};
    my($xspan, $yspan) = ($xmax - $xmin, $ymax - $ymin);
    my($s, $t);
    my $pt = $self->point;
    if ($xspan <= $yspan) {
        $s = $self->disp($pt);
        my $tx = $self->invx($pt, $xmin, $xmax);

        $t = $self->disp($tx);
        $s = $t if $s gt $t;

        my $tr = $self->revx($pt, $ymin, $ymax);
        $t = $self->disp($tr);
        $s = $t if $s gt $t;

        my $trx = $self->invx($tr, $xmin, $xmax);
        $t = $self->disp($trx);
        $s = $t if $s gt $t;
    }

    if ($yspan <= $xspan) {
        local @$self{qw{ xmin xmax ymin ymax }} = ($ymin, $ymax, $xmin, $xmax);
        my $ty = $self->trans($pt, $xmin, $xmax, $ymin, $ymax);
        $t = $self->disp($ty);
        $s = $t if !defined($s) || $s gt $t;

        my $tyx = $self->invx($ty, $ymin, $ymax);
        $t = $self->disp($tyx);
        $s = $t if $s gt $t;

        my $tyr = $self->revx($ty, $xmin, $xmax);
        $t = $self->disp($tyr);
        $s = $t if $s gt $t;

        my $tyrx = $self->invx($tyr, $ymin, $ymax);
        $t = $self->disp($tyrx);
        $s = $t if $s gt $t;
    }
    return $s;
}

sub add_new {
    my($self, $x, $y) = @_;
    $self = $self->clone;
    if ($self->{count} == 0) {
        @$self{qw{ xmin xmax ymin ymax }} = ($x, $x, $y, $y);
    } else {
        $self->{xmin} = $x if $self->{xmin} > $x;
        $self->{xmax} = $x if $self->{xmax} < $x;
        $self->{ymin} = $y if $self->{ymin} > $y;
        $self->{ymax} = $y if $self->{ymax} < $y;
    }
    $self->setv($self->point, $x, $y);
    my $disallowed = $self->disallowed;
    $self->setv($disallowed, $x, $y);
    my $nb = $self->neighbour;
    $self->unsetv($nb, $x, $y);
    for ([-1, 0], [1, 0], [0, -1], [0, 1]) {
        my($xn, $yn) = ($x + $_->[0], $y + $_->[1]);
        $self->setv($nb, $xn, $yn)
                unless $self->hasv($disallowed, $xn, $yn);
    }
    ++$self->{count};
    return $self;
}

sub expand {
    my($self) = @_;
    return if $self->count == 0;
    while ($self->{xmin} < 2) {
        unshift @{ $self->point }, '';
        unshift @{ $self->disallowed }, '';
        unshift @{ $self->neighbour }, '';
        ++$self->{xmin};
        ++$self->{xmax};
    }
    if ($self->{ymin} < 2) {
        for ($self->point, $self->disallowed, $self->neighbour) {
            @$_ = map "\x{0}" . $_, @$_;
        }
        $self->{ymin} += 8;
        $self->{ymax} += 8;
    }
    return;
}

sub colinear {
    my($self, $x, $y) = @_;
    my %seen;
    my $high = 0;
    for ($self->allv($self->point)) {
        my($x2, $y2) = @$_;
        my $v = ($y == $y2) ? "oo" : Math::BigRat->new($x - $x2, $y - $y2);
        ++$seen{$v};
        $high = $seen{$v} if $high < $seen{$v};
    }
# $x -= $self->{xmin}; $y -= $self->{ymin}; warn "col $high for [$x,$y] in @{[ $self->disp ]}\n";
    return $high + 1;
}

sub Dump {
    my($self) = @_;
    my @r = @$self{qw{ xmin xmax ymin ymax }};
    my @s = ($self->rdp, $self->rdd, $self->rdn);
    return sprintf <<EOS, $self->{count}, @r, @s;
count %s; xmin %s; xmax %s; ymin %s; ymax %s
pt: %s
ds: %s
nb: %s
EOS
}

sub rdisp {
    my($self, $data, $xmin, $xmax, $ymin, $ymax) = @_;
    return join '/', map {
        my $x = $_;
        join '', map {
            $self->hasv($data, $x, $_) ? '*' : '.'
        } $ymin .. $ymax;
    } $xmin .. $xmax;
}

sub rdp {
    my($self) = @_;
    return $self->rdisp($self->point, @$self{qw{ xmin xmax ymin ymax }});
}

sub rdd {
    my($self) = @_;
    return $self->rdisp($self->disallowed,
        $self->{xmin} - 1, $self->{xmax} + 1,
        $self->{ymin} - 1, $self->{ymax} + 1
    );
}

sub rdn {
    my($self) = @_;
    return $self->rdisp($self->neighbour,
        $self->{xmin} - 1, $self->{xmax} + 1,
        $self->{ymin} - 1, $self->{ymax} + 1
    );
}

sub disp {
    my($self, $data) = @_;
    $data //= $self->{point};
    my($xmin, $xmax, $ymin, $ymax) = @$self{qw{ xmin xmax ymin ymax }};
    my $xspan = $xmax - $xmin + 1;
    my $yspan = $ymax - $ymin + 1;
    return sprintf '%d:%d:%s',
            $xmax - $xmin + 1, $ymax - $ymin + 1,
            $self->rdisp($data, $xmin, $xmax, $ymin, $ymax);
}

sub hasv {
    my($self, $v, $x, $y) = @_;
    return 0 if $#$v < $x;
    return vec($v->[$x], $y, 1) ? 1 : 0;
}

sub setv {
    my($self, $v, $x, $y) = @_;
    push @$v, '' while @$v < $x;
    vec($v->[$x], $y, 1) = 1;
    return;
}

sub unsetv {
    my($self, $v, $x, $y) = @_;
    return if $#$v < $x;
    vec($v->[$x], $y, 1) = 0;
    if ($#$v == $x) {
        pop @$v while @$v && $v->[-1] !~ /[^\x{0}]/;
    }
    return;
}

sub allv {
    my($self, $v) = @_;
    return map {
        my $i = $_;
        my $vi = $v->[$i];
        my $len = 8 * length($vi);
        map [ $i, $_ ], grep vec($vi, $_, 1), 0 .. $len - 1;
    } 0 .. $#$v;
}

sub invx {
    my($self, $v, $xmin, $xmax) = @_;
    return [ ('') x $xmin, reverse @$v[ $xmin .. $xmax ] ];
}

sub revx {
    my($self, $v, $ymin, $ymax) = @_;
    my $w = [ ('') x @$v ];
    for my $i (0 .. $#$v) {
        my $src = $v->[$i];
        next if $src =~ /^\x{0}*\z/;
        my $dest = $w->[$i];
        for my $j (0 .. $ymax - $ymin) {
            vec($dest, $ymin + $j, 1) = vec($src, $ymax - $j, 1);
        }
        $w->[$i] = $dest;
    }
    return $w;
}

sub trans {
    my($self, $v, $xmin, $xmax, $ymin, $ymax) = @_;
    my $w = [ ('') x ($ymax + 1) ];
    for my $i ($xmin .. $xmax) {
        for my $j ($ymin .. $ymax) {
            vec($w->[$j], $i, 1) = vec($v->[$i] //= '', $j, 1);
        }
    }
    return $w;
}

1;
