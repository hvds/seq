package Prove::Poly;
use strict;
use warnings;

use Math::Prime::Util qw{ gcd lcm };
use Math::BigRat;
sub MBR { Math::BigRat->new(@_) }
my $zero = MBR(0);

=head1 NAME

Prove::Poly - representation of a multivariate polynomial

=head1 DESCRIPTION

Top-level object is [ [ m_i, c_i ]* ] where c_i is rational and
m_i are as below, representing sum(c_i m_i).

m_i are hashrefs of the form (a_i => x_i) where x_i are integer and
a_i are as below, representing prod(a_i^x_i).

a_i are packed strings of the form [length, n_i] all integers represented
as signed 8-bit characters, representing the linear combination sum(n_i v_i)
where v_i is the i'th variable, and v_0 is integer 1. For convenience,
the n_i are stored in reverse order, and are normalised such that there is
no factor common to all elements and the first element (ie the coefficient
of the highest indexed variable) is positive.

=cut

sub new {
    my($class, $data) = @_;
    return bless $data, $class;
}

# Get the constant value; it is a fatal error if we are not constant.
sub const {
    my($self) = @_;
    die "panic" if @$self > 1;
    my($mi, $ci) = @{ $self->[0] };
    die "panic" if keys %$mi;
    return $ci;
}

# assume normalized, and in signed 8-bit range
sub _packlc {
    my($lc) = @_;
    return pack "c*", $#$lc, reverse @$lc;
}
sub _unpacklc {
    my($sig) = @_;
    return [ reverse unpack "c*", substr $sig, 1 ];
}
sub _packmi {
    my($mi) = @_;
    return join '', map $_ . chr($mi->{$_}), sort keys %$mi;
}

# Initialize an object from a simple arrayref representing a linear
# combination of the variables and 1.
sub from_sc {
    my($class, $sc) = @_;
    return $class->new([ [ { _packlc($sc) => 1 }, 1 ] ]);
}

{
    my %comb;
    # return comb(a, b), assuming it fits in a perl integer
    sub _comb {
        my($a, $b) = @_;
        return $comb{"$a.$b"} //=
            ($b == 0 || $b == $a) ? 1 : (($a * _comb($a - 1, $b - 1)) / $b);
    }
}
        
# Given (m_i, c_i) where m_i has at least 2 a_i containing the variable
# at level, return an arrayref of [m_i, c_i]* representing the same
# expression distributed such that each m_i has at most 1 a_i containing
# that variable.
# Of the affected a_i, we pick the one at the highest power as our base a_0
# (breaking ties with the asciibetically highest packed LC), and split
# the others such that a_i = (k_i a_0 + a_i') with a_i' independent of
# the variable at level.
sub _distrib {
    my($mi, $ci, $level) = @_;
    my %mj = %$mi;
    my @key = sort { $mj{$b} <=> $mj{$a} || $b cmp $a }
            grep ord($_) >= $level, keys %mj;   # ref _unpacklc
    my $pa0 = shift @key;
    my $ua0 = _unpacklc($pa0);
    my $a0h = $ua0->[$level];
    my $q = 1;
    my @comp = ([ { }, 1 ]);
    for (@key) {
        my $ua1 = _unpacklc($_);
        my $a1h = $ua1->[$level];
        my $pow = delete $mj{$_};
        if ($a1h % $a0h) {
            $_ *= $a0h for ($a1h, @$ua1);
            $q *= MBR(1, $a0h);
        }
        my $cj = $a1h / $a0h;
        my($pa2, $ck) = _normal(
            map $ua1->[$_] - $cj * $ua0->[$_], 0 .. $level - 1
        );
use DDP; do { p [ $ua0, $ua1, $a0h, $a1h, $pow, $cj ]; die } if !$ck;
        my @newcomp = map [
            {
                ($_ ? ($pa0 => $_) : ()),
                (($pa2 && ($pow - $_)) ? ($pa2 => $pow - $_) : ()),
            },
            _comb($pow, $_) * ($cj ** $_) * ($ck ** ($pow - $_)),
        ], 0 .. $pow;
        @comp = map {
            my($m1, $c1) = @$_;
            map {
                my($m2, $c2) = @$_;
                my %m3 = %$m1;
                $m3{$_} += $m2->{$_} for keys %$m2;
                [ \%m3, $c1 * $c2 * $q ];
            } @newcomp;
        } @comp;
    }
    return [ map {
        my($mk, $ck) = @$_;
        my %ml = %mj;
        $ml{$_} += $mk->{$_} for keys %$mk;
        [ \%ml, $ci * $ck ];
    } @comp ];
}

# Integrate the supplied m_i with respect to the variable at level,
# returning an (m_i, c_i) pair.
sub _integrate {
    my($mi, $level) = @_;
    # ref _unpacklc
    my($key) = List::Util::first(sub { ord($_) >= $level }, keys %$mi);
    my %copy = %$mi;
    if (defined $key) {
        my $pow = ++$copy{$key};
        my($mul) = unpack('c', substr $key, 1, 1);   # ref _unpacklc
        return +(\%copy, MBR(1, $mul * $pow));
    }
    # cannot already exist
    ++$copy{_packlc([ (0) x $level, 1 ])};
    return +(\%copy, 1);
}

# Return a new object representing the expression obtained by integrating
# ourselves with respect to the variable at level.
sub integrate {
    my($self, $level) = @_;
    my %group;
    for (@$self) {
        my($mi, $ci) = @$_;
        my @topvar = grep ord($_) >= $level, keys %$mi;     # ref _unpacklc
        my $npair = (@topvar <= 1)
                ? [[$mi, $ci]] : _distrib($mi, $ci, $level);
        for (@$npair) {
            my($mj, $cj) = @$_;
            ($group{_packmi($mj)} //= [ $mj, $zero ])->[1] += $cj;
        }
    }
    return ref($self)->new([ map {
        my($mi, $ci) = @$_;
        my($mj, $cj) = _integrate($mi, $level);
        [ $mj, $ci * $cj ];
    } values %group ]);
}

# Return a new object representing the evaluation after integrating with
# respect to the variable at level within the range, expressed as an
# arrayref [low, high], each being in the form of an unpacked LC.
sub inteval {
    my($self, $level, $range) = @_;
    my %group;
    for (@$self) {
        my($mi, $ci) = @$_;
        # there must be exactly one; ref _unpacklc
        my $key = List::Util::first(sub { ord($_) >= $level }, keys %$mi)
                // die "panic";
        my $pow = $mi->{$key};
        for my $rx (0, 1) {
            my $ri = $range->[$rx];
            my $lc = _unpacklc($key);
            my $mult = $lc->[$level];
            my($pa, $cj) = _normal(
                map $lc->[$_] + ($ri->[$_] // 0) * $mult, 0 .. $level - 1
            );
            next unless $cj;
            $cj **= $pow;
            $cj = -$cj if $rx == 0;
            my %mj = do {
                local $mi->{$pa} = ($mi->{$pa} // 0) + $pow if defined $pa;
                delete local $mi->{$key};
                %$mi;
            };
            my $pmj = _packmi(\%mj);
            ($group{$pmj} //= [ \%mj, 0 ])->[1] += $ci * $cj;
        }
    }
    return ref($self)->new([ values %group ]);
}

# return (packed, const); packed may be undef
sub _normal {
    my(@a) = @_;
    pop @a while @a && !$a[-1];
    return +(undef, $a[0] // 0) if @a < 2;
    my $c = 1;
    if ($a[-1] < 0) {
        $_ = -$_ for (@a, $c);
    }
    my $g = gcd(@a);
    if ($g != 1) {
        $c *= $g;
        $_ /= $g for @a;
    }
    return +(_packlc(\@a), $c);
}

# for debugging: return our expression as a string
sub expand {
    my($self) = @_;
    my $s = join '+', map {
        my($m, $c) = @$_;
        (($c == 1) ? '' : $c) . join '', map {
            my($lc, $pow) = (_unpacklc($_), $m->{$_});
            my $lcs = join '+', map {
                my($n, $v) = ($lc->[$_], $_ ? chr(ord('a') - 1 + $_) : '');
                my $sgn = ($n < 0) ? '-' : '';
                (abs($n) == 1 && $v) ? "$sgn$v" : "$n$v";
            } grep $lc->[$_], 0 .. $#$lc;
            $lcs = "($lcs)" if $lcs =~ /\+/;
            ($pow == 1) ? $lcs : "$lcs^$pow";
        } sort keys %$m;
    } @$self;
    return $s =~ s/\+-/-/gr;
}
sub to_string { 
    my($self) = @_;
    return $self->expand;
}

1;
