package Constraint::Power;
use strict;
our @ISA = qw/ Constraint /;

use ModFunc qw/ quadvec mod_combine gcd /;
use warnings;
no warnings qw/ recursion /;

#
# Constraint::Power->new($c, $k, $x, $z, $opt_mpow)
# - specialization of Constraint, to find values matching the constraints
#   specified for the Constraint object $c with the additional requirement
#   that n + kd = xy^z, so for a given y we have d = (xy^z - n) / k.
#
sub new {
    my($class, $c, $k, $x, $z, $opt_mpow) = @_;
    my $self = $class->SUPER::new(
        'n' => $c->n(),
        'f' => $c->f(),
        'tell_count' => $c->tell_count(),
        't0' => $c->t0(),
        'min' => $c->min(),
        'max' => $c->max(),
        'check' => $c->check(),
        'tau' => $c->tau(),
        'parent' => $c,
    );

    ($z & 1) == 0 or die "Constraint::Power->new: \$z must be even (not $z)";
    @$self{qw{ pow_k pow_x pow_z }} = ($k, $x, $z);
    $self->{'min'} = $self->_dtoceily($c->min());
    $self->{'max'} = $self->_dtoceily($c->max());

    $_ = $self->convert_mod_override($_) for grep ref($_), @$opt_mpow;
    $self->mod_override($_) for @$opt_mpow;

    # Copy over all the constraints
    my %v;
    for my $mod (1 .. $self->check()) {
        for my $v (0 .. ($mod - 1) / 2) {
            my($subv, $submod) = $self->_mod_ytod($v, $mod);
            my $vec = $submod ? ($v{$submod} //= $c->c($submod)->[2]) : undef;
            if (!$submod || vec($vec, $subv, 1)) {
                $self->power_suppress($mod, $v);
                $self->power_suppress($mod, $mod - $v) if $v;
            }
        }
    }

    # Need to copy the pend list, but convert the trigger values. How to
    # do that without calculating a bunch of square roots?
#   $self->{'pend'} = [ map [
#       $self->_convert_trigger($_->[0]), @$_[1..$#$_]
#   ], @{ $c->{'pend'} } ];
    printf(
        "304 Throwing away %s pending values triggering in range [%s, %s]\n",
        0 + @$_, $_->[0][0], $_->[$#$_][0],
    ) for grep @$_, $c->{'pend'};

    return $self;
}

#
# Calculate floor(y) given d: floor(y) = floor(((n + kd) / x) ^ (1/z))
#
sub _dtoy {
    my($self, $val) = @_;
    my $base = $self->{n} + $self->{pow_k} * $val;
    return +($base / $self->{pow_x})->broot($self->{pow_z});
}

sub _dtoceily {
    my($self, $val) = @_;
    return 1 + $self->_dtoy($val - 1);
}

#
# Calculate d given y: d = (xy^z - n) / k
# We expect k should always divide the expression exactly, or raise an error;
# however in rare cases that's ok, these are marked by setting $::LOOSE_CUR.
#
sub _ytod {
    my($self, $val) = @_;
    my $base = $self->{pow_x} * $val ** $self->{pow_z} - $self->{n};
    my($div, $rem) = $base->bdiv($self->{pow_k});
    if ($rem != 0) {
        use Carp;
        confess sprintf(
            "_ytod(k = %s, x = %s, z = %s => %s) not divisible by k",
            @$self{qw{ pow_k pow_x pow_z }}, $val,
        ) unless $::LOOSE_CUR;
        # if loose is good enough, return the floor but mark it
        $div .= ':LOOSE';
    }
    return $div;
}

#
# Given y == y_m (mod m) and d = (xy^z - n) / k, return (d_s, s) as the
# value and modulus of the corresponding constraint on d, d == d_s (mod s).
# If no valid d is possible, returns s == 0.
#
sub _mod_ytod {
    my($self, $val, $mod) = @_;
    my($n, $k, $x, $z) = @$self{qw{ n pow_k pow_x pow_z }};
    my $base = $x * $val ** $z - $n;
    my $g = gcd($k, $x);
    my $gbase = $base / $g;
    my $gk = $k / $g;
    my $g2 = gcd($gk, $mod);
    if ($gk == 1) {
        return ($gbase % $mod, $mod);
    } elsif ($g2 == 1) {
        my $inv = $gk->bmodinv($mod);
        return (($gbase * $inv) % $mod, $mod);
    } elsif (($gbase % $g2) == 0) {
        my $gmod = $mod / $g2;
        my $g2base = $gbase / $g2;
        my $g2k = $gk / $g2;
        my $inv = $g2k->bmodinv($gmod);
        return (($g2base * $inv) % $gmod, $gmod);
    } else {
        return (1, 0);
    }
}

sub cur {
    my $self = shift;
    return $self->_ytod($self->{cur});
}

sub next {
    my($self) = @_;
    my $mult = $self->{'mult'};
    my $cur = $self->{'cur'} + $mult;
    my $rebuild = 0;
    do {
        my $trigger = ($self->{'pend'}[0] || [0])->[0];
        if ($trigger && $trigger < $cur) {
            $cur = $self->catchup($cur);
            $mult = $self->{'mult'};
            $rebuild = 1;
        }
    };
    my $sc = $self->{'sc'};
    my($t, $u) = (0, 0);
    $cur = Constraint::cnext($cur, $mult, $sc, $rebuild, $self->max);
    $self->{'cur'} = $cur;
    $self->{'tests'} += $t;
    $self->{'skipped'} += $u;
    $self->{'kept'} += 1;
    return undef if $cur > $self->{'max'};
    return $self->_ytod($cur);
}

sub convert_mod_override {
    my($self, $array) = @_;
    my($override, $which) = @$array;
    ++$which;
    my($mod, $op, $val) = ($override =~ m{ ^ (\d+) (=) (\d+) \z }x)
            or die "Invalid power mod override '$override'";
    my($found, $coval) = (0, undef);
    VAL: for (0 .. $mod - 1) {
        my($subval, $submod) = $self->_mod_ytod($_, $mod);
        next VAL unless $submod;
        if ($mod != $submod) {
            my $cosubval = ($val % $submod);
            next VAL unless $cosubval == $subval;
            printf <<WARN, $subval, $submod, $cosubval, $mod;
318 Ambiguous mod fix: y==%s (mod %s) yields d==%s (mod %s)
WARN
            $subval = $val; # treat as valid
        }
        next VAL unless $subval == $val;
        $coval = $_;
        last VAL if ++$found == $which;
    }
    die sprintf <<DIE, $found, $which, $val, $mod unless $found == $which;
519 Found only %s of %s values matching %s (mod %s)
DIE
    return "$mod=$coval";
}

sub mod_override {
    my($self, $override) = @_;
    my($mod, $op, $val) = ($override =~ m{ ^ (\d+) ([=!]) (\d+) \z }x)
            or die "Invalid power mod override '$override'";
    if ($op eq '=') {
        $self->power_require($mod, $val);
    } else {
        $self->power_suppress($mod, $val);
    }
}

sub power_require {
    my($self, $p, $v, $min) = @_;
    no warnings qw/ redefine /;
    local *suppress = sub { shift->SUPER::suppress(@_) };
    $self->SUPER::require($p, $v, $min);
}

sub power_suppress {
    my($self, $p, $v, $min, $depend) = @_;
    no warnings qw/ redefine /;
    local *suppress = sub { shift->SUPER::suppress(@_) };
    $self->SUPER::suppress($p, $v, $min, $depend);
}

sub suppress {
    my($self, $p, $v, $min, $depend) = @_;
    # we require pow_z even
    for my $w (0 .. ($p - 1) / 2) {
        my $subw = $self->_ytod($w) % $p;
        next if $subw != $v;
        $self->power_suppress($p, $w, $min, $depend);
        $self->power_suppress($p, $p - $w, $min, $depend) if $w;
    }
}

1;
