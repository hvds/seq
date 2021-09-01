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
#   that the k'th target is xy^z, so far a given y we have:
#     d = (xy^z - n) / k    (tauseq)
#     d = xy^z - kn         (addseq)
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
        'parent' => $c,
    );
    my $type = $c->type;
    $self->set_type($type);

    ($z & 1) == 0 or die "Constraint::Power->new: \$z must be even (not $z)";
    @$self{qw{ pow_k pow_x pow_z pow_g }} = ($k, $x, $z, gcd($k, $x));
    $self->{'min'} = $type->dtoceily($self, $c->min);
    $self->{'max'} = $type->dtoceily($self, $c->max);

    $_ = $self->convert_mod_override($_) for grep ref($_), @$opt_mpow;
    $self->mod_override($_) for @$opt_mpow;

    # Copy over all the constraints: we need only look at "uniquely disallowed"
    # since the rest will be duplicates.
    # FIXME: the mod_ytod calls can be really slow; consider moving the loop
    # body to a $type method to allow inlining.
    my @v = ('', map $c->c($_)->[1], 1 .. $self->check);
    for my $mod (1 .. $self->check()) {
        next if $c->c($mod)->[8] == 0;
        for my $v (0 .. ($mod - 1) / 2) {
            my($subv, $submod) = $type->mod_ytod($self, $v, $mod);
            if (!$submod || vec($v[$submod], $subv, 1)) {
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

sub pow_k { shift->{pow_k} }
sub pow_x { shift->{pow_x} }
sub pow_z { shift->{pow_z} }
sub pow_g { shift->{pow_g} }

sub cur {
    my $self = shift;
    return $self->type->ytod($self, $self->{cur});
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
    return $self->type->ytod($self, $cur);
}

sub convert_mod_override {
    my($self, $array) = @_;
    my($override, $which) = @$array;
    ++$which;
    my($mod, $op, $val) = ($override =~ m{ ^ (\d+) (=) (\d+) \z }x)
            or die "Invalid power mod override '$override'";
    my($found, $coval) = (0, undef);
    my $type = $self->type;
    VAL: for (0 .. $mod - 1) {
        my($subval, $submod) = $type->mod_ytod($self, $_, $mod);
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
    unless ($found == $which) {
        printf <<DIE, $found, $which, $val, $mod;
519 Found only %s of %s values matching %s (mod %s)
DIE
        exit 1;
    }
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
    my $type = $self->type;
    # we require pow_z even
    for my $w (0 .. ($p - 1) / 2) {
        my $subw = $type->ytod($self, $w) % $p;
        next if $subw != $v;
        $self->power_suppress($p, $w, $min, $depend);
        $self->power_suppress($p, $p - $w, $min, $depend) if $w;
    }
}

1;
