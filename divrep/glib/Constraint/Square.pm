package Constraint::Square;
use strict;
our @ISA = qw/ Constraint /;

use ModFunc qw/ quadvec mod_combine /;
use warnings;
no warnings qw/ recursion /;

#
# Constraint::Square->new($c, $k, $ty2)
# - specialization of Constraint, to find values matching the constraints
#   specified for the Constraint object $c with the additional requirement
#   that the value is of the form y^2 - $c->n() / $k for some y, with
#   tau(y^2) = ty2.
#
sub new {
	my($class, $c, $k, $ty2) = @_;
	my $self = $class->SUPER::new(
		'n' => $c->n(),
		'f' => $c->f(),
		'tell_count' => $c->tell_count(),
		't0' => $c->t0(),
		'min' => $c->min(),
		'max' => $c->max(),
		'check' => $c->check(),
		'tau' => $c->tau(),
	);
	my $diff = $self->{'sq_diff'} = $self->n() / $k;	# must divide
	$self->{'sq_tau'} = $ty2;

	$self->{'min'} = $self->_convert($c->min() - 1) + 1;
	$self->{'max'} = $self->_convert($c->max() - 1) + 1;

	# Copy over all the constraints
	for my $mod (1 .. $self->check()) {
		my $vec = $c->c($mod)->[2];
		for my $v (0 .. ($mod - 1) / 2) {
			my $subv = ($v * $v - $diff) % $mod;
			if (vec($vec, $subv, 1)) {
				$self->square_suppress($mod, $v);
				$self->square_suppress($mod, $mod - $v) if $v;
			}
		}
	}

	# Need to copy the pend list, but convert the trigger values. How to
	# do that without calculating a bunch of square roots?
#	$self->{'pend'} = [ map [
#		$self->_convert_trigger($_->[0]), @$_[1..$#$_]
#	], @{ $c->{'pend'} } ];
warn sprintf "Throwing away %s pending values triggering in range [%s, %s]\n", 0+@$_, $_->[0][0], $_->[$#$_][0] for grep @$_, $c->{'pend'};

	return $self;
}

sub _convert {
	my($self, $val) = @_;
	return +($val + $self->{'sq_diff'})->bsqrt;
}

sub cur {
	my $self = shift;
	$self->{'cur'} * $self->{'cur'} - $self->{'sq_diff'};
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
	$cur = Constraint::cnext($cur, $mult, $sc, $rebuild);
	$self->{'cur'} = $cur;
	$self->{'tests'} += $t;
	$self->{'skipped'} += $u;
	$self->{'kept'} += 1;
	return undef if $cur > $self->{'max'};
	return $cur * $cur - $self->{'sq_diff'};
}

sub square_suppress {
	my($self, $p, $v, $min, $depend) = @_;
	no warnings qw/ redefine /;
	local *suppress = sub { shift->SUPER::suppress(@_) };
	$self->SUPER::suppress($p, $v, $min, $depend);
}

sub suppress {
	my($self, $p, $v, $min, $depend) = @_;
	my $diff = $self->{'sq_diff'};
	for my $w (0 .. ($p - 1) / 2) {
		my $subw = ($w * $w - $diff) % $p;
		next if $subw != $v;
		$self->square_suppress($p, $w, $min, $depend);
		$self->square_suppress($p, $p - $w, $min, $depend) if $w;
	}
}

1;
