package Quadratic;
use strict;
use Math::Pari qw/ gcd /;

my %q;

sub quadres {
	my($class, $n, $mod) = @_;
	return 0 if gcd($n, $mod) > 1;
	my $v = $class->quadvec($mod);
	return vec($v, ($n % $mod), 1) ? 1 : -1;
}

sub is_residue {
	my($class, $n, $mod) = @_;
	my $v = $class->quadvec($mod);
	return vec($v, ($n % $mod), 1);
}

sub quadvec {
	my($class, $mod) = @_;
	${ $q{$mod} ||= do {
		my $v = "";
		vec($v, (($_ * $_) % $mod), 1) = 1 for 1 .. $mod - 1;
		\$v;
	} };
}
	
1;
