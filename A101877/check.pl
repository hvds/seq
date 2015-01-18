#!/opt/maths/bin/perl -w
use strict;
use Math::BigRat;

my($n, $k) = @ARGV;
open(my $pipe, "./C_int9 $n $k | tail -1 |") or die "./C_int9: $!";
my $set = <$pipe>;
close $pipe;

if ($set =~ /no solution/) {
	print "Error: a($n) = $k found no solution\n";
} elsif ($set =~ /^Possible/) {
	print "Dubious: a($n) = $k found $set";
} else {
	my $total = Math::BigRat->new(0);
	for (grep length, split /\s+/, $set) {
		$total += Math::BigRat->new("1 / $_");
	}
	if ($total == $n) {
		print "Verified: a($n) = $k is a valid candidate\n";
	} else {
		print "Error: a($n) = $k gave set with sum $total\n";
	}
}
