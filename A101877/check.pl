#!/opt/maths/bin/perl -w
use strict;
use Math::BigRat;

my $total = Math::BigRat->new(0);
for (grep length, map split(/\s+/), @ARGV) {
	$total += Math::BigRat->new("1 / $_");
}
print "$total\n";
