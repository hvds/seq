#!/usr/bin/perl -w
use strict;
use Math::GMP;

=head kernel of 2

For a number n, let f(n) be the set of numbers gotten by splitting n^2 at the 0 
digits.  For example

29648^2 = 879003904

so f(29648) = { 4, 39, 879 }

Let S be the smallest set of numbers containing 2 and fixed by f.  What is the 
largest element of S?
--------------------------------
- David Wilson 

A113917(n) is the largest element; A113918(n) is the number of elements in
the fixed set.

=cut

$| = 1;
my $base = shift(@ARGV) || 10;
my @pend = ('2');
my %seen = ('2' => 1);
my $max = 2;

while (@pend) {
    my $bn = shift @pend;
    my $n = Math::GMP::new_from_scalar_with_base($bn, $base);
    $max = $n if $max < $n;
    my $next = ($n * $n)->get_str_gmp($base);
# print "$bn -> $n -> $next\n";
    for (split /0+/, $next) {
        push @pend, $_ unless $seen{$_}++;
    }
}
printf "base %s: max %s, count %s\n", $base, $max, scalar keys %seen;
exit 0;

__END__
Extend to base b, then a(2) = 2, a(3) = 5575, a(4) = 2
1-11 13-21 31 34
121 141 221 224 243 311 432 441 443-444
1323 2131 2211 2221 2421 3342 4131 4341
11143 23424 31431 43341
112341 221241 414224
13132241
2121234311

2422434 24112114 33331211 42331331 432414341 4224424331 13243332331
