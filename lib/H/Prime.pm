package H::Prime;
use strict;
use warnings;

=head1 NAME

H::Prime - Useful functions relating to primes.

=head1 DESCRIPTION

Calculate basic functions relating to prime numbers and prime factorization.

=head1 FUNCTIONS

=over 4

=item prime_i ( i )

Returns the I<i>th prime number, such that prime_i(0) returns 2.

=back

=head1 TODO

This should detect the presence of C<Math::BigInt> or a specific bigint
library, and adapt.

For the small integer case, tighten up with either overflow detection
or explicit limits.

=cut

{
	my @p = (2, 3);

	sub _nextPrime {
		my $p = $p[$#p];
		DIV: while (1) {
			$p += 2;
			my $pc = 0;
			while (1) {
				my $d = $p[$pc++];
				last if $d * $d > $p;
				next DIV unless $p % $d;
			}
			$p[@p] = $p;
			return $p;
		}
	}

	sub prime_i {
		my($i) = @_;
		_nextPrime() while @p <= $i;
		return $p[$i];
	}
}

1;
