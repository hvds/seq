package H::Partition;
use strict;
use warnings;

=head1 NAME

H::Partition - Calculates the partitions of an integer

=head1 DESCRIPTION

Calculates the partitions of an integer I<n>, ie the collections of integers
that sum to I<n>, subject to various constraints, and in various styles.

Primary mechanisms cache the list of partitions for every I<n> up to the
greatest asked for; P(60) takes about 20s on my machine to return the list
(length 966,467), and consumes about 330MB.

Throughout, "array" style means each partition is represented as an arrayref
of integers; "string" style means the same partition is packed into a string
in which C<ord()> of each character is the corresponding integer.

=head1 METHODS

=over 4

=item strings ( n )

Returns the partitions of I<n> as an arrayref of strings.

This is the actual arrayref from the cache; the contents should not be
modified.

=item strings_max ( n, max )

Returns the partitions of I<n> as a tuple C<(arrayref, length)>; the arrayref
is the same returned by C<strings()>, but within that arrayref the strings
at indexes C<0 .. length - 1> represent those partitions that have no
component exceeding I<max>.

This is the actual arrayref from the cache; the contents should not be
modified.

=back

=head1 TO DO

Add non-caching iterator support; add array-style support; add other
constraints; find a better naming scheme.

=cut

{
	# initialize the caches for P(0)
	my @strings;             BEGIN { @strings = (['']) }
	my @strings_max_offsets; BEGIN { @strings_max_offsets = ([1]) }

	sub strings {
		my($class, $n) = @_;
		return $strings[$n] ||= do {
			my @partition;
			my @offset = (0);
			for my $start (1 .. $n) {
				my $rest = $n - $start;
				my $tail = $class->strings($rest);
				my $max = $start;
				$max = $rest if $max > $rest;
				my $offset = $strings_max_offsets[$rest][$max] || 0;
				my $start_char = chr($start);
				push @partition, map "$start_char$_", @$tail[0 .. $offset - 1];
				push @offset, scalar @partition;
			}
			$strings_max_offsets[$n] = \@offset;
			\@partition;
		};
	}

	sub strings_max {
		my($class, $n, $max) = @_;
		my $partition = $class->strings($n);
		my $offset = $strings_max_offsets[$max];
		return +($partition, $offset);
	}
}

1;
