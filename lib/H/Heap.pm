package H::Heap;
use strict;
use warnings;

my $index = 0;
my %names;

=head1 NAME

H::Heap - Simple configurable binary heap

=head1 DESCRIPTION

Provides a "minimum" heap, using the supplied comparator.

=head1 METHODS

=over 4

=item new ( comparator )

Returns a new empty heap object for the supplied comparator.

The comparator should be an eval-able string that compares two items
in C<$a> and C<$b>, and returns values corresponding to those of C<< <=> >>.

=item size ( )

Returns the number of items in the heap.

=item insert ( object )

Inserts the specified I<object> into the heap. Returns nothing.

=item fetch ( )

Removes and returns the least item from the heap. If the heap is empty,
returns C<undef>.

=item peek ( )

Returns the least item from the heap, without removing it. If the heap
is empty, returns C<undef>.

=item peek2 ( )

Returns the second least item from the heap. If the heap does not have
at least 2 items, returns C<undef>.

=back

=cut

sub new {
	my($class, $cmp_s) = @_;
	my $subclass = $names{$cmp_s} ||= do {
		my $c = 'H::Heap::Heap_' . $index++;
		eval sprintf '
			package %s;
			our @ISA = ($class);
			sub cmp_cv {
				my($self, $a, $b) = @_;
				%s;
			}
		', $c, $cmp_s;
		$c;
	};
	bless [], $subclass;
}

sub size {
    my($self) = @_;
    return scalar @$self;
}

sub insert {
	my($self, $new) = @_;
	my $node = @$self;
	push @$self, $new;
	while ($node) {
		my $parent = ($node - 1) >> 1;
		last if $self->cmp_cv(@$self[$parent, $node]) <= 0;
		@$self[$parent, $node] = @$self[$node, $parent];
		$node = $parent;
	}
	return;
}

sub peek {
    my($self) = @_;
    return $self->[0];
}

sub peek2 {
    my($self) = @_;
    return $self->[1] if @$self < 3;
    return $self->cmp_cv(@$self[1, 2]) < 0 ? $self->[1] : $self->[2];
}

sub fetch {
	my($self) = @_;
	return undef unless @$self;
	my $value = $self->[0];
	my $switch = pop @$self;
	$self->[0] = $switch if @$self;
	my $node = 0;
	my $size = @$self;
	while (1) {
		my $child = ($node << 1) + 1;
		last if $child >= $size;
		if ($child + 1 < $size
			&& $self->cmp_cv(@$self[$child, $child + 1]) > 0
		) {
			++$child;
		}
		last if $self->cmp_cv(@$self[$child, $node]) >= 0;
		@$self[$node, $child] = @$self[$child, $node];
		$node = $child;
	}
	return $value;
}

1;
