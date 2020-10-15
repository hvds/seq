package Graph::Hamilton;
use strict;
use warnings;

use List::Util qw{ max };

=head1 findLoop ( Graph $graph, method $method )

Invoked as a class method with a L<Graph::Undirected> object C<$graph> and
a reference to a cycle-finding method (either as a string method name or
as a subref), invokes the method to find a Hamilton cycle as a graph, then
returns the result post-processed into a simple arrayref of vertices in the
order they appear in the cycle, or C<undef> if no cycle was found.

=head1 find ( Graph $graph )

Invoked as a class method with a L<Graph::Undirected> object C<$graph>, finds
and returns a subgraph whose edges form a Hamilton cycle in C<$graph>,
if one exists, else returns C<undef>.

If the whole graph itself is a Hamilton cycle, C<$graph> itself is returned.

Note that this is an NP-complete algorithm, so it may take some time.
It works by deleting edges until only a cycle is left, so is likely to
be most efficient on sparse graphs.

=head1 findBFF ( Graph $graph )

Same as L</find>, but using the probabilistic algorithm described by
Bollobas, Fenner and Frieze in "An algorithm for finding Hamilton paths
and cycles in random graphs" (1985)
L<https://www.math.cmu.edu/~af1p/Texfiles/AFFHCIRG.pdf>. This may return
C<undef> even though a cycle exists, but the paper guarantees that for
B<random> graphs the failures asymptotically approach zero as the number
of nodes in the graph increases.

Note that this is a polynomial-time algorithm, but may use a lot of memory.
It works by extending a chain of links however possible, so is likely to
be most efficient on dense graphs.

The published algorithm does not include logic to deduplicate elements of Q;
this needs to be added to avoid explosive memory use - eg on a graph with
31 vertices and 44 edges, it generated over 1e6 elements of which only 346
were distinct.

The algorithm also does not appear to work: it is easy for it to get into
a situation in which the current list has
  [ a, b, ... c, d, e, ... f, g, h, ... i, j ]
.. where vertices a, d, j, g have only the edges a-b, a-e, d-c, d-e, j-i,
j-f, g-f, g-h. Then it will find only the four candidates involving zero
to two of the rotations [a .. d] and [g .. j] and then give up. One way
to break that deadlock would be to find a set [ ... k, l, ... m, n, ... ]
where we have edges k-m and l-n to permit an internal rotation that
traverses one or both of the [c, d, e] and [f, g, h] chunks.

=head1 findB ( Graph $graph )

Same as L</find>, a hybrid of it and L</findBFF>: we use an initial pass
to verify certain constraints (connected, every vertex having at least
two edges), and to identify edges that must be or must not be part of any
cycle. We then attempt to extend a list of connected edges as in L</findBFF>,
but with an explicit search when the ends of the list have only interior
edges, rather than stacking all possible permutations.

=cut

sub findLoop {
    my($class, $graph, $method) = @_;
    my $g = ref($method) ? $method->($graph) : $class->$method($graph);
    return undef unless $g;
    return [ $g->find_a_cycle ];
}

sub find {
    no warnings qw{ recursion };

    my($class, $g) = @_;
    return undef unless $g->vertices > 2 && $g->is_biconnected;
    my @d;
    while (1) {
        # Classify the vertices by their degree: if less than 2, we cannot have
        # a cycle; if exactly two, this vertex has no choices; if greater than
        # 2 it does have choices.
        @d = ();
        for my $v ($g->vertices) {
            my $d = $g->degree($v);
            return undef if $d < 2;
            $d[$d]{$v} = 1;
        }

        if (@d == 3) {
            # Every vertex has exactly two edges; we have a cycle iff the
            # remaining graph is connected.
            return undef unless $g->is_connected;
            return $g;
        }

        # For each >2 vertex: if it has two =2 neighbours, it must connect
        # to them, so we can disallow its other edges; if has has more than
        # two, no cycle is possible.
        my $h = $g->copy;
        my $improve = 0;
        my $eq2 = $d[2];
        for my $dv (grep defined, @d[3 .. $#d]) {
            for my $v (keys %$dv) {
                my($req, @non) = 0;
                for my $u ($g->neighbours($v)) {
                    if ($eq2->{$u}) {
                        ++$req;
                    } else {
                        push @non, $u;
                    }
                }
                return undef if $req > 2;
                next if $req < 2;
                $h->delete_edge($v, $_) for @non;
                $improve += @non;
            }
        }

        last unless $improve;
        $g = $h;
    }

    # Could not reduce the graph entirely deterministically, so pick a vertex
    # and try each possible way to reduce it.
    my($dv) = grep defined, @d[3 .. $#d];
    my($v) = keys %$dv;
    my @u = $g->neighbours($v);
    my($eq2) = grep $d[2]{$_}, @u;

    if (defined $eq2) {
        # We already have one fixed edge, so need only pick the other one
        my %suppress = ($eq2 => 1);
        for my $u (@u) {
            next if $u == $eq2;
            local $suppress{$u} = 1;
            my $h = $g->copy;
            $h->delete_edge($v, $_) for grep !$suppress{$_}, @u;
            my $result = $class->find($h);
            return $result if $result;
        }
        return undef;
    }

    my %suppress;
    for my $u1 (@u) {
        local $suppress{$u1} = 1;
        for my $u2 (@u) {
            next if $u2 == $u1;
            local $suppress{$u2} = 1;
            my $h = $g->copy;
            $h->delete_edge($v, $_) for grep !$suppress{$_}, @u;
            my $result = $class->find($h);
            return $result if $result;
        }
    }
    return undef;
}

=for algorithm

Input: a connected graph G = (V_n, E) of minimum degree at least 2.

let P be the path (1, w) where w = min {v: {1, v} \in E};
k := 1
L1: {      # stage k begins here
    Q_1 := P_k; s := 1; t := 1; \delta(Q_1) := 0;
    # \delta(Q_s) is the number of rotations in the sequence constructing
    # Q_s from Q_1
    repeat {
        let path Q_s have endpoints w_0, w_1 where w_0 < w_1;
        for i = 0, 1 {
            # Suppose that the edges incident with w_i and not contained in Q_s
            # are {w_i, x_1}, ..., {w_i, x_p} where x_1 < x_2 < ... < x_p;
            for j = 1 to p {
                if (x_j is not on Q_s) {
                    P_{k+1} := Q_s + {w_i, x_j};    # extension
                    k := k + 1;
                    next L1
                } elsif (x_j = w_{1-i}) {
                    let C be the cycle Q_s + {w_0, w_1};
                    if (C is a Hamiltonian cycle) {
                        return C
                    } else {
                        # starting from w_0, let u be the first vertex along
                        # Q_s which is adjacent to some vertex not in C; let v
                        # be the lowest numbered neighbour of u not in C and
                        # let u_1 and u_2 be the neighbours of u on C with
                        # u_1 < u_2, then:
                        P_{k+1} := C + {u,v} - {u,u_1}; k := k + 1;
                        next L1   # cycle extension
                    }
                } else {
                    t := t + 1;
                    Q_t := ROTATE(Q_s, {w_i, x_j});
                    \delta(Q_t) := \delta(Q_s) + 1
                }
            }
        }
        s := s + 1;
    } until \delta(Q_s) >= 2 * T + 1;
    # where T = ceil(log(n) / (log(d) - log(log(d)))) + 1 and d = 2 * m / n
    return undef
}

=cut

my $maxQ = 1;

sub findBFF {
    my($class, $g) = @_;
    my(@P, @Q, @delta, %Q);
    my $n = $g->vertices;
    my $m = $g->edges;
    return undef unless $n > 2 && $m > 2 && $g->is_connected;

    my $d = 2 * $m / $n;
    my $T = int(log($n) / (log($d) - log(log($d)))) + 2;
    my $lim_delta = 2 * $T + 1;

    @P = @{ _first_edge($g) };
    STAGE: {
print "init: [@P]\n";
        my($s, $t) = (0, 0);
        @Q = join '', map chr($_), @P;
        %Q = map +($_ => 1), @Q;
        $delta[$s] = 0;
        while (defined $Q[$s] && $delta[$s] <= $lim_delta) {
            my $q = [ map ord($_), split //, $Q[$s] ];
print "try: [@$q]\n";
            my %known = map +($q->[$_] => $_), 0 .. $#$q;
            my @w = ($q->[0], $q->[-1]);
print "  ends: @w\n";
            for my $i (0, 1) {
                my $adj = ($i == 0) ? $q->[1] : $q->[-2];
                for my $x_j (sort $g->neighbours($w[$i])) {
print "  try $w[$i]-$x_j\n";
print("  .. is adjacent\n"),
                    next if $x_j eq $adj;
                    if (!defined $known{$x_j}) {
                        @P = ($i == 0) ? ($x_j, @$q) : (@$q, $x_j);
                        @P = reverse @P if $P[0] gt $P[-1];
print "  .. is new, restart with new list\n";
                        redo STAGE;
                    }
                    if ($x_j eq $w[1 - $i]) {
print "  .. loops to other end\n";
                        if (@$q == $n) {
print "@{[ 0 + @Q ]} orderings tested, of which @{[ 0 + keys %Q ]} distinct\n";
                            return _result($g, $q);
                        }
                        for my $qi (0 .. $#$q) {
                            my $node = $q->[$qi];
                            my($new) = sort grep !defined $known{$_},
                                    $g->neighbours($node);
                            next unless defined $new;
                            @P = (@$q[$qi + 1 .. $#$q], @$q[0 .. $qi], $new);
                            @P = reverse @P if $P[0] gt $P[-1];
print "  ... so we cycle to new node at $node allowing new $new\n";
                            redo STAGE;
                        }
print "  ... but we find no new node\n";
                        next;
                    }
print "  .. so try for rotation\n";
print("  ... but stop, delta $delta[$s] reaches limit $lim_delta\n"),
                    next if $delta[$s] >= $lim_delta;
                    my $xi = $known{$x_j};
                    my @newQ = (($i == 0)
                        ? (reverse(@$q[0 .. $xi - 1]), @$q[$xi .. $#$q])
                        : (@$q[0 .. $xi], reverse(@$q[$xi + 1 .. $#$q]))
                    );
                    @newQ = reverse @newQ if $newQ[0] gt $newQ[-1];
                    my $packed = join '', map chr($_), @newQ;
print "  ... have rotation to index $xi\n";
print("  ... but we've seen it before\n"),
                    next if $Q{$packed}++;
print "  ... so add it\n";
                    $Q[++$t] = $packed;
                    $delta[$t] = $delta[$s] + 1;
                }
            }
            ++$s;
        }
    }
print "@{[ 0 + @Q ]} orderings tested, of which @{[ 0 + keys %Q ]} distinct\n";
print "delta $delta[-1] of $lim_delta\n";
    return undef;
}

sub _result {
    my($g, $path) = @_;
    return $g if $g->edges == @$path;
    my $h = $g->subgraph([]);
    for (0 .. $#$path - 1) {
        $h->add_edge($path->[$_], $path->[$_ + 1]);
    }
    $h->add_edge($path->[0], $path->[-1]);
    return $h;
}

sub _first_edge {
    my($g) = @_;
    my $u;
    for my $w ($g->vertices) {
        $u = $w if !defined($u) || $w lt $u;
    }
    my $v;
    for my $w ($g->neighbours($u)) {
        $v = $w if !defined($v) || $w lt $v;
    }
    return [ $u, $v ];
}

=for results

findB prepass: v=32, e=46 => v=13, e=21, req=6
first extension finds 13 vertices
first extension finds 13 vertices

findB prepass: v=473, e=713 => v=210, e=346, req=67
first extension finds 118 vertices
first extension finds 150 vertices

findB prepass: v=9641, e=14622 => v=3775, e=6524, req=1416
first extension finds 974 vertices
first extension finds 1304 vertices

=cut

sub findB {
    my($class, $g) = @_;
    my $s = {
        g => $g->copy,
        v => {},
        required => $g->subgraph([]),
        elided => [],
        list => [],
        inlist => {},
    };
    return undef unless _findB_prepass($s);
    warn sprintf "findB prepass: v=%d, e=%d => v=%d, e=%d, req=%d\n",
            0 + $g->vertices, 0 + $g->edges,
            0 + $s->{g}->vertices, 0 + $s->{g}->edges,
            0 + $s->{required}->edges;

    _findB_extend($s, 0, max $s->{g}->vertices);
    EXTEND: while (1) {
        my($first, $second) = @{ $s->{list} }[0, -1];
        my @order = sort { $b <=> $a } ($first, $second);
        for my $u (@order) {
            my $v = $s->{v}{$u}{out} // next;
            _findB_extend($s, ($u == $first) ? 0 : 1, $v);
            next EXTEND;
        }

        for my $u (@order) {
            for my $v (sort { $b <=> $a } $s->{g}->neighbours($u)) {
                next if ($u == $first && $v == $s->{list}[1])
                        || ($u == $second && $v == $s->{list}[-2]);
                my $iv = _index_of($s, $v);
                my $w = $s->{list}[ $u == $first ? $iv - 1 : $iv + 1 ];
                next if $s->{required}->has_edge($v, $w);
                next unless $s->{v}{$w}{out};
                _findB_rotate($s, $u == $first ? 0 : 1, $iv);
                next EXTEND;
            }
        }
        last;
    }
    warn sprintf "first extension finds %d vertices\n", 0 + @{ $s->{list} };

    local @$s{qw{ g required }} = map "$_", @$s{qw{ g required }};
    use Data::Dumper; local $Data::Dumper::Indent = 1; print Dumper($s);
    return "";
}

sub _findB_rotate {
    my($s, $end, $i) = @_;
    my $list = $s->{list};
    if ($end == 0) {
        @$list[0 .. $i - 1] = reverse @$list[0 .. $i - 1];
    } else {
        @$list[$i + 1 .. $#$list] = reverse @$list[$i + 1 .. $#$list];
    }
    return;
}

sub _index_of {
    my($s, $v) = @_;
    my $list = $s->{list};
    for (0 .. $#$list) {
        return $_ if $list->[$_] == $v;
    }
    die "cannot find $v in [@$list]";
}

sub _findB_extend {
    my($s, $end, $v) = @_;
    my $u = $s->{list}[$end ? -1 : 0];
    my $w = $s->{v}{$v}{req};

    if ($end) {
        push @{ $s->{list} }, $v;
    } else {
        unshift @{ $s->{list} }, $v;
    }
    $s->{inlist}{$v} = 1;
    for my $x ($s->{g}->neighbours($v)) {
        if ($v == ($s->{v}{$x}{out} // -1)) {
            $s->{v}{$x}{out} = max grep !$s->{inlist}{$_},
                    $s->{g}->neighbours($x);
        }
    }

    _findB_extend($s, $end, $w)
            if defined($w) && (!defined($u) || $w != $u);
    return;
}

sub _findB_prepass {
    my($s) = @_;
    my $g = $s->{g};
    my $req = $s->{required};
    my $el = $s->{elided};

    my @d;
    my $eg = $g->edges;
    my $er = $req->edges;
    while (1) {
        # Any vertex with only two edges requires both of them. We delete
        # that vertex, and replace it with a mandatory edge between its
        # two neighbours
        for my $v ($g->vertices) {
            my $d = $g->degree($v);
            return 0 if $d < 2;
            next if $d > 2;
            my($n1, $n2) = $g->neighbours($v);
            $_->delete_vertex($v) for ($g, $req);
            $_->add_edge($n1, $n2) for ($g, $req);
            push @$el, [ $n1, $v, $n2 ];
            # if a vertex has two required edges, any remaining edges
            # are disallowed
            for my $u ($n1, $n2) {
                my $du = $req->neighbours($u);
                next if $du < 2;
                die "panic: ($u) in $g and $req" if $du > 2;
                for my $w ($g->neighbours($u)) {
                    $g->delete_edge($u, $w) unless $req->has_edge($u, $w);
                }
            }
        }

        my $eg2 = $g->edges;
        my $er2 = $req->edges;
        last if $eg == $eg2 && $er == $er2;
        ($eg, $er) = ($eg2, $er2);
    }

    my $vh = $s->{v};
    for my $v ($g->vertices) {
        my $vs = $vh->{$v} = {};
        my @h = $req->neighbours($v);
        if (@h > 1) {
            die "logic error: unelided vertex with multiple required edges",
                ": ($v) $g $req";
        }
        $vs->{req} = $h[0] if @h;
        my %in = map +($_ => 1), @h;
        $vs->{out} = max grep !$in{$_}, $g->neighbours($v);
    }
    return 1;
}

1;
