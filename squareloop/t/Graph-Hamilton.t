use Test::More;
use strict;
use warnings;

use Graph;
my $class = 'Graph::Hamilton';

require_ok($class);

for my $method (qw{ find findBFF }) {
    for my $t ([
        '', 0, 'empty'
    ], [
        'a', 0, '1 node'
    ], [
        'a=b', 0, '2 nodes'
    ], [
        'a=b,a=c', 0, '3 nodes, 2 edges'
    ], [
        'a=b,a=c,b=c', 1, '3 nodes, complete'
    ], [
        'a=b,a=c,a=d,b=c,c=d', 2, '4 nodes, 5 edges'
    ], [
        'a=b,a=c,b=c,d=e,d=f,e=f', 0, '6 nodes, 2 disjoint loops'
    ], [
        'a=b,a=c,b=c,c=d,d=e,d=f,e=f', 0, '6 nodes, 2 loops with bridge'
    ]) {
        my($graph, $result, $legend) = @$t;
        my $text = "graph '$graph' ($legend)";
        my $g = _graphFromString($graph);
        my $actual = $class->$method($g);

        is("$g", $graph,
                "find() does not modify input $text");

        if ($result == 0) {
            # no cycle
            is($actual, undef, "$text returns undef");
        } elsif ($result == 1) {
            # identity
            is($actual, $g, "$text returns same graph");
        } elsif ($result == 2) {
            # subgraph
            isa_ok($actual, 'Graph::Undirected',
                    "$text returns a Graph::Undirected object");
            is_deeply([ sort $g->vertices ], [ sort $actual->vertices ],
                    "result for $text has the same nodes");
            my @bad_degree = map sprintf('%s: %s', $actual->degree($_), $_),
                grep $actual->degree($_) != 2,
                    $actual->vertices;
            is_deeply(\@bad_degree, [],
                    "each node in result for $text has 2 edges");
            ok($actual->is_connected,
                    "result for $text is connected");
            my @bad_edge = grep !$g->has_edge(@$_), $actual->edges;
            is_deeply(\@bad_edge, [],
                    "each edge in result for $text is in the original");
        } else { die "logic error" }
    }
}

done_testing();

sub _graphFromString {
    my($t) = @_;
    my $g = Graph::Undirected->new;
    for (split /,/, $t) {
        if (/^([^=]+)=([^=]+)\z/) {
            $g->add_edge($1, $2);
        } elsif (!/[=]/) {
            $g->add_vertex($_);
        } else { die "logic error" }
    }
    return $g;
}
