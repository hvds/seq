package Prove::Integral;
use strict;
use warnings;

use lib 'lib';
use Prove::LinCom;
use Prove::Poly;

=head1 NAME

Prove::Integral - representation of a stack of nested integrals

=head1 DESCRIPTION

Uses an arrayref of the form C<< [ <Prove::LinCom>* <Axiom::Expr> ] >>.

=cut

my $tag = 1;

sub new {
    my($class, $data) = @_;
    my $self = bless $data, $class;
#   $self->{tag} = join ':', grep defined($_), ($self->tag, $tag++);
    $self->{tag} = $tag++;
    $self->{expr} //= $self->istack->[-1];
    return $self;
}

sub expr { shift->{expr} }
sub scexpr { shift->{scexpr} }
sub loc { shift->{loc} }
sub istack { shift->{istack} }
sub iloc { shift->{iloc} }
sub ofork { shift->{ofork} }
sub tag { shift->{tag} }

sub clone {
    my($self) = @_;
    return bless {
        expr => $self->expr->copy,
        loc => [ @{ $self->loc } ],
        istack => [ map [ map $_->copy, @$_ ], @{ $self->istack } ],
        iloc => [ map [ @$_ ], @{ $self->iloc } ],
        ofork => $self->ofork,
#       tag => join(':', grep defined($_), ($self->tag, $tag++)),
        tag => $tag++,
    }, ref($self);
}

sub _expr {
    my($lc) = @_;
    my $e = join '+', map {
        my($i, $c) = @$_;
        $c . ($i ? chr(ord('a') + $i - 1) : '')
    } $lc->all_pair;
    $e =~ s/\+-/-/g;
    $e =~ s/(^|[-+])1([a-z])/$1$2/g;
    $e;
}

sub debugs {
    my($self) = @_;
    my $stack = $self->istack;
    my @part = map {
        my $v = chr(ord('a') + $_);
        sprintf '%s=%s^%s', $v, map _expr($_), @{ $stack->[$_] };
    } 0 .. $#$stack;
    return "@{[ $self->tag ]}: [@{[ join q{, }, @part ]}]";
}

sub debug {
    my($self, $check) = @_;
    my($stack, $iloc, $fork) = @$self{qw{ istack iloc ofork }};
    my @part = map {
        my $v = chr(ord('a') + $_);
        my $int = sprintf '%s=%s^%s', $v, map _expr($_), @{ $stack->[$_] };
        my $plus = ($fork & (1 << ($_+1))) ? '+' : '';
        my $loc = "[@{ $iloc->[$_] }]";
        "$loc $int $plus";
    } 0 .. $#$stack;
    warn "@{[ $self->tag ]}: [\n";
    warn "  $_\n" for @part;
    warn "  [@{ $self->loc }] @{[ $self->expr->str ]}\n";
    warn "]\n";
    return $self;
}

sub set_stack {
    my($self, $i, $which, $repl) = @_;
    $self->istack->[$i - 1][$which] = $repl;
    return;
}

sub i {
    my($self, $i) = @_;
    my $istack = $self->istack;
    die if $i < 1 || $i > @$istack;
    return $istack->[$i - 1];
}

sub loci {
    my($self, $i) = @_;
    my $iloc = $self->iloc;
    die if $i < 1 || $i > @$iloc;
    return $iloc->[$i - 1];
}

sub match {
    my($self, $match) = @_;
    my $matched = 0;
    my $iter = $self->expr->iter_tree;
    while (my $e = $iter->()) {
        # match names only
        return 1 unless $match->diff($e, 1, 2);
    }
    return 0;
}

sub range {
    my($self, $lc) = @_;
    my $istack = $self->istack;
    my($min, $max) = ([ @$lc ], [ @$lc ]);
    while (@$min > 1) {
        pop(@$min), next if $min->[-1] == 0;
        my($i, $c) = ($#$min, $min->[-1]);
        my $targ = $istack->[$i - 1][ ($c < 0) ? 1 : 0 ];
        $min->[$_] += $c * $targ->[$_] for 0 .. $#$targ;
        pop @$min;
    }
    while (@$max > 1) {
        pop(@$max), next if $max->[-1] == 0;
        my($i, $c) = ($#$max, $max->[-1]);
        my $targ = $istack->[$i - 1][ ($c < 0) ? 0 : 1 ];
        $max->[$_] += $c * $targ->[$_] for 0 .. $#$targ;
        pop @$max;
    }
    return [ $min->[0] // 0, $max->[0] // 0 ];
}

sub integrate {
    my($self) = @_;
    my $poly = Prove::Poly->from_sc($self->scexpr);
    my $stack = $self->istack;
    for (reverse (0 .. $#$stack)) {
        my $range = $stack->[$_];
        my $eval = $poly->integrate($_ + 1);
        $poly = $eval->inteval($_ + 1, $range);
    }
    return $poly->const;
}

1;
