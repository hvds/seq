package Prove::Construct;
use strict;
use warnings;

use lib 'lib';
use Axiom::Expr;
use Prove::LinCom;
use Prove::Line;
use Prove::Construct::Split;
use Prove::Construct::Resolve;

=head1 NAME

Prove::Construct - construct the set of proof lines needed

=head1 DESCRIPTION

Given a target expression of the form C< \min(a+b+...,c+d+...) >,
constructs the lines needed to resolve occurrences of that expression.
The strategy is as follows:

1) apply unarydistrib where needed to split pluslists within integrals
to isolate the expression;

2) for each occurrence of the target expression (iteratively), find
a suitable split point for an enclosing integral;

3) generate additional lemmas as needed to justify the itersplits;

4) generate lemmas for resolution of the C< \min() >, up to integration
of the most deeply nested integral if relevant;

The resulting information is supplied to the L<Prove::State> object.

=cut

sub new {
    my($class, $state, $plus, $minus) = @_;
    my $vare = $state->vare;
    my $vari = $state->vari;
    my @plusv = @$vare{@$plus};
    my @minusv = @$vare{@$minus};

    my %seen;
    ++$seen{$_} for (@$plus, @$minus);
    my @both = grep $seen{$_} > 1, keys %seen;
    my @allv = @$vare{sort keys %seen};
    %seen = map +($_ => 1), @both;
    my @plus = grep !$seen{$_}, @$plus;
    my @minus = grep !$seen{$_}, @$minus;
    my @vars = @$vari{ @plus, @minus };
    my @need = ((map 1, @plus), (map -1, @minus));
    my $pivot = Prove::LinCom->new_splice(\@vars, \@need);
    my $ce = _construct_expr($state->vare, $plus, $minus);
    my $self = bless {
        state => $state,
        needvar => \@vars,
        need => \@need,
        allv => \@allv,
        bothv => [@$vare{@both}],
        posneg => [ \@plusv, \@minusv ],
        pivot => $pivot,
        ce => $ce,
        reqs => [],
        lines => [ ],
    };
    return $self;
}

sub state { shift->{state} }
sub needvar { shift->{needvar} }
sub need { shift->{need} }
sub allv { shift->{allv} }
sub bothv { shift->{bothv} }
sub posneg { shift->{posneg} }
sub pivot { shift->{pivot} }
sub ce { shift->{ce} }
sub reqs { shift->{reqs} }
sub lines { shift->{lines} }
sub lastexpr { shift->lines->[-1]->expr }

sub build {
    my($self, $expr) = @_;
    push @{ $self->lines }, Prove::Line->new({
        type => 'condstart',
        expr => $expr,
    });
    while (1) {
        my $int = $self->find_next or last;
        my $lc = $self->pivot;
        my($what, $args) = $self->dispose($int, $lc);
        if ($what eq 'resolve') {
            $self->resolve($int, $args);
            next;
        }

        my($i, $lcs) = @$args;
        my $fork = $int->ofork;
        if ($fork & ~((1 << $i) - 1)) {
            $self->distrib($int, $i);
            next;
        }
        $self->itersplit($int, $i, $lcs);
    }

    my $tail = $self->gen_tail;
    push @{ $self->lines }, Prove::Line->new({
        type => 'condend',
        expr => $tail,
    });
    push @{ $self->lines }, Prove::Line->new({
        type => 'theorem',
        label => 'P',
        expr => $tail,
    });
    return;
}

sub find_next {
    my($self) = @_;
    my $expr = $self->lastexpr;
    my $ce = $self->ce;
    my $vari = $self->state->vari;
    my $loc;
    my $iter = $expr->iter_locn;
    while (1) {
        my($e, $l) = $iter->() or return;
        next if $e->diff($ce, 0, 2);
        $loc = $l;
        last;
    }
    my $ancestry = $expr->ancestry($loc);
    my(@stack, @iloc);
    my $fork = 0;
    for (0 .. $#$ancestry) {
        my $e = $ancestry->[$_];
        my $t = $e->type;
        if ($t eq 'integral') {
            my($f, $t) = @{ $e->args }[1, 2];
            push @iloc, [@$loc[0 .. $_ - 1]];
            push @stack, [
                Prove::LinCom->parse($f, $vari),
                Prove::LinCom->parse($t, $vari),
            ];
        } elsif ($t eq 'pluslist') {
            $fork |= (1 << $#stack);
        }
    }
    return Prove::Integral->new({
        loc => $loc,
        expr => $ce,
        istack => \@stack,
        iloc => \@iloc,
        ofork => $fork,
    });
}

sub emit {
    my($self) = @_;
    print <<PRE;
*import int
*var A

PRE
    Prove::Construct::Split->emit_lemmas;
    Prove::Construct::Resolve->emit_lemmas;
    Prove::Line->emit_set($self->lines);
}

sub dispose {
    my($self, $int, $lc) = @_;
    my $range = $int->range($lc);
#warn "range @{[ $self->state->mkexpr($lc)->str ]} is [@$range]\n";
    die if $range->[1] <= $range->[0];
    return +('resolve', 1) if $range->[1] <= 0;
    return +('resolve', 2) if $range->[0] >= 0;
    return +('split', $self->find_split($int, $lc));
}

sub find_split {
    my($self, $int, $lc) = @_;
    my($i, $c) = $lc->top;
    my $lc2 = $lc->popone;

    # for lc=b+c-d with d in [0,1-a], we need b+c>0 and b+c<1-a => -1+a+b+c<0
    # lege=0
    my $lc3 = $lc2->addmul($int->i($i)->[($c < 0) ? 0 : 1], $c);
    my $range = $int->range($lc3);
#warn "for $c of var $i at lege 0 try @{[ $self->state->mkexpr($lc3)->str ]} giving range [@$range]\n";
    return $self->find_split($int, $lc3) unless $range->[0] >= 0;

    # lege=1
    $lc3 = $lc2->addmul($int->i($i)->[($c < 0) ? 1 : 0], $c);
    $range = $int->range($lc3);
#warn "for $c of var $i at lege 1 try @{[ $self->state->mkexpr($lc3)->str ]} giving range [@$range]\n";
    return $self->find_split($int, $lc3) unless $range->[1] <= 0;

#warn "will split $i over @{[ $self->state->mkexpr(($c < 0) ? $lc2 : $lc2->negate)->str ]}\n";
    return [ $i, ($c < 0) ? $lc2 : $lc2->negate ];
}

sub itersplit {
    my($self, $int, $i, $lcs) = @_;
    return Prove::Construct::Split->gen_split($self, $int, $i, $lcs);
}

sub resolve {
    my($self, $int, $resolve) = @_;
    return Prove::Construct::Resolve->gen_resolve($self, $int, $resolve);
}

sub distrib {
    my($self, $int, $i) = @_;
    my $expr = $self->lastexpr;
    my $loc = $int->loc;
#warn "try distrib over $i at @{[ $expr->locate($loc)->str ]}\n";
    my $ancestry = $expr->ancestry($loc);

    my $plus_depth = (List::Util::first(sub {
        $ancestry->[$_]->type eq 'pluslist'
    }, reverse 0 .. $#$ancestry - 1)) // die "panic: no pluslist for distrib";
    my $int_loc = [ @$loc[0 .. $plus_depth - 2] ];
    my $int_expr = $expr->locate($int_loc);
    $int_expr->type eq 'integral' or die "panic: type @{[ $int_expr->type ]}";
    my($var, $from, $to, $plus_expr) = @{ $int_expr->args };
    $plus_expr->type eq 'pluslist' or die "panic: type @{[ $plus_expr->type ]}";

    # CHECKME: do we ever want to keep some plus exprs grouped?
    my @repl_args = map Axiom::Expr->new({
        type => 'integral',
        args => [ $var->copy, $from->copy, $to->copy, $_->copy ],
    }), @{ $plus_expr->args };
    my $repl = Axiom::Expr->new({
        type => 'pluslist',
        args => \@repl_args,
    });

    my $parent = $expr->parent($int_loc);
    my $start = 0;
    if ($parent && $parent->type eq 'pluslist') {
        my $off = pop(@$int_loc) - 1;
        my $args = $parent->args;
        $repl = Axiom::Expr->new({
            type => 'pluslist',
            args => [
                (map $_->copy, @$args[0 .. $off - 1]),
                @repl_args,
                (map $_->copy, @$args[$off + 1 .. $#$args]),
            ],
        });
        $start = $off + 1;
    }
    my $new_expr = $expr->substitute($int_loc, $repl);
    push @{ $self->lines }, Prove::Line->new({
        type => 'unarydistrib',
        expr => $new_expr,
    });
    return;
}

sub gen_tail {
    my($self) = @_;
    my $result = $self->lastexpr;
    my $ante = $self->state->expr;
    return Axiom::Expr->new({
        type => 'implies',
        args => [ $ante->copy, $result->copy ],
    });
}

sub _construct_expr {
    my($ve, $plus, $minus) = @_;
    my @pe = @$ve{@$plus};
    my @me = @$ve{@$minus};
    return Axiom::Expr->new({
        type => 'min',
        args => [
            (@pe > 1 ? Axiom::Expr->new({
                type => 'pluslist',
                args => [ map $_->copy, @pe ],
            }) : $pe[0]->copy),
            (@me > 1 ? Axiom::Expr->new({
                type => 'pluslist',
                args => [ map $_->copy, @me ],
            }) : $me[0]->copy),
        ],
    });
}

1;
