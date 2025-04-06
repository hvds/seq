package Prove::Construct::Resolve;
use strict;
use warnings;

use List::Util qw{};

use lib 'lib';
use Axiom::Expr;
use Prove::Line;
use Prove::Lemma;
Prove::Lemma->register('rgive');
Prove::Lemma->register('requate');

=head1 NAME

Prove::Construct::Resolve - construct proof lines to resolve C< \min(x, y) >

=head1 DESCRIPTION

Given a request to resolve a particular instance of C< \min(x, y) >,
records the give/equate pair to apply the resolution, along with proof
of the accompanying lemma to justify the give if needed.

=cut

sub new {
    my($class, $construct, $int, $resolve) = @_;
    my $state = $construct->state;
    return bless {
        state => $state,
        construct => $construct,
        oint => $int,
        resolve => $resolve,
    }, $class;
}

sub construct { shift->{construct} }
sub oint { shift->{oint} }
sub resolve { shift->{resolve} }
sub state { shift->construct->state }
sub lines { shift->construct->lines }
sub need { shift->{need} }
sub top {
    my($self) = @_;
    return $self->{top} //= ($self->construct->pivot->top)[0];
}
sub topn {
    my($self) = @_;
    return $self->{topn} //= $self->state->vars->[$self->top];
}
sub topvar {
    my($self) = @_;
    return $self->{topvar} //= $self->state->vare->{$self->topn};
}

sub gen_resolve {
    my($class, $construct, $int, $resolve) = @_;
    my $self = $class->new($construct, $int, $resolve);
    my $lemma = $self->find_lemma;
    push @{ $construct->lines }, Prove::Line->new({
        type => 'give',
        args => [ $lemma ? ("*with $lemma") : () ],
        expr => $self->gen_give,
    });
    my $equate = $self->find_equate;
    push @{ $construct->lines }, Prove::Line->new({
        type => 'equate',
        args => [ $equate ],
        expr => $self->gen_equate,
    });
    return;
}

sub gen_give {
    my($self) = @_;
    my $expr = $self->construct->lastexpr;
    my $iter = $expr->iter_locn($self->oint->loc);
    my $ce = $self->construct->ce;
    my $loc;
    while (1) {
        my($e, $l) = $iter->();
        next if $ce->diff($e, 0, 2);
        $loc = $l;
        last;
    }
    my $repl = Axiom::Expr->new({
        type => 'given',
        args => [
            $self->need->copy,
            $ce->copy,
        ],
    });
    return $expr->substitute($loc, $repl);
}

sub llr {
    my($self) = @_;
    return $self->{llr} //= do {
        my $lege = ($self->resolve == 1) ? 0 : 1;
        my $lc = $self->construct->pivot;   # a+d-e-f
        my($i, $c) = $lc->top;              # -f
        $c == -1 or die;
        my $left = Prove::LinCom->new([])->addi($i, -$c);   # +f
        my $right = $lc->popone;            # a+d-e
        [ $left, $lege, $right ];
    };
}

sub find_lemma {
    my($self) = @_;
    my $int = $self->oint;
    my($left, $lege, $right) = @{ $self->llr };
    my $olemma = Prove::Lemma->new($self->state, 'rgive', 'GM', $int,
            $left, $lege, $right);
    $self->{need} = $olemma->need;
    return $olemma->lemma;
}

sub gen_equate {
    my($self) = @_;
    my $expr = $self->construct->lastexpr;
    my $iter = $expr->iter_locn($self->oint->loc);
    my $loc;
    while (1) {
        my($e, $l) = $iter->();
        next unless $e->type eq 'given';
        $loc = $l;
        last;
    }
    my(undef, $lege, undef) = @{ $self->llr };
    my $posneg = $self->eposneg;
    return $expr->substitute($loc, $posneg->[$lege ? 1 : 0]->copy);
}

sub find_equate {
    my($self) = @_;
    my $name = ($self->resolve == 1) ? 'Min1' : 'Min2';
    Prove::Lemma->lemmas('requate')->{$name} //= $self->make_min($name);
    return $name;
}

sub emit_lemmas {
    my($class) = @_;
    my $lemmas = Prove::Lemma->lemmas('rgive');
    for my $sig (map $_->[0], sort {
        $a->[1] <=> $b->[1]
    } map [ $_, /^GM(\d+)$/ ], keys %$lemmas) {
        Prove::Line->emit_set($lemmas->{$sig});
    }
    $lemmas = Prove::Lemma->lemmas('requate');
    for my $sig (map $_->[0], sort {
        $a->[1] <=> $b->[1]
    } map [ $_, /^Min(\d+)$/ ], keys %$lemmas) {
        Prove::Line->emit_set($lemmas->{$sig});
    }
    return;
}

sub eposneg {
    my($self) = @_;
    return $self->{eposneg} //= do {
        my($pos, $neg) = @{ $self->construct->posneg };
        [
            ((@$pos > 1) ? Axiom::Expr->new({
                type => 'pluslist',
                args => [ map $_->copy, @$pos ],
            }) : $pos->[0]->copy),
            ((@$neg > 1) ? Axiom::Expr->new({
                type => 'pluslist',
                args => [ map $_->copy, @$neg ],
            }) : $neg->[0]->copy),
        ];
    };
}

sub _wrap {
    my($e, $va) = @_;
    $e = Axiom::Expr->new({
        type => 'forall',
        args => [ $_->copy, $e->copy ],
    }) for reverse @$va;
    return $e;
}

sub make_min {
    my($self, $name) = @_;
    my @line;
    my $line = sub {
        my($type, $args, $label, $expr) = @_;
        push @line, Prove::Line->new({
            type => $type,
            args => $args,
            label => $label,
            expr => $expr,
        });
    };
    my($left, $lege, $right) = @{ $self->llr };
    $left = $self->state->mkexpr($left);
    $right = $self->state->mkexpr($right);
    my $cond = Axiom::Expr->new({
        type => ($lege ? 'rle' : 'rge'),
        args => [ $left->copy, $right->copy ],
    });
    $line->('condstart', undef, undef, $cond);
    $line->('negate', undef, undef, Axiom::Expr->new({
        type => ($lege ? 'rge' : 'rle'),
        args => [ $left->negate->clean, $right->negate->clean ],
    }));
    my($pos, $neg) = @{ $self->eposneg };
    my $bothv = $self->construct->bothv;
    my $allv = $self->construct->allv;
    if ($left->negate->diff($pos, 0, 2)) {
        my $add = _wrap(Axiom::Expr->new({
            type => +($lege ? 'rge' : 'rle'),
            args => [ $pos->copy, $neg->copy ],
        }), $bothv);
        $line->('add', undef, undef, $add);
    }
    my $min = Axiom::Expr->new({
        type => 'min',
        args => [ $pos->copy, $neg->copy ],
    });
    my $result = Axiom::Expr->new({
        type => 'req',
        args => [
            $min->copy,
            ($lege ? $neg->copy : $pos->copy),
        ],
    });
    my $args = [ $lege ? 'int.M2' : 'int.M1' ];
    $line->('ponens', $args, undef, _wrap($result, $bothv));
    my $condend = _wrap(Axiom::Expr->new({
        type => 'implies',
        args => [
            $cond->copy,
            $result->copy,
        ],
    }), $allv);
    $line->('condend', undef, undef, $condend);
    my $given = _wrap(Axiom::Expr->new({
        type => 'req',
        args => [
            Axiom::Expr->new({
                type => 'given',
                args => [ $cond->copy, $min->copy ],
            }),
            ($lege ? $neg->copy : $pos->copy),
        ],
    }), $allv);
    $line->('given', undef, undef, $given->copy);
    $line->('theorem', undef, $name, $given->copy);
    return \@line;
}

1;
