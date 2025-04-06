package Prove::Construct::Split;
use strict;
use warnings;

use List::Util qw{};

use lib 'lib';
use Axiom::Expr;
use Prove::Line;
use Prove::Lemma;
Prove::Lemma->register('split');

=head1 NAME

Prove::Construct::Split - construct proof lines to split an integral

=head1 DESCRIPTION

Given a request to split an integral in the expression in a particular
way, records the itersplit needed to split at that point, along with proof
of the accompanying lemma to justify it if needed.

=cut

sub new {
    my($class, $construct, $int, $who, $split) = @_;
    my $state = $construct->state;
    return bless {
        state => $state,
        construct => $construct,
        oint => $int,
        who => $who,
        which => $state->vare->{ $state->vars->[$who] },
        osplit => $split,
        esplit => $state->mkexpr($split),
    }, $class;
}

sub construct { shift->{construct} }
sub oint { shift->{oint} }
sub who { shift->{who} }    # index of variable to split over
sub which { shift->{which} }  # expr of variable to split over
sub osplit { shift->{osplit} }
sub esplit { shift->{esplit} }
sub state { shift->construct->state }
sub lines { shift->construct->lines }

sub gen_split {
    my($class, $construct, $int, $who, $split) = @_;
#print STDERR "split $who via @{[ $construct->state->mkexpr($split)->str ]}: "; $int->debug;
    my $self = $class->new($construct, $int, $who, $split);
    my $lemma = $self->find_lemma;
    push @{ $construct->lines }, Prove::Line->new({
        type => 'itersplit',
        args => [ $lemma ? ("*with $lemma") : () ],
        expr => $self->gen_expr,
    });
    return;
}

sub gen_expr {
    my($self) = @_;
    my $int = $self->oint;
    my $esplit = $self->esplit;
    my @loc = @{ $int->loci($self->who) };
    my $off = pop @loc;
    my $base = $self->construct->lastexpr;
    my $expr = $base->locate(\@loc);
    my $old = $expr->args->[$off - 1];
    $old->type eq 'integral' or die "split target is not integral for $off at : @{[ $expr->str ]}\n";;
    my($var, $from, $to, $iexpr) = @{ $old->args };
    my $newl = Axiom::Expr->new({
        type => 'integral',
        args => [ map $_->copy, ($var, $from, $esplit, $iexpr) ],
    });
    my $newr = Axiom::Expr->new({
        type => 'integral',
        args => [ map $_->copy, ($var, $esplit, $to, $iexpr) ],
    });
    my $repl;
    if ($expr->type eq 'pluslist') {
        my $args = $expr->args;
        $repl = Axiom::Expr->new({
            type => 'pluslist',
            args => [
                (map $_->copy, @$args[0 .. $off - 2]),
                $newl, $newr,
                (map $_->copy, @$args[$off .. $#$args]),
            ],
        });
    } else {
        push @loc, $off;
        $repl = Axiom::Expr->new({
            type => 'pluslist',
            args => [ $newl, $newr ],
        });
    }
    my $result = $base->substitute(\@loc, $repl);
    return $result;
}

sub find_lemma {
    my($self) = @_;
    my $int = $self->oint;
    my $split = $self->osplit;
    my $range = $int->i($self->who);
    my $var = $self->which->name;

    my $olemma = Prove::Lemma->new($self->state, 'split', "S$var", $int,
            $split, 0, $range->[0]);
    $olemma->extend($split, 1, $range->[1]);
    return $olemma->lemma;
}

sub emit_lemmas {
    my($class) = @_;
    my $lemmas = Prove::Lemma->lemmas('split');
    for my $sig (map $_->[0], sort {
        $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]
    } map [ $_, /^S([a-z])(\d+)$/ ], keys %$lemmas) {
        Prove::Line->emit_set($lemmas->{$sig});
    }
    return;
}

1;
