package Prove::State;
use strict;
use warnings;

use lib 'lib';
use Axiom::Expr;
use Prove::Construct;
use Prove::Integral;
use Prove::LinCom;

=head1 NAME

Prove::State - status of proof

=head1 DESCRIPTION

=cut

sub new {
    my($class, $expr) = @_;
    my($vars, $vari, $vare) = _find_vars($expr);
    my $self = bless {
        expr => $expr,
        vars => $vars,
        vari => $vari,
        vare => $vare,
    };
    return $self;
}

sub expr { shift->{expr} }
sub vars { shift->{vars} }  # [ undef, name* ]
sub vari { shift->{vari} }  # { name => index }
sub vare { shift->{vare} }  # { name => varexpr }

sub name {
    my($self, $i) = @_;
    return $self->vars->[$i - 1];
}

sub _find_vars {
    my($expr) = @_;
    my @vars = (undef);     # const goes first
    my(%vari, %vare);
    my $capture = sub {
        my($v) = @_;
        my $n = $v->name;
        $vari{$n} //= do {
            push @vars, $n;
            $vare{$n} = $v;
            $#vars;
        };
        return;
    };
    my $recurse;
    $recurse = sub {
        my($e) = @_;
        my $t = $e->type;
        if ($t eq 'integral') {
            $capture->($e->args->[0]);
            return $recurse->($e->args->[3]);
        } elsif ($t eq 'pluslist') {
            for (@{ $e->args }) {
                return 1 if $recurse->($_);
            }
            return 1;
        } elsif ($t eq 'name') {
            $capture->($e);
            return 0;
        } elsif ($e->is_const) {
            return 0;
        } else {
            for (@{ $e->args }) {
                return 1 if $recurse->($_)
            }
            return 0;
        }
    };
    $recurse->($expr->args->[1]);   # skip 'A = ...'
    return +(\@vars, \%vari, \%vare);
}

sub construct {
    my($self, $plus, $minus) = @_;
    my $c = Prove::Construct->new($self, $plus, $minus);
    $c->build($self->expr);
    $c->emit;
    $self->{expr} = $c->lastexpr;
    return;
}

sub mkexpr {
    my($self, $lc) = @_;
    my($vars, $vare) = ($self->vars, $self->vare);
    my $args = [ map {
        my($i, $c) = @$_;
        ($i == 0) ? Axiom::Expr->new_const($c)
            : ($c == 1) ? $vare->{$vars->[$i]}->copy
            : Axiom::Expr->new({
                type => 'mullist',
                args => [
                    Axiom::Expr->new_const($c),
                    $vare->{$vars->[$i]}->copy,
                ],
            });
    } $lc->all_pair ];
    return +(@$args > 1) ? Axiom::Expr->new({
        type => 'pluslist',
        args => $args,
    }) : $args->[0];
}

sub mkcond {
    my($self, $int, $var, $lege) = @_;
    my($vars, $vare) = ($self->vars, $self->vare);
    my $lc = $int->i($var)->[$lege];
    my $sum = $self->mkexpr($lc);
    return Axiom::Expr->new({
        type => ($lege ? 'rle' : 'rge'),
        args => [ $vare->{$vars->[$var]}->copy, $sum ],
    });
}

1;
