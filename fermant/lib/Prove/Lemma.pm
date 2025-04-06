package Prove::Lemma;
use strict;
use warnings;

use List::Util qw{};

use lib 'lib';
use Axiom::Expr;
use Prove::Line;

=head1 NAME

Prove::Lemma - prove a lemma

=head1 DESCRIPTION

Given a request to show that C< left ge|le right > in the context of C<int>,
where C<left> and C<right> are L<Prove::LinCom> objects and C<int> is a
L<Prove::Integral>, one of a) determines that no lemma is needed, b) finds
an existing lemma that proves what is needed, or c) constructs a new lemma
to prove what is needed given the constraints from the nested integrals in
C<int>.

Classes must register with a typename, and provide a callback to construct
names for lemmas of that type. Constructed lemmas will be saved until
explicitly cleared.

=cut

my %cache;
my %registry;

sub new {
    my($class, $state, $type, $name, $int, $left, $lege, $right) = @_;
#warn "lemma for $name: @{[ $state->mkexpr($left)->str ]} @{[ $lege ? q{<=} : q{>=} ]} @{[ $state->mkexpr($right)->str ]}\n";
    return bless {
        state => $state,
        type => $type,
        %{ $registry{$type} // die "unknown type" },
        name => $name,
        oint => $int,
        scneed => [ [ $lege, $left, $right ] ],
        notneed => {},
    }, $class;
}

sub extend {
    my($self, $left, $lege, $right) = @_;
#my $state = $self->state; warn ".. extend with: @{[ $state->mkexpr($left)->str ]} @{[ $lege ? q{<=} : q{>=} ]} @{[ $state->mkexpr($right)->str ]}\n";
    push @{ $self->scneed }, [ $lege, $left, $right ];
    return;
}

sub register {
    my($class, $type) = @_;
    $registry{$type} //= {
        lemmas => {},       # name => [ Prove::Line ]
        lookup => {         # expr->str => name
            '' => ''        # we give an empty name if nothing needed
        },
    };
    return;
}

sub state { shift->{state} }
sub type { shift->{type} }
sub name { shift->{name} }
sub lemmas {
    my($self, $type) = @_;
    return $type ? $registry{$type}{lemmas} : $self->{lemmas};
}
sub lookup { shift->{lookup} }
sub oint { shift->{oint} }
sub scneed { shift->{scneed} }
sub notneed { shift->{notneed} }

sub lemma {
    my($self) = @_;
    return $self->{lemma} //= do {
        my $targ = $self->target;
        my $sig = $targ ? $targ->str : '';
        $self->lookup->{$sig} //= $self->make_lemma($targ);
    };
}

#
# Nothing is needed to prove a particular target if it is of the form
# var op expr, and we have exactly that provided by an integral constraint.
# In that case we don't need it as an antecedent or postcedent in the lemma,
# but may still need to provide the full list/expr of things proved (eg to
# provide to a 'give').
#
sub _notneed {
    my($left, $neg, $pos) = @_;
    return 0 if @$neg || @$pos != 1;
    my($i, $c) = @{ $pos->[0] };
    return 0 if $c != 1;
    my $diff = $left->addi($i, -$c);
    return 0 unless $diff->is_const;
    return $diff->const == 0 ? 1 : 0;
}

sub req {
    my($self) = @_;
    return $self->{req} //= do {
        my $int = $self->oint;
        my $need = $self->scneed;
        my @req;
        for my $ni (0 .. $#$need) {
            my($lege, $left, $right) = @{ $need->[$ni] };
#warn "find req for @{[ $self->state->mkexpr($left)->str ]} @{[ $lege ? q{<=} : q{>=} ]} @{[ $self->state->mkexpr($right)->str ]}\n";
            my $lc = $left->addmul($right->negate, 1);
            my(@neg, @pos);
            while (!$lc->is_const) {
#warn "tester @{[ $self->state->mkexpr($lc)->str ]}\n";
                my($i, $c) = $lc->top;
                push @{ ($c < 0) ? \@neg : \@pos }, [ $i, $c ];
                my $efflege = $lege;
                $efflege = 1 - $efflege if $c < 0;
#warn "stash [ $i, $c ] given @{[ $self->state->mkcond($int, $i, $efflege)->str ]}\n";
                $lc = $lc->addi($i, -$c);
                $lc = $lc->addmul($int->i($i)->[$efflege], $c);
            }
#warn "done with const @{[ $lc->const ]}\n";
die unless $lc->const == 0;
            if (_notneed($left, \@neg, \@pos)) {
                $self->notneed->{$ni} = 1;
            } else {
                push @req, [ $lege, [ reverse @neg ], [ reverse @pos ] ];
            }
        };
        \@req;
    }
}

sub aneed {
    my($self) = @_;
    return $self->{aneed} //= do {
        my $state = $self->state;
        [ map {
            my($lege, $left, $right) = @$_;
            Axiom::Expr->new({
                type => ($lege ? 'rle' : 'rge'),
                args => [ map $state->mkexpr($_), ($left, $right) ],
            });
        } @{ $self->scneed } ];
    };
}

sub arneed {
    my($self) = @_;
    return $self->{arneed} //= do {
        my $state = $self->state;
        my $need = $self->scneed;
        $self->req; # make sure we've checked what is really needed
        my $notneed = $self->notneed;
        [ map {
            my($lege, $left, $right) = @$_;
            Axiom::Expr->new({
                type => ($lege ? 'rle' : 'rge'),
                args => [ map +($state // die)->mkexpr($_), ($left, $right) ],
            });
        } @$need[grep !$notneed->{$_}, 0 .. $#$need] ];
    };
}

sub need {
    my($self) = @_;
    return $self->{need} //= do {
        my $aneed = $self->aneed;
        (@$aneed > 1) ? Axiom::Expr->new({
            type => 'andlist',
            args => $aneed,
        }) : $aneed->[0];
    };
}

sub rneed {
    my($self) = @_;
    return $self->{rneed} //= do {
        my $arneed = $self->arneed;
        (@$arneed > 1) ? Axiom::Expr->new({
            type => 'andlist',
            args => $arneed,
        }) : $arneed->[0];
    };
}

sub hcond {
    my($self) = @_;
    return $self->{hcond} //= do {
        my $state = $self->state;
        my $int = $self->oint;
        my %cond;
        for (@{ $self->req }) {
            my($lege, $neg, $pos) = @$_;
            for ([ -1, $neg ], [ 1, $pos ]) {
                my($sign, $subreqs) = @$_;
                my $eff = ($sign < 0) ? (1 - $lege) : $lege;
                for my $subreq (@$subreqs) {
                    my $var = $subreq->[0];
                    my $sig = "$var:$eff";
                    $cond{$sig} //= $state->mkcond($int, $var, $eff);
                }
            }
        }
        \%cond;
    };
}

sub icond {
    my($self, $i, $eff) = @_;
    return $self->hcond->{"$i:$eff"};
}

sub cond {
    my($self) = @_;
    return $self->{cond} //= do {
        my $cond = $self->hcond;
        my @cond = map $cond->{$_->[0]}, sort {
            $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]
        } map [ $_, split /:/, $_ ], keys %$cond;
        (@cond > 1) ? Axiom::Expr->new({
            type => 'andlist',
            args => \@cond,
        }) : $cond[0];
    };
}

sub var {
    my($self) = @_;
    return $self->{var} //= do {
        my %var;
        $self->cond->walk_tree(sub {
            my($e) = @_;
            # FIXME: make sure it's a local variable
            $var{$e->name} //= $e if $e->type eq 'name';
        });
        [ @var{sort keys %var} ];
    };
}

sub target {
    my($self) = @_;
    return $self->{target} //= ($self->rneed) ? do {
        my $targ = Axiom::Expr->new({
            type => 'implies',
            args => [ $self->cond->copy, $self->rneed->copy ],
        });
        $targ = Axiom::Expr->new({
            type => 'forall',
            args => [ $_->copy, $targ ],
        }) for reverse @{ $self->var };
        $targ;
    } : undef;
}

sub name_lemma {
    my($self) = @_;
    my $known = $self->lemmas;
    my $base = $self->name;
    my $i = 1;
    ++$i while $known->{"$base$i"};
    return "$base$i";
}

sub make_lemma {
    my($self, $targ) = @_;
    my $name = $self->name_lemma;
    my $req = $self->req;
    return '' unless @$req;
    my $arneed = $self->arneed;
    my $label = 'a';
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

    my $lcond = $name . $label++;
    $line->('condstart', undef, $lcond, $self->cond);

    my $index = 0;
    my $conjlege;
    for my $ri (0 .. $#$req) {
        my($lege, $neg, $pos) = @{ $req->[$ri] };
        if ($ri >= 1) {
            # the end of the previous chunk is the line to conjoin, so
            # label it and remember the label and the expr
            $conjlege = [ $name . $label++, $line[-1]->expr ];
            $line[-1]->set_label($conjlege->[0]);
        }
        my $conjsign;
        for ([ -1, $neg ], [ 1, $pos ]) {
            my($sign, $subreqs) = @$_;
            next unless @$subreqs;
            my $eff = ($sign < 0) ? (1 - $lege) : $lege;
            if (@$subreqs == 1) {
                my $args = ($index++) ? [ $lcond ] : undef;
                my($i, $c) = @{ $subreqs->[0] };
                my $cond = $self->icond($i, $eff);
                $line->('simplify', $args, undef, $cond);
            } else {
                my $args = ($index++) ? [ $lcond ] : undef;
                my $mpargs = [ ($eff) ? 'int.Q8' : 'int.Q7' ];
                my($i1, $c1) = @{ $subreqs->[0] };
                my $cond1 = $self->icond($i1, $eff);
                my($i2, $c2) = @{ $subreqs->[1] };
                my $cond2 = $self->icond($i2, $eff);
                my $cond = Axiom::Expr->new({
                    type => 'andlist',
                    args => [ $cond1, $cond2 ],
                });
                $line->('simplify', $args, undef, $cond);
                my $mpsub = sub {
                    my($c1, $c2) = @_;
                    my($r1, $r2) = map Axiom::Expr->new({
                        type => 'pluslist',
                        args => [
                            $c1->args->[$_]->copy,
                            $c2->args->[$_]->copy,
                        ],
                    })->clean, (0, 1);
                    return Axiom::Expr->new({
                        type => $cond1->type,
                        args => [ $r1, $r2 ],
                    });
                };
                my $mp = $mpsub->($cond1, $cond2);
                $line->('ponens', $mpargs, undef, $mp);
                for (2 .. $#$subreqs) {
                    my($i, $c) = @{ $subreqs->[$_] };
                    my $cond = $self->icond($i, $eff);
                    my $conjmore = $name . $label++;
                    $line[-1]->set_label($conjmore);
                    $line->('simplify', [ $lcond ], undef, $cond);
                    my $ce = Axiom::Expr->new({
                        type => 'andlist',
                        args => [ $mp->copy, $cond->copy ],
                    });
                    $line->('conjoin', [ $conjmore ], undef, $ce);
                    my $mp2 = $mpsub->($mp, $cond);
                    $line->('ponens', $mpargs, undef, $mp2);
                    $mp = $mp2;
                }
            }
            if ($sign < 0) {
                my $prev = $line[-1]->expr;
                my $expr = Axiom::Expr->new({
                    type => $prev->inverse_type,
                    args => [ map $_->negate, @{ $prev->args } ],
                })->clean;
                my $refer;
                if (@$pos) {
                    # there is + to follow
                    $refer = $name . $label++;
                    $conjsign = [ $expr, $refer ];
                }
                $line->('negate', undef, $refer, $expr);
            }
        }
        if ($conjsign) {
            my($neg, $refer) = @$conjsign;
            my $pos = $line[-1]->expr;
            my $expr = Axiom::Expr->new({
                type => 'andlist',
                args => [ $neg->copy, $pos->copy ],
            });
            $line->('conjoin', [ $refer ], undef, $expr);
            my $args = [ $lege ? 'int.Q8' : 'int.Q7' ];
            my($left, $right) = map Axiom::Expr->new({
                type => 'pluslist',
                args => [ $neg->args->[$_]->copy, $pos->args->[$_]->copy ],
            })->clean, (0, 1);
            $expr = Axiom::Expr->new({
                type => ($lege ? 'rle' : 'rge'),
                args => [ $left, $right ],
            });
            $line->('ponens', $args, undef, $expr);
        }
        my $type = $line[-1]->expr->type;
        my $thistarg = List::Util::first(sub { $_->type eq $type }, @$arneed);
        $line->('add', undef, undef, $thistarg->copy)
                if $thistarg->diff($line[-1]->expr, 1, 2);
        if ($ri >= 1) {
            my($label, $expr1) = @$conjlege;
            my $expr2 = $line[-1]->expr;
            $line->('conjoin', [ $label ], undef, Axiom::Expr->new({
                type => 'andlist',
                args => [ $expr1->copy, $expr2->copy ],
            }));
        }
    }
    $line->('condend', undef, undef, $targ->copy);
    $line->('theorem', undef, $name, $targ->copy);

    $self->lemmas->{$name} = \@line;
#Prove::Line->emit_set(\@line);
    return $name;
}

1;
