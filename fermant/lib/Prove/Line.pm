package Prove::Line;
use strict;
use warnings;

=head1 NAME

Prove::Line - represent a line of a proof

=cut

sub new {
    my($class, $data) = @_;
    return bless $data, $class;
}

for my $attr (qw{ type args label expr }) {
    my $method = sub {
        my($self) = @_;
        return $self->{$attr};
    };
    no strict 'refs';
    *$attr = $method;
}

sub set_label {
    my($self, $new) = @_;
    $self->{label} = $new;
    return $new;
}

sub emit {
    my($self, $indent) = @_;
    my $e = $self->expr->str;
    # FIXME: horrible hack until Axiom::Expr->str works better
    $e =~ s{\((-?\d*(?:[a-z]+(?:\^\d+|\{\d+\})?)*)\)}{$1}g;
    $e =~ s{\+-}{-}g;
    $e =~ s{([-+])1([a-z])}{$1$2}g;
    printf "%s%s%s%s: %s\n",
        ' ' x $indent,
        $self->type,
        $self->textargs,
        $self->textlabel,
        $e;
    return;
}

sub textargs {
    my $args = shift->args;
    return '' unless $args && @$args;
    return sprintf '(%s)', join ', ', @$args;
}

sub textlabel {
    my $label = shift->label;
    return '' unless defined($label) && length($label);
    return " $label";
}

sub emit_set {
    my($class, $set, $indent) = @_;
    $indent //= 2;
    for my $line (@$set) {
        my $type = $line->type;
        my($pre, $post) = ($type eq 'condstart') ? (0, 2)
            : ($type eq 'condend') ? (-2, 0)
            : ($type eq 'theorem' || $type eq 'axiom') ? (-2, 0)
            : (0, 0);
        $indent += $pre;
        $line->emit($indent);
        $indent += $post;
    }
    print "\n";
    return;
}

1;
