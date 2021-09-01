package ModFunc;
use strict;

use Math::GMP ();

require Exporter;
our @ISA = qw/ Exporter /;
our @EXPORT_OK = qw/
    quadres is_residue quadvec
    mod_combine
    gcd
/;

{
    my %q;

    #
    # Returns Legendre ( n / m ), which is 0 if gcd(n, m) > 1,
    # 1 if n is a quadratic residue (mod m), -1 otherwise.
    #
    sub quadres {
        my($n, $m) = @_;
        return 0 if gcd($n, $m) > 1;
        my $v = quadvec($m);
        return vec($v, ($n % $m), 1) ? 1 : -1;
    }

    #
    # Returns true if Legendre ( n / m ) is 0 or 1.
    #
    sub is_residue {
        my($n, $m) = @_;
        my $v = quadvec($m);
        return vec($v, ($n % $m), 1);
    }

    #
    # Returns an m-bit vector, such that bit i is true iff is_residue(i, m).
    #
    sub quadvec {
        my($m) = @_;
        ${ $q{$m} ||= do {
            my $v = "";
            vec($v, (($_ * $_) % $m), 1) = 1 for 0 .. $m - 1;
            \$v;
        } };
    }
}


#
# Given constraints (a, b) and (c, d) representing n=a(mod b), n=c(mod d),
# return (e, f) such that n=e(mod f) encapsulates both constraints.
#
sub mod_combine {
    my($a, $b, $c, $d) = @_;
    my $g = gcd($b, $d);
    if ($g > 1) {
        die <<DIE unless ($a % $g) == ($c % $g);
Inconsistent values: can't combine $a(mod $b) with $c(mod $d)
DIE
        while ($g > 1) {
            my $gb = gcd($b / $g, $g);
            my $gd = gcd($d / $g, $g);
            if ($gb > 1) {
                # p|gb constrained to higher power in (a, b) than in (c, d)
                $d /= $gb;
                $c %= $d;
            } elsif ($gd > 1) {
                # p|gd constrained to higher power in (c, d) than in (a, b)
                $b /= $gd;
                $a %= $b;
            } else {
                # p|g constrained to same power in both; pick one
                $b /= $g;
                $a %= $b;
            }
            $g = gcd($b, $d);
        }
    }
    # with the gcd reduced, we can safely take inverses
    my $bi = $b->bmodinv($d);
    my $di = $d->bmodinv($b);
    my $ab = (($a * $di) % $b) * $d;
    my $cd = (($c * $bi) % $d) * $b;
    my $bd = $b * $d;
    return (($ab + $cd) % $bd, $bd);
}

sub gcd {
    my($a, $b) = @_;
    $a = Math::GMP->new($a) unless ref $a;
    return $a->bgcd($b);
}

1;
