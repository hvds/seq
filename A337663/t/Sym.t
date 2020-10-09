use strict;
use warnings;
use Test::More;

use lib 'lib';
use Sym;

my @all = Sym->all;
is(8, 0 + @all, "there are 8 syms");

is(join('-', sort Sym->all_bits(0xff)), join('-', sort @all),
        "all_bits(0xff) finds them all");

# make a vals arrayref of arrayrefs from a string
sub _v {
    my($str) = @_;
    return [ map [ split /\s+/ ], split /;\s*/, $str ];
}

# make a string from a vals arrayref
sub _u {
    my($vals) = @_;
    return join '; ', map join(' ', @$_), @$vals;
}

for (
    [ 'Sym::xy', '1 2; 3 4' ],
    [ 'Sym::xy', '1 3 5 7' ],
    [ 'Sym::xy', '1; 3; 5; 7' ],
    [ 'Sym::xY', '1; 2; 3; 4' ],
    [ 'Sym::xY', '1 1; 2 2; 3 3' ],
    [ 'Sym::xY', '1 3 1; 2 4 2' ],
    [ 'Sym::Xy', '1 2 3 4' ],
    [ 'Sym::Xy', '1 2 3; 1 2 3' ],
    [ 'Sym::Xy', '1 2; 3 4; 1 2' ],
    [ 'Sym::XY', '1 1 0; 0 1 1' ],
    [ 'Sym::XY', '1 1 0; 0 0 0; 0 1 1' ],
    [ 'Sym::XY', '1 0; 1 1; 0 1' ],
    [ 'Sym::yx', '1 2; 2 3' ],
    [ 'Sym::yx', '1 3 5; 3 9 13; 5 13 21' ],
    [ 'Sym::yX', '1 1; 1 1' ],
    [ 'Sym::yX', '1 0 1; 0 2 0; 1 0 1' ],
    [ 'Sym::Yx', '1 1; 1 1' ],
    [ 'Sym::Yx', '1 0 1; 0 2 0; 1 0 1' ],
    [ 'Sym::YX', '0 1; 1 0' ],
    [ 'Sym::YX', '0 1 1; 3 0 1; 5 3 0' ],
) {
    my($class, $str) = @$_;
    my $v = _v($str);
    is($class->check($v), 1, "$class check $str ok");
    is(_u($class->transform($v)), $str, "$class transform $str ok");
}

for (
    [ 'Sym::xY', '1 2; 3 4', '2 1; 4 3' ],
    [ 'Sym::xY', '1 2 3 4', '4 3 2 1' ],
    [ 'Sym::Xy', '1 2; 3 4', '3 4; 1 2' ],
    [ 'Sym::Xy', '1; 2; 3; 4', '4; 3; 2; 1' ],
    [ 'Sym::XY', '1 2; 3 4', '4 3; 2 1' ],
    [ 'Sym::XY', '1; 2; 3; 4', '4; 3; 2; 1' ],
    [ 'Sym::XY', '1 2 3 4', '4 3 2 1' ],
    [ 'Sym::yx', '1 2; 3 4', '1 3; 2 4' ],
    [ 'Sym::yx', '1; 2; 3; 4', '1 2 3 4' ],
    [ 'Sym::yx', '1 2 3 4', '1; 2; 3; 4' ],
    [ 'Sym::yX', '1 2; 3 4', '3 1; 4 2' ],
    [ 'Sym::yX', '1; 2; 3; 4', '4 3 2 1' ],
    [ 'Sym::yX', '1 2 3 4', '1; 2; 3; 4' ],
    [ 'Sym::Yx', '1 2; 3 4', '2 4; 1 3' ],
    [ 'Sym::Yx', '1; 2; 3; 4', '1 2 3 4' ],
    [ 'Sym::Yx', '1 2 3 4', '4; 3; 2; 1' ],
    [ 'Sym::YX', '1 2; 3 4', '4 2; 3 1' ],
    [ 'Sym::YX', '1; 2; 3; 4', '4 3 2 1' ],
    [ 'Sym::YX', '1 2 3 4', '4; 3; 2; 1' ],
) {
    my($class, $str, $trans) = @$_;
    my $v = _v($str);
    is($class->check($v), 0, "$class check $str ok");
    is(_u($class->transform($v)), $trans, "$class transform $str ok");
}

for (
    [ 'Sym::xy', 3, 3, 0, 0, 0, 0 ],
    [ 'Sym::xy', 9, 1, 2, 2, 2, 2 ],
    [ 'Sym::xy', 9, 1, -2, -2, -2, -2 ],
    [ 'Sym::xY', 3, 3, 0, 0, 0, 2 ],
    [ 'Sym::xY', 9, 1, 2, 2, 2, -2 ],
    [ 'Sym::xY', 9, 1, -2, -2, -2, 2 ],
    [ 'Sym::Xy', 3, 3, 0, 0, 2, 0 ],
    [ 'Sym::Xy', 9, 1, 2, 2, 6, 2 ],
    [ 'Sym::Xy', 9, 1, -2, -2, 10, -2 ],
    [ 'Sym::XY', 3, 3, 0, 0, 2, 2 ],
    [ 'Sym::XY', 9, 1, 2, 2, 6, -2 ],
    [ 'Sym::XY', 9, 1, -2, -2, 10, 2 ],
    [ 'Sym::yX', 3, 3, 0, 0, 0, 2 ],
    [ 'Sym::yX', 9, 1, 2, 2, 2, 6 ],
    [ 'Sym::yX', 9, 1, -2, -2, -2, 10 ],
    [ 'Sym::Yx', 3, 3, 0, 0, 2, 0 ],
    [ 'Sym::Yx', 9, 1, 2, 2, -2, 2 ],
    [ 'Sym::Yx', 9, 1, -2, -2, 2, -2 ],
    [ 'Sym::YX', 3, 3, 0, 0, 2, 2 ],
    [ 'Sym::YX', 9, 1, 2, 2, -2, 6 ],
    [ 'Sym::YX', 9, 1, -2, -2, 2, 10 ],
) {
    my($class, $x, $y, $lx, $ly, $ex, $ey) = @$_;
    my $g = G->new($x, $y);
    my($tx, $ty) = @{ $class->transform_loc($g, [ $lx, $ly ]) };
    is("$tx.$ty", "$ex.$ey", "$class transform_loc ok");
}

done_testing();

# transform_loc needs something Group-like enough to respond to x and y.
package G {
    sub new { bless [ @_[1, 2] ], $_[0] }
    sub x { shift->[0] }
    sub y { shift->[1] }
};
