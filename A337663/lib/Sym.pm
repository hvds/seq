package Sym;
use strict;
use warnings;

use List::Util;

my %by_bits = map +($_->bit => $_), Sym->all;
my @ordered_bits = map $_->bit, Sym->all;

sub all {
    return qw{
        Sym::xy Sym::xY Sym::Xy Sym::XY
        Sym::yx Sym::yX Sym::Yx Sym::YX
    };
}

sub all_bits {
    my($class, $bits) = @_;
    return map $by_bits{$_}, grep !$_ || ($bits & $_), @ordered_bits;
}

# null symmetry
package Sym::xy {
    sub bit { 0 }
    sub is_transpose { 0 }

    sub check {
        my($class, $vals) = @_;
        return 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ @$loc ];
    }

    sub transform {
        my($class, $vals) = @_;
        return [ map [ @$_ ], @$vals ];
    }
};

# x -> x, y -> -y
package Sym::xY {
    sub bit { 1 }
    sub is_transpose { 0 }

    sub check {
        my($class, $vals) = @_;
        return 0 + List::Util::all sub {
            join('.', @$_) eq join('.', reverse @$_)
        }, @$vals;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $loc->[0], $group->y - 1 - $loc->[1] ];
    }

    sub transform {
        my($class, $vals) = @_;
        return [ map [ reverse @$_ ], @$vals ];
    }
};

# x -> -x, y -> y
package Sym::Xy {
    sub bit { 2 }
    sub is_transpose { 0 }

    sub check {
        my($class, $vals) = @_;
        my $x = $#$vals;
        return 0 + List::Util::all sub {
            join('.', @{ $vals->[$_] }) eq join('.', @{ $vals->[$x - $_] })
        }, 0 .. ($x + 1) >> 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $group->x - 1 - $loc->[0], $loc->[1] ];
    }

    sub transform {
        my($class, $vals) = @_;
        return [ map [ @$_ ], reverse @$vals ];
    }
};

# x -> -x, y -> -y
package Sym::XY {
    sub bit { 4 }
    sub is_transpose { 0 }

    sub check {
        my($class, $vals) = @_;
        my $x = $#$vals;
        return 0 + List::Util::all sub {
            join('.', @{ $vals->[$_] })
                    eq join('.', reverse @{ $vals->[$x - $_] })
        }, 0 .. $x;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $group->x - 1 - $loc->[0], $group->y - 1 - $loc->[1] ];
    }

    sub transform {
        my($class, $vals) = @_;
        return [ map [ reverse @$_ ], reverse @$vals ];
    }
};

# x -> y, y -> x
package Sym::yx {
    sub bit { 8 }
    sub is_transpose { 1 }

    sub check {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return 0 unless $x == $y;
        for my $i (0 .. $x) {
            for my $j ($i + 1 .. $y) {
                return 0 unless $vals->[$i][$j] == $vals->[$j][$i];
            }
        }
        return 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $loc->[1], $loc->[0] ];
    }

    sub transform {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return [ map {
            my $i = $_; [ map $vals->[$_][$i], 0 .. $x ]
        } 0 .. $y ];
    }
};

# x -> y, y -> -x
package Sym::yX {
    sub bit { 16 }
    sub is_transpose { 1 }

    sub check {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return 0 unless $x == $y;
        for my $i (0 .. $x) {
            for my $j (0 .. $y) {
                return 0 unless $vals->[$i][$j] == $vals->[$j][$x - $i];
            }
        }
        return 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $loc->[1], $group->x - 1 - $loc->[0] ];
    }

    sub transform {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return [ map {
            my $i = $_; [ map $vals->[$_][$i], reverse 0 .. $x ]
        } 0 .. $y ];
    }
};

# x -> -y, y -> x
package Sym::Yx {
    sub bit { 32 }
    sub is_transpose { 1 }

    sub check {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return 0 unless $x == $y;
        for my $i (0 .. $x) {
            for my $j (0 .. $y) {
                return 0 unless $vals->[$i][$j] == $vals->[$y - $j][$i];
            }
        }
        return 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $group->y - 1 - $loc->[1], $loc->[0] ];
    }

    sub transform {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return [ map {
            my $i = $_; [ map $vals->[$_][$i], 0 .. $x ]
        } reverse 0 .. $y ];
    }
};

# x -> -y, y -> -x
package Sym::YX {
    sub bit { 64 }
    sub is_transpose { 1 }

    sub check {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return 0 unless $x == $y;
        for my $i (0 .. $x) {
            for my $j (0 .. $y) {
                return 0 unless $vals->[$i][$j] == $vals->[$y - $j][$x - $i];
            }
        }
        return 1;
    }

    sub transform_loc {
        my($class, $group, $loc) = @_;
        return [ $group->y - 1 - $loc->[1], $group->x - 1 - $loc->[0] ];
    }

    sub transform {
        my($class, $vals) = @_;
        my($x, $y) = ($#$vals, $#{ $vals->[0] });
        return [ map {
            my $i = $_; [ map $vals->[$_][$i], reverse 0 .. $x ]
        } reverse 0 .. $y ];
    }
};

1;
