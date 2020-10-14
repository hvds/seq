use strict;
use warnings;
use Test::More;

use lib 'lib';
use Sym;

my @all = Sym->all;
is(8, 0 + @all, "there are 8 syms");

is(join('-', sort Sym->all_bits(0xff)), join('-', sort @all),
        "all_bits(0xff) finds them all");

for my $s (@all) {
    my $t = $s->is_transpose;
    my $match = ($s =~ /::xy/i) ? 0 : ($s =~ /::yx/i) ? 1
            : die "unexpected Sym class '$s'";
    is(!$t, !$match, "$s is_transpose");
}

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

# make a class name from a short name
sub _c {
    my($short) = @_;
    return "Sym::$short";
}

# make a dimension string from a vals arrayref
sub _d {
    my($vals) = @_;
    return join('-', 0 + @$vals, 0 + @{ $vals->[0] });
}

open(my $f, '<', 't/test_sym.h')
        or die "t/test_sym.h: $!";
my %data;
while (<$f>) {
    my($type, $data) = m{
        ^ \s* \(
            test_(\w+)_t
        \)\{ \s*
            (.*?)
        \s* \} ,? $
    }x or next;
    my @args = map s{^"(.*)"$}{$1}r, split /,\s*/, $data;
    push @{ $data{$type} }, \@args;
}
close $f;
for my $type (qw{ sym asym loc }) {
    my $tests = delete $data{$type}
            // die "Failed to find tests of type '$type'";
    if ($type eq 'sym') {
        for (@$tests) {
            my($s, $x, $y, $d) = @$_;
            my $class = _c($s);
            my $vals = _v($d);
            _d($vals) eq "$x-$y" or die "unwrapped '$d'";
            is($class->check($vals), 1, "$class check $d");
            my $trans = $class->transform($vals);
            _d($trans) eq "$x-$y" or die "$class transform $d";
            is(_u($trans), $d, "$class transform $d");
        }
    } elsif ($type eq 'asym') {
        for (@$tests) {
            my($s, $x, $y, $d, $e) = @$_;
            my $class = _c($s);
            my $vals = _v($d);
            _d($vals) eq "$x-$y" or die "unwrapped '$d'";
            is($class->check($vals), 0, "$class check $d");
            my $trans = $class->transform($vals);
            _d($trans) eq ($class->is_transpose ? "$y-$x" : "$x-$y")
                    or die "$class transform $d dims";
            is(_u($trans), $e, "$class transform $d");
        }
    } elsif ($type eq 'loc') {
        for (@$tests) {
            my($s, $x, $y, $lx, $ly, $ex, $ey) = @$_;
            my $class = _c($s);
            my $same = ($lx == $ex && $ly == $ey);
            my($tx, $ty) = @{ $class->transform_loc($x, $y, [ $lx, $ly ]) };
            is("$tx-$ty", "$ex-$ey", "$class transform_loc $lx.$ly");
            is($class->check_loc($x, $y, [ $lx, $ly ]), $same,
                    "$class check_loc $lx.$ly");
        }
    }
}
if (keys %data) {
    die "Unexpected test data types: @{[ join ', ', keys %data ]}\n";
}

done_testing();
