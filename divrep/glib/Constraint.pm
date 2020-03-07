package Constraint;
use strict;
use ModFunc qw/ quadvec mod_combine /;
use Math::Prime::Util qw();
use warnings;
no warnings qw/ recursion /;

sub MBI { return Math::GMP->new(@_) }

my $debug = 0;
my $BIT32 = 2 * (1 << 31);

=head

Usage:
  my $c = Constraint->new(
    'n' => $n, 'f' => $f,            # we are searching for A165500(n)==f
    'min' => $min, 'max' => $max,    # setting d in [min, max]
    'tell_count' => $tell_count, 't0' => $t0, 'tau' => $tau,
  );

  # first apply all constraints
  $c->suppress(3, 1);            # disallow d==1 (mod 3)
  $c->suppress(3, 1, 27);        # disallow d==1 (mod 3) for d>27
  $c->require(5, 0);             # require d==0 (mod 5)
  $c->require(5, 0, 25);         # require d==0 (mod 5) for d>25
  $c->series(13, 169);           # disallow n+kd==0 (mod 13) for 0<k<f, d>169

  # then initialize
  $c->init;

  # then use
  while (my $d = $c->next) {
    # use $d
  }
  # done

=cut

sub debug_more { ++$debug }
sub debug { $debug }

sub new {
    my($class) = shift;
    my $self = bless {
        c => {},            # hash lookup of constraints by modulus
        allc => [],         # list of all constraints sorted by modulus
        sc => [],           # list of active constraints sorted by ?density?
        pend => [],         # list of not-yet-active constraints
        mult => MBI(1),     # d = mod_mult (mod mult) is required by constraints
        mod_mult => MBI(0),
        tests => 0,
        skipped => 0,
        kept => 0,
        max => undef,       # set in init()
        cur => undef,       # set in init()
        @_,                 # n, f, tell_count, t0, min, max, check, tau
    }, $class;
    $_ = MBI($_ || 0) for @$self{qw/ n f min max tau /};
    return $self;
}

for my $method (qw/
    cur n f tell_count t0 min max check tau min_potency mult mod_mult
/) {
    my $sub = sub { shift->{$method} };
    no strict 'refs';
    *$method = $sub;
}

sub elapsed {
    my $self = shift;
    my $t1 = times();
    return $t1 - $self->{'t0'};
}

sub c {
    my($self, $n) = @_;
    $self->{'c'}{$n} ||= do {
        my $nb = ref($n) ? $n : MBI($n);
        my $div = _highfactors($nb);    # x: x>1, px=n
        my $c = [
            $nb,    # [0] the modulus
            '',     # [1] bit vector of uniquely disallowed values
            '',     # [2] bit vector of disallowed values
            0,      # [3] count of uniquely disallowed values
            0,      # [4] TRUE when incorporated in mult/mod_mult
            0,      # [5] TRUE if new disallowed since last recalc
            $div,   # [6] d: d | p
            [],     # [7] d: p | d
        ];
        # apply dependencies for my divisors
        for my $d (@$div) {
            my $cd =  $self->c($d);
            push @{ $cd->[7] }, $n;
            my $q = $n / $d;
            for my $v (0 .. $d - 1) {
                next unless vec($cd->[2], $v, 1);
                for my $m (0 .. $q - 1) {
                    vec($c->[2], $m * $d + $v, 1) = 1;
                    $c->[5] = 1
                }
            }
        }
        ($debug > 2) && warn "init $n with vec [@{[ unpack 'b*', $c->[2] ]}]\n";
        _insert($self->{'allc'}, $c);
        $c;
    };
}

sub init {
    my($self) = @_;
    my $min = $self->min();

    # We may get a very large number of pending values: just recording the
    # least (and maybe storing the rest in a different structure) may be
    # preferable
    @{ $self->{pend} } = sort { $a->[0] <=> $b->[0] } @{ $self->{pend} };

    # apply any pending constraints that are immediately active
    # and force a recalc
    $self->{'cur'} = $self->catchup($min, 1) - $self->{'mult'};

    ($debug > 3) && warn "inited: ", $self->Dump();
}

sub recalc {
    my($self, $cur) = @_;
    while (1) {
        #
        # find the current list of active constraints
        #
        $self->find_active();

        #
        # that may change the fixed mod - if $cur is still valid, we're done
        #
        my $diff = _diff_to_mod($cur, @$self{qw/ mod_mult mult /}) or last;
        $cur += $diff;
($diff > 0 && $debug > 2) && warn "recalc: bump cur by $diff to $cur\n";

        #
        # increasing $cur may mean new pending constraints get activated.
        # if not, we're done, else we must re-check what's active.
        #
        $self->catchup($cur) or last;
    }
    $self->pack_sc();
    return $cur;
}

#
# Find all active constraints. Moduli permitting only a single value are
# combined into <mod_mult, mult>; constraints for other moduli are included
# in the list $self->{'sc'}. Optimisation of the list is done as a separate
# step.
#
sub find_active {
    my $self = shift;
    my @sc;
    my($mod_mult, $mult) = @$self{qw/ mod_mult mult /};
    for my $c (@{ $self->{'allc'} }) {
        # must be in order to ensure p is incorporated before p^2
        my $d = $c->[0];
        next if $c->[4];    # already incorporated in mult
        if ($c->[5]) {
            # touched since we last looked
            $c->[5] = 0;
            if (unpack('%32b*', $c->[2]) >= $d - 1) {
                my $v = (grep !vec($c->[2], $_, 1), 0 .. $d - 1)[0];
                unless (defined $v) {
                    printf <<OUT, $d, $self->elapsed;
402 Error: all values (mod %s) disallowed (%.3fs)
OUT
                    exit 1;
                }
                ($debug > 2)
                        && warn "fix $v(mod $d) in $mod_mult(mod $mult)\n";
                ($mod_mult, $mult) = mod_combine($mod_mult, $mult, $v, $d);
                ($debug > 1) && warn "now fixed: $mod_mult(mod $mult)\n";
                @$self{qw/ mod_mult mult /} = ($mod_mult, $mult);
                $c->[4] = 1;
                next;
            }
        }
        push @sc, $c if $c->[3];
    }
    $self->{'sc'} = \@sc;
    printf "302 checking %s(mod %s): %s\n",
            $mod_mult, $mult, join ' ', map "$_->[3]/$_->[0]", @sc;
}

#
# Average number of values tested against modular constraints before
# a value gets through to the full testers regime.
#
sub frequency {
    my($self) = @_;
    return $self->{freq} //= do {
        $self->find_active unless $self->{sc};
        my $f = 1;
        for (@{ $self->{sc} }) {
            my $m = numify($_->[0]);
            my $d = numify($_->[3]);
            $f *= $m / ($m - $d);
        }
        $f;
    };
}

#
# Optimise the list of active constraints: where practical, combine multiple
# moduli into a single test; then sort the remaining list in order of potency.
#
# TODO: further combine active constraints into inactive moduli as long as
# the target modulus is within range.
#
sub pack_sc {
    my $self = shift;
    # Cache this before packing, since we don't maintain the unique values
    # count in the packed list.
    $self->frequency;
    my @aux = map +{
        'sc' => $_,
        'n' => $_->[0],
        'nn' => numify($_->[0]),
        'ncheck' => $_->[3],
        'vcheck' => $_->[1],
    }, @{ $self->{'sc'} };

    # where n_i divides n_j, merge a_i into a_j
    AI: for (my $i = 0; $i < @aux; ++$i) {
        my $ai = $aux[$i];
        my $ni = $ai->{'nn'};
        for my $j ($i + 1 .. $#aux) {
            my $aj = $aux[$j];
            next if ($aj->{'nn'} % $ni) > 0;
            # will merge
            splice @aux, $i, 1;
            _aux_merge($aj, $ai);
            redo AI;
        }
    }

    # Sort by potency ((number of uniquely disallowed values (mod m)) / m).
    $_->{'potency'} = numify($_->{'ncheck'}) / $_->{'nn'} for @aux;
    @aux = sort { $b->{'potency'} <=> $a->{'potency'} } @aux;

    # where lcm(n_i, n_j) is within the requested range, merge a_j into a_i
    my $maxsize = $self->check;
    AI2: for (my $i = 0; $i < @aux; ++$i) {
        my $ai = $aux[$i];
        my $ni = $ai->{'n'};
        for my $j ($i + 1 .. $#aux) {
            my $aj = $aux[$j];
            my $nij = $ni->blcm($aj->{'n'});
            next if $nij > $maxsize;
            splice(@aux, $j, 1);
            $ai->{'sc'} = $self->c($nij);
            _aux_merge2($ai, $aj, $nij);
            $ai->{'potency'} = numify($ai->{'ncheck'}) / numify($nij);
            redo AI2;
        }
    }

    # Sort by potency again
    @aux = sort { $b->{'potency'} <=> $a->{'potency'} } @aux;

    # now suppress anything insufficiently potent
    my $min = $self->min_potency;
    if ($min) {
        for (my $i = 0; $i < @aux; ++$i) {
            next if $aux[$i]->{'potency'} * $min >= 1;
            $debug >= 1 and warn "Splice sc[$i..$#aux] for min potency $min\n";
            splice @aux, $i;
            last;
        }
    }

    printf "303 packed %s(mod %s): %s\n",
            @$self{qw/ mod_mult mult /},
            join ' ', map "$_->{'ncheck'}/$_->{'n'}", @aux;
    $self->{'sc'} = [ map $_->{'sc'}, @aux ];
}

sub _aux_merge {
    my($dest, $src) = @_;
    my($nd, $ns) = map $_->{'n'}, ($dest, $src);
    my $q = $nd / $ns;
    for my $x (0 .. $ns - 1) {
        next unless vec($src->{'vcheck'}, $x, 1);
        for my $y (0 .. $q - 1) {
            my $z = $y * $ns + $x;
            next if vec($dest->{'vcheck'}, $z, 1);
            vec($dest->{'vcheck'}, $z, 1) = 1;
            ++$dest->{'ncheck'};
        }
    }
}

sub _aux_merge2 {
    my($dest, $src, $size) = @_;
    my $nd = $dest->{'n'};
    my $gcd = $size / $nd;
    my $vd = $dest->{'vcheck'};
    for my $x (0 .. $nd - 1) {
        next unless vec($vd, $x, 1);
        for my $y (1 .. $gcd - 1) {
            vec($vd, $y * $nd + $x, 1) = 1;
        }
    }
    $dest->{'n'} = $size;
    $dest->{'ncheck'} *= $gcd;
    $dest->{'vcheck'} = $vd;
    _aux_merge($dest, $src);
}

sub catchup {
    my($self, $cur, $recalc) = @_;
    my $pend = $self->{pend};
    while (@$pend && $pend->[0][0] < $cur) {
        my($effective, $action) = @{ shift @$pend };
($debug > 2) && warn "catchup: trigger pend($effective) at cur=$cur\n";
        $recalc |= $action->();
    }
    return $recalc ? $self->recalc($cur) : $cur;
}

sub Dump {
    my $self = { %{ +shift } };
    defined($self->{$_}) && ($self->{$_} =
        "$self->{$_}"
    ) for qw/
        mult mod_mult tests skipped kept max cur check f min n tau tell_count
    /;
    defined($self->{$_}) && ($self->{$_} = [
        map "c_$_->[0]", @{ $self->{$_} }
    ]) for qw/ sc allc /;
    defined($self->{$_}) && ($self->{$_} = [
        map [ map "$_", @$_ ], @{ $self ->{$_} }
    ]) for qw/ pend /;
    defined($self->{$_}) && ($self->{$_} = do {
        my $s = $self->{$_};
        my $d = {};
        for (keys %$s) {
            $d->{$_} = [
                "$s->{$_}[0]",                # p
                unpack("b*", $s->{$_}[1]),    # v1
                unpack("b*", $s->{$_}[2]),    # v2
                "$s->{$_}[3]",                # unique
                "$s->{$_}[4]",                # fixed
                "$s->{$_}[5]",                # changed
            ];
        }
        $d;
    }) for qw/ c /;
    use Data::Dumper;
    local $Data::Dumper::Indent=1;
    local $Data::Dumper::Sortkeys = sub {
        my $h = shift;
        my @keys = keys %$h;
        defined($keys[0]) && ($keys[0] =~ /^\d+$/)
            ? [ sort { $a <=> $b } @keys ]
            : [ sort { $a cmp $b } @keys ]
    };
    return Dumper($self);
}

use Inline C => <<'INLINE', NAME => 'Constraint'; no AutoLoader;
SV* _av_fetch_actual(AV* av, I32 off) {
    SV** svp = av_fetch(av, off, FALSE);
    if (!svp || !*svp) {
        return &PL_sv_undef;
    } else {
        return *svp;
    }
}

int _bool(SV* sv) {
    int v = SvOK(sv) ? 1 : 0;
    SvREFCNT_dec(sv);
    return v;
}

#define MAXUINT 0xffffffff

typedef struct scmod_t {
    SV* svmod;
    UV mod;
    int veclen;
    char* vec;
    UV basemod;
    UV multmod;
} scmod;

typedef struct sclist_t {
    int len;
    scmod* sc;
    SV* base;
    SV* mult;
    UV iter;
    UV limit;
} sclist;

sclist* scl = (sclist *)NULL; /* values extracted from arrayref sc (malloced) */

void _free_sc() {
    if (scl) {
        if (scl->sc) {
            int i;
            for (i = 0; i < scl->len; ++i) {
                if (scl->sc[i].mod) SvREFCNT_dec(scl->sc[i].svmod);
            }
            free(scl->sc);
        }
        if (scl->base) SvREFCNT_dec(scl->base);
        if (scl->mult) SvREFCNT_dec(scl->mult);
        free(scl);
    }
}

void _init_sc(SV* cur, SV* mult, AV* sc) {
    /* sc is an arrayref of [ modulus, any, bitvector, ... ] sets, sorted
       in order of descending usefulness. Any value of cur should be skipped
       if the (cur % modulus)'th bit of the bitvector is set for any of the
       listed constraints.
     */
    SV* entry;
    AV* aentry;
    SV* vec;
    I32 i;
    UV limit;
    UV mod;

    SAVETMPS;
    _free_sc();
    scl = (sclist *)malloc(sizeof(sclist));
    scl->base = newSVsv(cur);
    scl->mult = mult;
    SvREFCNT_inc(mult);
    scl->iter = 0;
    scl->limit = MAXUINT;

    scl->len = (sc) ? av_len(sc) + 1 : 0;
    if (scl->len == 0) {
        scl->sc = (scmod *)NULL;
    } else {
        scl->sc = (scmod *)malloc(scl->len * sizeof(scmod));
    }
    for (i = 0; i < scl->len; ++i) {
        scmod* scm = &scl->sc[i];
        entry = _av_fetch_actual(sc, i);
        if (!SvROK(entry)) croak("Invalid entry in sc[%d]\n", i);
        aentry = (AV*)SvRV(entry);
        scm->mod = mod = sv_2uv(_av_fetch_actual(aentry, 0));
        scm->svmod = newSVuv(mod);
        vec = _av_fetch_actual(aentry, 2);
        scm->vec = SvPVX(vec);
        scm->veclen = SvCUR(vec);
        scm->basemod = sv_2uv(amagic_call(cur, scm->svmod, modulo_amg, 0));
        scm->multmod = sv_2uv(amagic_call(mult, scm->svmod, modulo_amg, 0));
        /* more precisely, limit = (MAXUINT - basemod) / multmod, but this
         * slightly looser calculation saves recalculating the limit when
         * we reset the base
         */
        limit = (MAXUINT - mod) / scm->multmod;
        if (scl->limit > limit) scl->limit = limit;
    }
    FREETMPS;
}

SV* _calc_result() {
    SV* sviter = newSVuv(scl->iter);
    SV* prod = amagic_call(scl->mult, sviter, mult_amg, 0);
    SV* result = amagic_call(scl->base, prod, add_amg, 0);
    SvREFCNT_dec(sviter);
    SvREFCNT_inc(result);
    return result;
}



#define MAX_TMPS 1000
unsigned long long int count_tests = 0;
unsigned long long int count_skipped = 0;
unsigned long long int count_kept = 0;

U32 htests() { return count_tests >> 32; }
U32 ltests() { return count_tests & (0xffffffff); }
U32 hskipped() { return count_skipped >> 32; }
U32 lskipped() { return count_skipped & (0xffffffff); }
U32 hkept() { return count_kept >> 32; }
U32 lkept() { return count_kept & (0xffffffff); }

/* Find the next value worth testing.
 * cur: GMP BigInt, the current value to try
 * mult: GMP BigInt, the step to increment cur by
 * sc: arrayref of [mod, ?, vec, ...] arrayrefs
 *
 * Repeatedly increase cur by mult until a value is found such that
 * vec(sc[i].vec, (cur % sc[i].mod), 1) == 0 for all i.
 *
 * Returns the first qualifying value of cur.
 *
 * Note: optimizations here assume we'll get to do many tests before
 * having to recalculate at (roughly) floor(MAXUINT / mult); if mult is
 * high enough (plausible for high prime runs) that may become invalid.
 *
 * TODO:
 *   further optimize main loop: rearrange the bit-vectors so that we need
 * only calculate (modval = iter % scm->mod), and replace the in-loop
 * increments of count_tests and count_skipped with out-loop additions by
 * a calculated value. The rearranged bit-vectors should need recalculating
 * only around every MAX_INT iterations.
 */
SV* cnext(SV* cur, SV* mult, AV* sc, I32 rebuild, SV* max) {
    int found;
    I32 i;
    UV off;
    UV modval;

    if (!scl || rebuild) _init_sc(cur, mult, sc);

    while (1) {
        for (; scl->iter < scl->limit; ++scl->iter) {
            found = 1;
            for (i = 0; i < scl->len; ++i) {
                scmod* scm = &scl->sc[i];
                ++count_tests;
                modval = ((scl->iter % scm->mod) * scm->multmod + scm->basemod) % scm->mod;
                off = modval >> 3;
                if ((off < scl->sc[i].veclen)
                        && (scl->sc[i].vec[off] & (1 << (modval & 7)))) {
                    ++count_skipped;
                    found = 0;
                    break;
                }
            }
            if (found) {
                /* return base + iter * mult */
                SV* result = _calc_result();
                ++scl->iter;    /* ready for next time round */
                ++count_kept;
                return result;
            }
        }
        /* rebase sc */
        {
            SV *newbase, *cmp;
            int i;
            SAVETMPS;
            newbase = _calc_result();
            SvREFCNT_dec(scl->base);
            scl->base = newbase;
            scl->iter = 0;
            for (i = 0; i < scl->len; ++i) {
                scmod* scm = &scl->sc[i];
                scm->basemod = sv_2uv(
                        amagic_call(newbase, scm->svmod, modulo_amg, 0));
            }
            /* if $cur > $max, return $cur so we can abort */
            cmp = amagic_call(newbase, max, gt_amg, 0);
            if (SvTRUE_nomg(cmp)) {
                FREETMPS;
                return newbase;
            }
            FREETMPS;
        }
    }
    /* NOTREACHED */
}
INLINE

sub tests {
    (Constraint::htests() * $BIT32) + Constraint::ltests();
}
sub kept {
    (Constraint::hkept() * $BIT32) + Constraint::lkept();
}
sub skipped {
    (Constraint::hskipped() * $BIT32) + Constraint::lskipped();
}

sub next {
    my($self) = @_;
    my $mult = $self->{'mult'};
    my $cur = $self->{'cur'} + $mult;
    my $rebuild = 0;
    do {
        my $trigger = ($self->{'pend'}[0] || [0])->[0];
        if ($trigger && $trigger < $cur) {
            $cur = $self->catchup($cur);
            $mult = $self->{'mult'};
            $rebuild = 1;
        }
    };
    my $sc = $self->{'sc'};
    my($t, $u) = (0, 0);
    $cur = Constraint::cnext($cur, $mult, $sc, $rebuild, $self->max);
    $self->{'cur'} = $cur;
    $self->{'tests'} += $t;
    $self->{'skipped'} += $u;
    $self->{'kept'} += 1;
    return undef if $cur > $self->{'max'};
    return $cur;
}

#
# If n+kd = xy^z for fixed x and z, we can just search for appropriate y.
#
sub fix_power {
    my($self, $k, $x, $z, $opt_mq) = @_;
    require Constraint::Power;
    return Constraint::Power->new($self, $k, $x, $z, $opt_mq);
}

#
# If n+ad forces a factorization of n+bd, we have some special cases.
#
sub fix_fact {
    my($self) = @_;
    require Constraint::Fact;
    return Constraint::Fact->new($self);
}

sub disallowed {
    my($self, $n, $v) = @_;
    my $c = $self->{'c'}{$n};
    return vec($c->[2], $v, 1) ? 1 : 0 if $c;
    for my $fp (factor_exp($n)) {
        my $pp = $fp->[0] ** $fp->[1];
        my $c = $self->{'c'}{$pp} or return 0;
        return 0 unless vec($c->[2], ($v % $pp), 1);
    }
    return 1;
}

sub require {
    my($self, $p, $v, $min) = @_;
    for (grep $_ != $v, 0 .. $p - 1) {
        $self->suppress($p, $_, $min);
    }
}

sub suppress {
    my($self, $p, $v, $min, $depend) = @_;
    $min ||= 0; $depend ||= 0;
    ($debug > 3) && warn "s: [ $p, $v, $min, $depend ]\n";

    my $c = $self->c($p);
    if ($min > $self->{'min'}) {
        (($debug > 2) && warn "suppress @{[ $depend ? '' : 'in' ]}dependent $v(mod $p): ignore pend\n"),
                return 0 if vec($c->[2], $v, 1);
        ($debug > 1) && warn "suppress (pend) $v(mod $p), d>$min\n";
        push @{ $self->{'pend'} }, [ $min, sub { $self->suppress($p, $v) } ];
        return 0;
    }
    if ($depend) {
        my $effect = 0;
        if (vec($c->[1], $v, 1)) {
            ($debug > 2) && warn "suppress $v(mod $p): now dependent\n";
            vec($c->[1], $v, 1) = 0;
            --$c->[3];
            $effect = 1;
        }
        (($debug > 2) && warn "suppress $v(mod $p): ignore\n"),
                return $effect if vec($c->[2], $v, 1);
        ($debug > 2) && warn "suppress dependent $v(mod $p): apply\n";
        vec($c->[2], $v, 1) = 1;
        $c->[5] = 1;
    } else {
        (($debug > 2) && warn "suppress $v(mod $p): ignore\n"),
                return 0 if vec($c->[2], $v, 1);
        ($debug > 1) && warn "suppress independent $v(mod $p): apply\n";
# FIXME: check for "all suppressed" here to shortcircuit the cascade storm
        vec($c->[2], $v, 1) = 1;
        vec($c->[1], $v, 1) = 1;
        ++$c->[3];
        $c->[5] = 1;
    }

    for my $d (@{ $c->[7] }) {
        # p | d
        for (0 .. $d / $p - 1) {
            $self->suppress($d, $_ * $p + $v, 0, 1);
        }
    }

    DP: for my $d (@{ $c->[6] }) {
        # d | p
        my $dv = $v % $d;
        for (0 .. $p / $d - 1) {
            next DP unless vec($c->[2], $_ * $d + $dv, 1);
        }
        for (0 .. $p / $d - 1) {
            my $mv = $_ * $d + $dv;
            next unless vec($c->[1], $mv, 1);
            ($debug > 2) && warn "suppress $mv(mod $p): self-dependent\n";
            vec($c->[1], $mv, 1) = 0;
            --$c->[3];
        }
        $self->suppress($d, $dv, 0, 0);
    }
    return 1;
}

#
# Disallow value | n+kd whenever n+kd > min.
#
sub series {
    my($self, $value, $min) = @_;
    $min ||= 0;
    ($debug > 1) && warn "series (0, $value, $min)\n";
    my($n, $f) = @$self{qw/ n f /};
    for my $k (1 .. $f - 1) {
        my $g = $k->bgcd($value);
        my($gn, $gnr) = $n->bdiv($g);
        next if $gnr != 0;    # No constraint: the suppressed modulus can't happen
        my($gk, $gmod) = ($k / $g, $value / $g);
        my $kmin = int(($min - $n) / $k) + 1;
        $self->suppress($gmod, (-$gn * $gk->bmodinv($gmod)) % $gmod, $kmin);
    }
    ($debug > 1) && warn "series (0, $value, $min): applied\n";
}

#
# array insert by binary chop, using key $value->[0]
# it is known the element to be inserted is not already present
#
sub _insert {
    my($aref, $val) = @_;
    # we know A[low] < val < A[high]
    my($low, $high) = (-1, @$aref + 0);
    while ($low + 1 < $high) {
        my $mid = ($low + $high) >> 1;
        if ($aref->[$mid][0] < $val->[0]) {
            $low = $mid;
        } else {
            $high = $mid;
        }
    }
    splice(@$aref, $low + 1, 0, $val);
}

#
# Return arrayref of factors f of n: f>1, fp=n -> p prime.
#
sub _highfactors {
    my $n = shift;
    return [ grep $_ > 1, map $n / $_->[0], factor_exp($n) ];
}

#
# Given n, val, mod: return min(m >= 0: m+n == val (mod mod))
#
sub _diff_to_mod {
    my($n, $val, $mod) = @_;
    my $offset = $val - ($n % $mod);
    $offset < 0 ? $offset + $mod : $offset;
}

sub factor_exp {
    my($n) = @_;
    return () if $n == 1;
    return Math::Prime::Util::factor_exp($n);
}

sub numify {
    my($n) = @_;
    my $v = ref($n) ? $n->intify : $n;
    return $v if $v;
    # Math::GMP intifies to 0 if out of range of IV/UV
    return "$n" + 0;
}

1;
