package Constraint::Fact;
use parent qw{ Constraint };
use strict;
use warnings;
use Math::Prime::Util qw{ factor_exp is_prime };
use ModFunc qw{ gcd };

=head1 Constraint::Fact

In some cases we can get additional powerful constraints by showing that
two n + kd values are linked in a way that constrains the factorization
of both.

=head2 n = p^5

For example, given:
  n = p^5, p an odd prime
  p ~| d
  a = (p + 1) / 2
  b = (p - 1) / 2 = a - 1
  k = a^2
then:
  1) we will have n + pd = pq^2, so d = q^2 - p^4 for some prime q
  2) n + kd = (aq - bp^2)(aq + bp^2)
  3) So long as aq - bp^2 > 1, we therefore require one of those two
     factors to be prime and the other a prime square.

For p in (3, 5, 7) that gives us:
  n = 3^5, n + 3d = 3q^2 => a = 2, b = 1, n + 4d = (2q - 9)(2q + 9)
  n = 5^5, n + 5d = 5q^2 => a = 3, b = 2, n + 9d = (3q - 50)(3q + 50)
  n = 7^5, n + 7d = 7q^2 => a = 4, b = 3, n + 16d = (4q - 147)(4q + 147)

So we should search by choosing successive primes r, and for each one:
  (try r^2 = aq + bp^2)
    - reject if r^2 - bp^2 is not divisible by a
    - construct q = (r^2 - bp^2) / a
    - reject if q is not prime
    - reject if aq - bp^2 (== r^2 - 2bp^2) is not prime
    - construct d = q^2 - p^4
    - n + pd, n + kd are now satisfied, try the rest
  (try r^2 = aq - bp^2)
    - reject if r^2 + bp^2 is not divisible by a
    - construct q = (r^2 + bp^2) / a
    - reject if q is not prime
    - reject if aq + bp^2 (== r^2 + 2bp^2) is not prime
    - construct d = q^2 - p^4
    - n + pd, n + kd are now satisfied, try the rest

If r' is the next prime after r, we need to be sure that the higher
value from r will not be accepted in preference to the lower value from r'
if the latter is the smaller. For that we need r'^2 - bp^2 > r^2 + bp^2, or
(r'+r)(r'-r) > 2bp^2. Since r' >= r + 2, it is sufficient to require
r > (bp^2 - 2) / 2 which means requiring d > ((b^2.p^4/4 + 1) / a)^2 - p^4
before allowing this strategy to be used. For p in (3, 5, 7) that gives:
  p = 3: a = 2, b = 1, require d > (85/8)^2 - 81 = 31
  p = 5: a = 3, b = 2, require d > (626/3)^2 - 625 = 42,916
  p = 7: a = 4, b = 3, require d > (21613/16)^2 - 2401 = 1,822,293

We also need aq - bp^2 > 1, so d > ((1 + bp^2)^2)/a^2 - p^4, giving:
  p = 3: a = 2, b = 1, require d > 25 - 81 = -56
  p = 5: a = 3, b = 2, require d > 289 - 625 = -336
  p = 7: a = 4, b = 3, require d > 1369 - 2401 = -1032
.. so no actual constraint there.

Specifically when p=3, there are further connections which lead to
a contradiction:
  a = 2, b = 1, so n + 4d = (2q - 9)(2q + 9);
  r^2 cannot be 2q + 9, since that would require q to be even, so r^2 = 2q - 9;
  so q == 2 mod 3;
  now n + d = (q + 9 - 3r)(q + 9 + 3r); since both factors are == 2 mod 3,
  neither can be a square, so we cannot have tau(n + d) = 6.
This gives us the C<[243, 4]> exception in gtauseq:check_exceptions.

=head2 198

For g(198, 12), find_fixed() finds that n + 10 d = 8 p^2; we then find
5(n + 11d) = 11 . (2p + 3)(2p - 3) = 11st.

Assume m_5 is the principal value of p (mod 5), and m_11 the same (mod 11).
We then have a variety of possibilities: m_5 is constrained to be 1 or 4;
if m_11 \in (4, 7) then one or other factor is also divisible by 11; in that
case we need st = 55qr, with q and r prime, and m_5 and m_11 will tell us
which factor is divisible by 5 and which by 11 (possibly both the same).

For m_11 \notin (4, 7) we have st = 5q^2r, with q and r prime; since
s and t differ by 6, clearly they are both odd, and whichever is 1 (mod 4)
must be the square. Again, m_5 will tell us where the factor of 5 is,
and m_5 and m_11 together will tell us which must include the square.

=head2 2401

For g(2401, 3), secondary fix_power finds that we must satisfy the Pell
equation C< 4x^2 - y^2 = 2401 >; but considering possibilities C<(mod 2)>
shows immediately that no solution is possible.

(Future handling of Pell equations may make it possible to discover this
sort of thing automatically.)

=cut

sub new {
    my($class, $c) = @_;
    my $n = $c->n;

    my @fact = factor_exp($n);
    if (@fact == 1 && $fact[0][1] == 5 && ($fact[0][0] & 1)) {
        return Constraint::Fact::P5->new($c, \@fact);
    }

    my $f = $c->f;
    if ($n == 198 && $f > 11) {
        for my $m (2 .. $c->check) {
            # We'd prefer to do this, for consistency, but we don't have
            # $opt_cp available.
            #next if $opt_cp && cp_skip($opt_cp, $n, \@fact);

            next unless is_prime($m);
            next if gcd($m, $n) > 1;

            # n + 11d = 11/5 (2y - 3)(2y + 3) must be 11^2 pq or 11p^2 q
            # with y^2 = (5d + 99) / 4; so we need the smaller of p and q
            # to satisfy p^2 >= (2y - 3) / 5 => 5p^4 + 6p^2 - 18 >= d.
            # Any smaller p, then, must not divide n + 11d.
            my $m2 = $m * $m;
            my $lim = $m2 * $m2 * 5 + $m2 * 6 - 18;
            $c->suppress($m, -18 % $m, $lim);
        }
        return;
    }

    die "Constraint::Fact->new: no special factors known for f($n, $f)";
}

package Constraint::Fact::P5 {
    use parent qw{ Constraint };
    use strict;
    use warnings;
    no warnings qw/ recursion /;
    use Math::Prime::Util qw{ is_prime next_prime };

    sub new {
        my($class, $c, $fact) = @_;
        my $self = $class->SUPER::new(
            'n' => $c->n(),
            'f' => $c->f(),
            'tell_count' => $c->tell_count(),
            't0' => $c->t0(),
            'min' => $c->min(),
            'max' => $c->max(),
            'check' => $c->check(),
            'tau' => $c->tau(),
            'parent' => $c,
        );

        @$fact == 1 && $fact->[0][1] == 5 && ($fact->[0][0] & 1)
                or die "Constraint::Fact::P5->new: expect n = p^5, odd p";
        my $p = Math::GMP->new($fact->[0][0]);
        my $minf = $p + ($p / 2) ** 2;
        $minf < $self->{'f'}
                or die "Constraint::Fact->new: for n=p^5 expect f > $minf";
        my $a = ($p + 1) / 2;
        my $b = $a - 1;
        my $k = $a * $a;
        $self->{'p5'} = {
            p => $p, p2 => $p * $p, p4 => $p * $p * $p * $p,
            a => $a, b => $b, k => $k,
            bp2 => $b * $p * $p,
        };

        $self->{'min'} = $self->_dtoceilr($c->min);
        $self->{'max'} = $self->_dtoceilr($c->max);

        # todo: copy over constraints
        # todo: copy over pend list

        return $self;
    }

    #
    # Calculate floor(r) given d
    #
    sub _dtor {
        my($self, $d) = @_;
        my $p5 = $self->{'p5'};
        # d = q^2 - p^4 => q = sqrt(d + p^4)
        my $q = ($d + $p5->{'p4'})->bsqrt;
        # r^2 = aq +/- bp^2
        my $sq = $p5->{'a'} * $q - $p5->{'bp2'};
        return $sq < 0 ? Math::GMP->new(0) : $sq->bsqrt;
    }

    sub _dtoceilr {
        my($self, $d) = @_;
        my $p5 = $self->{'p5'};
        my $q = 1 + ($d + $p5->{'p4'} - 1)->bsqrt;
        my $sq = $p5->{'a'} * $q + $p5->{'bp2'};
        return $sq < 1 ? Math::GMP->new(0) : 1 + ($sq - 1)->bsqrt;
    }

    #
    # Calculate d given r and direction (-1 or +1)
    #
    sub _rtod {
        my($self, $r, $dir) = @_;
        my $p5 = $self->{'p5'};
        # r^2 = aq +/- bp^2 => q = (r^2 -/+ bp^2) / a
        my $q = ($r * $r - $dir * ($p5->{'bp2'})) / $p5->{'a'};
        # d = q^2 - p^4
        return $q * $q - $p5->{'p4'};
    }

    sub cur {
        my($self) = @_;
        return $self->_rtod($self->{cur}, $self->{dir} // -1);
    }

    sub next {
        my($self) = @_;
        my($a, $bp2, $p4) = @{ $self->{'p5'} }{qw{ a bp2 p4 }};
        my $r = $self->{'cur'};
        my $dir = $self->{'dir'} // +1;
        my($q, $aq, $rem, $alt);
        my $t = 0;

        while (1) {
            ++$t;
            if ($dir > 0) {
                $r = next_prime($r);
                $dir = -1;
                return undef if $r > $self->{'max'};
            } else {
                $dir = +1;
            }
            if ($dir < 0) {
                # (try r^2 = aq + bp^2)
                $aq = $r * $r - $bp2;
                ($q, $rem) = $aq->bdiv($a);
                next if $rem != 0;
                next unless is_prime($q);
                $alt = $aq - $bp2;
                next unless is_prime($alt);
            } else {
                # (try r^2 = aq - bp^2)
                $aq = $r * $r + $bp2;
                ($q, $rem) = $aq->bdiv($a);
                next if $rem != 0;
                next unless is_prime($q);
                $alt = $aq + $bp2;
                next unless is_prime($alt);
            }
            last;
        }
        $self->{'cur'} = $r;
        $self->{'dir'} = $dir;
        $self->{'tests'} += $t;
        $self->{'skipped'} += $t - 1;
        $self->{'kept'} += 1;
        return $q * $q - $p4;
    }
};

1;
