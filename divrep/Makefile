MPUGMP = /src/perl/Math-Prime-Util-GMP
PERL = /opt/maths-GMP/lib/perl5/5.34.0/x86_64-linux/CORE
COUL = coul.c coulfact.c diag.c rootmod.c
HOUL = coulfact.h diag.h rootmod.h
CFACTOR = ${MPUGMP}/factor.c ${MPUGMP}/prime_iterator.c ${MPUGMP}/ecm.c ${MPUGMP}/pbrent63.c ${MPUGMP}/isaac.c ${MPUGMP}/tinyqs.c ${MPUGMP}/squfof126.c ${MPUGMP}/simpqs.c ${MPUGMP}/primality.c ${MPUGMP}/utility.c ${MPUGMP}/gmp_main.c ${MPUGMP}/bls75.c ${MPUGMP}/real.c ${MPUGMP}/ecpp.c
HFACTOR = ${MPUGMP}/factor.h

coul: Makefile ${COUL} ${HOUL} ${CFACTOR} ${HFACTOR}
	gcc -o coul -g -DSTANDALONE ${COUL} ${CFACTOR} -I${MPUGMP} -I${PERL} -lgmp -lm
