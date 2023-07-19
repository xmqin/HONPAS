# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
SIESTA_ARCH=itanium2
#
# Serial compilation
#
FC=efc -Vaxlib
#
FC_ASIS=$(FC)
#
FFLAGS= -O2 -tpp2 -W0
FFLAGS_DEBUG= -g -O0
LDFLAG=#-Vaxlib 
RANLIB=echo
LAPACK=-L/opt/intel/mkl72/lib/64 -lmkl_lapack64 -lmkl -lguide -lpthread
LIBS=$(LAPACK)
SYS=bsd
FPPFLAGS= -DWXML_INIT_FIX -DALLOC_SAVE_BUG
COMP_LIBS=
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#



