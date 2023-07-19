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
SIESTA_ARCH=intel-checks
#
# Serial compilation without the need of any installed libraries.
# (That would be too simple: ifc needs the -Vaxlib option at link time)
#
# Following line for Pentium 4 , with static alloc
#FC=ifc  -tpp7 -O2  -Vaxlib -static -xW
#
# More conservative
FC=ifort
#
FC_ASIS=$(FC)
#
FFLAGS= -O1 -ftrapuv -CB -u
FFLAGS_DEBUG= -g 
LDFLAGS=
RANLIB=echo
LIBS=  
SYS=bsd
FPPFLAGS=
COMP_LIBS=linalg.a
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

