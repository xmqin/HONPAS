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
SIESTA_ARCH=gfortran-nolibs
#
# Serial compilation without the need of any installed libraries.
#
FC=gfortran
#
FC_ASIS=$(FC)
#
FFLAGS=-O2
FFLAGS_DEBUG= -g -Wall -Wextra
LDFLAGS=
RANLIB=echo
LIBS=  
SYS=nag
FPPFLAGS=-DGFORTRAN -DFC_HAVE_FLUSH -DFC_HAVE_ABORT          # Note this !!
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

