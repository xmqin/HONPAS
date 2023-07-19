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
SIESTA_ARCH=darwin
#
FC=f90
FC_ASIS=$(FC)
#
FFLAGS=  -O
FFLAGS_DEBUG= -g
RANLIB=echo "... do we need to ranlib this? :"
LDFLAGS=
COMP_LIBS=linalg.a
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
LIBS= -lblas -lU77 $(NETCDF_LIBS)
SYS=bsd
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
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
