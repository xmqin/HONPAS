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
SIESTA_ARCH=pgf95-macosx32
# The only thing you should change is the location of the PGI libraries
# on your computer
#
FC=pgf95
FC_ASIS=$(FC)
#
FFLAGS= -tp core2             # Other options might be useful
FFLAGS_CHECKS= -O0 -g -Mchkptr -Mbounds -traceback
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
LIBS= -L/opt/pgi/osx86/7.2-5/lib -llapack -lblas \
       $(MPI_LIBS)  $(NETCDF_LIBS)
SYS=cpu_time
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
# Compile atom.f and electrostatic.f without optimization.
# Make sure that the dependency is explicit, so that these
# lines work with VPATH
#
atom.o: atom.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
electrostatic.o: electrostatic.f
	$(FC) -c $(FFLAGS_DEBUG) $<
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
