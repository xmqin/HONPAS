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
SIESTA_ARCH=gfortran-nolibs-netcdf
#
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
#
# --- Edit the location of your netcdf files
#
NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2-gfortran
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
#
FPPFLAGS=-DGFORTRAN -DCDF -DFC_HAVE_FLUSH -DFC_HAVE_ABORT     # Note this !!
COMP_LIBS=linalg.a
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
LIBS=$(NETCDF_LIBS)
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

