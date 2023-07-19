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
SIESTA_ARCH=hpcx

FC=mpxlf90_r -bmaxdata:1000000000
FC_ASIS=$(FC)

FFLAGS= -qfixed -O3 -qarch=pwr4 -qtune=pwr4 -qstrict
FFLAGS_DEBUG= -g -C -qflag=W:W -qinitauto -qsave -qmaxmem=16000
# -qipa gives speed optimization of a few percent, but takes long to link

LIBS= lapack.o -L/usr/local/lib -lscalapack -lblacsF77init -lblacs -lessl
#LIBS= lapack.o -L/hpcx/home/e05/e05/jdg/libs -lscalapack -lblacsF77init -lblacs -lessl
COMP_LIBS=
#
SYS=ibm
RANLIB=ranlib
MPI_INCLUDE=/usr/lpp/ppe.poe/include
MPI_INTERFACE=libmpi_f90.a
#
FPPFLAGS=-WF,-DMPI,-DNODAT
FREE_F90=-qsuffix=f=f90 -qfree=f90
FREE_F90_CPP=-qsuffix=cpp=F90 -qfree=f90
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FREE_F90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FREE_F90_CPP) $<
