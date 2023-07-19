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
SIESTA_ARCH=pathf90-2.0-mpi
#
# For an opteron cluster with Pathscale compiler at UniCan -- parallel
#
FC=/localhome/jjunquer/ATC/mpich_gm64/bin/mpif90
FC_ASIS=$(FC)
#
FFLAGS= -m64 -ipa -Ofast -Wl,-R/opt/PathScale/lib/2.0/
FFLAGS_DEBUG= -g 
LDFLAGS= -m64 -ipa -Ofast
RANLIB=/opt/PathScale/x86_64-pathscale-linux/bin/ranlib
COMP_LIBS=
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=            #  -DCDF
#
KINDS=4 8
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=
MPI_LIBS=
FPPFLAGS_MPI=-DMPI
#
LIBS=/localhome/jjunquer/Libraries/SCALAPACK/LIB/libscalapack.a \
	/localhome/jjunquer/Libraries/BLACS/LIB/blacsF77init_MPI-LINUX-0.a \
	/localhome/jjunquer/Libraries/BLACS/LIB/blacs_MPI-LINUX-0.a \
	/localhome/jjunquer/Libraries/BLACS/LIB/blacsCinit_MPI-LINUX-0.a \
	/localhome/jjunquer/acml2.5.1/pathscale64/lib/libacml.a 
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) # -DGRID_DP
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
