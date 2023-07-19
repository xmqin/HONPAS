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
#
# File for regatta power4 processors, parallel version.
# Source: Fedwa El-Mellouhi <f.el.mellouhi@umontreal.ca>
# Note the use of ibm_pessl.f as auxiliary file.
# 
SIESTA_ARCH=ibmp4-mpi

FC=mpxlf
FC_ASIS=$(FC)

FFLAGS=-O3 -qarch=pwr4 -qtune=auto -qcache=auto  -qstrict

FFLAGS_DEBUG= -g -C -qflag=W:W -qinitauto -qsave  -qnolm -qstrict
# -qipa gives speed optimization of a few percent, but takes long to
link
# -qmaxmem=16000 set the maximum memory to 16000, by default it is set
to 2048ko

LIBS= -Wl,-bD:2000000000 -qsave \
       -lpessl -lblacs  \
       /usr/lpp/ppe.poe/lib/libmpi.a \
       -qnolm

#      /home/root/SCA/SCALAPACK/scalapack_RS6K.a \
#      /home/root/SCA/SCALAPACK/pblas_RS6K.a \
#      /home/root/SCA/SCALAPACK/tools_RS6K.a \
#      -lpessl \
#      -lblacs \
#      /usr/lpp/ppe.poe/lib/libmpi.a \
#      -qnolm

# essl lacks zhegv of lapack used for scalar diagonalization.
# It is not called in mpi runs, but it is within a not preprocessed
# if, and requires to be resolved. The following line is thus added
# to link the scalar lapack and blas routines provided within the siesta
# package, even though it has no influence on the run
COMP_LIBS=linalg.a
#
# Add code for descinit, not in pessl
#
SYS=ibm_pessl
RANLIB=ranlib
MPI_INCLUDE=/usr/lpp/ppe.poe/include
MPI_INTERFACE=libmpi_f90.a
#
FPPFLAGS=-WF,-DMPI
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

#
