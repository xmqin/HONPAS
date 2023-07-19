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
# Makefile for WFS utilities
#
# The VPATH directive below helps to re-use f2kcli.F90 from the top Src.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
all: readwf readwfx info_wfsx wfs2wfsx wfsx2wfs wfsnc2wfsx
#
OBJDIR=Obj
VPATH=../../Src
#
include ../../$(OBJDIR)/arch.make
#
# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
INCFLAGS:=$(NETCDF_INCFLAGS) $(INCFLAGS) # -I../../$(OBJDIR)
#
# We assume that the top arch.make defines the relevant NetCDF symbols,
# including the NETCDF_INCFLAGS line. If not, uncomment the following lines, 
# using appropriate information
#
# NETCDF_ROOT=/some/path
# INCFLAGS=-I$(NETCDF_ROOT)/include
# NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
#
readwf: readwf.o
	$(FC) $(LDFLAGS) -o $@ $<  
readwfx: readwfx.o
	$(FC) $(LDFLAGS) -o $@ $<  
info_wfsx: info_wfsx.o
	$(FC) $(LDFLAGS) -o $@ $<  
wfs2wfsx: wfs2wfsx.o
	$(FC) $(LDFLAGS) -o $@ $<  
wfsx2wfs: wfsx2wfs.o
	$(FC) $(LDFLAGS) -o $@ $<  
wfsnc2wfsx: wfsnc2wfsx.o
	$(FC) $(LDFLAGS) -o $@ $<  $(NETCDF_LIBS)
#
clean: 
	rm -f *.o  *.*d readwf readwfx wfs2wfsx wfsx2wfs wfsnc2wfsx
#










