# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for LWF to netCDF converter
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
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: lwf2cdf
all: lwf2cdf
#
OBJDIR=Obj
include ../../$(OBJDIR)/arch.make
#
# We assume that the top arch.make defines the relevant NetCDF symbols,
# including the INCFLAGS line. If not, uncomment the following lines, 
# using appropriate information
#
# NETCDF_ROOT=/some/path
# INCFLAGS=-I$(NETCDF_ROOT)/include
# NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf
#
#
INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
#
lwf2cdf: lwf2cdf.o
	$(FC) $(LDFLAGS) -o $@ lwf2cdf.o  $(NETCDF_LIBS)
#
clean: 
	rm -f *.o lwf2cdf  *.o *.*d
#





