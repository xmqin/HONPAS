# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for PLSTM 
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

OBJDIR=Obj

all: plstm plsts

include ../../../$(OBJDIR)/arch.make
vpath %.F90 ../../Grid

#
# The m_gridfunc module (a dependency) has netCDF support, which is currently
# not used by plstm.
# Compiling with netCDF in some systems forces loading the parallel libraries, which
# sometimes prevents interactive use of the program.
# For now we deactivate netCDF in m_gridfunc in this Makefile.
# When netCDF support is added to plstm, the proper symbol below can be defined, if needed.
#
USE_NETCDF_IN_PLSTM=
ifdef USE_NETCDF_IN_PLSTM
   INCFLAGS += $(NETCDF_INCFLAGS)
   LIBS = $(NETCDF_LIBS)
else
   # Remove netCDF functionality from m_gridfunc
   FPPFLAGS += $(DEFS_PREFIX)-UCDF
endif

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

plstm.o: m_gridfunc.o options.o
plsts.o: m_gridfunc.o
plstm: options.o plstm.o m_gridfunc.o
	$(FC) $(LDFLAGS) -o $@  plstm.o m_gridfunc.o options.o $(LIBS)


plsts: plsts.o m_gridfunc.o 
	$(FC) $(LDFLAGS) -o $@  plsts.o m_gridfunc.o $(LIBS)


clean: 
	rm -f *.o *.mod plstm plsts

