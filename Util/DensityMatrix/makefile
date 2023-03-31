# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for DM to DM-netCDF converter
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
default: dm2cdf cdf2dm dmbs2dm dmUnblock
all: dm2cdf cdf2dm dmfilter dmbs2dm dmUnblock dm_noncol_sign_flip4
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
OPT_OBJS=m_getopts.o f2kcli.o 
#
dmbs2dm.o: $(OPT_OBJS)
dmbs2dm: dmbs2dm.o $(OPT_OBJS)
	$(FC) $(LDFLAGS) -o $@ dmbs2dm.o $(OPT_OBJS)
#
dmUnblock.o: $(OPT_OBJS)
dmUnblock: dmUnblock.o $(OPT_OBJS)
	$(FC) $(LDFLAGS) -o $@ dmUnblock.o $(OPT_OBJS)
#
dm2cdf: dm2cdf.o
	$(FC) $(LDFLAGS) -o $@ dm2cdf.o  $(NETCDF_LIBS)
#
cdf2dm: cdf2dm.o
	$(FC) $(LDFLAGS) -o $@ cdf2dm.o  $(NETCDF_LIBS)
#
dmfilter: dmfilter.o
	$(FC) $(LDFLAGS) -o $@ dmfilter.o
#
dm_noncol_sign_flip4: dm_noncol_sign_flip4.o
	$(FC) $(LDFLAGS) -o $@ dm_noncol_sign_flip4.o
#
clean: 
	rm -f *.o dm2cdf cdf2dm dmfilter dmbs2dm dmUnblock *.o *.*d
	rm -f dm_noncol_sign_flip4
#
m_getopts.o: f2kcli.o






