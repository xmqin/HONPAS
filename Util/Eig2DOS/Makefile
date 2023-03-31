# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
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
.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
#
all: Eig2DOS
#
include ../../$(OBJDIR)/arch.make
#
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
INCFLAGS:=$(INCFLAGS) $(NETCDF_INCFLAGS)
#
OBJS=m_getopts.o f2kcli.o 
#
default: Eig2DOS
#
Eig2DOS: $(OBJS) Eig2DOS.o
	$(FC) -o $@ $(LDFLAGS) $(OBJS) Eig2DOS.o
#
clean:
	rm -f *.o *.mod Eig2DOS

# Dependecies
m_getopts.o: f2kcli.o
Eig2DOS.o: m_getopts.o f2kcli.o
