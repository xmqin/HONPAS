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

.SUFFIXES:
.SUFFIXES: .f .f90 .F .F90 .o

VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj

all: eigfat2plot gnubands

include ../../$(OBJDIR)/arch.make


FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI 
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI 
INCFLAGS:=$(INCFLAGS) $(NETCDF_INCFLAGS)

FAT_OBJS=m_getopts.o f2kcli.o 

gnubands: $(OBJS) gnubands.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) gnubands.o

eigfat2plot: $(OBJS) eigfat2plot.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) eigfat2plot.o

clean:
	rm -f *.o *.mod eigfat2plot gnubands

dep:
	-sfmakedepend --depend=obj --modext=o \
		*.f90 *.F90 *.f *.F

# DO NOT DELETE THIS LINE - used by make depend
eigfat2plot.o: f2kcli.o m_getopts.o
gnubands.o: f2kcli.o m_getopts.o
m_getopts.o: f2kcli.o
