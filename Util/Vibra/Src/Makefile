# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Vibra package
#
# Pablo Ordejon, August 1998
# Modified by Alberto Garcia, February 2006
#
# Points to note, until we switch to a better building system:
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
.SUFFIXES: .f .o .a 
#
default: fcbuild vibra
#
# You should have compiled Siesta in the main Src directory (or
# in an architecture-specific object directory) before attempting 
# to compile the Vibra package.
#
OBJDIR=Obj
ARCH_MAKE=../../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
# this path is as seen from the fdf subdir...
TOPSRC:=../../../../Src
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF $(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF $(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
#
FDF=libfdf.a
FDF_MAKEFILE=$(TOPSRC)/fdf/makefile
FDF_INCFLAGS:=-I $(TOPSRC)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC_SERIAL)" \
                          "VPATH=$(TOPSRC)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)" \
                          "MPI_INTERFACE= " \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
BUILD-SRCS= fcbuild.f chkdim.f recoor.f io.f
BUILD-OBJS=$(BUILD-SRCS:.f=.o)

VIBRA-SRCS= vibra.f chkdim.f recoor.f hermdp.f klines.f outbands.f \
            reclat.f parse.f io.f
VIBRA-OBJS=$(VIBRA-SRCS:.f=.o)

#
fcbuild.o vibra.o:  vibra.h
#
fcbuild: $(FDF) $(BUILD-OBJS)
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD-OBJS) $(FDF)

vibra : $(FDF) $(VIBRA-OBJS) 
	$(FC) -o vibra \
	       $(LDFLAGS) $(VIBRA-OBJS) $(FDF)
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f fcbuild vibra *.o *.mod libfdf.a
	rm -rf ./fdf




