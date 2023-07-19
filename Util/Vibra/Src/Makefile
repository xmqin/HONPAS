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
default: fcbuild vibrator
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
TOPSRC:=$(shell pwd)/../../../Src
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
#
FDF=libfdf.a
FDF_MAKEFILE=$(TOPSRC)/fdf/makefile
FDF_INCFLAGS:=-I $(TOPSRC)/fdf $(INCFLAGS)
$(FDF): 
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC_SERIAL)" \
                          "VPATH=$(TOPSRC)/fdf" \
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
BUILD-SRCS= fcbuild.f chkdim.f recoor.f io.f
BUILD-OBJS=$(BUILD-SRCS:.f=.o)

VIBRA-SRCS= vibrator.f chkdim.f recoor.f hermdp.f paste.f klines.f outbands.f \
            reclat.f parse.f io.f
VIBRA-OBJS=$(VIBRA-SRCS:.f=.o)

#
fcbuild.o vibrator.o:  vibra.h
#
fcbuild: $(FDF) $(BUILD-OBJS)
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD-OBJS) $(FDF)

vibrator : $(FDF) $(VIBRA-OBJS) 
	$(FC) -o vibrator \
	       $(LDFLAGS) $(VIBRA-OBJS) $(FDF)
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f fcbuild vibrator *.o *.mod libfdf.a
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)



