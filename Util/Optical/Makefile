# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#
# Makefile for Optical utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a .f90 .F90

# No VPATH should be used here.
# This utility is self-contained and will interfere
# with the optical.F code (in Src)
OBJDIR=Obj

all: optical_input optical

include ../../$(OBJDIR)/arch.make

# This is needed on some systems to avoid loading the parallel libraries, which
# sometimes force running on queuing systems

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive

##FFLAGS=$(FFLAGS_DEBUG)

optical_input: optical_input.o
	$(FC) $(LDFLAGS) -o $@ optical_input.o

optical: optical.o
	$(FC) $(LDFLAGS) -o $@ optical.o

clean: 
	rm -f *.o optical_input optical
