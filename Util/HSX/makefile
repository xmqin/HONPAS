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
# Makefile for HSX utilities
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
OBJDIR=Obj
#
all: hsx2hs hs2hsx
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
FFLAGS=$(FFLAGS_DEBUG)
#
#------------------
hsx2hs: hsx_m.o hsx2hs.o
	$(FC) $(LDFLAGS) -o $@ hsx_m.o  hsx2hs.o
#
hs2hsx: hsx_m.o hs2hsx.o
	$(FC) $(LDFLAGS) -o $@ hsx_m.o  hs2hsx.o
#
#------------------
clean: 
	rm -f *.o hsx2hs hs2hsx *.o *.*d
#










