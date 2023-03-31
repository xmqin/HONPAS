# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
# Makefile for Client-Server example
#
# JM Soler, J Nara, A Garcia, E Anglada
#
WORK_DIR=$(shell pwd)
SIESTA_DIR=$(WORK_DIR)/../../../..
SRC_DIR=$(SIESTA_DIR)/Src
OBJDIR=Obj
OBJ_DIR=$(SIESTA_DIR)/$(OBJDIR)
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH=$(SRC_DIR)

all: default
default: protoNEB

ARCH_MAKE=$(OBJ_DIR)/arch.make
include $(ARCH_MAKE)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
# Routines used by pipes version (VPATH knows where they are)
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
WXML=libwxml.a
XC=libSiestaXC.a
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(LIBS)
#
libs_collected:
	(cd $(OBJ_DIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WORK_DIR) )
	touch libs_collected
#
protoNEB: libs_collected protoNEB.o 
	$(FC) $(LDFLAGS) -o protoNEB \
	      protoNEB.o $(ALL_LIBS)
#
intermediate: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make clean to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: clean
clean:
	@echo "==> Cleaning all intermediate and executable files"
	rm -f protoNEB
	rm -f *.a *.o *.mod libs_collected

