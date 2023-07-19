# Makefile for macroave
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: macroave
#
OBJDIR=Obj
VPATH:=$(shell pwd)/../../../Src
include ../../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
# Note that machine-specific files are now in top Src directory.
#
OBJS = io.o iorho.o paste.o macroave.o thetaft.o \
	surpla.o volcel.o recipes.o hdr_io.o
#
#
# This is crude but will have to do for now.
# Note : precision must be the first module
#
MOD_OBJS=defs_basis.o defs_common.o
COM_OBJS=$(OBJS) 
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)
#
$(COM_OBJS): 
#
macroave:  $(ALL_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o macroave $(ALL_OBJS) 
#
permute.o: f2kcli.o
permute:   f2kcli.o permute.o 
	$(FC) $(FFLAGS) $(LDFLAGS) -o permute f2kcli.o permute.o
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f macroave permute *.o  *.a *.pcl *.pc *.mod
	rm -f *.mod
#



