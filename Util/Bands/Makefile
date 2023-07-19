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
all: eigfat2plot new.gnubands gnubands
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
FAT_OBJS=m_getopts.o f2kcli.o 
#
new.gnubands: $(OBJS) new.gnubands.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) new.gnubands.o
#
eigfat2plot: $(OBJS) eigfat2plot.o
	$(FC) -o $@ $(LDFLAGS) $(FAT_OBJS) eigfat2plot.o
#
gnubands: gnubands.o
	$(FC) -o $@ $(LDFLAGS)  gnubands.o
#
clean:
	rm -f *.o *.mod new.gnubands eigfat2plot gnubands

# DO NOT DELETE THIS LINE - used by make depend
m_getopts.o: f2kcli.o
new.gnubands.o: m_getopts.o f2kcli.o
eigfat2plot.o: m_getopts.o f2kcli.o
