#
# Modern Makefile for ATM
# Supports the standard mechanisms in the SIESTA tree
# For traditional compilation, needing only an f77 compiler,
# and no CML support, see makefile-old.
#
# ---------------------------------------------------------------
#
.SUFFIXES:
.SUFFIXES: .F90 .f .o .a 
#
default: atm
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Src, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#
OBJDIR=Obj
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
#----------------------------------------------$$$$$
# EDIT the following lines to suit your system
#
AUX_OBJS= auxf95.o
LIBS= $(EXTRA_LIBS)
INCFLAGS=
#--------------------------------------------------------------
# FoX whatnot
#
# FoX build targets:
FoX_configured=../../$(OBJDIR)/FoX/.config
FoX_built=../../$(OBJDIR)/FoX/.FoX
#
# This is how we pick up modules and libraries for FoX:
FoX_FCFLAGS=`../../$(OBJDIR)/FoX/FoX-config --fcflags`
FoX_LIBS=`../../$(OBJDIR)/FoX/FoX-config --libs --wxml`
#
# And add them to global compiler flags:
INCFLAGS:=$(FoX_FCFLAGS) $(INCFLAGS)
#
# First, it needs to be configured. This may have been done
# by the main SIESTA configure, but in case not:
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
#
# Note that we hardwire the need to use the configure script in the main Src directory.
#
TOP_LEVEL=$(shell pwd)/../../Src
$(FoX_configured):
	(cd ../../$(OBJDIR)/FoX; touch arch.make ; \
         CONFIGURE="$(TOP_LEVEL)/FoX/configure"; \
         $$CONFIGURE VPATH="$(TOP_LEVEL)/FoX" \
         FC="$(FC_SERIAL)" FCFLAGS="$(FFLAGS)" \
         --enable-wcml $(DUMMY_FOX) || false )
#
# Note we have to make sure to use the same compiler as SIESTA,
# and pick up all the same FFLAGS, and also remember to maybe ask for a dummy library.
# Note also that the "false" clause will stop the 'make' process if the configuration fails.
#
# Then we want to make FoX itself. Like so:
$(FoX_built): $(FoX_configured)
	(cd ../../$(OBJDIR)/FoX; $(MAKE))
#
# Now ensure that building FoX is inserted into the dependency chain.
# Routines using FoX:
#
pseudoXML.o: $(FoX_built)
#
#--------------------------------------------------------------

#----------------------------------------------$$$$$
#
INCFLAGS:=$(INCFLAGS)  ####         -I../../$(OBJDIR)
#
#
XC_OBJS= xc.o
#
# Replace definition above with line below to use the 
# old Soler-Balbas XC routines.
#
#XC_OBJS= atomxc.o  ldaxc.o ggaxc.o pbexc.o pw92xc.o pzxc.o \
#        pw92c.o exchng.o
#
MAIN_OBJS= atm.o  nucl_z.o   difnrl.o   difrel.o   dmixp.o   dsolv1.o   \
        pseudo.o  wf.o  hsc.o  excorr.o \
        tm2.o  ker.o  wrapup.o   splift.o   spliq.o  \
        velect.o   vionic.o   wtrans.o   v0pp.o  chg_mism.o  \
        dsolv2.o   etotal.o   ext.o  input.o  header.o orban.o   \
        potran.o   potrv.o potrvs.o  potrw.o   prdiff.o  string.o  root.o  \
        genrhs.o  logder.o ode.o denplot.o scilib.o  recipes.o \
	compat_params.o  pcc_exp.o  change_valence.o pswf_store.o \
	coreq.o get_unit.o gauleg.o divdif.o symbols.o \
        prversion.o pseudoXML.o pxf.o force_underflow.o
#
OBJS= $(MAIN_OBJS) $(XC_OBJS) 
#
provide_xml:
	@echo "Make sure that the wxml library is already compiled."
	@echo "Otherwise the following command might fail"
	@echo
#
pseudoXML.o: provide_xml
#
atm:	$(FoX_built) $(OBJS) $(AUX_OBJS)
	$(FC) $(LDFLAGS) -o atm $(OBJS) $(AUX_OBJS) $(LIBS) $(FoX_LIBS)

#
clean:	
	rm -f *.o atm
#
atm.o: radial.h orbital.h param.h charge.h elecpot.h energy.h
#
chg_mism.o: radial.h nonlinear.h linear.h coeffs.h
compat_params.o: compat.h
denplot.o:  radial.h charge.h
difnrl.o: radial.h param.h ion.h elecpot.h
difrel.o: radial.h param.h ion.h elecpot.h
dsolv1.o: radial.h orbital.h ion.h charge.h elecpot.h energy.h
dsolv2.o: radial.h orbital.h param.h ion.h charge.h elecpot.h energy.h
etotal.o: orbital.h param.h energy.h
excorr.o: radial.h param.h
genrhs.o: nonlinear.h
input.o:  radial.h orbital.h param.h compat.h
logder.o: radial.h orbital.h param.h ion.h elecpot.h energy.h ode_blk.h \
          ode_path.h
ode.o	: radial.h ode_blk.h
orban.o : radial.h param.h ion.h elecpot.h
prversion.o: version.h
pseudo.o: radial.h orbital.h param.h elecpot.h
symbols.o: set2.h
tm2.o   : radial.h orbital.h ion.h param.h elecpot.h energy.h tm2_blk.h 
v0pp.o  : radial.h tm2_blk.h
velect.o: radial.h orbital.h param.h charge.h elecpot.h energy.h compat.h
vionic.o: radial.h orbital.h param.h ion.h charge.h 
change_valence.o: radial.h param.h  charge.h ion.h
wf.o hsc.o ker.o : radial.h orbital.h ion.h param.h elecpot.h energy.h
wrapup.o: radial.h orbital.h param.h ion.h charge.h elecpot.h energy.h \
          compat.h
wtrans.o potran.o: plot.h
ninput.o header.o: input.h
pseudoXML.o      : param.h radial.h ion.h orbital.h charge.h
pswf_store.o pseudoXML.o : pseudowave.h
#









