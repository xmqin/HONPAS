#
.SUFFIXES:
.SUFFIXES: .f .o .f90
#
# This can also be achieved with a "FORCE" empty target.      
# xml: FORCE
#	(cd ../../Src/xmlparser ; make)
# FORCE:
#
#
.PHONY: xml clean
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
default: xml2psf
#
OBJDIR=Obj
TOPSRC=../../Src
#
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
INCFLAGS:=$(INCFLAGS) -I../../$(OBJDIR)
LIBS= ../../$(OBJDIR)/libxmlparser.a
#
OBJS= sys.o io.o pseudopotential.o m_pseudo_types.o m_pseudo.o xml2psf_helper.o  xml2psf.o
#
xml2psf: xml $(OBJS)
	$(FC) -o xml2psf $(LDFLAGS) $(OBJS) $(LIBS)
#
xml:
	@echo "Make sure that the xmlparser library is already compiled"
	@echo "The following command might fail if you are using VPATH"
	(cd ../../$(OBJDIR)/xmlparser ; $(MAKE) "VPATH=$(TOPSRC)/xmlparser")
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
clean:
	rm -f *.o *.*d xml2psf
#
# DO NOT DELETE THIS LINE - used by make depend
io.o: sys.o
m_pseudo.o: m_pseudo_types.o
xml2psf.o: m_pseudo.o m_pseudo_types.o pseudopotential.o
xml2psf_helper.o: m_pseudo_types.o pseudopotential.o
