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
# Makefile for test-xml CML-file diffing tool
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
#
.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
.PHONY: wxml xmlparser clean
#
default: test-xml
OBJDIR=Obj
VPATH=../../Src
#
include ../../$(OBJDIR)/arch.make
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
#
INCFLAGS:=$(INCFLAGS) -I../../$(OBJDIR)/xmlparser # To pick up mod files
INCFLAGS:=$(INCFLAGS) -I../../$(OBJDIR)/wxml  # To pick up mod files
#
LIBS= ../../$(OBJDIR)/xmlparser/libxmlparser.a ../../$(OBJDIR)/wxml/libwxml.a
#
DOM_OBJ_FILES= flib_dom.o \
          m_dom_attribute.o m_dom_element.o m_dom_node.o m_dom_types.o \
          m_dom_debug.o m_dom_error.o m_dom_nodelist.o m_dom_utils.o \
          m_dom_document.o m_dom_namednodemap.o m_dom_parse.o \
          m_strings.o
#
test.o: f2kcli.o m_getopts.o m_compare.o
m_compare.o: m_compare_tol.o corresponding_node.o string_utilities.o
corresponding_node.o: string_utilities.o
string_utilities.o: m_compare_tol.o	
m_getopts.o: f2kcli.o
m_compare_tol.o: $(DOM_OBJ_FILES)
corresponding_node.o: $(DOM_OBJ_FILES)
#
OTHER_OBJS=pxf.o
#
OBJS=$(DOM_OBJ_FILES) test.o f2kcli.o m_getopts.o \
     corresponding_node.o string_utilities.o \
     compare_m.o m_compare_tol.o $(OTHER_OBJS)
#
test-xml: wxml xmlparser $(OBJS)
	@echo
	@echo  "-- Building test-xml program. "
	@echo  "-- Some compilers might have trouble with the DOM subsystem. "
	@echo  "-- Make sure also that 'INCFLAGS' is part of the default compilation rules in arch.make"
	@echo
	$(FC) $(FFLAGS) $(LDFLAGS) -o test-xml $(OBJS) $(LIBS)

LIST_OBJS= list_nodes.o $(DOM_OBJ_FILES) $(OTHER_OBJS)
list:   
	$(FC) $(FFLAGS) $(LDFLAGS)  -o list $(LIST_OBJS) $(LIBS)

clean:
	rm -f *.o *.*d  test-xml list
#
#
wxml:
	@echo "Making sure that the wxml library is compiled..."
	(cd ../../$(OBJDIR)/wxml ; make "VPATH=$(VPATH)/wxml")
xmlparser:
	@echo "Making sure that the xmlparser library is compiled..."
	(cd ../../$(OBJDIR)/xmlparser ; make "VPATH=$(VPATH)/xmlparser")
#
dep:
	sfmakedepend --depend=obj --modext=o *.f *.f90 *.F *.F90 
#
# DO NOT DELETE THIS LINE - used by make depend
compare_m.o: corresponding_node.o flib_dom.o m_compare_tol.o m_strings.o
compare_m.o: string_utilities.o
corresponding_node.o: flib_dom.o string_utilities.o
flib_dom.o: m_dom_attribute.o m_dom_document.o m_dom_element.o
flib_dom.o: m_dom_namednodemap.o m_dom_node.o m_dom_nodelist.o m_dom_parse.o
flib_dom.o: m_dom_types.o m_dom_utils.o m_strings.o
m_compare_tol.o: m_strings.o
m_dom_attribute.o: m_dom_node.o m_dom_types.o m_strings.o
m_dom_document.o: m_dom_types.o m_strings.o
m_dom_element.o: m_dom_attribute.o m_dom_debug.o m_dom_document.o
m_dom_element.o: m_dom_namednodemap.o m_dom_node.o m_dom_nodelist.o
m_dom_element.o: m_dom_types.o m_strings.o
m_dom_namednodemap.o: m_dom_types.o m_strings.o
m_dom_node.o: m_dom_debug.o m_dom_error.o m_dom_namednodemap.o m_dom_nodelist.o
m_dom_node.o: m_dom_types.o m_strings.o
m_dom_nodelist.o: m_dom_types.o
m_dom_parse.o: m_dom_debug.o m_dom_document.o m_dom_element.o m_dom_node.o
m_dom_parse.o: m_dom_types.o
m_dom_types.o: m_strings.o
m_dom_utils.FoX.o: m_dom_debug.o m_dom_document.o m_dom_element.o
m_dom_utils.FoX.o: m_dom_namednodemap.o m_dom_node.o m_dom_types.o m_strings.o
m_dom_utils.o: m_dom_debug.o m_dom_document.o m_dom_element.o
m_dom_utils.o: m_dom_namednodemap.o m_dom_node.o m_dom_types.o m_strings.o
string_utilities.o: m_compare_tol.o m_strings.o
test.o: compare_m.o flib_dom.o m_compare_tol.o
compare_tol_m.o: m_compare_tol.o
