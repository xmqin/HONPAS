.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libxmlparser.a 
include ../arch.make
#
OBJFILES= m_debug.o m_buffer.o m_elstack.o m_dictionary.o m_io.o \
          m_reader.o m_converters.o m_charset.o  m_entities.o \
          m_xml_error.o  m_fsm.o  m_xml_parser.o flib_sax.o

#
libxmlparser.a: $(OBJFILES) 
	$(AR) $(ARFLAGS_EXTRA) cru libxmlparser.a $(OBJFILES)
	-$(RANLIB) libxmlparser.a
	cp libxmlparser.a ..
	@cp -p *.*d ..
	@echo "xmlparser/work.pc" >> ../work.pcl
#
#         (Note trick for Intel Fortran Compiler)
#
INCFLAGS:= $(INCFLAGS) -I../wxml
#
#
PARENT_VPATH=`grep '^VPATH' ../Makefile`

wxml: ../wxml/flib_wcml.f90
	(cd ../wxml ; make "$(PARENT_VPATH)/wxml")

clean:
	rm -f *.o *.*d *.a 
#
# DO NOT DELETE THIS LINE - used by make depend
m_reader.o: m_io.o m_debug.o
m_elstack.o: m_buffer.o
m_dictionary.o: m_buffer.o
m_entities.o: m_buffer.o m_xml_error.o
m_converters.o: m_debug.o
m_fsm.o: m_buffer.o m_dictionary.o m_charset.o m_entities.o m_elstack.o
m_xml_error.o: m_elstack.o
m_xml_parser.o: m_buffer.o m_elstack.o m_reader.o m_fsm.o m_dictionary.o \
                m_debug.o m_xml_error.o
flib_sax.o:  m_dictionary.o m_xml_parser.o m_converters.o m_xml_error.o

