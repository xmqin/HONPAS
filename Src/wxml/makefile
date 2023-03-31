.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
module:  libwxml.a
include ../arch.make
#
OBJFILES= m_wxml_buffer.o m_wxml_array_str.o m_wxml_dictionary.o m_wxml_elstack.o \
          m_wxml_text.o m_wxml_escape.o  m_wxml_core.o \
          m_wxml_overloads.o flib_wxml.o flib_wstml.o flib_wcml.o \
	  m_wxml_error.o m_wcml_coma.o
#
libwxml.a: $(OBJFILES)
	$(AR) $(ARFLAGS_EXTRA) cru libwxml.a $(OBJFILES)
	-$(RANLIB) libwxml.a
	cp libwxml.a ..
	@cp -p *.*d ..
	@echo "wxml/work.pc" >> ../work.pcl
#
#         (Note trick for Intel Fortran Compiler)
#
clean:
	rm -f *.o *.*d *.a
# DO NOT DELETE THIS LINE - used by make depend
flib_wcml.o: flib_wstml.o flib_wxml.o m_wcml_coma.o
flib_wstml.o: flib_wxml.o
flib_wxml.o: m_wxml_core.o m_wxml_overloads.o m_wxml_text.o
m_wcml_coma.o: flib_wstml.o flib_wxml.o
m_wxml_buffer.o: m_wxml_error.o m_wxml_escape.o
m_wxml_core.o: m_wxml_buffer.o m_wxml_dictionary.o m_wxml_elstack.o
m_wxml_core.o: m_wxml_escape.o m_wxml_array_str.o
m_wxml_dictionary.o: m_wxml_array_str.o m_wxml_escape.o
m_wxml_elstack.o: m_wxml_error.o
m_wxml_overloads.o: m_wxml_core.o m_wxml_text.o
