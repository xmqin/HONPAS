program xml2psf

  use pseudopotential, only: pseudopotential_t, pseudo_write_formatted
  use m_pseudo_types,  only: pseudo_t
  use m_pseudo,        only: begin_element, end_element, pcdata_chunk
  use m_pseudo,        only: pseudo    ! Shared structure
  use flib_sax,        only: xml_t, open_xmlfile, xml_parse

 implicit none 

 type(pseudopotential_t)         :: p
 type(pseudo_t), pointer         :: psxml
 type(xml_t)                     :: fxml

 integer :: iostat

 call open_xmlfile("PSXML",fxml,iostat)
 if (iostat /=0) stop "Cannot open file"
                                                                         
 call xml_parse(fxml, begin_element,end_element,pcdata_chunk,verbose=.false.)

 psxml => pseudo

 call xml2psf_helper( psxml, p )
 call pseudo_write_formatted("PSF",p)

end program xml2psf
