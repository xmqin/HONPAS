module m_dom_document

use m_dom_types
use m_strings

private

  !-------------------------------------------------------  
  ! METHODS FOR DOCUMENT NODES
  !-------------------------------------------------------   
  public :: createDocumentNode
  public :: createDocumentFragment
  public :: createTextNode
  public :: createAttribute
  public :: createElement
  public :: createComment
  public :: createCdataSection

CONTAINS

  !-----------------------------------------------------------
  ! METHODS FOR DOCUMENT NODES
  !-----------------------------------------------------------

  function createDocumentNode()

    type(fnode), pointer :: createDocumentNode
    
    createDocumentNode => createNode()
    createDocumentNode % nodeType = DOCUMENT_NODE
    createDocumentNode % nodeName = "#document"
    
  end function createDocumentNode
!-------------------------------------------------------------------
  function createDocumentFragment()

    type(fnode), pointer :: createDocumentFragment
    
    createDocumentFragment => createNode()
    createDocumentFragment % nodeType = DOCUMENT_FRAGMENT_NODE
    createDocumentFragment % nodeName = "#document-fragment"
    
  end function createDocumentFragment
!-------------------------------------------------------------------
  function createTextNode(data)

    character(len=*), intent(in) :: data
    type(fnode), pointer :: createTextNode
    
    createTextNode => createNode()
    createTextNode % nodeType = TEXT_NODE
    createTextNode % nodeName = "#text"
    createTextNode % nodeValue = data  ! NB need to split this string 
                                       ! across several nodes 
    
  end function createTextNode

  !-----------------------------------------------------------

  function createAttribute(name)
    
    character(len=*), intent(in) :: name
    type(fnode), pointer :: createAttribute
  
    createAttribute => createNode()
    createAttribute % nodeName = name
    createAttribute % nodeType = ATTRIBUTE_NODE
  
  end function createAttribute

  !-----------------------------------------------------------

  function createElement(tagName)
    
    character(len=*), intent(in) :: tagName
    type(fnode), pointer :: createElement
  
    createElement => createNode()
    createElement % nodeName = tagName
    createElement % nodeType = ELEMENT_NODE
  
  end function createElement

  !-----------------------------------------------------------

  function createComment(data)
    
    character(len=*), intent(in) :: data
    type(fnode), pointer :: createComment
  
    createComment => createNode()
    createComment % nodeName = "#comment"
    createComment % nodeValue = data
    createComment % nodeType = COMMENT_NODE
  
  end function createComment

  !-----------------------------------------------------------

  function createCdataSection(data)
    
    character(len=*), intent(in) :: data
    type(fnode), pointer :: createCdataSection
  
    createCdataSection => createNode()
    createCdataSection % nodeName = "#cdata-section"
    createCdataSection % nodeValue = data
    createCdataSection % nodeType = CDATA_SECTION_NODE
  
  end function createCdataSection

  
end module m_dom_document
