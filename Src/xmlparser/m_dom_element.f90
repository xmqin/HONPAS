module m_dom_element

use m_dom_types
use m_dom_namednodemap
use m_dom_nodelist
use m_dom_attribute
use m_dom_document
use m_dom_debug
use m_dom_node
use m_strings

private

  !-------------------------------------------------------   
  ! METHODS FOR ELEMENT NODES
  !-------------------------------------------------------   
  public :: getTagName
  public :: getElementsByTagName
  public :: getElementsByTagAttrName
  public :: getAttribute
  public :: getAttributeNode
  public :: setAttribute
  public :: setAttributeNode
  public :: removeAttribute
  public :: normalize       !--- combines adjacent text nodes ---!

CONTAINS

  !-----------------------------------------------------------
  !  METHODS FOR ELEMENT NODES
  !-----------------------------------------------------------
  function getTagName(element)

    type(fnode), intent(in) :: element   
    type(string)            :: getTagName

    if (element % nodeType == ELEMENT_NODE) then
       getTagName = element % nodeName 
    else
       getTagName = ''
    endif

  end function getTagName

  !-----------------------------------------------------------
 

 function getElementsByTagName(element, tag) result(nodelist)
    type(fnode), pointer         :: element
    character(len=*), intent(in) :: tag
    type(fnodeList), pointer     :: nodelist 

    type(fnode), pointer        :: np

    nodelist => null()

    np => element
    if (dom_debug) print *, "Going into search for tag: ", trim(tag)
    call search(np)

    CONTAINS

    recursive subroutine search(np)
    type(fnode), pointer        :: np

    type(string)                :: name

    !
    ! Could replace the calls to helper methods by direct lookups of node 
    ! components to make it faster.
    ! 
    
    do
       if (.not. associated(np)) exit
       select case(np%nodeType)

          case(DOCUMENT_NODE) 
             ! special case ... search its children 
             if (hasChildNodes(np)) call search(getFirstChild(np))
             ! will exit for lack of siblings
          case(ELEMENT_NODE)

             name = getNodeName(np)
             if (dom_debug) print *, "exploring node: ", char(name)
             if ((tag == "*") .or. (tag == name)) then
                call append(nodelist,np)
                if (dom_debug) print *, "found match ", nodelist%length
             endif
             if (hasChildNodes(np)) call search(getFirstChild(np))

          case default
             
             ! do nothing

        end select

        if (associated(np,element)) exit  ! no siblings of element...
        np => getNextSibling(np)

     enddo
    
    end subroutine search

  end function getElementsByTagName

  !-----------------------------------------------------------

 function getElementsByTagAttrName(element, tag, attr,value) result(nodelist)
    type(fnode), pointer         :: element
    character(len=*), intent(in) :: tag,attr
    character(len=*), intent(in) :: value

    type(fnodeList), pointer     :: nodelist, sametag
    integer                      :: i, length
    type(fnode), pointer         :: np, attr_n
    logical                      :: dom_debug_local
    type(string)                 :: value_s,name_s

    nodelist => null()

    dom_debug_local = .false.

    sameTag => getElementsByTagName(element,tag)
    if (dom_debug .or. dom_debug_local) then
       print *, "Going into search for nodes with tag: ", trim(tag)
       print *, "Going into search for node with attr: ", trim(attr)
    endif

    length=getLength(sameTag)

    do i=0,length-1
       np => item(sameTag,i)
       attr_n => getAttributeNode(np,attr)
       if (associated(attr_n)) then
          name_s = getNodeName(attr_n)
          value_s = getNodeValue(attr_n)
          if (len(value_s) > 0) then
             if (trim(char(value_s)) == trim(value))then
                call append(nodelist,np)
             endif
          endif
       endif
    enddo

  end function getElementsByTagAttrName

!--------------------------------------------------------------

  function getAttribute(element, name)
    
    type(fnode), intent(in) :: element
    character(len=*), intent(in) :: name
    type(string)                 :: getAttribute

    type(fnode), pointer :: nn

    getAttribute = ""  ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    nn => getNamedItem(element%attributes,name)
    if (.not. associated(nn)) RETURN
    
    getAttribute = nn%nodeValue

        
  end function getAttribute

  !-----------------------------------------------------------

  function getAttributeNode(element, name)
    
    type(fnode), intent(in) :: element
    type(fnode), pointer    :: getAttributeNode
    character(len=*), intent(in) :: name

    getAttributeNode => null()     ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    getAttributeNode => getNamedItem(element%attributes,name)

  end function getAttributeNode
  
  !-----------------------------------------------------------

  subroutine setAttributeNode(element, newattr)
    type(fnode), pointer :: element
    type(fnode), pointer :: newattr

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) then
       if (dom_debug) print *, "not an element node in setAttributeNode..."
       RETURN
    endif

    dummy => setNamedItem(element%attributes,newattr)
     
  end subroutine setAttributeNode

!-------------------------------------------------------------------
  subroutine setAttribute(element, name, value)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    type(fnode), pointer      :: newattr

    newattr => createAttribute(name)
    call setValue(newattr,value)
    call setAttributeNode(element,newattr)

  end subroutine setAttribute

  !-----------------------------------------------------------

  subroutine removeAttribute(element, name)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) RETURN
    if (.not. associated(element%attributes)) RETURN

    dummy => removeNamedItem(element%attributes,name)
     
  end subroutine removeAttribute

  !-----------------------------------------------------------
  recursive subroutine normalize(element)
    type(fnode), pointer         :: element

    type(fnode), pointer        :: np, ghost
    logical                     :: first

    type(fnode), pointer        :: head

    first = .true.  ! next Text node will be first

    if (dom_debug) print *, "Normalizing: ", trim(element%nodeName)
    np => element%firstChild
    ! 
    do
       if (.not. associated(np)) exit
       select case(np%nodeType)

          case(TEXT_NODE) 
             if (first) then
                if (dom_debug) print *, "normalize: found first in chain"
                head => np
                first = .false.
                np => getNextSibling(np)
             else                    ! a contiguous text node
                if (dom_debug) print *, "normalize: found second in chain"
                head%nodeValue = head%nodeValue // np%nodeValue
                head%nextSibling => np%nextSibling
                if (associated(np,np%parentNode%lastChild)) then
                   np%parentNode%lastChild => head
                   head%nextSibling => null()
                else
                   np%nextSibling%previousSibling => head
                endif
                ghost => np
                np => getNextSibling(np)
                call destroyNode(ghost)
             endif

          case(ELEMENT_NODE)

             first = .true.
             if (dom_debug) print *, "element sibling: ", trim(np%nodeName)
             if (hasChildNodes(np)) call normalize(np)
             np => getNextSibling(np)

          case default
             
             ! do nothing, just mark that we break the chain of text nodes
             if (dom_debug) print *, "other sibling: ", trim(np%nodeName)
             first = .true.
             np => getNextSibling(np)

        end select

     enddo

    end subroutine normalize


end module m_dom_element
