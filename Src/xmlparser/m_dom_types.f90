Module m_dom_types

  use m_strings

  implicit none

  private

  !-------------------------------------------------------   
  ! A GENERIC NODE
  !-------------------------------------------------------   
  type, public :: fnode
     type(string)         :: nodeName    
     type(string)         :: nodeValue   
!!!     character(len=200)    :: nodeName  = ""
!!!     character(len=200)    :: nodeValue = ""
     integer              :: nc              = 0 
     integer              :: nodeType        = 0
     type(fnode), pointer :: parentNode      => null()
     type(fnodeList), pointer :: childNodes  => null()  ! New
     type(fnode), pointer :: firstChild      => null()
     type(fnode), pointer :: lastChild       => null()
     type(fnode), pointer :: previousSibling => null()
     type(fnode), pointer :: nextSibling     => null()
     type(fnode), pointer :: ownerDocument   => null()
     type(fnamedNodeMap), pointer :: attributes => null()
  end type fnode

  !-----------------------------------------------------------
  !  ONE WAY TO IMPLEMENT A NAMEDNODEMAP  (dictionary)
  !-----------------------------------------------------------

  ! Linked list of name/node pairs, with overall length variable

  type, public :: fnamedNode
     type(string)                   :: name
!!!     character(len=100)            :: name
     type(fnode), pointer          :: node => null()
     type(fnamedNode), pointer     :: next => null()
  end type fnamedNode

  type, public :: fnamedNodeMap
     integer :: length = 0
     type(fnamedNode), pointer  :: head => null()
     type(fnamedNode), pointer  :: tail => null()
  end type fnamedNodeMap

  !-----------------------------------------------------------
  !  ONE WAY TO IMPLEMENT A NODELIST 
  !-----------------------------------------------------------

  type, public :: flistNode
     type(fnode), pointer          :: node => null()
     type(flistNode), pointer      :: next => null()
  end type flistNode

  type, public :: fnodeList
     integer                      :: length = 0
     type(flistNode), pointer     :: head => null()
     type(flistNode), pointer     :: tail => null()
  end type fnodeList

!========================================================================
  integer, save, private          :: allocated_nodes = 0
!========================================================================

  !-------------------------------------------------------   
  ! NODETYPES
  !-------------------------------------------------------   
  integer, parameter, public :: ELEMENT_NODE                = 1
  integer, parameter, public :: ATTRIBUTE_NODE              = 2
  integer, parameter, public :: TEXT_NODE                   = 3
  integer, parameter, public :: CDATA_SECTION_NODE          = 4
  integer, parameter, public :: ENTITY_REFERENCE_NODE       = 5
  integer, parameter, public :: ENTITY_NODE                 = 6
  integer, parameter, public :: PROCESSING_INSTRUCTION_NODE = 7
  integer, parameter, public :: COMMENT_NODE                = 8
  integer, parameter, public :: DOCUMENT_NODE               = 9
  integer, parameter, public :: DOCUMENT_TYPE_NODE          = 10
  integer, parameter, public :: DOCUMENT_FRAGMENT_NODE      = 11
  integer, parameter, public :: NOTATION_NODE               = 12

  public :: node_class
  public :: createNode
  public :: destroyNode
  public :: destroyNamedNodeMap
  public :: destroyNodeList
  public :: getNumberofAllocatedNodes

CONTAINS

  function getNumberofAllocatedNodes() result(n)
    integer   :: n

    n = allocated_nodes
  end function getNumberofAllocatedNodes

!--------------------------------------------------------------
  function createNode() result(node)
    type(fnode), pointer  :: node

    allocate(node)
    allocated_nodes = allocated_nodes + 1

  end function createNode
!--------------------------------------------------------------

  function node_class(nodetype) result(class)
    integer, intent(in) :: nodetype
    character(len=10)  ::     class

    select case(nodetype)
    case(ELEMENT_NODE)
       class = "element"
    case(ATTRIBUTE_NODE)
       class = "attribute"
    case(TEXT_NODE)
       class = "text"
    case(COMMENT_NODE)
       class = "comment"
    case(DOCUMENT_NODE)
       class = "document"
    end select
  end function node_class

  subroutine destroyNamedNodeMap(nodemap)
    type(fnamedNodeMap), pointer :: nodemap

    type(fnamednode), pointer  :: nnp
    type(fnode), pointer       :: ghost
    
    if (.not. associated(nodemap)) return
    nnp => nodemap%head
    do while (associated(nnp))
       call unstring(nnp%name)
       ghost => nnp%node
       nnp => nnp%next
       call destroyNode(ghost)      ! We might not want to really destroy
    enddo
  end subroutine destroyNamedNodeMap

  subroutine destroyNodeList(nodelist)
    type(fnodeList), pointer :: nodelist

    type(flistnode), pointer   :: p
    type(fnode), pointer       :: ghost
    
    if (.not. associated(nodelist)) return
    p => nodelist%head
    do while (associated(p))
       ghost => p%node
       p => p%next
       call destroyNode(ghost)      ! We might not want to really destroy
    enddo
  end subroutine destroyNodeList

  recursive subroutine destroyNode(node)
    type(fnode), pointer  :: node
    
    type(fnode), pointer  :: np, ghost
    
    np => node
    do while (associated(np))
       if (associated(np%firstChild)) then
          call destroyNode(np%firstChild)
       endif
       if (associated(np%attributes)) call destroyNamedNodeMap(np%attributes)
       call unstring(np%nodeName)
       call unstring(np%nodeValue)
       if (associated(np%previousSibling)) & 
                np%previousSibling%nextSibling => np%nextSibling
       if (associated(np%nextSibling)) & 
                np%nextSibling%previousSibling => np%previousSibling
       if (associated(np%parentNode)) then
          if (associated(np%parentNode%firstChild,np)) &
               np%parentNode%firstChild => null()
          if (associated(np%parentNode%lastChild,np)) &
               np%parentNode%lastChild => null()
       endif
       if (associated(np,node)) then    
          deallocate(np)
          allocated_nodes = allocated_nodes - 1
          EXIT                           ! do not destroy siblings
       else
          ghost => np
          np => np%nextSibling
          deallocate(ghost)
          allocated_nodes = allocated_nodes - 1
       endif
    enddo
    node => null()     ! superfluous ?
  end subroutine destroyNode

end module m_dom_types

