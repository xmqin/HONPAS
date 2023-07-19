module corresponding_node
  use flib_dom
  use string_utilities, only: clean_string, only_numbers

  logical :: corresponding_debug = .false.

  public :: find_corresponding_node
  type(fnode),pointer, public :: modified_corresponding_last => null()

  private
contains

  !--------------------------------------------------------------------------

  recursive function find_corresponding_node(ref,mod_node) result (corresponding)
    ! Given a reference node this subroutine finds the corresponding node
    ! in the modified list.
    type(fnode), pointer     :: ref  !The reference node 
    !type(fnode), pointer     :: mod  !The list of modified nodes.
    type(fnode), pointer     :: mod_node ! The modified node corresponding to the reference node.

    !Internal vars.
    integer                  :: i,n_elem
    type(fnode), pointer     :: node
    type(fnodelist),pointer  :: mod_list
    logical                  :: corresponding 

    corresponding = .false.

    if (.not. associated(modified_corresponding_last) )then
       stop "It was impossible to find a corresponding node"
    endif
    nullify(mod_node)

    !Get the tag and print nodes properties
 
    if(corresponding_debug)then
       print*, "********************************************************"
       print*, "Looking for node corresponding to:"
       call node_properties(ref)
    endif
    
    corresponding = is_corresponding_node(ref,modified_corresponding_last)
    if (corresponding) then
       mod_node => modified_corresponding_last       
       return
    endif
    
    mod_list=>getChildNodes(modified_corresponding_last)
    n_elem=getLength(mod_list)
    
    !Loop over the elements of the list and compare the nodes
    !until we found the corresponding one.
    if (.not. corresponding) then
       do i=0,n_elem-1
          node => item(mod_list,i)
          corresponding = is_corresponding_node(ref,node)  
          if (corresponding) then
             mod_node => node
             return
          endif
       enddo
    endif

    if (.not. corresponding) then
       !We couldn't find a corresponding node, so
       !Set the parent of modified_corresponding_last as modified_corresponding_last
       if(corresponding_debug) then
          print *,"Couldn't find corresponding, updating modified_corresponding_last"
       endif
       modified_corresponding_last => getParentNode(modified_corresponding_last)
       !Look for corresponding again!
       corresponding = find_corresponding_node(ref, mod_node)
    endif

  end function find_corresponding_node

  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------

  recursive function is_corresponding_node(ref,mod) result (same)
    use flib_dom

    type(fnode),pointer :: ref
    type(fnode),pointer :: mod
    logical :: same

    integer :: i,n_type_ref, n_type_mod,length_ref,length_mod, length
    type(string) :: sn_ref,sn_mod,sv_ref, sv_mod
    type(string) :: attr_name_ref, attr_value_ref, attr_value_mod, attr_name_mod
    type(fnamedNodeMap),pointer :: attr_ref,attr_mod
    type(fnode), pointer ::  c_node_ref, c_node_mod

    same = .true.

    !Check name
    sn_ref= getNodeName(ref)
    sn_mod= getNodeName(mod)

    if (char(sn_ref) /= "" .and. char(sn_mod) /= "") then
       sn_ref = clean_string(sn_ref)
       sn_mod = clean_string(sn_mod)
    endif

    if (char(sn_ref) /= char(sn_mod)) same = .false.
    if (corresponding_debug)  print *, "   ///////////////////////////////////is_corresponding-begin"
    if (corresponding_debug)  print *, "     is_corresponding: checking nodes with names=",char(sn_ref),"|",char(sn_mod)


    !Check node type
    n_type_ref = getNodeType(ref)
    n_type_mod = getNodeType(mod)

    if (n_type_ref /= n_type_mod ) same = .false.
    if (corresponding_debug) print *,"     is_corresponding: node types=",n_type_ref,"|",n_type_mod

    !Compare values if they aren't numeric
    sv_ref = getNodeValue(ref)
    sv_mod = getNodeValue(mod)

    if (len(sv_ref) > 0 .and. len(sv_mod) > 0)then

       if ( .not. only_numbers(sv_ref) .and. .not. only_numbers(sv_mod))then
          if(corresponding_debug) print *,"    is_corresponding: values=",char(sv_ref),char(sv_mod)
          if (sv_ref /= sv_mod)  same = .false.
       endif
    endif


    if (same .and. hasAttributes(ref) .and. hasAttributes(mod)) then

       attr_ref => getAttributes(ref)
       attr_mod => getAttributes(mod)

       length_ref = getlength(attr_ref)
       length_mod = getlength(attr_mod)
       length     = min(length_ref,length_mod)

       if (corresponding_debug) print *, "     is_corresponding: legths=",length_ref,length_mod

       if (length_ref /= length_mod) same = .false.

       if (same .and. length_ref > 0)then
          do i=0,length-1
             if(corresponding_debug) print *,"     is_corresponding: i attr=",i
             c_node_ref => item(attr_ref,i)
             c_node_mod => item(attr_mod,i)           

             !Check attributes names and values.
             attr_name_ref = getNodeName(c_node_ref)
             attr_name_mod = getNodeName(c_node_mod)

             if(corresponding_debug) print *,"       is_corresponding: attr names=",&
        trim(char(attr_name_ref)),"|",trim(char(attr_name_mod))

             if (attr_name_ref /= attr_name_mod) then
                same = .false.
                exit
             endif

             attr_value_ref = getNodeValue(c_node_ref)
             attr_value_mod = getNodeValue(c_node_mod)

             if(corresponding_debug) print *,"       is_corresponding: attr values=",&
         char(attr_value_ref),"|",char(attr_value_mod)

             if (attr_value_ref /= attr_value_mod) then
                same = .false.
                exit
             endif

          enddo
       endif
    endif
   
    if(corresponding_debug)  print *, "   /////////////////////////////////////////////is_corresponding-end"

    if(corresponding_debug)print *,"     is_corresponding_node: same?",same

  end function is_corresponding_node

  !----------------------------------------------------------------------

  subroutine node_properties(node)
    type(fnode),pointer :: node

    type(string) :: name, value
    type(string) :: attr_name, attr_value
    type(fnamedNodeMap),pointer :: attr
    type(fnode), pointer ::  sub_node
    integer      :: i,length

    name = getNodeName(node)
    if (len(name)>0) print *,"   Node_properties: name=",trim(char(name))
    value = getNodeValue(node)
    if (len(value)>0) print *,"   Node_properties: value=",trim(char(value))

    attr => getAttributes(node)
    length = getLength(attr)

    do i=0,length-1
       print *,"   Node_properties: subnode i=",i
       sub_node => item(attr,i)
       attr_name = getNodeName(sub_node)
       if (len(attr_name)>0) print *,"   Node_properties: subnode name=",trim(char(attr_name))
       attr_value = getNodeValue(sub_node)
       if (len(attr_value)>0) print *,"   Node_properties: subnode value=",trim(char(attr_value))
    enddo 

  end subroutine node_properties
end module corresponding_node
