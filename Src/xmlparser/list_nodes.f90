!Test the xml output of siesta 
!A xml file is compared against a  xml

program list_nodes 
use flib_dom
implicit none

type(fnode), pointer     :: f_parsed
type(fnodeList),pointer  :: lattice
type(fnode),pointer      :: sub_node
integer                  :: n_children,i


!read(5,*) input
f_parsed => parsefile("h2o.xml",verbose=.true.)

call list(f_parsed)

!Example with the latticevectors.
lattice => getElementsByTagName(f_parsed,"latticeVector")

n_children = getLength(lattice)
do i=0,n_children-1
        sub_node => item(lattice,i)
        call list(sub_node) 
enddo   



contains

recursive subroutine list(node)
use flib_dom
! Given two nodes lists this subroutine compares all the nodes, one by one.
type(fnode), pointer, intent(in) :: node

!Internal vars.
logical                  :: has_children
integer                  :: i,n_children
type(fnodelist),pointer  :: sub_list
type(fnode),pointer      :: sub_node

has_children = hasChildNodes(Node)

call print_node(node)

if ( has_children) then
        sub_list => getchildNodes(node)
        n_children = getLength(sub_list)
        do i=0,n_children-1
                sub_node => item(sub_list,i)
                call list(sub_node) 
        enddo   
endif

end subroutine list

subroutine print_node(node)
type(fnode),pointer, intent(in)  :: node

integer :: j,length
type(string) :: sn,sv
type(fnamedNodeMap),pointer :: attr
type(fnode), pointer :: c_node

sn = getNodeName(node)
print *, " Name:",char(sn)

j = getNodeType(node)
if ( j == TEXT_NODE) then
        print *, "   Text node"
elseif( j ==  ELEMENT_NODe) then
        print *, "   Element_node"
elseif( j == DOCUMENT_NODE)then
        print *, "   Document_node"
else
       print*, " other kind!"
       print*, " complete the list"
       stop
endif

sv = getNodevalue(node)
print *, "   Value:",char(sv)

if (hasAttributes(node))then
        print *, "   Attributes:"
        attr => getAttributes(node)
        length = getlength(attr)
        do j=0,length-1
                c_node => item(attr,j)
                sn = getNodeName(c_node)
                print ("(a,a20)"), "    attr. name:",char(sn)
                sv = getNodeValue(c_node)
                print ("(a,a20)") , "    attr. value:",char(sv)
                print *, "   ..............." 
        enddo
endif

print *,"-----------------------"

end subroutine print_node
end program list_nodes 



