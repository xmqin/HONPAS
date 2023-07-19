! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
program pdosxml
!
! Driver for Siesta PDOS processing
!
use f2kcli      ! Remove if it gives you trouble (see below)
use flib_sax
use m_pdos      ! Defines begin_element, end_element, pcdata_chunk

integer      :: iostat, i
type(xml_t)  :: fxml

!
! If you have trouble with the command-line processor f2kcli, remove the
! lines from here...
!
integer             :: narg
character(len=50)   :: filein

narg = command_argument_count()
if (narg == 1) then
   call get_command_argument(1,filein)
else
   stop "Usage: pdos filename"
endif
!  ... to here  AND  change the following line
!
call  open_xmlfile(filein,fxml,iostat)
! to this:
!call  open_xmlfile("PDOS.INPUT",fxml,iostat)
!
! and make sure your input file is renamed (or linked) to PDOS.INPUT

if (iostat /=0) stop "Cannot open file"

call xml_parse(fxml, &
                begin_element,end_element,pcdata_chunk,verbose=.false.)

do i=1, n_energies
   write(6,"(3f12.5)") energies(i), dos1(i), dos2(i)
enddo


end program pdosxml














