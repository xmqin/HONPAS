! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
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














