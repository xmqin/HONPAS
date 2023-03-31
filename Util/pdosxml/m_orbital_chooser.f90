! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_orbital_chooser
!
! Determines which orbitals to consider when processing PDOS data
!
type, public :: orbital_id_t
   integer  :: n
   integer  :: l
   integer  :: m
   integer  :: z
   integer  :: index
   integer  :: atom_index
   character(len=40)  :: species
end type orbital_id_t

public :: want_orbital

CONTAINS

function want_orbital(orbid) result(wantit)
type(orbital_id_t), intent(in)   :: orbid
logical                          :: wantit

!
!  Examples
!  
!  1. Want only s-orbitals
!
!     wantit = ( orbid%l == 0 )
!
!  2. Want only n=3 orbitals
!
!     wantit = ( orbid%n == 3 )
!
!  2. Want 3p orbitals
!
!     wantit = ( ( orbid%n == 3 ) .and. (orbid%l == 0 ) )
!
!  3. Want Oxygen orbitals
!
!     wantit = ( orbid%species == "O" )
!
!wantit = .true.
!
      wantit = ( orbid%species == "O" )

end function want_orbital

end module m_orbital_chooser
