! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module iofa_m
  implicit none

  public :: iofa
contains
  subroutine iofa( na, fa , suffix )
    ! *******************************************************************
    ! Writes forces in eV/Ang
    ! Emilio Artacho, Feb. 1999
    ! Nick Papior, Jun 2015
    ! Updated precision Nick Papior, Jun 2018
    ! Changed for explicit suffix Nick Papior, September 2019
    ! ********** INPUT **************************************************
    ! integer na           : Number atoms
    ! real*8  fa(3,na)     : Forces on the atoms
    ! character(*) suffix  : suffix for the file
    ! *******************************************************************

    use files,     only : slabel, label_length
    use precision, only : dp
    use units,     only : Ang, eV

    implicit none

    integer,  intent(in) :: na
    real(dp), intent(in) :: fa(3,na)
    character(len=*),  intent(in) :: suffix

    external :: io_assign, io_close

    ! Internal 
    character(len=label_length+4) :: fname
    integer :: ia, iu, ix
    real(dp) :: tmp

    fname = trim(slabel) // '.' // trim(suffix)

    call io_assign( iu )
    open( iu, file=fname, form='formatted', status='unknown', position='rewind')      

    write(iu,'(i6)') na
    tmp = Ang / eV
    do ia = 1, na
      write(iu,'(i6,3(tr1,e17.9))') ia, fa(:,ia) * tmp
    end do

    call io_close( iu )

  end subroutine iofa

end module iofa_m
