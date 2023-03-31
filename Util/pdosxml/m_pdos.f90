! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_pdos
!
use flib_sax
use m_orbital_chooser

private

!
! It defines the routines that are called from xml_parser in response
! to particular events.
!
public  :: begin_element, end_element, pcdata_chunk
private :: die

logical, private  :: in_energy_values = .false. 
logical, private  :: in_orbital = .false. , in_data = .false.
logical, private  :: in_nspin = .false.

logical, private  :: spin_polarized = .false.

integer, public, save  :: ndata, n_energies

integer, parameter, public    :: NMAX  = 4000

real, dimension(NMAX), public :: energies
real, dimension(NMAX), public :: dos1, dos2  , data 

type(orbital_id_t), private, save     ::  orbid

CONTAINS  !===========================================================

!----------------------------------------------------------------------
subroutine begin_element(name,attributes)
character(len=*), intent(in)    :: name
type(dictionary_t), intent(in)  :: attributes

character(len=100)  :: value
integer             :: status


select case(name)

      case ("pdos")

         dos1 = 0.0
         dos2 = 0.0

      case ("energy_values")
         in_energy_values = .true.

         n_energies = 0             ! To start the build up
         
         call get_value(attributes,"units",value,status)
         if (status /= 0 ) call die("Cannot determine energy units")

      case ("nspin")
         write(unit=0,fmt="(a)") "Found nspin element"
         in_nspin = .true.

      case ("orbital")
         in_orbital = .true.

         call get_value(attributes,"l",value,status)
         if (status /= 0 ) call die("Cannot determine l for orbital")
         read(unit=value,fmt=*) orbid%l

         call get_value(attributes,"n",value,status)
         if (status /= 0 ) call die("Cannot determine n for orbital")
         read(unit=value,fmt=*) orbid%n

         call get_value(attributes,"m",value,status)
         if (status /= 0 ) call die("Cannot determine m for orbital")
         read(unit=value,fmt=*) orbid%m

         call get_value(attributes,"z",value,status)
         if (status /= 0 ) call die("Cannot determine z for orbital")
         read(unit=value,fmt=*) orbid%z

         call get_value(attributes,"index",value,status)
         if (status /= 0 ) call die("Cannot determine orb index for orbital")
         read(unit=value,fmt=*) orbid%index

         call get_value(attributes,"atom_index",value,status)
         if (status /= 0 ) call die("Cannot determine atom index for orbital")
         read(unit=value,fmt=*) orbid%atom_index

         call get_value(attributes,"species",orbid%species,status)
         if (status /= 0 )   &
                 call die("Cannot determine atomic species for orbital")

      case ("data")
         in_data = .true.
         data = 0.0
         ndata = 0             ! To start the build up

end select

end subroutine begin_element
!----------------------------------------------------------------------

subroutine end_element(name)
character(len=*), intent(in)     :: name

select case(name)

      case ("orbital")
         in_orbital = .false.

      case ("energy_values")
         in_energy_values = .false.

      case ("nspin")
         in_nspin = .false.

      case ("data")
      !
      ! We are done filling up the data arrays
      ! Check that we got the right number of items
      !
         in_data = .false.
         if (spin_polarized) then
            if (ndata /= 2*n_energies) then
               print *, "ndata, n_energies", ndata, n_energies
               STOP "npts mismatch"
            endif
         else
            if (ndata /= n_energies) then
               print *, "ndata, n_energies", ndata, n_energies
               STOP "npts mismatch"
            endif
         endif

         if (want_orbital(orbid))  then
            write(unit=0,fmt="(a,6(i4),1x,a)")  &
               "Orbital:(n.l.m.z.index.atom_index.species):",  &
            orbid%n, orbid%l, orbid%m, orbid%z, orbid%index,   &
            orbid%atom_index,  trim(orbid%species)

            if (spin_polarized) then
               dos1(1:n_energies) = dos1(1:n_energies) + data(1:ndata:2)
               dos2(1:n_energies) = dos2(1:n_energies) + data(2:ndata:2)
            else
               dos1(1:n_energies) = dos1(1:n_energies) + data(1:ndata)
            endif
         endif

end select

end subroutine end_element
!----------------------------------------------------------------------

subroutine pcdata_chunk(chunk)
character(len=*), intent(in) :: chunk

integer :: nspin, ndata_spin
integer, dimension(10) :: nspin_dummy

if (len_trim(chunk) == 0) RETURN     ! skip empty chunk

if (in_data) then
!
   call build_data_array(chunk,data,ndata)
endif
if (in_energy_values) then
!
      call build_data_array(chunk,energies,n_energies)
endif

ndata_spin = 0
if (in_nspin) then
!
! This is the only safe way to read the numeric element content
!
   call build_data_array(chunk,nspin_dummy,ndata_spin)
   if (ndata_spin /= 1) &
      write(unit=0,fmt="(a)") &
        "Too many numbers in nspin element. Taking first one."
   nspin = nspin_dummy(1)
   spin_polarized = (nspin == 2)
   if (spin_polarized) then
      write(unit=0,fmt="(a)") "NOTE: Spin polarized system"
      write(unit=0,fmt="(a)") &
               "There will be two PDOS columns in the output file"
   endif

endif

end subroutine pcdata_chunk
!----------------------------------------------------------------------

      subroutine die(str)
      character(len=*), intent(in), optional   :: str
      if (present(str)) then
         write(unit=0,fmt="(a)") trim(str)
      endif
      write(unit=0,fmt="(a)") "Stopping Program"
      stop
      end subroutine die


end module m_pdos












