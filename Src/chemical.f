! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module chemical

      use sys
      use precision, only: dp

      implicit none

      private

      public :: atomic_number
      public :: number_of_species, species_label
      public :: is_floating, is_bessel, is_synthetic
      public :: read_chemical_types, print_chemical_type

      ! Public due to Bcast routines
      public :: chemical_types, chemical_list

      ! Species information
      type chemical_types
         integer                    :: no_of_species
         character(len=20), pointer :: spec_label(:)
         integer, pointer           :: z(:)
      end type chemical_types

      type(chemical_types), save :: chemical_list


      CONTAINS

      subroutine check(i)
      integer, intent(in) :: i
      if (i.lt.0 .or. i.gt.chemical_list%no_of_species)
     $     call die("Wrong species number requested")
      end subroutine check

      function number_of_species()
      integer number_of_species
      number_of_species = chemical_list%no_of_species
      end function number_of_species

      function species_label(i)
      character(len=20) species_label
      integer, intent(in)  :: i

      call check(i)
      species_label = chemical_list%spec_label(i)
      end function species_label

      function atomic_number(i)
      integer atomic_number
      integer, intent(in)  :: i

      call check(i)
      atomic_number = chemical_list%z(i)
      end function atomic_number
! -------
      function is_floating(i)
      logical is_floating
      integer, intent(in)  :: i

      call check(i)
      is_floating = (chemical_list%z(i) .le. 0)
      end function is_floating
! -------

      function is_bessel(i)
      logical is_bessel
      integer, intent(in)  :: i

      call check(i)
      is_bessel = (chemical_list%z(i) .eq. -100)
      end function is_bessel
! -------
!     Checks whether we are dealing with a synthetic atom
!
      function is_synthetic(i)
      logical is_synthetic
      integer, intent(in)  :: i

      call check(i)
      ! Note that we could have a synthetic ghost atom, with
      ! z <= -200
      is_synthetic = (abs(chemical_list%z(i)) .gt. 200)
      end function is_synthetic
!---
      subroutine read_chemical_types(silent)

      use parallel,    only : Node
      use fdf

      logical, intent(in), optional :: silent ! default .false.

      integer nsp, isp
      integer ns_read

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      character(len=20) :: label
      character(len=256) :: msg

      integer :: z, is
      logical :: found
      logical :: lsilent

      ! Determine whether we should be silent
      lsilent = .false.
      if ( present(silent) ) lsilent = silent
      if ( Node /= 0 ) lsilent = .true.

      ! Default to 0
      nsp = fdf_integer('Number_of_species',0)

      ! The most important thing to find is
      ! the block containing the species
      found = fdf_block('Chemical_species_label', bfdf)
      if (.not. found )
     $     call die("Block Chemical_species_label does not exist.")

      if ( nsp == 0 ) then
         ns_read = fdf_block_linecount('Chemical_species_label', 'iin')
      else
         ns_read = nsp
      end if
      ! If they are not equal we notify the user
      if ( nsp /= ns_read ) then
         nsp = ns_read
      end if
      if ( nsp == 0 ) call die("No species found!!!")

      allocate(chemical_list%spec_label(nsp))
      allocate(chemical_list%z(nsp))
      chemical_list%no_of_species = nsp

      ns_read = 0
      do while( fdf_bline(bfdf,pline) )
         if ( .not. fdf_bmatch(pline,'iin') ) cycle

         ns_read = ns_read + 1
         
         ! Get species information
         isp = fdf_bintegers(pline,1)
         label = fdf_bnames(pline,1)
         z = fdf_bintegers(pline,2)

         ! We cannot test label names in this
         ! loop as isp may be non-linear
         if ( isp < 1 .or. nsp < isp )
     $ call die("Wrong specnum in Chemical_species_label")

         chemical_list%z(isp) = z
         chemical_list%spec_label(isp) = label
        
      end do
      if ( ns_read /= nsp )
     &     call die("Not enough species in block")

      if ( .not. lsilent ) then
         ! Align output, always
         do isp = 1 , nsp
            call print_chemical_type(isp)
         end do
         write(*,*) ! new-line
      end if

      ! Check that none of the chemical species are the
      ! same
      if ( nsp > 1 ) then
       do z = 1 , nsp - 1
        do is = z+1, nsp
         if (trim(species_label(z))==trim(species_label(is))) then
          write(msg,'(2(a,i0,2a),2a)')
     & "Specie index/label = ",
     & z,'/',trim(species_label(z)), " has same label as ",
     & is,'/',trim(species_label(is)),". ",
     & " Use a different one for hygienic reasons."
          call die(trim(msg))
         end if
        end do
       end do
      end if

      end subroutine read_chemical_types

      subroutine print_chemical_type(isp)
      integer, intent(in)  :: isp
      
      character(len=256) :: label
      integer :: z
      logical :: floating, bessel, synthetic

      label = species_label(isp)
      z = atomic_number(isp)

      floating  =     z <= 0
      bessel    =     z == -100
      synthetic = abs(z) > 200

      ! Bessel *must* be checked first!
      if ( bessel ) then
         write(6,'(a,i3,3a)') 'Species number: ',isp,
     $        ' Label: ', trim(label),
     $        ' (floating Bessel functions)'
      else if ( floating ) then
         write(6,'(a,i3,a,i4,3a)') 'Species number: ',isp,
     $        ' Atomic number: ', z,
     $        ' Label: ', trim(label),
     $        ' (floating PAOs)'
      else
         write(6,'(a,i3,a,i4,2a)') 'Species number: ',isp,
     $        ' Atomic number: ', z,
     $        ' Label: ', trim(label)
      end if
      end subroutine print_chemical_type

      end module chemical
