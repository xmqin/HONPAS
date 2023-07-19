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
      module chemical

      use sys
      use precision, only: dp

      implicit none

      private
      public number_of_species, atomic_number, species_label
      public is_floating, is_bessel, is_synthetic
      public read_chemical_types, print_chemical_type
!
!     Species information
!
      type chemical_types
         integer                       :: no_of_species
         character(len=20), pointer    :: spec_label(:)
         integer, pointer              :: z(:)
      end type chemical_types

      type(chemical_types), save       :: chemical_list


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
      is_synthetic = (chemical_list%z(i) .gt. 200)
      end function is_synthetic
!---
      subroutine read_chemical_types(silent)

      use parse
      use fdf

      logical, intent(in), optional :: silent

      integer nsp, isp
      integer ns_read

      type(block), pointer  :: bp
      type(parsed_line), pointer  :: p
      character(len=132) line

      character(len=20) label
      integer  z, is
      logical floating, bessel, found, synthetic
      logical printing_allowed

      printing_allowed = .true.
      if (present(silent)) then
         if (silent) printing_allowed = .false.
      endif

      nsp = fdf_integer('Number_of_species',0)
      if (nsp.eq.0) call die("No species found!!!")

      allocate(chemical_list%spec_label(nsp))
      allocate(chemical_list%z(nsp))
      chemical_list%no_of_species = nsp
     
      nullify(bp)
      found = fdf_block('Chemical_species_label',bp)
      if (.not. found )
     $     call die("Block Chemical_species_label does not exist.")

      ns_read = 0
      loop: DO
        if (.not. fdf_bline(bp,line)) exit loop
        ns_read = ns_read + 1
        p => digest(line)
        if (.not. match(p,"iin"))
     $       call die("Wrong format in Chemical_species_label")
        isp = integers(p,1)
        if (isp .gt. nsp .or. isp .lt. 1)
     $       call die("Wrong specnum in Chemical_species_label")
        label = names(p,1)
        z = integers(p,2)
        floating = (z.le.0)
        bessel = (z.eq.-100)
        synthetic = (z.gt.200)

        if (printing_allowed) then
           if (.not.floating) then
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' Atomic number:',  z
           elseif (bessel) then
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' (floating Bessel functions)'
           else
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' Atomic number:',  z,
     $             ' (floating PAOs)'
           endif
        endif
!
        if (ns_read > 1) then
          do is = 1, ns_read-1
            if (trim(chemical_list%spec_label(is)) == trim(label)) then
               call die("Species label " // trim(label) // " repeated."
     $              // " Use a different one for hygienic reasons.")
            endif
         enddo
        endif

        chemical_list%z(isp) = z
        chemical_list%spec_label(isp) = label
 
       call destroy(p)
      enddo loop
      if (ns_read .ne. nsp) call die("Not enough species in block")
      call destroy(bp)

      end subroutine read_chemical_types

      subroutine print_chemical_type(isp)
      integer, intent(in)  :: isp
      
      character(len=256) :: label
      real(dp)           :: z
      logical            :: floating, bessel, synthetic

      label = species_label(isp)
      z = atomic_number(isp)

        floating = (z.le.0)
        bessel = (z.eq.-100)
        synthetic = (z.gt.200)

        if (.not.floating) then
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' Atomic number:',  z
           elseif (bessel) then
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' (floating Bessel functions)'
           else
              write(6,*) 'Species number: ', isp,
     $             ' Label: ', trim(label),
     $             ' Atomic number:',  z,
     $             ' (floating PAOs)'
       endif
      end subroutine print_chemical_type

      end module chemical






