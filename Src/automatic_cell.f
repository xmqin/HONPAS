! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      subroutine automatic_cell(ucell,scell,na_u,xa,isa,charnet)

      use precision, only: dp
      use atmfuncs,  only: rcut, nofis
      use parallel,  only: IONode
      use units,     only: Ang
      use atm_types,    only : nspecies

      implicit none

      real(dp), dimension(3,3), intent(out)   :: ucell
      real(dp), dimension(3,3), intent(out)   :: scell
      integer, intent(in)                     :: na_u
      real(dp), dimension(3,*), intent(in)    :: xa
      integer, dimension(*), intent(in)       :: isa
      real(dp), intent(in)                    :: charnet


      integer  :: ix, ia, is, iv, io
      real(dp) :: rc, xmin, xmax, rcmax

      ! We need a measure of the "extent of the atom" in terms of orbitals
      ! and Vna

      real(dp), allocatable  :: rcatom(:)

      allocate(rcatom(nspecies))

      ! Find the largest size for each species
      ! In previous versions rc(Vna) was used, but if the basis set
      ! contains large orbitals which are not occupied in the ground
      ! state, this is incorrect

      do is = 1, nspecies
         rcmax = rcut(is,0)   ! rc(Vna)
         do io = 1, nofis(is)
            rcmax = max(rcut(is,io),rcmax)
         enddo
         rcatom(is) = rcmax
      enddo

        ucell(1:3,1:3) = 0.0_dp
        scell(1:3,1:3) = 0.0_dp
        do ix = 1,3
          xmin =  huge(1._dp)
          xmax = -xmin
          do ia = 1,na_u
            is = isa(ia)
            rc = rcatom(is)
            xmin = min( xmin, xa(ix,ia)-rc )
            xmax = max( xmax, xa(ix,ia)+rc )
          enddo
!         Use a 10% margin for atomic movements
          ucell(ix,ix) = 1.10_dp * (xmax - xmin)
          scell(ix,ix) = ucell(ix,ix)
        enddo

        ! If the system is charged we build a cubic cell, 
        ! since that is the only case we know how to handle

        if (charnet .ne. 0.0_dp) then
          xmax = -huge(1._dp)
          do ix = 1,3
            if (ucell(ix,ix) .gt. xmax) xmax = ucell(ix,ix)
          enddo
          do ix = 1,3
            ucell(ix,ix) = xmax
            scell(ix,ix) = xmax
          enddo
        endif

        if (IOnode) then
          write(6,'(/,a,3(/,a,3f12.6))')
     .      'siesta: Automatic unit cell vectors (Ang):',
     .      ('siesta:', (ucell(ix,iv)/Ang,ix=1,3), iv =1,3)
        endif

        deallocate(rcatom)

        end subroutine automatic_cell

