! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_proximity_check
      private
      public :: proximity_check
      CONTAINS

      subroutine proximity_check( rmax )
      USE precision, only:dp
      USE siesta_options
      use units, only: Ang
      use siesta_geom
      use chemical, only: is_floating
      use parallel,     only : IOnode
      use neighbour,    only : jna=>jan, r2ij, mneighb,
     &                         reset_neighbour_arrays
      use m_ts_global_vars, only : onlyS
      implicit none
      real(dp), intent(in)  :: rmax
      integer :: ii, jj, jamin, idxneighJ, nneigbI, isel, nna
      real(dp):: rmin, r2min 

!     Check if any two atoms are unreasonably close
      isel = 0   ! Get back all neighbors, regardless of whether ja>ia

!     Initialize neighb
      call mneighb( scell, rmax, na_s, xa, 0, isel, nna )
      do ii = 1,na_s
        ! Check whether the atom is a floating orbital
        ! If that is the case, simply skip checking...
        if ( is_floating(isa(ii)) ) cycle
        r2min = huge(1._dp)
        jamin = 0
        call mneighb( scell, rmax, na_s, xa, ii, isel, nneigbI )
        do jj = 1,nneigbI
          idxneighJ = jna(jj)
          if ( is_floating(isa(idxneighJ)) ) cycle
          if ( r2ij(jj).lt.r2min .and. idxneighJ.ge.ii ) then
!           Check that it is not the same atom
            if ( idxneighJ.ne.ii .or. r2ij(jj).gt.1.d-12 ) then
              r2min = r2ij(jj)
              jamin = idxneighJ
            endif
          endif
        enddo
        rmin = sqrt( r2min )
        if (IOnode) then
          if (( rmin .lt. rijmin ) .and. (.not. onlyS))
     &      write(6,'(a,2i6,a,f12.6,a)') 'siesta: WARNING: Atoms', ii,
     &                    jamin, ' too close: rij =', rmin/Ang, ' Ang'
        endif
      enddo

      call reset_neighbour_arrays( )

      end subroutine proximity_check
      end module m_proximity_check
