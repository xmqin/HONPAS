! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      module m_rhofft

      use precision

      integer, parameter, public  :: FORWARD = -1
      integer, parameter, public  :: BACKWARD = +1

      public :: rhofft

      CONTAINS

      subroutine rhofft( CELL, N1, N2, N3, Mesh, nspin, RHO, rhog,
     $                   direction)

C *********************************************************************
C WRITTEN BY J.M.SOLER, JUNE 1995.
C **************** INPUT **********************************************
C REAL*8  CELL(3,3)     : UNIT CELL VECTORS
C INTEGER N1,N2,N3      : NUMBER OF LOCAL MESH DIVISIONS IN EACH CELL VECTOR
C INTEGER Mesh(3)       : Number of global mesh divisions
C REAL*4  RHO(N1,N2,N3) : DENSITIY AT MESH POINTS
C *********************************************************************

C     Modules
      use precision,   only : dp, grid_p
      use sys,         only : die
      use m_fft,       only : fft     ! 3-D fast Fourier transform
      use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
      use cellsubs,    only : volcel  ! Finds unit cell volume

      implicit          none

C     Input/output variables
      integer               :: N1, N2, N3, Mesh(3), nspin
      integer, intent(in)   :: direction
      real(grid_p)          :: RHO(N1*N2*N3,nspin)
      real(dp)              :: CELL(3,3)
      real(grid_p)          :: rhog(2,N1*N2*N3,nspin)

C     Local variables
      integer               :: ispin, ng
      real(dp)              :: VOLUME, b(3,3)

#ifdef DEBUG
      call write_debug( '    PRE rhofft' )
#endif

C     Start time counter
      call timer( 'rhofft', 1 )

C     Find unit cell volume
      VOLUME = VOLCEL( CELL )

C     Find reciprocal lattice vectors
      call reclat(CELL, B, 1 )

      NG = mesh(1) * mesh(2) * mesh(3)

      if (direction == FORWARD) then
         rhog(1,:,:) = RHO(:,:)
         rhog(2,:,:) = 0.0_grid_p

         do ispin = 1, nspin
            call fft( rhog(1,1,ispin), Mesh, -1 )
         enddo
         ! Units: electrons per unit cell
         rhog(:,:,:) = rhog(:,:,:) * volume / dble(ng)
      else if (direction == BACKWARD) then
         do ispin = 1, nspin
            call fft( rhog(1,1,ispin), Mesh, +1 )
         enddo
         rho(:,:) = rhog(1,:,:) * ng / volume
         !  print *, "sum(rho)*dvol", sum(rho)*volume/dble(ng)
      else
         call die("wrong fft direction")
      endif

C     Stop time counter
      call timer( 'rhofft', 2 )

#ifdef DEBUG
      call write_debug( '    POS rhofft' )
#endif

      end subroutine rhofft

      end module m_rhofft
