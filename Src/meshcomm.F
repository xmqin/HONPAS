! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module MeshComm
C     Modules
      use precision
      implicit none

      integer,          save :: MCsize          ! Number of communications
      integer, pointer, save :: MCBuffer(:)     ! Buffer that contains all
                                                ! the data
      integer, pointer, save :: MCsrc(:)        ! Source Process ID
      integer, pointer, save :: MCdst(:)        ! Destiny Process ID
      integer, pointer, save :: MCinZ(:)        ! initial index on Z
      integer, pointer, save :: MCfiZ(:)        ! Final index on Z
      integer, pointer, save :: MCinY(:)        ! initial index on Y
      integer, pointer, save :: MCfiY(:)        ! Final index on Y
      integer,          save :: udMeshLim(2,3)  ! Lim. of uniform dist.
      integer,          save :: nddMeshLim(2,3) ! Lim. of nested disec. dist.

      CONTAINS

C ==================================================================
C Sets the mesh limits of the mesh and all the related mesh sizes
C to the Uniform distribution
C ==================================================================
C subroutine setUdMeshLim( nsm, nsp, nml, nmpl, ntml, ntpl )
C
C INPUT:
C integer nsm          : Number of mesh sub-divisions in each direction
C integer nsp          : Number of sub-points of each mesh point
C
C OUTPUT:
C integer nml(3)       : Local Number of Mesh divisions of each cell vector
C integer nmpl         : Number of mesh points in unit cell locally
C integer ntml(3)      : Total number of mesh points stored locally
C integer ntpl         : Number of mesh Total Points in unit cell
C                      : (including subpoints) locally
C
C BEHAVIOR:
C Sets the uniform distribution as the one to be used. We also compute
C some mesh sizes.
C
C ==================================================================
      subroutine setUdMeshLim( nsm, nsp, nml, nmpl, ntml, ntpl )
      use mesh, only: meshLim
      implicit none
      integer, intent(in)  :: nsm, nsp
      integer, intent(out) :: nml(3), nmpl, ntml(3), ntpl
!------------------------------------------------------------------------- BEGIN
      MeshLim = nddMeshLim
      nml(1) = (MeshLim(2,1)-MeshLim(1,1)) + 1
      nml(2) = (MeshLim(2,2)-MeshLim(1,2)) + 1
      nml(3) = (MeshLim(2,3)-MeshLim(1,3)) + 1
      nmpl   = nml(1)*nml(2)*nml(3)
      ntml   = nml*nsm
      ntpl   = nmpl*nsp
!--------------------------------------------------------------------------- END
      end subroutine setUdMeshLim

C ==================================================================
C Sets the mesh limits of the mesh and all the related mesh sizes
C to the Nested Disection distribution
C ==================================================================
C subroutine setNddMeshLim( nsm, nsp, nml, nmpl, ntml, ntpl )
C
C INPUT:
C integer nsm          : Number of mesh sub-divisions in each direction
C integer nsp          : Number of sub-points of each mesh point
C
C OUTPUT:
C integer nml(3)       : Local Number of Mesh divisions of each cell vector
C integer nmpl         : Number of mesh points in unit cell locally
C integer ntml(3)      : Total number of mesh points stored locally
C integer ntpl         : Number of mesh Total Points in unit cell
C                      : (including subpoints) locally
C
C BEHAVIOR:
C Sets the nested disection distribution as the one to be used. We
C also compute some mesh sizes.
C
C ==================================================================
      subroutine setNddMeshLim( nsm, nsp, nml, nmpl, ntml, ntpl )
      use mesh, only: meshLim
      implicit none
      integer, intent(in)  :: nsm, nsp
      integer, intent(out) :: nml(3), nmpl, ntml(3), ntpl
!------------------------------------------------------------------------- BEGIN
      MeshLim = udMeshLim
      nml(1) = (MeshLim(2,1)-MeshLim(1,1)) + 1
      nml(2) = (MeshLim(2,2)-MeshLim(1,2)) + 1
      nml(3) = (MeshLim(2,3)-MeshLim(1,3)) + 1
      nmpl   = nml(1)*nml(2)*nml(3)
      ntml   = nml*nsm
      ntpl   = nmpl*nsp
!--------------------------------------------------------------------------- END
      end subroutine setNddMeshLim

      end module MeshComm
