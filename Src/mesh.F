! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module mesh
!
! Stores quantities that are connected with the mesh
!
      use precision, only : dp
      implicit none

! Mesh variables
! ----------------------------------------------------------------------
! real*8  cmesh(3,3)    : Mesh-cell vectors
! real*8  dxa(3,na)     : Atom position within mesh-cell
! integer idop(mop)     : Extended-mesh-index displacement of points
!                       : within a sphere of radius rmax
! integer indexp(nep)   : Translation from extended to normal mesh index
! integer ipa(na)       : Mesh cell in which atom is
! integer mop           : Maximum number of non-zero orbital points
! integer nmeshg(3)     : Total number of mesh points in each direction
! integer ne(3)         : Number of mesh-Extension intervals 
!                       : in each direction
! integer nem(3)        : Extended-mesh divisions in each direction
! integer nmsc(3)       : Mesh divisions of each supercell vector
! integer nmuc(3)       : Mesh divisions of each unit cell vector
! integer nsm           : Number of mesh sub-divisions in each direction
! integer nsp           : Number of sub-points of each mesh point
! integer nusc          : Number of unit cells in each supercell dir.
! real*8  rcmesh(3,3)   : Reciprocal mesh-cell vectors
!                       : (WITHOUT 2*pi factor)
! real*8  xdop(3,mop)   : Vector to mesh points within rmax
! real*8  xdsp(3,nsp)   : Vector to mesh sub-points
! integer meshLim(2,3)  : Upper an lower limits of the mesh in every process
! integer iatfold(3,na) : Supercell vector that keeps track of the 
!                         of the folding of the atomic coordinates in the mesh 
C ----------------------------------------------------------------------
      integer,  pointer, save :: idop(:)
      integer,  pointer, save :: ipa(:) => null()
      real(dp), pointer, save :: dxa(:,:) => null()
      real(dp), pointer, save :: xdop(:,:)
      real(dp), pointer, save :: xdsp(:,:) ! Allocated in dhscf
      integer,           save :: mop
      integer,           save :: ne(3)
      integer,           save :: nem(3)
      integer,           save :: nmsc(3)
      integer,save :: nmuc(3)       ! Mesh points in each unit cell vector
      integer,save :: nusc(3)       ! Unit cells in each supercell dir.
      integer,           save :: meshLim(2,3)  ! My processor's box 
                                               ! of mesh points:
                                               !   myBox(1,:) : lower bounds
                                               !   myBox(2,:) : upper bounds
      integer,           save :: nmeshg(3)
      integer,           save :: nsm
      integer,           save :: nsp
      real(dp),          save :: cmesh(3,3)
      real(dp),          save :: rcmesh(3,3)
      integer,  pointer, save :: indexp(:) => null()
      integer,  pointer, save :: iatfold(:,:) => null()
      end module mesh
