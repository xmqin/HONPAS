! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module planed
      use precision
      public :: plane
      private

      CONTAINS

      SUBROUTINE PLANE( NA, NPLAMAX, IDIMEN, OPTION, 
     .                  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, 
     .                  NPX, NPY, NPZ, COORPO, NORMALP, 
     .                  DIRVER1, DIRVER2,
     .                  XA, NAPLA, INDICES, ISCALE,
     .                  LATPOINT, PLAPOINT, XAPLANE )

C **********************************************************************
C This subroutine calculates coordinates of the points of a plane 
C or in the 3D-grid in two reference frames: in the lattice reference 
C frame (lrf) and in another one centered on the plane in which the 
C third coordiante is always cero (prf)
C **********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .   NPLAMAX, OPTION, NPX, NPY, NPZ, ISCALE, NA, NAPLA, INDICES(NA),
     .   IDIMEN

      real(dp), INTENT(IN) ::
     .   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .   NORMALP(3), COORPO(3,3), XA(3,NA)

      real(dp), INTENT(IN) ::  DIRVER1(3), DIRVER2(3)

      real(dp), INTENT(OUT) ::
     .   PLAPOINT(NPLAMAX,3), LATPOINT(NPLAMAX,3),
     .   XAPLANE(3,NA)

C **** INPUTS **********************************************************
C INTEGER NA         : Number of atoms in supercell
C INTEGER NPLAMAX    : Maximum number of points in the plane
C INTEGER IDIMEN     : Specify if the run is to plot quantities
C                      in a plane or in a 3D grid (2 or 3, respect)
C INTEGER OPTION     : Choose one of these options to define a plane:
C                      1 = components of the normal vector
C                      2 = equation of two straigth lines
C                      3 = coordinates of three points
C                      4 = indices of three atoms
C INTEGER NPX,NPY,NPZ: Number of points generated along x and y
C                      directions (and z for 3D-grids) 
C REAL(DP) XMIN, XMAX  : Limits of the plane in the PRF for x-direction
C REAL(DP) YMIN, YMAX  : Limits of the plane in the PRF for y-direction
C REAL(DP) ZMIN, ZMAX  : Limits of the plane in the PRF for z-direction
C REAL(DP)  NORMAL(3)  : Components of the normal vector used to define
C                      the plane
C REAL(DP)  DIRVER1(3) : Components of the first vector contained
C                      in the plane
C                      (Only used if ioption = 2)
C REAL(DP)  DIRVER2(3) : Components of the second vector contained
C                       in the plane
C                       (Only used if ioption = 2)
C REAL(DP)  COORPO(3,3): Coordinates of the three points used to define
C                      the plane (Only used if ioption = 3)
C INTEGER NAPLA      : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA): Indices of the atoms whose coordinates will
C                      be rotated from the lattice reference frame
C                      to the in-plane reference frame
C INTEGER ISCALE     : Unit of the points of the plane
C REAL(DP)  XA(3,NA)   : Atomic positions in cartesian coordinates
C                      (in bohr)
C **** OUTPUT **********************************************************
C REAL(DP)  LATPOINT(NPLAMAX,3) : Coordinates of the points 
C                               (lattice ref. frame)
C REAL(DP)  PLAPOINT(NPLAMAX,3) : Coordinates of the points 
C                               (plane refer. frame)
C REAL(DP)  XAPLANE(3,NA): Atomic coordinates in plane reference frame
C **********************************************************************

C Internal variables ---------------------------------------------------

      integer i,j
      real(dp) length,modnor,moddi1
      real(dp) normal(3)
      real(dp) xplane(3),yplane(3)
      real(dp) origin(3)
      real(dp) mrot(3,3),inmrot(3,3)
      real(dp) :: tmpvec1(3), tmpvec2(3)

      external length, atompla


C **********************************************************************
C REAL(DP) MODNOR, MODDI1       : Length of the vector normal and dirver1
C REAL(DP) XPLANE(3), YPLANE(3) 
C        ZPLANE(3)            : Plane reference frame
C REAL(DP) ORIGIN(3)            : Coordinates of a point of the plane 
C                               that will be considered the origin of the 
C                               reference frame fixed in the plane. 
C REAL(DP) MROT(3,3)            : matrix that relates the in plane
C                               reference frame with the lattice reference
C                               frame.
C REAL(DP) INMROT(3,3)          : Inverse of the previous matrix
C **********************************************************************

      if(option .eq. 1) then

         do i = 1,3
            xplane(i) = coorpo(2,i) - coorpo(1,i)
         enddo

         modnor = length(normalp)
         moddi1 = length(xplane)

         do i = 1,3 
            normal(i) = (1/modnor)*normalp(i)
            xplane(i) = (1/moddi1)*xplane(i)
         enddo

c         write(*,*)(normal(i),i=1,3)
c         write(*,*)(xplane(i),i=1,3)
*=======================================================================
      else if(option .eq. 2) then     

         call crossv(dirver1,dirver2,normal)
         
         modnor = length(normal)
         moddi1 = length(dirver1)

         do i = 1,3
            normal(i) = (1/modnor) * normal(i)
            xplane(i) = (1/moddi1) * dirver1(i)
         enddo

*=======================================================================
      else if( (option .eq. 3) .or. (option .eq. 4) ) then

         do i = 1,3
            tmpvec1(i) = coorpo(2,i) - coorpo(1,i)
            tmpvec2(i) = coorpo(3,i) - coorpo(1,i)
         enddo
        
         call crossv(tmpvec1,tmpvec2,normal)
        
         modnor = length(normal)
         moddi1 = length(tmpvec1)
        
         do i = 1,3
            normal(i) = (1/modnor)*normal(i)
            xplane(i) = (1/moddi1)*tmpvec1(i)
         enddo
          
      endif

      call crossv(normal,xplane,yplane)

      do i = 1,3
         origin(i) = coorpo(1,i)
      enddo


*
*       We define now the matrix that describe the rotation.
*
      do j = 1,3
         mrot(j,1) = xplane(j) 
         mrot(j,2) = yplane(j)
         mrot(j,3) = normal(j)
      enddo
        
      do i = 1,3
         do j = 1,3
            inmrot(i,j) = mrot(j,i)
         enddo
      enddo

      call popla(nplamax,npx,npy,npz,xmin,xmax,ymin,ymax,zmin,zmax,
     .          idimen,plapoint)
        
      call rotation(nplamax,npx,npy,npz,mrot,inmrot,origin,
     .              plapoint,latpoint)

      call atompla( na, origin, xa, inmrot, idimen, iscale, 
     .              napla, indices, xaplane )

      end subroutine plane


*=======================================================================
*=======================================================================
      subroutine popla(nplamax,npx,npy,npz,xmin,xmax,ymin,ymax,
     .           zmin,zmax,idimen,plapoint)
*
*     This subroutine generates the coordinates of the points of the plane in
*    a reference frame fixed in the plane and where the z-axis is the normal
*    vector (in this way the third coordinate of the points is always zero).
*
*-----------------------------------------------------------------------
*     VARIABLES
*
        implicit none
    
        integer nplamax
        integer npx,npy,npz,i,k,l,m,idimen
        real(dp) xmin,xmax,ymin,ymax,zmin,zmax,deltax,deltay,deltaz
        real(dp) plapoint(nplamax,3)
*-----------------------------------------------------------------------
*     INICIALIZATION
*
        deltax=(xmax-xmin)/(npx - 1)
        deltay=(ymax-ymin)/(npy - 1)
        if (npz .ne. 1) then
          deltaz=(zmax-zmin)/(npz - 1)
        else
          deltaz=0.0d0
        endif

*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCK
*

C Do the loops differently if 2D or 3D. 
C This is necessary to comply with the ordering
C expected by Cube, and to maintain the former 
C ordering in the 2D case from former versions

        i = 1
        if (idimen .eq. 2) then
          do k = 1,npx
            do l = 1,npy
               plapoint(i,1) = xmin + (k-1)*deltax
               plapoint(i,2) = ymin + (l-1)*deltay
               plapoint(i,3) = 0.0d0
               i = i+1
            enddo
          enddo
        else if (idimen .eq. 3) then
          do m = 1,npz
             do l = 1,npy
               do k = 1,npx
                  plapoint(i,1) = xmin + (k-1)*deltax
                  plapoint(i,2) = ymin + (l-1)*deltay
                  plapoint(i,3) = zmin + (m-1)*deltaz
                  i = i+1
                enddo
            enddo
          enddo
        endif

        end subroutine popla
*=======================================================================
*=======================================================================       

      subroutine rotation(nplamax,npx,npy,npz,mrot,inmrot,origin,
     .                    plapoint,latpoint)
*
*     This subroutine makes the transformation of the coordinates of the points
*    of the plane in the plane reference frame to the lattice reference frame.
*    
*-----------------------------------------------------------------------
*     VARIABLES
*
      implicit none

      integer nplamax
      integer npx,npy,npz,i,j
      real(dp) mrot(3,3),inmrot(3,3)
      real(dp) origin(3),vaux1(3),vaux2(3)
      real(dp) plapoint(nplamax,3),latpoint(nplamax,3)

      external matvect
*    
*     vaux1,vaux2 = auxiliar vectors needed in the computation of latpoint.
*-----------------------------------------------------------------------
*     COMPUTATIONAL BLOCK
*

      do i = 1, npx*npy*npz
         do j = 1,3
           vaux1(j) = plapoint(i,j)
         enddo

         call matvect(mrot,vaux1,vaux2)

         do j = 1,3
            latpoint(i,j) = vaux2(j) + origin(j)
         enddo
      enddo

      end subroutine rotation
*=======================================================================
*=======================================================================
      subroutine dotv(a,b,c)
*
*     This subroutine calculates the dot product of two vectors defined
*    in cartesian coordinates.
*
      real(dp) a(3),b(3)
      real(dp) c
*
*     a(3), b(3) : cartesian coordintes of the two vectors.
*     c : dot product of the two vectors.
*

      c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      end subroutine dotv
*=======================================================================
*=======================================================================
      SUBROUTINE CROSSV( A, B, C )

C **********************************************************************
C This subroutine calculates the cross product of two vectors defined 
C in cartesian coordinates.
C **********************************************************************

      real(dp), INTENT(IN) ::
     .   A(3), B(3)

      real(dp), INTENT(OUT) ::
     .   C(3)

C **** INPUT ***********************************************************
C REAL(DP) A(3)  : Cartesian components of the first vector
C REAL(DP) B(3)  : Cartesian components of the second vector
C **** OUTPUT **********************************************************
C REAL(DP) C(3)  : Cartesian components of cross product vector
C **********************************************************************

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      END SUBROUTINE CROSSV

      end module planed
