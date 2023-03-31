! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
!**********************************************************************
! module sorting
!
! Provides the following public procedures:
!   iorder  ! Orders elements of an integer array according to some index
!   ordix   ! Finds an order index of increasing array elements
!   order   ! Orders elements of a real array according to some index
!   ordvec  ! Orders a set of vectors in a well defined order
!**********************************************************************

      MODULE sorting

      CONTAINS

      SUBROUTINE ORDVEC( TOL, NX, NV, V, INDEX )
      
C **********************************************************************
C Orders a set of vectors in a well defined order: by last coordinate;
C if last coordinate is equal, by before-last corrdinate, etc.
C Written by J.Soler. October 1997.
C ************ INPUT ***************************************************
C REAL*8  TOL      : Tolerance to consider two coordinates equal.
C INTEGER NX       : Number of vector coordinates
C INTEGER NV       : Number of vectors
C REAL*8  V(NX,NV) : Vectors to be ordered
C ************ OUTPUT **************************************************
C REAL*8  V(NX,NV) : Ordered vectors
C INTEGER INDEX(NV): Order index such that Vout(IX,IV)=Vin(IX,INDEX(IV))
C **********************************************************************

      IMPLICIT          NONE
      INTEGER           NV, NX
      INTEGER           INDEX(NV)
      DOUBLE PRECISION  TOL, V(NX,NV)

      INTEGER           IV, IV0, IV1, IX, JX

      INTEGER :: IAUX(NV)           ! automatic array

      DO IV = 1,NV
        INDEX(IV) = IV
      ENDDO
      DO IX = NX,1,-1
        IV0 = 0
  30    CONTINUE
          DO IV1 = IV0+1,NV-1
            DO JX = IX+1,NX
              IF (ABS(V(JX,IV1+1)-V(JX,IV1)) .GT. TOL) GOTO 60
            ENDDO
          ENDDO
          IV1 = NV
  60      CONTINUE
          IF (IV1 .GT. IV0+1) THEN
!JMS 2008_03_02
!            CALL ORDIX(  V(IX:IX,IV0+1:), NX, IV1-IV0, IAUX )
            CALL ORDIX(  V(IX,IV0+1), NX, IV1-IV0, IAUX )
            CALL ORDER(  V(1,IV0+1),  NX, IV1-IV0, IAUX )
            CALL IORDER( INDEX(IV0+1), 1, IV1-IV0, IAUX )
          ENDIF
          IV0 = IV1
        IF (IV0 .LT. NV-1) GO TO 30
      ENDDO

      END SUBROUTINE ORDVEC


      SUBROUTINE ordix( x, m, n, indx )
! *******************************************************************
! SUBROUTINE ORDIX( X, M, N, INDX )
! Makes an index table of increasing array X, with size N and
! stride M using the heapsort algorithm.
! Written by J.M.Soler, May.2015
! *************** INPUT *********************************************
! REAL*8  X(M,N)   : Array with the values to be ordered
! INTEGER M, N     : Dimensions of array X
! *************** OUTPUT ********************************************
! INTEGER INDX(N)  : Array which gives the increasing order of X(1,I):
!                    X(1,INDX(I)) .LE. X(1,INDX(I+1)) )
! *************** USAGE *********************************************
! Example to order atomic positions X(I,IA), I=1,3, IA=1,NA by
! increasing z coordinate:
!    CALL ORDIX( X(3,1), 3, NA, INDEX )
!    CALL ORDER( X(1,1), 3, NA, INDEX )
! *************** ALGORITHM *****************************************
! A hierarchical 'family tree' (heap) is generated, with each parent
! k older than its two children (2*k and 2*k+1). Persons k with k>np
! are children (np = highest power of 2 smaller than n). Persons with
! np/2<k<=np are parents (but only those with k<=n/2 have actually
! one or two children). Persons np/4<k<=np/2 are grandparents, etc.
! Then, the person with k=1 (the oldest) is removed and the tree is
! reconstructed. This is iterated until all members are picked.
! Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
! *******************************************************************

      implicit none
      integer, parameter  :: dp = kind(1.d0)
      real(dp),parameter  :: tol = 1.e-12_dp  ! tolerance for value comparisons

      integer, intent(in) :: m, n     ! Dimensions of array x
      real(dp),intent(in) :: x(m,n)   ! Array with the values to be ordered
                                      ! will order against the first element (m)
      integer, intent(out):: indx(n)  ! Increasing order of x(1,:)

      integer:: k, nFamily, parent
      real(dp):: ageTol

      ! Construct the heap (family tree)
      indx = (/(k,k=1,n)/)            ! initial array order, to be modified

      ageTol = tol*maxval(x(1,:)) - tol*minval(x(1,:))  ! tolerance for age comparisons

      ! Swap to create the actual parent tree

      nFamily = n                     ! number of persons in the family tree
      do parent = n/2,1,-1            ! sift 'parents' down the tree
        call siftDown(parent)         ! siftDown inherits age and indx arrays
      enddo

      ! Reduce the tree size, retiring its succesive patriarchs (first element)
      do nFamily = n-1,1,-1           ! nFamily is the new size of the tree
        call swap( indx(1), indx(nFamily+1) ) ! swap patriarch and youngest child
        call siftDown(1)              ! now recolocate child in tree
      enddo

      contains

      subroutine siftDown( person )   ! place person in family tree

      implicit none
      integer, intent(in) :: person

      ! Inherited from ordix: age(:), ageTol, indx(:), nFamily
      integer:: child, sw, parent

      ! initialize
      parent = person                    ! assume person is a parent
      child = 2 * parent                 ! first child of parent
      sw = parent                        ! sorted swap child
      do while ( child <= nFamily )      ! iterate the sift-down process
         ! check current child
         if ( x(1,indx(sw)) < x(1,indx(child)) - ageTol ) then
            sw = child
         end if
         ! Check neighbouring child
         if ( child < nFamily ) then
            if ( x(1,indx(sw)) < x(1,indx(child+1)) - ageTol ) then
               sw = child + 1
            end if
         end if
         if ( sw == parent ) then
            exit ! break
         else
            call swap( indx(parent) , indx(sw) )
            ! update for next sift
            parent = sw
            child = sw * 2
         end if
      end do
      
      end subroutine siftDown

      subroutine swap(i,j)            ! exchange integers i and j
      integer:: i,j,k
      k = i
      i = j
      j = k
      end subroutine swap

      END SUBROUTINE ordix


      SUBROUTINE ORDER( X, M, N, INDX )
C *******************************************************************
C Orders array X(M,N) according to array INDEX(N), which may be
C   generated by routine ORDIX.
C Written by J.M.Soler. May'96.
C *************** INPUT *********************************************
C INTEGER M, N     : Dimensions of array X
C INTEGER INDX(N)  : Array which gives the desired order
C *************** INPUT AND OUTPUT **********************************
C REAL*8  X(M,N) : Array(s) to be ordered: Xout(I,J) = Xin(I,INDEX(J))
C *******************************************************************
      IMPLICIT          NONE
      INTEGER           M, N, INDX(N)
      DOUBLE PRECISION  X(M,N)
      X = X(:,INDX)
      END SUBROUTINE ORDER


      SUBROUTINE IORDER( IA, M, N, INDX )
C *******************************************************************
C Orders integer array IA(M,N) according to array INDEX(N), which
C   may be generated by routine ORDIX.
C Written by J.M.Soler. May'96 and Oct'97.
C *************** INPUT *********************************************
C INTEGER M, N     : Dimensions of array IA
C INTEGER INDX(N)  : Array which gives the desired order
C *************** INPUT AND OUTPUT **********************************
C REAL*8  IA(M,N): Array(s) to be ordered: IAout(I,J) = IAin(I,INDEX(J))
C *******************************************************************
      IMPLICIT NONE
      INTEGER  M, N
      INTEGER  IA(M,N), INDX(N)
      IA = IA(:,INDX)
      END SUBROUTINE IORDER

      END MODULE sorting
