!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! These routines should ALWAYS be kept separable from the SIESTA package.
! No routines must ever be dependent on any SIESTA MODULES.

! A module that supplements routines that generally are missing in 
! The FORTRAN standard.
! Some are also more "exotic" in terms of their use and the FORTRAN standard.

! DEVELOPER REMARK:
!   Many of these functions are not GNU compliant for nested calls.
!   I.e. SORT(UNIQ(arr)) is NOT allowed and will result in errorneous
!   results.

! The following routines are supplied:
!  - VNORM: perform euclidean norm calculations on vectors or matrices.
!           for complex numbers it is equilvalent to VNORM(abs(z))
!  - SORT : allows sorting of an array. It currently only exists
!           for an integer array. However, this is also the most stringent
!           case as there has not to be any margin of error (EPS)
!  - SORT_QUICK : allows sorting of an array (using the quick-sort algorithm)
!  - INDEX_SORT_HEAP : allows returning an index array that sorts the input array (using the heap-sort algorithm)
!  - UNIQ : returns the unique integer values of an array
!  - UNIQC: returns the number of unique integer vaules of an array
!  - SFIND: will return the index of an integer in a sorted array (by logical
!           fast stepping). It has optional parameters to search:
!             * NEAREST=+1: if not found, returns the first index HIGHER
!                           than the sought value (first index if first
!                           value is larger than the searched value)
!             * NEAREST=-1: if not found, returns the first index LOWER
!                           than the sought value (last index if last
!                           value is larger than the searched value)
!             * NEAREST= 0: if not found, returns -1 (NORMAL behaviour)
!  - MODP : Returns "FORTRAN" indexed MOD, i.e. if MOD would return 0,
!           MODP returns the divisor.
!  - MODULOP : Returns "FORTRAN" indexed MODULO, i.e. if MODULO would return 0,
!           MODULOP returns the divisor.
!  - EYE  : a routine to generate identity matrices
!           this is a subroutine as I think a function would produce
!           a large memory foot-print before it actually saves to the array.
!           When called on a complex array, only the REAL diagonal is 1.
!  - TRANSPOSE : transpose matrices
!  - TRACE : transpose matrices
!  - ROTATE : returns a vector rotated theta angles (in radians).
!             For 2D points this rotation is implicitly directioned (in the plane)
!             For 3D points a direction is [1,2,3] is needed to specify the rotation
!             direction.
!  - PROJ : A projection method to project a vector onto a certain subspace
!           Exists only for 3D spaces.
!  - VEC_PROJ : A projection method to project a vector onto another vector

module intrinsic_missing

  implicit none

  private

  integer, parameter :: sp = selected_real_kind(5,10)
  integer, parameter :: dp = selected_real_kind(10,100)

  ! Naming scheme must be "different" than what could be supplied by
  ! the standard, for instance, NORM2 is part of F08STD!
  ! VNORM calculates the euclidiean norm of a vector or matrix
  interface VNORM
     module procedure VNORM_sp_1 ,VNORM_sp_2
     module procedure VNORM_dp_1 ,VNORM_dp_2
     module procedure VNORM_cp_1 ,VNORM_cp_2
     module procedure VNORM_zp_1 ,VNORM_zp_2
  end interface

! Pure functions.. (i.e. callable in interface declarations..)
  public :: VNORM
  public :: FSORT, SORT, SORT_QUICK
  public :: INDEX_SORT_HEAP
  interface INDEX_SORT_HEAP
    module procedure INDEX_SORT_HEAP_I4
  end interface INDEX_SORT_HEAP
  public :: UNIQ, UNIQC
  public :: SFIND
  interface SFIND
     module procedure SORTED_FIND, SORTED_CLOSE
  end interface SFIND
#ifdef INTRINSIC_MISSING_TEST
  public :: SORTED_FIND
  public :: SORTED_REC_FIND
  public :: SORTED_BINARY_FIND
#endif


! Elemental functions (can be called on arrays)
  public :: MODP, MODULOP

! Missing matrix stuff
  public :: EYE
  interface EYE
     module procedure EYE_i_2D
     module procedure EYE_sp_2D, EYE_dp_2D
     module procedure EYE_cp_2D, EYE_zp_2D
     module procedure EYE_i_1D
     module procedure EYE_sp_1D, EYE_dp_1D
     module procedure EYE_cp_1D, EYE_zp_1D
  end interface

  public :: TRANSPOSE
  interface TRANSPOSE
     module procedure TRANSPOSE_dp_2D
     module procedure TRANSPOSE_zp_2D
     module procedure TRANSPOSE_dp_1D
     module procedure TRANSPOSE_zp_1D
  end interface

  public :: TRACE
  interface TRACE
     module procedure TRACE_dp_2D
     module procedure TRACE_zp_2D
     module procedure TRACE_dp_1D
     module procedure TRACE_zp_1D
  end interface

  ! Calculate the angle between two vectors
  public :: ANGLE
  interface ANGLE
     module procedure ANGLE_sp
     module procedure ANGLE_dp
  end interface

! Rotation of points
  public :: ROTATE
  interface ROTATE
     module procedure ROTATE_2D
     module procedure ROTATE_3D
  end interface
  
  ! Projection of N-D vector to N-D space
  public :: SPC_PROJ
  interface SPC_PROJ
     module procedure SPC_PROJ_sp
     module procedure SPC_PROJ_dp
  end interface

  ! Index of projection of N-D vector to N-D space
  public :: IDX_SPC_PROJ
  interface IDX_SPC_PROJ
     module procedure IDX_SPC_PROJ_sp
     module procedure IDX_SPC_PROJ_dp
  end interface

  ! Projection of a N-D vector onto N-D vector
  public :: VEC_PROJ
  interface VEC_PROJ
     module procedure VEC_PROJ_sp
     module procedure VEC_PROJ_dp
  end interface

  ! Projection of a N-D vector onto N-D vector
  ! by keeping the scaling of the onto projected vector
  public :: VEC_PROJ_SCA
  interface VEC_PROJ_SCA
     module procedure VEC_PROJ_SCA_sp
     module procedure VEC_PROJ_SCA_dp
  end interface

contains


! ROTATE a point around origo with some angle \theta
! This is a purely plane rotation
  pure subroutine ROTATE_2D(v,theta)
    real(dp), intent(inout) :: v(2)
    real(dp), intent(in) :: theta
    real(dp) :: rmT(2,2), vv(2)
    rmT(1,1) = cos(theta)
    rmT(1,2) = sin(theta)
    rmT(2,1) = -rmT(1,2)
    rmT(2,2) = rmT(1,1)

    vv(1) = sum(rmT(:,1) * v)
    vv(2) = sum(rmT(:,2) * v)
    v = vv

  end subroutine ROTATE_2D

  pure subroutine ROTATE_3D(v,theta,dir)
    real(dp), intent(inout) :: v(3)
    real(dp), intent(in) :: theta
    integer, intent(in) :: dir
    real(dp) :: rmT(3,3), vv(3)
    rmT(:,:) = 0._dp
    if ( dir == 3 ) then
       rmT(1,1) = cos(theta)
       rmT(1,2) = sin(theta)
       rmT(2,1) = -rmT(1,2)
       rmT(2,2) = rmT(1,1)
       rmT(3,3) = 1._dp
    else if ( dir == 2 ) then
       rmT(1,1) = cos(theta)
       rmT(1,3) = sin(theta)
       rmT(2,2) = 1._dp
       rmT(3,1) = -rmT(1,3)
       rmT(3,3) = rmT(1,1)
    else if ( dir == 1 ) then
       rmT(1,1) = 1._dp
       rmT(2,2) = cos(theta)
       rmT(2,3) = sin(theta)
       rmT(3,2) = -rmT(2,3)
       rmT(3,3) = rmT(2,2)
    end if

    vv(1) = sum(rmT(:,1) * v)
    vv(2) = sum(rmT(:,2) * v)
    vv(3) = sum(rmT(:,3) * v)
    v = vv

  end subroutine ROTATE_3D


! A MOD function which behaves differently on the edge.
! It will NEVER return 0, but instead return the divisor.
! This makes it useful in FORTRAN do-loops which are not zero-based
! but 1-based.
! Thus we have the following scheme:
!  MODP(x',x) == x         for x' == x
!  MODP(x',x) == MOD(x',x) for x' /= x
  elemental function MODP(a,p)
    integer, intent(in) :: a,p
    integer :: MODP
    MODP = MOD(a-1,p) + 1
  end function MODP

! A MODULO function which behaves differently on the edge.
! It will NEVER return 0, but instead return the divisor.
! This makes it useful in FORTRAN do-loops which are not zero-based
! but 1-based.
! Thus we have the following scheme:
!  MODULOP(x',x) == x         for x' == x
!  MODULOP(x',x) == MODULOP(x',x) for x' /= x
  elemental function MODULOP(a,p)
    integer, intent(in) :: a,p
    integer :: MODULOP
    MODULOP = MODULO(a-1,p) + 1
  end function MODULOP

! Function to return the unique COUNT of an integer array.
! Thus will return how many DIFFERENT entries there exists.
  pure function UNIQC(array)
    ! We dont need to see the size of the array
    integer, intent(in) :: array(:)
    integer :: dummy(ubound(array,dim=1))
    integer :: UNIQC
    integer :: i, DA
    
    DA = ubound(array,dim=1)

    ! The first entry is of course not found anywhere else
    if ( DA == 0 ) then
       UNIQC = 0
       return
    end if

    ! Everything else has to be checked...
    UNIQC = 1
    ! Save the first entry
    dummy(1) = array(1)
    do i = 2 , DA
       ! Apparently, inserting a do-loop is worse than count?
       ! This makes me worried.... ???
!       do ii = 1 , UNIQC
!          if ( dummy(ii) == array(i) ) cycle find_num
!       end do
!       UNIQC = UNIQC + 1
!       dummy(UNIQC) = array(i)

       if ( count(dummy(1:UNIQC) == array(i)) == 0 ) then
          UNIQC = UNIQC + 1
          dummy(UNIQC) = array(i)
       end if
    end do

  end function UNIQC

! Function to return the unique COUNT of an integer array (THIS ARRAY HAS
! TO BE SORTED).
! Thus will return how many DIFFERENT entries there exists.
  pure function UNIQC_SORTED(array) result(UNIQC)
    ! We dont need to see the size of the array
    integer, intent(in) :: array(:)
    integer :: UNIQC
    integer :: i, DA
    
    DA = ubound(array,dim=1)

    ! The first entry is of course not found anywhere else
    if ( DA == 0 ) then
       UNIQC = 0
       return
    end if

    ! Everything else has to be checked...
    UNIQC = 1
    ! Save the first entry
    do i = 2 , DA
       if ( array(i-1) /= array(i) ) then
          UNIQC = UNIQC + 1
       end if
    end do

  end function UNIQC_SORTED

  pure function UNIQ(array)
    integer, intent(in) :: array(:)
    integer :: UNIQ(UNIQC(array)) ! we cannot know the size ( it will be up to the user to 
    ! Have the correct size available
    integer :: i,j,DA

    DA = ubound(array,dim=1)

    ! Everything else has to be checked...
    UNIQ(1) = array(1)
    j = 1
    do i = 2 , DA
       ! Apparently, inserting a do-loop is worse than count?
       ! This makes me worried.... ???
!       do ii = 1 , j
!          if ( UNIQ(ii) == array(i) ) cycle find_num
!       end do
!       j = j + 1
!       UNIQ(j) = array(i)

       if ( count(UNIQ(1:j) == array(i)) == 0 ) then
          j = j + 1
          UNIQ(j) = array(i)
       end if

    end do

  end function UNIQ

! Returns an integer array sorted
! It will contain dublicates if encountered. And hence sorting will amount to
! arrays (/1,2,2,3,4,4,5/) for example.
! This sorting routine has been optimized for consecutive segments
! in the original array, hence, sorting on arrays with:
!   [1,2,3,4,25,26,27,28,15,16,17] are VERY fast!
  pure function FSORT(array) result(SO)
    integer, intent(in) :: array(:)
    integer :: SO(ubound(array,dim=1))
    integer :: i,j,DA,h,FM

    DA = ubound(array,dim=1)
    
    ! The first entry is of course not found anywhere else
    if ( DA == 0 ) then
       SO = array
       return
    else if ( DA == 1 ) then
       SO = array
       return
    end if
    
    ! Everything else has to be checked...
    ! This sorting routine is very efficient, and it works
    ! Consider to change this function to a quick-sort algorithm...
    SO(1) = array(1)
    i = 2
    sort_loop: do while ( i <= DA )
       if ( array(i) <= SO(1) ) then
          !print '(a,1000(tr1,i0))','F1',SO(1:i-1),array(i)
          !  put in front of the sorted array
          call insert_front(DA,SO,array,i,FM)
          i = i + FM
          !print '(a,1000(tr1,i0))','F2',SO(1:i-1)
       else if ( SO(i-1) <= array(i) ) then
          !print '(a,1000(tr1,i0))','B1',SO(1:i-1),array(i)
          call insert_back(DA,SO,array,i,FM)
          i = i + FM
          !print '(a,1000(tr1,i0))','B2',SO(1:i-1)
       else if ( i-1 < 35 ) then
          !print '(a,1000(tr1,i0))','M1',SO(1:i-1),array(i)
          ! We assume that for array segments below 35 elements
          ! it will be faster to do array search
          expSearch: do j = 2 , i - 1
             ! It will always be better to search for the end 
             ! of the overlapping region. i.e. less elements to move
             if ( array(i) <= SO(j) ) then
                h = j
                call insert_mid(DA,SO,array,h,i,FM)
                i = i + FM
                exit expSearch ! exit the loop
             end if
          end do expSearch
          !print '(a,1000(tr1,i0))','M2',SO(1:i-1)

       else
          !print '(a,1000(tr1,i0))','S1',SO(1:i-1),array(i)

          ! search using SFIND,
          ! We are taking advantage that both SO(1) and SO(i-1)
          ! has been checked, hence the -1
          j = SFIND(SO(2:i-1),array(i),NEAREST=+1) + 1

          ! Insert directly, we have found what we searched for
          call insert_mid(DA,SO,array,j,i,FM)
          i = i + FM

          !print '(a,1000(tr1,i0))','S2',SO(1:i-1)
       end if
    end do sort_loop

    !do i = 1 , DA -1 
    !   if ( so(i) > so(i+1) )then
    !      print *,i,DA,so(i),so(i+1)
    !      print *,sum(so),sum(array)
    !      call die('wrong sort')
    !   end if
    !end do

    !if ( sum(so) /= sum(array) ) then
    !   print *,sum(so),sum(array)
    !   call die('wrong sort')
    !end if

  contains

    pure subroutine insert_mid(DA,SO,array,sF,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      ! The place where we found a SO(sF) <= array(sA)
      integer, intent(in out) :: sF
      integer, intent(in) :: sA ! The current reached iteration in array
      integer, intent(out) :: P ! the number of inserted values

      ! The last insertion point
      integer :: lA
      
      ! First we will skip to the last non-SAME value
      ! I.e. where we can do the insert in SO
      do while ( sF < sA - 1 .and. SO(sF) == array(sA) )
         sF = sF + 1
      end do
      
      ! Now SO(sF) < array(sA)
      if ( sF >= sA )  then
         call insert_back(DA,SO,array,sA,P)
         return
      end if

!******* OLD
      ! We wish to find the last element of array                               
      ! we can move into the sort'ed array                                      
!      lA = sA + 1                                                               
      ! We know we are in the middle of the sort array                          
      ! hence, we can exploit the next element in the sorted                    
      ! array                                                                   
!      do while ( SO(sF-1) <= array(lA-1) .and. &                                
!           array(lA) < SO(sF) .and. &                                           
!           array(lA-1) <= array(lA) .and. &                                     
!           lA <= DA  ) ! must be in the range [sF;sA-1]                         
!         lA = lA + 1                                                            
!      end do                                                                    
!      do while ( SO(sF) == array(lA) .and. lA <= DA  )                          
         ! must be in the range [sF;sA-1] 
!         lA = lA + 1                                                            
!      end do                       


      lA = sA + 1
      ! We know we are in the middle of the sort array
      ! hence, we can exploit the next element in the sorted
      ! array
      do while ( lA <= DA )
         ! If the previous SO value is larger than the insertion
         if ( SO(sF-1)    >  array(lA - 1) ) exit
         ! If the array is not consecutive
         if ( array(lA-1) >  array(lA) )     exit
         ! If the insertion point array is not consecutive
         if ( array(lA)   >  SO(sF) )        exit
         ! We need to ensure an overcount of 1
         lA = lA + 1
      end do
      !if ( lA <= DA ) then
      !   do while ( SO(sF) == array(lA) .and. lA <= DA  )
           ! must be in the range [sF;sA-1]
      !      lA = lA + 1
      !end do

      ! The number of elements we are pushing in
      P = lA - sA
      ! We have "overcounted"
      !lA = lA - 1

      !print '(6(tr1,a,tr1,i4))', &
      !     'SO |-1| =',SO(sF-1), &
      !     'SO |0| =',SO(sF), &
      !     'SO |+1| =',SO(sF+1)
      !print '(6(tr1,a,tr1,i4))','P===',P,&
      !     'LH1',sA+P-1-(sF+P),&
      !     'RH1',sA-1-sF,&
      !     'LH2',sF+P-1-sF,&
      !     'RH2',sA+P-1-sA
      !print '(6(tr1,i4))',array(sA:sA+5)

      ! Copy the mid to the front of the SO
      SO(sF+P:sA+P-1) = SO(sF:sA-1)
      SO(sF:sF+P-1)   = array(sA:sA+P-1)
      
    end subroutine insert_mid

    pure subroutine insert_front(DA,SO,array,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      integer, intent(in)     :: sA
      integer, intent(out)    :: P ! The number of Pasted values

      ! The last insertion point
      integer :: lA, i

      i = sA + 1
      do lA = i , DA
         ! if the previous element is larger than the current element
         if ( array(lA-1) >= array(lA) ) exit
         ! if the checked point is larger than the insertion point
         if ( array(lA) >= SO(1) ) exit
      end do
      i = lA
      do lA = i , DA
         if ( array(lA) /= SO(1) ) exit
      end do

      
!      lA = sA + 1
      ! Take all the values which are smaller than SO(1)
!      do while ( array(lA-1) < array(lA) .and. &
!           array(lA) < SO(1) .and. &
!           lA <= DA ) 
!         lA = lA + 1
!      end do
      ! Take all the values which are EQUAL to SO(1)
!      do while ( array(lA) == SO(1) .and. lA <= DA ) 
!         lA = lA + 1
!      end do

      ! Number of points found in the sort routine
      P = lA - sA
      ! Correct the overstepping (remark the above line counts correctly!)
      !lA = lA - 1
      !print '(6(tr1,a,tr1,i4))',&
      !     'LH1',P+sA-1-(1+P),&
      !     'RH1',sA-1-1,&
      !     'LH2',P-1,&
      !     'RH2',lA-sA

      ! Copy over the values
      SO(P+1:P+sA-1) = SO(1:sA-1)
      SO(1:P)        = array(sA:lA-1)
      
    end subroutine insert_front

    pure subroutine insert_back(DA,SO,array,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      integer, intent(in)     :: sA
      integer, intent(out)    :: P ! 

      ! The last insertion point
      integer :: lA,i
      
      lA = sA + 1
      ! Step until SO(sA-1) /= array(lA)
      do lA = sA + 1 , DA
         if ( SO(sA-1) /= array(lA-1) ) exit
      end do
! OLD
!      do while ( SO(sA-1) == array(lA-1) .and. lA <= DA )
!         lA = lA + 1
!      end do
      
      ! Step until the last element of SO, SO(sA-1), is not
      ! smaller than the array value
      i = lA
      do lA = i , DA
         if ( SO(sA-1) >= array(lA-1) ) exit
         if ( array(lA-1) >= array(lA) ) exit
      end do
!      do while ( SO(sA-1) < array(lA-1) .and. &
!           array(lA-1) < array(lA) .and. & ! asserts the elements we choose are consecutive
!           lA <= DA )
!         lA = lA + 1
!      end do
      
      ! The number of elements we are pushing in
      P = lA - sA
      !print '(6(tr1,a,tr1,i4))',&
      !     'L/R H',sA+P-1-sA      
      SO(sA:sA+P-1) = array(sA:sA+P-1)

    end subroutine insert_back

  end function FSORT

  ! Returns an integer array sorted
! It will contain dublicates if encountered. And hence sorting will amount to
! arrays (/1,2,2,3,4,4,5/) for example.
! This sorting routine has been optimized for consecutive segments
! in the original array, hence, sorting on arrays with:
!   [1,2,3,4,25,26,27,28,15,16,17] are VERY fast!
  pure subroutine SORT(DA, array, so)
    integer, intent(in) :: DA
    integer, intent(in) :: array(DA)
    integer, intent(out) :: SO(DA)
    integer :: i, j, h, FM

    ! No sorting
    if ( DA == 0 ) return

    SO(1) = array(1)
    if ( DA == 1 ) return
    
    ! Everything else has to be checked...
    ! This sorting routine is very efficient, and it works
    ! Consider to change this function to a quick-sort algorithm...
    i = 2
    sort_loop: do while ( i <= DA )
       if ( array(i) <= SO(1) ) then
          !print '(a,1000(tr1,i0))','F1',SO(1:i-1),array(i)
          !  put in front of the sorted array
          call insert_front(DA,SO,array,i,FM)
          i = i + FM
          !print '(a,1000(tr1,i0))','F2',SO(1:i-1)
       else if ( SO(i-1) <= array(i) ) then
          !print '(a,1000(tr1,i0))','B1',SO(1:i-1),array(i)
          call insert_back(DA,SO,array,i,FM)
          i = i + FM
          !print '(a,1000(tr1,i0))','B2',SO(1:i-1)
       else if ( i-1 < 35 ) then
          !print '(a,1000(tr1,i0))','M1',SO(1:i-1),array(i)
          ! We assume that for array segments below 35 elements
          ! it will be faster to do array search
          expSearch: do j = 2 , i - 1
             ! It will always be better to search for the end 
             ! of the overlapping region. i.e. less elements to move
             if ( array(i) <= SO(j) ) then
                h = j
                call insert_mid(DA,SO,array,h,i,FM)
                i = i + FM
                exit expSearch ! exit the loop
             end if
          end do expSearch
          !print '(a,1000(tr1,i0))','M2',SO(1:i-1)

       else
          !print '(a,1000(tr1,i0))','S1',SO(1:i-1),array(i)

          ! search using SFIND,
          ! We are taking advantage that both SO(1) and SO(i-1)
          ! has been checked, hence the -1
          j = SFIND(SO(2:i-1),array(i),NEAREST=+1) + 1

          ! Insert directly, we have found what we searched for
          call insert_mid(DA,SO,array,j,i,FM)
          i = i + FM

          !print '(a,1000(tr1,i0))','S2',SO(1:i-1)
       end if
    end do sort_loop

    !do i = 1 , DA -1 
    !   if ( so(i) > so(i+1) )then
    !      print *,i,DA,so(i),so(i+1)
    !      print *,sum(so),sum(array)
    !      call die('wrong sort')
    !   end if
    !end do

    !if ( sum(so) /= sum(array) ) then
    !   print *,sum(so),sum(array)
    !   call die('wrong sort')
    !end if

  contains

    pure subroutine insert_mid(DA,SO,array,sF,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      ! The place where we found a SO(sF) <= array(sA)
      integer, intent(in out) :: sF
      integer, intent(in) :: sA ! The current reached iteration in array
      integer, intent(out) :: P ! the number of inserted values

      ! The last insertion point
      integer :: lA, i
      
      ! First we will skip to the last non-SAME value
      ! I.e. where we can do the insert in SO
      do while ( sF < sA - 1 .and. SO(sF) == array(sA) )
         sF = sF + 1
      end do
      
      ! Now SO(sF) < array(sA)
      if ( sF >= sA )  then
         call insert_back(DA,SO,array,sA,P)
         return
      end if

!******* OLD
      ! We wish to find the last element of array                               
      ! we can move into the sort'ed array                                      
!      lA = sA + 1                                                               
      ! We know we are in the middle of the sort array                          
      ! hence, we can exploit the next element in the sorted                    
      ! array                                                                   
!      do while ( SO(sF-1) <= array(lA-1) .and. &                                
!           array(lA) < SO(sF) .and. &                                           
!           array(lA-1) <= array(lA) .and. &                                     
!           lA <= DA  ) ! must be in the range [sF;sA-1]                         
!         lA = lA + 1                                                            
!      end do                                                                    
!      do while ( SO(sF) == array(lA) .and. lA <= DA  )                          
         ! must be in the range [sF;sA-1] 
!         lA = lA + 1                                                            
!      end do                       


      lA = sA + 1
      ! We know we are in the middle of the sort array
      ! hence, we can exploit the next element in the sorted
      ! array
      do while ( lA <= DA )
         ! If the previous SO value is larger than the insertion
         if ( SO(sF-1)    >  array(lA - 1) ) exit
         ! If the array is not consecutive
         if ( array(lA-1) >  array(lA) )     exit
         ! If the insertion point array is not consecutive
         if ( array(lA)   >  SO(sF) )        exit
         ! We need to ensure an overcount of 1
         lA = lA + 1
      end do
      !if ( lA <= DA ) then
      !   do while ( SO(sF) == array(lA) .and. lA <= DA  )
           ! must be in the range [sF;sA-1]
      !      lA = lA + 1
      !end do

      ! The number of elements we are pushing in
      P = lA - sA
      ! We have "overcounted"
      !lA = lA - 1

      !print '(6(tr1,a,tr1,i4))', &
      !     'SO |-1| =',SO(sF-1), &
      !     'SO |0| =',SO(sF), &
      !     'SO |+1| =',SO(sF+1)
      !print '(6(tr1,a,tr1,i4))','P===',P,&
      !     'LH1',sA+P-1-(sF+P),&
      !     'RH1',sA-1-sF,&
      !     'LH2',sF+P-1-sF,&
      !     'RH2',sA+P-1-sA
      !print '(6(tr1,i4))',array(sA:sA+5)

      ! Copy the mid to the front of the SO
      if ( sF + P <= sA - 1 ) then
         SO(sF+P:sA+P-1) = SO(sF:sA-1)
      else
         do i = sF , sA - 1
            SO(i+P) = SO(i)
         end do
      end if
      !SO(sF:sF+P-1)   = array(sA:sA+P-1)
      do i = 0 , P - 1
         SO(sF+i) = array(sA+i)
      end do
      
    end subroutine insert_mid

    pure subroutine insert_front(DA,SO,array,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      integer, intent(in)     :: sA
      integer, intent(out)    :: P ! The number of Pasted values

      ! The last insertion point
      integer :: lA, i

      i = sA + 1
      do lA = i , DA
         ! if the previous element is larger than the current element
         if ( array(lA-1) >= array(lA) ) exit
         ! if the checked point is larger than the insertion point
         if ( array(lA) >= SO(1) ) exit
      end do
      i = lA
      do lA = i , DA
         if ( array(lA) /= SO(1) ) exit
      end do

      
!      lA = sA + 1
      ! Take all the values which are smaller than SO(1)
!      do while ( array(lA-1) < array(lA) .and. &
!           array(lA) < SO(1) .and. &
!           lA <= DA ) 
!         lA = lA + 1
!      end do
      ! Take all the values which are EQUAL to SO(1)
!      do while ( array(lA) == SO(1) .and. lA <= DA ) 
!         lA = lA + 1
!      end do

      ! Number of points found in the sort routine
      P = lA - sA
      ! Correct the overstepping (remark the above line counts correctly!)
      !lA = lA - 1
      !print '(6(tr1,a,tr1,i4))',&
      !     'LH1',P+sA-1-(1+P),&
      !     'RH1',sA-1-1,&
      !     'LH2',P-1,&
      !     'RH2',lA-sA

      ! Copy over the values
      if ( P + 1 <= sA - 1 ) then
         SO(P+1:P+sA-1) = SO(1:sA-1)
      else
         do i = 1 , sA - 1
            SO(P+i) = SO(i)
         end do
      end if
      !SO(1:P)        = array(sA:lA-1)
      do i = 1 , P
         SO(i) = array(sA+i-1)
      end do
      
    end subroutine insert_front

    pure subroutine insert_back(DA,SO,array,sA,P)
      integer, intent(in)     :: DA
      integer, intent(in out) :: SO(DA)
      integer, intent(in)     :: array(DA)
      integer, intent(in)     :: sA
      integer, intent(out)    :: P ! 

      ! The last insertion point
      integer :: lA,i
      
      lA = sA + 1
      ! Step until SO(sA-1) /= array(lA)
      do lA = sA + 1 , DA
         if ( SO(sA-1) /= array(lA-1) ) exit
      end do
! OLD
!      do while ( SO(sA-1) == array(lA-1) .and. lA <= DA )
!         lA = lA + 1
!      end do
      
      ! Step until the last element of SO, SO(sA-1), is not
      ! smaller than the array value
      i = lA
      do lA = i , DA
         if ( SO(sA-1) >= array(lA-1) ) exit
         if ( array(lA-1) >= array(lA) ) exit
      end do
!      do while ( SO(sA-1) < array(lA-1) .and. &
!           array(lA-1) < array(lA) .and. & ! asserts the elements we choose are consecutive
!           lA <= DA )
!         lA = lA + 1
!      end do
      
      ! The number of elements we are pushing in
      P = lA - sA
      !print '(6(tr1,a,tr1,i4))',&
      !     'L/R H',sA+P-1-sA
      do i = 0 , P - 1
         SO(sA+i) = array(sA+i)
      end do

    end subroutine insert_back

  end subroutine SORT

  subroutine index_sort_heap_i4( N, A, IDX )
! *******************************************************************
! SUBROUTINE index_sort_heap_i4( N, A, IDX )
! Makes an index table of increasing array A, with size N
! using the heapsort algorithm.
! Written by J.M.Soler, May.2015
! Re-written by N.R.Papior, Mar.2018, taken from ORDIX (sorting.f)
! *************** INPUT *********************************************
! INTEGER N     : Dimensions of array A
! integer A(N)  : Array with the values to be ordered
! *************** OUTPUT ********************************************
! INTEGER IDX(N)  : Array which gives the increasing order of A(I):
!                    A(IDX(I)) .LE. A(IDX(I+1))
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
    
    integer, intent(in) :: N      ! Dimensions of array x
    integer, intent(in) :: A(n)   ! Array with the values to be ordered
                                  ! will order against the first element (m)
    integer, intent(inout) :: idx(n) ! Increasing order of A(:)
    
    integer:: k, nFamily, parent
    
    ! Construct the heap (family tree)
    idx = (/(k,k=1,n)/)            ! initial array order, to be modified
    
    
    ! Swap to create the actual parent tree
    nFamily = n                     ! number of persons in the family tree
    do parent = n/2,1,-1            ! sift 'parents' down the tree
      call siftDown(parent)         ! siftDown inherits age and indx arrays
    end do

    ! Reduce the tree size, retiring its succesive patriarchs (first element)
    do nFamily = n-1,1,-1           ! nFamily is the new size of the tree
      call swap( idx(1), idx(nFamily+1) ) ! swap patriarch and youngest child
      call siftDown(1)              ! now recolocate child in tree
    end do
    
  contains
    
    subroutine siftDown( person )   ! place person in family tree
      integer, intent(in) :: person
      
      ! Inherited from ordix: age(:), ageTol, idx(:), nFamily
      integer:: child, sw, parent
      
      ! initialize
      parent = person                    ! assume person is a parent
      child = 2 * parent                 ! first child of parent
      sw = parent                        ! sorted swap child
      do while ( child <= nFamily )      ! iterate the sift-down process
        ! check current child
        if ( A(idx(sw)) < A(idx(child)) ) then
          sw = child
        end if
        ! Check neighbouring child
        if ( child < nFamily ) then
          if ( A(idx(sw)) < A(idx(child+1)) ) then
            sw = child + 1
          end if
        end if
        if ( sw == parent ) then
          exit ! break
        else
          call swap( idx(parent) , idx(sw) )
          ! update for next sift
          parent = sw
          child = sw * 2
        end if
      end do
      
    end subroutine siftDown
    
    subroutine swap(i,j) ! exchange integers i and j
      integer:: i,j,k
      k = i
      i = j
      j = k
    end subroutine swap
    
  end subroutine index_sort_heap_i4

  
  recursive pure subroutine sort_quick(n, array)
    integer, intent(in) :: n
    integer, intent(inout) :: array(n)

    integer :: div

    if ( n <= 1 ) return

    ! Retrieve the partition ID
    call partition(n, array, div)
    call sort_quick(div-1, array(1))
    call sort_quick(n-div+1, array(div))

  contains

    ! Partition an array by swapping
    pure subroutine partition(n, array, mark)
      integer, intent(in) :: n
      integer, intent(inout) :: array(n)
      integer, intent(out) :: mark
      
      integer :: i, j, x, t

      ! Simple check for n == 2
      ! Note that n == 1 will NEVER happen because
      ! of the sort_quick check
      if ( n == 2 ) then
         if ( array(2) < array(1) ) then
            t = array(2)
            array(2) = array(1)
            array(1) = t
         end if
         mark = 2
         return
      end if

      i = n / 2
      ! Find the median of the searched elements
      if ( array(1) < array(n) ) then
         if ( array(n) < array(i) ) then
            x = array(n)
         else if ( array(1) < array(i) ) then
            x = array(i)
         else
            x = array(1)
         end if
      else !if ( array(n) < array(1) )
         if ( array(1) < array(i) ) then
            x = array(1)
         else if ( array(n) < array(i) ) then
            x = array(i)
         else
            x = array(n)
         end if
      end if

      ! Starting counters..
      i = 0
      j = n + 1

      do
         ! find point from below which is
         ! above median
         j = j - 1
         do while ( j > 0 )
            if ( array(j) <= x ) exit
            j = j - 1
         end do

         i = i + 1
         do while ( i < n )
            if ( array(i) >= x ) exit
            i = i + 1
         end do

         if ( i < j ) then
            ! exchange array(i) and array(j)
            t = array(i)
            array(i) = array(j)
            array(j) = t
         else if ( i == j ) then
            mark = i + 1
            return
         else
            mark = i
            return
         end if
      end do

    end subroutine partition

  end subroutine sort_quick

! ***************** old CORRECT version ******************
! Returns an integer array sorted
! It will contain dublicates if encountered. And hence sorting will amount to
! arrays (/1,2,2,3,4,4,5/) for example.
!   function SORT(array)
!    integer, intent(in) :: array(:)
!    integer :: SORT(ubound(array,dim=1))
!    integer :: i,j,sa
!
!    sa = size(array)
!
!    ! The first entry is of course not found anywhere else
!    if ( sa == 0 ) then
!       SORT = array
!       return
!    else if ( sa == 1 ) then
!       SORT = array
!       return
!    end if
!    SORT(1) = array(1)
!    sort_loop: do i = 2 , sa
!       do j = 1, i-1
! It will always be better to search for the end 
! of the overlapping region. i.e. less elements to move
!          if ( array(i) < SORT(j) ) then
!             SORT(j+1:i) = SORT(j:i-1)
!             SORT(j) = array(i)
!             cycle sort_loop
!          end if
!       end do
!       SORT(i) = array(i)
!    end do sort_loop
!  end function SORT


! Functions for returning the Euclediean norm of a vector
! This takes two function interfaces:
!  VNORM(vec(:,:)) which will compute the length of each column and return 
!  a vector having length ubound(vec(:,:),dim=2), i.e. 
!  vector length of fastests index
  pure function VNORM_sp_1(vec) result(norm)
    real(sp), intent(in) :: vec(:)
    real(sp) :: norm
    integer :: i
    norm = vec(1) * vec(1)
    do i = 2 , ubound(vec,dim=1)
       norm = norm + vec(i) * vec(i)
    end do
    norm = sqrt(norm)
  end function VNORM_sp_1
  pure function VNORM_sp_2(vec) result(norm)
    real(sp), intent(in) :: vec(:,:)
    real(sp) :: norm(ubound(vec,dim=2))
    integer :: i
    do i = 1 , ubound(vec,dim=2)
       norm(i) = VNORM(vec(:,i))
    end do
  end function VNORM_sp_2
  pure function VNORM_dp_1(vec) result(norm)
    real(dp), intent(in) :: vec(:)
    real(dp) :: norm
    integer :: i
    norm = vec(1) * vec(1)
    do i = 2 , ubound(vec,dim=1)
       norm = norm + vec(i) * vec(i)
    end do
    norm = dsqrt(norm)
  end function VNORM_dp_1
  pure function VNORM_dp_2(vec) result(norm)
    real(dp), intent(in) :: vec(:,:)
    real(dp) :: norm(ubound(vec,dim=2))
    integer :: i
    do i = 1 , ubound(vec,dim=2)
       norm(i) = VNORM(vec(:,i))
    end do
  end function VNORM_dp_2
  pure function VNORM_cp_1(vec) result(norm)
    complex(sp), intent(in) :: vec(:)
    real(sp) :: norm
    integer :: i
    norm = conjg(vec(1)) * vec(1)
    do i = 2 , ubound(vec,dim=1)
       norm = norm + conjg(vec(i)) * vec(i)
    end do
    norm = sqrt(norm)
  end function VNORM_cp_1
  pure function VNORM_cp_2(vec) result(norm)
    complex(sp), intent(in) :: vec(:,:)
    real(sp) :: norm(ubound(vec,dim=2))
    integer :: i
    do i = 1 , ubound(vec,dim=2)
       norm(i) = VNORM(vec(:,i))
    end do
  end function VNORM_cp_2
  pure function VNORM_zp_1(vec) result(norm)
    complex(dp), intent(in) :: vec(:)
    real(dp) :: norm
    integer :: i
    norm = dconjg(vec(1)) * vec(1)
    do i = 2 , ubound(vec,dim=1)
       norm = norm + dconjg(vec(i)) * vec(i)
    end do
    norm = dsqrt(norm)
  end function VNORM_zp_1
  pure function VNORM_zp_2(vec) result(norm)
    complex(dp), intent(in) :: vec(:,:)
    real(dp) :: norm(ubound(vec,dim=2))
    integer :: i
    do i = 1 , ubound(vec,dim=2)
       norm(i) = VNORM(vec(:,i))
    end do
  end function VNORM_zp_2

  pure function ANGLE_sp(v1,v2) result(ang)
    real(sp), intent(in) :: v1(:), v2(:)
    real(sp) :: ang
    ang = acos( sum(v1 * v2) / (VNORM(v1)+VNORM(v2)) )
  end function ANGLE_sp
  pure function ANGLE_dp(v1,v2) result(ang)
    real(dp), intent(in) :: v1(:), v2(:)
    real(dp) :: ang
    ang = acos( sum(v1 * v2) / (VNORM(v1)+VNORM(v2)) )
  end function ANGLE_dp


! We add an algorithm for search for SORTED entries
! This will relieve a lot of programming in other segments
! of the codes.
! What it does is the following:
!  1. Get a sorted array of some size
!  2. Get a value which should be found in the array
!  3. Returns the index of the value found in the array
!  4. If the value is not found, return 0
!  5. If something went wrong, return -1
  pure function SORTED_FIND(array,val) result(idx)
    integer, intent(in) :: array(:)
    integer, intent(in) :: val
    ! Used internal variables
    integer :: DA,h,FM, idx
    
    ! Retrieve the size of the array we search in
    DA = ubound(array,dim=1)

    ! Initialize to default value
    idx = 0
    if ( DA == 0 ) return

    ! The two easiest cases, i.e. they are not in the array...
    if ( val < array(1) ) then
       return
    else if ( val == array(1) ) then
       idx = 1
       return
    else if ( array(DA) < val ) then
       return
    else if ( val == array(DA) ) then
       idx = DA
       return
    end if
    
    ! If it is size 1 or 2, we can immediately return,
    ! we have already checked the first/last index
    if ( DA <= 2 ) return

    ! An *advanced* search algorithm...

    ! Search the sorted array for an entry
    ! We know it *must* have one
    h = DA / 2
    idx = h ! Start in the middle
    ! The integer correction (due to round of errors when 
    ! calculating the new half...
    FM = MOD(h,2)
    do while ( h > 0 ) ! While we are still searching...

       if ( h >= 3 ) then
          h = h + FM
          ! This will only add 1 every other time 
          FM = MOD(h,2)
       end if

       ! integer division is fast
       h = h / 2

       if ( val < array(idx) ) then
          ! the value we search for is smaller than 
          ! the current checked value, hence we step back
          !print *,'stepping down',i,h
          idx = idx - h
       else if ( array(idx) < val ) then
          ! the value we search for is larger than 
          ! the current checked value, hence we step forward
          !print *,'stepping up',i,h
          idx = idx + h
       else
          !print *,'found',i
          ! We know EXACTLY where we are...
          return
          ! We found it!!!
       end if
    end do
          
    ! We need to ensure a range to search in
    ! The missing integers are *only* necesseary when 
    ! the search pattern is in the same direction.
    ! This can easily be verified...
    h  = 1 + FM

    ! The missing integer count ensures the correct range
    FM = max(idx - h, 1 )
    h  = min(idx + h, DA)

    ! The index will *most* likely be close to 'i'
    ! Hence we start by searching around it
    ! However, we know that val /= array(i)

    do idx = FM,  h 
       if ( val == array(idx) ) return
    end do
    
    ! Default value is *not found*
    idx = 0

  end function SORTED_FIND

  pure function SORTED_BINARY_FIND(array,val) result(idx)
    integer, intent(in) :: array(:)
    integer, intent(in) :: val

    ! Used internal variables
    integer :: idx
    integer :: small, large, d

    ! Retrieve the size of the array we search in
    d = ubound(array,dim=1)

    ! Initialize to default value
    idx = 0
    if ( d == 0 ) return

    ! The two easiest cases, i.e. they are not in the array...
    if ( val < array(1) ) then
       return
    else if ( val == array(1) ) then
       idx = 1
       return
    else if ( array(d) < val ) then
       return
    else if ( val == array(d) ) then
       idx = d
       return
    end if
    
    ! If it is size 1 or 2, we can immediately return,
    ! we have already checked the first/last index
    if ( d <= 2 ) return

    small = 1
    large = d
    do while ( small + 1 < large )
      idx = (small + large) / 2
      if ( array(idx) < val ) then
        small = idx
      else
        large = idx
      end if
    end do
    if ( array(idx) == val ) return
    do idx = small, large
      if ( array(idx) == val ) return
    end do
    idx = 0

  end function SORTED_BINARY_FIND

  pure recursive function SORTED_REC_FIND(array,val) result(idx)
    integer, intent(in) :: array(:)
    integer, intent(in) :: val
    ! Used internal variables
    integer :: DA, idx, tmp
    
    ! Retrieve the size of the array we search in
    DA = ubound(array,dim=1)

    ! Initialize to default value
    idx = 0
    if ( DA == 0 ) return

    ! The two easiest cases, i.e. they are not in the array...
    if ( val < array(1) ) then
       return
    else if ( val == array(1) ) then
       idx = 1
       return
    else if ( array(DA) < val ) then
       return
    else if ( val == array(DA) ) then
       idx = DA
       return
    end if
    
    ! If it is size 1, we can immediately return,
    ! we have already checked the first index
    if ( DA <= 2 ) then
       idx = 0
       return
    end if

    idx = DA / 2
    if ( val <= array(idx) ) then
       idx = SORTED_REC_FIND(array(1:idx),val)
    else
       tmp = SORTED_REC_FIND(array(idx+1:),val)
       if ( tmp == 0 ) then
          idx = 0
       else
          idx = idx + tmp
       end if
    end if

  end function SORTED_REC_FIND

! We add an algorithm for search for SORTED entries
! This will relieve a lot of programming in other segments
! of the codes.
! What it does is the following:
!  1. Get a sorted array of some size
!  2. Get a value which should be found in the array
!  3. Returns the index of the value found in the array
!  4. If the value is not found, return 0
!  5. If something went wrong, return -1
  pure function SORTED_CLOSE(array,val,NEAREST) result(idx)
    integer, intent(in) :: array(:)
    integer, intent(in) :: val
    integer, intent(in) :: NEAREST
    ! Used internal variables
    integer :: DA,h,FM, idx
    
    ! Retrieve the size of the array we search in
    DA = ubound(array,dim=1)

    ! Initialize to default value
    idx = 0
    if ( DA == 0 ) return

    ! The two easiest cases, i.e. they are not in the array...
    if ( val < array(1) ) then
       ! We need only handle if the user requests
       ! *closest above*
       if ( NEAREST > 0 ) idx = 1
       return
    else if ( val == array(1) ) then
       idx = 1
       return
    else if ( array(DA) < val ) then
       ! We need only handle if the user requests
       ! *closest below*
       if ( NEAREST < 0 ) idx = DA
       return
    else if ( val == array(DA) ) then
       idx = DA
       return
    end if
    
    ! If it is size 1, we can immediately return,
    ! we have already checked the first index
    if ( DA == 1 ) return

    ! An *advanced* search algorithm...

    ! Search the sorted array for an entry
    ! We know it *must* have one
    h = DA / 2
    if ( h < 1 ) h = DA ! This ensures a correct handling
    idx = h ! Start in the middle
    ! The integer correction (due to round of errors when 
    ! calculating the new half...
    FM = MOD(h,2)
    do while ( h > 1 ) ! While we are still searching...

       if ( h >= 3 ) then
          h = h + FM
          ! This will only add 1 every other time 
          FM = MOD(h,2)
       end if
       ! integer division is faster. :)
       h = h / 2

       ! We must ensure to never step UNDER zero
       ! This we can only do by using "floor"
       !h = int(h*0.5_dp)

       ! This makes sure that we 'track' the potential
       ! missing integer while performing 'floor'
       ! This missing integer, will only arise
       ! when the number is non-divisable by two
       !FM = FM + MOD(h,2)

       if ( val < array(idx) ) then
          ! the value we search for is smaller than 
          ! the current checked value, hence we step back
          !print *,'stepping down',i,h
          idx = idx - h
       else if ( array(idx) < val ) then
          ! the value we search for is larger than 
          ! the current checked value, hence we step forward
          !print *,'stepping up',i,h
          idx = idx + h
       else
          !print *,'found',i
          ! We know EXACTLY where we are...
          return
          ! We found it!!!
       end if
    end do
          
    ! We need to ensure a range to search in
    ! The missing integers are *only* necesseary when 
    ! the search pattern is in the same direction.
    ! This can easily be verified...
    h  = max(h,1) + FM + 1

    ! The missing integer count ensures the correct range
    FM = max(idx - h, 1 )
    h  = min(idx + h, DA)

    ! The index will *most* likely be close to 'i'
    ! Hence we start by searching around it
    ! However, we know that val /= array(i)

    ! We need to determine the method of indexing
    if ( NEAREST < 0 ) then
       do idx = FM, h - 1
          if ( val == array(idx) ) then
             return
          else if ( val < array(idx+1) ) then
             return
          end if
       end do
       ! This checks for array(h)
       if ( val == array(h) ) then
          idx = h
          return
       end if
    else if ( NEAREST == 0 ) then
       do idx = FM,  h 
          if ( val == array(idx) ) return
       end do
    else if ( NEAREST > 0 ) then
       do idx = FM, h 
          if ( val == array(idx) ) then
             return
          else if ( val < array(idx) .and. idx > 1 ) then
             if ( val > array(idx-1) ) return
          end if
       end do
    else
       ! ERROR
       idx = -1
    end if
    
    ! Default value is *not found*
    idx = 0

  end function SORTED_CLOSE

  
  ! Matrix operations missing
  ! The reason for choosing a subroutine for these
  ! are the direct impact on memory for very large matrices.
  ! This ensures direct writes, instead of temporary eye-arrays...
  subroutine EYE_i_2D(size,array,I)
    integer, intent(in) :: size
    integer, intent(out) :: array(size,size)
    integer, intent(in), optional :: I
    integer :: j, k, lI
    lI = 1
    if ( present(I) ) lI = I
!$OMP parallel do default(shared), private(k,j)
    do k = 1 , size
       do j = 1 , size
          array(j,k) = 0
       end do
       array(k,k) = lI
    end do
!$OMP end parallel do
  end subroutine EYE_i_2D
  subroutine EYE_sp_2D(size,array,I)
    integer, intent(in) :: size
    real(sp), intent(out) :: array(size,size)
    real(sp), intent(in), optional :: I
    real(sp) :: lI
    integer :: j, k
    lI = 1._sp
    if ( present(I) ) lI = I
!$OMP parallel do default(shared), private(k,j)
    do k = 1 , size
       do j = 1 , size
          array(j,k) = 0._sp
       end do
       array(k,k) = lI
    end do
!$OMP end parallel do
  end subroutine EYE_sp_2D
  subroutine EYE_dp_2D(size,array,I)
    integer, intent(in) :: size
    real(dp), intent(out) :: array(size,size)
    real(dp), intent(in), optional :: I
    real(dp) :: lI
    integer :: j, k
    lI = 1._dp
    if ( present(I) ) lI = I
!$OMP parallel do default(shared), private(k,j)
    do k = 1 , size
       do j = 1 , size
          array(j,k) = 0._dp
       end do
       array(k,k) = lI
    end do
!$OMP end parallel do
  end subroutine EYE_dp_2D
  subroutine EYE_cp_2D(size,array,I)
    integer, intent(in) :: size
    complex(sp), intent(out) :: array(size,size)
    complex(sp), intent(in), optional :: I
    complex(sp) :: lI
    integer :: j, k
    lI = cmplx(1._sp,0._sp)
    if ( present(I) ) lI = I
!$OMP parallel do default(shared), private(k,j)
    do k = 1 , size
       do j = 1 , size
          array(j,k) = 0._sp
       end do
       array(k,k) = lI
    end do
!$OMP end parallel do
  end subroutine EYE_cp_2D
  subroutine EYE_zp_2D(size,array,I)
    integer, intent(in) :: size
    complex(dp), intent(out) :: array(size,size)
    complex(dp), intent(in), optional :: I
    complex(dp) :: lI
    integer :: j, k
    lI = dcmplx(1._dp,0._dp)
    if ( present(I) ) lI = I
!$OMP parallel do default(shared), private(k,j)
    do k = 1 , size
       do j = 1 , size
          array(j,k) = dcmplx(0._dp,0._dp)
       end do
       array(k,k) = lI
    end do
!$OMP end parallel do
  end subroutine EYE_zp_2D

  subroutine EYE_i_1D(size,array)
    integer, intent(in) :: size
    integer, intent(out) :: array(size*size)
    call EYE_i_2D(size,array)
  end subroutine EYE_i_1D
  subroutine EYE_sp_1D(size,array)
    integer, intent(in) :: size
    real(sp), intent(out) :: array(size*size)
    call EYE_sp_2D(size,array)
  end subroutine EYE_sp_1D
  subroutine EYE_dp_1D(size,array)
    integer, intent(in) :: size
    real(dp), intent(out) :: array(size*size)
    call EYE_dp_2D(size,array)
  end subroutine EYE_dp_1D
  subroutine EYE_cp_1D(size,array)
    integer, intent(in) :: size
    complex(sp), intent(out) :: array(size*size)
    call EYE_cp_2D(size,array)
  end subroutine EYE_cp_1D
  subroutine EYE_zp_1D(size,array)
    integer, intent(in) :: size
    complex(dp), intent(out) :: array(size*size)
    call EYE_zp_2D(size,array)
  end subroutine EYE_zp_1D

  
  subroutine TRANSPOSE_dp_2D(size,array)
    integer, intent(in) :: size
    real(dp), intent(inout) :: array(size,size)
    integer :: i, j
    real(dp) :: d
!$OMP parallel do default(shared), private(i,j,d)
    do i = 1 , size - 1
       do j = i + 1 , size
          d = array(j,i)
          array(j,i) = array(i,j)
          array(i,j) = d
       end do
    end do
!$OMP end parallel do
  end subroutine TRANSPOSE_dp_2D
  subroutine TRANSPOSE_zp_2D(size,array)
    integer, intent(in) :: size
    complex(dp), intent(inout) :: array(size,size)
    integer :: i, j
    complex(dp) :: z
!$OMP parallel do default(shared), private(i,j,z)
    do i = 1 , size - 1
       do j = i + 1 , size
          z = array(j,i)
          array(j,i) = array(i,j)
          array(i,j) = z
       end do
    end do
!$OMP end parallel do
  end subroutine TRANSPOSE_zp_2D
  
  subroutine TRANSPOSE_dp_1D(size,array)
    integer, intent(in) :: size
    real(dp), intent(inout) :: array(size*size)
    call TRANSPOSE_dp_2D(size,array)
  end subroutine TRANSPOSE_dp_1D
  subroutine TRANSPOSE_zp_1D(size,array)
    integer, intent(in) :: size
    complex(dp), intent(inout) :: array(size*size)
    call TRANSPOSE_zp_2D(size,array)
  end subroutine TRANSPOSE_zp_1D


  function TRACE_dp_2D(size,array) result(T)
    integer, intent(in) :: size
    real(dp), intent(in) :: array(size,size)
    real(dp) :: T
    integer :: i
    T = 0._dp
!$OMP parallel do default(shared), private(i), reduction(+:T)
    do i = 1 , size
       T = T + array(i,i)
    end do
!$OMP end parallel do
  end function TRACE_dp_2D 
  function TRACE_zp_2D(size,array) result(T)
    integer, intent(in) :: size
    complex(dp), intent(in) :: array(size,size)
    complex(dp) :: T
    integer :: i
    T = dcmplx(0._dp,0._dp)
!$OMP parallel do default(shared), private(i), reduction(+:T)
    do i = 1 , size
       T = T + array(i,i)
    end do
!$OMP end parallel do
  end function TRACE_zp_2D
  
  function TRACE_dp_1D(size,array) result(T)
    integer, intent(in) :: size
    real(dp), intent(in) :: array(size*size)
    real(dp) :: T
    T = TRACE_dp_2D(size,array)
  end function TRACE_dp_1D
  function TRACE_zp_1D(size,array) result(T)
    integer, intent(in) :: size
    complex(dp), intent(in) :: array(size*size)
    complex(dp) :: T
    T = TRACE_zp_2D(size,array)
  end function TRACE_zp_1D


  ! Projections from one direction onto ND space
  pure function SPC_PROJ_sp(space,vin) result(vout)
    real(sp), intent(in) :: space(:,:), vin(:)
    real(sp) :: vout(size(vin)), tmp(size(vin))
    integer :: i
    do i = 1 , size(vin)
       tmp(:)  = space(:,i) / VNORM(space(:,i))
       vout(i) = sum( vin(:) * tmp(:) )
    end do
  end function SPC_PROJ_sp
  pure function SPC_PROJ_dp(space,vin) result(vout)
    real(dp), intent(in) :: space(:,:), vin(:)
    real(dp) :: vout(size(vin)), tmp(size(vin))
    integer :: i
    do i = 1 , size(vin)
       tmp(:)  = space(:,i) / VNORM(space(:,i))
       vout(i) = sum( vin(:) * tmp(:) )
    end do
  end function SPC_PROJ_dp

  pure function IDX_SPC_PROJ_sp(space,vin,mag) result(idx)
    real(sp), intent(in) :: space(:,:), vin(:)
    ! If mag is false (default) it takes the largest
    ! POSITIVE projection.
    ! If mag is true it takes the largest MAGNITUDE projection
    logical, intent(in), optional :: mag
    integer :: idx
    real(sp) :: spc(size(vin))
    integer :: i
    idx = 1
    ! Find the largest contributing space-vector
    ! We _only_ compare against the projected
    ! vector onto the normal length of the vector
    if ( present(mag) ) then
       if ( mag ) then
          spc(:) = abs( SPC_PROJ(space,vin) )
       else
          spc(:) = SPC_PROJ(space,vin)
       end if
    else
       spc(:) = SPC_PROJ(space,vin)
    end if
    ! Find the largest contributing one
    ! it does not matter whether it is plus 
    ! or minus, so long as it its magnitude
    ! is larger (one would _never_ have two
    ! vectors of opposite direction (i.e. a 2D surface)
    do i = 2 , size(vin)
       if ( spc(idx) < spc(i) ) then
          idx = i
       end if
    end do
  end function IDX_SPC_PROJ_sp

  pure function IDX_SPC_PROJ_dp(space,vin,mag) result(idx)
    real(dp), intent(in) :: space(:,:), vin(:)
    ! If mag is false (default) it takes the largest
    ! POSITIVE projection.
    ! If mag is true it takes the largest MAGNITUDE projection
    logical, intent(in), optional :: mag
    integer :: idx
    real(dp) :: spc(size(vin))
    integer :: i
    idx = 1
    if ( present(mag) ) then
       if ( mag ) then
          spc(:) = abs( SPC_PROJ(space,vin) )
       else
          spc(:) = SPC_PROJ(space,vin)
       end if
    else
       spc(:) = SPC_PROJ(space,vin)
    end if
    do i = 2 , size(vin)
       if ( spc(idx) < spc(i) ) then
          idx = i
       end if
    end do
  end function IDX_SPC_PROJ_dp


  ! Scalar projection of 'vin=a' onto 'vec=b'.
  !   a . b / |b|
  ! I.e. for these vectors:
  !   vec = [2, 0, 0]
  !   vin = [3, 1, 1]
  ! a = [3, 0, 0]
  pure function VEC_PROJ_SCA_sp(vec,vin) result(a)
    real(sp), intent(in) :: vec(:), vin(:)
    real(sp) :: a
    a = sum(vec*vin) / VNORM(vec)
  end function VEC_PROJ_SCA_sp
  pure function VEC_PROJ_SCA_dp(vec,vin) result(a)
    real(dp), intent(in) :: vec(:), vin(:)
    real(dp) :: a
    a = sum(vec*vin) / VNORM(vec)
  end function VEC_PROJ_SCA_dp

  ! Projection of 'vin=a' onto 'vec=b'
  !   a . b / ( b . b ) * b
  pure function VEC_PROJ_sp(vec,vin) result(vout)
    real(sp), intent(in) :: vec(:), vin(:)
    real(sp) :: vout(size(vec))
    vout = sum(vec * vin) / sum(vec * vec) * vec
  end function VEC_PROJ_sp
  pure function VEC_PROJ_dp(vec,vin) result(vout)
    real(dp), intent(in) :: vec(:), vin(:)
    real(dp) :: vout(size(vec))
    vout = sum(vec * vin) / sum(vec * vec) * vec
  end function VEC_PROJ_dp

end module intrinsic_missing

#ifdef INTRINSIC_MISSING_TEST
program test
  use intrinsic_missing

  integer, parameter :: sp = selected_real_kind(5,10)
  integer, parameter :: dp = selected_real_kind(10,100)
  real(dp) :: cell(3,3), v(3)
  real :: t1 , t2

  integer :: i, j
  integer, parameter :: N = 21355
  integer :: list(N)
  

  cell(:,:) = 0._dp
  cell(1,1) = 10._dp
  cell(1,2) = 1._dp
  cell(2,2) = 5._dp
  cell(1,3) = 1._dp
  cell(2,3) = 5._dp
  cell(3,3) = 10._dp

  v(:) = 0._dp
  v(1) = 2._dp
  print *,1,'==',IDX_SPC_PROJ(cell,v)
  v(2) = 5._dp
  print *,2,'==',IDX_SPC_PROJ(cell,v)
  v(3) = 10._dp
  print *,3,'==',IDX_SPC_PROJ(cell,v)
  v(2) = v(1)
  v(3) = 0._dp
  print *,2,'==',IDX_SPC_PROJ(cell,v)
  v(2) = v(1) / 10._dp
  v(3) = 0._dp
  print *,1,'==',IDX_SPC_PROJ(cell,v)

  do i = 1 , N
     list(i) = i
  end do

  call cpu_time(t1)
  do j = 2 , N
     do i = 1 , j
        if ( SORTED_FIND(list(1:j),i) /= i ) print *,'Error',i
     end do
     do i = 1 , j
        if ( SORTED_FIND(list(1:j),i+j) /= 0 ) print *,'Error',i
     end do
  end do
  call cpu_time(t2)
  print '(a,f10.5,a)','Timing of SORTED_FIND: ',t2-t1,' secs'

  call cpu_time(t1)
  do j = 2 , N
     do i = 1 , j
        if ( SORTED_REC_FIND(list(1:j),i) /= i ) print *,'Error',i
     end do
     do i = 1 , j
        if ( SORTED_REC_FIND(list(1:j),i+j) /= 0 ) print *,'Error',i
     end do
  end do
  call cpu_time(t2)
  print '(a,f10.5,a)','Timing of SORTED_REC_FIND: ',t2-t1,' secs'

  call cpu_time(t1)
  do j = 2 , N
     do i = 1 , j
        if ( SORTED_BINARY_FIND(list(1:j),i) /= i ) print *,'Error',i
     end do
     do i = 1 , j
        if ( SORTED_BINARY_FIND(list(1:j),i+j) /= 0 ) print *,'Error',i
     end do
  end do
  call cpu_time(t2)
  print '(a,f10.5,a)','Timing of SORTED_BINARY_FIND: ',t2-t1,' secs'

end program test
#endif
