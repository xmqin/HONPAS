! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module class_Sparsity
  !! Implements a reference-counted derived type ([bud](|page|/datastructures/1-buds.html))
  !! to hold the [sparsity pattern](|page|/datastructures/2-sparse.html)  of the program's matrices.
  !! Objects of this type (notably [[sparse_matrices:sparse_pattern]] can be passed around
  !! and put in container types.
  !! The legacy indexing arrays can be linked to this bud's data ([[class_Sparsity:attach]])
  !! when necessary.
  use alloc, only: re_alloc, de_alloc
  
  implicit none
  
  public :: newSparsity
  public :: print_type
  public :: attach
  public :: nrows, nrows_g
  public :: ncols, ncols_g
  public :: n_row, n_col
  public :: nnzs, list_ptr, list_col
  public :: equivalent

  character(len=*), parameter :: mod_name= "class_Sparsity"

  ! This is the "meat" of the type
  !
  type Sparsity_
    integer :: refCount = 0
    character(len=36)  :: id = "null_id"
    !----------------------
    character(len=256) :: name = "null_sparsity"
    integer            :: nrows = 0             ! Local number of rows
    integer            :: nrows_g = 0           ! Global number or rows
    integer            :: ncols = 0             ! Local number of columns
    integer            :: ncols_g = 0           ! Global number or columns
    integer            :: nnzs  = 0             ! Local number of non-zeros
    integer, pointer   :: n_col(:)     =>null() ! Nonzero cols of each row
    integer, pointer   :: list_col(:)  =>null() ! Index of nonzero columns
    integer, pointer   :: list_ptr(:)  =>null() ! First element of each row
  end type Sparsity_

  ! This is a wrapper type to be passed around
  type Sparsity
     type(Sparsity_), pointer :: data => null()
  end type Sparsity

  interface nrows
     module procedure nrowsSparsity
  end interface

  interface nrows_g
     module procedure nrows_gSparsity
  end interface

  interface ncols
     module procedure ncolsSparsity
  end interface

  interface ncols_g
     module procedure ncols_gSparsity
  end interface

! ******************
! Specific routines to retrieve direct information
! about the sparsity entries.

  interface attach
     module procedure attachSparsity
  end interface

  interface equivalent
     module procedure equivalentSparsity
  end interface

  interface n_col
     module procedure n_colSparsity
     module procedure n_colSparsityI
  end interface

  interface n_row
     ! This interface ensures that one can retrieve
     ! all information from the sparse class.
     module procedure n_rowSparsityI
  end interface

  interface nnzs
     module procedure nnzsSparsity
  end interface

  interface list_ptr
     module procedure list_ptrSparsity
     module procedure list_ptrSparsityI
  end interface
  interface list_col
     module procedure list_colSparsity
     module procedure list_colSparsityI
  end interface

  interface print_type
     module procedure printSparsity
  end interface

!===========================
#define TYPE_NAME Sparsity
#include "basic_type.inc"
!===========================

  subroutine delete_Data(spdata)
    type(Sparsity_) :: spdata
    call de_alloc( spdata%n_col,   &
         name="n_col " // trim(spdata%name),routine="Sparsity")	
    call de_alloc( spdata%list_ptr,   &
         name="list_ptr " // trim(spdata%name),routine="Sparsity")	
    call de_alloc( spdata%list_col,   &
         name="list_col " // trim(spdata%name),routine="Sparsity") 
  end subroutine delete_Data

!--------------------------------------------------------------------    
! For easy compatability we have added ncols and ncols_g as
! optional arguments.
! In case, not specified, values of nrows_g are used. (block-cyclic)
  subroutine newSparsity(sp,nrows,nrows_g,nnzs,num,listptr,list,name, &
       ncols,ncols_g)

    type(Sparsity), intent(inout) :: sp

    integer, intent(in)           :: nrows, nrows_g, nnzs
    integer, intent(in)           :: num(:), listptr(:)
    integer, intent(in), optional :: list(:)
    character(len=*), intent(in), optional  :: name
    integer, intent(in), optional :: ncols, ncols_g

   ! We release the previous incarnation
   ! This means that we relinquish access to the previous
   ! memory location. It will be deallocated when nobody
   ! else is using it.

    call init(sp)
    
    sp%data%name = trim(name)
    
    call re_alloc(sp%data%n_col,1,nrows, &
         name="n_col " // trim(sp%data%name),routine="Sparsity")	
    call re_alloc(sp%data%list_ptr,1,nrows, &
         name="list_ptr " // trim(sp%data%name),routine="Sparsity")	

    sp%data%nrows = nrows
    sp%data%nrows_g = nrows_g
    if ( present(ncols_g) ) then
       sp%data%ncols_g = ncols_g
    else
       ! We default it to a Block-cyclic distribution, hence
       ! the number of columns is the same as the number of 
       ! global rows (we cannot guess super-cells, that would indeed be amazing)
       sp%data%ncols_g = nrows_g
    end if
    if ( present(ncols) ) then
       sp%data%ncols = ncols
    else
       ! Again, block cyclic has the maximum number of columns
       ! in each block
       sp%data%ncols = sp%data%ncols_g
    end if
    sp%data%nnzs  = nnzs
    sp%data%n_col(1:nrows) = num(1:nrows)
    sp%data%list_ptr(1:nrows) = listptr(1:nrows)

    if (nnzs /= sum(num(1:nrows))) then
       call die("nnzs mismatch in new_sparsity")
    endif

    call re_alloc(sp%data%list_col,1,nnzs, &
         name="list_col " // trim(sp%data%name),routine="Sparsity")
    if ( present(list) ) then
       sp%data%list_col(1:nnzs) = list(1:nnzs)
    else
       sp%data%list_col(1:nnzs) = 0
    end if

    call tag_new_object(sp)

  end subroutine newSparsity

! This works beautifully with gfortran, but intel creeps out
! when destroying the object. :(
!  subroutine pntSparsity(sp,nrows,nrows_g,num,listptr,list,name, &
!       ncols,ncols_g)
!    
!    use alloc, only : alloc_count
!
!    type(Sparsity), intent(inout) :: sp
!
!    integer, intent(in)           :: nrows, nrows_g
!    integer, pointer              :: num(:), listptr(:), list(:)
!    character(len=*), intent(in)  :: name
!    integer, intent(in), optional :: ncols, ncols_g
!
!   ! We release the previous incarnation
!   ! This means that we relinquish access to the previous
!   ! memory location. It will be deallocated when nobody
!   ! else is using it.
!
!    call init(sp)
!    
!    sp%data%name = trim(name)
!
!    if ( size(num) /= nrows ) then
!       call die('pntSparsity: n_col size does not match nrows')
!    end if
!    if ( size(listptr) /= nrows ) then
!       call die('pntSparsity: list_ptr size does not match nrows')
!    end if
!
!    ! Make pointers
!    sp%data%n_col    => num(:)
!    sp%data%list_ptr => listptr(:)
!    sp%data%list_col => list(:)
!    
!    sp%data%nrows   = nrows
!    sp%data%nrows_g = nrows_g
!    if ( present(ncols_g) ) then
!       sp%data%ncols_g = ncols_g
!    else
!       ! We default it to a Block-cyclic distribution, hence
!       ! the number of columns is the same as the number of 
!       ! global rows (we cannot guess super-cells, that would indeed be amazing)
!       sp%data%ncols_g = nrows_g
!    end if
!    if ( present(ncols) ) then
!       sp%data%ncols = ncols
!    else
!       ! Again, block cyclic has the maximum number of columns
!       ! in each block
!       sp%data%ncols = sp%data%ncols_g
!    end if
!    sp%data%nnzs = sum(num(1:nrows))
!    if ( size(list) /= sp%data%nnzs ) then
!       call die('pntSparsity: list size does not match sum(n_col)')
!    end if
!
!    ! Track the sparsity count, this *might* produce wrong
!    ! count of memory!
!    ! TODO CHECK MEMORY COUNT
!    call alloc_count( nrows, 'I', 'n_col '// trim(sp%data%name), "Sparsity" )
!    call alloc_count( nrows, 'I', 'list_ptr '// trim(sp%data%name), "Sparsity" )
!    call alloc_count( sp%data%nnzs, 'I', 'list_col '// trim(sp%data%name), "Sparsity" )
!
!    call tag_new_object(sp)
!
!  end subroutine pntSparsity
  
  pure function nrowsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows
  end function nrowsSparsity

  elemental function n_rowSparsityI(this,col) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: col
    integer                    :: p
    p = count(this%data%list_col == col)
  end function n_rowSparsityI

  pure function nrows_gSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%nrows_g
  end function nrows_gSparsity

  pure function ncolsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%ncols
  end function ncolsSparsity

  pure function ncols_gSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    n = this%data%ncols_g
  end function ncols_gSparsity


  pure function nnzsSparsity(this) result (n)
    type(Sparsity), intent(in) :: this
    integer                    :: n
    if ( initialized(this) ) then
       n = this%data%nnzs
    else
       n = 0
    end if
  end function nnzsSparsity

  function n_colSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%n_col
  end function n_colSparsity
  elemental function n_colSparsityI(this,row) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: row
    integer                    :: p
    p = this%data%n_col(row)
  end function n_colSparsityI

  function list_ptrSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%list_ptr
  end function list_ptrSparsity
  elemental function list_ptrSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%list_ptr(i)
  end function list_ptrSparsityI

  function list_colSparsity(this) result (p)
    type(Sparsity), intent(in) :: this
    integer, pointer           :: p(:)
    p => this%data%list_col
  end function list_colSparsity
  elemental function list_colSparsityI(this,i) result (p)
    type(Sparsity), intent(in) :: this
    integer, intent(in)        :: i
    integer                    :: p
    p = this%data%list_col(i)
  end function list_colSparsityI

  ! Generic routine for retrieval of information in one line...
  subroutine attachSparsity(this,D,n_col,list_col,list_ptr, &
       nrows,nrows_g,ncols,ncols_g, &
       nnzs)
    type(Sparsity), intent(inout) :: this
    logical, intent(in) , optional :: D ! DUMMY, force named arguments
    ! Optional arguments to retrieve all information with named arguments
    integer, pointer    , optional :: n_col(:), list_col(:), list_ptr(:)
    integer, intent(out), optional :: nrows, nrows_g
    integer, intent(out), optional :: ncols, ncols_g
    integer, intent(out), optional :: nnzs
    if ( present(D) ) call die('PROGRAMMING ERROR, named args please')

    ! the arrays
    if ( present(n_col) ) n_col => n_colSparsity(this)
    if ( present(list_col) ) list_col => list_colSparsity(this)
    if ( present(list_ptr) ) list_ptr => list_ptrSparsity(this)
    ! the integers
    if ( present(nrows) ) nrows = nrowsSparsity(this)
    if ( present(nrows_g) ) nrows_g = nrows_gSparsity(this)
    if ( present(ncols) ) ncols = ncolsSparsity(this)
    if ( present(ncols_g) ) ncols_g = ncols_gSparsity(this)
    if ( present(nnzs) ) nnzs = nnzsSparsity(this)

  end subroutine attachSparsity

  function equivalentSparsity(sp1,sp2) result(equivalent)
    type(Sparsity), intent(inout) :: sp1, sp2
    logical :: equivalent
    integer, pointer :: ncol1(:), l_col1(:), ncol2(:), l_col2(:)
    integer :: lno1, lno2, no1, no2

    equivalent = same(sp1,sp2)
    if ( equivalent ) return
    if ( initialized(sp1) ) then
       if (.not. initialized(sp2) ) return
    else if ( initialized(sp2) ) then
       if (.not. initialized(sp1) ) return
    end if
    ! If they are the same object, immediately return
    ! with true
    equivalent = sp1%data%id == sp2%data%id
    if ( equivalent ) return

    ! Both are initialized
    call attach(sp1,nrows=lno1,nrows_g=no1, &
         n_col = ncol1, list_col = l_col1 )
    call attach(sp2,nrows=lno2,nrows_g=no2, &
         n_col = ncol2, list_col = l_col2 )
    equivalent = lno1 == lno2
    if ( .not. equivalent ) return
    equivalent = no1 == no2
    if ( .not. equivalent ) return
    equivalent = all(ncol1 == ncol2)
    if ( .not. equivalent ) return
    equivalent = all(l_col1 == l_col2)

  end function equivalentSparsity

  subroutine printSparsity(sp)
    type(Sparsity), intent(in) :: sp

    if (.not. initialized(sp) ) then
       print "(a)", "Sparsity Not Associated"
       RETURN
    endif

    print "(a,/,t4,2(a,i0),a,f0.4,2(a,i0),a)", &
                "  <sparsity:"//trim(sp%data%name), &
                " nrows_g=", sp%data%nrows_g, &
                " nrows=", sp%data%nrows, &
                " sparsity=",real(sp%data%nnzs)/&
                real(sp%data%nrows_g)/real(sp%data%ncols_g), &
                " nnzs=",sp%data%nnzs,", refcount: ",  &
                refcount(sp), ">"
  end subroutine printSparsity

end module class_Sparsity
