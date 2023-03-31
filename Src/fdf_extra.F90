
module fdf_extra

  ! Module which extends the fdf routines
  ! by a little...
  use fdf
  use m_region
  
  implicit none
  
  private
  
  public :: fdf_bnext
  public :: fdf_brange
  public :: fdf_bregions

contains

  function fdf_bnext(bfdf,pline) result(has)
    type(block_fdf), intent(inout) :: bfdf
    type(parsed_line), pointer :: pline
    logical :: has
    do 
      has = fdf_bline(bfdf,pline)
      if ( .not. has ) return
      if ( fdf_bntokens(pline) > 0 ) return
    end do
    
  end function fdf_bnext

  !> Read a range of integers from a parsed line
  !!
  !! It returns the list of integers in a `tRng` object.
  subroutine fdf_brange(pline,r, low, high)
    type(parsed_line), pointer :: pline
    type(tRgn), intent(out) :: r
    ! Practically min/max for the region
    ! in question, however, it allows for 
    ! specifying negative numbers.
    integer, intent(in) :: low, high

    integer :: j, n, i, i1, i2, il, nl, step
    integer, allocatable :: list(:)
    character(len=32) :: g

    ! We allocate a temporary list
    allocate(list(high-low+1))
    
    ! Number of read elements
    n = 0

    ! If lists exists, we use those
    if ( fdf_bnnames(pline) == 1 ) then
      
      ! Get number of lists on current line
      nl = fdf_bnlists(pline)

      if ( nl == 0 ) then

        ! The input line is a simple list of integers
        do j = 1 , fdf_bnintegers(pline)
          i = fdf_bintegers(pline,j)
          n = n + 1 
          list(n) = correct(i,low,high)
        end do

      else

        do il = 1, nl

          ! Read in number of items
          i1 = -1
          call fdf_blists(pline,il,i1,list(n+1:))
          if ( i1 + n > size(list) ) then
            print *,'Parsed line: ',trim(pline%line)
            call die('fdf_brange: number of elements in block list &
                &is too many to fit the maximal range of the &
                &list. Please correct.')
          end if
          if ( i1 == 0 ) then
            print *,'Parsed line: ',trim(pline%line)
            call die('fdf_brange: a block list with zero elements is not &
                &allowed, please correct input.')
          end if

          ! Read in actual list
          call fdf_blists(pline,il,i1,list(n+1:n+i1))

          do i = 1 , i1
            list(n+i) = correct(list(n+i),low,high)
          end do

          ! update n
          n = n + i1

        end do

      end if

    else if ( fdf_bnnames(pline) == 3 ) then

      g = fdf_bnames(pline,2)
      if ( .not. leqi(g,'from') ) then
        print *,'Parsed line: ',trim(pline%line)
        call die('fdf_brange: error in range block: &
            &from <int> to/plus/minus <int> is ill formatted')
      end if

      g = fdf_bnames(pline,3)
      if ( fdf_bnintegers(pline) < 2 ) then
        print *,'Parsed line: ',trim(pline%line)
        call die('fdf_brange: error in range block &
            &from <int> to/plus/minus <int> is ill formatted')
      end if

      ! Read in arguments
      i1 = fdf_bintegers(pline,1)
      i2 = fdf_bintegers(pline,2)
      if ( fdf_bnintegers(pline) > 2 ) then
        step = fdf_bnintegers(pline,3)
        if ( step == 0 ) &
            call die('fdf_extra: stepping MUST be different from 0')
      else
        step = 0
      end if

      ! Check how the arguments are read
      if ( leqi(g,'to') ) then
        ! do nothing....
      else if ( leqi(g,'plus') ) then
        i2 = i1 + i2 - 1
      else if ( leqi(g,'minus') ) then
        ! possibly correct step
        if ( step > 0 ) step = - step
        i2 = i1 - i2 + 1
      else
        print *,'Parsed line: ',trim(pline%line)
        call die('fdf_brange: unrecognized designator of ending range, &
            &[to, plus, minus] accepted.')
      end if

      ! Figure out the step-level
      if ( step == 0 ) then
        if ( i1 <= i2 ) then
          step = 1
        else
          step = -1
        end if
      end if
      
      ! Check input for list creation
      if ( (i1 < i2 .and. step < 0) .or. &
          (i1 > i2 .and. step > 0) ) then
        print *,'Parsed line: ',trim(pline%line)
        print *,i1,i2,step
        call die('fdf_brange: block range is not consecutive')
      end if

      ! Create list
      do i = i1 , i2 , step
        n = n + 1
        ! correct for wrap-arounds
        list(n) = correct(i,low,high)
      end do

    else

      print *,'Parsed line: ',trim(pline%line)
      call die('fdf_brange: error in range block, input not recognized')

    end if

    do j = 1 , n
      if ( list(j) < low .or. high < list(j) ) then
        print *,'Parsed line: ',trim(pline%line)
        call die('fdf_brange: error in range block. Input is beyond range')
      end if
    end do

    call rgn_list(r,n,list)
    deallocate(list)

  contains

    pure function correct(in,low,high) result(out)
      integer, intent(in) :: in, low, high
      integer :: out
      if ( low == 1 .and. high >= low ) then
        ! Special case, we use wrap-arounds
        out = in
        if ( out == 0 ) then
          ! this is not allowed
          return
        end if
        if ( out < 0 ) then
          ! correct for the zero-base
          out = out + 1
          do while ( out < 1 )
            out = out + high
          end do
        end if
        if ( out > high ) then
          do while ( high < out )
            out = out - high
          end do
        end if
      else
        if ( in < low ) then
          out = high + in + 1
        else if ( high < in ) then
          out = in - high + low - 1
        else
          out = in
        end if
      end if

    end function correct
    
  end subroutine fdf_brange
  
  subroutine fdf_bregions(bName, n, n_r, rgns)
    
    ! Name of block
    character(len=*), intent(in) :: bName
    ! Wrapper counter
    integer, intent(in) :: n
    ! Number of regions found
    integer, intent(out) :: n_r
    type(tRgn), intent(inout), allocatable :: rgns(:)

    ! ** local variables
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    type(tRgn) :: r1
    integer :: i, il, ic
    character(len=64) :: g
    character(len=64), allocatable :: rlist(:)
    logical :: found

    n_r = 0
    if ( allocated(rgns) ) deallocate(rgns)

    ! Get number of regions
    il = fdf_block_linecount(bName)
    if ( il == 0 ) return

    if ( .not. fdf_block(bName,bfdf) ) then
      call die('fdf_bregions: failed implementation.')
    end if

    ! allocate
    allocate(rlist(il))

    ! first count number of differently named regions
    do while ( fdf_bnext(bfdf,pline) ) 

      found = .false.
      if ( n_r > 0 ) then
        g = fdf_bnames(pline,1)
        do i = 1 , n_r
          if ( leqi(g,rlist(i)) ) then
            found = .true.
            exit
          end if
        end do
      end if
      if ( .not. found ) then
        n_r = n_r + 1
        rlist(n_r) = g
      end if

    end do

    ! Clean-up
    deallocate(rlist)
    call fdf_brewind(bfdf)
    
    allocate(rgns(0:n_r))
    
    il = 0
    do while ( fdf_bnext(bfdf,pline) ) 

      g = fdf_bnames(pline,1)

      ! Check if the name already has been read (then
      ! we accumulate the atoms)
      found = .false.
      ic = il + 1
      if ( il > 0 ) then
        do i = 1 , il
          if ( leqi(g,rgns(i)%name) ) then
            ic = i
            exit
          end if
        end do
      end if
      if ( ic == il + 1 ) then
        ! we have a new name
        il = ic

        ! We can read in a range
        call fdf_brange(pline,r1,1,n)
        if ( r1%n == 0 ) then
          print *,'Region: ',trim(g)
          call die('fdf_bregions: Could not read in anything in region!')
        end if
        call rgn_union(rgns(il),r1,rgns(il))
        rgns(il)%name = trim(g)

      else

        call fdf_brange(pline,r1,1,n)
        if ( r1%n == 0 ) then
          print *,'Region: ',trim(g)
          call die('fdf_bregions: Could not read in anything in region!')
        end if
        call rgn_union(rgns(ic),r1,rgns(ic))
        rgns(ic)%name = trim(g)
          
      end if
      
    end do
    
    call rgn_delete(r1)
    
  end subroutine fdf_bregions

end module fdf_extra
