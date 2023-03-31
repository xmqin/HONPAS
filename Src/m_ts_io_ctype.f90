!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_io_ctype
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use precision, only : dp

  implicit none

  integer, parameter :: C_N_NAME_LEN = 20

  integer, parameter :: c_N = 128
  type :: ts_c_io
     ! the name of the integration segment
     character(len=C_N_NAME_LEN) :: name = ' '
     ! The integral bounds
     real(dp) :: a, b
     ! (text-form)
     character(len=c_N) :: ca = ' ', cb = ' '
     ! The division of the grid
     real(dp) :: d
     ! (text-form)
     character(len=c_N) :: cd = ' '
     ! Number of points in the grid
     integer :: N = 0
     ! (text-form)
     character(len=c_N) :: cN = ' '
     ! the integration method (g-legendre etc.)
     character(len=c_N) :: method = ' '
     ! the type of integration (eq|neq|tran)
     character(len=4) :: type = ' '
     ! the integration part (circle|square|line|tail)
     character(len=c_N) :: part = ' '
     ! the options attached to the integration method
     ! (this will take the form of a linked list)
     type(ts_c_opt_ll), pointer :: opt => null()
  end type ts_c_io

  type ts_c_opt_ll
     character(len=c_N) :: opt = ' '
     character(len=c_N) :: val = ' '
     type(ts_c_opt_ll), pointer :: next => null()
  end type ts_c_opt_ll

  interface copy
     module procedure copy_
  end interface copy
  private :: copy_

  interface operator(.eq.)
     module procedure equal_
  end interface 
  private :: equal_

contains

  function equal_(a,b) result(eq)
    type(ts_c_io), intent(in) :: a, b
    logical :: eq

    eq = a%name == b%name

  end function equal_


  subroutine copy_(from,to)
    type(ts_c_io), intent(in) :: from
    type(ts_c_io), intent(out) :: to
    to%name = from%name
    to%a    = from%a
    to%b    = from%b
    to%ca   = from%ca
    to%cb   = from%cb
    to%d    = from%d
    to%cd   = from%cd
    to%N    = from%N
    to%cN   = from%cN
    to%method = from%method
    to%type = from%type
    to%part = from%part
    to%opt  => from%opt
  end subroutine copy_

  function fdf_nc_iotype(prefix,suffix) result(n)
    use fdf

    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: suffix
    integer :: n

    if ( len_trim(suffix) == 0 ) then
      n = fdf_block_linecount(trim(prefix)//'.Contours', 'n')
    else
      n = fdf_block_linecount(trim(prefix)//'.Contours'//trim(suffix), 'n')
    end if

  end function fdf_nc_iotype

  function fdf_name_c_iotype(prefix,suffix,i) result(name)
    use fdf
    character(len=*), intent(in) :: prefix, suffix
    integer, intent(in) :: i
    character(len=C_N_NAME_LEN) :: name

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: n
    logical :: found

    name = ' '

    n = 0
    if ( len_trim(suffix) == 0 ) then
       found = fdf_block(trim(prefix)//'.Contours',bfdf)
    else
       found = fdf_block(trim(prefix)//'.Contours.'//trim(suffix),bfdf)
    end if
    if ( .not. found ) return

    ! Find the name of the contour
    n = 0
    do while ( fdf_bline(bfdf,pline) )
      if ( fdf_bnnames(pline) == 0 ) cycle
      n = n + 1 
      if ( n == i ) then
        name = fdf_bnames(pline,1)
        call fdf_bclose(bfdf)
        return
      end if
    end do

  end function fdf_name_c_iotype

  subroutine c_io_add_opt(this,opt,val)
    use fdf, only : leqi
    type(ts_c_io), intent(in out) :: this
    character(len=*), intent(in) :: opt, val
    type(ts_c_opt_ll), pointer :: copt, new_opt

    if ( c_io_has_opt(this,opt) ) return
    
    nullify(new_opt)
    allocate(new_opt)
    new_opt%opt = trim(opt)
    new_opt%val = trim(val)

    if ( .not. associated(this%opt) ) then
       this%opt => new_opt
    else
       copt => this%opt
       do while ( associated(copt%next) )
          copt => copt%next
       end do
       copt%next => new_opt
    end if
    
  end subroutine c_io_add_opt

  function c_io_has_opt(this,opt) result(has)
    use fdf, only : leqi
    type(ts_c_io), intent(in) :: this
    character(len=*), intent(in) :: opt
    logical :: has
    type(ts_c_opt_ll), pointer :: copt

    copt => this%opt

    has = .false.
    do while ( associated(copt) )
       has = has .or. leqi(copt%opt,opt)
       copt => copt%next
    end do

  end function c_io_has_opt

  function c_io_get_opt(this,opt) result(res)
    use fdf, only : leqi
    type(ts_c_io), intent(in) :: this
    character(len=*), intent(in) :: opt
    character(len=c_N) :: res
    type(ts_c_opt_ll), pointer :: copt

    copt => this%opt

    res = ' '
    do while ( associated(copt) )
       if ( leqi(copt%opt,opt) ) then
          res = copt%val
          return
       end if
       copt => copt%next
    end do

  end function c_io_get_opt

  function ts_exists_contour_block(prefix,suffix,bName,bfdf) result(exist)
    use fdf

    character(len=*), intent(in) :: prefix, suffix
    character(len=C_N_NAME_LEN), intent(in) :: bName
    type(block_fdf), optional :: bfdf
    logical :: exist

    character(len=c_N) :: g

    ! if the block does not exist, return
    if ( len_trim(suffix) > 0 ) then
       g = trim(prefix)//'.Contour.'//trim(suffix)
    else
       g = trim(prefix)//'.Contour'
    end if
    g = trim(g) // '.' // trim(bName)
    if ( present(bfdf) ) then
       exist = fdf_block(trim(g),bfdf)
       if ( .not. exist ) then
          bfdf%label = trim(g)
       end if
    else
       exist = fdf_isblock(trim(g))
    end if

  end function ts_exists_contour_block

  subroutine ts_read_contour_block(prefix,suffix,bName,c, kT, V) 

    use units, only : eV
    use fdf
    use parse, only : search_fun, characters, ntokens

    character(len=*), intent(in) :: prefix, suffix
    character(len=C_N_NAME_LEN), intent(in) :: bName
    type(ts_c_io), intent(inout) :: c
    real(dp), intent(in) :: kT
    real(dp), intent(in), optional :: V

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()

    character(len=c_N) :: opt, val
    integer :: iS, iE
    logical :: is_file
    
    ! if the block does not exist, return
    if ( .not. ts_exists_contour_block(prefix,suffix,bName,bfdf=bfdf) ) then
       call die('Could not find block: '//trim(bfdf%label))
    end if

    ! the contour has already been assigned values
    if ( len_trim(c%name) /= 0 ) return

    ! Read the information in block
    c%name = trim(bName)
    
    ! We must ensure the block be organized as this:
    !  part circle|line|tail|user
    !  from <a> to <b> 
    !    points <N> / separation <d>
    !      method <method>
    !       opt <opt-1>
    !       opt <opt-2>
    !       opt ...
    !       opt <opt-N>

    ! { "part <part>" circle|square|line|tail
    ! This defaults to line
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find part segment in contour.')
    end if

    ! get the 'part'
    iS = search_fun('part',pline)
    if ( iS < 0 ) iS = search_fun('p',pline)
    if ( iS < 0 ) then
       c%part = 'line'
    else if ( fdf_bnnames(pline) < 2 ) then
       call die('Block: '//trim(bName)//' has not described the part properly. &
            &Must have part <part>')
    else
       c%part = fdf_bnames(pline,2)
    end if
    if ( leqi(c%part,'circle') ) then
       c%part = 'circle' ! ensures capitalization!!!! DON'T EDIT
    else if ( leqi(c%part,'square') ) then
       c%part = 'square'
    else if ( leqi(c%part,'line') ) then
       c%part = 'line'
    else if ( leqi(c%part,'tail') ) then
       c%part = 'tail'
    else if ( leqi(c%part,'user') ) then
       c%part = 'user'
    else
       call die('Part of the contour could not be recognized as circle|square|line|tail|user')
    end if
   

    ! } "part

    ! { "from <a> to <b>"
    if ( iS > 0 ) then
       if ( .not. move2names() ) then
          call die('Block: '//trim(bName)//'. &
               &Could not find from <a> to <b> segment in contour')
       end if
    end if

    ! get the 'from'
    iS = search_fun('from',pline)
    if ( iS <= 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &from <a> to <b> line cannot find "from".')
    end if
    
    iE = search_fun('to',pline)
    if ( iE <= 0 ) then
       call die('Block: '//trim(bName)//' is not build correctly. &
            &from <a> to <b> line cannot find "to".')
    end if
    ! get the line as input
    call pline_E_parse(pline,iS,c%ca,val=c%a,V=V,kT=kT, &
         before=iE)
    if ( leqi(c%ca,'next') ) call die('Block: '//trim(bName)//' can not &
         &have a==next')

    ! get <b>
    iS = iE
    call pline_E_parse(pline,iS,c%cb,val=c%b,V=V,kT=kT)
    if ( leqi(c%cb,'previous') .or. & 
         leqi(c%cb,'prev') ) call die('Block: '//trim(bName)//' can not &
         &have b==previous')

    ! } "from"

    ! { "points <N>"
    if ( .not. move2names() ) then
       call die('Block: '//trim(bName)//'. &
            &Could not find points segment in contour')
    end if
    
    ! we now read the points or separation
    is_file = .false.
    iS = search_fun('points',pline)
    if ( iS < 0 ) iS = search_fun('p',pline)
    iE = search_fun('delta',pline)
    if ( iE < 0 ) iE = search_fun('d',pline)
    if ( iS < 0 .and. iE < 0 ) then
      iS = search_fun('file',pline)
      if ( iS < 0 ) iS = search_fun('user',pline)
      if ( iS < 0 ) then
        call die('Block: '//trim(bName)//' is not build correctly. &
            &Could not decipher points/delta/file')
      end if
      is_file = .true.
    end if
    
    ! if we have points we simply read in the number
    if ( is_file ) then
      
      ! The file contains all the relevant information.
      c%cN = fdf_bnames(pline, 2) ! file name
      c%method = 'user' ! user-defined method
      c%N = file_energy_points(c%cN) ! read number of energy points

    else if ( 0 <= iS ) then
       c%N  = fdf_bintegers(pline,1,after=iS) ! first integer
       c%cN = characters(pline,1,-1,after=iS)
       if ( c%N < 1 ) then
          call die('Block: '//trim(bName)//' is not valid. &
               &A negative amount of integration points is not &
               &a valid input.')
       end if
    else ! we have a delta designation
       ! notice that we can actually use kT here
       call pline_E_parse(pline,iE,c%cd,val=c%d,kT=kT)
       if ( c%d <= 0._dp ) then
          call die('Block: '//trim(bName)//' is not valid. &
               &The dE designator is negative or zero.')
       end if
       ! Ensure that "inf" is not one of the end-points
       if ( index(c%ca,'inf') > 0 .or. &
            index(c%cb,'inf') > 0 ) then
          call die('Block: '//trim(bName)//' is not valid. &
               &The dE designator is not allowed with an "inf" boundary.')
       end if
    end if

    ! } "points"

    ! { "method <method>"
    if ( move2names() ) then
      val = fdf_bnames(pline, 1)
      if ( leqi(val, 'method') ) then
        if ( fdf_bnnames(pline) < 2 ) then
          call die('Block: '//trim(bName)//' has not described the method properly. &
              &Must have method <method>')
        end if

        ! the method should be a one-name thing
        c%method = fdf_bnames(pline,2)
        
      else if ( .not. is_file ) then
        ! The method segment has not been found.
        ! Error out
        call die('Block: '//trim(bName)//'. &
            &Could not find method <method> segment in contour')
      else if ( is_file ) then
        ! The name is most probably an option, step back in the block
        if ( .not. fdf_bbackspace(bfdf) ) then
          call die('Block: '//trim(bName)//' parsing went wrong!')
        end if
      end if
    end if
    if ( is_file ) c%method = 'user'

    ! } "method"

    ! { "opt <option>"
    nullify(c%opt)
    do 
       ! if we don't find anything simply exit the optional reading
       if ( .not. move2names() ) exit
       
       opt = trim(characters(pline,1,1,after=1))
       if ( ntokens(pline,after=1) > 1 ) then
          val = trim(characters(pline,1,1,after=2))
          do iE = 2 , ntokens(pline,after=2)
             val = trim(val)//' '//trim(characters(pline,iE,iE,after=2))
          end do
       else
          val = ' '
       end if

       call c_io_add_opt(c,opt,val)

    end do
    ! } "opt"
    
  contains

    function file_energy_points(file) result(ne)
      use m_io, only: io_assign, io_close
      
      character(len=*), intent(in) :: file
      integer :: iu
      character(len=256) :: line
      integer :: ne, iostat
      logical :: exist

      call io_assign(iu)

      ! Open the contour file
      inquire(file=trim(file), exist=exist)
      if ( .not. exist ) then
        call die('Block: '//trim(bName)//' requested an external contour file. &
            &The file: '//trim(file)//' could not be found!')
      end if

      ! Open the file
      open(iu, file=trim(file), form='formatted', status='old')

      ne = 0
      do
        read(iu, '(a)', iostat=iostat) line
        if ( iostat /= 0 ) exit
        if ( len_trim(line) == 0 ) cycle
        line = trim(adjustl(line))
        if ( line(1:1) == '#' ) cycle
        
        ne = ne + 1
      end do

      call io_close(iu)
      
    end function file_energy_points

    function move2names() result(found)
      logical :: found
      found = .false.
      do while ( fdf_bline(bfdf,pline) ) 
       
         ! if no names exist we can loop to the next line
         ! thus commenting out a contour works very good!
         if ( fdf_bnnames(pline) == 0 ) then
            cycle
         else
            found = .true.
            exit
         end if
      end do
      
    end function move2names

  end subroutine ts_read_contour_block

  function ts_c_bisphysical(pline,ind)
    use parse
    implicit none
!------------------------------------------------- Input Variables
    type(parsed_line), pointer        :: pline
    integer(ip), intent(in)           :: ind

!------------------------------------------------ Output Variables
    logical                           :: ts_c_bisphysical
    
    ts_c_bisphysical = match(pline,'vn',after=ind)

  end function ts_c_bisphysical
    
  function ts_c_bphysical(pline,ind,defunit)
    use fdf
    use parse, only : die
    
    implicit none
!------------------------------------------------- Input Variables
    type(parsed_line), pointer        :: pline
    integer,     intent(in)           :: ind
    character(len=*), intent(in)      :: defunit

!------------------------------------------------ Output Variables
    real(dp)                          :: ts_c_bphysical
    character(len=10)                 :: unitstr

    if (.not. ts_c_bisphysical(pline,ind) ) then
       call die('PARSE module: ts_c_bphysical', 'Not enough values and names in line',   &
            'm_ts_io_ctype', -1)
    endif
    
    ! Label with value
    ts_c_bphysical = fdf_bvalues(pline, 1, ind)
    
    unitstr = fdf_bnames(pline, 1, ind+1)
    if (.not. leqi(unitstr, defunit)) &
         ts_c_bphysical = ts_c_bphysical * fdf_convfac(unitstr, defunit)
    
  end function ts_c_bphysical

  subroutine pline_E_parse(pline,after,c, &
       val,V,kT, &
       before)
    use fdf
    use parse, only : characters, ntokens
    use units, only : eV
    
    ! This routine will parse a pline from an index and collect the values from units
    ! into one value and return the value and the string which can build it again
    type(parsed_line), pointer :: pline
    integer, intent(in) :: after ! search physical values from after this index
    character(len=c_N), intent(out) :: c ! the string that can construct this physical value in the shortest sense
    ! value of the energy
    real(dp), intent(out), optional :: val
    real(dp), intent(in),  optional :: V, kT
    integer, intent(in), optional   :: before

    ! Local parameters
    real(dp) :: tmp
    character(len=c_N) :: g
    logical :: add, get_val, has_V, has_kT, absolute
    integer :: i, j, offset

    ! initialize the output string
    c = ' '
    ! whether we should collect the value
    get_val = present(val)
    if ( get_val ) then
       val = 0._dp  ! initialize
       has_V = present(V)
       has_kT = present(kT)
    end if

    i = after - 1
!    if ( present(before) ) then
!       write(*,*)'CHARS: ',trim(characters(pline,after+1,before-1))
!    else
!       write(*,*)'CHARS: ',trim(characters(pline,1,-1,after=after))
!    end if
    add = .true.

    do while ( ntokens(pline,after=i) > 1 ) 
       i = i + 1
       if ( present(before) ) then
          !print *,'Before: ',i,before
          if ( i >= before ) then
             c = adjustl(c)
             return
          end if
       end if
         
       ! first we need to figure out if we are dealing with 
       ! a single plus or negative sign...
       g = fdf_bnames(pline,1,after=i)
!       print *,'CHECK: ',trim(g),' ',trim(characters(pline,1,1,after=i))
       offset = 2
       if ( trim(characters(pline,1,1,after=i)) == trim(g) ) then
          offset = 1 ! for units the unit is at the index
       else if ( ntokens(pline,after=i) > 0 ) then
!          print *,'CHECK2: ',trim(g),' ',trim(characters(pline,1,1,after=i+1))
          if ( trim(characters(pline,1,1,after=i+1)) == trim(g) ) then
             offset = 0
          end if
       end if
!print *,'PE: ',trim(g),' offset:',offset
       if ( leqi(g,'+') ) then
          add = .true.
          c = trim(c)//' +'
          cycle
          
       else if ( leqi(g,'-') ) then
          add = .false.
          c = trim(c)//' -'
          cycle
          
       end if

       ! else we need to figure out if it is something which is not a specific unit...
         
       !print *,'Tokens: ',trim(tokens(pline,1,after=i)),' text: ',match(pline,'v',after=i),i,trim(g)
       ! Parse compressed elements
       if ( leqi(g(1:1),'+') ) then
          add = .true.
          c = trim(c)//' +'
          g = g(2:)

       else if ( leqi(g(1:1),'-') ) then
          add = .false.
          c = trim(c)//' -'
          g = g(2:)

       end if

       ! we need to be sure to capture everything...
       if ( leqi(g,'prev') .or. &
            leqi(g,'previous') .or. &
            leqi(g,'next') ) then
          c = trim(c) // ' ' //trim(g)

       else if ( leqi(g,'inf') ) then
          ! note: we have already added sign

          if ( add ) then
             if ( get_val ) val = huge(1._dp)
             c = ' + inf'
          else
             if ( get_val ) val = -huge(1._dp)
             c = ' - inf'
          end if

          ! you can't add/subtract anything meaningful to
          ! "inf" without getting inf again...
          c = adjustl(c)
          return

       else if ( leqi(g(1:1),'v') .or. leqi(g(1:2),'|v') ) then

          ! distinguish between the absolute part or the regular V
          absolute = leqi(g(1:1),'|')
          if ( absolute ) then
             if ( .not. leqi(g(3:3),'|') ) then
                call die('Requesting the absolute bias point requires it to &
                     &be formatted as: |V|/<fraction>')
             end if
             ! remove the | signs
             g = g(2:2)//trim(g(4:))
          end if

          ! Get fraction
          j = parse_frac(g)
          if ( j > 1 ) then
             if ( absolute ) then
                write(g,'(a,i0)') '|V|/',j
             else
                write(g,'(a,i0)') 'V/',j
             end if
          else
             if ( absolute ) then
                g = '|V|'
             else
                g =  'V'
             end if
          end if

          ! get the V designation
          c = trim(c)//' '//trim(g)
          
          if ( get_val ) then
             if ( .not. has_V ) call die('You cannot request &
                  &V value in this segment')
             
             ! Notice that this is with respect to the sign of the bias...
             ! the sign of the bias is important, ONLY when dealing with
             ! 2 electrodes
             if ( absolute ) then
                tmp = abs(V) / real(j,dp)
             else
                tmp = V / real(j,dp)
             end if
             !print*,'Found : '//trim(g)//' at value: ',tmp/eV

             if ( add ) then
                val = val + tmp
             else
                val = val - tmp
             end if

          end if

       else

          if ( offset > 1 ) then
             c = adjustl(c)
             return ! we have too many things unrecognizable
          end if

          ! We are dealing with a proper unit...
          ! check that we have a physical quantity that can be described...
          if ( .not. ts_c_bisphysical(pline,i-offset) ) then
             !print *,'NOT PHY: ',characters(pline,1,2,after=i-offset)
             c = adjustl(c)
             return ! if we can't find
          end if
          
          if ( leqi(g(1:2),'kt') .or. leqi(g(1:3),'kbt') ) then
             
             ! Get fraction and parse to correct format
             j = parse_frac(g)
             if ( j > 1 ) then
                write(g,'(a,i0)') ' kT/',j
             else
                g = ' kT'
             end if

             if ( get_val .and. .not. has_kT ) call die('You cannot request &
                  &kT value in this segment')

             ! Read in the former value
             ! We should not remove pm as that will not be read in
             g = trim(characters(pline,1,1,after=i-offset))//trim(g)
             
             c = trim(c)//' '//trim(g)
             
             ! get the value...
             tmp = fdf_bvalues(pline,1,after=i-offset) * kT / j
          else

             ! this will make the code break if the unit is wrongly assigned...
             tmp = ts_c_bphysical(pline,i-offset,'Ry')

             c = trim(c)//' '//trim(characters(pline,1,1,after=i-offset))//' '&
                  //trim(characters(pline,2,2,after=i-offset))

          end if

          if ( get_val ) then
             if ( add ) then
                val = val + tmp
             else
                if ( tmp > 0._dp ) then
                   val = val - tmp
                else
                   val = val + tmp
                end if
             end if
          end if


          ! skip both the value and unit...
          if ( offset == 0 ) then
             i = i + 1
          end if

       end if

       ! we must have the resetting here, otherwise single +- statements
       ! will not be caught
       add = .true.
    end do

    ! be sure to adjustl c
    c = adjustl(c)
    
  contains
    
    function remove_pm(s)
      character(len=*), intent(in) :: s
      character(len=len_trim(s)) :: remove_pm
      if ( leqi(s(1:1),'+') .or. &
           leqi(s(1:1),'-') ) then
         remove_pm = trim(s(2:))
      else
         remove_pm = trim(s)
      end if
    end function remove_pm

    function parse_frac(s) result(frac)
      character(len=*), intent(in) :: s
      integer :: frac
      integer :: i, stat

      ! Always divide by 1 (no fraction)
      frac = 1
      ! parenthesis position
      i = index(s,'/')

      if ( i < 1 ) return

      read(s(i+1:),'(i9)',iostat=stat) frac
      if ( stat /= 0 ) then
         call die('Fractional parameter chosen cannot be distinguished: '// &
              trim(s)//' expecting '//trim(s(i:)))
      end if
      
    end function parse_frac

  end subroutine pline_E_parse

  subroutine ts_print_contour_block(prefix,c)

    use parallel, only : IONode

    character(len=*), intent(in) :: prefix
    type(ts_c_io), intent(in) :: c
    type(ts_c_opt_ll), pointer :: opt
    character(len=300) :: name

    if ( .not. IONode ) return

    if ( c%name(1:1) == '*' ) then
       name = trim(prefix) // trim(c%name(2:))
    else
       name = trim(prefix) // trim(c%name)
    end if

    ! Start by writing out the block beginning
    write(*,'(a,a)') '%block ',trim(name)

    ! For transport contours we are ALWAYS
    ! using the LINE contour
    if ( trim(c%type) /= 'tran' ) then
       write(*,'(2a)') '  part ',trim(c%part)
    end if

    ! move to third column...
    write(*,'(t3)',advance='no')
      
    write(*,'(4a)') '   from ',trim(c%ca),' to ',trim(c%cb)

    ! Ensure we correctly get the user-definition
    if ( c%method == 'user' ) then
      write(*,'(t7,a,tr1,a)') 'file', trim(c%cN)
      
    else if ( len_trim(c%cd) /= 0 ) then
      ! we have delta designation
      write(*,'(t7,a,tr1,a)') 'delta', trim(c%cd)
      
    else       
      ! Print the number of points...
      write(*,'(t7,a,tr1,i0)') 'points', c%N
      
    end if

    ! Print the method
    write(*,'(t9,a,tr1,a)') 'method', trim(c%method)

    opt => c%opt
    do while ( associated(opt) )
       if ( len_trim(opt%val) > 0 ) then
          write(*,'(t10,a,2(tr1,a))') 'opt', trim(opt%opt),trim(opt%val)
       else
          write(*,'(t10,a,tr1,a)') 'opt', trim(opt%opt)
       end if
       opt => opt%next
    end do
    
    write(*,'(a,a)') '%endblock ',trim(name)
    
  end subroutine ts_print_contour_block


  ! *****
  ! This routine fixes the inputs for the contours according to those given by 
  ! the input electrode
  ! 1.) It fixes the bounds next to each other if they have designated
  !     'next' or 'previous'
  ! 2.) Set the method to be equilibrium/non-equilibrium/transport
  ! 3.) Calculate number of points if dE specified
  ! 4.) Checks whether the contours are connected
  subroutine ts_fix_contour(cur,next,prev,connected)
    use fdf, only : leqi
    type(ts_c_io), intent(inout) :: cur
    type(ts_c_io), intent(inout), optional :: next, prev
    logical, intent(out), optional :: connected
    integer :: same_ref

    if ( present(connected) ) connected = .true.

    ! we need this to "look ahead"
    if ( present(next) ) then
       if ( leqi(next%ca,'prev') .or. leqi(next%ca,'previous') ) then
          next%a = cur%b
          if ( leqi(cur%cb,'next') ) then
             write(*,'(a)') 'Current contour have next: '//trim(cur%name)
             write(*,'(a)') 'Next contour have prev: '//trim(next%name)
             call die('Connecting two contours by next and prev is invalid. &
                  &An explicit value is needed in one of the places.')
          end if
       end if
    end if

    if ( leqi(cur%cb,'next') ) then
       if ( present(next) ) then
          cur%b = next%a
       else
          write(*,'(a)') 'Erroneous contour: '//trim(cur%name)
          call die('The contour segment is not attached &
               &to a following segment (next does not exist).')
       end if
    end if

    if ( leqi(cur%ca,'prev') .or. leqi(cur%ca,'previous') ) then
       if ( present(prev) ) then
          cur%a = prev%b
       else
          write(*,'(a)') 'Erroneous contour: '//trim(cur%name)
          call die('The contour segment is not attached &
               &to a previous segment (prev does not exist).')
       end if
    end if

    ! Correct if 'b' connects 'a' through prev, or 'a' to 'b' 
    ! through next
    same_ref = 0
    if ( index(cur%ca,'next') > 0 ) then
       cur%a = cur%a + cur%b
       same_ref = same_ref + 1
    end if
    if ( index(cur%cb,'prev') > 0 .or. index(cur%cb,'previous') > 0 ) then
       cur%b = cur%b + cur%a
       same_ref = same_ref + 1
    end if
    if ( same_ref > 1 ) then
       write(*,'(a)') 'Erroneous contour: '//trim(cur%name)
       call die('The contour segment must not reference it self, &
            &next to prev is NOT allowed.')
    end if

    ! we can compare bounds
    if ( present(connected) ) then
       connected = c_io_conn(cur,prev=prev,next=next)
    end if
    
    ! at this point both boundaries MUST exist
    if ( len_trim(cur%cd) > 0 .and. len_trim(cur%cN) == 0 ) then
       cur%N = nint(abs(cur%b - cur%a)/cur%d)
    end if
         
  end subroutine ts_fix_contour

  function c_io_conn(cur,next,prev) result(conn)
    type(ts_c_io), intent(inout) :: cur
    type(ts_c_io), intent(inout), optional :: next, prev
    logical :: conn

    conn = .true.
    
    if ( present(prev) ) then
       if ( abs(cur%a - prev%b) > 1.e-8_dp ) then
          conn = .false.
       end if
    end if
    if ( present(next) ) then
       if ( abs(next%a - cur%b) > 1.e-8_dp ) then
          conn = .false.
       end if
    end if
    
    ! We do not allow single points
    if ( abs(cur%a - cur%b) < 1.e-8_dp ) then
       conn = .false.
    end if
    
  end function c_io_conn

  subroutine ts_fix_contours(N,c_io,Ni,idx)

    use parallel, only: IONode
    use units, only : eV
    
    integer, intent(in) :: N
    type(ts_c_io), intent(inout) :: c_io(N)
    integer, intent(in) :: Ni
    integer, intent(in) :: idx(Ni)

    integer :: i, start, j
    integer :: cur, next, prev
    logical :: done

    if ( Ni == 1 ) then
       call ts_fix_contour( c_io(idx(1)) )
       return
    end if
    
    do i = 1 , Ni
       cur = idx(i)
       if ( 1 < i .and. i < Ni ) then
          prev = idx(i-1)
          next = idx(i+1)
          call ts_fix_contour( c_io(cur), &
               prev=c_io(prev), next=c_io(next) )
       else if ( 1 < i ) then
          prev = idx(i-1)
          call ts_fix_contour( c_io(cur), &
               prev=c_io(prev) )
       else if ( i < Ni ) then
          next = idx(i+1)
          call ts_fix_contour( c_io(cur), &
               next=c_io(next) )
       else
          call die('error')
       end if
    end do

    call update_done()
    
    if ( .not. done ) then

       ! We try to track next
       start = 1
       do while ( start <= Ni )

          i = start
          do
             cur = idx(i)
             if ( index(c_io(cur)%cb,'next') == 0 ) exit
             i = i + 1
          end do

          do j = i , start , -1
             cur = idx(j)
             ! only track forward (we are following)
             if ( j < Ni ) then
                next = idx(j+1)
                call ts_fix_contour( c_io(cur), &
                     next=c_io(next) )
             else if ( 1 < j ) then
                prev = idx(j-1)
                call ts_fix_contour( c_io(cur) , &
                     prev=c_io(prev) )
             else
                call ts_fix_contour( c_io(cur) )
             end if
             
          end do
          
          start = i + 1
       end do

    end if
       
    call update_done()
    
    if ( .not. done ) then
       if ( IONode ) then
          ! printing the contour segments
          do i = 1 , Ni
             call ts_print_contour_block('ERROR.TS.',c_io(idx(i)))
          end do
          do i = 1 , Ni
             write(*,*)c_io(idx(i))%a/eV,c_io(idx(i))%b/eV
          end do
       end if
       call die('Error in contour continuity.')
    end if

  contains

    subroutine update_done()
      
      ! check connection
      done = .true.

      cur = idx(1)
      if ( Ni > 1 ) then
         next = idx(2)
         done = done .and. &
              c_io_conn(c_io(cur), next=c_io(next))
      end if
      do i = 2 , Ni - 1
         prev = idx(i-1)
         cur  = idx(i)
         next = idx(i+1)
         done = done .and. &
              c_io_conn(c_io(cur), &
              prev=c_io(prev), next=c_io(next))
      end do
      cur = idx(Ni)
      if ( Ni > 1 ) then
         prev = idx(Ni-1)
         done = done .and. &
              c_io_conn(c_io(cur), prev=c_io(prev))
      end if

    end subroutine update_done
    
  end subroutine ts_fix_contours
  
end module m_ts_io_ctype
