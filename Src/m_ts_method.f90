module m_ts_method

  use m_region

  implicit none

  public
  save

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA.
  integer, parameter :: TS_FULL = 1

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the heavily optimized tri-diagonalization
  integer, parameter :: TS_BTD = 2

  ! This is the transiesta version utilizing the 
  ! full sparsity pattern of SIESTA as well 
  ! as the MUMPS library
  integer, parameter :: TS_MUMPS = 3

  ! The default solution method (it will be reset
  ! after option reading)
  integer :: ts_method = TS_BTD

  ! The BTD spectral function calcultation method
  ! Currently we only do propagation and column methods
  ! This may be abstracted later for other variants
  integer, parameter :: TS_BTD_A_PROPAGATION = 0
  integer, parameter :: TS_BTD_A_COLUMN = 1

  ! The actual method used 
  integer :: ts_A_method = TS_BTD_A_PROPAGATION

  ! The buffer atoms have type = -1, dev = 0, Electrodes = E_idx
  integer, parameter :: TYP_BUFFER = -1
  integer, parameter :: TYP_DEVICE = 0

  ! Containers for which atom/orbital is what
  integer, private :: no_u_TS = 0
  integer, allocatable, private :: a_type(:)
  integer, allocatable, private :: o_type(:)
  integer, allocatable, private :: a_offset(:)
  integer, allocatable, private :: o_offset(:)

  ! Containers for r_[ao]Buf%n for easier handling.
  integer :: na_Buf = 0
  integer :: no_Buf = 0
  type(tRgn) :: r_aBuf, r_oBuf ! the buffer region
  type(tRgn) :: r_aC,   r_oC   ! the entire calculation region

  ! We create a table for converting transiesta orbitals to 
  ! siesta orbitals:
  !   s_io = r_pvt%r(ts_io)
  ! This allows to easily convert from one to the other
  type(tRgn) :: r_pvt

  ! Local routine
  private :: set_type

contains

  subroutine ts_init_regions(prefix, na_u, lasto)

    use alloc
    use fdf
    use fdf_extra, only : fdf_brange

    character(len=*), intent(in) :: prefix
    integer, intent(in) :: na_u, lasto(0:na_u)

    integer :: i, ia

    ! prepare to read in the data...
    character(len=len_trim(prefix)+13) :: bName
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    character(len=32) :: g
    type(tRgn) :: r_tmp, r_diff

    no_u_TS = lasto(na_u)

    ! allocate regions
    if ( allocated(a_type) ) then

      deallocate(a_type,a_offset,o_type,o_offset)

      ! Clean regions
      call rgn_delete(r_aBuf,r_oBuf,r_aC,r_oC,r_pvt)

    end if

    allocate(a_type(na_u)   ,a_offset(na_u)    )
    allocate(o_type(no_u_TS),o_offset(no_u_TS) )
    a_type(:)   = TYP_DEVICE
    a_offset(:) = 0
    o_type(:)   = TYP_DEVICE
    o_offset(:) = 0

    ! old options for buffer atoms
    call fdf_obsolete('TS.BufferAtomsLeft')
    call fdf_obsolete('TS.BufferAtomsRight')

    ! Read in TS.Atoms.Buffer
    bName = trim(prefix)//'.Atoms.Buffer'
    ! Read in buffer region via the block or list
    if ( fdf_islist(bName) ) then

      ! Query size of list
      i = -1
      call rgn_init(r_aBuf, 1)
      call fdf_list(bName, i, r_aBuf%r)
      call rgn_init(r_aBuf, i)
      call fdf_list(bName, r_aBuf%n, r_aBuf%r)
      do i = 1 , r_aBuf%n
        if ( r_aBuf%r(i) < 0 ) then
          r_aBuf%r(i) = r_aBuf%r(i) + na_u + 1
        end if
      end do

    else if ( fdf_block(bName,bfdf) ) then

      ! read by line and set them to be buffer atoms
      do while ( fdf_bline(bfdf,pline) ) 

        ! empty line
        if ( fdf_bnnames(pline) == 0 ) cycle

        g = fdf_bnames(pline,1)
        if ( leqi(g,'atom') .or. leqi(g,'position') ) then

          call fdf_brange(pline, r_tmp, 1, na_u)
          if ( r_tmp%n == 0 ) &
              call die('Could not read in any atoms &
              &in line of TS.Atoms.Buffer')
          call rgn_union(r_aBuf,r_tmp,r_aBuf)

        end if

        ! Atoms NOT in the buffer region...
        if ( leqi(g,'not-atom') .or. leqi(g,'not-position') .or. &
            leqi(g,'-atom') .or. leqi(g,'-position') ) then
          call fdf_brange(pline, r_tmp, 1, na_u)
          if ( r_tmp%n == 0 ) &
              call die('Could not read in any atoms &
              &in line of TS.Atoms.Buffer')
          call rgn_union(r_diff,r_tmp,r_diff)
          
        end if

      end do

    end if

    if ( r_diff%n > 0 ) then
      ! Remove any atoms explicitly forced to not be part of it
      call rgn_complement(r_diff, r_aBuf, r_aBuf)
    end if
    call rgn_delete(r_diff)
    
    if ( r_aBuf%n > 0 ) then
      ! Remove all 0 entries
      ! This is a "helper" because the user may
      ! request -10 -- 10
      call rgn_init(r_tmp,1,val=0)
      call rgn_complement(r_tmp,r_aBuf,r_aBuf)
      call rgn_uniq(r_aBuf)
      do i = 1 , r_aBuf%n
        call set_type(TYP_BUFFER,r_aBuf%r(i),na_u,lasto)
      end do
    end if
    call rgn_delete(r_tmp)

    ! Sort the atoms
    call rgn_sort(r_aBuf)

    ! Create the calculation region
    call rgn_range(r_aC, 1, na_u)
    call rgn_complement(r_aBuf, r_aC, r_aC)

    ! Convert atom regions to orbital regions
    call rgn_Atom2Orb(r_aBuf,na_u,lasto,r_oBuf)
    call rgn_Atom2Orb(r_aC,na_u,lasto,r_oC)

    ! Update counting buffers
    na_Buf = r_aBuf%n
    no_Buf = r_oBuf%n

    ! Create the "pivoting" array
    call rgn_init(r_pvt,lasto(na_u)-r_oBuf%n)
    ia = 0
    do i = 1 , lasto(na_u)
      if ( orb_type(i) == TYP_BUFFER ) cycle
      ia = ia + 1
      if ( ia > r_pvt%n ) call die('Error in programming!')
      r_pvt%r(ia) = i
    end do

    ! Name the regions
    r_aBuf%name = '[A]-buffer'
    r_oBuf%name = '[O]-buffer'

    r_aC%name = '[A]-calculation'
    r_oC%name = '[O]-calculation'

    r_pvt%name = '[O]-pivot'

  end subroutine ts_init_regions

  subroutine ts_init_electrodes(na_u,lasto,N_Elec,Elecs)
    use m_ts_electype

    integer,    intent(in) :: na_u, lasto(0:na_u)
    integer,    intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer :: iE, ia

    do iE = 1 , N_Elec
      do ia = 1 , TotUsedAtoms(Elecs(iE))
        call set_type(iE,Elecs(iE)%idx_a - 1 + ia,na_u,lasto)
      end do
    end do

  end subroutine ts_init_electrodes

  subroutine set_type(typ,ia,na_u,lasto)
    integer, intent(in) :: typ, ia, na_u,lasto(0:na_u)
    integer :: i, no
    if ( ia > na_u ) then
      call die('Error in specifying the type of an atom!. &
          &Atoms specified is above the total number of atoms!')
    end if
    if ( a_type(ia) /= TYP_DEVICE ) then
      write(*,'(2(a,i0))') 'Trying to set atom ',ia,' to type: ',typ
      write(*,'(2(a,i0))') 'Atom ',ia,' is already: ',a_type(ia)

      call die('Error in setup. Atoms are having two types, check for &
          &electrode and buffer atom overlap...')
    end if
    a_type(ia) = typ
    o_type(lasto(ia-1)+1:lasto(ia)) = typ
    if ( typ == TYP_BUFFER ) then
      do i = ia , na_u
        a_offset(i) = a_offset(i) + 1
      end do
      no = lasto(ia) - lasto(ia-1)
      do i = lasto(ia-1) + 1 , lasto(na_u)
        o_offset(i) = o_offset(i) + no
      end do
    end if
  end subroutine set_type

  elemental function orb_type(io) result(typ)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: typ
    typ = o_type(ucorb(io,no_u_TS))
  end function orb_type

  elemental function orb_offset(io) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: io
    integer :: off
    off = o_offset(ucorb(io,no_u_TS))
  end function orb_offset

  elemental function ts2s_orb(io) result(off)
    integer, intent(in) :: io
    integer :: off
    do off = io + o_offset(io) , no_u_TS
      if ( o_type(off) == TYP_BUFFER ) cycle ! the buffer atoms are NOT transiesta
      if ( off - o_offset(off) == io ) return
    end do
  end function ts2s_orb

  elemental function atom_offset(ia) result(off)
    use geom_helper, only : UCORB
    integer, intent(in) :: ia
    integer :: off
    off = a_offset(ia)
  end function atom_offset

  elemental function atom_type(ia) result(typ)
    integer, intent(in) :: ia
    integer :: typ
    typ = a_type(ia)
  end function atom_type

  elemental function a_isBuffer(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) == TYP_BUFFER
  end function a_isBuffer

  elemental function a_isElec(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) > 0
  end function a_isElec

  elemental function a_isDev(ia) result(typ)
    integer, intent(in) :: ia
    logical :: typ
    typ = a_type(ia) == TYP_DEVICE
  end function a_isDev

end module m_ts_method

