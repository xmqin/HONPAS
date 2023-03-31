! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This program has been fully implemented by:
!  Nick Papior, 2015

! This pivoting program will pivot atoms in the sparsity
! pattern by reading in a sparsity pattern from
! a given SIESTA file.
! Subsequently the sparsity pattern will be printed
! out in any of the provided formats.
! Currently the output formats are:
! 1. METIS graph
! 2. GRAPHVIZ input
program pvtsp

  use precision, only : dp

  use class_OrbitalDistribution
  use class_Sparsity

  use geom_helper, only : UCORB
  use m_io_s
  use m_os, only : file_exist
  use create_sparsity_SC

  use m_region

  use m_pivot
  use m_pivot_methods

  ! Get the transfer matrix possibility
  use create_sparsity_SC

  use m_sparsity_handling

  character(len=264) :: fname, arg
  
  ! The sparsity
  type(OrbitalDistribution) :: fdit
  type(Sparsity) :: sp_full, sp_uc

  integer :: no_u, na_u
  integer, pointer :: ncol(:), l_ptr(:), l_col(:)

  ! Parameters for the unit-cell, if the reading file
  ! is a TSHS file
  real(dp) :: uc(3,3)
  integer, allocatable :: lasto(:), isc_off(:,:)

  integer :: N_arg, i_arg, i, o1, o2
  real :: r

  character(len=64) :: fmethod
  integer :: method
  type(tRgn) :: pvt
  
  logical :: is_graphviz, is_metis, is_atom
  logical :: is_unit_cell
  logical :: has_weight

  integer :: tm(3), graph_method
  integer :: n_nzs

  method = -1
  graph_method = 1 ! regular graph
  is_unit_cell = .false.
  is_atom = .false.
  is_metis = .false.
  is_graphviz = .false.
  has_weight = .false.
  tm(:) = TM_ALL

  N_arg = command_argument_count()
  ! Get the file-name from the command-line
  if ( N_arg < 1 ) then
     stop "Call this program with a file-name *.DM, *.TSDE, *.TSHS"
  end if
  i_arg = 1
  do while ( i_arg <= N_arg )
    
     call get_command_argument(i_arg, value = arg)
     arg = trim_em(arg)

     select case ( arg )
     case ( '-out', '-pvt' )
       ! Specify the pivoting method (equivalent -out gps to -gps for instance)
       i_arg = i_arg + 1
       call get_command_argument(i_arg, value = arg)
       arg = '-' // trim(arg)
     end select

     
     select case ( arg )

       ! Output type:
     case ( '-graphviz', '-graph', '-viz' )
       is_graphviz = .true.
     case ( '-metis-stdout' )
       is_metis = .true.
       is_unit_cell = .true.

     case ( '-atom' )
       ! Check if it is atomic sparsity pattern
       is_atom = .true.

     case ( '-di', '-digraph' )
       ! Check if it is an digraph
       is_graphviz = .true.
       graph_method = 2

       ! Check the method for pivoting
     case ( '-cm' )
       method = PVT_CUTHILL_MCKEE
        fmethod = 'CuthillMckee'
      case ( '-rev-cm' )
        method = PVT_REV_CUTHILL_MCKEE
        fmethod = 'revCuthillMckee'
      case ( '-gps' )
        method = PVT_GPS
        fmethod = 'GibbsPooleStockmeyer'
      case ( '-rev-gps' )
        method = PVT_REV_GPS
        fmethod = 'revGibbsPooleStockmeyer'
      case ( '-pcg' )
        method = PVT_PCG
        fmethod = 'PeripheralConnectGraph'
      case ( '-rev-pcg' )
        method = PVT_REV_PCG
        fmethod = 'revPeripheralConnectGraph'
      case ( '-ggps' )
        method = PVT_GGPS
        fmethod = 'GeneralGibbsPooleStockmeyer'
      case ( '-rev-ggps' )
        method = PVT_REV_GGPS
        fmethod = 'revGeneralGibbsPooleStockmeyer'
      case ( '-scramble' )
        method = 0 ! 0 signals scrambling
        fmethod = 'Scramble'
        
#ifdef SIESTA__METIS
      case ( '-metis', '-nodend' )
        method = PVT_METIS_NODEND
        fmethod = 'metis-NodeND'
      case ( '-partgraphkway' )
        method = PVT_METIS_PARTGRAPHKWAY
        fmethod = 'metis-PartGraphKway'
      case ( '-partgraphrecursive' )
        method = PVT_METIS_PARTGRAPHRECURSIVE
        fmethod = 'metis-PartGraphRecursive'
#endif

      case ( '-w', '-weight' )
        ! force metis (the user wants weights...)
        is_metis = .true.
        ! Check for weight in the metis printing
        has_weight = .true.
        is_unit_cell = .true.

      case ( '-uc', '-unit-cell' )
        ! Check for the unit-cell region
        ! All tm are 0
        is_unit_cell = .true.
        
      case ( '-a' )
        i_arg = i_arg + 1
        call get_command_argument(i_arg, value = arg)
        read(arg,'(i16)') tm(1)
      case ( '-b' )          
        i_arg = i_arg + 1
        call get_command_argument(i_arg, value = arg)
        read(arg,'(i16)') tm(2)
      case ( '-c' )
        i_arg = i_arg + 1
        call get_command_argument(i_arg, value = arg)
        read(arg,'(i16)') tm(3)

      case ( '-h', '-help' )
        call help()
        stop

      end select

     i_arg = i_arg + 1
  end do

  ! Get the file-name for the sparsity pattern
  call get_command_argument(N_arg, value = fname)

  ! open file
  if ( .not. file_exist(fname) ) then
     print *,trim(fname)
     stop "File does not exist. Please provide file that exists..."
   else
     write(*,'(2a)')'Reading file: ', trim(fname)
  end if

  ! Read in the sparsity pattern
  call populate_Sp()

  ! First convert transfer elements
  if ( any(tm /= TM_ALL) ) then
    call crtSparsity_SC(sp_full,sp_uc, TM = tm, &
        ucell = uc, isc_off = isc_off)
    sp_full = sp_uc
    deallocate(isc_off)
  end if

  if ( is_unit_cell ) then
    call delete(sp_uc)
    ! Convert to UC
    call crtSparsity_SC(sp_full,sp_uc, UC = .true. )
  end if
  call delete(sp_full)

  if ( is_atom ) then
     call newDistribution(no_u,-1,fdit, name='fake dist')
     call SpOrb_to_SpAtom(fdit,sp_uc,na_u,lasto,sp_full)
     deallocate(lasto)
     call delete(fdit)
     sp_uc = sp_full
     call delete(sp_full)
  end if

  ! Convert sparsity pattern to graph
  call attach(sp_uc, nrows_g = no_u, &
       n_col = ncol, list_ptr = l_ptr, &
       list_col = l_col, nnzs = n_nzs )

  ! If the method is existing, do the pivoting
  if ( method >= 0 ) then
     if ( method > 0 ) then
        call sp_pvt(no_u,sp_uc,pvt,method)
     else

        ! Method is SCRAMBLE
        call rgn_range(pvt,1,no_u)
        pvt%name = 'Scramble'
        pvt%sorted = .false.

        ! Scramble 50 times the number of entries
        do i_arg = 1 , no_u * 50
           call random_number(r)
           o1 = 1 + floor(r * real(no_u))
           call random_number(r)
           o2 = 1 + floor(r * real(no_u))
           i = pvt%r(o1)
           pvt%r(o1) = pvt%r(o2)
           pvt%r(o2) = i
        end do

     end if

     ! Save the graphviz file
     fmethod = trim(fmethod) // '.gv'
     call sp2graphviz(fmethod,no_u,n_nzs,ncol,l_ptr,l_col, &
          method = graph_method, pvt = pvt )

  end if

  if ( is_metis ) call sp2metis()

  call delete(sp_uc)

contains

  subroutine sp2metis()
    integer :: io, o, o2, jo, j
    integer, pointer :: cur_col(:), col(:)
    integer, allocatable :: cw(:)
    character(len=20) :: fmt

    ! Start writing it out
    if ( has_weight ) then
       write(*,'(2(tr1,i0),tr1,a)') no_u , (n_nzs - no_u)/2,'001'
    else
       write(*,'(2(tr1,i0))') no_u , (n_nzs - no_u)/2
    end if

    ! Write out each graph point
    write(fmt,'(i0,a)') 2*no_u,'(tr1,i0)'
    do io = 1 , no_u
       
       cur_col => l_col(l_ptr(io)+1:l_ptr(io)+ncol(io))

       ! allocate graph and weight
       allocate(cw((ncol(io)-1) * 2))
     
       ! Create the weights graph
       o2 = 0
       do o = 1 , ncol(io)
          
          jo = cur_col(o)
          if ( jo == io ) cycle

          ! step position
          o2 = o2 + 2

          ! copy over position
          cw(o2-1) = jo

          ! create weights
          col => l_col(l_ptr(jo)+1:l_ptr(jo)+ncol(jo))
          cw(o2) = 0
          
          do j = 1 , ncol(jo)
             if ( any(col(j) == cur_col) ) cw(o2) = cw(o2) + 1
          end do
          
       end do
       
       if ( has_weight ) then
          write(*,'('//trim(fmt)//')') cw
       else
          write(*,'('//trim(fmt)//')') pack(cur_col, cur_col /= io)
       end if
       
       deallocate(cw)
    end do

  end subroutine sp2metis

  subroutine populate_Sp()

    use m_ts_io, only : tshs_version

    integer, parameter :: iu = 212
    integer :: nl, five(5), i, nsc(3)
    logical :: is_simple

    is_simple = .true.
    nl = len_trim(fname)
    if ( nl > 5 ) then
       if ( fname(nl-3:nl) == 'TSHS' ) then
          if ( tshs_version(fname) /= 1 ) then
             stop 'We require the latest TSHS version [1]'
          end if
          is_simple = .false.
       end if
    end if

    open(iu, file = trim(fname) , status = 'old', form = 'unformatted' )

    if ( is_simple ) then
       read(iu) no_u, i
       
       call io_read_Sp(iu, no_u, sp_full, 'sp')

    else
       ! Read in the full information
       read(iu) ! version
       read(iu) five
       na_u = five(1)
       no_u = five(2)
       n_nzs = five(5)
       read(iu) nsc
       read(iu) uc
       read(iu) ! Gamma
       read(iu) ! kscell
       read(iu) ! Ef
       read(iu) ! istep
       allocate(lasto(0:na_u))
       read(iu) lasto
       call io_read_Sp(iu,no_u,sp_full, 'sp')
       do i = 1 , no_u * (five(4) + 1)
          read(iu) ! S and H
       end do
       ! Now read in the isc_off
       i = product(nsc)
       allocate(isc_off(3,i))
       read(iu) isc_off
    end if

    close(iu)
    
  end subroutine populate_Sp

  function trim_em(s) result(f)
    character(len=*), intent(in) :: s
    character(len=250) :: f
    
    if ( s(1:2) == '--' ) then
       f = s(2:)
    else
       f = s
    end if

  end function trim_em

  subroutine help()
    character(len=20), parameter :: gf = '(tr3,a,'':'',/,tr8,a)'
    character(len=*), parameter :: nf = '(tr8,a)'

    character(len=*), parameter :: fm = '(tr11,a18,": ",a)'

    write(*,'(a)') 'The following options are available for pvtsp:'
    write(*,'(a)') 
    write(*,gf) '--help|-h','show this help menu'
    write(*,gf) '--atom','pivot in the atomic sparsity pattern, instead of the orbital(only for TSHS)'
    write(*,gf) '--graphviz|--graph','make graphviz output'
    write(*,gf) '--digraph|-di','create a directed graph'
    write(*,gf) '--metis-stdout','make METIS output (on STDOUT)'
    write(*,gf) '--pvt <method>','pivot according to a specific method'
    write(*,nf) '--<method> can be one of the following:'
    write(*,fm) 'cm', 'Cuthill-Mckee'
    write(*,fm) 'rev-cm', 'reverse Cuthill-Mckee'
    write(*,fm) 'gps', 'Gibbs-Poole-Stockmeyer'
    write(*,fm) 'rev-gps', 'reverse Gibbs-Poole-Stockmeyer'
    write(*,fm) 'pcg', 'Peripheral connectivity graph'
    write(*,fm) 'rev-pcg', 'reverse Peripheral connectivity graph'
    write(*,fm) 'ggps', 'General Gibbs-Poole-Stockmeyer'
    write(*,fm) 'rev-ggps', 'reverse General Gibbs-Poole-Stockmeyer'
#ifdef SIESTA__METIS
    write(*,fm) 'nodend', 'Metis NodeND pivoting'
    write(*,fm) 'partgraphkway', 'Metis PartGraphKway pivoting'
    write(*,fm) 'partgraphrecursive', 'Metis PartGraphRecursive pivoting'
#endif
    write(*,fm) 'scramble', 'Scramble the sparsity pattern'
    write(*,'(a)')
    write(*,gf) '--unit-cell|-uc','convert to unit-cell sparsity pattern'
    write(*,gf) '--a|-a <i>','use only ith supercell as connectivity graph (in A direction)'
    write(*,gf) '--b|-b <i>','use only ith supercell as connectivity graph (in B direction)'
    write(*,gf) '--c|-c <i>','use only ith supercell as connectivity graph (in C direction)'

  end subroutine help

end program pvtsp
