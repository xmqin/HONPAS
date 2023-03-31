! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This program has been fully implemented by:
!  Nick Papior, 2015

! This program enables the conversion of a TSHS file
! from an old format to a newer format.
! The new format is easy to extend and is less memory consuming.
program tshs2tshs
  
  use units
  use parallel
  use precision, only : dp
  use m_sparse
  use m_ts_io
  use m_ts_io_version
  use geom_helper, only : ucorb
  use class_Sparsity
  use class_OrbitalDistribution
  use class_dSpData1D
  use class_dSpData2D

  implicit none

  ! strings used to print out energies...
  character(len=500) :: filein, fileout, arg
  integer :: vin, vout

  integer :: iarg, narg
  logical :: exists, force

  ! ******************* TSHS ********************
  logical :: onlyS
  logical :: Gamma, TSGamma
  real(dp) :: ucell(3,3)
  integer :: na_u, no_l, no_u, no_s, n_nzs, nspin, nsc(3)
  real(dp), pointer :: xa(:,:) ! (3,na_u)
  integer, pointer :: lasto(:) ! (0:na_u) 
  integer, pointer :: ncol(:), l_ptr(:), l_col(:) ! (no_u)
  real(dp), pointer :: xij(:,:) => null()! (3,n_nzs)
  integer, pointer :: isc_off(:,:) => null()
  type(Sparsity) :: sp
  type(dSpData2D) :: dH, dxij
  type(dSpData1D) :: dS
  real(dp), pointer :: H(:,:), S(:) !(n_nzs,nspin),(n_nzs)
  integer :: kscell(3,3)
  real(dp) :: kdispl(3)
  real(dp) :: Ef, Qtot, Temp
  integer :: istep, ia1
  ! not really part of TSHS anymore
  integer, allocatable :: iza(:)
  ! *********************************************
  type(OrbitalDistribution), pointer :: dit
  integer :: i, n_s
  integer, allocatable :: indxuo(:) ! (no_s) 
  
  force = .false.
  IONode = .true.
  Node = 0
  Nodes = 1

  ! Default version out file to be 0 (old format)
  vout = 0

  ! Here we start the routines
  filein  = 'none'
  fileout = 'none'
  narg = command_argument_count()
  iarg = 1
  do while( iarg <= narg )
     arg = ' '
     call get_command_argument(iarg,arg)
     select case ( arg )
     case ( '-v' , '--version' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,'(i10)') vout
     case ( '-f' , '--force' )
        force = .true.
     case ( '-o', '--out', '-out' )
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        fileout = arg
     case ( '-h', '--help', '-help' )
        call help
     case default
        if ( arg(1:1) == '-' ) then
           write(0,'(a)') 'Either of the two errors has been encountered for the option "'//trim(arg)//'":'

           write(0,'(a)') ' 1) The option is not recognised'
           write(0,'(a)') ' 2) Input fdf cannot start with a hyphen "-"'
           call nl(0)
           call help
        end if
        if ( filein /= 'none' ) then
           fileout = arg
        else
           filein = arg
        end if
     end select
     iarg = iarg + 1
  end do
  if ( filein == 'none' ) then
     write(0,'(a)') 'Could not find input file on the command line'
     call nl(0)
     call help
  end if
  if ( fileout == 'none' ) then
     write(0,'(a)') 'Could not find output file on the command line'
     call nl(0)
     call help
  end if

  ! check whether the file exists
  inquire(file=filein,exist=exists)
  if (.not. exists ) then
     write(0,'(a)') 'Input file does not exist...'
     call nl(0)
     call help
  end if
  inquire(file=fileout,exist=exists)
  if ( exists ) then
     write(0,'(a)') 'Out put file already exist... You cannot overwrite...'
     call nl(0)
     call help
  end if

  ! In case the file is a NetCDF file
  i = len_trim(filein)
  if ( filein(i-1:i) == 'nc' ) then
     vin = -1 ! signals NCDF
  else
     ! get the input version
     vin = TSHS_version(filein)
  end if

  select case ( vin )
  case ( -1, 0, 1 )
     ! Do nothing
  case default
     write(0,'(a)') 'Version not recognized by this program, please update...'
     write(0,'(a)') 'Versions allowed are: 0,1'
     write(0,'(a,i0,a)') 'Found version: ',vin,' in file: '//trim(filein)
     call nl(0)
     call help
  end select

  select case ( vout )
  case ( 0, 1 )
     ! Do nothing
  case default
     write(0,'(a)') 'Version not recognized by this program, please update...'
     write(0,'(a)') 'Versions allowed are: 0,1'
     write(0,'(a,i0,a)') 'Requested version: ',vin,' for file: '//trim(fileout)
     call nl(0)
     call help
  end select

  if ( vin == vout .and. .not. force ) then
     write(0,'(a)') 'There is no need to convert the file.'
     write(0,'(a)') 'The version is the same...'
     call nl(0)
     call help
  end if

  write(*,'(a)') 'Reading in '//trim(filein)

  ! Read in TSHS
  call ts_read_tshs(filein, &
       onlyS, Gamma, TSGamma, &
       ucell, nsc, na_u, no_u, nspin,  &
       kscell, kdispl, &
       xa, lasto, &
       sp, dH, dS, isc_off, &
       Ef, Qtot, Temp, &
       istep, ia1)

  ! calculate known sets
  no_l = no_u
  n_s = product(nsc)
  no_s = n_s * no_u
  n_nzs = nnzs(sp)
  dit => dist(dS) ! S always exists
  if ( .not. onlyS ) H => val(dH)
  S => val(dS)

  call attach(sp,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

  write(*,'(a)') 'Writing to '//trim(fileout)

  if ( vout == 0 ) then

     ! Create dxij
     call newdSpData2D(sp,3,dit,dxij,'xij',sparsity_dim=2)
     
     ! Create the xij array
     call offset_xij(ucell,n_s,isc_off,na_u,xa,lasto,dxij)

     allocate(iza(na_u))
     iza(:) = 0
     xij => val(dxij)
     call write_TSHS_0(fileout, &
          onlyS, Gamma, TSGamma, &
          ucell, na_u, no_l, no_u, no_s, n_nzs, nspin,  &
          kscell, kdispl, &
          xa, iza, lasto, &
          ncol, l_ptr, l_col, xij, &
          H, S, Ef, Qtot, Temp, istep, ia1)
     deallocate(iza)

     call delete(dxij)

  else if ( vout == 1 ) then

     allocate(indxuo(no_s))
     do i = 1 , no_s
        indxuo(i) = ucorb(i,no_u)
     end do

     call ts_write_tshs(fileout, onlyS, Gamma, TSGamma, &
          ucell, nsc, isc_off, na_u, no_s, nspin, &
          kscell, kdispl, &
          xa, lasto, &
          dH, dS, indxuo, &
          Ef, Qtot, Temp, istep, ia1)

     deallocate(indxuo)

  end if

  if ( .not. onlyS ) call delete(dH)
  call delete(dS)
  call delete(sp)

  deallocate(xa,lasto)

  write(*,'(a)') trim(fileout)//' written.'

contains
  
  subroutine nl(u)
    integer, intent(in), optional :: u
    if ( present(u) ) then
       write(u,*)
    else
       write(*,*)
    end if
  end subroutine nl

  subroutine help()
    write(0,'(a)') 'Helps converting an old TranSIESTA input to the new format'
    write(0,'(a)') 'Can also read SIESTA.nc files (if compiled with NCDF_4)'
    write(0,'(a)') 'Options:'
    write(0,'(a)') '  -v|--version [0,1]:'
    write(0,'(a)') '          defines the output version of the TSHS file:'
    write(0,'(a)') '            -v 0 is the old TSHS format'
    write(0,'(a)') '            -v >0 is a newer TSHS format'
    write(0,'(a)') '  -o|--out <filename>:'
    write(0,'(a)') '          The output file name (must not exist)'
    write(0,'(a)') ' -h|--help : this help'
    stop
  end subroutine help

end program tshs2tshs 
