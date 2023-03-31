! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior, 2014

! Module for saving data
module m_tbt_save

  use precision, only : dp

  implicit none

  private

  integer, save :: cmp_lvl  = 0
  logical, save :: save_parallel = .false.
  
  ! Optional directory
  character(len=128), save :: save_dir = ' '
  public :: save_dir

  public :: init_save_options, print_save_options
  public :: name_save
#ifdef NCDF_4
  interface tbt_cdf_precision
     module procedure cdf_precision_real
     module procedure cdf_precision_cmplx
  end interface tbt_cdf_precision
  public :: tbt_cdf_precision
  public :: open_cdf_save
  public :: init_cdf_save
  public :: init_cdf_E_check
  public :: cdf_get_E_idx
  public :: cdf_get_kpt_idx
  public :: cdf_save_E
  public :: state_cdf_save
  public :: state_cdf_save_Elec
  public :: state_cdf_save_sp_dev
  public :: state_cdf2ascii
  public :: local_save_DOS
#else
  public :: init_save
  public :: init_save_Elec
  public :: step_kpt_save
  public :: state_save
  public :: state_save_Elec
  public :: end_save
#endif

#ifdef MPI
  public :: save_attach_buffer
#endif

  ! Type to control the energies that is contained
  ! This lets every processor know what the other processors have
  type :: tNodeE
     integer, allocatable :: iE(:)
     real(dp), allocatable :: E(:)
  end type tNodeE
  public :: tNodeE
  public :: MPI_BcastNode
  public :: save_parallel

#ifdef MPI
  ! Local variable to provide common routines
  ! for saving DOS
  real(dp), pointer :: rbuff1d(:)
#endif
  
contains

  subroutine MPI_BcastNode(iE, cE, nE)
#ifdef MPI
    use mpi_siesta, only : MPI_AllGather, MPI_Comm_World
    use mpi_siesta, only : MPI_Integer, MPI_Double_Precision
#endif
    integer, intent(in) :: iE
    complex(dp), intent(in) :: cE
    type(tNodeE), intent(inout) :: nE

#ifdef MPI
    real(dp) :: rE
    integer :: ierr

    rE = real(cE,dp)
    call MPI_AllGather(iE, 1, MPI_Integer, &
         nE%iE(0), 1, MPI_Integer, &
         MPI_Comm_World, ierr)
    call MPI_AllGather(rE, 1, MPI_Double_Precision, &
         nE%E(0), 1, MPI_Double_Precision, &
         MPI_Comm_World, ierr)
#else
    nE%iE(0) = iE
    nE%E(0) = real(cE,dp)
#endif

  end subroutine MPI_BcastNode

#ifdef NCDF_4

  ! Opens the save file accordingly to the setup parameters
  subroutine open_cdf_save(fname,ncdf)

    use netcdf_ncdf, ncdf_parallel => parallel
    
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif

    ! The file-name to be opened
    character(len=*), intent(in) :: fname
    type(hNCDF), intent(inout) :: ncdf

#ifdef MPI
    ! Open the netcdf file
    if ( save_parallel ) then
       call ncdf_open(ncdf,fname, mode=ior(NF90_WRITE,NF90_MPIIO), &
            comm = MPI_COMM_WORLD )
    else
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
    end if
#else
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#endif

  end subroutine open_cdf_save

#endif
  
  subroutine init_save_options()
    use m_verbosity, only: verbosity
    use parallel, only : IONode
#ifdef NCDF_PARALLEL
    use parallel, only : Nodes
#endif
    use fdf
    use m_os, only : dir_exist
#ifdef MPI
    use mpi_siesta, only : MPI_Barrier, MPI_Comm_World
#endif
    integer :: ldir

    cmp_lvl = fdf_get('CDF.Compress',0)
    cmp_lvl = fdf_get('TBT.CDF.Compress',cmp_lvl)
    if ( cmp_lvl < 0 ) cmp_lvl = 0
    if ( cmp_lvl > 9 ) cmp_lvl = 9
#ifdef NCDF_PARALLEL
    save_parallel = fdf_get('CDF.MPI',.false.)
    save_parallel = fdf_get('TBT.CDF.MPI',save_parallel)
    if ( Nodes == 1 ) save_parallel = .false.
    if ( save_parallel ) then
       cmp_lvl = 0
    end if
#else
    save_parallel = .false.
#endif

    save_dir = fdf_get('TBT.Directory',' ')
    ! Correct with suffix
    ldir = len_trim(save_dir)
    if ( ldir == 0 ) then
       ! do nothing, no save directory
    else if ( ldir == 1 ) then
       ! do nothing, we need only check for '/'
    else
       if ( save_dir(ldir-1:ldir) == '/.' ) save_dir = save_dir(1:ldir-1)
    end if
    ldir = len_trim(save_dir)
    if ( ldir > 0 ) then
       if ( save_dir(ldir:ldir) /= '/' ) save_dir = trim(save_dir)//'/'
    end if

    ! First try and create the directory
    if ( .not. dir_exist(save_dir, Bcast = .true. ) ) then
       if ( IONode ) then
          if ( verbosity > 3 ) &
               write(*,'(2a)') '*** Trying to create non-existing directory: ', &
               trim(save_dir)
          ! TODO OS call
          call system('mkdir -p '//trim(save_dir))
       end if
    end if
#ifdef MPI
    ! It may exist on other nodes
    call MPI_Barrier(MPI_Comm_World,ldir)
#endif

    if ( save_parallel ) then
       if ( .not. dir_exist(save_dir, all = .true. ) ) then
          write(*,'(a)') 'tbt: Parallel IO is not allowed by your &
               &file-system. Please remove TBT.CDF.MPI from FDF'
          call die('Directory: '//trim(save_dir)//' not visible &
               &to all processors, or simply does not exist.')
       end if
    else if ( .not. dir_exist(save_dir, Bcast = .true. ) ) then
       call die('Directory: '//trim(save_dir)//' does not exist.')
    end if
    
  end subroutine init_save_options

  subroutine print_save_options()

    use parallel, only: IONode
#ifdef NCDF_4
    use netcdf_ncdf, only : NF90_FLOAT, NF90_DOUBLE
#endif

    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'
    character(len=*), parameter :: f10='(''tbt: '',a,t53,''='',tr4,a)'
    character(len=*), parameter :: f11='(''tbt: '',a)'
    character(len=*), parameter :: f12='(''tbt: '',a,t53,''='',tr2,i0)'
#ifdef NCDF_4
    integer :: prec
#endif

    if ( .not. IONode ) return

    if ( len_trim(save_dir) > 0 ) then
       write(*,f10)'Data files stored in folder',trim(save_dir)
    else
       write(*,f11)'Data files stored in current folder'
    end if

#ifdef NCDF_4
    if ( cmp_lvl > 0 ) then
       write(*,f12) 'Compression level of TBT.nc files',cmp_lvl
    else
       write(*,f11)'No compression of TBT.nc files'
    end if
    call cdf_precision_real('none', 'single', prec)
    if ( prec == NF90_FLOAT ) then
       write(*,f10) 'Default NetCDF precision','single'
    else
       write(*,f10) 'Default NetCDF precision','double'
    end if
#ifdef NCDF_PARALLEL
    write(*,f1)'Use parallel MPI-IO for NetCDF file',save_parallel
#else
    write(*,f11)'Parallel MPI-IO not possible'
#endif
#endif
    
  end subroutine print_save_options

#ifdef NCDF_4
  subroutine cdf_precision_real(name,default,prec)

    use fdf, only : fdf_get, leqi
    use parallel, only : IONode
    use netcdf_ncdf, only : NF90_FLOAT, NF90_DOUBLE

    character(len=*), intent(in) :: name, default
    integer, intent(out) :: prec

    character(len=20) :: tmp

    ! Default, unless otherwise stated
    prec = NF90_FLOAT

    tmp = fdf_get('TBT.CDF.Precision',default)
    if ( leqi(tmp,'double') ) then
       prec = NF90_DOUBLE
    else if ( leqi(tmp,'single') &
         .or. leqi(tmp,'float') ) then
       prec = NF90_FLOAT
    else if ( IONode ) then
       write(*,'(a)')'WARNING: Could not recognize TBT.CDF.Precision, &
            &must be single|float|double, will use single.'
    end if

    if ( leqi(name, 'none') ) return
    
    tmp = fdf_get('TBT.CDF.'//trim(name)//'.Precision',tmp)
    if ( leqi(tmp,'double') ) then
       prec = NF90_DOUBLE
    else if ( leqi(tmp,'single') &
         .or. leqi(tmp,'float') ) then
       prec = NF90_FLOAT
    end if
    
  end subroutine cdf_precision_real

  subroutine cdf_precision_cmplx(name,default,prec)

    use netcdf_ncdf, only : NF90_FLOAT, NF90_DOUBLE
    use netcdf_ncdf, only : NF90_FLOAT_COMPLEX, NF90_DOUBLE_COMPLEX

    character(len=*), intent(in) :: name, default
    logical, intent(out) :: prec

    integer :: iprec

    ! Retrieve the precision as if it where a real
    call tbt_cdf_precision(name,default,iprec)
    select case ( iprec )
    case ( NF90_FLOAT )
       prec = NF90_FLOAT_COMPLEX
    case ( NF90_DOUBLE )
       prec = NF90_DOUBLE_COMPLEX
    case default
       call die('Unrecognized precision for complex value')
    end select

  end subroutine cdf_precision_cmplx

  subroutine init_cdf_save(fname,TSHS,r,btd,ispin, &
      N_Elec, Elecs, raEl, roElpd, btd_El, &
      nkpt, kpt, wkpt, NE, Eta, &
      a_Dev, a_Buf, sp_dev_sc, &
      save_DATA )

    use parallel, only : Node

    use m_os, only : file_exist

    use dictionary, assign_int => assign
    use netcdf_ncdf, ncdf_parallel => parallel
    use m_ncdf_io, only : cdf_w_Sp
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Integer, MPI_Logical
    use mpi_siesta, only : MPI_Comm_Self, MPI_Barrier
#endif
    use m_tbt_hs, only : tTSHS
    use m_ts_electype
    use m_region
    use class_OrbitalDistribution
    use class_Sparsity

    ! The file-name
    character(len=*), intent(in) :: fname
    ! The full Hamiltonian and system at present investigation.
    ! Note the H have been shifted to zero energy
    type(tTSHS), intent(in) :: TSHS
    ! The device region that we are checking
    ! This is the device regions pivot-table!
    ! Btd is the blocks in the BTD
    type(tRgn), intent(in) :: r, btd
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(tRgn), intent(in) :: raEl(N_Elec), roElpd(N_Elec), btd_El(N_Elec)
    integer, intent(in) :: nkpt
    real(dp), intent(in), target :: kpt(3,nkpt), wkpt(nkpt)
    integer, intent(in) :: NE
    real(dp), intent(in) :: Eta
    type(tRgn), intent(in) :: a_Dev
    ! In case the system has some buffer atoms.
    type(tRgn), intent(in) :: a_Buf
    ! The device sparsity pattern
    type(Sparsity), intent(inout) :: sp_dev_sc
    ! Options read from tbt_options
    type(dictionary_t), intent(inout) :: save_DATA

    character(len=50) :: tmp
    type(hNCDF) :: ncdf, grp
    type(dictionary_t) :: dic
    logical :: exist, sme, isGamma
    integer :: iEl, jEl, i, nnzs_dev, N_eigen, no_e
    integer :: prec_DOS, prec_T, prec_Teig, prec_J, prec_COOP, prec_DM
    type(OrbitalDistribution) :: fdit
    real(dp) :: mem
    character(len=2) :: unit
    real(dp), allocatable :: r2(:,:)
    type(tRgn) :: a_Dev_sort, r_tmp
#ifdef TBT_PHONON
    character(len=*), parameter :: T_unit = 'g0'
    character(len=*), parameter :: COHP_unit = 'Ry'
#else
    character(len=*), parameter :: T_unit = 'G0'
    character(len=*), parameter :: COHP_unit = 'Ry/Ry'
#endif
#ifdef MPI
    integer :: MPIerror
#endif

    ! We have to sort the device atoms.
    ! We however, know that a_Buf *is* sorted.
    call rgn_copy(a_Dev, a_Dev_sort)
    call rgn_sort(a_Dev_sort)

    ! In case the user thinks the double precision is too much
    call tbt_cdf_precision('DOS','single',prec_DOS)
    call tbt_cdf_precision('T','single',prec_T)
    call tbt_cdf_precision('T.Eig','single',prec_Teig)
    call tbt_cdf_precision('Current','single',prec_J)
    call tbt_cdf_precision('COOP','single',prec_COOP)
    call tbt_cdf_precision('DM','single',prec_DM)

    isGamma = all(TSHS%nsc(:) == 1)

    mem = 0._dp

    if ( 'T-eig' .in. save_DATA ) then
       call assign_int(N_eigen,save_DATA,'T-eig')
    else
       N_eigen = 0
    end if

    ! If compiled with net-cdf we ALWAYS save in this format
    exist = file_exist(fname, Bcast = .true. )

    ! in case it already exists...
    if ( exist ) then

       ! We check the content, and if it is the same,
       ! we allow continuation.
       call ncdf_open(ncdf,fname,mode=NF90_NOWRITE)

       dic = ('no_u'.kv. TSHS%no_u) // ('na_u'.kv. TSHS%na_u ) &
            // ('no_d'.kv. r%n ) // ('nkpt'.kv. nkpt )
       dic = dic // ('na_d'.kv. a_Dev_sort%n) // ('n_btd'.kv.btd%n)
       if ( a_Buf%n > 0 ) then
          dic = dic // ('na_b'.kv. a_Buf%n)
       end if
       call ncdf_assert(ncdf,sme,dims=dic)
       call delete(dic)
#ifdef MPI
       call MPI_Bcast(sme,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. sme ) then
          call die('Dimensions in the '//trim(fname)//' file &
               &does not conform to the current simulation.')
       end if

       ! Check the variables
       dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) ) // &
            ('pivot'.kvp. r%r ) // ('btd'.kvp.btd%r)
       dic = dic // ('xa'.kvp. TSHS%xa) // ('a_dev'.kvp.a_Dev_sort%r )
       dic = dic // ('nsc'.kvp. TSHS%nsc)
       if ( a_Buf%n > 0 )then
          dic = dic // ('a_buf'.kvp.a_Buf%r )
       end if
       call ncdf_assert(ncdf,sme,vars=dic, d_EPS = 1.e-4_dp )
       call delete(dic,dealloc=.false.) ! we have them pointing...
#ifdef MPI
       call MPI_Bcast(sme,1,MPI_Logical,0, &
            MPI_Comm_World,MPIerror)
#endif
       if ( .not. sme ) then
          call die('pivot, lasto, xa or a_buf in the '//trim(fname)//' &
               &file does not conform to the current simulation.')
       end if

       ! Check the k-points
       allocate(r2(3,nkpt))
       do i = 1 , nkpt
          call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
       end do
       dic = ('kpt'.kvp.r2) // ('wkpt'.kvp. wkpt)
       call ncdf_assert(ncdf,sme,vars=dic, d_EPS = 1.e-7_dp )
       if ( .not. sme ) then
          call die('k-points or k-weights are not the same')
       end if
       call delete(dic,dealloc = .false. )
       deallocate(r2)

       call ncdf_close(ncdf)
       
       ! The file has exactly the same content..
       if ( sme ) then

          ! We just need to set the weights to ensure
          ! unity weight.
          if ( Node == 0 .and. .not. isGamma ) then
             !call ncdf_open(ncdf,fname,mode=NF90_WRITE)
             !call ncdf_put_var(ncdf,'wkpt',wkpt,start=(/1/))
             !call ncdf_close(ncdf)
             write(*,'(a)') 'tbt: Continuation run on old TBT.nc file'
             write(*,'(a)') 'tbt: *** WARNING ***'
             write(*,'(a)') 'tbt: *** You need to make sure all &
                  &energy-points are fully contained in your current &
                  &energy range ***'
          end if

          call die('Currently the '//trim(fname)//' file exists, &
               &we do not currently implement a continuation scheme.')
          
       end if

       ! We complain to the user about it and DIE
       call die('The file content in '//trim(fname)//' &
            &is not consistent with this setup. Please delete the &
            &file.')

    else
       
       if ( Node == 0 ) then
          write(*,'(2a)')'tbt: Initializing data file: ',trim(fname)
       end if

    end if

    ! We need to create the file
    call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)

    ! Save the current system size
    call ncdf_def_dim(ncdf,'no_u',TSHS%no_u)
    call ncdf_def_dim(ncdf,'na_u',TSHS%na_u)
    ! Even for Gamma, it makes files unified
    !call ncdf_def_dim(ncdf,'nkpt',NF90_UNLIMITED) ! Parallel does not work
    call ncdf_def_dim(ncdf,'nkpt',nkpt)
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)
    call ncdf_def_dim(ncdf,'na_d',a_Dev_sort%n)
    call ncdf_def_dim(ncdf,'no_d',r%n)
    !call ncdf_def_dim(ncdf,'ne',NF90_UNLIMITED) ! Parallel does not work
    call ncdf_def_dim(ncdf,'ne',NE)
    call ncdf_def_dim(ncdf,'n_s',product(TSHS%nsc))
    call ncdf_def_dim(ncdf,'n_btd',btd%n)

    ! Create eigenvalue dimension, if needed
    if ( N_eigen > 0 ) then
       call ncdf_def_dim(ncdf,'neig',N_eigen)
    end if
    if ( a_Buf%n > 0 ) then
       call ncdf_def_dim(ncdf,'na_b',a_Buf%n) ! number of buffer-atoms
    end if

#ifdef TBT_PHONON
    dic = ('source'.kv.'PHtrans')
#else
    dic = ('source'.kv.'TBtrans')
#endif

    tmp = datestring()
    dic = dic//('date'.kv.tmp(1:10))
    if ( isGamma ) then
       dic = dic // ('Gamma'.kv.'true')
    else
       dic = dic // ('Gamma'.kv.'false')
    end if
    if ( TSHS%nspin > 1 ) then
       if ( ispin == 1 ) then
          dic = dic // ('spin'.kv.'UP')
       else
          dic = dic // ('spin'.kv.'DOWN')
       end if
    end if

    ! We initialize the counter for the current reached
    ! k-point and energy-point
    !if ( isGamma ) then
    !   dic = dic // ('k_idx_cur'.kv.1)
    !   dic = dic // ('k_idx_cur'.kv.0)
    !end if
    
    call ncdf_put_gatt(ncdf, atts = dic )
    call delete(dic)

    ! Create all the variables needed to save the states
    dic = ('info'.kv.'Last orbitals of the equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, TSHS%na_u)
    
    dic = dic//('info'.kv.'Unit cell')//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, 3, 3)
    
    dic = dic//('info'.kv.'Atomic coordinates')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
         atts = dic , chunks = (/3, TSHS%na_u/) )
    mem = mem + calc_mem(NF90_DOUBLE, 3, TSHS%na_u)
    call delete(dic)

    dic = ('info'.kv.'Supercell offsets')
    call ncdf_def_var(ncdf,'isc_off',NF90_INT,(/'xyz', 'n_s'/), &
         atts = dic)
    mem = mem + calc_mem(NF90_INT, 3, product(TSHS%nsc))
    
    dic = dic // ('info'.kv.'Number of supercells in each direction')
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
         atts = dic)
    mem = mem + calc_mem(NF90_INT, 3)

    dic = dic // ('info'.kv.'Device region orbital pivot table')
    call ncdf_def_var(ncdf,'pivot',NF90_INT,(/'no_d'/), &
         atts = dic)
    mem = mem + calc_mem(NF90_INT, r%n)

    dic = dic // ('info'.kv.'Blocks in BTD for the pivot table')
    call ncdf_def_var(ncdf,'btd',NF90_INT,(/'n_btd'/), &
         atts = dic)
    mem = mem + calc_mem(NF90_INT, btd%n)

    dic = dic // ('info'.kv.'Index of device atoms')
    call ncdf_def_var(ncdf,'a_dev',NF90_INT,(/'na_d'/), &
         atts = dic)
    mem = mem + calc_mem(NF90_INT, a_Dev_sort%n)

    if ( a_Buf%n > 0 ) then
       dic = dic // ('info'.kv.'Index of buffer atoms')
       call ncdf_def_var(ncdf,'a_buf',NF90_INT,(/'na_b'/), &
            atts = dic)
       mem = mem + calc_mem(NF90_INT, a_Buf%n)
    end if

    dic = dic // ('info'.kv.'k point')//('unit'.kv.'b')
    call ncdf_def_var(ncdf,'kpt',NF90_DOUBLE,(/'xyz ','nkpt'/), &
         atts = dic)
    call delete(dic)
    dic = ('info'.kv.'k point weights')
    call ncdf_def_var(ncdf,'wkpt',NF90_DOUBLE,(/'nkpt'/), &
         atts = dic , chunks = (/1/) )
    mem = mem + calc_mem(NF90_DOUBLE, 4, nkpt) ! kpt and wkpt

#ifdef TBT_PHONON
    dic = dic//('info'.kv.'Frequency')//('unit'.kv.'Ry')
#else
    dic = dic//('info'.kv.'Energy')//('unit'.kv.'Ry')
#endif
    call ncdf_def_var(ncdf,'E',NF90_DOUBLE,(/'ne'/), atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, NE)

    dic = dic//('info'.kv.'Imaginary part for device')
#ifdef TBT_PHONON
    dic = dic//('unit'.kv.'Ry**2')
#endif
    call ncdf_def_var(ncdf,'eta',NF90_DOUBLE,(/'one'/), atts = dic)

    if ( 'DOS-Gf' .in. save_DATA ) then
       dic = dic // ('info'.kv.'Density of states')//('unit'.kv.'1/Ry')
       call ncdf_def_var(ncdf,'DOS',prec_DOS,(/'no_d','ne  ','nkpt'/), &
            atts = dic , chunks = (/r%n,1,1/) , compress_lvl = cmp_lvl )
       mem = mem + calc_mem(prec_DOS, r%n, NE, nkpt)
    end if

    ! Clean-up dictionary
    call delete(dic)

    call ncdf_put_var(ncdf,'nsc',TSHS%nsc)
    call ncdf_put_var(ncdf,'isc_off',TSHS%isc_off)
    call ncdf_put_var(ncdf,'pivot',r%r)
    call ncdf_put_var(ncdf,'cell',TSHS%cell)
    call ncdf_put_var(ncdf,'xa',TSHS%xa)
    call ncdf_put_var(ncdf,'lasto',TSHS%lasto(1:TSHS%na_u))
    call ncdf_put_var(ncdf,'a_dev',a_Dev_sort%r)
    call ncdf_put_var(ncdf,'btd',btd%r)
    if ( a_Buf%n > 0 )then
       call ncdf_put_var(ncdf,'a_buf',a_Buf%r)
    end if

    ! We are now done with a_Dev_sort
    call rgn_delete(a_Dev_sort)

    ! Save all k-points
    ! Even though they are in an unlimited dimension,
    ! we save them instantly.
    ! This ensures that a continuation can check for 
    ! the same k-points in the same go.
    allocate(r2(3,nkpt))
    do i = 1 , nkpt
       call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
    end do
    call ncdf_put_var(ncdf,'kpt',r2)
    call ncdf_put_var(ncdf,'wkpt',wkpt)
    deallocate(r2)

    call ncdf_put_var(ncdf,'eta',Eta)

    sme = 'orb-current' .in. save_DATA
    sme = sme .or. ('COOP-Gf' .in. save_DATA)
    sme = sme .or. ('COOP-A' .in. save_DATA)
    sme = sme .or. ('COHP-Gf' .in. save_DATA)
    sme = sme .or. ('COHP-A' .in. save_DATA)
    sme = sme .or. ('DM-Gf' .in. save_DATA)
    sme = sme .or. ('DM-A' .in. save_DATA)
    if ( sme ) then
       
       ! In case we need to save the device sparsity pattern
       ! Create dimensions
       nnzs_dev = nnzs(sp_dev_sc)
       call ncdf_def_dim(ncdf,'nnzs',nnzs_dev)

       call delete(dic)

       dic = ('info'.kv.'Number of non-zero elements per row')
       call ncdf_def_var(ncdf,'n_col',NF90_INT,(/'no_u'/), &
            atts=dic)
       mem = mem + calc_mem(NF90_INT, TSHS%no_u)

       dic = dic//('info'.kv. &
            'Supercell column indices in the sparse format ')
       call ncdf_def_var(ncdf,'list_col',NF90_INT,(/'nnzs'/), &
            compress_lvl=cmp_lvl,atts=dic, chunks = (/nnzs_dev/) )
       mem = mem + calc_mem(NF90_INT, nnzs_dev)

#ifdef MPI
       call newDistribution(TSHS%no_u,MPI_Comm_Self,fdit,name='TBT-fake dist')
#else
       call newDistribution(TSHS%no_u,-1           ,fdit,name='TBT-fake dist')
#endif

       call cdf_w_Sp(ncdf,fdit,sp_dev_sc)
       call delete(fdit)
       call delete(dic)

    end if

    dic = dic // ('unit'.kv.'1/Ry')
    if ( 'DM-Gf' .in. save_DATA ) then
       dic = dic // ('info'.kv.'Green function density matrix')
       call ncdf_def_var(ncdf,'DM',prec_DM,(/'nnzs','ne  ','nkpt'/), &
           atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
       mem = mem + calc_mem(prec_DM, nnzs_dev, NE, nkpt)
    end if
    if ( 'COOP-Gf' .in. save_DATA ) then
       dic = dic // ('info'.kv.'Crystal orbital overlap population')
       call ncdf_def_var(ncdf,'COOP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
            atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
       mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
    end if
    if ( 'COHP-Gf' .in. save_DATA ) then
       dic = dic // ('info'.kv.'Crystal orbital Hamilton population')//('unit'.kv.COHP_unit)
       call ncdf_def_var(ncdf,'COHP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
            atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
       mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
    end if

    
    do iEl = 1 , N_Elec

       call delete(dic)

       call ncdf_def_grp(ncdf,trim(Elecs(iEl)%name),grp)

       ! Define atoms etc.
       i = TotUsedAtoms(Elecs(iEl))
       call ncdf_def_dim(grp,'na',i)

       dic = ('info'.kv.'Electrode atoms')
       call rgn_range(r_tmp, ELecs(iEl)%idx_a, ELecs(iEl)%idx_a + i - 1)
       call ncdf_def_var(grp,'a',NF90_INT,(/'na'/), atts = dic)
       call ncdf_put_var(grp,'a',r_tmp%r)
       mem = mem + calc_mem(NF90_INT, r_tmp%n)
       call rgn_delete(r_tmp)

       call ncdf_def_dim(grp,'na_down',raEl(iEl)%n)
       dic = dic//('info'.kv.'Electrode + downfolding atoms')
       call ncdf_def_var(grp,'a_down',NF90_INT,(/'na_down'/), atts = dic)
       call ncdf_put_var(grp,'a_down',raEl(iEl)%r)
       mem = mem + calc_mem(NF90_INT, raEl(iEl)%n)

       ! Save generic information about electrode
       dic = dic//('info'.kv.'Bloch expansion')
       call ncdf_def_var(grp,'bloch',NF90_INT,(/'xyz'/), atts = dic)
       call ncdf_put_var(grp,'bloch',Elecs(iEl)%Bloch%B)
       mem = mem + calc_mem(NF90_INT, 3)

       call ncdf_def_dim(grp,'no_down',roElpd(iEl)%n)

       dic = dic//('info'.kv.'Downfolding region orbital pivot table')
       call ncdf_def_var(grp,'pivot_down',NF90_INT,(/'no_down'/), atts = dic)
       call ncdf_put_var(grp,'pivot_down',roElpd(iEl)%r)
       mem = mem + calc_mem(NF90_INT, roElpd(iEl)%n)

       call ncdf_def_dim(grp,'n_btd',btd_El(iEl)%n)

       dic = dic//('info'.kv.'Blocks in BTD downfolding for the pivot table')
       call ncdf_def_var(grp,'btd',NF90_INT,(/'n_btd'/), atts = dic)
       call ncdf_put_var(grp,'btd',btd_El(iEl)%r)
       mem = mem + calc_mem(NF90_INT, btd_El(iEl)%n)

       no_e = Elecs(iEl)%o_inD%n
       call ncdf_def_dim(grp,'no_e',no_e)

       dic = dic//('info'.kv.'Orbital pivot table for self-energy')
       call ncdf_def_var(grp,'pivot',NF90_INT,(/'no_e'/), atts = dic)
       call ncdf_put_var(grp,'pivot',Elecs(iEl)%o_inD%r)
       mem = mem + calc_mem(NF90_INT, no_e)

       dic = dic//('info'.kv.'Chemical potential')//('unit'.kv.'Ry')
       call ncdf_def_var(grp,'mu',NF90_DOUBLE,(/'one'/), atts = dic)
       call ncdf_put_var(grp,'mu',Elecs(iEl)%mu%mu)

#ifdef TBT_PHONON
       dic = dic//('info'.kv.'Phonon temperature')
#else
       dic = dic//('info'.kv.'Electronic temperature')
#endif
       call ncdf_def_var(grp,'kT',NF90_DOUBLE,(/'one'/), atts = dic)
       call ncdf_put_var(grp,'kT',Elecs(iEl)%mu%kT)

       dic = dic//('info'.kv.'Imaginary part for self-energies')
#ifdef TBT_PHONON
       dic = dic//('unit'.kv.'Ry**2')
#endif
       call ncdf_def_var(grp,'eta',NF90_DOUBLE,(/'one'/), atts = dic)
       call ncdf_put_var(grp,'eta',Elecs(iEl)%Eta)

       dic = dic//('info'.kv.'Accuracy of the self-energy')//('unit'.kv.'Ry')
      call ncdf_def_var(grp,'Accuracy',NF90_DOUBLE,(/'one'/), atts = dic)
      call ncdf_put_var(grp,'Accuracy',Elecs(iEl)%accu)
      call delete(dic)

       mem = mem + calc_mem(NF90_DOUBLE, 4) ! mu, kT, eta, Accuracy

       call delete(dic)

       if ( ('DOS-Elecs' .in. save_DATA) .and. .not. Elecs(iEl)%out_of_core ) then

          call ncdf_def_dim(grp,'no_u',Elecs(iEl)%no_u)
          call ncdf_def_dim(grp,'na_u',Elecs(iEl)%na_u)

          dic = ('info'.kv.'Last orbitals of the equivalent atom')
          call ncdf_def_var(grp,'lasto',NF90_INT,(/'na_u'/), &
               atts = dic)
          mem = mem + calc_mem(NF90_INT, Elecs(iEl)%na_u)
          
          dic = dic//('info'.kv.'Bulk transmission')//('unit'.kv.T_unit)
          call ncdf_def_var(grp,'T',prec_T,(/'ne  ','nkpt'/), &
               atts = dic)
          mem = mem + calc_mem(prec_T, NE, nkpt)

          dic = dic//('info'.kv.'Unit cell')//('unit'.kv.'Bohr')
          call ncdf_def_var(grp,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
               atts = dic)
          mem = mem + calc_mem(NF90_DOUBLE, 3, 3)
          
          dic = dic//('info'.kv.'Atomic coordinates')
          call ncdf_def_var(grp,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
               atts = dic , chunks = (/3, Elecs(iEl)%na_u/) )
          mem = mem + calc_mem(NF90_DOUBLE, 3, Elecs(iEl)%na_u)

          ! Create electrode constant variables
          call ncdf_put_var(grp,'cell',Elecs(iEl)%cell)
          call ncdf_put_var(grp,'xa',Elecs(iEl)%xa)
          call ncdf_put_var(grp,'lasto',Elecs(iEl)%lasto(1:))

          dic = dic//('info'.kv.'Bulk density of states of electrode')
          dic = dic//('unit'.kv.'1/Ry')
          call ncdf_def_var(grp,'DOS',prec_DOS,(/'no_u','ne  ','nkpt'/), &
               atts = dic, chunks = (/Elecs(iEl)%no_u,1,1/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_DOS, Elecs(iEl)%no_u, NE, nkpt)

       end if
       call delete(dic)
       
       ! Now we will only add information that is calculated
       if ( iEl == N_Elec ) then
         ! check if all are calculated
         if ( ('DOS-A-all' .nin. save_DATA) .and. &
             ('T-all'.nin. save_DATA) ) cycle
       end if

       ! All have this unit
       dic = ('unit'.kv.'1/Ry')

       if ( 'DM-A' .in. save_DATA ) then
          dic = dic//('info'.kv.'Spectral function density matrix')
          call ncdf_def_var(grp,'DM',prec_DM,(/'nnzs','ne  ','nkpt'/), &
              atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_DM, nnzs_dev, NE, nkpt)
       end if

       if ( 'DOS-A' .in. save_DATA ) then
          dic = dic//('info'.kv.'Spectral function density of states')
          call ncdf_def_var(grp,'ADOS',prec_DOS,(/'no_d','ne  ','nkpt'/), &
               atts = dic, chunks = (/r%n,1,1/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_DOS, r%n, NE, nkpt)
       end if

       if ( 'COOP-A' .in. save_DATA ) then
          dic = dic//('info'.kv.'Crystal orbital overlap population')
          call ncdf_def_var(grp,'COOP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
               atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
       end if

       if ( 'COHP-A' .in. save_DATA ) then
          dic = dic//('info'.kv.'Crystal orbital Hamilton population')//('unit'.kv.COHP_unit)
          call ncdf_def_var(grp,'COHP',prec_COOP,(/'nnzs','ne  ','nkpt'/), &
               atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_COOP, nnzs_dev, NE, nkpt)
       end if

       ! All quantities here are transmissions.
       dic = dic//('unit'.kv.T_unit)

       if ( 'orb-current' .in. save_DATA ) then
          dic = dic//('info'.kv.'Orbital transmission')
          call ncdf_def_var(grp,'J',prec_J,(/'nnzs','ne  ','nkpt'/), &
               atts = dic , chunks = (/nnzs_dev/) , compress_lvl=cmp_lvl)
          mem = mem + calc_mem(prec_J, nnzs_dev, NE, nkpt)
       end if

       tmp = trim(Elecs(iEl)%name)
       do jEl = 1 , N_Elec
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-sum-out' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( iEl /= jEl ) then

             dic = dic//('info'.kv.'Transmission')
             call ncdf_def_var(grp,trim(Elecs(jEl)%name)//'.T',prec_T,(/'ne  ','nkpt'/), &
                  atts = dic )
             mem = mem + calc_mem(prec_T, NE, nkpt)

             if ( N_eigen > 0 ) then
                dic = dic//('info'.kv.'Transmission eigenvalues')
                call ncdf_def_var(grp,trim(Elecs(jEl)%name)//'.T.Eig',prec_Teig, &
                     (/'neig','ne  ','nkpt'/), &
                     atts = dic )
                mem = mem + calc_mem(prec_Teig, N_eigen, NE, nkpt)

             end if
             
          else

             ! For the same electrode we retain the group
             ! and utilise this for saving the reflection.
             dic = dic//('info'.kv.'Out transmission correction')
             call ncdf_def_var(grp,trim(tmp)//'.C',prec_T,(/'ne  ','nkpt'/), &
                 atts = dic )
             mem = mem + calc_mem(prec_T, NE, nkpt)

             if ( N_eigen > 0 ) then
                dic = dic//('info'.kv.'Out transmission eigenvalues')
                call ncdf_def_var(grp,trim(tmp)//'.C.Eig',prec_Teig, &
                     (/'neig','ne  ','nkpt'/), &
                     atts = dic )
                mem = mem + calc_mem(prec_Teig, N_eigen, NE, nkpt)
             end if

             dic = dic//('info'.kv.'Gf transmission')
             call ncdf_def_var(grp,trim(tmp)//'.T',prec_T,(/'ne  ','nkpt'/), &
                  atts = dic )
             mem = mem + calc_mem(prec_T, NE, nkpt)

          end if
          
       end do

    end do

    call delete(dic)

    call ncdf_close(ncdf)

#ifdef MPI
    ! Ensure that the processors are aligned
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

    if ( Node == 0 ) then
      call pretty_memory(mem, unit)
      write(*,'(3a,f8.3,tr1,a/)') 'tbt: Estimated file size of ', trim(fname), ':', &
          mem, unit
    end if

  contains

    pure function calc_mem(prec_nf90, n1, n2, n3) result(kb)
      use precision, only: dp
      integer, intent(in) :: prec_nf90, n1
      integer, intent(in), optional :: n2, n3
      real(dp) :: kb

      kb = real(n1, dp) / 1024._dp
      if ( present(n2) ) kb = kb * real(n2, dp)
      if ( present(n3) ) kb = kb * real(n3, dp)
      
      select case ( prec_nf90 )
      case ( NF90_INT, NF90_FLOAT )
        kb = kb * 4
      case ( NF90_DOUBLE )
        kb = kb * 8
      end select
      
    end function calc_mem

    pure subroutine pretty_memory(mem, unit)
      use precision, only: dp
      real(dp), intent(inout) :: mem
      character(len=2), intent(out) :: unit

      unit = 'KB'
      if ( mem > 1024._dp ) then
        mem = mem / 1024._dp
        unit = 'MB'
        if ( mem > 1024._dp ) then
          mem = mem / 1024._dp
          unit = 'GB'
          if ( mem > 1024._dp ) then
            mem = mem / 1024._dp
            unit = 'TB'
          end if
        end if
      end if

    end subroutine pretty_memory
    
  end subroutine init_cdf_save

  subroutine init_cdf_E_check(fname,E,NE)

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast
    use mpi_siesta, only : MPI_Integer
    use mpi_siesta, only : Mpi_double_precision
#endif

    character(len=*), intent(in) :: fname
    real(dp), intent(inout), allocatable :: E(:)
    integer, intent(in) :: NE

    type(hNCDF) :: ncdf
    integer :: cur_NE
#ifdef MPI
    integer :: MPIerror
#endif

    ! Open the netcdf file
    call ncdf_open(ncdf,trim(fname), mode=NF90_NOWRITE)
    
    call ncdf_inq_dim(ncdf,'ne',len=cur_NE)
    
#ifdef MPI
    call MPI_BCast(cur_NE,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

    ! Allocate all previous E
    allocate(E(NE+cur_NE))
    E(:) = huge(1._dp) ! Initialize to something no sane person would sample

    call ncdf_get_var(ncdf,'E',E,count=(/cur_NE/))
#ifdef MPI
    call MPI_Bcast(E(1),cur_NE,MPI_DOUBLE_PRECISION,0, &
         MPI_Comm_World, MPIerror)
#endif

    call ncdf_close(ncdf)
    
  end subroutine init_cdf_E_check

  subroutine cdf_get_E_idx(E,NE,zE,iE,have)
    ! the last NE energy points are, "possibly" dublicate
    real(dp), intent(in) :: E(:)
    ! The current number of searching energy-points
    integer, intent(in) :: NE
    ! The current energy-point
    complex(dp), intent(in) :: zE
    ! The current index in the local pattern
    ! If this is 0 or negative, we know we
    ! are dealing with a fake energy-point.
    ! Then we will immediately return
    integer, intent(inout) :: iE
    logical, intent(out) :: have

    real(dp), parameter :: Eta = 7.349806700083788e-06_dp ! 0.0001 eV
    real(dp) :: rE
    integer :: i

    have = .false.

    ! It is already a fake energy point
    if ( iE <= 0 ) return

    rE = real(zE,dp)
    i = size(E)
    do iE = 1 , i
       if ( abs(rE - E(iE)) < Eta ) then
          ! We have a match
          ! if the index is not in the end NE points
          ! then we know that it is already in the
          ! file.
          have = iE <= i - NE
          return
       end if
    end do

  end subroutine cdf_get_E_idx

  subroutine cdf_get_kpt_idx(fname,bkpt,ikpt)

    use parallel, only : Node
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast
    use mpi_siesta, only : MPI_Integer
#endif

    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: bkpt(3)
    integer, intent(inout) :: ikpt

    type(hNCDF) :: ncdf
    real(dp), allocatable :: rkpt(:,:)
    ! This is limiting us at 1000000 k-points, in each direction
    ! Sigh...
    real(dp), parameter :: Eta = 1.e-6_dp 
    integer :: nk, ik
#ifdef MPI
    integer :: MPIerror
#endif

    ! First we check that the kpoint exists
    ! If it does not, then definetely we do not have the 
    ! energy point.

    if ( Node == 0 ) then
       
       call ncdf_open(ncdf,trim(fname), mode=NF90_WRITE)

       ! Retrieve the attributes 
       call ncdf_inq_dim(ncdf,'nkpt',len=nk)

       ! Initialize to signal that we are going to extend
       ! the position
       ikpt = nk + 1
       
       allocate(rkpt(3,nk))
       call ncdf_get_var(ncdf,'kpt',rkpt)
       
       ! We check the k-point
       ! When they are equal we have processed all energy points
       ! currently in it.
       ! This puts a restriction on the continuation
       ! process, we only allow increasing the energy density
       ! at will, increasing the k-points at will is
       ! only allowed by using the 
       ! user input k-point file and have the SAME k-points
       ! in the beginning (the USER HAS TO DO THIS!!!!!)

       do ik = 1 , nk
          if ( all(abs(rkpt(:,ik) - bkpt(:)) < Eta) ) then
             ikpt = ik
             exit
          end if
       end do
       
       if ( ikpt == nk + 1 ) then

          ! We need to add the k-point to the list
          call ncdf_put_var(ncdf,'kpt',bkpt,start=(/1,nk+1/))

       else

          ! we signal that the k-point 
          ! is already present in the save-file.
          ikpt = - ikpt

       end if

       call ncdf_close(ncdf)

    end if

#ifdef MPI
    call MPI_Bcast(ikpt,1,MPI_Integer,0,MPI_Comm_World,MPIerror)
#endif

  end subroutine cdf_get_kpt_idx

  subroutine cdf_save_E(ncdf,nE)
    use parallel, only : Node, Nodes
    use netcdf_ncdf, ncdf_parallel => parallel

    type(hNCDF), intent(inout) :: ncdf
    type(tNodeE), intent(in) :: nE

    integer :: iN, idx

    if ( save_parallel ) then
       
       ! Create count
       idx = nE%iE(Node)
       if ( idx <= 0 ) then
          idx = 1
          iN = 0
       else
          iN = 1
       end if
       
       call ncdf_put_var(ncdf,'E',nE%E(Node), &
            start = (/idx/), count = (/iN/) )
       
    else

       ! We save the energy
#ifdef MPI
       do iN = 0 , Nodes - 1
          if ( nE%iE(iN) <= 0 ) cycle
          call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
       end do
#else
       call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
#endif
       
    end if

  end subroutine cdf_save_E

  subroutine state_cdf_save(ncdf, ikpt, nE, N_Elec, Elecs, DOS, T, &
       N_eigen, Teig, save_DATA)
    
    use parallel, only : Node, Nodes

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in) :: DOS(:,:)
    real(dp), intent(in) :: T(N_Elec+1,N_Elec)
    integer, intent(in) :: N_eigen
    real(dp), intent(in) :: Teig(N_eigen,N_Elec,N_Elec)
    type(dictionary_t), intent(in) :: save_DATA

    type(hNCDF) :: grp
    integer :: iEl, jEl, NDOS, iN, idx(2), cnt(2)
    character(len=30) :: tmp, tmp2
#ifdef MPI
    integer :: NT
    real(dp), allocatable :: thisDOS(:)
    real(dp), allocatable :: rT(:,:,:)
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-DOS-T',1)
#endif

    NDOS = size(DOS,dim=1)
           
#ifdef MPI
    if ( .not. save_parallel ) then
       if ( N_eigen > NDOS ) then
          allocate(thisDOS(N_eigen))
       else
          allocate(thisDOS(NDOS))
       end if
       call save_attach_buffer(thisDOS)
    end if
    NT = ( N_Elec + 1 ) * N_Elec
#endif

    if ( save_parallel ) then

       idx(1) = nE%iE(Node)
       cnt(1) = 1
       if ( idx(1) <= 0 ) then
          idx(1) = 1
          cnt(1) = 0
       end if
       call ncdf_put_var(ncdf,'E',nE%E(Node),start=idx(1:1), &
            count = cnt(1:1) )

    else

       ! Save the different options given to this routine
       ! We need to save the energy
#ifdef MPI
       do iN = 0 , Nodes - 1
          if ( nE%iE(iN) <= 0 ) cycle
          call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
       end do
#else
       call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
#endif

    end if

    if ( 'DOS-Gf' .in. save_DATA ) then

       ! We save the DOS
       ! This is the DOS from the Green's function
       call local_save_DOS(ncdf,'DOS',ikpt,nE,NDOS,DOS(1:NDOS,1))

    end if

    if ( 'DOS-A' .in. save_DATA ) then

       ! We save the DOS calculated from the spectral function

       do iEl = 1 , N_Elec
          if ( iEl == N_Elec ) then
             ! check if all are calculated
             if ( ('DOS-A-all' .nin. save_DATA) .and. &
                  ('T-all'.nin. save_DATA) ) cycle
          end if
          
          call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)
          
          call local_save_DOS(grp,'ADOS',ikpt,nE,NDOS,DOS(1:NDOS,1+iEl))

       end do
       
    end if

#ifdef MPI
    if ( Node == 0 .and. .not. save_parallel ) then
       allocate(rT(N_Elec+1,N_Elec,Nodes-1))
       do iN = 1 , Nodes - 1
          call MPI_Recv(rT(1,1,iN),NT,Mpi_double_precision, &
               iN, iN, Mpi_comm_world,status,MPIerror)
       end do
    else if ( Node /= 0 .and. .not. save_parallel ) then
       call MPI_Send(T(1,1),NT,Mpi_double_precision, &
            0, Node, Mpi_comm_world,MPIerror)
    end if
#endif

    ! Save transmission function
    do iEl = 1 , N_Elec
       if ( iEl == N_Elec .and. ('T-all' .nin. save_DATA) ) cycle

       ! Open group of electrode
       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       do jEl = 1 , N_Elec
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-sum-out' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( jEl == iEl ) then
             tmp  = trim(Elecs(jEl)%name)//'.C'
             tmp2 = trim(Elecs(jEl)%name)//'.T'
          else
             tmp  = trim(Elecs(jEl)%name)//'.T'
          end if

          idx = (/nE%iE(Node),ikpt/)
          cnt(:) = 1
          if ( idx(1) <= 0 ) then
             cnt = 0
             idx = 1
          end if
          ! Save data
          call ncdf_put_var(grp,tmp,T(jEl,iEl),start = idx, &
               count = cnt )
          if ( iEl == jEl ) then
             call ncdf_put_var(grp,tmp2,T(N_Elec+1,iEl), &
                  start = idx, count = cnt )
          end if
       
#ifdef MPI
          if ( Node == 0 .and. .not. save_parallel ) then
             do iN = 1 , Nodes - 1
                if ( nE%iE(iN) > 0 ) then
                   call ncdf_put_var(grp,tmp,rT(jEl,iEl,iN), &
                        start = (/nE%iE(iN),ikpt/) )
                   if ( iEl == jEl ) then
                      call ncdf_put_var(grp,tmp2,rT(N_Elec+1,iEl,iN), &
                           start = (/nE%iE(iN),ikpt/) )
                   end if
                end if
             end do
          end if
#endif

          if ( N_eigen > 0 ) then
             call local_save_DOS(grp,trim(tmp)//'.Eig',ikpt,nE,&
                  N_eigen,Teig(:,jEl,iEl))
          end if

       end do
    end do

#ifdef MPI
    if ( allocated(thisDOS) ) deallocate(thisDOS)
    if ( allocated(rT) ) deallocate(rT)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-DOS-T',2)
#endif

  end subroutine state_cdf_save

  subroutine state_cdf_save_Elec(ncdf, ikpt, nE, N_Elec, Elecs, &
       DOS, T, &
       save_DATA)
    
    use parallel, only : Node, Nodes

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Double_Precision
#endif
    use m_ts_electype

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in) :: DOS(:,:), T(N_Elec)
    type(dictionary_t), intent(in) :: save_DATA

    type(hNCDF) :: grp
    integer :: iEl
    integer :: N, idx(2), cnt(2)
#ifdef MPI
    integer :: iN
    real(dp), allocatable, target :: thisDOS(:)
    real(dp), allocatable :: rT(:,:)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-El',1)
#endif

    if ( 'DOS-Elecs' .in. save_DATA ) then

#ifdef MPI
       if ( .not. save_parallel ) then
          allocate(thisDOS(size(DOS,1)))
          call save_attach_buffer(thisDOS)
          allocate(rT(N_Elec,0:Nodes-1))
          call MPI_Gather(T(1),N_Elec,Mpi_Double_Precision, &
               rT(1,0),N_Elec,Mpi_Double_Precision, &
               0,MPI_COMM_WORLD,iEl)
       end if
#endif

       do iEl = 1 , N_Elec
          ! Skip electrodes which uses the out-of-core ability
          if ( Elecs(iEl)%out_of_core ) cycle
          
          call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

          idx = (/nE%iE(Node),ikpt/)
          cnt(:) = 1
          if ( idx(1) <= 0 ) then
             cnt = 0
             idx = 1
          end if

          ! DOS-elecs also calculates the transport
          call ncdf_put_var(grp,'T',T(iEl),start = idx, &
               count = cnt )
#ifdef MPI
          if ( .not. save_parallel ) then
             do iN = 1 , Nodes - 1
                if ( nE%iE(iN) <= 0 ) cycle
                call ncdf_put_var(grp,'T',rT(iEl,iN), &
                     start = (/nE%iE(iN),ikpt/) )
             end do
          end if
#endif
          
          N = Elecs(iEl)%no_u
          call local_save_DOS(grp,'DOS',ikpt,nE,N,DOS(1:N,iEl))
          
       end do

#ifdef MPI
       if ( allocated(thisDOS) ) deallocate(thisDOS)
       if ( allocated(rT) ) deallocate(rT)
#endif

    end if

#ifdef TBTRANS_TIMING
    call timer('cdf-w-El',2)
#endif

  end subroutine state_cdf_save_Elec

  subroutine local_save_DOS(grp,var,ikpt,nE,N,DOS)

    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif

    use netcdf_ncdf, ncdf_parallel => parallel

    type(hNCDF), intent(inout) :: grp
    character(len=*), intent(in) :: var
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N
    real(dp), intent(in) :: DOS(N)

    integer :: iN, cnt(3), idx(3)
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    if ( save_parallel ) then

       idx = (/1,nE%iE(Node),ikpt/)
       cnt(1) = N
       cnt(2) = 1
       cnt(3) = 1
       if ( idx(2) <= 0 ) then
          cnt = 0
          idx = 1
       end if
       call ncdf_put_var(grp,var,DOS,start = idx, &
            count = cnt )

       return

    end if

    if ( nE%iE(Node) > 0 ) then
       call ncdf_put_var(grp,var,DOS,start = (/1,nE%iE(Node),ikpt/) )
    end if
    
#ifdef MPI
    if ( .not. save_parallel ) then
       if ( Node == 0 ) then
          do iN = 1 , Nodes - 1
             if ( nE%iE(iN) <= 0 ) cycle
             call MPI_Recv(rbuff1d(1),N,MPI_double_precision,iN,iN, &
                  Mpi_comm_world,status,MPIerror)
             call ncdf_put_var(grp,var,rbuff1d(1:N),start = (/1,nE%iE(iN),ikpt/) )
          end do
       else if ( nE%iE(Node) > 0 ) then
          call MPI_Send(DOS(1),N,MPI_double_precision,0,Node, &
               Mpi_comm_world,MPIerror)
       end if
    end if
#endif
    
  end subroutine local_save_DOS
  
  subroutine state_cdf_save_sp_dev(ncdf, ikpt, nE, var_name, dat, El)
    
    use parallel, only : Node, Nodes
    use class_dSpData1D

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
    use mpi_siesta, only : MPI_Send, MPI_Recv
    use mpi_siesta, only : MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    character(len=*), intent(in) :: var_name
    type(dSpData1D), intent(inout) :: dat
    type(Elec), intent(inout), optional :: El

    type(hNCDF) :: grp
    integer :: nnzs_dev, iN, cnt(3), idx(3)
    real(dp), pointer :: D(:)
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-sp-dev',1)
#endif

    if ( present(El) ) then
       call ncdf_open_grp(ncdf,trim(El%name),grp)
    else
       ! Copy information to grp (no opening)
       grp = ncdf
    end if

    ! Get data and size of data
    D => val(dat)
    nnzs_dev = size(D)
    
    ! Save the current
    idx = (/1,nE%iE(Node),ikpt/)
    cnt(1) = nnzs_dev
    cnt(2) = 1
    cnt(3) = 1
    if ( idx(2) <= 0 ) then
       cnt = 0
       idx = 1
    end if
    call ncdf_put_var(grp,var_name,D,start = idx, count = cnt )

#ifdef MPI
    if ( .not. save_parallel ) then
       if ( Node == 0 ) then
          do iN = 1 , Nodes - 1
             if ( nE%iE(iN) > 0 ) then
                call MPI_Recv(D(1),nnzs_dev,Mpi_double_precision, &
                     iN, iN, Mpi_comm_world,status,MPIerror)
                call ncdf_put_var(grp,var_name,D,start = (/1,nE%iE(iN),ikpt/) )
             end if
          end do
       else if ( nE%iE(Node) > 0 ) then
          call MPI_Send(D(1),nnzs_dev,Mpi_double_precision, &
               0, Node, Mpi_comm_world,MPIerror)
       end if
    end if
#endif

#ifdef TBTRANS_TIMING
    call timer('cdf-w-sp-dev',2)
#endif

  end subroutine state_cdf_save_sp_dev


  ! Routine for reading in the TBT.nc file
  ! and convert it to regular transmission files.
  subroutine state_cdf2ascii(fname,nspin,ispin,N_Elec,Elecs,N_E,rW,save_DATA)

    use parallel, only : Node
    use units, only : eV
#ifdef TBT_PHONON
    use units, only : Kelvin
#endif

    use variable
    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel

    use m_interpolate, only : crt_pivot

    use m_timestamp, only : datestring
    use m_ts_electype

    character(len=*), intent(in) :: fname
    integer, intent(in) :: nspin, ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_E
    ! This quantity is the dE weight, with dE in Ry units
    real(dp), intent(in) :: rW(N_E)
    type(dictionary_t), intent(in) :: save_DATA

    character(len=256) :: ascii_file, tmp
    type(hNCDF) :: ncdf, grp
    logical :: exist
    integer :: iEl, jEl, i, N_eigen
    integer :: NE, nkpt, no_d, no_e
    real(dp), allocatable :: rkpt(:,:), rwkpt(:)
    real(dp), allocatable :: rE(:)
    real(dp), allocatable :: r2(:,:), r3(:,:,:)

    real(dp), parameter :: Coulomb = 1.6021766208e-19_dp
    real(dp), parameter :: eV2J = Coulomb
    real(dp), parameter :: h_eVs = 4.135667662e-15_dp
    real(dp), parameter :: h_eVfs = 4.135667662_dp
    real(dp), parameter :: hbar_eVs = 0.6582119514e-15_dp
    real(dp), parameter :: hbar_eVfs = 0.6582119514_dp
    real(dp) :: Current, eRy

#ifdef TBT_PHONON
    character(len=*), parameter :: T_unit = ' [g0]'
    real(dp) :: dT, kappa
#else
    character(len=*), parameter :: T_unit = ' [G0]'
    real(dp) :: Power, V, dd
#endif
    integer, allocatable :: pvt(:)

    ! In case we are doing something parallel, 
    ! we simply read in and write them in text based formats
    if ( Node /= 0 ) return

    call timer('cdf2ascii',1)

    tmp = datestring()
    tmp = tmp(1:10)

    ! Open the netcdf file
    call ncdf_open(ncdf,fname, mode=NF90_NOWRITE)

    ! First we read in all dimensions
    call ncdf_inq_dim(ncdf,'ne',len=NE)
    if ( NE /= N_E ) call die('Error when re-reading the number of &
        &energy-points')
    call ncdf_inq_dim(ncdf,'nkpt',len=nkpt)
    call ncdf_inq_dim(ncdf,'no_d',len=no_d)
    call ncdf_inq_dim(ncdf,'neig',exist=exist,len=N_eigen)
    if ( .not. exist ) then
      N_eigen = 0
    end if

    ! Allocate space
    allocate(rE(NE),pvt(NE))
    allocate(rkpt(3,nkpt),rwkpt(nkpt))
    ! Nearly all quantities are like this...
    allocate(r2(NE,nkpt))

    ! Read in common information
    call ncdf_get_var(ncdf,'E',rE)
    ! Convert energy to eV
    rE = rE / eV
    ! Create pivot table
    call crt_pivot(NE,rE,pvt)

    call ncdf_get_var(ncdf,'kpt',rkpt)
    call ncdf_get_var(ncdf,'wkpt',rwkpt)

    if ( 'DOS-Gf' .in. save_DATA ) then

      ! Get (orbital summed) DOS (in /eV)
      call get_DOS(ncdf, 'DOS', no_d, NE, nkpt, r2)

      if ( nkpt > 1 ) then
        call name_save(ispin,nspin,ascii_file,end='DOS')
        call save_DAT(ascii_file,nkpt,rkpt,rwkpt,no_d,NE,rE,pvt,r2,'DOS [1/eV]', &
            '# DOS calculated from the Green function, k-resolved')
      end if
      call name_save(ispin,nspin,ascii_file,end='AVDOS')
      call save_DAT(ascii_file,1,rkpt,rwkpt,no_d,NE,rE,pvt,r2,'DOS [1/eV]', &
          '# DOS calculated from the Green function, k-averaged')

    end if

#ifdef TBT_PHONON
    if ( Node == 0 .and. N_Elec > 1 ) then
      write(*,'(/,a)')'Heatflow (ensure frequency range covers temperature tails):'
    end if
#else
    if ( Node == 0 .and. N_Elec > 1 ) then
      write(*,'(/,a)')'Currents (ensure entire Fermi function window):'
    end if
#endif

    ! We should now be able to create all the files
    do iEl = 1 , N_Elec

      ! Always open the group...
      call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

      if ( ('DOS-Elecs' .in. save_DATA) .and. .not. Elecs(iEl)%out_of_core ) then

        ! Get bulk-transmission
        call ncdf_get_var(grp,'T',r2)
        if ( nkpt > 1 ) then
          call name_save(ispin,nspin,ascii_file,end='BTRANS',El1=Elecs(iEl))
          call save_DAT(ascii_file,nkpt,rkpt,rwkpt,0,NE,rE,pvt,r2,'T'//T_unit,&
              '# Bulk transmission, k-resolved')
        end if
        call name_save(ispin,nspin,ascii_file,end='AVBTRANS',El1=Elecs(iEl))
        call save_DAT(ascii_file,1,rkpt,rwkpt,0,NE,rE,pvt,r2,'T'//T_unit, &
            '# Bulk transmission, k-averaged')

        ! Bulk DOS
        call get_DOS(grp, 'DOS', Elecs(iEl)%no_u, NE, nkpt, r2)
        if ( nkpt > 1 ) then
          call name_save(ispin,nspin,ascii_file,end='BDOS',El1=Elecs(iEl))
          call save_DAT(ascii_file,nkpt,rkpt,rwkpt, &
              Elecs(iEl)%no_u,NE,rE,pvt,r2,'DOS [1/eV]',&
              '# Bulk DOS, k-resolved')
        end if
        call name_save(ispin,nspin,ascii_file,end='AVBDOS',El1=Elecs(iEl))
        call save_DAT(ascii_file,1,rkpt,rwkpt, &
            Elecs(iEl)%no_u,NE,rE,pvt,r2,'DOS [1/eV]',&
            '# Bulk DOS, k-averaged')

      end if

      ! We do not calculate the last electrode
      ! unless requested
      if ( iEl == N_Elec ) then
        ! check if all are calculated
        if ( ('DOS-A-all' .nin. save_DATA) .and. &
            ('T-all'.nin. save_DATA) ) cycle
      end if

      if ( 'DOS-A' .in. save_DATA ) then

        ! Spectral DOS
        call get_DOS(grp, 'ADOS', no_d, NE, nkpt, r2)
        if ( nkpt > 1 ) then
          call name_save(ispin,nspin,ascii_file,end='ADOS',El1=Elecs(iEl))
          call save_DAT(ascii_file,nkpt,rkpt,rwkpt,no_d,NE,rE,pvt,r2,'DOS [1/eV]',&
              '# DOS calculated from the spectral function, k-resolved')
        end if
        call name_save(ispin,nspin,ascii_file,end='AVADOS',El1=Elecs(iEl))
        call save_DAT(ascii_file,1,rkpt,rwkpt,no_d,NE,rE,pvt,r2,'DOS [1/eV]', &
            '# DOS calculated from the spectral function, k-averaged')

      end if

      if ( N_eigen > 0 ) allocate(r3(N_eigen,NE,nkpt))

      do jEl = 1 , N_Elec
        ! Calculating iEl -> jEl is the
        ! same as calculating jEl -> iEl, hence if we
        ! do not wish to assert this is true, we do not
        ! calculate this.
        if ( ('T-all' .nin. save_DATA ) .and. &
            jEl < iEl ) cycle
        if ( ('T-sum-out' .nin. save_DATA ) .and. &
            iEl == jEl ) cycle

        if ( iEl == jEl ) then
          call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.C',r2)
          ! Save the variable to ensure the correct sum in the transmission
          if ( nkpt > 1 ) then
            call name_save(ispin,nspin,ascii_file,end='CORR', El1=Elecs(iEl))
            call save_DAT(ascii_file,nkpt,rkpt,rwkpt,0,NE,rE,pvt,r2,'TC'//T_unit,&
                '# Out transmission correction, k-resolved')
          end if
          call name_save(ispin,nspin,ascii_file,end='AVCORR', El1=Elecs(iEl))
          call save_DAT(ascii_file,1,rkpt,rwkpt,0,NE,rE,pvt,r2,'TC'//T_unit,&
              '# Out transmission correction, k-averaged')

          if ( N_eigen > 0 ) then
            call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.C.Eig',r3)
            if ( nkpt > 1 ) then
              call name_save(ispin,nspin,ascii_file,end='CEIG', El1=Elecs(iEl) )
              call save_EIG(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,N_eigen,r3,'TCeig'//T_unit,&
                  '# Out transmission correction eigenvalues, k-resolved')
            end if
            call name_save(ispin,nspin,ascii_file,end='AVCEIG', El1=Elecs(iEl) )
            call save_EIG(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,N_eigen,r3,'TCeig'//T_unit,&
                '# Out transmission correction eigenvalues, k-averaged')
          end if

          ! The transmission is now the total incoming wave 
          ! [G-G^\dagger].\Gamma
          call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.T',r2)
        else
          call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.T',r2)
          if ( N_eigen > 0 ) then
            call ncdf_get_var(grp,trim(Elecs(jEl)%name)//'.T.Eig',r3)
            if ( nkpt > 1 ) then
              call name_save(ispin,nspin,ascii_file,end='TEIG', &
                  El1=Elecs(iEl), El2=Elecs(jEl))
              call save_EIG(ascii_file,nkpt,rkpt,rwkpt,NE,rE,pvt,N_eigen,r3,'Teig'//T_unit,&
                  '# Transmission eigenvalues, k-resolved')
            end if
            call name_save(ispin,nspin,ascii_file,end='AVTEIG', &
                El1=Elecs(iEl), El2=Elecs(jEl))
            call save_EIG(ascii_file,1,rkpt,rwkpt,NE,rE,pvt,N_eigen,r3,'Teig'//T_unit,&
                '# Transmission eigenvalues, k-averaged')
          end if
        end if

        ! Save transmission
        if ( nkpt > 1 ) then
          call name_save(ispin,nspin,ascii_file,end='TRANS', &
              El1=Elecs(iEl), El2=Elecs(jEl))
          call save_DAT(ascii_file,nkpt,rkpt,rwkpt,0,NE,rE,pvt,r2,'T'//T_unit,&
              '# Transmission, k-resolved')
        end if
        call name_save(ispin,nspin,ascii_file,end='AVTRANS', &
            El1=Elecs(iEl), El2=Elecs(jEl))
        call save_DAT(ascii_file,1,rkpt,rwkpt,0,NE,rE,pvt,r2,'T'//T_unit,&
            '# Transmission, k-averaged')

        ! The array r2 now contains the k-averaged transmission.
#ifdef TBT_PHONON
        ! Now we calculate the thermal current
        ! nb function is: nb(E-E1) - nb(E-E2) IMPORTANT
        !   \int d\omega 1/(2\pi) \hbar \omega [nb_1-nb_2] T(\omega)
        ! ==\int d\omega h \omega [nb_1-nb_2] T(\omega)
        Current = 0._dp
        ! Thermal conductance:
        !   \kappa = hbar ^2 / (2\pi kB T ^2)
        !            \int d \omega  \omega^2 e^{kB T \omega / \hbar}
        !              / ( e^{kB T \omega / \hbar} - 1 )^2
        kappa = 0._dp
!$OMP parallel do default(shared), private(i,eRy), &
!$OMP&  reduction(+:Current,kappa)
        do i = 1 , NE
          ! We have rE in eV, hence the conversion
          eRy = rE(i) * eV
          Current = Current + r2(i,1) * rW(i) * rE(i) * nb2(eRy, &
              Elecs(iEl)%mu%mu, Elecs(iEl)%mu%kT, &
              Elecs(jEl)%mu%mu, Elecs(jEl)%mu%kT )
          kappa = kappa + r2(i,1) * rW(i) * rE(i) ** 2 * &
              dnb(eRy, Elecs(iEl)%mu%mu, Elecs(iEl)%mu%kT)
        end do
!$OMP end parallel do

        ! 'Current' is now in [Ry] [eV]

        ! rE is already in eV, r2 and nb are unit-less.
        ! rW is in Ry => / eV
        Current = Current / eV
        ! Change from \omega -> \hbar\omega yields:
        !   \hbar d\omega -> d (\hbar\omega) = d E
        ! Hence we need to divide by:
        !   \hbar == 2pi h in [eV fs]
        Current = Current / h_eVfs

        ! 'kappa' is now in [Ry] [eV]**2
        ! (also the conversion from \omega->\hbar\omega introduces division by hbar
        !    hbar = 2\pi h
        kappa = kappa / eV / h_eVfs
        ! 'kappa' in [eV] ** 2 / [fs]
        kappa = kappa / (Elecs(iEl)%mu%kT/eV * Elecs(iEl)%mu%kT/Kelvin)
        dT = ( Elecs(iEl)%mu%kT - Elecs(jEl)%mu%kT ) / Kelvin

        if ( Node == 0 ) then
          write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
              ' -> ',trim(Elecs(jEl)%name),', dT [K] / J [eV/fs]: ', &
              dT, ' K / ',Current,' eV/fs'
          if ( abs(dT) < 10._dp ) then
            ! Only calculate kappa for small dT < 10 Kelvin
            ! Possibly we should calculate kappa(T)
            ! and save to a file.
            ! We advocate this to be calculated off-site.
            write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
                ' -> ',trim(Elecs(jEl)%name),', T+dT [K] / kappa [eV/(fs K)]: ', &
                Elecs(iEl)%mu%kT / Kelvin, ' K / ',kappa,' eV/(fs K)'
          end if
        end if
#else
        ! Now we calculate the current
        ! nf function is: nF(E-E1) - nF(E-E2) IMPORTANT
        Current = 0._dp
        Power = 0._dp
        if ( iEl == jEl ) then
          ! Do nothing
        else
!$OMP parallel do default(shared), private(i,dd,eRy), &
!$OMP&  reduction(+:Current,Power)
          do i = 1 , NE
            eRy = rE(i) * eV
            ! We have rE in eV, hence the conversion
            dd = r2(i,1) * rW(i) * nf2(eRy, &
                Elecs(iEl)%mu%mu, Elecs(iEl)%mu%kT, &
                Elecs(jEl)%mu%mu, Elecs(jEl)%mu%kT )
            Current = Current + dd
            ! rE is in eV, mu is in Ry
            Power = Power + dd * ( rE(i) - Elecs(iEl)%mu%mu / eV )
          end do
!$OMP end parallel do
        end if

        ! 'Current' is now in [Ry]

        ! rW is in Ry => / eV
        Current = Current / eV
        ! G0 = e^2 / h        ! no spin degeneracy
        !   e = 1 C = 1.602...e-19
        !   h = [eV s]
        Current = Current * Coulomb ** 2 / (h_eVs * eV2J)

        ! 'Power' is now in [Ry] [eV]

        ! rW is in Ry => / eV
        !  and apply 1 / h [eV s]
        Power = Power / eV / h_eVs
        ! Power is now in [eV] / [s]
        !   [eV] => [J] 
        Power = Power * eV2J

        ! Calculate applied bias
        V = ( Elecs(iEl)%mu%mu - Elecs(jEl)%mu%mu ) / eV

        if ( Node == 0 .and. iEl /= jEl ) then
          write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
              ' -> ',trim(Elecs(jEl)%name),', V [V] / I [A]: ', &
              V, ' V / ',Current,' A'
          write(*,'(4a,2(g12.6,a))') trim(Elecs(iEl)%name), &
              ' -> ',trim(Elecs(jEl)%name),', V [V] / P [W]: ', &
              V, ' V / ',Power,' W'
        end if
#endif

      end do

      ! Clean-up TEIG variable
      if ( allocated(r3) ) deallocate(r3)

    end do

    if ( Node == 0 ) write(*,*) ! new-line

    ! Clean-up
    deallocate(rE,rkpt,rwkpt,pvt)
    deallocate(r2)

    call ncdf_close(ncdf)

    call timer('cdf2ascii',2)

  contains

    subroutine get_DOS(grp, var, no, NE, nkpt, r2)
      type(hNCDF), intent(inout) :: grp
      character(len=*), intent(in) :: var
      integer, intent(in) :: no, NE, nkpt
      real(dp), intent(out) :: r2(NE, nkpt)

      integer :: io, ie, ik
      real(dp) :: DOS

      allocate(r3(no,NE,nkpt))

      call ncdf_get_var(grp, var, r3)

      ! Immediately sum all orbitals and convert to 1/eV
!$OMP parallel do default(shared), private(ik,ie,io,DOS)
      do ik = 1, nkpt
        do ie = 1, NE
          DOS = 0._dp
          do io = 1, no
            DOS = DOS + r3(io,ie,ik)
          end do
          r2(ie,ik) = DOS * eV
        end do
      end do
!$OMP end parallel do

      deallocate(r3)

    end subroutine get_DOS

    subroutine save_DAT(fname,nkpt,kpt,wkpt,N,NE,E,ipiv,DAT,value,header)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nkpt, NE, ipiv(NE), N
      real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt), E(NE)
      real(dp), intent(inout) :: DAT(NE,nkpt)
      character(len=*), intent(in) :: value, header

      integer :: iu, ik, i

      call io_assign(iu)
      open( iu, file=trim(fname), form='formatted', status='unknown' ) 

      write(iu,'(a)') trim(header)
      write(iu,'(a)') '# Date: '//trim(tmp)
#ifdef TBT_PHONON
      write(iu,'(a,a9,tr1,a16)')"#","Omega [eV]", value
#else
      write(iu,'(a,a9,tr1,a16)')"#","E [eV]", value
#endif

      do ik = 1 , nkpt 
        if ( nkpt > 1 ) then
          write(iu,'(/,a6,3(e16.8,'' ''),a,e15.8)') '# kb= ',kpt(:,ik) ,'w= ',wkpt(ik)
        end if
        do i = 1 , NE
#ifdef TBT_PHONON
          write(iu,'(f10.6,tr1,e16.8)') E(ipiv(i)), DAT(ipiv(i),ik)
#else
          write(iu,'(f10.5,tr1,e16.8)') E(ipiv(i)), DAT(ipiv(i),ik)
#endif
        end do
        if ( nkpt > 1 ) then
          ! Update the average values in the first entry
          if ( ik == 1 ) then
            call dscal(NE, wkpt(1), DAT(1,1), 1)
          else
            call daxpy(NE, wkpt(ik), DAT(1,ik), 1, DAT(1,1), 1)
          end if
        end if
      end do

      call io_close(iu)

    end subroutine save_DAT

    subroutine save_EIG(fname,nkpt,kpt,wkpt,NE,E,ipiv,neig,EIG,value,header)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: nkpt, NE, neig, ipiv(NE)
      real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt), E(NE)
      real(dp), intent(inout) :: EIG(neig,NE,nkpt)
      character(len=*), intent(in) :: value, header

      integer :: iu, ik, i, ie
      character(len=20) :: fmt

      ! Create format
      write(fmt,'(a,i0,a)')'(f10.5,tr1,',neig,'e16.8)'
      call io_assign(iu)
      open( iu, file=trim(fname), form='formatted', status='unknown' ) 

      write(iu,'(a)') trim(header)
      write(iu,'(a)') '# Date: '//trim(tmp)
#ifdef TBT_PHONON
      write(iu,'(a,a9,tr1,a16)')"#","Omega [eV]", value
#else
      write(iu,'(a,a9,tr1,a16)')"#","E [eV]", value
#endif
      do ik = 1 , nkpt 
        if ( nkpt > 1 ) then
          write(iu,'(/,a6,3(e16.8,'' ''),a,e15.8)') &
              '# kb= ',kpt(:,ik) ,'w= ',wkpt(ik)
        end if
        do i = 1 , NE
          write(iu,fmt) E(ipiv(i)), EIG(:,ipiv(i),ik)
        end do
        if ( nkpt > 1 ) then
          ! Update the average values in the first entry
          if ( ik == 1 ) then
            call dscal(NE*neig, wkpt(1), EIG(1,1,1), 1)
          else
            call daxpy(NE*neig, wkpt(ik), EIG(1,1,ik), 1, EIG(1,1,1), 1)
          end if
        end if
      end do

      call io_close(iu)

    end subroutine save_EIG

#ifdef TBT_PHONON
    elemental function nb2(E,E1,kT1,E2,kT2)
      real(dp), intent(in) :: E,E1,kT1,E2,kT2
      real(dp) :: nb2
      nb2 = nb(E,E1,kT1) - nb(E,E2,kT2)
    end function nb2
    elemental function nb(E,Ef,kT)
      real(dp), intent(in) :: E,Ef,kT
      real(dp) :: nb
      nb = 1._dp/(exp((E-Ef)/kT)-1._dp)
    end function nb
    elemental function dnb(E,Ef,kT)
      real(dp), intent(in) :: E,Ef,kT
      real(dp) :: dnb
      dnb = exp((E-Ef)/kT)
      dnb = dnb / ( dnb - 1 ) ** 2
    end function dnb
#else
    elemental function nf2(E,E1,kT1,E2,kT2)
      real(dp), intent(in) :: E,E1,kT1,E2,kT2
      real(dp) :: nf2
      nf2 = nf(E,E1,kT1) - nf(E,E2,kT2)
    end function nf2
    elemental function nf(E,Ef,kT)
      real(dp), intent(in) :: E,Ef,kT
      real(dp) :: nf
      nf = 1._dp/(exp((E-Ef)/kT)+1._dp)
    end function nf
#endif

  end subroutine state_cdf2ascii

#endif


#ifdef MPI
  subroutine save_attach_buffer(array)
    real(dp), intent(inout), target :: array(:)
    rbuff1d => array(:)
  end subroutine save_attach_buffer
#endif

  ! Get the file name
  subroutine name_save(ispin,nspin,fname,end,El1,El2)
    use files, only : slabel
    use m_ts_electype
    integer, intent(in) :: ispin, nspin
    character(len=*), intent(out) :: fname
    character(len=*), intent(in), optional :: end ! designator of the file
    type(Elec), intent(in), optional :: El1, El2

    fname = ' '
#ifdef TBT_PHONON
    fname = trim(save_dir)//trim(slabel)//'.PHT'
#else
    fname = trim(save_dir)//trim(slabel)//'.TBT'
#endif

    ! Now figure out the file name
    if ( nspin > 1 ) then
       if( ispin .eq. 1 ) fname = trim(fname)//"_UP"
       if( ispin .eq. 2 ) fname = trim(fname)//"_DN"
    end if

    ! Add the designator
    if ( present(end) ) then
       fname = trim(fname)//'.'//trim(end)
    end if

    if ( present(El1) ) then
       fname = trim(fname)//'_'//trim(El1%name)
       if ( present(El2) ) then
          fname = trim(fname)//'-'//trim(El2%name)
       end if
    end if

  end subroutine name_save

#ifndef NCDF_4
  ! These routines are for creating the output data in
  ! pure ASCII format.
  

  ! This routine prepares the files for saving ASCII format
  ! data.
  ! NOTE that ASCII data will only be created in case
  ! of Netcdf not being compiled in
  subroutine init_save(iounits,ispin,nspin,no_d, N_Elec, Elecs, &
       N_eigen, save_DATA)
    
    use parallel, only : Node

    use dictionary
    use m_timestamp, only : datestring
    use m_ts_electype

    integer, intent(inout) :: iounits(:)
    integer, intent(in)    :: ispin, nspin
    integer, intent(in)    :: no_d
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: N_eigen
    type(dictionary_t), intent(in) :: save_DATA

    character(len=128) :: ascii_file, tmp
    integer :: iu, cu

    integer :: iEl, jEl

    ! Only the IO-Node can prepare the output data.
    if ( Node /= 0 ) return

    tmp = datestring()
    tmp = tmp(1:10)

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call name_save(ispin,nspin,ascii_file,end='DOS')

       call io_assign(iu)
       open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
       write(iu,'(a)') '# DOS calculated from the Green function, k-resolved'
       write(iu,'(a)') '# Date: '//trim(tmp)
       write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'DOS [1/eV]'
       
       iounits(cu) = iu

       cu = cu + 1
       
    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then

          call name_save(ispin,nspin,ascii_file,end='ADOS',El1=Elecs(iEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
          write(iu,'(a)') '# DOS calculated from the spectral function, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'DOS [1/eV]'

          iounits(cu) = iu
          
          cu = cu + 1

       end if

       do jEl = 1 , N_Elec

          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-sum-out' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle


          if ( iEl == jEl ) then

             call name_save(ispin,nspin,ascii_file,end='CORR', &
                  El1=Elecs(iEl))
             
             call io_assign(iu)
             open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
             write(iu,'(a)') '# Out transmission correction, k-resolved'
             write(iu,'(a)') '# Date: '//trim(tmp)
             write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Correction'

             iounits(cu) = iu
             
             cu = cu + 1

             if ( N_eigen > 0 ) then
                call name_save(ispin,nspin,ascii_file,end='CEIG', &
                     El1=Elecs(iEl))
                
                call io_assign(iu)
                open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
                write(iu,'(a)') '# Out transmission correction eigenvalues, k-resolved'
                write(iu,'(a)') '# Date: '//trim(tmp)
                write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Eigenvalues'
                
                iounits(cu) = iu
                
                cu = cu + 1
             end if
             
          end if
             
          call name_save(ispin,nspin,ascii_file,end='TRANS', &
               El1=Elecs(iEl), El2=Elecs(jEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown')
          write(iu,'(a)') '# Transmission, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Transmission'

          iounits(cu) = iu
          
          cu = cu + 1

          if ( jEl /= iEl .and. N_eigen > 0 ) then
             call name_save(ispin,nspin,ascii_file,end='TEIG', &
                  El1=Elecs(iEl), El2=Elecs(jEl))
             
             call io_assign(iu)
             open( iu, file=trim(ascii_file), form='formatted', status='unknown')
             write(iu,'(a)') '# Transmission eigenvalues, k-resolved'
             write(iu,'(a)') '# Date: '//trim(tmp)
             write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'Eigenvalues'
             
             iounits(cu) = iu
             
             cu = cu + 1
          end if
          
       end do
       
    end do

  end subroutine init_save

  subroutine init_save_Elec(iounits,ispin,nspin,N_Elec,Elecs,save_DATA)
    
    use parallel, only : Node

    use dictionary
    use m_timestamp, only : datestring
    use m_ts_electype

    integer, intent(inout) :: iounits(:)
    integer, intent(in)    :: ispin, nspin
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(dictionary_t), intent(in) :: save_DATA

    character(len=128) :: ascii_file, tmp
    integer :: iu, cu

    integer :: iEl

    ! Only the IO-Node can prepare the output data.
    if ( Node /= 0 ) return

    tmp = datestring()
    tmp = tmp(1:10)

    cu = 1

    if ( 'DOS-Elecs' .in. save_DATA ) then

       do iEl = 1 , N_Elec
          ! Skip electrodes which uses GF files (they don't contain the information)
          if ( Elecs(iEl)%out_of_core ) cycle

          call name_save(ispin,nspin,ascii_file,end='BDOS',El1=Elecs(iEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
          write(iu,'(a)') '# Bulk DOS, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'DOS [1/eV]'

          iounits(cu) = iu
          
          cu = cu + 1

          call name_save(ispin,nspin,ascii_file,end='BTRANS',El1=Elecs(iEl))

          call io_assign(iu)
          open( iu, file=trim(ascii_file), form='formatted', status='unknown' ) 
          write(iu,'(a)') '# Bulk transmission, k-resolved'
          write(iu,'(a)') '# Date: '//trim(tmp)
          write(iu,'(a,a9,tr1,a16)')'#','E [eV]', 'T'

          iounits(cu) = iu
          
          cu = cu + 1

       end do
       
    end if
    
  end subroutine init_save_Elec

  subroutine step_kpt_save(iounits,nkpt,bkpt,wkpt)
    
    use parallel, only : Node

    integer, intent(in) :: iounits(:), nkpt
    real(dp), intent(in) :: bkpt(3), wkpt

    integer :: cu, nu
    logical :: is_open

    if ( Node /= 0 ) return
    if ( nkpt == 1 ) return

    cu = 1
    nu = size(iounits)

    do while ( cu <= nu )

       inquire(iounits(cu), opened = is_open)
       if ( .not. is_open ) exit

       call wrt_k(iounits(cu),bkpt,wkpt)
       cu = cu + 1
       
    end do

  contains

    subroutine wrt_k(iu,bkpt,wkpt)
      integer, intent(in) :: iu
      real(dp) :: bkpt(3), wkpt
      write(iu,'(/,a6,3(e16.8,'' ''),a,e15.8)') '# kb= ',bkpt(:) ,'w= ',wkpt
    end subroutine wrt_k
    
  end subroutine step_kpt_save


  subroutine end_save(iounits)
    
    use parallel, only : Node

    integer, intent(in) :: iounits(:)

    integer :: cu, nu
    logical :: is_open

    if ( Node /= 0 ) return

    cu = 1
    nu = size(iounits)

    do while ( cu <= nu )

       inquire(iounits(cu), opened = is_open)
       if ( .not. is_open ) exit

       call io_close(iounits(cu))
       cu = cu + 1
       
    end do

  end subroutine end_save

  ! This routine prepares the files for saving ASCII format
  ! data.
  ! NOTE that ASCII data will only be created in case
  ! of Netcdf not being compiled in
  subroutine state_save(iounits,no_d, nE,N_Elec,Elecs, DOS, T, &
       N_eigen, Teig, &
       save_DATA )
    
    use parallel, only : Nodes
    use units, only : eV
    use m_interpolate, only : crt_pivot

    use dictionary
    use m_ts_electype

    integer, intent(in)    :: iounits(:)
    integer, intent(in)    :: no_d
    type(tNodeE), intent(in) :: nE
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in)   :: DOS(:,:), T(:,:)
    integer, intent(in)    :: N_eigen
    real(dp), intent(in)   :: Teig(N_eigen,N_Elec,N_Elec)
    type(dictionary_t), intent(in) :: save_DATA

    integer :: cu
    integer :: iEl, jEl, N
    real(dp), allocatable, target :: thisDOS(:)
    integer, allocatable :: ipvt(:)

    allocate(ipvt(Nodes))
    ipvt = 1

    N = size(DOS,dim=1)

#ifdef MPI
    if ( Nodes > 1 ) then
       call crt_pivot(Nodes,nE%E,ipvt)
       if ( N_eigen > N ) then
          allocate(thisDOS(N_eigen))
       else
          allocate(thisDOS(N))
       end if
       call save_attach_buffer(thisDOS)
    end if
#endif

    ! Correct for rank indices
    ipvt(:) = ipvt(:) - 1

    cu = 1

    if ( 'DOS-Gf' .in. save_DATA ) then

       call local_save_DAT(iounits(cu),nE,ipvt,no_d,DOS(1:no_d,1),fact=eV)

       cu = cu + 1

    end if

    do iEl = 1 , N_Elec

       ! We do not calculate the last electrode
       ! unless requested
       if ( iEl == N_Elec ) then
          ! check if all are calculated
          if ( ('DOS-A-all' .nin. save_DATA) .and. &
               ('T-all'.nin. save_DATA) ) cycle
       end if
       
       if ( 'DOS-A' .in. save_DATA ) then

          call local_save_DAT(iounits(cu),nE,ipvt,no_d,DOS(1:no_d,1+iEl),fact=eV)
          
          cu = cu + 1
          
       end if

       do jEl = 1 , N_Elec

          ! Calculating iEl -> jEl is the
          ! same as calculating jEl -> iEl, hence if we
          ! do not wish to assert this is true, we do not
          ! calculate this.
          if ( ('T-all' .nin. save_DATA ) .and. &
               jEl < iEl ) cycle
          if ( ('T-sum-out' .nin. save_DATA ) .and. &
               iEl == jEl ) cycle

          if ( jEl == iEl ) then
             ! Note this is reversed according to the 
             ! creation of the arrays.
             ! This is because the reflection is 1->1
             ! and the transmission is the G.\Gamma
             ! flux. Hence we simply reverse the print-outs.
             call local_save_DAT(iounits(cu),nE,ipvt,1,T(N_Elec+1:N_Elec+1,iEl))
             cu = cu + 1

             if ( N_eigen > 0 ) then
                call local_save_EIG(iounits(cu),nE,ipvt,N_eigen,Teig(:,jEl,iEl))
                cu = cu + 1
             end if

          end if

          call local_save_DAT(iounits(cu),nE,ipvt,1,T(jEl:jEl,iEl))
          
          cu = cu + 1
          
          if ( jEl /= iEl .and. N_eigen > 0 ) then
             call local_save_EIG(iounits(cu),nE,ipvt,N_eigen,Teig(:,jEl,iEl))
             cu = cu + 1
          end if

       end do
       
    end do

    deallocate(ipvt)
#ifdef MPI
    if ( allocated(thisDOS) ) deallocate(thisDOS)
#endif

  end subroutine state_save

  ! This routine prepares the files for saving ASCII format
  ! data.
  ! NOTE that ASCII data will only be created in case
  ! of Netcdf not being compiled in
  subroutine state_save_Elec(iounits,nE,N_Elec,Elecs, &
       DOS, T, &
       save_DATA )
    
    use parallel, only : Nodes
    use units, only : eV
    use m_interpolate, only : crt_pivot

    use dictionary
    use m_ts_electype

    integer, intent(in)    :: iounits(:)
    type(tNodeE), intent(in) :: nE
    integer, intent(in)    :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    real(dp), intent(in)   :: DOS(:,:), T(N_Elec)
    type(dictionary_t), intent(in) :: save_DATA

    integer :: cu, N
    integer :: iEl
    integer, allocatable :: ipvt(:)
#ifdef MPI
    real(dp), allocatable, target :: thisDOS(:)
#endif
    
    allocate(ipvt(Nodes))
    ipvt = 1

#ifdef MPI
    if ( Nodes > 1 ) then
       call crt_pivot(Nodes,nE%E,ipvt)
       allocate(thisDOS(size(DOS,dim=1)))
       call save_attach_buffer(thisDOS)
    end if
#endif

    ! Correct for rank indices
    ipvt(:) = ipvt(:) - 1

    cu = 1

    if ( 'DOS-Elecs' .in. save_DATA ) then

       do iEl = 1 , N_Elec
          ! Skip electrodes which uses GF files (they don't contain the information)
          if ( Elecs(iEl)%out_of_core ) cycle

          N = Elecs(iEl)%no_u
          call local_save_DAT(iounits(cu),nE,ipvt,N,DOS(1:N,iEl),fact=eV)
          
          cu = cu + 1

          ! save bulk transmission
          call local_save_DAT(iounits(cu),nE,ipvt,1,T(iEl))

          cu = cu + 1
          
       end do

    end if

    deallocate(ipvt)
#ifdef MPI
    if ( allocated(thisDOS) ) deallocate(thisDOS)
#endif

  end subroutine state_save_Elec
  
  subroutine local_save_DAT(iu,nE,ipvt,N,DATA,fact)
    use parallel, only : Node, Nodes
    use units, only : eV

#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif

    integer, intent(in) :: iu
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: ipvt(:), N
    real(dp), intent(in) :: DATA(N)
    real(dp), intent(in), optional :: fact

    integer :: iN, i
    real(dp) :: rnd
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    if ( present(fact) ) then
      rnd = fact
    else
      rnd = 1._dp
    end if
    
    if ( Node == 0 ) then
       do iN = 0 , Nodes - 1
          i = ipvt(iN+1) ! sorting E
#ifdef MPI
          if ( nE%iE(i) <= 0 ) cycle ! if the energy point is fake, discard
#endif
          if ( i == 0 ) then ! local node
             write(iu,'(f10.5,tr1,e16.8)') nE%E(i) / eV, sum(DATA(1:N)) * rnd
          else
#ifdef MPI
             call MPI_Recv(rbuff1d(1),1,MPI_double_precision,i,i, &
                  Mpi_comm_world,status,MPIerror)
             write(iu,'(f10.5,tr1,e16.8)') nE%E(i) / eV, rbuff1d(1)
#else
             call die('Error')
#endif
          end if
       end do
    else if ( nE%iE(Node) > 0 ) then
#ifdef MPI
      rnd = sum(DATA(1:N)) * rnd
      call MPI_Send(rnd,1,MPI_double_precision,0,Node, &
          Mpi_comm_world,MPIerror)
#endif
    end if

  end subroutine local_save_DAT

  subroutine local_save_EIG(iu,nE,ipvt,N,EIG)
    use parallel, only : Node, Nodes
    use units, only : eV

#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif

    integer, intent(in) :: iu
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: ipvt(:), N
    real(dp), intent(in) :: EIG(N)

    integer :: iN, i
    character(len=20) :: fmt
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    write(fmt,'(a,i0,a)')'(f10.5,tr1,',N,'e16.8)'

    if ( Node == 0 ) then
       do iN = 0 , Nodes - 1
          i = ipvt(iN+1) ! sorting E
#ifdef MPI
          if ( nE%iE(i) <= 0 ) cycle ! if the energy point is fake, discard
#endif
          if ( i == 0 ) then ! local node
             write(iu,fmt) nE%E(i) / eV, EIG(:)
          else
#ifdef MPI
             call MPI_Recv(rbuff1d(1),N,MPI_double_precision,i,i, &
                  Mpi_comm_world,status,MPIerror)
             write(iu,fmt) nE%E(i) / eV, rbuff1d(1:N)
#else
             call die('Error')
#endif
          end if
       end do
    else if ( nE%iE(Node) > 0 ) then
#ifdef MPI
       call MPI_Send(EIG(1),N,MPI_double_precision,0,Node, &
            Mpi_comm_world,MPIerror)
#endif
    end if

  end subroutine local_save_EIG

#endif

end module m_tbt_save
