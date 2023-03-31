! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---


! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication

! Creation of the hssigma files.
module m_tbt_sigma_save

  use units, only : dp

  use m_tbt_hs, only : tTSHS
  use m_tbt_save, only : tNodeE
#ifdef NCDF_4
  use m_tbt_save, only : tbt_cdf_precision
#endif
  
  implicit none

  private 

  public :: init_Sigma_options, print_Sigma_options

#ifdef NCDF_4
  logical, save :: sigma_save      = .false.
  logical, save :: sigma_mean_save = .false.
  logical, save :: sigma_parallel  = .false.
  integer, save :: cmp_lvl    = 0

  public :: open_cdf_Sigma
  public :: init_Sigma_save
  public :: state_Sigma_save
  public :: state_Sigma2mean
#endif

contains


  subroutine init_Sigma_options(save_DATA)

    use dictionary
    use fdf

    type(dictionary_t), intent(inout) :: save_DATA

#ifdef NCDF_4

    sigma_save = fdf_get('TBT.CDF.SelfEnergy.Save',.false.)
    sigma_save = fdf_get('TBT.SelfEnergy.Save',sigma_save)
    if ( sigma_save ) then
      sigma_mean_save = fdf_get('TBT.CDF.SelfEnergy.Save.Mean',.false.)
      sigma_mean_save = fdf_get('TBT.SelfEnergy.Save.Mean',sigma_mean_save)
    end if
    cmp_lvl = fdf_get('CDF.Compress',0)
    cmp_lvl = fdf_get('TBT.CDF.Compress',cmp_lvl)
    cmp_lvl = fdf_get('TBT.CDF.SelfEnergy.Compress',cmp_lvl)
    if ( cmp_lvl < 0 ) cmp_lvl = 0
    if ( cmp_lvl > 9 ) cmp_lvl = 9
#ifdef NCDF_PARALLEL
    sigma_parallel = fdf_get('TBT.CDF.MPI',.false.)
    sigma_parallel = fdf_get('TBT.CDF.SelfEnergy.MPI',sigma_parallel)
    if ( sigma_parallel ) then
       cmp_lvl = 0
    end if
#endif

    if ( sigma_save .and. fdf_get('TBT.SelfEnergy.Only',.false.) ) then
       save_DATA = save_DATA // ('Sigma-only'.kv.1)
    end if
#endif
    
  end subroutine init_Sigma_options

  subroutine print_Sigma_options( save_DATA )

    use parallel, only: IONode
    use dictionary

    type(dictionary_t), intent(inout) :: save_DATA

    character(len=*), parameter :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'
    character(len=*), parameter :: f12='(''tbt: '',a,t53,''='',tr2,i0)'
    character(len=*), parameter :: f11='(''tbt: '',a)'

    if ( .not. IONode ) return
    
#ifdef NCDF_4
    write(*,f1)'Saving downfolded self-energies',sigma_save
    if ( .not. sigma_save ) return
#ifdef NCDF_PARALLEL
    write(*,f1)'Use parallel MPI-IO for self-energy file', sigma_parallel
#endif

    write(*,f1)'Only calc downfolded self-energies', &
         ('Sigma-only'.in.save_DATA)
    if ( cmp_lvl > 0 ) then
       write(*,f12)'Compression level of TBT.SE.nc files',cmp_lvl
    else
       write(*,f11)'No compression level of TBT.SE.nc files'
    end if
    write(*,f1)'k-average downfolded self-energies',sigma_mean_save
#else
    write(*,f11)'Saving downfolded self-energies not enabled (NetCDF4)'
#endif

  end subroutine print_Sigma_options

#ifdef NCDF_4

  subroutine open_cdf_Sigma(fname, ncdf)
    use netcdf_ncdf, ncdf_parallel => parallel

#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD
#endif

    character(len=*), intent(in) :: fname
    type(hNCDF), intent(inout) :: ncdf

    if ( .not. sigma_save ) return

#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
       call ncdf_open(ncdf,fname, mode=ior(NF90_WRITE,NF90_MPIIO), &
            comm = MPI_COMM_WORLD )
    else
#endif
       call ncdf_open(ncdf,fname, mode=NF90_WRITE)
#ifdef NCDF_PARALLEL
    end if
#endif
    
  end subroutine open_cdf_Sigma

  ! Save the self-energies of the electrodes and
  subroutine init_Sigma_save(fname, TSHS, r, btd, ispin, &
      N_Elec, Elecs, raEl, roElpd, btd_El, &
      nkpt, kpt, wkpt, NE, Eta, &
      a_Dev, a_Buf)

    use parallel, only : IONode
    use m_os, only : file_exist

    use netcdf_ncdf, ncdf_parallel => parallel
    use m_timestamp, only : datestring
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Bcast, MPI_Logical, MPI_Barrier
#endif
    use m_ts_electype
    use m_region
    use dictionary

    ! The file name that we save in
    character(len=*), intent(in) :: fname
    ! The full Hamiltonian and system at present investigation.
    ! Note the H have been shifted to zero energy
    type(tTSHS), intent(in) :: TSHS
    ! The device region that we are checking
    ! This is the device regions pivot-table!
    type(tRgn), intent(in) :: r, btd
    integer, intent(in) :: ispin
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(tRgn), intent(in) :: raEl(N_Elec), roElpd(N_Elec), btd_El(N_Elec)
    integer, intent(in) :: nkpt
    real(dp), intent(in) :: kpt(3,nkpt), wkpt(nkpt)
    integer, intent(in) :: NE
    real(dp), intent(in) :: Eta

    ! Device atoms
    type(tRgn), intent(in) :: a_Dev
    ! Buffer atoms
    type(tRgn), intent(in) :: a_Buf

    type(hNCDF) :: ncdf, grp
    type(dictionary_t) :: dic
    type(tRgn) :: r_tmp

    logical :: prec_Sigma
    logical :: exist, same
    character(len=256) :: char
    real(dp) :: mem
    character(len=2) :: unit
    integer :: i, iEl, no_e
    real(dp), allocatable :: r2(:,:)
#ifdef MPI
    integer :: MPIerror
#endif

    if ( .not. sigma_save ) return

    exist = file_exist(fname, Bcast = .true. )

    mem = 0._dp

    call tbt_cdf_precision('SelfEnergy','single',prec_Sigma)

    ! in case it already exists...
    if ( exist ) then

      ! Create a dictionary to check that the sigma file is the
      ! same
      call ncdf_open(ncdf,fname)

      dic = ('no_u'.kv.TSHS%no_u) // ('na_u'.kv.TSHS%na_u) // &
          ('nkpt'.kv.nkpt ) // ('no_d'.kv.r%n) // &
          ('ne'.kv. NE ) // ('n_btd'.kv.btd%n)
      dic = dic // ('na_d'.kv. a_Dev%n)
      if ( a_Buf%n > 0 ) then
        dic = dic // ('na_b'.kv.a_Buf%n)
      end if
      call ncdf_assert(ncdf,exist,dims=dic)
      call delete(dic)
#ifdef MPI
      call MPI_Bcast(same,1,MPI_Logical,0, &
          MPI_Comm_World,MPIerror)
#endif
      if ( .not. same ) then
        call die('Dimensions in the TBT.SE.nc file does not conform &
            &to the current simulation.')
      end if

      do iEl = 1 , N_Elec
        call ncdf_open_grp(ncdf,Elecs(iEl)%name,grp)
        dic = dic // ('no_e'.kv.Elecs(iEl)%o_inD%n)
        call ncdf_assert(grp,exist,dims=dic)
        if ( .not. exist ) then
          write(*,*) 'Assertion of dimensions in file: '//trim(fname)//' failed.'

          call die('We could not assert the dimensions TBT.SE.nc file.')
        end if
      end do
      call delete(dic)

      ! Check the variables
      ! Check the variables
      dic = ('lasto'.kvp. TSHS%lasto(1:TSHS%na_u) ) // &
          ('pivot'.kvp. r%r ) // ('btd'.kvp.btd%r)
      call rgn_copy(a_Dev, r_tmp)
      call rgn_sort(r_tmp)
      dic = dic // ('a_dev'.kvp.r_tmp%r )
      dic = dic // ('xa'.kvp. TSHS%xa)
      if ( a_Buf%n > 0 )then
        dic = dic // ('a_buf'.kvp.a_Buf%r )
      end if
      call ncdf_assert(ncdf,same,vars=dic, d_EPS = 1.e-4_dp )
      call delete(dic,dealloc=.false.) ! we have them pointing...
#ifdef MPI
      call MPI_Bcast(same,1,MPI_Logical,0, &
          MPI_Comm_World,MPIerror)
#endif
      if ( .not. same ) then
        call die('pivot, lasto, xa or a_buf in the TBT.nc file does &
            &not conform to the current simulation.')
      end if
      call rgn_delete(r_tmp)

      ! Check the k-points
      allocate(r2(3,nkpt))
      do i = 1 , nkpt
        call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
      end do
      dic = ('kpt'.kvp. r2) // ('wkpt'.kvp. wkpt)
      call ncdf_assert(ncdf,same,vars=dic, d_EPS = 1.e-7_dp )
      if ( .not. same ) then
        call die('k-points or k-weights are not the same')
      end if
      call delete(dic,dealloc = .false. )
      deallocate(r2)

      call die('Currently the TBT.SE.nc file exists, &
          &we do not currently implement a continuation scheme.')

      ! We currently overwrite the Sigma-file
      if ( IONode ) then
        write(*,'(2a)')'tbt: Overwriting self-energy file: ',trim(fname)
      end if

    else

      if ( IONode ) then
        write(*,'(2a)')'tbt: Initializing self-energy file: ',trim(fname)
      end if

    end if

    ! We need to create the file
#ifdef NCDF_PARALLEL
    if ( sigma_parallel ) then
      call ncdf_create(ncdf,fname, mode=ior(NF90_NETCDF4,NF90_MPIIO), overwrite=.true., &
          comm = MPI_COMM_WORLD, &
          parallel = .true. )
    else
      call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
    end if
#else
    call ncdf_create(ncdf,fname, mode=NF90_NETCDF4, overwrite=.true.)
#endif

    ! Save the current system size
    call ncdf_def_dim(ncdf,'no_u',TSHS%no_u)
    call ncdf_def_dim(ncdf,'na_u',TSHS%na_u)
    call ncdf_def_dim(ncdf,'nkpt',nkpt) ! Even for Gamma, it makes files unified
    call ncdf_def_dim(ncdf,'xyz',3)
    call ncdf_def_dim(ncdf,'one',1)
    call ncdf_def_dim(ncdf,'na_d',a_Dev%n)
    call ncdf_def_dim(ncdf,'no_d',r%n)
    call ncdf_def_dim(ncdf,'ne',NE)
    call ncdf_def_dim(ncdf,'n_s',product(TSHS%nsc))
    call ncdf_def_dim(ncdf,'n_btd',btd%n)
    if ( a_Buf%n > 0 ) then
      call ncdf_def_dim(ncdf,'na_b',a_Buf%n) ! number of buffer-atoms
    end if

#ifdef TBT_PHONON
    dic = ('source'.kv.'PHtrans-SE')
#else
    dic = ('source'.kv.'TBtrans-SE')
#endif

    char = datestring()
    dic = dic//('date'.kv.char(1:10))
    if ( all(TSHS%nsc(:) == 1) ) then
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
    call ncdf_put_gatt(ncdf, atts = dic )
    call delete(dic)

    ! Create all the variables needed to save the states
    dic = ('info'.kv.'Last orbitals of the equivalent atom')
    call ncdf_def_var(ncdf,'lasto',NF90_INT,(/'na_u'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, TSHS%na_u)

    dic = dic//('info'.kv.'Unit cell')
    dic = dic//('unit'.kv.'Bohr')
    call ncdf_def_var(ncdf,'cell',NF90_DOUBLE,(/'xyz','xyz'/), &
        atts = dic)
    dic = dic//('info'.kv.'Atomic coordinates')
    call ncdf_def_var(ncdf,'xa',NF90_DOUBLE,(/'xyz ','na_u'/), &
        atts = dic)
    call delete(dic)
    mem = mem + calc_mem(NF90_DOUBLE, 3, TSHS%na_u + 3)

    dic = ('info'.kv.'Supercell offsets')
    call ncdf_def_var(ncdf,'isc_off',NF90_INT,(/'xyz', 'n_s'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, 3, product(TSHS%nsc))
    
    dic = dic//('info'.kv.'Number of supercells in each direction')
    call ncdf_def_var(ncdf,'nsc',NF90_INT,(/'xyz'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, 3)
    
    dic = ('info'.kv.'Device region orbital pivot table')
    call ncdf_def_var(ncdf,'pivot',NF90_INT,(/'no_d'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, r%n)

    dic = dic // ('info'.kv.'Blocks in BTD for the pivot table')
    call ncdf_def_var(ncdf,'btd',NF90_INT,(/'n_btd'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, btd%n)

    dic = dic//('info'.kv.'Index of device atoms')
    call ncdf_def_var(ncdf,'a_dev',NF90_INT,(/'na_d'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_INT, a_Dev%n)

    if ( a_Buf%n > 0 ) then
      dic = dic//('info'.kv.'Index of buffer atoms')
      call ncdf_def_var(ncdf,'a_buf',NF90_INT,(/'na_b'/), &
          atts = dic)
      mem = mem + calc_mem(NF90_INT, a_Buf%n)
    end if

    dic = dic//('info'.kv.'k point')//('unit'.kv.'b')
    call ncdf_def_var(ncdf,'kpt',NF90_DOUBLE,(/'xyz ','nkpt'/), &
        atts = dic)
    call delete(dic)
    dic = dic//('info'.kv.'k point weights')
    call ncdf_def_var(ncdf,'wkpt',NF90_DOUBLE,(/'nkpt'/), &
        atts = dic)
    mem = mem + calc_mem(NF90_DOUBLE, 4, nkpt)

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
    mem = mem + calc_mem(NF90_DOUBLE, 1)

    call delete(dic)

    call ncdf_put_var(ncdf,'nsc',TSHS%nsc)
    call ncdf_put_var(ncdf,'isc_off',TSHS%isc_off)
    call ncdf_put_var(ncdf,'pivot',r%r)
    call ncdf_put_var(ncdf,'cell',TSHS%cell)
    call ncdf_put_var(ncdf,'xa',TSHS%xa)
    call ncdf_put_var(ncdf,'lasto',TSHS%lasto(1:TSHS%na_u))
    call rgn_copy(a_Dev, r_tmp)
    call rgn_sort(r_tmp)
    call ncdf_put_var(ncdf,'a_dev',r_tmp%r)
    call ncdf_put_var(ncdf,'btd',btd%r)
    call rgn_delete(r_tmp)
    if ( a_Buf%n > 0 ) then
      call ncdf_put_var(ncdf,'a_buf',a_Buf%r)
      mem = mem + calc_mem(NF90_INT, a_Buf%n)
    end if

    ! Save all k-points
    allocate(r2(3,nkpt))
    do i = 1 , nkpt
      call kpoint_convert(TSHS%cell,kpt(:,i),r2(:,i),1)
    end do
    call ncdf_put_var(ncdf,'kpt',r2)
    call ncdf_put_var(ncdf,'wkpt',wkpt)
    deallocate(r2)

    call ncdf_put_var(ncdf,'eta',Eta)

    do iEl = 1 , N_Elec

      call ncdf_def_grp(ncdf,trim(Elecs(iEl)%name),grp)

      ! Define atoms etc.
      i = TotUsedAtoms(Elecs(iEl))
      call ncdf_def_dim(grp,'na',i)

      dic = dic//('info'.kv.'Electrode atoms')
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
      
      dic = dic//('info'.kv.'Blocks in BTD downfolding for the pivot_down table')
      call ncdf_def_var(grp,'btd',NF90_INT,(/'n_btd'/), atts = dic)
      call ncdf_put_var(grp,'btd',btd_El(iEl)%r)
      mem = mem + calc_mem(NF90_INT, btd_El(iEl)%n)

      no_e = Elecs(iEl)%o_inD%n
      call ncdf_def_dim(grp,'no_e',no_e)

      dic = ('info'.kv.'Orbital pivot table for self-energy')
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

      dic = dic//('info'.kv.'Imaginary part of self-energy')
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

      dic = dic//('info'.kv.'Downfolded self-energy')
#ifdef TBT_PHONON
      dic = dic//('unit'.kv.'Ry**2')
#else
      dic = dic//('unit'.kv.'Ry')
#endif
      ! Chunking greatly reduces IO cost
      call ncdf_def_var(grp,'SelfEnergy', prec_Sigma, &
          (/'no_e','no_e','ne  ','nkpt'/), compress_lvl = cmp_lvl, &
          atts = dic , chunks = (/no_e,no_e,1,1/) )
      call delete(dic)

      if ( prec_Sigma ) then
        ! real + imag
        mem = mem + calc_mem(NF90_DOUBLE, no_e, no_e, NE, nkpt) * 2
      else
        mem = mem + calc_mem(NF90_DOUBLE, no_e, no_e, NE, nkpt)
      end if

    end do

    call ncdf_close(ncdf)

#ifdef MPI
    ! Ensure that the processors are aligned
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

    if ( IONode ) then
      call pretty_memory(mem, unit)
      write(*,'(3a,f8.3,tr1,a/)') 'tbt: Estimated file size of ', trim(fname), ':', &
          mem, unit
    end if

  contains

    pure function calc_mem(prec_nf90, n1, n2, n3, n4) result(kb)
      use precision, only: dp
      integer, intent(in) :: prec_nf90, n1
      integer, intent(in), optional :: n2, n3, n4
      real(dp) :: kb

      kb = real(n1, dp) / 1024._dp
      if ( present(n2) ) kb = kb * real(n2, dp)
      if ( present(n3) ) kb = kb * real(n3, dp)
      if ( present(n4) ) kb = kb * real(n4, dp)

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

  end subroutine init_Sigma_save

  subroutine state_Sigma_save(ncdf, ikpt, nE, N_Elec, Elecs, nzwork,zwork)

    use parallel, only : Node, Nodes

    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
#endif
    use m_ts_electype

    ! The file name we save too
    type(hNCDF), intent(inout) :: ncdf
    integer, intent(in) :: ikpt
    type(tNodeE), intent(in) :: nE
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    integer, intent(in) :: nzwork
    complex(dp), intent(inout), target :: zwork(nzwork)

    type(hNCDF) :: grp
    integer :: iEl, iN, no_e, n_e
    complex(dp), pointer :: Sigma2D(:,:)
#ifdef MPI
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif

    if ( .not. sigma_save ) return

    ! Save the energy-point
    if ( parallel_io(ncdf) ) then
       if ( nE%iE(Node) > 0 ) then
          call ncdf_put_var(ncdf,'E',nE%E(Node),start = (/nE%iE(Node)/) )
       end if
    else
       do iN = 0 , Nodes - 1
          if ( nE%iE(iN) <= 0 ) cycle
          call ncdf_put_var(ncdf,'E',nE%E(iN),start = (/nE%iE(iN)/) )
       end do
    end if

#ifdef MPI
    if ( .not. sigma_parallel .and. Nodes > 1 ) then
       no_e = 0
       do iEl = 1 , N_Elec
          no_e = max(no_e,Elecs(iEl)%o_inD%n)
       end do
       n_e = no_e ** 2
       if ( n_e > nzwork ) then
          call die('Could not re-use the work array for Sigma &
               &communication.')
       end if
    end if
#endif

    do iEl = 1 , N_Elec
       
       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       no_e = Elecs(iEl)%o_inD%n

       ! Create new pointer to make the below things much easier
       call pass2pnt(no_e, Elecs(iEl)%Sigma, Sigma2D)

       if ( nE%iE(Node) > 0 ) then
          call ncdf_put_var(grp,'SelfEnergy', Sigma2D, &
               start = (/1,1,nE%iE(Node),ikpt/) )
       end if

#ifdef MPI
       if ( .not. sigma_parallel .and. Nodes > 1 ) then
          n_e = no_e ** 2
          if ( Node == 0 ) then
             ! Because we are using a work-array to retrieve data
             call pass2pnt(no_e, zwork, Sigma2D)
             do iN = 1 , Nodes - 1
                if ( nE%iE(iN) <= 0 ) cycle
                call MPI_Recv(Sigma2D(1,1),n_e,MPI_Double_Complex,iN,iN, &
                     Mpi_comm_world,status,MPIerror)
                call ncdf_put_var(grp,'SelfEnergy',Sigma2D, &
                     start = (/1,1,nE%iE(iN),ikpt/) )
             end do
          else if ( nE%iE(Node) > 0 ) then
             call MPI_Send(Sigma2D(1,1),n_e,MPI_Double_Complex,0,Node, &
                  Mpi_comm_world,MPIerror)
          end if
       end if
#endif

    end do

  end subroutine state_Sigma_save

  subroutine state_Sigma2mean(fname, N_Elec, Elecs)

    use parallel, only : IONode

    use dictionary
    use netcdf_ncdf, ncdf_parallel => parallel
#ifdef MPI
    use mpi_siesta, only : MPI_COMM_WORLD, MPI_Gather
    use mpi_siesta, only : MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX
    use mpi_siesta, only : MPI_Integer, MPI_STATUS_SIZE
    use mpi_siesta, only : Mpi_double_precision
    use mpi_siesta, only : MPI_Barrier
#endif
    use m_ts_electype

    ! The file name we save too
    character(len=*), intent(in) :: fname
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)

    type(dictionary_t) :: dic
    type(hNCDF) :: ncdf, grp
    integer :: iEl, iE, ikpt
    integer :: NE, nkpt, no_e
    real(dp), allocatable :: rwkpt(:)
    complex(dp), allocatable :: c2(:,:)
    complex(dp), pointer :: Sigma(:,:)

#ifdef MPI
    integer :: MPIerror
#endif

    ! If we should not save the mean, we return immediately.
    if ( .not. sigma_save ) return
    if ( .not. sigma_mean_save ) return

    if ( .not. IONode ) then
#ifdef MPI
       call MPI_Barrier(Mpi_comm_world,MPIerror)
#endif
       return
    end if

    call timer('SE-mean', 1)

    ! We do this on one processor
    call ncdf_open(ncdf,fname, mode=NF90_WRITE)

    ! We read in the dimensions
    call ncdf_inq_dim(ncdf,'ne',len=NE)
    call ncdf_inq_dim(ncdf,'nkpt',len=nkpt)

    ! Allocate space
    allocate(rwkpt(nkpt))
    call ncdf_get_var(ncdf,'wkpt',rwkpt)

    ! When taking the mean of self-energies
    ! we need the transpose, hence we need half the
    ! contribution from Sigma and Sigma^T
    rwkpt(:) = 0.5_dp * rwkpt(:)

    ! Loop over all electrodes
    do iEl = 1 , N_Elec

      call delete(dic)

       ! We need to extend the netcdf file with the SigmaMean
       ! variable

       call ncdf_open_grp(ncdf,trim(Elecs(iEl)%name),grp)

       ! Get size of Sigma
       call ncdf_inq_dim(grp,'no_e',len=no_e)

       dic = ('info'.kv.'Downfolded self-energy, k-averaged')
       dic = dic//('unit'.kv.'Ry')
       ! Chunking greatly reduces IO cost
       call ncdf_def_var(grp,'SelfEnergyMean',NF90_DOUBLE_COMPLEX, &
            (/'no_e','no_e','ne  '/), chunks = (/no_e,no_e,1/) , &
            atts = dic ,compress_lvl = cmp_lvl )
       call delete(dic)

       ! Allocate space for the self-energy mean
       allocate(c2(no_e,no_e))

       ! Point the sigma
       ! This is a hack to ease the processing
       call pass2pnt(no_e,Elecs(iEl)%Sigma(1:no_e**2),Sigma)

       ! loop over all energy points
       do iE = 1 , NE

         ! Loop over k-points to average
         call ncdf_get_var(grp,'SelfEnergy',Sigma, &
             start=(/1,1,iE,1/) )
         
         c2(:,:) = rwkpt(1) * ( Sigma + transpose(Sigma) )

         do ikpt = 2 , nkpt
           
           ! Loop over k-points to average
           call ncdf_get_var(grp,'SelfEnergy',Sigma, &
               start=(/1,1,iE,ikpt/) )
           
           c2(:,:) = c2(:,:) + rwkpt(ikpt) * ( Sigma + transpose(Sigma) )
           
         end do

         call ncdf_put_var(grp,'SelfEnergyMean',c2, start=(/1,1,iE/) )

       end do

       deallocate(c2)

    end do

    deallocate(rwkpt)
    
    call ncdf_close(ncdf)

    call timer('SE-mean', 2)
        
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

  end subroutine state_Sigma2mean

  subroutine pass2pnt(no,Sigma,new_pnt)
    integer :: no
    complex(dp), target :: Sigma(no,no)
    complex(dp), pointer :: new_pnt(:,:)
    new_pnt => Sigma(:,:)
  end subroutine pass2pnt

#endif

end module m_tbt_sigma_save

  

