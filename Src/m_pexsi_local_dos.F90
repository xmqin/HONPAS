
! Tangled code
  MODULE m_pexsi_local_dos
#ifdef SIESTA__PEXSI
  private
  public :: pexsi_local_dos

  CONTAINS
  subroutine pexsi_local_dos( )
    use m_energies
  
    use sparse_matrices
    USE siesta_options
    use siesta_geom
    use atomlist,       only: indxuo, indxua           
    use atomlist,       only: qtot, no_u, no_l
    use atomlist,       only: iphorb                   
    use atomlist,       only: datm, no_s, iaorb        
    use fdf
    use files,          only : slabel     
    use files,          only : filesOut_t ! type for output file names
    use parallel,       only:  SIESTA_worker
    use m_ntm
    use m_forces,       only: fa
    use m_spin,         only: nspin
    use atomlist,       only: qs => qtots
    use m_energies,     only: efs
    use m_dhscf,        only: dhscf
  
    implicit none
  
    integer   :: dummy_iscf = 1
    real(dp)  :: dummy_str(3,3), dummy_strl(3,3)  ! for dhscf call
    real(dp)  :: dummy_dipol(3)
  
    real(dp)  :: factor, g2max
    real(dp)  :: energy, broadening
  
    type(filesOut_t)     :: filesOut  ! blank output file names
  
    energy = fdf_get('PEXSI.LDOS.Energy',0.0_dp,"Ry")
    broadening = fdf_get('PEXSI.LDOS.Broadening',0.01_dp,"Ry")
  
    ! Note that we re-use Dscf, so it will be obliterated
    call get_LDOS_SI( no_u, no_l, nspin,  &
         maxnh, numh, listh, H, S,  &
         Dscf, energy, broadening)
  
    if (SIESTA_worker) then
       !Find the LDOS in the real space mesh
       filesOut%rho = trim(slabel) // '.LDSI'
       g2max = g2cut
  
       ! There is too much clutter here, because dhscf is
       ! a "kitchen-sink" routine that does too many things
  
       call dhscf( nspin, no_s, iaorb, iphorb, no_l, &
            no_u, na_u, na_s, isa, xa_last, indxua,  &
            ntm, 0, 0, 0, filesOut,                  &
            maxnh, numh, listhptr, listh, Dscf, Datm, maxnh, H, &
            Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc, &
            dummy_dipol, dummy_str, fa, dummy_strl )
    endif
  
  END subroutine pexsi_local_dos
  subroutine get_LDOS_SI( no_u, no_l, nspin_in,  &
       maxnh, numh, listh, H, S,  &
       LDOS_DM, energy, broadening)
  
  use precision, only  : dp
  use fdf
  use units,       only: eV, pi
  use sys,         only: die
  use m_mpi_utils, only: broadcast
  use parallel, only   : SIESTA_worker, BlockSize
  use parallel, only   : SIESTA_Group, SIESTA_Comm
  use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
  use class_Distribution
  use alloc,             only: re_alloc, de_alloc
  use mpi_siesta
  use iso_c_binding
  use f_ppexsi_interface, only: f_ppexsi_options
  use f_ppexsi_interface, only: f_ppexsi_plan_finalize
  use f_ppexsi_interface, only: f_ppexsi_plan_initialize
  use f_ppexsi_interface, only: f_ppexsi_selinv_complex_symmetric_matrix
  use f_ppexsi_interface, only: f_ppexsi_load_real_symmetric_hs_matrix
  use f_ppexsi_interface, only: f_ppexsi_set_default_options
  use f_ppexsi_interface, &
        only: f_ppexsi_symbolic_factorize_complex_symmetric_matrix
  
  implicit          none
  
  integer, intent(in)          :: maxnh, no_u, no_l, nspin_in
  integer, intent(in), target  :: listh(maxnh), numh(no_l)
  real(dp), intent(in), target :: H(maxnh,nspin_in), S(maxnh)
  real(dp), intent(in)         :: energy, broadening
  real(dp), intent(out)        :: LDOS_DM(maxnh,nspin_in)
  integer        :: ih, i
  integer        :: info
  logical        :: write_ok
  !------------
  external         :: timer
  integer          :: World_Comm, mpirank, ierr
  !
  integer   :: norbs
  integer   :: nspin
  integer :: PEXSI_Pole_Group, PEXSI_Spatial_Group, World_Group
  integer, allocatable :: pexsi_pole_ranks_in_world(:)
  integer  :: nworkers_SIESTA
  integer, allocatable :: siesta_ranks_in_world(:)
  integer :: PEXSI_Pole_Group_in_World
  integer, allocatable :: PEXSI_Pole_ranks_in_World_Spin(:,:)
  integer :: PEXSI_Pole_Comm, PEXSI_Spatial_Comm, PEXSI_Spin_Comm
  integer :: numNodesTotal
  integer :: npPerPole
  logical  :: PEXSI_worker
  !
  type(Distribution)   :: dist1
  type(Distribution), allocatable, target   :: dist2_spin(:)
  type(Distribution), pointer :: dist2
  integer  :: pbs, color, spatial_rank, spin_rank
  type(aux_matrix), allocatable, target :: m1_spin(:)
  type(aux_matrix) :: m2
  type(aux_matrix), pointer :: m1
  integer :: nrows, nnz, nnzLocal, numColLocal
  integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
  !
  real(dp), pointer, dimension(:) :: &
          HnzvalLocal=>null(), SnzvalLocal=>null(),  &
          DMnzvalLocal => null()
  !
  integer :: ispin, pexsi_spin
  type(f_ppexsi_options) :: options
  !
  integer                :: isSIdentity
  integer                :: verbosity
  integer(c_intptr_t)    :: plan
    integer :: numProcRow, numProcCol
    integer :: outputFileIndex
  real(dp), pointer, dimension(:) :: AnzvalLocal => null()
  real(dp), pointer, dimension(:) :: AinvnzvalLocal => null()
  integer :: loc
  !  --------  for serial compilation
  !
  ! Our global communicator is a duplicate of the passed communicator
  !
  call MPI_Comm_Dup(true_MPI_Comm_World, World_Comm, ierr)
  call mpi_comm_rank( World_Comm, mpirank, ierr )
  
  call timer("pexsi-ldos", 1)  
  
  if (SIESTA_worker) then
     ! rename some intent(in) variables, which are only
     ! defined for the Siesta subset
     norbs = no_u
     nspin = nspin_in
  endif
  !
  call broadcast(norbs,comm=World_Comm)
  call broadcast(nspin,World_Comm)
  call MPI_Comm_Group(World_Comm,World_Group, ierr)
  call MPI_Group_Size(SIESTA_Group, nworkers_SIESTA, ierr)
  allocate(siesta_ranks_in_world(nworkers_SIESTA))
  call MPI_Group_translate_ranks( SIESTA_Group, nworkers_SIESTA, &
       (/ (i,i=0,nworkers_SIESTA-1) /), &
       World_Group, siesta_ranks_in_world, ierr )
  call newDistribution(dist1,World_Comm,siesta_ranks_in_world, &
                       TYPE_BLOCK_CYCLIC,BlockSize,"bc dist")
  deallocate(siesta_ranks_in_world)
  call MPI_Barrier(World_Comm,ierr)
  
  call mpi_comm_size( World_Comm, numNodesTotal, ierr )
  
  npPerPole  = fdf_get("PEXSI.np-per-pole",4)
  npPerPole  = fdf_get("PEXSI.LDOS.np-per-pole",npPerPole)
  if (nspin*npPerPole > numNodesTotal) &
       call die("PEXSI.np-per-pole is too big for MPI size")
  
  ! "Row" communicator for independent PEXSI operations on each spin
  ! The name refers to "spatial" degrees of freedom.
  color = mod(mpirank,nspin)    ! {0,1} for nspin = 2, or {0} for nspin = 1
  call MPI_Comm_Split(World_Comm, color, mpirank, PEXSI_Spatial_Comm, ierr)
  
  ! "Column" communicator for spin reductions
  color = mpirank/nspin       
  call MPI_Comm_Split(World_Comm, color, mpirank, PEXSI_Spin_Comm, ierr)
  
  ! Group and Communicator for first-pole team of PEXSI workers
  !
  call MPI_Comm_Group(PEXSI_Spatial_Comm, PEXSI_Spatial_Group, Ierr)
  call MPI_Group_incl(PEXSI_Spatial_Group, npPerPole,   &
       (/ (i,i=0,npPerPole-1) /),&
       PEXSI_Pole_Group, Ierr)
  call MPI_Comm_create(PEXSI_Spatial_Comm, PEXSI_Pole_Group,&
       PEXSI_Pole_Comm, Ierr)
  
  
  call mpi_comm_rank( PEXSI_Spatial_Comm, spatial_rank, ierr )
  call mpi_comm_rank( PEXSI_Spin_Comm, spin_rank, ierr )
  PEXSI_worker = (spatial_rank < npPerPole)   ! Could be spin up or spin down
  
  ! PEXSI blocksize 
  pbs = norbs/npPerPole
  
  ! Careful with this. For the purposes of matrix transfers,
  ! we need the ranks of the Pole group
  ! in the "bridge" communicator/group (World)
  
  allocate(pexsi_pole_ranks_in_world(npPerPole))
  call MPI_Comm_Group(World_Comm, World_Group, Ierr)
  
  call MPI_Group_translate_ranks( PEXSI_Pole_Group, npPerPole, &
       (/ (i,i=0,npPerPole-1) /), &
       World_Group, pexsi_pole_ranks_in_world, ierr )
  
  ! What we need is to include the actual world ranks
  ! in the distribution object
  allocate (PEXSI_Pole_ranks_in_World_Spin(npPerPole,nspin))
  call MPI_AllGather(pexsi_pole_ranks_in_world,npPerPole,MPI_integer,&
       PEXSI_Pole_Ranks_in_World_Spin(1,1),npPerPole, &
       MPI_integer,PEXSI_Spin_Comm,ierr)
  
  ! Create distributions known to all nodes
  allocate(dist2_spin(nspin))
  do ispin = 1, nspin
     call newDistribution(dist2_spin(ispin), World_Comm, &
          PEXSI_Pole_Ranks_in_World_Spin(:,ispin),  &
          TYPE_PEXSI, pbs, "px dist")
  enddo
  deallocate(pexsi_pole_ranks_in_world,PEXSI_Pole_Ranks_in_World_Spin)
  call MPI_Barrier(World_Comm,ierr)
  
  
  pexsi_spin = spin_rank+1  ! {1,2}
  ! This is done serially on the Siesta side, each time
  ! filling in the structures in one PEXSI set
  
  allocate(m1_spin(nspin))
  do ispin = 1, nspin
  
     m1 => m1_spin(ispin)
  
     if (SIESTA_worker) then
        m1%norbs = norbs
        m1%no_l  = no_l
        m1%nnzl  = sum(numH(1:no_l))
        m1%numcols => numH
        m1%cols    => listH
        allocate(m1%vals(2))
        m1%vals(1)%data => S(:)
        m1%vals(2)%data => H(:,ispin)
  
     endif  ! SIESTA_worker
  
     call timer("redist_orbs_fwd", 1)
  
     ! Note that we cannot simply wrap this in a pexsi_spin test, as
     ! there are Siesta nodes in both spin sets.
     ! We must discriminate the PEXSI workers by the distribution info
     dist2 => dist2_spin(ispin)
     call redistribute_spmatrix(norbs,m1,dist1,m2,dist2,World_Comm)
     
     call timer("redist_orbs_fwd", 2)
  
     if (PEXSI_worker .and. (pexsi_spin == ispin) ) then
  
        nrows = m2%norbs          ! or simply 'norbs'
        numColLocal = m2%no_l
        nnzLocal    = m2%nnzl
        call MPI_AllReduce(nnzLocal,nnz,1,MPI_integer,MPI_sum,PEXSI_Pole_Comm,ierr)
  
        call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_ldos")
        colptrLocal(1) = 1
        do ih = 1,numColLocal
           colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
        enddo
  
        rowindLocal => m2%cols
        SnzvalLocal => m2%vals(1)%data
        HnzvalLocal => m2%vals(2)%data
  
        call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","pexsi_ldos")
  
        call memory_all("after transferring H+S for PEXSI-LDOS",PEXSI_Pole_Comm)
  
     endif ! PEXSI worker
  enddo
  
  ! Make these available to all
  ! (Note that the values are those on process 0, which is in the spin=1 set
  ! In fact, they are only needed for calls to the interface, so the broadcast
  ! could be over PEXSI_Spatial_Comm only.
  
  call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
  call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)
  
  call memory_all("after setting up H+S for PEXSI LDOS",World_comm)
  
  ! -- 
    isSIdentity = 0
  !
    call f_ppexsi_set_default_options( options )
    ! Ordering flag:
    !   1: Use METIS
    !   0: Use PARMETIS/PTSCOTCH
    options%ordering = fdf_get("PEXSI.ordering",1)
    ! Number of processors for symbolic factorization
    ! Only relevant for PARMETIS/PT_SCOTCH
    options%npSymbFact = fdf_get("PEXSI.np-symbfact",1)
    options%verbosity = fdf_get("PEXSI.verbosity",1)
    call get_row_col(npPerPole,numProcRow,numProcCol)
  
    ! Set the outputFileIndex to be the pole index.
    ! Starting from PEXSI v0.8.0, the first processor for each pole outputs
    ! information
  
    if( mod( mpirank, npPerPole ) .eq. 0 ) then
      outputFileIndex = mpirank / npPerPole;
    else
      outputFileIndex = -1;
    endif
  !
  ! Note that even though we only use one pole's worth of processors, we
  ! still use the full spatial PEXSI communicator in the plan.
  ! Failing to do so leads to an error. This is not sufficiently documented.
  !
    plan = f_ppexsi_plan_initialize(&
        PEXSI_Spatial_Comm,&
        numProcRow,&
        numProcCol,&
        outputFileIndex,&
        info) 
  
  call check_info(info,"plan_initialize in LDOS")
  call f_ppexsi_load_real_symmetric_hs_matrix(&
        plan,&
        options,&
        nrows,&
        nnz,&
        nnzLocal,&
        numColLocal,&
        colptrLocal,&
        rowindLocal,&
        HnzvalLocal,&
        isSIdentity,&
        SnzvalLocal,&
        info) 
  
  call check_info(info,"load_real_sym_hs_matrix in LDOS")
  
    call f_ppexsi_symbolic_factorize_complex_symmetric_matrix(&
         plan, &
         options,&
         info)
    call check_info(info,"factorize complex matrix in LDOS")
  
  options%isSymbolicFactorize = 0 ! We do not need it anymore
  if (PEXSI_worker) then
  
     if(mpirank == 0) then
        write(6,"(a,f16.5,f10.5)") &
             'Calling PEXSI LDOS routine. Energy and broadening (eV) ', &
             energy/eV, broadening/eV
        write(6,"(a,i4)") &
             'Processors working on selected inversion: ', npPerPole
     endif
  
     call timer("pexsi-ldos-selinv", 1)
  
     ! Build AnzvalLocal as H-zS, where z=(E,broadening), and treat
     ! it as a one-dimensional real array with 2*nnzlocal entries
  
     call re_alloc(AnzvalLocal,1,2*nnzLocal,"AnzvalLocal","pexsi_ldos")
     call re_alloc(AinvnzvalLocal,1,2*nnzLocal,"AinvnzvalLocal","pexsi_ldos")
  
     loc = 1
     do i = 1, nnzLocal
        AnzvalLocal(loc) = Hnzvallocal(i) - energy*Snzvallocal(i)
        AnzvalLocal(loc+1) =  - broadening*Snzvallocal(i)
        loc = loc + 2
     enddo
  
     call f_ppexsi_selinv_complex_symmetric_matrix(&
          plan,&
          options,&
          AnzvalLocal,&
          AinvnzvalLocal,&
          info) 
  
     call check_info(info,"selinv complex matrix in LDOS")
  
     ! Get DMnzvalLocal as 1/pi * Imag(Ainv...)
  
     loc = 1
     do i = 1, nnzLocal
        DMnzvalLocal(i) = (1.0_dp/pi) * AinvnzvalLocal(loc+1)
        loc = loc + 2
     enddo
     call de_alloc(AnzvalLocal,"AnzvalLocal","pexsi_ldos")
     call de_alloc(AinvnzvalLocal,"AinvnzvalLocal","pexsi_ldos")
  
     call timer("pexsi-ldos-selinv", 2)
     !
  endif ! PEXSI_worker
  
  do ispin = 1, nspin
  
     m1 => m1_spin(ispin)
  
     if (PEXSI_worker .and. (pexsi_spin == ispin)) then
        ! Prepare m2 to transfer
  
        call de_alloc(colPtrLocal,"colPtrLocal","pexsi_ldos")
  
        call de_alloc(m2%vals(1)%data,"m2%vals(1)%data","pexsi_ldos")
        call de_alloc(m2%vals(2)%data,"m2%vals(2)%data","pexsi_ldos")
  
       deallocate(m2%vals)
       allocate(m2%vals(1))
       m2%vals(1)%data => DMnzvalLocal(1:nnzLocal)
  
     endif
  
     ! Prepare m1 to receive the results
     if (SIESTA_worker) then
        nullify(m1%vals(1)%data)    ! formerly pointing to S
        nullify(m1%vals(2)%data)    ! formerly pointing to H
        deallocate(m1%vals)
        nullify(m1%numcols)         ! formerly pointing to numH
        nullify(m1%cols)            ! formerly pointing to listH
     endif
  
     call timer("redist_orbs_bck", 1)
     dist2 => dist2_spin(ispin)
     call redistribute_spmatrix(norbs,m2,dist2,m1,dist1,World_Comm)
     call timer("redist_orbs_bck", 2)
  
     if (PEXSI_worker .and. (pexsi_spin == ispin)) then
        call de_alloc(DMnzvalLocal, "DMnzvalLocal", "pexsi_ldos")
  
        nullify(m2%vals(1)%data)    ! formerly pointing to DM
        deallocate(m2%vals)
        ! allocated in the direct transfer
        call de_alloc(m2%numcols,"m2%numcols","pexsi_ldos")
        call de_alloc(m2%cols,   "m2%cols",   "pexsi_ldos")
     endif
  
     if (SIESTA_worker) then
  
        LDOS_DM(:,ispin)  = m1%vals(1)%data(:)    
        ! Check no_l
        if (no_l /= m1%no_l) then
           call die("Mismatch in no_l")
        endif
        ! Check listH
        if (any(listH(:) /= m1%cols(:))) then
           call die("Mismatch in listH")
        endif
  
        call de_alloc(m1%vals(1)%data,"m1%vals(1)%data","pexsi_ldos")
        deallocate(m1%vals)
        ! allocated in the direct transfer
        call de_alloc(m1%numcols,"m1%numcols","pexsi_ldos") 
        call de_alloc(m1%cols,   "m1%cols",   "pexsi_ldos")
  
     endif
  enddo
  call timer("pexsi-ldos", 2)
  
  
  call delete(dist1)
  do ispin = 1, nspin
     call delete(dist2_spin(ispin))
  enddo
  deallocate(dist2_spin)
  deallocate(m1_spin)
  
  call MPI_Comm_Free(PEXSI_Spatial_Comm, ierr)
  call MPI_Comm_Free(PEXSI_Spin_Comm, ierr)
  call MPI_Comm_Free(World_Comm, ierr)
  
  ! This communicator was created from a subgroup.
  ! As such, it is MPI_COMM_NULL for those processes
  ! not in the subgroup (non PEXSI_workers). Only
  ! defined communicators can be freed
  if (PEXSI_worker) then
     call MPI_Comm_Free(PEXSI_Pole_Comm, ierr)
  endif
  
  ! We finalize the plan here
  call f_ppexsi_plan_finalize( plan, info )
  
  call MPI_Group_Free(PEXSI_Spatial_Group, ierr)
  call MPI_Group_Free(PEXSI_Pole_Group, ierr)
  call MPI_Group_Free(World_Group, ierr)
  
  CONTAINS
      
  subroutine get_row_col(np,nrow,ncol)
  integer, intent(in)  :: np
  integer, intent(out) :: nrow, ncol
  !
  ! Finds the factors nrow and ncol such that nrow*ncol=np,
  ! are as similar as possible, and nrow>=ncol.
  ! For prime np, ncol=1, nrow=np.
  
  ncol  = floor(sqrt(dble(np)))
  do
    nrow = np/ncol
    if (nrow*ncol == np) exit
    ncol = ncol - 1
  enddo
  end subroutine get_row_col
  
  subroutine check_info(info,str)
  integer, intent(in) :: info
  character(len=*), intent(in) :: str
  
      if(mpirank == 0) then
         if (info /= 0) then
            write(6,*) trim(str) // " info : ", info
            call die("Error exit from " // trim(str) // " routine")
         endif
        call pxfflush(6)
      endif       
  end subroutine check_info 
  
  end subroutine get_LDOS_SI
#endif
  end module m_pexsi_local_dos
! End of tangled code
