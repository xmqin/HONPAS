
! Tangled code
module m_pexsi_dos
#ifdef SIESTA__PEXSI
    use precision, only  : dp

  public :: pexsi_dos

CONTAINS

subroutine pexsi_dos(no_u, no_l, nspin_in,  &
     maxnh, numh, listhptr, listh, H, S, qtot, ef_in)

    use fdf
    use parallel, only   : SIESTA_worker, BlockSize
    use parallel, only   : SIESTA_Group, SIESTA_Comm
    use m_mpi_utils, only: globalize_max
    use m_mpi_utils, only: broadcast
    use units,       only: eV
    use m_redist_spmatrix, only: aux_matrix, redistribute_spmatrix
    use class_Distribution
    use alloc,             only: re_alloc, de_alloc
#ifdef MPI
    use mpi_siesta
#endif

    use iso_c_binding
    use f_ppexsi_interface, only: f_ppexsi_options
    use f_ppexsi_interface, only: f_ppexsi_plan_finalize
    use f_ppexsi_interface, only: f_ppexsi_plan_initialize
    use f_ppexsi_interface, only: f_ppexsi_inertia_count_real_symmetric_matrix
    use f_ppexsi_interface, only: f_ppexsi_load_real_symmetric_hs_matrix
    use f_ppexsi_interface, only: f_ppexsi_set_default_options
    use f_ppexsi_interface, &
          only: f_ppexsi_symbolic_factorize_real_symmetric_matrix

  implicit          none

  integer, intent(in)  :: maxnh, no_u, no_l, nspin_in
  integer, intent(in), target  :: listh(maxnh), numh(no_l), listhptr(no_l)
  real(dp), intent(in), target :: H(maxnh,nspin_in), S(maxnh)
  real(dp), intent(in)        :: qtot ! Total number of electrons
  real(dp), intent(in)        :: ef_in  ! Fermi energy
integer        :: ih, i
integer        :: info
logical        :: write_ok
!------------
external         :: timer
integer          :: World_Comm, mpirank, ierr
!
real(dp)  :: numElectronExact, ef
integer   :: norbs, scf_step
!
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
integer(c_intptr_t)    :: plan
  integer :: numProcRow, numProcCol
  integer :: outputFileIndex
type(aux_matrix), allocatable, target :: m1_spin(:)
type(aux_matrix) :: m2
type(aux_matrix), pointer :: m1
integer :: nrows, nnz, nnzLocal, numColLocal
integer, pointer, dimension(:) ::  colptrLocal=> null(), rowindLocal=>null()
!
real(dp), pointer, dimension(:) :: &
        HnzvalLocal=>null(), SnzvalLocal=>null(),  &
        DMnzvalLocal => null() , EDMnzvalLocal => null(), &
        FDMnzvalLocal => null()
!
integer :: ispin, pexsi_spin
type(f_ppexsi_options) :: options
!
integer                :: isSIdentity
integer                :: verbosity
integer  :: npoints, nShifts, numNodesSpatial
real(dp), allocatable :: edos(:), intdos(:)
real(dp) :: emin, emax, delta
logical  :: ef_reference
real(dp), allocatable :: intdos_spin(:,:)
integer :: j, lun
!  --------  for serial compilation
#ifndef MPI
    call die("PEXSI needs MPI")
#else
!
! Our global communicator is a duplicate of the passed communicator
!
call MPI_Comm_Dup(true_MPI_Comm_World, World_Comm, ierr)
call mpi_comm_rank( World_Comm, mpirank, ierr )

call timer("pexsi_dos", 1)  

if (SIESTA_worker) then

   ! rename and broadcast some intent(in) variables, which are only
   ! defined for the Siesta subset

   norbs = no_u
   nspin = nspin_in
   numElectronExact = qtot 
   ef = ef_in

endif
!
call broadcast(norbs,comm=World_Comm)
call broadcast(numElectronExact,World_Comm)
call broadcast(nspin,World_Comm)
call broadcast(ef,World_Comm)
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

! -- Prepare plan
call get_row_col(npPerPole,numProcRow,numProcCol)

  ! Set the outputFileIndex to be the pole index.
  ! Starting from PEXSI v0.8.0, the first processor for each pole outputs
  ! information

  if( mod( mpirank, npPerPole ) .eq. 0 ) then
    outputFileIndex = mpirank / npPerPole;
  else
    outputFileIndex = -1;
  endif

   plan = f_ppexsi_plan_initialize(&
      PEXSI_Spatial_Comm,&
      numProcRow,&
      numProcCol,&
      outputFileIndex,&
      info) 
   call check_info(info,"plan_initialize in DOS")

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

      call re_alloc(colptrLocal,1,numColLocal+1,"colptrLocal","pexsi_solver")
      colptrLocal(1) = 1
      do ih = 1,numColLocal
         colptrLocal(ih+1) = colptrLocal(ih) + m2%numcols(ih)
      enddo

      rowindLocal => m2%cols
      SnzvalLocal => m2%vals(1)%data
      HnzvalLocal => m2%vals(2)%data

      call re_alloc(DMnzvalLocal,1,nnzLocal,"DMnzvalLocal","pexsi_solver")
      call re_alloc(EDMnzvalLocal,1,nnzLocal,"EDMnzvalLocal","pexsi_solver")
      call re_alloc(FDMnzvalLocal,1,nnzLocal,"FDMnzvalLocal","pexsi_solver")

      call memory_all("after setting up H+S for PEXSI (PEXSI_workers)",PEXSI_Pole_Comm)

   endif ! PEXSI worker
enddo

! Make these available to all
! (Note that the values are those on process 0, which is in the spin=1 set
! In fact, they are only needed for calls to the interface, so the broadcast
! could be over PEXSI_Spatial_Comm only.

call MPI_Bcast(nrows,1,MPI_integer,0,World_Comm,ierr)
call MPI_Bcast(nnz,1,MPI_integer,0,World_Comm,ierr)

call memory_all("after setting up H+S for PEXSI",World_comm)

!
call f_ppexsi_set_default_options( options )

isSIdentity = 0

! Ordering flag:
!   1: Use METIS
!   0: Use PARMETIS/PTSCOTCH
options%ordering = fdf_get("PEXSI.ordering",1)

! Number of processors for symbolic factorization
! Only relevant for PARMETIS/PT_SCOTCH
options%npSymbFact = fdf_get("PEXSI.np-symbfact",1)

verbosity = fdf_get("PEXSI.verbosity",1)
options%verbosity = verbosity

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

call check_info(info,"load_real_sym_hs_matrix")


call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
     plan, &
     options,&
     info)
call check_info(info,"symbolic_factorize_real_symmetric_matrix")

options%isSymbolicFactorize = 0 ! We do not need it anymore
call mpi_comm_size( PEXSI_Spatial_Comm, numNodesSpatial, ierr )
nShifts = numNodesSpatial/npPerPole

npoints = fdf_get("PEXSI.DOS.npoints",200)

npoints = ceiling(dble(npoints)/nShifts) * nShifts

allocate (edos(npoints),intdos(npoints))
emin = fdf_get("PEXSI.DOS.emin",-1.0_dp,"Ry")
emax = fdf_get("PEXSI.DOS.emax",+1.0_dp,"Ry")

ef_reference = fdf_get("PEXSI.DOS.Ef.Reference",.true.)
if (ef_reference) then
   emin = emin + ef
   emax = emax + ef
endif
  delta = (emax-emin)/(npoints-1)
  do j = 1, npoints
     edos(j) = emin + (j-1)*delta
  enddo

call timer("pexsi-raw-inertia-ct", 1)

if(mpirank == 0) then
   write (6,"(a,f12.5,a,f12.5,a,a,i4)")  &
                'Calling inertia_count for DOS: [', &
                emin/eV, ",", emax/eV, "] (eV)", &
                " Nshifts: ", npoints
endif

call f_ppexsi_inertia_count_real_symmetric_matrix(&
  plan,&
  options,&
  npoints,&
  edos,&
  intdos,&
  info) 
call check_info(info,"inertia_count_real_symmetric_matrix in DOS")

call timer("pexsi-raw-inertia-ct", 2)

allocate(intdos_spin(npoints,nspin))
call MPI_Gather( intdos, npoints, MPI_double_precision, &
     intdos_spin(1,1), npoints, MPI_double_precision, &
     0, PEXSI_Spin_Comm, ierr )

if (mpirank == 0) then
   call io_assign(lun)
   open(unit=lun,file="PEXSI_INTDOS",form="formatted",status="unknown", &
        position="rewind",action="write")
   write(lun,"(2f15.6,i6,i2,a)") ef/eV, qtot, npoints, nspin, &
        "# (Ef, qtot, npoints, nspin) / npoints lines: E(eV), IntDos(E)"
   ! Note that the inertia count procedure counts the number of eigenvalues,
   ! so for nspin=1 we have to multiply by two to get the number of states
   do j=1,npoints
      if (nspin == 1) then
         write(lun,"(f15.6,f15.2)") edos(j)/eV, 2*intdos(j)
      else
         write(lun,"(f15.6,2f15.2)") edos(j)/eV, intdos_spin(j,1), intdos_spin(j,2)
      endif
   enddo
endif

call timer("pexsi_dos", 2)


deallocate(edos,intdos,intdos_spin)

if (SIESTA_worker) then
   ! Its pointers were not actually allocated
   deallocate(m1%vals)
endif

call delete(dist1)
do ispin = 1, nspin
   call delete(dist2_spin(ispin))
enddo
deallocate(dist2_spin)
deallocate(m1_spin)

if (PEXSI_worker) then

   call de_alloc(colptrLocal,"colptrLocal","pexsi_dos")

   call de_alloc(m2%numcols,"m2%numcols","m_pexsi_dos")
   call de_alloc(m2%cols,"m2%cols","m_pexsi_dos")
   do j=1,size(m2%vals)
      call de_alloc(m2%vals(j)%data,"m2%vals(j)%data","m_pexsi_dos")
   enddo
   deallocate(m2%vals)

endif

call f_ppexsi_plan_finalize( plan, info )

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

call MPI_Group_Free(PEXSI_Spatial_Group, ierr)
call MPI_Group_Free(PEXSI_Pole_Group, ierr)
call MPI_Group_Free(World_Group, ierr)
#endif

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

end subroutine pexsi_dos
#endif
end module m_pexsi_dos
! End of tangled code
