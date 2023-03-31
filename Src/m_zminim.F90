! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_zminim

use fdf,            only : fdf_boolean, fdf_integer, fdf_get, fdf_physical
use files,          only : slabel
use parallel,       only : ProcessorY, BlockSize, Node, Nodes
use precision,      only : dp
use siesta_options, only : fixspin
use sys,            only : die
#ifdef MPI
use mpi_siesta,     only : mpi_integer, mpi_double_precision, mpi_double_complex, mpi_comm_world, mpi_sum, mpi_status_size
use parallelsubs,   only : GetNodeOrbs, GlobalToLocalOrb, WhichNodeOrb
use parallelsubs,   only : set_BlockSizeDefault
#endif

implicit none

external :: io_assign, io_close

!**** PRIVATE ***********************************!

private

type multispin
  complex(dp), allocatable :: mtrx(:,:)
end type multispin

real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp

complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
complex(dp), parameter :: cmplx_2=(2.0_dp,0.0_dp)
complex(dp), parameter :: cmplx_half=(0.5_dp,0.0_dp)
complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

logical, save :: LongOut                     ! print detailed output?
logical, save :: UseCholesky                 ! use Cholesky factorization?
logical, save :: UseSparse                   ! use sparse algebra?
logical, save :: WriteCoeffs                 ! write WF coeffs. at the end of each MD iteration?
logical, save :: ReadCoeffs                  ! read WF coeffs. at the start of the calculation?
logical, allocatable, save :: FirstCall(:,:) ! first time minim_cg is called? (k point, up/down spin)
#ifdef MPI
logical, save :: Use2D                       ! use 2D data decomposition?
#endif

integer, save :: N_occ_loc_1D(1:2)         ! num. of WFs (local 1D) (up/down spin)
integer, save :: h_dim_loc_1D              ! num. of AOs (local 1D)
integer, save :: N_occ_loc(1:2,1:2)        ! num. of WFs (local 1D or 2D) (up/down spin)
integer, save :: h_dim_loc(1:2)            ! num. of AOs (local 1D or 2D) (up/down spin)
integer, save :: N_occ_diff                ! difference between num. of WFs for up/down spin
#ifdef MPI
integer, save :: BlockSize_c               ! ScaLAPACK blocking factor for WFs
integer, save :: desc1(1:9)                ! descriptor for operator matrix in AO basis (1D or 2D)
integer, save :: desc2(1:9,1:2)            ! descriptor for operator matrix in WF basis (1D or 2D) (up/down spin)
integer, save :: desc3(1:9,1:2)            ! descriptor for WF coeffs. matrix (1D or 2D) (up/down spin)
integer, allocatable :: numh_recv(:)       ! num. of nonzero elements of each row of sparse matrices
integer, allocatable :: listhptr_recv(:)   ! pointer to start of row in listh
integer, allocatable :: listh_recv(:)      ! list of nonzero elements of each row of sparse matrices
integer, allocatable :: h_dim_l2g_recv(:)  ! local-to-global index transform for AOs
integer, allocatable, save :: h_dim_l2g(:) ! local-to-global index transform for AOs
#endif

#ifdef MPI
complex(dp), allocatable :: As_recv(:) ! work matrix
#endif

type(multispin), allocatable, save :: c(:,:)      ! WF coeffs. matrix (after Cholesky fact.) (k point, up/down spin)
type(multispin), allocatable, save :: c_orig(:,:) ! WF coeffs. matrix (before Cholesky fact.) (k point, up/down spin)

!**** PUBLIC ************************************!

public :: zminim

!************************************************!

contains

!================================================!
! use the orbital minimization method (OMM) to   !
! solve the eigenvalue problem (double complex   !
! routine, for extended k grid calculations)     !
!================================================!
subroutine zminim(CalcE,PreviousCallDiagon,iscf,istp,nbasis,nspin,h_dim,nhmax,numh,listhptr,listh,d_sparse,eta,qs,no,xij,indxuo,nk,&
                  kpoint,wk,h_sparse,s_sparse,t_sparse)
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! Calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?

  integer, intent(in) :: iscf               ! SCF iteration num.
  integer, intent(in) :: istp               ! MD iteration num.
  integer, intent(in) :: nbasis             ! dimension of numh and listhptr
  integer, intent(in) :: nspin              ! num. of spins
  integer, intent(in) :: h_dim              ! num. of AOs (global)
  integer, intent(in) :: nhmax              ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(1:nbasis)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(1:nbasis) ! pointer to start of row in listh
  integer, intent(in) :: listh(1:nhmax)     ! list of nonzero elements of each row of sparse matrices
  integer, intent(in) :: no                 ! num. of AOs in supercell (global)
  integer, intent(in) :: indxuo(1:no)       ! index of equivalent AO in unit cell
  integer, intent(in) :: nk                 ! num. of k points

  real(dp), intent(in) :: qs(1:2)                             ! num. of electrons per spin
  real(dp), intent(in) :: eta(1:2)                            ! chemical potential for Kim functional
  real(dp), intent(in) :: xij(1:3,1:nhmax)                    ! vectors between AO centres (sparse)
  real(dp), intent(in) :: kpoint(1:3,1:nk)                    ! k point vectors
  real(dp), intent(in) :: wk(1:nk)                            ! k point weights
  real(dp), intent(in), optional :: h_sparse(1:nhmax,1:nspin) ! hamiltonian matrix (sparse)
  real(dp), intent(in), optional :: t_sparse(1:nhmax)         ! kinetic energy matrix (sparse)
  real(dp), intent(in), optional :: s_sparse(1:nhmax)         ! overlap matrix (sparse)

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: d_sparse(1:nhmax,1:nspin) ! (energy-)density matrix (sparse)

  !**** LOCAL ***********************************!

  logical :: UpdatePrecon     ! update the preconditioner?
  logical :: UsePrecon        ! use the preconditioner?
  logical :: UpdateSparseComm ! update nhmax_max?

  integer :: ik
  integer :: ispin
  integer :: io
  integer :: juo
  integer :: jo
  integer :: j
  integer :: k
  integer :: ind
  integer, save :: N_occ(1:2)                         ! num. of WFs (global) (up/down spin)
  integer, save :: precon_default                     ! default number of SCF steps for which to use preconditioning
  integer, save :: precon_first_step                  ! number of SCF steps for which to use preconditioning (first MD step)
  integer, save :: last_call(1:2)=0                   ! iscf and istp value of last module call
  integer, allocatable, save :: last_precon_update(:) ! istp value of last preconditioner update (k point)
#ifdef MPI
  integer :: BlockSize_c_default
#endif

  real(dp), save :: sum_wk
  real(dp), allocatable, save :: kxij(:,:,:)

  complex(dp), allocatable :: hk_sparse(:)   ! hamiltonian matrix at k point (sparse)
  complex(dp), allocatable :: sk_sparse(:)   ! overlap matrix at k point (sparse)
  complex(dp), allocatable :: dk_sparse(:)   ! (energy-)density matrix at k point (sparse)
  complex(dp), allocatable :: h_dense1D(:,:) ! hamiltonian matrix (dense 1D)
  complex(dp), allocatable :: t_dense1D(:,:) ! kinetic energy matrix (dense 1D)
  complex(dp), allocatable :: s_dense1D(:,:) ! overlap matrix (dense 1D)
  complex(dp), allocatable :: d_dense1D(:,:) ! (energy-)density matrix (dense 1D)

  !**********************************************!

  call timer('zmin',1)

  if (.not. allocated(FirstCall)) then
    allocate(FirstCall(1:nk,1:2))
    FirstCall=.true.
  end if
  if (.not. allocated(c)) allocate(c(1:nk,1:2))
  if (.not. allocated(c_orig)) allocate(c_orig(1:nk,1:2))

  if (all(FirstCall(1:nk,1:2))) then

    LongOut=fdf_boolean('OMM.LongOutput',.false.)

    UseCholesky=fdf_boolean('OMM.UseCholesky',.true.)
    UseSparse=fdf_boolean('OMM.UseSparse',.false.)
    if (UseSparse) call die('ERROR: sparse algebra not yet implemented for k points!')
    if (UseSparse .and. UseCholesky) call die('ERROR: sparse algebra not compatible with Cholesky factorization!')

    WriteCoeffs=fdf_boolean('OMM.WriteCoeffs',.false.)
    ReadCoeffs=fdf_boolean('OMM.ReadCoeffs',.false.)

    if (.not. UseCholesky) then
      precon_default=fdf_integer('OMM.Precon',-1)
      precon_first_step=fdf_integer('OMM.PreconFirstStep',precon_default)
    end if

    if (nspin==1) then
      N_occ(1)=nint(0.5_dp*qs(1))
      if (abs(N_occ(1)-0.5_dp*qs(1))>1.0d-10) call die('ERROR: only integer number of electrons per spin allowed!')
    else
      if (.not. fixspin) call die('ERROR: OMM for spin unpolarized calculations only supports fixed spin!')
      do ispin=1,nspin
        N_occ(ispin)=nint(qs(ispin))
        if (abs(N_occ(ispin)-qs(ispin))>1.0d-10) call die('ERROR: only integer number of electrons per spin allowed!')
      end do
      N_occ_diff=N_occ(1)-N_occ(2)
    end if

    h_dim_loc_1D=nbasis

    sum_wk=0.0_dp
    do ik=1,nk
      sum_wk=sum_wk+wk(ik)
    end do

    allocate(last_precon_update(1:nk))
    last_precon_update=0

#ifdef MPI
    Use2D=fdf_boolean('OMM.Use2D',.true.)
    if (UseSparse .and. Use2D) call die('ERROR: sparse algebra not compatible with 2D data decomposition!')

    ! calculate the ScaLAPACK blocking factor for distributing the WF coeffs. matrix
    call set_blocksizedefault(Nodes,N_occ(1),BlockSize_c_default)
    BlockSize_c=fdf_integer('OMM.BlockSize',BlockSize_c_default)
#endif

  end if

  if (.not. allocated(kxij)) then
    allocate(kxij(1:2,1:nhmax,1:nk))
    do ik=1,nk
      do ind=1,nhmax
        kxij(1,ind,ik)=cos(kpoint(1,ik)*xij(1,ind)+&
                           kpoint(2,ik)*xij(2,ind)+&
                           kpoint(3,ik)*xij(3,ind))
        kxij(2,ind,ik)=sin(kpoint(1,ik)*xij(1,ind)+&
                           kpoint(2,ik)*xij(2,ind)+&
                           kpoint(3,ik)*xij(3,ind))
      end do
    end do
  end if

  if (UseSparse) then
    allocate(dk_sparse(1:nhmax))
  else
    allocate(d_dense1D(1:h_dim,1:h_dim_loc_1D))
  end if

  if (CalcE) then

    UsePrecon=.false.
    UpdatePrecon=.false.
    UpdateSparseComm=.false.

  else

    if (UseSparse) then
      allocate(hk_sparse(1:nhmax))
    else
      allocate(h_dense1D(1:h_dim,1:h_dim_loc_1D))
    end if

    ! decide whether preconditioning should be used in this minimization call
    if (UseCholesky) then
      UsePrecon=.false.
      UpdatePrecon=.false.
    else
      if (istp==1) then
        if ((iscf<=precon_first_step) .or. (precon_first_step<0)) then
          UsePrecon=.true.
        else
          UsePrecon=.false.
        end if
      else
        if ((iscf<=precon_default) .or. (precon_default<0)) then
          UsePrecon=.true.
        else
          UsePrecon=.false.
        end if
      end if
    end if

  end if

  do ik=1,nk

    if (.not. calcE) then

      ! if this is the first time the module is called for this MD step, convert the new overlap
      ! matrix from sparse to dense
      if (UseSparse) then
        if (istp/=last_call(2)) then
          UpdateSparseComm=.true.
          if (.not. allocated(sk_sparse)) allocate(sk_sparse(1:nhmax))
          !! convert s_sparse to sk_sparse !!
        else
          UpdateSparseComm=.false.
        end if
        if (UsePrecon .and. (istp/=last_precon_update(ik))) then
          if (.not. allocated(s_dense1D)) allocate(s_dense1D(1:h_dim_loc_1D,1:h_dim))
          s_dense1D=cmplx_0
          do io=1,h_dim_loc_1D
            do j=1,numh(io)
              ind=listhptr(io)+j
              jo=listh(ind)
              juo=indxuo(jo)
              s_dense1D(io,juo)=s_dense1D(io,juo)+cmplx(+s_sparse(ind)*kxij(1,ind,ik),&
                                                        -s_sparse(ind)*kxij(2,ind,ik),dp)
            end do
          end do
        end if
      else
        if (istp/=last_call(2)) then
          if (.not. allocated(s_dense1D)) allocate(s_dense1D(1:h_dim,1:h_dim_loc_1D))
          s_dense1D=cmplx_0
          do io=1,h_dim_loc_1D
            do j=1,numh(io)
              ind=listhptr(io)+j
              jo=listh(ind)
              juo=indxuo(jo)
              s_dense1D(juo,io)=s_dense1D(juo,io)+cmplx(+s_sparse(ind)*kxij(1,ind,ik),&
                                                        -s_sparse(ind)*kxij(2,ind,ik),dp)
            end do
          end do
        end if
      end if

      ! if this is the first time we are using preconditioning for this MD step, convert also the
      ! new kinetic energy matrix from sparse to dense
      if (UsePrecon) then
        if (istp/=last_precon_update(ik)) then
          if (UseSparse) then
            if (.not. allocated(t_dense1D)) allocate(t_dense1D(1:h_dim_loc_1D,1:h_dim))
            t_dense1D=cmplx_0
            do io=1,h_dim_loc_1D
              do j=1,numh(io)
                ind=listhptr(io)+j
                jo=listh(ind)
                juo=indxuo(jo)
                t_dense1D(io,juo)=t_dense1D(io,juo)+cmplx(+t_sparse(ind)*kxij(1,ind,ik),&
                                                          -t_sparse(ind)*kxij(2,ind,ik),dp)
              end do
            end do
          else
            if (.not. allocated(t_dense1D)) allocate(t_dense1D(1:h_dim,1:h_dim_loc_1D))
            t_dense1D=cmplx_0
            do io=1,h_dim_loc_1D
              do j=1,numh(io)
                ind=listhptr(io)+j
                jo=listh(ind)
                juo=indxuo(jo)
                t_dense1D(juo,io)=t_dense1D(juo,io)+cmplx(+t_sparse(ind)*kxij(1,ind,ik),&
                                                          -t_sparse(ind)*kxij(2,ind,ik),dp)
              end do
            end do
          end if
          UpdatePrecon=.true.
          last_precon_update(ik)=istp
        else
          UpdatePrecon=.false.
        end if
      else
        UpdatePrecon=.false.
      end if

    end if

    do ispin=1,nspin

      if (.not. calcE) then
        ! convert the hamiltonian matrix from sparse to dense (and shift the eingevalue spectrum
        ! w.r.t. the chemical potential reference)
        if (UseSparse) then
          !! convert h_sparse to hk_sparse !!
        else
          h_dense1D=cmplx_0
          do io=1,h_dim_loc_1D
            do j=1,numh(io)
              ind=listhptr(io)+j
              jo=listh(ind)
              juo=indxuo(jo)
              h_dense1D(juo,io)=h_dense1D(juo,io)+cmplx(+(h_sparse(ind,ispin)-eta(ispin)*s_sparse(ind))*kxij(1,ind,ik),&
                                                        -(h_sparse(ind,ispin)-eta(ispin)*s_sparse(ind))*kxij(2,ind,ik),dp)
            end do
          end do
        end if
      end if

      ! call the routine to perform the energy minimization
      if (UseSparse) then
        if (CalcE) then
          call minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),&
                               nspin,ispin,UpdatePrecon,UsePrecon,UpdateSparseComm,nk,ik,kpoint,dk_sparse)
        else
          if (allocated(s_dense1D)) then
            if (allocated(t_dense1D)) then
              call minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),&
                                   nspin,ispin,UpdatePrecon,UsePrecon,UpdateSparseComm,nk,ik,kpoint,dk_sparse,hk_sparse,&
                                   sk_sparse,s_dense1D,t_dense1D)
            else
              call minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),&
                                   nspin,ispin,UpdatePrecon,UsePrecon,UpdateSparseComm,nk,ik,kpoint,dk_sparse,hk_sparse,&
                                   sk_sparse,s_dense1D)
            end if
          else
            if (allocated(t_dense1D)) then
              call minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),&
                                   nspin,ispin,UpdatePrecon,UsePrecon,UpdateSparseComm,nk,ik,kpoint,dk_sparse,hk_sparse,&
                                   sk_sparse,t_dense1D=t_dense1D)
            else
              call minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),&
                                   nspin,ispin,UpdatePrecon,UsePrecon,UpdateSparseComm,nk,ik,kpoint,dk_sparse,hk_sparse,sk_sparse)
            end if
          end if
        end if
      else
        if (CalcE) then
          call minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),nspin,ispin,UpdatePrecon,UsePrecon,nk,ik,&
                        kpoint,d_dense1D)
        else
          if (allocated(s_dense1D)) then
            if (allocated(t_dense1D)) then
              call minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),nspin,ispin,UpdatePrecon,UsePrecon,nk,&
                            ik,kpoint,d_dense1D,h_dense1D,s_dense1D,t_dense1D)
            else
              call minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),nspin,ispin,UpdatePrecon,UsePrecon,nk,&
                            ik,kpoint,d_dense1D,h_dense1D,s_dense1D)
            end if
          else
            if (allocated(t_dense1D)) then
              call minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),nspin,ispin,UpdatePrecon,UsePrecon,nk,&
                            ik,kpoint,d_dense1D,h_dense1D,t_dense1D=t_dense1D)
            else
              call minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ(ispin),eta(ispin),nspin,ispin,UpdatePrecon,UsePrecon,nk,&
                            ik,kpoint,d_dense1D,h_dense1D)
            end if
          end if
        end if
      end if

      ! convert the (energy-)density matrix from dense to sparse
      if (UseSparse) then
        !! convert dk_sparse to d_sparse !!
      else
        if ((ik==1) .and. (ispin==1)) d_sparse=0.0_dp
        do io=1,h_dim_loc_1D
          do j=1,numh(io)
            ind=listhptr(io)+j
            jo=listh(ind)
            juo=indxuo(jo)
            d_sparse(ind,ispin)=d_sparse(ind,ispin)+(wk(ik)/sum_wk)*(real(d_dense1D(juo,io),dp)*kxij(1,ind,ik)-&
                                                                     aimag(d_dense1D(juo,io))*kxij(2,ind,ik))
          end do
        end do
      end if

    end do
  end do

  if (.not. CalcE) then
    if (allocated(t_dense1D)) deallocate(t_dense1D)
    if (allocated(s_dense1D)) deallocate(s_dense1D)
    if (UseSparse) then
      if (allocated(sk_sparse)) deallocate(sk_sparse)
      deallocate(hk_sparse)
    else
      deallocate(h_dense1D)
    end if
  end if
  if (UseSparse) then
    deallocate(dk_sparse)
  else
    deallocate(d_dense1D)
  end if
  if (CalcE) deallocate(kxij)

  last_call(1:2)=(/iscf,istp/)

  call timer('zmin',2)

  end subroutine zminim

!================================================!
! minimize the Kim functional by conjugate       !
! gradients (dense routine)                      !
!================================================!
subroutine minim_cg(CalcE,PreviousCallDiagon,iscf,h_dim,N_occ,eta,nspin,ispin,UpdatePrecon,UsePrecon,nk,ik,kpoint,d_dense1D,&
                    h_dense1D,s_dense1D,t_dense1D)
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?
  logical, intent(in) :: UpdatePrecon       ! update the preconditioner?
  logical, intent(in) :: UsePrecon          ! use the preconditioner?

  integer, intent(in) :: iscf  ! SCF iteration num.
  integer, intent(in) :: h_dim ! num. of AOs (global)
  integer, intent(in) :: N_occ ! num. of WFs (global)
  integer, intent(in) :: nspin ! num. of spins
  integer, intent(in) :: ispin ! up/down spin
  integer, intent(in) :: nk    ! num. of k points
  integer, intent(in) :: ik    ! k point num.

  real(dp), intent(in) :: eta         ! chemical potential for Kim functional
  real(dp), intent(in) :: kpoint(:,:) ! k point vectors

  complex(dp), intent(in), optional :: h_dense1D(:,:) ! hamiltonian matrix in AO basis (dense 1D)
  complex(dp), intent(in), optional :: s_dense1D(:,:) ! overlap matrix in AO basis (dense 1D)
  complex(dp), intent(in), optional :: t_dense1D(:,:) ! kinetic energy matrix in AO basis (dense 1D)

  !**** OUTPUT **********************************!

  complex(dp), intent(out) :: d_dense1D(:,:) ! (energy-)density matrix in AO basis (dense 1D)

  !**** LOCAL ***********************************!

  character(len=6) :: k_point_name
  character(len=100) :: WF_COEFFS_filename
#ifdef MPI
  character(len=5) :: Node_name
#endif

  logical :: new_s
  logical :: conv
  logical :: ls_conv
  logical :: ls_fail
  logical :: ReadCoeffs2

  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: m
  integer :: n
  integer :: info
  integer :: icg                           ! CG step num.
  integer :: n_step_max=100                ! max. num. steps for CG minimization
  integer :: lwork
  integer, allocatable :: ipiv(:)
#ifdef MPI
  integer :: liwork
  integer :: mpi_status(1:mpi_status_size) ! MPI status
  integer, save :: ictxt                   ! handle for main BLACS context (1D or 2D)
  integer, save :: ictxt_1D                ! handle for additional BLACS context (1D)
  integer, save :: ictxt_1D_T              ! handle for additional BLACS context (1D transposed)
  integer, save :: desc1_1D(1:9)           ! descriptor for operator matrix in AO basis (1D)
  integer, save :: desc3_1D_T(1:9,1:2)     ! descriptor for WF coeffs. matrix (1D transposed)
  integer, allocatable :: iwork(:)
  integer, external :: numroc
#else
  integer, external :: ilaenv
#endif

  real(dp) :: rn
  real(dp) :: rn2
  real(dp) :: rn3
  real(dp) :: rn4
  real(dp) :: E_diff
  real(dp) :: TrQS
  real(dp) :: lambda
  real(dp) :: lambda_n
  real(dp) :: lambda_d
  real(dp) :: lambda_n_tot
  real(dp) :: lambda_d_tot
  real(dp) :: E_OMM                         ! OMM functional energy
  real(dp) :: E_OMM_old                     ! OMM functional energy at previous step
  real(dp) :: coeff(0:4)                    ! coeffs. of the quartic equation
  real(dp), save :: t_precon_scale          ! kinetic energy scale for the preconditioning
  real(dp), save :: cg_tol                  ! convergence tolerance of CG minimization
  real(dp), allocatable :: dwork1(:,:)      ! work matrix
  real(dp), allocatable, save :: x_min(:,:) ! position of minimum
#ifdef MPI
  real(dp) :: Tr_loc
  real(dp), save :: scale_Cholesky          ! Cholesky factorization eigenvalue scaling factor
  real(dp), allocatable :: dwork2(:,:)      ! work matrix
#endif

  complex(dp), allocatable :: work1(:,:)     ! work matrix
  complex(dp), allocatable :: Sd(:,:)        ! g^T*s*g
  complex(dp), allocatable :: Sdd(:,:)       ! g^T*s*g
  complex(dp), allocatable :: g(:,:)         ! gradient
  complex(dp), allocatable :: g_p(:,:)       ! gradient at previous step
  complex(dp), allocatable :: pg(:,:)        ! preconditioned gradient
  complex(dp), allocatable :: pg_p(:,:)      ! preconditioned gradient at previous step
  complex(dp), allocatable :: d(:,:)         ! conjugate search direction
  complex(dp), allocatable :: hc(:,:)        ! h*c
  complex(dp), allocatable :: hg(:,:)        ! h*g
  complex(dp), allocatable :: h_dense(:,:)   ! hamiltonian matrix in AO basis (dense 2D or Cholesky transformed)
#ifdef MPI
  complex(dp), allocatable :: d_dense2D(:,:) ! (energy-)density matrix in AO basis (dense 2D)
  complex(dp), allocatable :: work2(:,:)     ! work matrix
  complex(dp), allocatable :: work3(:,:)     ! work matrix
#endif

  type(multispin), save :: twoI(1:2)               ! identity matrix x2
  type(multispin), allocatable, save :: s_dense(:) ! overlap matrix in AO basis (dense 1D or 2D)
  type(multispin), allocatable, save :: p_dense(:) ! preconditioning matrix in AO basis (dense 1D or 2D)
  type(multispin), allocatable, save :: H(:,:)     ! hamiltonian matrix in WF basis
  type(multispin), allocatable, save :: S(:,:)     ! overlap matrix in WF basis
  type(multispin), allocatable, save :: Hd(:,:)    ! g^T*h*c
  type(multispin), allocatable, save :: Hdd(:,:)   ! g^T*h*g
  type(multispin), allocatable, save :: sc(:,:)    ! s*c
  type(multispin), allocatable, save :: sg(:,:)    ! s*g
  type(multispin), allocatable, save :: cd(:,:)    ! work matrix

  !**********************************************!

  call timer('m_cg',1)

  if (.not. allocated(x_min)) allocate(x_min(1:nk,1:2))
  if (.not. allocated(s_dense)) allocate(s_dense(1:nk))
  if (.not. allocated(p_dense)) allocate(p_dense(1:nk))
  if (.not. allocated(H)) allocate(H(1:nk,1:2))
  if (.not. allocated(S)) allocate(S(1:nk,1:2))
  if (.not. allocated(Hd)) allocate(Hd(1:nk,1:2))
  if (.not. allocated(Hdd)) allocate(Hdd(1:nk,1:2))
  if (.not. allocated(sc)) allocate(sc(1:nk,1:2))
  if (.not. allocated(sg)) allocate(sg(1:nk,1:2))
  if (.not. allocated(cd)) allocate(cd(1:nk,1:2))

  ! if this is the first time the minimization module is called, several things need to be done
  ! (detailed below)  
  if (FirstCall(ik,ispin)) then

#ifdef MPI
    if (all(FirstCall(1:nk,1:2))) then

      ! initialize the BLACS process grids
      if (Use2D) then
        call blacs_get(-1,0,ictxt_1D)
        call blacs_gridinit(ictxt_1D,'C',1,Nodes)
        call blacs_get(ictxt_1D,10,ictxt)
        call blacs_gridinit(ictxt,'C',processorY,Nodes/processorY)
        call blacs_get(ictxt,10,ictxt_1D_T)
        call blacs_gridinit(ictxt_1D_T,'C',Nodes,1)
      else
        call blacs_get(-1,0,ictxt)
        call blacs_gridinit(ictxt,'C',1,Nodes)
        call blacs_get(ictxt,10,ictxt_1D_T)
        call blacs_gridinit(ictxt_1D_T,'C',Nodes,1)
      end if

    end if

    if (all(FirstCall(1:nk,ispin))) then

      ! calculate the local dimensions of the AO and WF matrices
      call blacs_gridinfo(ictxt_1D_T,i,j,k,l)
      N_occ_loc_1D(ispin)=numroc(N_occ,BlockSize_c,k,0,Nodes)
      if (Use2D) then
        call blacs_gridinfo(ictxt,i,j,k,l)
        h_dim_loc(1)=numroc(h_dim,BlockSize,k,0,processorY)
        h_dim_loc(2)=numroc(h_dim,BlockSize,l,0,Nodes/processorY)
        N_occ_loc(1,ispin)=numroc(N_occ,BlockSize_c,k,0,processorY)
        N_occ_loc(2,ispin)=numroc(N_occ,BlockSize_c,l,0,Nodes/processorY)
      else
        h_dim_loc(1)=h_dim
        h_dim_loc(2)=h_dim_loc_1D
        N_occ_loc(1,ispin)=N_occ
        N_occ_loc(2,ispin)=N_occ_loc_1D(ispin)
      end if

      ! initialize the matrix descriptors
      if (Use2D) then
        call descinit(desc1_1D,h_dim,h_dim,BlockSize,BlockSize,0,0,ictxt_1D,h_dim,info)
        if (info/=0) call die('ERROR: desc1_1D setup has failed in minim!')
      end if
      call descinit(desc1,h_dim,h_dim,BlockSize,BlockSize,0,0,ictxt,h_dim_loc(1),info)
      if (info/=0) call die('ERROR: desc1 setup has failed in minim!')
      call descinit(desc2(1:9,ispin),N_occ,N_occ,BlockSize_c,BlockSize_c,0,0,ictxt,N_occ_loc(1,ispin),info)
      if (info/=0) call die('ERROR: desc2 setup has failed in minim!')
      call descinit(desc3(1:9,ispin),N_occ,h_dim,BlockSize_c,BlockSize,0,0,ictxt,N_occ_loc(1,ispin),info)
      if (info/=0) call die('ERROR: desc3 setup has failed in minim!')
      call descinit(desc3_1D_T(1:9,ispin),N_occ,h_dim,BlockSize_c,BlockSize,0,0,ictxt_1D_T,N_occ_loc_1D(ispin),info)
      if (info/=0) call die('ERROR: desc3_1D_T setup has failed in minim!')

    end if
#else
      N_occ_loc_1D(ispin)=N_occ
      h_dim_loc(1)=h_dim
      h_dim_loc(2)=h_dim
      N_occ_loc(1,ispin)=N_occ
      N_occ_loc(2,ispin)=N_occ
#endif

    ! allocate the WF coeffs. matrix
    allocate(c(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    if (UseCholesky) allocate(c_orig(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))

    ! if this is the first SCF step, then we need to initialize the WF coeffs. matrix with random
    ! numbers between -0.5 and 0.5 (normalize at the end to avoid instabilities), unless we are
    ! reading them from file
    if (ReadCoeffs) then
      write(k_point_name,'(i6)') ik
#ifdef MPI
      write(Node_name,'(i5)') Node
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        end if
      end if
#else
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))
        end if
      end if
#endif
      inquire(file=trim(WF_COEFFS_filename),exist=ReadCoeffs2)
    else
      ReadCoeffs2=.false.
    end if
    if (ReadCoeffs2) then
      call io_assign(i)
      open(i,file=trim(WF_COEFFS_filename),form='unformatted',status='old',action='read')
      if (UseCholesky) then
        read(i) c_orig(ik,ispin)%mtrx
      else
        read(i) c(ik,ispin)%mtrx
      end if
      call io_close(i)
    else
      if (ik==1) then
        if ((ispin==1) .or. (N_occ_diff/=0)) then
          call rand_init
          do i=1,h_dim_loc(2)
            do j=1,N_occ_loc(1,ispin)
              call random_number(rn)
              call random_number(rn2)
              call random_number(rn3)
              call random_number(rn4)
              c(ik,ispin)%mtrx(j,i)=cmplx(sign(0.5_dp*rn,rn2-0.5_dp),sign(0.5_dp*rn3,rn4-0.5_dp),dp)
            end do
          end do
          c(ik,ispin)%mtrx=1.0d-2*c(ik,ispin)%mtrx/sqrt(real(h_dim,dp))
        else
          c(ik,2)%mtrx=c(ik,1)%mtrx
        end if
      else
        rn2=999999.9_dp
        do i=1,ik-1
          rn=(kpoint(1,ik)-kpoint(1,i))**2+&
             (kpoint(2,ik)-kpoint(2,i))**2+&
             (kpoint(3,ik)-kpoint(3,i))**2
          if (rn<rn2) then
            rn2=rn
            j=i
          end if
        end do
        c(ik,ispin)%mtrx=c(j,ispin)%mtrx
      end if
    end if

    if (all(FirstCall(1:nk,1:2))) then
      t_precon_scale=fdf_physical('OMM.TPreconScale',10.0_dp,'Ry')
      cg_tol=fdf_get('OMM.RelTol',1.0d-9)
    end if

  end if

#ifdef MPI
  if (Use2D) allocate(d_dense2D(1:h_dim_loc(1),1:h_dim_loc(2)))
#endif

  if (CalcE) then

    ! calculate the energy-density matrix: e=c^T*[(2*I-S)*(H+eta*S)]*c
#ifdef MPI
    call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Hd(ik,ispin)%mtrx,1,1,desc2(1:9,ispin),cmplx_1,H(ik,ispin)%mtrx,&
                 1,1,desc2(1:9,ispin))
#else
    do k=1,N_occ
      do l=1,N_occ
        H(ik,ispin)%mtrx(k,l)=H(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Hd(ik,ispin)%mtrx(l,k))
      end do
    end do
#endif
    H(ik,ispin)%mtrx=H(ik,ispin)%mtrx+x_min(ik,ispin)*Hd(ik,ispin)%mtrx+x_min(ik,ispin)**2*Hdd(ik,ispin)%mtrx
    allocate(work1(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    if (UseCholesky) then
#ifdef MPI
      if (Use2D) then
        call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense2D,&
                          work1,cd(ik,ispin)%mtrx)
      else
        call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense1D,&
                          work1,cd(ik,ispin)%mtrx)
      end if
#else
      call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense1D,&
                        work1,cd(ik,ispin)%mtrx)
#endif
    else
#ifdef MPI
      if (Use2D) then
        call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense2D,&
                          work1,cd(ik,ispin)%mtrx)
      else
        call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense1D,&
                          work1,cd(ik,ispin)%mtrx)
      end if
#else
      call calc_densmat(h_dim,N_occ,ispin,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense1D,work1,&
                        cd(ik,ispin)%mtrx)
#endif
    end if
#ifdef MPI
    if (Use2D) call pzgemr2d(h_dim,h_dim,d_dense2D,1,1,desc1,d_dense1D,1,1,desc1_1D,ictxt)
#endif
    if (nspin==1) d_dense1D=cmplx_2*d_dense1D

    if (WriteCoeffs) then
      call io_assign(i)
      write(k_point_name,'(i6)') ik
#ifdef MPI
      write(Node_name,'(i5)') Node
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        end if
      end if
#else
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))
        end if
      end if
#endif
      open(i,file=trim(WF_COEFFS_filename),form='unformatted',status='replace',action='write')
      if (UseCholesky) then
        write(i) c_orig(ik,ispin)%mtrx
      else
        write(i) c(ik,ispin)%mtrx
      end if
      call io_close(i)
    end if

    deallocate(work1)
    deallocate(cd(ik,ispin)%mtrx)
    if (allocated(p_dense(ik)%mtrx)) deallocate(p_dense(ik)%mtrx)
    if (.not. UseCholesky) then
      deallocate(sg(ik,ispin)%mtrx)
      deallocate(sc(ik,ispin)%mtrx)
    end if
    if (allocated(s_dense(ik)%mtrx)) deallocate(s_dense(ik)%mtrx)
    deallocate(S(ik,ispin)%mtrx)
    deallocate(Hdd(ik,ispin)%mtrx)
    deallocate(Hd(ik,ispin)%mtrx)
    deallocate(H(ik,ispin)%mtrx)
#ifdef MPI
    if (Use2D) deallocate(d_dense2D)
#endif

    call timer('m_cg',2)

    return

  end if

  if (.not. allocated(H(ik,ispin)%mtrx)) allocate(H(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(Hd(ik,ispin)%mtrx)) allocate(Hd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(Hdd(ik,ispin)%mtrx)) allocate(Hdd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(S(ik,ispin)%mtrx)) then
    allocate(S(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
    if (.not. allocated(s_dense(ik)%mtrx)) then
      allocate(s_dense(ik)%mtrx(1:h_dim_loc(1),1:h_dim_loc(2)))
#ifdef MPI
      if (Use2D) then
        call pzgemr2d(h_dim,h_dim,s_dense1D,1,1,desc1_1D,s_dense(ik)%mtrx,1,1,desc1,ictxt)
      else
        s_dense(ik)%mtrx=s_dense1D
      end if
#else
      s_dense(ik)%mtrx=s_dense1D
#endif
      if (UseCholesky) then
#ifdef MPI
        call pzpotrf('U',h_dim,s_dense(ik)%mtrx,1,1,desc1,info)
        if (info/=0) call die('ERROR: pzpotrf has failed in minim!')
#else
        call zpotrf('U',h_dim,s_dense(ik)%mtrx,h_dim,info)
        if (info/=0) call die('ERROR: zpotrf has failed in minim!')
#endif
      end if
    end if
    new_s=.true.
  else
    new_s=.false.
  end if
  if (.not. UseCholesky) then
    if (.not. allocated(sc(ik,ispin)%mtrx)) allocate(sc(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    if (.not. allocated(sg(ik,ispin)%mtrx)) allocate(sg(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  end if

  if (UseCholesky) then
    if ((new_s .and. (.not. FirstCall(ik,ispin))) .or. &
        PreviousCallDiagon .or. &
        (ReadCoeffs .and. FirstCall(ik,ispin))) then
      if (.not. PreviousCallDiagon) c(ik,ispin)%mtrx=c_orig(ik,ispin)%mtrx
#ifdef MPI
      call pztrmm('R','U','C','N',N_occ,h_dim,cmplx_1,s_dense(ik)%mtrx,1,1,desc1,c(ik,ispin)%mtrx,1,1,desc3(1:9,ispin))
#else
      call ztrmm('R','U','C','N',N_occ,h_dim,cmplx_1,s_dense(ik)%mtrx,h_dim,c(ik,ispin)%mtrx,N_occ)
#endif
    end if
  end if

#ifdef MPI
  if (UseCholesky .or. Use2D) then
    allocate(h_dense(1:h_dim_loc(1),1:h_dim_loc(2)))
    if (Use2D) then
      call pzgemr2d(h_dim,h_dim,h_dense1D,1,1,desc1_1D,h_dense,1,1,desc1,ictxt)
    else
      h_dense=h_dense1D
    end if
  end if
#else
  if (UseCholesky) then
    allocate(h_dense(1:h_dim_loc(1),1:h_dim_loc(2)))
    h_dense=h_dense1D
  end if
#endif
  if (UseCholesky) then
#ifdef MPI
    call pzhegst(1,'U',h_dim,h_dense,1,1,desc1,s_dense(ik)%mtrx,1,1,desc1,scale_Cholesky,info)
    if (info/=0) call die('ERROR: pzhegst has failed in minim!')
    allocate(work1(1:h_dim_loc(1),1:h_dim_loc(2)))
    allocate(work2(1:h_dim_loc(1),1:h_dim_loc(2)))
    allocate(work3(1:h_dim_loc(1),1:h_dim_loc(2)))
    call pztranc(h_dim,h_dim,cmplx_1,h_dense,1,1,desc1,cmplx_0,work1,1,1,desc1)
    work2=cmplx_0
    call pzlaset('U',h_dim,h_dim,cmplx_1,cmplx_half,work2,1,1,desc1)
    work3=cmplx_0
    call pzlaset('L',h_dim,h_dim,cmplx_1,cmplx_half,work3,1,1,desc1)
    h_dense=work2*h_dense+work3*work1
    deallocate(work3)
    deallocate(work2)
    deallocate(work1)
#else
    call zhegst(1,'U',h_dim,h_dense,h_dim,s_dense(ik)%mtrx,h_dim,info)
    if (info/=0) call die('ERROR: zhegst has failed in minim!')
    do i=1,h_dim-1
      do j=i+1,h_dim
        h_dense(j,i)=conjg(h_dense(i,j))
      end do
    end do
#endif
  end if

  ! calculate the preconditioning matrix (s+t/tau)^-1
  if (UpdatePrecon .and. (ispin==1)) then
    if (.not. allocated(p_dense(ik)%mtrx)) allocate(p_dense(ik)%mtrx(1:h_dim_loc(1),1:h_dim_loc(2)))
#ifdef MPI
    if (Use2D) then
      call pzgemr2d(h_dim,h_dim,t_dense1D,1,1,desc1_1D,p_dense(ik)%mtrx,1,1,desc1,ictxt)
    else
      p_dense(ik)%mtrx=t_dense1D
    end if
#else
    p_dense(ik)%mtrx=t_dense1D
#endif
    p_dense(ik)%mtrx=s_dense(ik)%mtrx+p_dense(ik)%mtrx/t_precon_scale
#ifdef MPI
    allocate(ipiv(1:h_dim_loc(1)+BlockSize))
    call pzgetrf(h_dim,h_dim,p_dense(ik)%mtrx,1,1,desc1,ipiv,info)
    if (info/=0) call die('ERROR: pzgetrf has failed in minim!')
    allocate(work1(1:1,1:1))
    allocate(iwork(1:1))
    call pzgetri(h_dim,p_dense(ik)%mtrx,1,1,desc1,ipiv,work1,-1,iwork,-1,info)
    if (info/=0) call die('ERROR: pzgetri has failed in minim!')
    liwork=iwork(1)
    deallocate(iwork)
    lwork=work1(1,1)
    deallocate(work1)
    allocate(work1(1:lwork,1:1))
    allocate(iwork(1:liwork))
    call pzgetri(h_dim,p_dense(ik)%mtrx,1,1,desc1,ipiv,work1,lwork,iwork,liwork,info)
    if (info/=0) call die('ERROR: pzgetri has failed in minim!')
    deallocate(iwork)
    deallocate(work1)
    deallocate(ipiv)
#else
    allocate(ipiv(1:h_dim))
    lwork=h_dim*ilaenv(1,'zhetrf','U',h_dim,-1,-1,-1)
    allocate(work1(1:lwork,1:1))
    call zhetrf('U',h_dim,p_dense(ik)%mtrx,h_dim,ipiv,work1,lwork,info)
    if (info/=0) call die('ERROR: zhetrf has failed in minim!')
    deallocate(work1)
    allocate(work1(1:h_dim,1:1))
    call zhetri('U',h_dim,p_dense(ik)%mtrx,h_dim,ipiv,work1,info)
    if (info/=0) call die('ERROR: zhetri has failed in minim!')
    deallocate(work1)
    deallocate(ipiv)
    do i=1,h_dim-1
      do j=i+1,h_dim
        p_dense(ik)%mtrx(j,i)=conjg(p_dense(ik)%mtrx(i,j))
      end do
    end do
#endif
  end if

  allocate(Sd(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  allocate(Sdd(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  allocate(g(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(g_p(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  if (UsePrecon) then
    allocate(pg(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    allocate(pg_p(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  end if
  allocate(d(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(hc(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(hg(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(work1(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))

  ! first we calculate the energy and gradient for our initial guess, with the following steps:
  ! -calculate the hamiltonian in WF basis: H=c^T*h*c
  if (allocated(h_dense)) then
    call calc_A(h_dim,N_occ,ispin,h_dense,c(ik,ispin)%mtrx,H(ik,ispin)%mtrx,hc)
  else
    call calc_A(h_dim,N_occ,ispin,h_dense1D,c(ik,ispin)%mtrx,H(ik,ispin)%mtrx,hc)
  end if
  ! -calculate the overlap matrix in WF basis: S=c^T*s*c
  if (UseCholesky) then
#ifdef MPI
    if (new_s .or. PreviousCallDiagon) call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),&
                                                   c(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),cmplx_0,S(ik,ispin)%mtrx,1,1,&
                                                   desc2(1:9,ispin))
#else
    if (new_s .or. PreviousCallDiagon) call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,N_occ,c(ik,ispin)%mtrx,N_occ,&
                                                  cmplx_0,S(ik,ispin)%mtrx,N_occ)
#endif
  else
    if (new_s .or. PreviousCallDiagon) then
      call calc_A(h_dim,N_occ,ispin,s_dense(ik)%mtrx,c(ik,ispin)%mtrx,S(ik,ispin)%mtrx,sc(ik,ispin)%mtrx)
    else
      sc(ik,ispin)%mtrx=sc(ik,ispin)%mtrx+x_min(ik,ispin)*sg(ik,ispin)%mtrx
    end if
  end if
  ! -calculate the gradient: g=2*(2*h*c-s*c*H-h*c*S)
  !  (note that we *reuse* h*c and s*c contained in hc and sc from the previous call to calc_A)
  if (UseCholesky) then
    call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,c(ik,ispin)%mtrx)
  else
    call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,sc(ik,ispin)%mtrx)
  end if
  ! -calculate the preconditioned gradient by premultiplying g by (s+t/tau)^-1
  if (UsePrecon) then
#ifdef MPI
    call pzgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,1,1,desc3(1:9,ispin),p_dense(ik)%mtrx,1,1,desc1,cmplx_0,pg,1,1,desc3(1:9,ispin))
#else
    call zgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,N_occ,p_dense(ik)%mtrx,h_dim,cmplx_0,pg,N_occ)
#endif
  end if
  ! -calculate the additional matrices:
  !  Hd=g^T*h*c
  !  Sd=g^T*s*c
  !  Hdd=g^T*h*g
  !  Sdd=g^T*s*g
  !  (again, h*c has already been calculated, although h*g has not)
  !  and, finally, the coeffs. of the quartic line search equation in the direction g
  !  (the energy at c is given by the zeroth-order coeff. c(0))
  if (UsePrecon) then
#ifdef MPI
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),pg,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                desc2(1:9,ispin))
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),pg,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                desc2(1:9,ispin))
#else
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,pg,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,pg,N_occ,cmplx_0,Sd,N_occ)
#endif
  else
#ifdef MPI
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                desc2(1:9,ispin))
#else
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,g,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
#endif
    if (UseCholesky) then
#ifdef MPI
      call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                  desc2(1:9,ispin))
#else
      call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,N_occ,g,N_occ,cmplx_0,Sd,N_occ)
#endif
    else
#ifdef MPI
      call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                  desc2(1:9,ispin))
#else
      call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,g,N_occ,cmplx_0,Sd,N_occ)
#endif
    end if
  end if
  if (UsePrecon) then
    if (allocated(h_dense)) then
      call calc_A(h_dim,N_occ,ispin,h_dense,pg,Hdd(ik,ispin)%mtrx,hg)
    else
      call calc_A(h_dim,N_occ,ispin,h_dense1D,pg,Hdd(ik,ispin)%mtrx,hg)
    end if
    call calc_A(h_dim,N_occ,ispin,s_dense(ik)%mtrx,pg,Sdd,sg(ik,ispin)%mtrx)
  else
    if (allocated(h_dense)) then
      call calc_A(h_dim,N_occ,ispin,h_dense,g,Hdd(ik,ispin)%mtrx,hg)
    else
      call calc_A(h_dim,N_occ,ispin,h_dense1D,g,Hdd(ik,ispin)%mtrx,hg)
    end if
    if (UseCholesky) then
#ifdef MPI
      call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,g,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Sdd,1,1,desc2(1:9,ispin))
#else
      call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,g,N_occ,g,N_occ,cmplx_0,Sdd,N_occ)
#endif
    else
      call calc_A(h_dim,N_occ,ispin,s_dense(ik)%mtrx,g,Sdd,sg(ik,ispin)%mtrx)
    end if
  end if
  call calc_coeff(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,Hd(ik,ispin)%mtrx,Sd,Hdd(ik,ispin)%mtrx,Sdd,coeff,work1)
  E_OMM=coeff(0)

  ! this is the main loop of the CG algorithm. We perform a series of line minimizations, with the
  ! gradient g at each new step being modified to obtain the search direction d
  if (Node==0) then
    if ((ik==1) .and. (ispin==1)) then
      print'(a)', '+---------------------------------------------+'
      if (UseCholesky) then
        print'(a)', '| OMM (Cholesky factorization)                |'
      else if (UsePrecon) then
        print'(a)', '| OMM (preconditioning)                       |'
      else
        print'(a)', '| OMM                                         |'
      end if
      print'(a)',      '+---------------------------------------------+'
    end if
    if (nk>1) then
      print'(a,i6,a,i6,a)', '| k point ', ik, ' of ', nk, '                    |'
    end if
    if (nspin==2) then
      if (ispin==1) then
        print'(a)', '| up spin                                     |'
      else
        print'(a)', '| down spin                                   |'
      end if
    end if
    if ((nk>1) .or. (nspin==2)) print'(a)', '+---------------------------------------------+'
    if (LongOut) print'(a)', '|             E_OMM            E_diff         |'
  end if
  conv=.false.
  d=cmplx_0
  icg=0
  do i=1,n_step_max
    lambda=0.0_dp
    do j=1,2*h_dim*N_occ-1
      if (UsePrecon) then
        d=pg+cmplx(lambda,0.0_dp,dp)*d
      else
        d=g+cmplx(lambda,0.0_dp,dp)*d
      end if
      g_p=g
      if (UsePrecon) pg_p=pg
      E_OMM_old=E_OMM
      ! if this is not the first CG step, we have to recalculate Hd, Sd, Hdd, Sdd, and the coeffs.
      if (icg>0) then
#ifdef MPI
        call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                    desc2(1:9,ispin))
#else
        call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,d,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
#endif
        if (UseCholesky) then
#ifdef MPI
          call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,&
                      1,desc2(1:9,ispin))
#else
          call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,c(ik,ispin)%mtrx,N_occ,d,N_occ,cmplx_0,Sd,N_occ)
#endif
        else
#ifdef MPI
          call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,&
                      1,desc2(1:9,ispin))
#else
          call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,d,N_occ,cmplx_0,Sd,N_occ)
#endif
        end if
        if (allocated(h_dense)) then
          call calc_A(h_dim,N_occ,ispin,h_dense,d,Hdd(ik,ispin)%mtrx,hg)
        else
          call calc_A(h_dim,N_occ,ispin,h_dense1D,d,Hdd(ik,ispin)%mtrx,hg)
        end if
        if (UseCholesky) then
#ifdef MPI
          call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,d,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Sdd,1,1,&
                      desc2(1:9,ispin))
#else
          call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,d,N_occ,d,N_occ,cmplx_0,Sdd,N_occ)
#endif
        else
          call calc_A(h_dim,N_occ,ispin,s_dense(ik)%mtrx,d,Sdd,sg(ik,ispin)%mtrx)
        end if
        call calc_coeff(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,Hd(ik,ispin)%mtrx,Sd,Hdd(ik,ispin)%mtrx,Sdd,coeff,work1)
      end if
      ! using the coeffs. calculated anlytically, we can find the minimum of the functional in the
      ! search direction, and calculate the energy at that minimum
      call solve_quartic(coeff(0:4),x_min(ik,ispin),ls_fail)
      ! in certain regions of the coeffs. space the line search gives no minimum--this occurs when there
      ! are positive eigenvalues in the eigenspecturm which are significantly occupied by our coeffs.
      ! matrix; the only known cure, unfortunately, is to scale down the entire matrix, thus returning to
      ! a  safe region of the coeffs. space.
      if (ls_fail) then
        if (Node==0) print'(a)', '| WARNING: Rescaling coefficients!            |'
        E_OMM=3.0*E_OMM
        c(ik,ispin)%mtrx=cmplx(0.5_dp,0.0_dp,dp)*c(ik,ispin)%mtrx
        ls_conv=.false.
      else
        ! if the line search is successful, move to the minimum
        E_OMM=coeff(4)*x_min(ik,ispin)**4+&
              coeff(3)*x_min(ik,ispin)**3+&
              coeff(2)*x_min(ik,ispin)**2+&
              coeff(1)*x_min(ik,ispin)+&
              coeff(0)
        c(ik,ispin)%mtrx=c(ik,ispin)%mtrx+x_min(ik,ispin)*d
        ls_conv=.true.
      end if
      ! recalculate S at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        S(ik,ispin)%mtrx=cmplx(0.25_dp,0.0_dp,dp)*S(ik,ispin)%mtrx
      else
#ifdef MPI
        call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Sd,1,1,desc2(1:9,ispin),cmplx_1,S(ik,ispin)%mtrx,1,1,&
                     desc2(1:9,ispin))
#else
        do k=1,N_occ
          do l=1,N_occ
            S(ik,ispin)%mtrx(k,l)=S(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Sd(l,k))
          end do
        end do
#endif
        S(ik,ispin)%mtrx=S(ik,ispin)%mtrx+x_min(ik,ispin)*Sd+x_min(ik,ispin)**2*Sdd
      end if
      E_diff=2.0_dp*abs((E_OMM-E_OMM_old)/(E_OMM+E_OMM_old))
      if ((Node==0) .and. LongOut) print'(a,2(1x,i5),2(1x,es15.7e3),1x,a)', '|', i, j, E_OMM, E_diff, '|'
      icg=icg+1
      if (E_diff<=cg_tol) then
        conv=.true.
        exit
      end if
      ! recalculate H at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        H(ik,ispin)%mtrx=cmplx(0.25_dp,0.0_dp,dp)*H(ik,ispin)%mtrx
      else
#ifdef MPI
        call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Hd(ik,ispin)%mtrx,1,1,desc2(1:9,ispin),cmplx_1,&
                     H(ik,ispin)%mtrx,1,1,desc2(1:9,ispin))
#else
        do k=1,N_occ
          do l=1,N_occ
            H(ik,ispin)%mtrx(k,l)=H(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Hd(ik,ispin)%mtrx(l,k))
          end do
        end do
#endif
        H(ik,ispin)%mtrx=H(ik,ispin)%mtrx+x_min(ik,ispin)*Hd(ik,ispin)%mtrx+x_min(ik,ispin)**2*Hdd(ik,ispin)%mtrx
      end if
      ! recalculate g at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        hc=cmplx(0.5_dp,0.0_dp,dp)*hc
        if (.not. UseCholesky) sc(ik,ispin)%mtrx=cmplx(0.5_dp,0.0_dp,dp)*sc(ik,ispin)%mtrx
        g=g_p+cmplx(1.5_dp,0.0_dp,dp)*hc
      else
        hc=hc+x_min(ik,ispin)*hg
        if (UseCholesky) then
          call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,c(ik,ispin)%mtrx)
        else
          sc(ik,ispin)%mtrx=sc(ik,ispin)%mtrx+x_min(ik,ispin)*sg(ik,ispin)%mtrx
          call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,sc(ik,ispin)%mtrx)
        end if
      end if
      if (UsePrecon) then
#ifdef MPI
        call pzgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,1,1,desc3(1:9,ispin),p_dense(ik)%mtrx,1,1,desc1,cmplx_0,pg,1,1,&
                    desc3(1:9,ispin))
#else
        call zgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,N_occ,p_dense(ik)%mtrx,h_dim,cmplx_0,pg,N_occ)
#endif
      end if
      if (ls_conv) then
        if (UsePrecon) then
          lambda_n=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_n=lambda_n+real(pg(k,l),dp)*(real(g(k,l),dp)-real(g_p(k,l),dp))
              lambda_n=lambda_n+aimag(pg(k,l))*(aimag(g(k,l))-aimag(g_p(k,l)))
            end do
          end do
          lambda_d=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_d=lambda_d+real(pg_p(k,l),dp)*real(g_p(k,l),dp)
              lambda_d=lambda_d+aimag(pg_p(k,l))*aimag(g_p(k,l))
            end do
          end do
        else
          lambda_n=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_n=lambda_n+real(g(k,l),dp)*(real(g(k,l),dp)-real(g_p(k,l),dp))
              lambda_n=lambda_n+aimag(g(k,l))*(aimag(g(k,l))-aimag(g_p(k,l)))
            end do
          end do
          lambda_d=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_d=lambda_d+real(g_p(k,l),dp)*real(g_p(k,l),dp)
              lambda_d=lambda_d+aimag(g_p(k,l))*aimag(g_p(k,l))
            end do
          end do
        end if
#ifdef MPI
        call mpi_allreduce(lambda_n,lambda_n_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        call mpi_allreduce(lambda_d,lambda_d_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        lambda=lambda_n_tot/lambda_d_tot
#else
        lambda=lambda_n/lambda_d
#endif
      else
        exit
      end if
    end do
    if (conv) exit
  end do
  if (i>n_step_max) then
    if (Node==0) print'(a)', '| WARNING: OMM failed to converge!            |'
  end if
  if ((Node==0) .and. LongOut) print'(a)', '+---------------------------------------------+'

  deallocate(work1)
  deallocate(hg)
  deallocate(hc)
  deallocate(d)
  if (UsePrecon) then
    deallocate(pg_p)
    deallocate(pg)
  end if
  deallocate(g_p)
  deallocate(g)
  deallocate(Sdd)
  deallocate(Sd)
  if (allocated(h_dense)) deallocate(h_dense)

  ! calculate the density matrix: d=c*(2*I-S)*c^T
  if (.not. allocated(twoI(ispin)%mtrx)) then
    allocate(twoI(ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
#ifdef MPI
    call pzlaset('A',N_occ,N_occ,cmplx_0,cmplx_2,twoI(ispin)%mtrx,1,1,desc2(1:9,ispin))
#else
    twoI(ispin)%mtrx=cmplx_0
    do i=1,N_occ
      twoI(ispin)%mtrx(i,i)=cmplx_2
    end do
#endif
  end if
  if (.not. allocated(cd(ik,ispin)%mtrx)) allocate(cd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  if (UseCholesky) then
    c_orig(ik,ispin)%mtrx=c(ik,ispin)%mtrx
#ifdef MPI
    call pztrsm('R','U','C','N',N_occ,h_dim,cmplx_1,s_dense(ik)%mtrx,1,1,desc1,c_orig(ik,ispin)%mtrx,1,1,desc3(1:9,ispin))
    if (Use2D) then
      call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense2D,cd(ik,ispin)%mtrx)
    else
      call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense1D,cd(ik,ispin)%mtrx)
    end if
#else
    call ztrsm('R','U','C','N',N_occ,h_dim,cmplx_1,s_dense(ik)%mtrx,h_dim,c_orig(ik,ispin)%mtrx,N_occ)
    call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c_orig(ik,ispin)%mtrx,d_dense1D,cd(ik,ispin)%mtrx)
#endif
  else
#ifdef MPI
    if (Use2D) then
      call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense2D,cd(ik,ispin)%mtrx)
    else
      call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense1D,cd(ik,ispin)%mtrx)
    end if
#else
    call calc_densmat(h_dim,N_occ,ispin,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,d_dense1D,cd(ik,ispin)%mtrx)
#endif
  end if
#ifdef MPI
  if (Use2D) then
    call pzgemr2d(h_dim,h_dim,d_dense2D,1,1,desc1,d_dense1D,1,1,desc1_1D,ictxt)
    deallocate(d_dense2D)
  end if
#endif
  if (nspin==1) d_dense1D=cmplx_2*d_dense1D

  ! calculate the trace of S to make sure we are occupying the right number of eigenstates in our
  ! solution
#ifdef MPI
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+(real(twoI(ispin)%mtrx(i,j),dp)-real(S(ik,ispin)%mtrx(i,j),dp))*real(S(ik,ispin)%mtrx(i,j),dp)+&
                    (aimag(twoI(ispin)%mtrx(i,j))-aimag(S(ik,ispin)%mtrx(i,j)))*aimag(S(ik,ispin)%mtrx(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrQS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  TrQS=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrQS=TrQS+(real(twoI(ispin)%mtrx(i,j),dp)-real(S(ik,ispin)%mtrx(i,j),dp))*real(S(ik,ispin)%mtrx(i,j),dp)+&
                (aimag(twoI(ispin)%mtrx(i,j))-aimag(S(ik,ispin)%mtrx(i,j)))*aimag(S(ik,ispin)%mtrx(i,j))
    end do
  end do
#endif
  if (Node==0) then
    if (nspin==1) then
      print'(a,i5,a)',    '| minim: icg             = ', icg, '              |'
      print'(a,f13.7,a)', '| minim: 2*Tr[(2*I-S)*S] = ', 2.0_dp*TrQS, '      |'
    else
      print'(a,i5,a)',    '| minim: icg           = ', icg, '                |'
      print'(a,f13.7,a)', '| minim: Tr[(2*I-S)*S] = ', TrQS, '        |'
    end if
    print'(a)',       '+---------------------------------------------+'
  end if

  if (FirstCall(ik,ispin)) FirstCall(ik,ispin)=.false.

  call timer('m_cg',2)

end subroutine minim_cg

!================================================!
! minimize the Kim functional by conjugate       !
! gradients (sparse routine)                     !
!================================================!
subroutine minim_cg_sparse(nhmax,numh,listhptr,listh,CalcE,PreviousCallDiagon,iscf,h_dim,N_occ,eta,nspin,ispin,UpdatePrecon,&
                           UsePrecon,UpdateSparseComm,nk,ik,kpoint,d_sparse,h_sparse,s_sparse,s_dense1D,t_dense1D)
  implicit none

  !**** INPUT ***********************************!

  logical, intent(in) :: CalcE              ! calculate the energy-density matrix from the existing coeffs.?
  logical, intent(in) :: PreviousCallDiagon ! Previous SCF iteration solved by diagonalization?
  logical, intent(in) :: UpdatePrecon       ! update the preconditioner?
  logical, intent(in) :: UsePrecon          ! use the preconditioner?
  logical, intent(in) :: UpdateSparseComm   ! update nhmax_max?

  integer, intent(in) :: nhmax       ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(:)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(:) ! pointer to start of row in listh
  integer, intent(in) :: listh(:)    ! list of nonzero elements of each row of sparse matrices
  integer, intent(in) :: iscf        ! SCF iteration num.
  integer, intent(in) :: h_dim       ! num. of AOs (global)
  integer, intent(in) :: N_occ       ! num. of WFs (global)
  integer, intent(in) :: nspin       ! num. of spins
  integer, intent(in) :: ispin       ! up/down spin
  integer, intent(in) :: nk          ! num. of k points
  integer, intent(in) :: ik          ! k point num.

  real(dp), intent(in) :: eta         ! chemical potential for Kim functional
  real(dp), intent(in) :: kpoint(:,:) ! k point vectors

  complex(dp), intent(in), optional :: h_sparse(:)    ! hamiltonian matrix in AO basis (sparse)
  complex(dp), intent(in), optional :: s_sparse(:)    ! overlap matrix in AO basis (sparse)
  complex(dp), intent(in), optional :: s_dense1D(:,:) ! overlap matrix in AO basis (dense 1D)
  complex(dp), intent(in), optional :: t_dense1D(:,:) ! kinetic energy matrix in AO basis (dense 1D)

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: d_sparse(:) ! (energy-)density matrix in AO basis (sparse)

  !**** LOCAL ***********************************!

  character(len=6) :: k_point_name
  character(len=100) :: WF_COEFFS_filename
#ifdef MPI
  character(len=5) :: Node_name
#endif

  logical :: new_s
  logical :: conv
  logical :: ls_conv
  logical :: ls_fail
  logical :: ReadCoeffs2

  integer :: i
  integer :: j
  integer :: k
  integer :: l
  integer :: m
  integer :: n
  integer :: info
  integer :: icg                           ! CG step num.
  integer :: n_step_max=100                ! max. num. steps for CG minimization
  integer :: lwork
  integer, allocatable :: ipiv(:)
#ifdef MPI
  integer :: liwork
  integer :: mpi_status(1:mpi_status_size) ! MPI status
  integer, save :: ictxt                   ! handle for main BLACS context (1D)
  integer, save :: nhmax_max
  integer, save :: h_dim_loc_max
  integer, allocatable :: iwork(:)
  integer, external :: numroc
#else
  integer, external :: ilaenv
#endif

  real(dp) :: rn
  real(dp) :: rn2
  real(dp) :: rn3
  real(dp) :: rn4
  real(dp) :: E_diff
  real(dp) :: TrQS
  real(dp) :: lambda
  real(dp) :: lambda_n
  real(dp) :: lambda_d
  real(dp) :: lambda_n_tot
  real(dp) :: lambda_d_tot
  real(dp) :: E_OMM                         ! OMM functional energy
  real(dp) :: E_OMM_old                     ! OMM functional energy at previous step
  real(dp) :: coeff(0:4)                    ! coeffs. of the quartic equation
  real(dp), save :: t_precon_scale          ! kinetic energy scale for the preconditioning
  real(dp), save :: cg_tol                  ! convergence tolerance of CG minimization
  real(dp), allocatable :: dwork2(:,:)      ! work matrix
  real(dp), allocatable, save :: x_min(:,:) ! position of minimum
#ifdef MPI
  real(dp) :: Tr_loc
#endif

  complex(dp), allocatable :: work1(:,:)           ! work matrix
  complex(dp), allocatable :: Sd(:,:)              ! g^T*s*g
  complex(dp), allocatable :: Sdd(:,:)             ! g^T*s*g
  complex(dp), allocatable :: g(:,:)               ! gradient
  complex(dp), allocatable :: g_p(:,:)             ! gradient at previous step
  complex(dp), allocatable :: pg(:,:)              ! preconditioned gradient
  complex(dp), allocatable :: pg_p(:,:)            ! preconditioned gradient at previous step
  complex(dp), allocatable :: d(:,:)               ! conjugate search direction
  complex(dp), allocatable :: hc(:,:)              ! h*c
  complex(dp), allocatable :: hg(:,:)              ! h*g
#ifdef MPI
  complex(dp), allocatable :: work2(:,:)           ! work matrix
  complex(dp), allocatable :: work3(:,:)           ! work matrix
#endif

  type(multispin), save :: twoI(1:2)                 ! identity matrix x2
  type(multispin), allocatable, save :: p_dense1D(:) ! preconditioning matrix in AO basis (dense 1D)
  type(multispin), allocatable, save :: H(:,:)       ! hamiltonian matrix in WF basis
  type(multispin), allocatable, save :: S(:,:)       ! overlap matrix in WF basis
  type(multispin), allocatable, save :: Hd(:,:)      ! g^T*h*c
  type(multispin), allocatable, save :: Hdd(:,:)     ! g^T*h*g
  type(multispin), allocatable, save :: sc(:,:)      ! s*c
  type(multispin), allocatable, save :: sg(:,:)      ! s*g
  type(multispin), allocatable, save :: cd(:,:)      ! work matrix

  !**********************************************!

  if (.not. allocated(x_min)) allocate(x_min(1:nk,1:2))
  if (.not. allocated(p_dense1D)) allocate(p_dense1D(1:nk))
  if (.not. allocated(H)) allocate(H(1:nk,1:2))
  if (.not. allocated(S)) allocate(S(1:nk,1:2))
  if (.not. allocated(Hd)) allocate(Hd(1:nk,1:2))
  if (.not. allocated(Hdd)) allocate(Hdd(1:nk,1:2))
  if (.not. allocated(sc)) allocate(sc(1:nk,1:2))
  if (.not. allocated(sg)) allocate(sg(1:nk,1:2))
  if (.not. allocated(cd)) allocate(cd(1:nk,1:2))

  call timer('m_cg',1)

  ! if this is the first time the minimization module is called, several things need to be done
  ! (detailed below)  
  if (FirstCall(ik,ispin)) then

#ifdef MPI
    if (all(FirstCall(1:nk,1:2))) then

      ! calculate the local-to-global index transformations needed for the sparse matrix-matrix
      ! multiplications
      allocate(h_dim_l2g(1:h_dim_loc_1D))
      j=0
      k=0
      l=0
      do i=1,h_dim
        k=k+1
        if (j==Node) then
          l=l+1
          h_dim_l2g(l)=i
        end if
        if (k==BlockSize) then
          k=0
          j=j+1
          if (j==Nodes) j=0
        end if
      end do

      ! find the largest value of h_dim_loc over all MPI processes needed for the sparse matrix
      ! multiplications
      h_dim_loc_max=0
      do i=0,Nodes-1
        if (i==Node) then
          k=h_dim_loc_1D
        end if
        call mpi_bcast(j,1,mpi_integer,i,mpi_comm_world,info)
        call mpi_bcast(k,1,mpi_integer,i,mpi_comm_world,info)
        if (k>h_dim_loc_max) h_dim_loc_max=k
      end do

      ! initialize the BLACS process grids
      call blacs_get(-1,0,ictxt)
      call blacs_gridinit(ictxt,'C',Nodes,1)

    end if

    if (all(FirstCall(1:nk,ispin))) then

      ! calculate the local dimensions of the AO and WF matrices
      call blacs_gridinfo(ictxt,i,j,k,l)
      N_occ_loc_1D(ispin)=numroc(N_occ,BlockSize_c,k,0,Nodes)
      h_dim_loc(1)=h_dim_loc_1D
      h_dim_loc(2)=h_dim
      N_occ_loc(1,ispin)=N_occ_loc_1D(ispin)
      N_occ_loc(2,ispin)=N_occ

      ! initialize the matrix descriptors
      call descinit(desc1,h_dim,h_dim,BlockSize,BlockSize,0,0,ictxt,h_dim_loc(1),info)
      if (info/=0) call die('ERROR: desc1 setup has failed in minim!')
      call descinit(desc2(1:9,ispin),N_occ,N_occ,BlockSize_c,BlockSize_c,0,0,ictxt,N_occ_loc(1,ispin),info)
      if (info/=0) call die('ERROR: desc2 setup has failed in minim!')
      call descinit(desc3(1:9,ispin),N_occ,h_dim,BlockSize_c,BlockSize,0,0,ictxt,N_occ_loc(1,ispin),info)
      if (info/=0) call die('ERROR: desc3 setup has failed in minim!')

    end if
#else
      N_occ_loc_1D(ispin)=N_occ
      h_dim_loc(1)=h_dim
      h_dim_loc(2)=h_dim
      N_occ_loc(1,ispin)=N_occ
      N_occ_loc(2,ispin)=N_occ
#endif

    ! allocate the WF coeffs. matrix
    allocate(c(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))

    ! if this is the first SCF step, then we need to initialize the WF coeffs. matrix with random
    ! numbers between -0.5 and 0.5 (normalize at the end to avoid instabilities), unless we are
    ! reading them from file
    if (ReadCoeffs) then
      write(k_point_name,'(i6)') ik
#ifdef MPI
      write(Node_name,'(i5)') Node
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        end if
      end if
#else
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))
        end if
      end if
#endif
      inquire(file=trim(WF_COEFFS_filename),exist=ReadCoeffs2)
    else
      ReadCoeffs2=.false.
    end if
    if (ReadCoeffs2) then
      call io_assign(i)
      open(i,file=trim(WF_COEFFS_filename),form='unformatted',status='old',action='read')
      read(i) c(ik,ispin)%mtrx
      call io_close(i)
    else
      if (ik==1) then
        if ((ispin==1) .or. (N_occ_diff/=0)) then
          call rand_init
          do i=1,h_dim_loc(2)
            do j=1,N_occ_loc(1,ispin)
              call random_number(rn)
              call random_number(rn2)
              call random_number(rn3)
              call random_number(rn4)
              c(ik,ispin)%mtrx(j,i)=cmplx(sign(0.5_dp*rn,rn2-0.5_dp),sign(0.5_dp*rn3,rn4-0.5_dp),dp)
            end do
          end do
          c(ik,ispin)%mtrx=1.0d-2*c(ik,ispin)%mtrx/sqrt(real(h_dim,dp))
        else
          c(ik,2)%mtrx=c(ik,1)%mtrx
        end if
      else
        rn2=999999.9_dp
        do i=1,ik-1
          rn=(kpoint(1,ik)-kpoint(1,i))**2+&
             (kpoint(2,ik)-kpoint(2,i))**2+&
             (kpoint(3,ik)-kpoint(3,i))**2
          if (rn<rn2) then
            rn2=rn
            j=i
          end if
        end do
        c(ik,ispin)%mtrx=c(j,ispin)%mtrx
      end if
    end if

    if (all(FirstCall(1:nk,1:2))) then
      t_precon_scale=fdf_physical('OMM.TPreconScale',10.0_dp,'Ry')
      cg_tol=fdf_get('OMM.RelTol',1.0d-9)
    end if

  end if

#ifdef MPI
  ! find the largest value of nhmax over all MPI processes needed for the sparse matrix
  ! multiplications
  if (UpdateSparseComm) then
    nhmax_max=0
    do i=0,Nodes-1
      if (i==Node) then
        j=nhmax
      end if
      call mpi_bcast(j,1,mpi_integer,i,mpi_comm_world,info)
      call mpi_bcast(k,1,mpi_integer,i,mpi_comm_world,info)
      if (j>nhmax_max) nhmax_max=j
    end do
  end if

  allocate(numh_recv(1:h_dim_loc_max))
  allocate(listhptr_recv(1:h_dim_loc_max))
  allocate(listh_recv(1:nhmax_max))
  allocate(h_dim_l2g_recv(1:h_dim_loc_max))
  allocate(As_recv(1:nhmax_max))
#endif

  if (CalcE) then

    ! calculate the energy-density matrix: e=c^T*[(2*I-S)*(H+eta*S)]*c
#ifdef MPI
    call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Hd(ik,ispin)%mtrx,1,1,desc2(1:9,ispin),cmplx_1,H(ik,ispin)%mtrx,&
                 1,1,desc2(1:9,ispin))
#else
    do k=1,N_occ
      do l=1,N_occ
        H(ik,ispin)%mtrx(k,l)=H(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Hd(ik,ispin)%mtrx(l,k))
      end do
    end do
#endif
    H(ik,ispin)%mtrx=H(ik,ispin)%mtrx+x_min(ik,ispin)*Hd(ik,ispin)%mtrx+x_min(ik,ispin)**2*Hdd(ik,ispin)%mtrx
    allocate(work1(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    call calc_densmat_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,H(ik,ispin)%mtrx+cmplx(eta,0.0_dp,dp)*S(ik,ispin)%mtrx,&
                             c(ik,ispin)%mtrx,d_sparse,work1,cd(ik,ispin)%mtrx)
    if (nspin==1) d_sparse=2.0_dp*d_sparse

    if (WriteCoeffs) then
      call io_assign(i)
      write(k_point_name,'(i6)') ik
#ifdef MPI
      write(Node_name,'(i5)') Node
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))//'.'//trim(adjustl(Node_name))
        end if
      end if
#else
      if (nspin==1) then
        WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS.'//trim(adjustl(k_point_name))
      else
        if (ispin==1) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_UP.'//trim(adjustl(k_point_name))
        else if (ispin==2) then
          WF_COEFFS_filename=trim(slabel)//'.WF_COEFFS_DOWN.'//trim(adjustl(k_point_name))
        end if
      end if
#endif
      open(i,file=trim(WF_COEFFS_filename),form='unformatted',status='replace',action='write')
      write(i) c(ik,ispin)%mtrx
      call io_close(i)
    end if

    deallocate(work1)
    deallocate(cd(ik,ispin)%mtrx)
    if (allocated(p_dense1D(ik)%mtrx)) deallocate(p_dense1D(ik)%mtrx)
    deallocate(sg(ik,ispin)%mtrx)
    deallocate(sc(ik,ispin)%mtrx)
    deallocate(S(ik,ispin)%mtrx)
    deallocate(Hdd(ik,ispin)%mtrx)
    deallocate(Hd(ik,ispin)%mtrx)
    deallocate(H(ik,ispin)%mtrx)

#ifdef MPI
    deallocate(As_recv)
    deallocate(h_dim_l2g_recv)
    deallocate(listh_recv)
    deallocate(listhptr_recv)
    deallocate(numh_recv)
#endif

    call timer('m_cg',2)

    return

  end if

  if (.not. allocated(H(ik,ispin)%mtrx)) allocate(H(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(Hd(ik,ispin)%mtrx)) allocate(Hd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(Hdd(ik,ispin)%mtrx)) allocate(Hdd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  if (.not. allocated(S(ik,ispin)%mtrx)) then
    allocate(S(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
    new_s=.true.
  else
    new_s=.false.
  end if
  if (.not. allocated(sc(ik,ispin)%mtrx)) allocate(sc(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  if (.not. allocated(sg(ik,ispin)%mtrx)) allocate(sg(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))

  ! calculate the preconditioning matrix (s+t/tau)^-1
  if (UpdatePrecon .and. (ispin==1)) then
    if (.not. allocated(p_dense1D(ik)%mtrx)) allocate(p_dense1D(ik)%mtrx(1:h_dim_loc(1),1:h_dim_loc(2)))
    p_dense1D(ik)%mtrx=s_dense1D+t_dense1D/t_precon_scale
#ifdef MPI
    allocate(ipiv(1:h_dim_loc(1)+BlockSize))
    call pzgetrf(h_dim,h_dim,p_dense1D(ik)%mtrx,1,1,desc1,ipiv,info)
    if (info/=0) call die('ERROR: pzgetrf has failed in minim!')
    allocate(work1(1:1,1:1))
    allocate(iwork(1:1))
    call pzgetri(h_dim,p_dense1D(ik)%mtrx,1,1,desc1,ipiv,work1,-1,iwork,-1,info)
    if (info/=0) call die('ERROR: pzgetri has failed in minim!')
    liwork=iwork(1)
    deallocate(iwork)
    lwork=work1(1,1)
    deallocate(work1)
    allocate(work1(1:lwork,1:1))
    allocate(iwork(1:liwork))
    call pzgetri(h_dim,p_dense1D(ik)%mtrx,1,1,desc1,ipiv,work1,lwork,iwork,liwork,info)
    if (info/=0) call die('ERROR: pzgetri has failed in minim!')
    deallocate(iwork)
    deallocate(work1)
    deallocate(ipiv)
#else
    allocate(ipiv(1:h_dim))
    lwork=h_dim*ilaenv(1,'zhetrf','U',h_dim,-1,-1,-1)
    allocate(work1(1:lwork,1:1))
    call zhetrf('U',h_dim,p_dense1D(ik)%mtrx,h_dim,ipiv,work1,lwork,info)
    if (info/=0) call die('ERROR: zhetrf has failed in minim!')
    deallocate(work1)
    allocate(work1(1:h_dim,1:1))
    call zhetri('U',h_dim,p_dense1D(ik)%mtrx,h_dim,ipiv,work1,info)
    if (info/=0) call die('ERROR: zhetri has failed in minim!')
    deallocate(work1)
    deallocate(ipiv)
    do i=1,h_dim-1
      do j=i+1,h_dim
        p_dense1D(ik)%mtrx(j,i)=p_dense1D(ik)%mtrx(i,j)
      end do
    end do
#endif
  end if

  allocate(Sd(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  allocate(Sdd(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
  allocate(g(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(g_p(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  if (UsePrecon) then
    allocate(pg(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
    allocate(pg_p(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  end if
  allocate(d(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(hc(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(hg(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  allocate(work1(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))

  ! first we calculate the energy and gradient for our initial guess, with the following steps:
  ! -calculate the hamiltonian in WF basis: H=c^T*h*c
  call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,h_sparse,c(ik,ispin)%mtrx,H(ik,ispin)%mtrx,hc)
  ! -calculate the overlap matrix in WF basis: S=c^T*s*c
  if (new_s .or. PreviousCallDiagon) then
    call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,s_sparse,c(ik,ispin)%mtrx,S(ik,ispin)%mtrx,sc(ik,ispin)%mtrx)
  else
    sc(ik,ispin)%mtrx=sc(ik,ispin)%mtrx+x_min(ik,ispin)*sg(ik,ispin)%mtrx
  end if
  ! -calculate the gradient: g=2*(2*h*c-s*c*H-h*c*S)
  !  (note that we *reuse* h*c and s*c contained in hc and sc from the previous call to calc_A)
  call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,sc(ik,ispin)%mtrx)
  ! -calculate the preconditioned gradient by premultiplying g by (s+t/tau)^-1
  if (UsePrecon) then
#ifdef MPI
    call pzgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,1,1,desc3(1:9,ispin),p_dense1D(ik)%mtrx,1,1,desc1,cmplx_0,pg,1,1,&
                desc3(1:9,ispin))
#else
    call zgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,N_occ,p_dense1D(ik)%mtrx,h_dim,cmplx_0,pg,N_occ)
#endif
  end if
  ! -calculate the additional matrices:
  !  Hd=g^T*h*c
  !  Sd=g^T*s*c
  !  Hdd=g^T*h*g
  !  Sdd=g^T*s*g
  !  (again, h*c has already been calculated, although h*g has not)
  !  and, finally, the coeffs. of the quartic line search equation in the direction g
  !  (the energy at c is given by the zeroth-order coeff. c(0))
  if (UsePrecon) then
#ifdef MPI
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),pg,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                desc2(1:9,ispin))
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),pg,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                desc2(1:9,ispin))
#else
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,pg,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,pg,N_occ,cmplx_0,Sd,N_occ)
#endif
  else
#ifdef MPI
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                desc2(1:9,ispin))
#else
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,g,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
#endif
#ifdef MPI
    call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),g,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                desc2(1:9,ispin))
#else
    call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,g,N_occ,cmplx_0,Sd,N_occ)
#endif
  end if
  if (UsePrecon) then
    call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,h_sparse,pg,Hdd(ik,ispin)%mtrx,hg)
    call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,s_sparse,pg,Sdd,sg(ik,ispin)%mtrx)
  else
    call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,h_sparse,g,Hdd(ik,ispin)%mtrx,hg)
    call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,s_sparse,g,Sdd,sg(ik,ispin)%mtrx)
  end if
  call calc_coeff(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,Hd(ik,ispin)%mtrx,Sd,Hdd(ik,ispin)%mtrx,Sdd,coeff,work1)
  E_OMM=coeff(0)

  ! this is the main loop of the CG algorithm. We perform a series of line minimizations, with the
  ! gradient g at each new step being modified to obtain the search direction d
  if (Node==0) then
    if ((ik==1) .and. (ispin==1)) then
      print'(a)', '+---------------------------------------------+'
      if (UsePrecon) then
        print'(a)', '| OMM (sparse algebra+preconditioning)        |'
      else
        print'(a)', '| OMM (sparse algebra)                        |'
      end if
      print'(a)', '+---------------------------------------------+'
    end if
    if (nk>1) then
      print'(a,i6,a,i6,a)', '| k point ', ik, ' of ', nk, '                    |'
    end if
    if (nspin==2) then
      if (ispin==1) then
        print'(a)', '| up spin                                     |'
      else
        print'(a)', '| down spin                                   |'
      end if
    end if
    if ((nk>1) .or. (nspin==2)) print'(a)', '+---------------------------------------------+'
    if (LongOut) print'(a)', '|             E_OMM            E_diff         |'
  end if
  conv=.false.
  d=cmplx_0
  icg=0
  do i=1,n_step_max
    lambda=0.0_dp
    do j=1,2*h_dim*N_occ-1
      if (UsePrecon) then
        d=pg+cmplx(lambda,0.0_dp,dp)*d
      else
        d=g+cmplx(lambda,0.0_dp,dp)*d
      end if
      g_p=g
      if (UsePrecon) pg_p=pg
      E_OMM_old=E_OMM
      ! if this is not the first CG step, we have to recalculate Hd, Sd, Hdd, Sdd, and the coeffs.
      if (icg>0) then
#ifdef MPI
        call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Hd(ik,ispin)%mtrx,1,1,&
                    desc2(1:9,ispin))
        call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,1,1,desc3(1:9,ispin),d,1,1,desc3(1:9,ispin),cmplx_0,Sd,1,1,&
                    desc2(1:9,ispin))
#else
        call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,hc,N_occ,d,N_occ,cmplx_0,Hd(ik,ispin)%mtrx,N_occ)
        call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,sc(ik,ispin)%mtrx,N_occ,d,N_occ,cmplx_0,Sd,N_occ)
#endif
        call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,h_sparse,d,Hdd(ik,ispin)%mtrx,hg)
        call calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,s_sparse,d,Sdd,sg(ik,ispin)%mtrx)
        call calc_coeff(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,Hd(ik,ispin)%mtrx,Sd,Hdd(ik,ispin)%mtrx,Sdd,coeff,work1)
      end if
      ! using the coeffs. calculated anlytically, we can find the minimum of the functional in the
      ! search direction, and calculate the energy at that minimum
      call solve_quartic(coeff(0:4),x_min(ik,ispin),ls_fail)
      ! in certain regions of the coeffs. space the line search gives no minimum--this occurs when there
      ! are positive eigenvalues in the eigenspecturm which are significantly occupied by our coeffs.
      ! matrix; the only known cure, unfortunately, is to scale down the entire matrix, thus returning to
      ! a  safe region of the coeffs. space.
      if (ls_fail) then
        if (Node==0) print'(a)', '| WARNING: Rescaling coefficients!            |'
        E_OMM=3.0*E_OMM
        c(ik,ispin)%mtrx=cmplx(0.5_dp,0.0_dp,dp)*c(ik,ispin)%mtrx
        ls_conv=.false.
      else
        ! if the line search is successful, move to the minimum
        E_OMM=coeff(4)*x_min(ik,ispin)**4+&
              coeff(3)*x_min(ik,ispin)**3+&
              coeff(2)*x_min(ik,ispin)**2+&
              coeff(1)*x_min(ik,ispin)+&
              coeff(0)
        c(ik,ispin)%mtrx=c(ik,ispin)%mtrx+x_min(ik,ispin)*d
        ls_conv=.true.
      end if
      ! recalculate S at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        S(ik,ispin)%mtrx=cmplx(0.25_dp,0.0_dp,dp)*S(ik,ispin)%mtrx
      else
#ifdef MPI
        call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Sd,1,1,desc2(1:9,ispin),cmplx_1,S(ik,ispin)%mtrx,1,1,&
                     desc2(1:9,ispin))
#else
        do k=1,N_occ
          do l=1,N_occ
            S(ik,ispin)%mtrx(k,l)=S(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Sd(l,k))
          end do
        end do
#endif
        S(ik,ispin)%mtrx=S(ik,ispin)%mtrx+x_min(ik,ispin)*Sd+x_min(ik,ispin)**2*Sdd
      end if
      E_diff=2.0_dp*abs((E_OMM-E_OMM_old)/(E_OMM+E_OMM_old))
      if ((Node==0) .and. LongOut) print'(a,2(1x,i5),2(1x,es15.7e3),1x,a)', '|', i, j, E_OMM, E_diff, '|'
      icg=icg+1
      if (E_diff<=cg_tol) then
        conv=.true.
        exit
      end if
      ! recalculate H at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        H(ik,ispin)%mtrx=cmplx(0.25_dp,0.0_dp,dp)*H(ik,ispin)%mtrx
      else
#ifdef MPI
        call pzgeadd('C',N_occ,N_occ,cmplx(x_min(ik,ispin),0.0_dp,dp),Hd(ik,ispin)%mtrx,1,1,desc2(1:9,ispin),cmplx_1,&
                     H(ik,ispin)%mtrx,1,1,desc2(1:9,ispin))
#else
        do k=1,N_occ
          do l=1,N_occ
            H(ik,ispin)%mtrx(k,l)=H(ik,ispin)%mtrx(k,l)+x_min(ik,ispin)*conjg(Hd(ik,ispin)%mtrx(l,k))
          end do
        end do
#endif
        H(ik,ispin)%mtrx=H(ik,ispin)%mtrx+x_min(ik,ispin)*Hd(ik,ispin)%mtrx+x_min(ik,ispin)**2*Hdd(ik,ispin)%mtrx
      end if
      ! recalculate g at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        hc=cmplx(0.5_dp,0.0_dp,dp)*hc
        sc(ik,ispin)%mtrx=cmplx(0.5_dp,0.0_dp,dp)*sc(ik,ispin)%mtrx
        g=g_p+cmplx(1.5_dp,0.0_dp,dp)*hc
      else
        hc=hc+x_min(ik,ispin)*hg
        sc(ik,ispin)%mtrx=sc(ik,ispin)%mtrx+x_min(ik,ispin)*sg(ik,ispin)%mtrx
        call calc_grad(h_dim,N_occ,ispin,H(ik,ispin)%mtrx,S(ik,ispin)%mtrx,g,hc,sc(ik,ispin)%mtrx)
      end if
      if (UsePrecon) then
#ifdef MPI
        call pzgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,1,1,desc3(1:9,ispin),p_dense1D(ik)%mtrx,1,1,desc1,cmplx_0,pg,1,1,&
                    desc3(1:9,ispin))
#else
        call zgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,g,N_occ,p_dense1D(ik)%mtrx,h_dim,cmplx_0,pg,N_occ)
#endif
      end if
      if (ls_conv) then
        if (UsePrecon) then
          lambda_n=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_n=lambda_n+real(pg(k,l),dp)*(real(g(k,l),dp)-real(g_p(k,l),dp))
              lambda_n=lambda_n+aimag(pg(k,l))*(aimag(g(k,l))-aimag(g_p(k,l)))
            end do
          end do
          lambda_d=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_d=lambda_d+real(pg_p(k,l),dp)*real(g_p(k,l),dp)
              lambda_d=lambda_d+aimag(pg_p(k,l))*aimag(g_p(k,l))
            end do
          end do
        else
          lambda_n=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_n=lambda_n+real(g(k,l),dp)*(real(g(k,l),dp)-real(g_p(k,l),dp))
              lambda_n=lambda_n+aimag(g(k,l))*(aimag(g(k,l))-aimag(g_p(k,l)))
            end do
          end do
          lambda_d=0.0_dp
          do k=1,N_occ_loc(1,ispin)
            do l=1,h_dim_loc(2)
              lambda_d=lambda_d+real(g_p(k,l),dp)*real(g_p(k,l),dp)
              lambda_d=lambda_d+aimag(g_p(k,l))*aimag(g_p(k,l))
            end do
          end do
        end if
#ifdef MPI
        call mpi_allreduce(lambda_n,lambda_n_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        call mpi_allreduce(lambda_d,lambda_d_tot,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
        lambda=lambda_n_tot/lambda_d_tot
#else
        lambda=lambda_n/lambda_d
#endif
      else
        exit
      end if
    end do
    if (conv) exit
  end do
  if (i>n_step_max) then
    if (Node==0) print'(a)', '| WARNING: OMM failed to converge!            |'
  end if
  if (Node==0 .and. LongOut) print'(a)', '+---------------------------------------------+'

  deallocate(work1)
  deallocate(hg)
  deallocate(hc)
  deallocate(d)
  if (UsePrecon) then
    deallocate(pg_p)
    deallocate(pg)
  end if
  deallocate(g_p)
  deallocate(g)
  deallocate(Sdd)
  deallocate(Sd)

  ! calculate the density matrix: d=c*(2*I-S)*c^T
  if (.not. allocated(twoI(ispin)%mtrx)) then
    allocate(twoI(ispin)%mtrx(1:N_occ_loc(1,ispin),1:N_occ_loc(2,ispin)))
#ifdef MPI
    call pzlaset('A',N_occ,N_occ,cmplx_0,cmplx_2,twoI(ispin)%mtrx,1,1,desc2(1:9,ispin))
#else
    twoI(ispin)%mtrx=cmplx_0
    do i=1,N_occ
      twoI(ispin)%mtrx(i,i)=cmplx_2
    end do
#endif
  end if
  if (.not. allocated(cd(ik,ispin)%mtrx)) allocate(cd(ik,ispin)%mtrx(1:N_occ_loc(1,ispin),1:h_dim_loc(2)))
  call calc_densmat_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,twoI(ispin)%mtrx-S(ik,ispin)%mtrx,c(ik,ispin)%mtrx,&
                           d_sparse,cd(ik,ispin)%mtrx)
  if (nspin==1) d_sparse=2.0_dp*d_sparse

#ifdef MPI
    deallocate(As_recv)
    deallocate(h_dim_l2g_recv)
    deallocate(listh_recv)
    deallocate(listhptr_recv)
    deallocate(numh_recv)
#endif

  ! calculate the trace of S to make sure we are occupying the right number of eigenstates in our
  ! solution
#ifdef MPI
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+(real(twoI(ispin)%mtrx(i,j),dp)-real(S(ik,ispin)%mtrx(i,j),dp))*real(S(ik,ispin)%mtrx(i,j),dp)+&
                    (aimag(twoI(ispin)%mtrx(i,j))-aimag(S(ik,ispin)%mtrx(i,j)))*aimag(S(ik,ispin)%mtrx(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrQS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  TrQS=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrQS=TrQS+(real(twoI(ispin)%mtrx(i,j),dp)-real(S(ik,ispin)%mtrx(i,j),dp))*real(S(ik,ispin)%mtrx(i,j),dp)+&
                (aimag(twoI(ispin)%mtrx(i,j))-aimag(S(ik,ispin)%mtrx(i,j)))*aimag(S(ik,ispin)%mtrx(i,j))
    end do
  end do
#endif
  if (Node==0) then
    if (nspin==1) then
      print'(a,i5,a)',    '| minim: icg             = ', icg, '              |'
      print'(a,f13.7,a)', '| minim: 2*Tr[(2*I-S)*S] = ', 2.0_dp*TrQS, '      |'
    else
      print'(a,i5,a)',    '| minim: icg           = ', icg, '                |'
      print'(a,f13.7,a)', '| minim: Tr[(2*I-S)*S] = ', TrQS, '        |'
    end if
    print'(a)',       '+---------------------------------------------+'
  end if

  if (FirstCall(ik,ispin)) FirstCall(ik,ispin)=.false.

  call timer('m_cg',2)

end subroutine minim_cg_sparse

!================================================!
! calculate the coeffs. of the quartic line      !
! search equation from three energy points and   !
! two gradient points                            !
!================================================!
subroutine fit_quartic(x,y,g,c)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: x(1:3) ! three x-points {x_i}
  real(dp), intent(in) :: y(1:3) ! y(x_i) at the three points
  real(dp), intent(in) :: g(1:2) ! (dy/dx)|x_i at the three points

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: c(0:4) ! coeffs. of the quartic equation

  !**********************************************!

  !call timer('m_solve_quartic',1)

  ! the following expressions for the coeffs. were produced automatically using Maple 12
  c(4)=(x(3)**3*x(2)*g(1)-3*x(1)*x(2)**2*y(1)+3*y(3)*x(1)*x(2)**2+x(1)**2*x(2)**2*g(1)+x(3)*x(2)**3*&
       g(1)+2*x(1)**2*x(3)**2*g(2)-3*x(2)*x(3)**2*y(1)+3*y(2)*x(1)**2*x(2)-x(3)**3*x(1)*g(1)+x(3)**3*&
       x(2)*g(2)-x(2)**2*x(3)**2*g(2)-x(1)**2*x(2)**2*g(2)-2*x(2)**2*x(3)**2*g(1)+3*x(2)*x(3)**2*&
       y(2)+x(1)**2*x(3)**2*g(1)-x(1)*x(2)**3*g(1)-3*x(1)*x(3)**2*y(1)-x(3)**3*x(1)*g(2)+3*x(1)*&
       x(3)**2*y(2)+x(2)*g(2)*x(1)**3-3*y(3)*x(1)**2*x(2)-x(3)*g(2)*x(1)**3+6*x(1)*x(3)*x(2)*y(1)+2*&
       x(1)*x(3)*x(2)**2*g(2)+x(1)*x(3)**2*x(2)*g(1)-x(1)*x(3)**2*x(2)*g(2)-2*x(1)**2*x(3)*x(2)*g(1)-&
       x(1)**2*x(3)*x(2)*g(2)+x(1)*x(3)*x(2)**2*g(1)-6*x(1)*x(3)*x(2)*y(2)+2*x(3)**3*y(1)-2*x(3)**3*&
       y(2)+x(2)**3*y(1)+y(3)*x(1)**3-y(3)*x(2)**3-y(2)*x(1)**3)/(-2*x(3)**3*x(1)**4+x(3)**4*x(1)**3-&
       x(1)**2*x(2)**5-x(3)**4*x(2)**3-x(2)**5*x(3)**2-3*x(1)**4*x(2)**3+2*x(3)**3*x(2)**4+x(1)**5*&
       x(3)**2+3*x(2)**4*x(1)**3+4*x(3)**3*x(2)*x(1)**3-4*x(3)**3*x(1)*x(2)**3+2*x(1)*x(3)*x(2)**5+4*&
       x(1)**4*x(3)*x(2)**2+8*x(1)**2*x(3)**2*x(2)**3+x(1)**4*x(3)**2*x(2)-x(1)*x(3)**2*x(2)**4-2*&
       x(1)**5*x(3)*x(2)-4*x(1)**2*x(3)*x(2)**4+x(1)**5*x(2)**2-8*x(2)**2*x(3)**2*x(1)**3-3*x(3)**4*&
       x(1)**2*x(2)+3*x(3)**4*x(1)*x(2)**2)

  c(3)=-(-x(1)*g(1)+2*c(4)*x(1)**4-x(1)*g(2)+4*x(1)*c(4)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3+x(2)*&
       g(2)-2*c(4)*x(2)**4+2*y(1)-2*y(2))/(x(1)**3+3*x(1)*x(2)**2-3*x(2)*x(1)**2-x(2)**3)

  c(2)=-(-y(2)+c(4)*x(2)**4+c(3)*x(2)**3+x(2)*g(1)-4*x(2)*c(4)*x(1)**3-3*x(2)*c(3)*x(1)**2+y(1)+3*&
       c(4)*x(1)**4+2*c(3)*x(1)**3-x(1)*g(1))/(x(1)**2-2*x(1)*x(2)+x(2)**2)

  c(1)=g(1)-4*c(4)*x(1)**3-3*c(3)*x(1)**2-2*c(2)*x(1)

  c(0)=y(1)-c(4)*x(1)**4-c(3)*x(1)**3-c(2)*x(1)**2-c(1)*x(1)

  !if (Node==0) print*, 'f(x)=',c(4),'*x**4+',c(3),'*x**3+',c(2),'*x**2+',c(1),'*x+',c(0)

  !call timer('m_fit_quartic',2)

end subroutine fit_quartic

!================================================!
! find the minimum for the quartic line search   !
! equation                                       !
!================================================!
subroutine solve_quartic(c,x_min,fail)
  implicit none

  !**** INPUT ***********************************!

  real(dp), intent(in) :: c(0:4) ! coeffs. of the quartic equation

  !**** OUTPUT **********************************!

  logical, intent(out) :: fail ! did we fail to find a minimum?

  real(dp), intent(out) :: x_min ! position of minimum

  !**** LOCAL ***********************************!

  integer :: i
  integer :: x_order(1:3)

  real(dp) :: t(1:3)
  real(dp) :: z(1:3)
  real(dp) :: a
  real(dp) :: b
  real(dp) :: d
  real(dp) :: Q
  real(dp) :: R
  real(dp) :: theta
  real(dp) :: S
  real(dp) :: U

  !**********************************************!

  !call timer('m_solve_quartic',1)

  fail=.false.

  !if (c(4)<0.0_dp) then
  !  if (Node==0) print*, '#WARNING: Function is unbounded!'
  !  !stop
  !end if

  ! in order to find the minimum of the quartic equation, we have to solve a cubic equation; the
  ! following method is taken from Numerical Recipes
  a=3.0_dp*c(3)/(4.0_dp*c(4))
  b=2.0_dp*c(2)/(4.0_dp*c(4))
  if ((abs(b)>=1.0d11) .or. (abs(c(4))<=1.0d-11)) then
    !if (Node==0) print*, '#WARNING: Function is quadratic!'
    x_min=-0.5_dp*c(1)/c(2)
    return
  end if
  d=c(1)/(4.0_dp*c(4))

  Q=(a**2-3.0_dp*b)/9.0_dp
  R=(2.0_dp*a**3-9.0_dp*a*b+27.0_dp*d)/54.0_dp
  if (R**2<Q**3) then
    theta=acos(R/sqrt(Q**3))
    t(1)=-2.0_dp*sqrt(Q)*cos(theta/3.0_dp)-a/3.0_dp
    t(2)=-2.0_dp*sqrt(Q)*cos((theta+2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    t(3)=-2.0_dp*sqrt(Q)*cos((theta-2.0_dp*Pi)/3.0_dp)-a/3.0_dp
    z(1:3)=c(4)*t(1:3)**4+c(3)*t(1:3)**3+c(2)*t(1:3)**2+c(1)*t(1:3)+c(0)
    if (c(4)>0.0_dp) then
      if (all(z(1)>=z(2:3))) then
        x_order(1:3)=(/1,2,3/)
      else if (z(2)>z(3)) then
        x_order(1:3)=(/2,3,1/)
      else
        x_order(1:3)=(/3,1,2/)
      end if
      if ((0.0_dp<=t(x_order(1))) .and. (t(x_order(2))<=t(x_order(1)))) then
        x_min=t(x_order(2))
      else
        x_min=t(x_order(3))
      end if
    else
      if (all(z(1)<=z(2:3))) then
        x_min=t(1)
      else if (z(2)<z(3)) then
        x_min=t(2)
      else
        x_min=t(3)
      end if
    end if
  else
    S=-sign(1.0_dp,R)*(abs(R)+sqrt(R**2-Q**3))**(1.0_dp/3.0_dp)
    if (S==0.0_dp) then
      U=0.0_dp
    else
      U=Q/S
    end if
    x_min=(S+U)-(a/3.0_dp)
    if (c(4)<0.0_dp) fail=.true.
  end if

  !call timer('m_solve_quartic',2)

end subroutine solve_quartic

!================================================!
! calculate the gradient of the Kim functional:  !
! g=2*(2*h*C-s*C*H-h*C*S)                        !
!================================================!
subroutine calc_grad(h_dim,N_occ,ispin,H,S,grad,hc,sc)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim ! num. of AOs (global)
  integer, intent(in) :: N_occ ! num. of WFs (global)
  integer, intent(in) :: ispin ! up/down spin (1/2)

  complex(dp), intent(in) :: H(:,:) ! hamiltonian matrix in WF basis
  complex(dp), intent(in) :: S(:,:) ! overlap matrix in WF basis

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: grad(:,:) ! gradient of Kim functional
  complex(dp), intent(inout) :: hc(:,:)   ! h*c
  complex(dp), intent(inout) :: sc(:,:)   ! s*c

  !**********************************************!

  call timer('m_calc_grad',1)

  grad=cmplx(4.0_dp,0.0_dp,dp)*hc

#ifdef MPI
  call pzgemm('N','N',N_occ,h_dim,N_occ,-cmplx_2,S,1,1,desc2(1:9,ispin),hc,1,1,desc3(1:9,ispin),cmplx_1,grad,1,1,&
              desc3(1:9,ispin))
  call pzgemm('N','N',N_occ,h_dim,N_occ,-cmplx_2,H,1,1,desc2(1:9,ispin),sc,1,1,desc3(1:9,ispin),cmplx_1,grad,1,1,&
              desc3(1:9,ispin))
#else
  call zgemm('N','N',N_occ,h_dim,N_occ,-cmplx_2,S,N_occ,hc,N_occ,cmplx_1,grad,N_occ)
  call zgemm('N','N',N_occ,h_dim,N_occ,-cmplx_2,H,N_occ,sc,N_occ,cmplx_1,grad,N_occ)
#endif

  call timer('m_calc_grad',2)

end subroutine calc_grad

!================================================!
! calculate operator matrix in WF basis (dense   !
! routine)                                       !
!================================================!
subroutine calc_A(h_dim,N_occ,ispin,Ap,c,A,Apc)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim ! num. of AOs (global)
  integer, intent(in) :: N_occ ! num. of WFs (global)
  integer, intent(in) :: ispin ! up/down spin (1/2)

  complex(dp), intent(in) :: Ap(:,:) ! operator matrix in AO basis
  complex(dp), intent(in) :: c(:,:)  ! WF coeffs. matrix

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: A(:,:)   ! operator matrix in WF basis
  complex(dp), intent(inout) :: Apc(:,:) ! work matrix

  !**********************************************!

  call timer('m_calc_A',1)

#ifdef MPI
  call pzgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,c,  1,1,desc3(1:9,ispin),Ap,1,1,desc1,           cmplx_0,Apc,1,1,desc3(1:9,ispin))
  call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,Apc,1,1,desc3(1:9,ispin),c, 1,1,desc3(1:9,ispin),cmplx_0,A,  1,1,desc2(1:9,ispin))
#else
  call zgemm('N','N',N_occ,h_dim,h_dim,cmplx_1,c,  N_occ,Ap,h_dim,cmplx_0,Apc,N_occ)
  call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,Apc,N_occ,c, N_occ,cmplx_0,A,  N_occ)
#endif

  call timer('m_calc_A',2)

end subroutine calc_A

!================================================!
! calculate operator matrix in WF basis (sparse  !
! routine)                                       !
!================================================!
subroutine calc_A_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,As,c,A,Asc)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim       ! num. of AOs (global)
  integer, intent(in) :: N_occ       ! num. of WFs (global)
  integer, intent(in) :: ispin       ! up/down spin (1/2)
  integer, intent(in) :: nhmax       ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(:)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(:) ! pointer to start of row in listh
  integer, intent(in) :: listh(:)    ! list of nonzero elements of each row of sparse matrices

  complex(dp), intent(in) :: As(:)  ! operator matrix in AO basis (sparse)
  complex(dp), intent(in) :: c(:,:) ! WF coeffs. matrix

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: A(:,:)   ! operator matrix in WF basis
  complex(dp), intent(inout) :: Asc(:,:) ! work matrix

  !**** LOCAL ***********************************!

  integer :: i
  integer :: j
  integer :: l
#ifdef MPI
  integer :: n_comm
  integer :: info
  integer :: h_dim_loc_recv
  integer :: nhmax_recv
#endif

  !**********************************************!

  call timer('m_calc_A',1)

#ifdef MPI
  Asc=0.0_dp
  do n_comm=0,Nodes-1
    if (n_comm==Node) then
      h_dim_loc_recv=h_dim_loc_1D
      nhmax_recv=nhmax
      numh_recv(1:h_dim_loc_1D)=numh(1:h_dim_loc_1D)
      listhptr_recv(1:h_dim_loc_1D)=listhptr(1:h_dim_loc_1D)
      listh_recv(1:nhmax)=listh(1:nhmax)
      As_recv(1:nhmax)=As(1:nhmax)
      h_dim_l2g_recv(1:h_dim_loc_1D)=h_dim_l2g(1:h_dim_loc_1D)
    end if
    call mpi_bcast(h_dim_loc_recv,   1,              mpi_integer,       n_comm,mpi_comm_world,info)
    call mpi_bcast(nhmax_recv,       1,              mpi_integer,       n_comm,mpi_comm_world,info)
    call mpi_bcast(numh_recv(1),     h_dim_loc_recv, mpi_integer,       n_comm,mpi_comm_world,info)
    call mpi_bcast(listhptr_recv(1), h_dim_loc_recv, mpi_integer,       n_comm,mpi_comm_world,info)
    call mpi_bcast(listh_recv(1),    nhmax_recv,     mpi_integer,       n_comm,mpi_comm_world,info)
    call mpi_bcast(As_recv(1),       nhmax_recv,     mpi_double_complex,n_comm,mpi_comm_world,info)
    call mpi_bcast(h_dim_l2g_recv(1),h_dim_loc_recv, mpi_integer,       n_comm,mpi_comm_world,info)
    if (n_comm==Node) then
      do i=1,h_dim_loc_1D
        do j=1,numh(i)
          l=listhptr(i)+j
          Asc(:,h_dim_l2g(i))=Asc(:,h_dim_l2g(i))+As(l)*c(:,listh(l))
        end do
      end do
    else
      do i=1,h_dim_loc_recv
        do j=1,numh_recv(i)
          l=listhptr_recv(i)+j
          Asc(:,h_dim_l2g_recv(i))=Asc(:,h_dim_l2g_recv(i))+As_recv(l)*c(:,listh_recv(l))
        end do
      end do
    end if
  end do
  call pzgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,Asc,1,1,desc3(1:9,ispin),c,1,1,desc3(1:9,ispin),cmplx_0,A,1,1,desc2(1:9,ispin))
#else
  Asc=cmplx_0
  do i=1,h_dim
    do j=1,numh(i)
      l=listhptr(i)+j
      Asc(:,i)=Asc(:,i)+As(l)*c(:,listh(l))
    end do
  end do
  call zgemm('N','C',N_occ,N_occ,h_dim,cmplx_1,Asc,N_occ,c,N_occ,cmplx_0,A,N_occ)
#endif

  call timer('m_calc_A',2)

end subroutine calc_A_sparse

!================================================!
! calculate operator matrix in AO basis (dense   !
! routine)                                       !
!================================================!
subroutine calc_densmat(h_dim,N_occ,ispin,A,c1,Ap,cA,c2)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim ! num. of AOs (global)
  integer, intent(in) :: N_occ ! num. of WFs (global)
  integer, intent(in) :: ispin ! up/down spin (1/2)

  complex(dp), intent(in) :: A(:,:)            ! Operator matrix in WF basis
  complex(dp), intent(in) :: c1(:,:)           ! WF coeffs. matrix
  complex(dp), intent(in), optional :: c2(:,:) ! pre-multiplied WF coeffs. matrix

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: Ap(:,:) ! operator matrix in orbital basis
  complex(dp), intent(inout) :: cA(:,:) ! work matrix

  !**********************************************!

  call timer('m_calc_densmat',1)

#ifdef MPI
  if (present(c2)) then
    call pzgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,1,1,desc2(1:9,ispin),c2,1,1,desc3(1:9,ispin),cmplx_0,cA,1,1,desc3(1:9,ispin))
  else
    call pzgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,1,1,desc2(1:9,ispin),c1,1,1,desc3(1:9,ispin),cmplx_0,cA,1,1,desc3(1:9,ispin))
  end if
  call pzgemm('C','N',h_dim,h_dim,N_occ,cmplx_1,c1,1,1,desc3(1:9,ispin),cA,1,1,desc3(1:9,ispin),cmplx_0,Ap,1,1,desc1)
#else
  if (present(c2)) then
    call zgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,N_occ,c2,N_occ,cmplx_0,cA,N_occ)
  else
    call zgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,N_occ,c1,N_occ,cmplx_0,cA,N_occ)
  end if
  call zgemm('C','N',h_dim,h_dim,N_occ,cmplx_1,c1,N_occ,cA,N_occ,cmplx_0,Ap,h_dim)
#endif

  call timer('m_calc_densmat',2)

end subroutine calc_densmat

!================================================!
! calculate operator matrix in AO basis (sparse  !
! routine)                                       !
!================================================!
subroutine calc_densmat_sparse(h_dim,N_occ,ispin,nhmax,numh,listhptr,listh,A,c1,As,cA,c2)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim       ! num. of AOs (global)
  integer, intent(in) :: N_occ       ! num. of WFs (global)
  integer, intent(in) :: ispin       ! up/down spin (1/2)
  integer, intent(in) :: nhmax       ! first dimension of listh and sparse matrices
  integer, intent(in) :: numh(:)     ! num. of nonzero elements of each row of sparse matrices
  integer, intent(in) :: listhptr(:) ! pointer to start of row in listh
  integer, intent(in) :: listh(:)    ! list of nonzero elements of each row of sparse matrices

  complex(dp), intent(in) :: A(:,:)            ! operator matrix in WF basis
  complex(dp), intent(in) :: c1(:,:)           ! WF coeffs. matrix
  complex(dp), intent(in), optional :: c2(:,:) ! pre-multiplied WF coeffs. matrix

  !**** OUTPUT **********************************!

  complex(dp), intent(out) :: As(:) ! operator matrix in AO basis (sparse)

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: cA(:,:) ! work matrix

  !**** LOCAL ***********************************!

  integer :: i
  integer :: j
  integer :: k
  integer :: l
#ifdef MPI
  integer :: n_comm
  integer :: info
  integer :: h_dim_loc_recv
  integer :: nhmax_recv
#endif

  !**********************************************!

  call timer('m_calc_densmat',1)

#ifdef MPI
  if (present(c2)) then
    call pzgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,1,1,desc2(1:9,ispin),c2,1,1,desc3(1:9,ispin),cmplx_0,cA,1,1,desc3(1:9,ispin))
  else
    call pzgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,1,1,desc2(1:9,ispin),c1,1,1,desc3(1:9,ispin),cmplx_0,cA,1,1,desc3(1:9,ispin))
  end if
  do n_comm=0,Nodes-1
    if (n_comm==Node) then
      h_dim_loc_recv=h_dim_loc_1D
      nhmax_recv=nhmax
      numh_recv(1:h_dim_loc_1D)=numh(1:h_dim_loc_1D)
      listhptr_recv(1:h_dim_loc_1D)=listhptr(1:h_dim_loc_1D)
      listh_recv(1:nhmax)=listh(1:nhmax)
      h_dim_l2g_recv(1:h_dim_loc_1D)=h_dim_l2g(1:h_dim_loc_1D)
    end if
    call mpi_bcast(h_dim_loc_recv,   1,              mpi_integer,n_comm,mpi_comm_world,info)
    call mpi_bcast(nhmax_recv,       1,              mpi_integer,n_comm,mpi_comm_world,info)
    call mpi_bcast(numh_recv(1),     h_dim_loc_recv, mpi_integer,n_comm,mpi_comm_world,info)
    call mpi_bcast(listhptr_recv(1), h_dim_loc_recv, mpi_integer,n_comm,mpi_comm_world,info)
    call mpi_bcast(listh_recv(1),    nhmax_recv,     mpi_integer,n_comm,mpi_comm_world,info)
    call mpi_bcast(h_dim_l2g_recv(1),h_dim_loc_recv, mpi_integer,n_comm,mpi_comm_world,info)
    As_recv=0.0_dp
    do i=1,h_dim_loc_recv
      do j=1,numh_recv(i)
        l=listhptr_recv(i)+j
        do k=1,N_occ_loc_1D(ispin)
          As_recv(l)=As_recv(l)+cA(k,listh_recv(l))*conjg(c1(k,h_dim_l2g_recv(i)))
        end do
      end do
    end do
    call mpi_reduce(As_recv(1),As(1),nhmax_recv,mpi_double_complex,mpi_sum,n_comm,mpi_comm_world,info)
  end do
#else
  if (present(c2)) then
    call zgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,N_occ,c2,N_occ,cmplx_0,cA,N_occ)
  else
    call zgemm('N','N',N_occ,h_dim,N_occ,cmplx_1,A,N_occ,c1,N_occ,cmplx_0,cA,N_occ)
  end if
  As=0.0_dp
  do i=1,h_dim
    do j=1,numh(i)
      l=listhptr(i)+j
      do k=1,N_occ
        As(l)=As(l)+cA(k,listh(l))*conjg(c1(k,i))
      end do
    end do
  end do
#endif

  call timer('m_calc_densmat',2)

end subroutine calc_densmat_sparse

!================================================!
! calculate coeffs. of the quartic line search   !
! equation using analytical expressions          !
!================================================!
subroutine calc_coeff(h_dim,N_occ,ispin,H,S,Hd,Sd,Hdd,Sdd,coeff,SdH)
  implicit none

  !**** INPUT ***********************************!

  integer, intent(in) :: h_dim ! num. of AOs (global)
  integer, intent(in) :: N_occ ! num. of WFs (global)
  integer, intent(in) :: ispin ! up/down spin (1/2)

  complex(dp), intent(in) :: H(:,:)   ! hamiltonian matrix in WF basis
  complex(dp), intent(in) :: S(:,:)   ! overlap matrix in WF basis
  complex(dp), intent(in) :: Hd(:,:)  ! g^T*h*c
  complex(dp), intent(in) :: Sd(:,:)  ! g^T*s*c
  complex(dp), intent(in) :: Hdd(:,:) ! g^T*h*g
  complex(dp), intent(in) :: Sdd(:,:) ! g^T*h*g

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: coeff(0:4) ! coeffs. of the quartic equation

  complex(dp), intent(inout) :: SdH(:,:) ! work matrix

  !**** LOCAL ***********************************!

  integer :: i, j
#ifdef MPI
  integer :: info
#endif

  real(dp) :: Tr_loc
  real(dp) :: TrH
  real(dp) :: TrHS
  real(dp) :: TrHd
  real(dp) :: TrHdS
  real(dp) :: TrHSd
  real(dp) :: TrHdd
  real(dp) :: TrHddS
  real(dp) :: TrHSdd
  real(dp) :: TrHdSd
  real(dp) :: TrHdSdH
  real(dp) :: TrHddSd
  real(dp) :: TrHdSdd
  real(dp) :: TrHddSdd

#ifdef MPI
  complex(dp), external :: pzlatra
#endif

  !**********************************************!

  call timer('m_calc_coeff',1)

#ifdef MPI
  TrH=pzlatra(N_occ,H,1,1,desc2(1:9,ispin))
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(H(i,j),dp)*real(S(i,j),dp)+aimag(H(i,j))*aimag(S(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  TrHd=pzlatra(N_occ,Hd,1,1,desc2(1:9,ispin))
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hd(i,j),dp)*real(S(i,j),dp)+aimag(Hd(i,j))*aimag(S(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHdS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(H(i,j),dp)*real(Sd(i,j),dp)+aimag(H(i,j))*aimag(Sd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  TrHdd=pzlatra(N_occ,Hdd,1,1,desc2(1:9,ispin))
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hdd(i,j),dp)*real(S(i,j),dp)+aimag(Hdd(i,j))*aimag(S(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHddS,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(H(i,j),dp)*real(Sdd(i,j),dp)+aimag(H(i,j))*aimag(Sdd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  call pztranc(N_occ,N_occ,cmplx_1,Sd,1,1,desc2(1:9,ispin),cmplx_0,SdH,1,1,desc2(1:9,ispin))
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hd(i,j),dp)*real(SdH(i,j),dp)+aimag(Hd(i,j))*aimag(SdH(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHdSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hd(i,j),dp)*real(Sd(i,j),dp)+aimag(Hd(i,j))*aimag(Sd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHdSdH,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hdd(i,j),dp)*real(Sd(i,j),dp)+aimag(Hdd(i,j))*aimag(Sd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHddSd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hd(i,j),dp)*real(Sdd(i,j),dp)+aimag(Hd(i,j))*aimag(Sdd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHdSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)

  Tr_loc=0.0_dp
  do i=1,N_occ_loc(1,ispin)
    do j=1,N_occ_loc(2,ispin)
      Tr_loc=Tr_loc+real(Hdd(i,j),dp)*real(Sdd(i,j),dp)+aimag(Hdd(i,j))*aimag(Sdd(i,j))
    end do
  end do
  call mpi_allreduce(Tr_loc,TrHddSdd,1,mpi_double_precision,mpi_sum,mpi_comm_world,info)
#else
  TrH=0.0_dp
  do i=1,N_occ
    TrH=TrH+real(H(i,i),dp)
  end do
  TrHS=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHS=TrHS+real(H(i,j),dp)*real(S(i,j),dp)+aimag(H(i,j))*aimag(S(i,j))
    end do
  end do

  TrHd=0.0_dp
  do i=1,N_occ
    TrHd=TrHd+real(Hd(i,i),dp)
  end do
  TrHdS=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHdS=TrHdS+real(Hd(i,j),dp)*real(S(i,j),dp)+aimag(Hd(i,j))*aimag(S(i,j))
    end do
  end do
  TrHSd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHSd=TrHSd+real(H(i,j),dp)*real(Sd(i,j),dp)+aimag(H(i,j))*aimag(Sd(i,j))
    end do
  end do

  TrHdd=0.0_dp
  do i=1,N_occ
    TrHdd=TrHdd+real(Hdd(i,i),dp)
  end do
  TrHddS=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHddS=TrHddS+real(Hdd(i,j),dp)*real(S(i,j),dp)+aimag(Hdd(i,j))*aimag(S(i,j))
    end do
  end do
  TrHSdd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHSdd=TrHSdd+real(H(i,j),dp)*real(Sdd(i,j),dp)+aimag(H(i,j))*aimag(Sdd(i,j))
    end do
  end do
  TrHdSd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHdSd=TrHdSd+real(Hd(i,j),dp)*real(Sd(j,i),dp)-aimag(Hd(i,j))*aimag(Sd(j,i))
    end do
  end do
  TrHdSdH=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHdSdH=TrHdSdH+real(Hd(i,j),dp)*real(Sd(i,j),dp)+aimag(Hd(i,j))*aimag(Sd(i,j))
    end do
  end do

  TrHddSd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHddSd=TrHddSd+real(Hdd(i,j),dp)*real(Sd(i,j),dp)+aimag(Hdd(i,j))*aimag(Sd(i,j))
    end do
  end do
  TrHdSdd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHdSdd=TrHdSdd+real(Hd(i,j),dp)*real(Sdd(i,j),dp)+aimag(Hd(i,j))*aimag(Sdd(i,j))
    end do
  end do

  TrHddSdd=0.0_dp
  do i=1,N_occ
    do j=1,N_occ
      TrHddSdd=TrHddSdd+real(Hdd(i,j),dp)*real(Sdd(i,j),dp)+aimag(Hdd(i,j))*aimag(Sdd(i,j))
    end do
  end do
#endif

  coeff(0)=2.0_dp*TrH-TrHS
  coeff(1)=2.0_dp*(2.0_dp*TrHd-TrHdS-TrHSd)
  coeff(2)=2.0_dp*(TrHdd-TrHdSd-TrHdSdH)-TrHddS-TrHSdd
  coeff(3)=-2.0_dp*(TrHddSd+TrHdSdd)
  coeff(4)=-TrHddSdd

  call timer('m_calc_coeff',2)

end subroutine calc_coeff

!================================================!
! random number generator                        !
! -initialize with:                              !
!  call rand_init()                              !
! -generate new number with:                     !
!  call random_numer(rn)                         !
!  where where rn is a real(dp) variable         !
!================================================!
subroutine rand_init
  implicit none

  !**** LOCAL ***********************************!

  character(10) :: system_time

  integer :: i
  integer :: rand_size
  integer, allocatable :: rand_seed(:)

  real(dp) :: rtime
  real(dp) :: rn

  !**********************************************!

  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size))
  call date_and_time(time=system_time)
  read (system_time,*) rtime
  rand_seed=(Node+1)*int(rtime*1000.0_dp)
  call random_seed(put=rand_seed)
  deallocate(rand_seed)

  do i=1,10000
    call random_number(rn)
  end do

end subroutine rand_init

end module m_zminim
