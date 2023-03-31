!
! Copyright (C) 1996-2020	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2020, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! Module for correcting the density matrix for retaining a constant charge density
! The idea is to introduce several different schemes of charge corrections.

module ts_dq_m

  use precision, only: dp
  
  implicit none

  public 

  ! Method parameters for the charge-correction
  integer, save :: TS_DQ_METHOD = 0
  integer, parameter :: TS_DQ_METHOD_BUFFER = 1
  integer, parameter :: TS_DQ_METHOD_FERMI = 2
  real(dp), save :: TS_DQ_FACTOR = 0.8_dp

  real(dp), save :: TS_DQ_FERMI_TOLERANCE = 0.01_dp
  real(dp), save :: TS_DQ_FERMI_MAX = 0.1102471_dp ! 1.5 eV
  real(dp), save :: TS_DQ_FERMI_ETA = 7.349806700083787e-5_dp ! 0.001 eV
  real(dp), save :: TS_DQ_FERMI_SCALE = 75._dp


  type :: ts_dq_t
    !< Whether we should run the Ef charge correction
    logical :: run = .false.

    !< Values per chemical potential
    type(ts_dq_t_mu_), pointer :: mus(:) => null()

  contains

    procedure, pass :: initialize => ts_dq_t_initialize
    procedure, pass :: initialize_dq => ts_dq_t_initialize_dq
    procedure, pass :: calculate_dEf => ts_dq_t_calculate_dEf
    procedure, pass :: get_index => ts_dq_t_get_index
    procedure, pass :: delete => ts_dq_t_delete

  end type ts_dq_t

  !< Container for dq and eta *per* mu
  type :: ts_dq_t_mu_
    !< Charge at chemical potentials (create a list)
    real(dp), pointer :: dq(:) => null()
    !< Eta value at corresponding  chemical potentials (create a list)
    real(dp), pointer :: eta(:) => null()
    !< The linear total equilibrium index of the point
    integer, pointer :: idxE(:) => null()
  end type ts_dq_t_mu_

  type(ts_dq_t), public, save :: ts_dq

  private :: dp

contains

  subroutine ts_dq_read( )

    use fdf, only : fdf_get, leqi
    use units, only: eV
    character(len=64) :: chars

    chars = fdf_get('TS.ChargeCorrection', 'none')
    chars = fdf_get('TS.dQ',chars)
    TS_DQ_METHOD = 0
    if ( leqi(chars,'none') ) then
      TS_DQ_METHOD = 0
    else if ( leqi(chars,'b') .or. leqi(chars,'buffer') ) then
      TS_DQ_METHOD = TS_DQ_METHOD_BUFFER
    else if ( leqi(chars,'fermi') ) then
      TS_DQ_METHOD = TS_DQ_METHOD_FERMI
    else
      call die('TS.dQ: Charge correction method unknown, &
          &only one of [none|buffer|fermi] are allowed.')
    end if

    TS_DQ_FERMI_TOLERANCE = &
        fdf_get('TS.ChargeCorrection.Fermi.Tolerance', 0.01_dp)
    TS_DQ_FERMI_TOLERANCE = &
        fdf_get('TS.dQ.Fermi.Tolerance', TS_DQ_FERMI_TOLERANCE)
    ! Truncation of fermi-level change (default 1.5 eV)
    TS_DQ_FERMI_MAX = fdf_get('TS.ChargeCorrection.Fermi.Max', 1.5_dp * eV,'Ry')
    TS_DQ_FERMI_MAX = fdf_get('TS.dQ.Fermi.Max', TS_DQ_FERMI_MAX,'Ry')
    ! The scale of the dDmax and dHmax required to run the charge correction.
    ! Once dDmax and dHmax are below d*criteria * Scale we allow running the charge
    ! correction
    TS_DQ_FERMI_SCALE = fdf_get('TS.dQ.Fermi.Scale', 50._dp)
    ! The eta value at which to extrapolate the charge density
    TS_DQ_FERMI_ETA = fdf_get('TS.dQ.Fermi.Eta', 0.001_dp * eV, 'Ry')

    ! Factor for charge-correction
    TS_DQ_FACTOR = fdf_get('TS.ChargeCorrection.Factor', 0.8_dp)
    TS_DQ_FACTOR = fdf_get('TS.dQ.Factor', TS_DQ_FACTOR)
    if ( TS_DQ_FACTOR <= 0.0_dp ) then
      call die("TS.dQ.Factor: Charge correction factor must be larger than 0")
    end if

  end subroutine ts_dq_read


  !< Correct total charge such that the cell becomes charge neutral.
  !!
  !! This correction method works by moving excess charge to the buffer
  !! region. Since this region should not influence the TS region it
  !! should be relatively stable.
  !!
  !! It will however, influence the XC and Hartree energies
  !! for these regions.
  !! Currently, we can't correct the total energies related
  !! to these *in-correct* terms.
  subroutine ts_dq_buffer(N_Elec, &
      dit, sp, nspin, n_nzs, DM, EDM, S, Qtot)
    
    use m_ts_method
    use parallel, only : Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use geom_helper, only : UCORB

    use ts_charge_m, only: ts_charge_get
    
! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: N_Elec
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix and energy density matrix
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    ! The overlap
    real(dp), intent(in) :: S(n_nzs)
    ! Total charge of the system
    real(dp), intent(in) :: Qtot

! **********************
! * LOCAL variables    *
! **********************
    ! The charge in the regions
    real(dp), allocatable :: Q(:,:)
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_lo, no_u, lio, io, ind, ispin
    real(dp) :: reD
!    real(dp) :: addQ(nspin)

    ! Reset to zero if not existing
    if ( no_Buf == 0 ) return

    allocate(Q(0:2+N_Elec*2,nspin))
    
    call ts_charge_get(N_Elec, dit, sp, nspin, n_nzs, DM, S, Q = Q)

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
         n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=no_lo,nrows_g=no_u)
    
    ! Calculate the density factor for obtaining the correct charge
    ! For the left buffer region
    reD = (Qtot-sum(Q(:,:))) / sum(Q(1,:))

    ! Apply charge-correction factor 
    ! This will reduce "heavy" charge fluctuations and
    ! should guard against this.
    reD = reD * TS_DQ_FACTOR + 1._dp

    ! immediately deallocate charge
    deallocate(Q)

!    addQ(:) = 0.0_dp
    do ispin = 1 , nspin
       do lio = 1 , no_lo
          
          ! obtain the global index of the orbital.
          io = index_local_to_global(dit,lio,Node)

          if ( orb_type(io) /= TYP_BUFFER ) cycle

          ! Loop number of entries in the row... (index frame)
          do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

             if ( orb_type(l_col(ind)) /= TYP_BUFFER ) cycle

!             addQ(ispin) = addQ(ispin) + DM(ind,ispin) * (reD - 1._dp)

             DM(ind,ispin) = DM(ind,ispin) * reD
             ! Energy density matrices for buffer regions are not necessary
             ! to change, but we do it anyways
             EDM(ind,ispin) = EDM(ind,ispin) * reD

          end do
       end do
       
    end do

  end subroutine ts_dq_buffer

  !< Populate the ts_dq type for the arrays that it should read
  !!
  !! This gets populated according to the MPI distribution of the
  !! energy-points.
  !! So if the order of energy-points is changed in m_ts_*[gk]
  !! routines, then it should also be adapted here.
  subroutine ts_dq_t_initialize(this, N_mu, mus)
    
    use sorting, only: ordix
    use m_ts_chem_pot, only: ts_mu

    use m_ts_cctype, only: ts_c_idx
    use m_ts_contour_eq,  only : Eq_E, c2energy, get_c_io_index, Eq_c, Eq_linear_index

    class(ts_dq_t), intent(inout) :: this
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in) :: mus(N_mu)

    real(dp), parameter :: TOLERANCE_E = 0.0000001_dp
    type(ts_c_idx) :: cE
    real(dp), allocatable :: E(:)
    integer, allocatable :: indx(:), perm(:)
    integer :: imu, ic, idx, ie, N

    ! Clean if necessary
    if ( associated(this%mus) ) call this%delete()

    ! Allocate
    allocate(this%mus(N_mu))

    do imu = 1, N_mu
      ! Loop equilibrium sections for this chemical potential
      N = 0
      do ic = 1, size(mus(imu)%Eq_seg)
        idx = get_c_io_index(mus(imu)%Eq_seg(ic))

        ! Now check the energies lying close to the
        do ie = 1, size(Eq_c(idx)%c)
          if ( abs(real(Eq_c(idx)%c(ie), dp) - mus(imu)%mu) < TOLERANCE_E ) N = N + 1
        end do
      end do

      ! Pre-allocate for all so we can sort
      allocate(E(N), indx(N), perm(N))
      N = 0
      do ic = 1, size(mus(imu)%Eq_seg)
        idx = get_c_io_index(mus(imu)%Eq_seg(ic))

        ! Now check the energies lying close to the
        do ie = 1, Eq_c(idx)%c_io%N
          if ( abs(real(Eq_c(idx)%c(ie), dp) - mus(imu)%mu) < TOLERANCE_E ) then
            N = N + 1
            ! This ensures we can sort according to ETA value
            E(N) = abs(aimag(Eq_c(idx)%c(ie)) - TS_DQ_FERMI_ETA)
            indx(N) = Eq_linear_index(idx, ie)
          end if
        end do
      end do
      ! Now we have the full energy spectrum for energies close to the chemical potentials
      call ordix(E, 1, N, perm)

      allocate(this%mus(imu)%dq(min(5,N)))
      allocate(this%mus(imu)%eta(min(5,N)))
      allocate(this%mus(imu)%idxE(min(5,N)))

      do ie = 1, min(5, N)
        idx = indx(perm(ie))
        cE = Eq_E(idx)
        this%mus(imu)%eta(ie) = aimag(cE%e)
        this%mus(imu)%idxE(ie) = idx
      end do

      deallocate(E, indx, perm)

    end do

    ! At this point we will *maximally* have allocated
    !   N_mu * 5 * 3 arrays
    ! which for 5 different chemical potentials will be ~ 0.6 kB ~ 0 ;)
      
  end subroutine ts_dq_t_initialize

  !< Zero out the dq variables contained in this type
  subroutine ts_dq_t_initialize_dq(this)
    
    class(ts_dq_t), intent(inout) :: this
    integer :: i

    if ( .not. associated(this%mus) ) return

    do i = 1, size(this%mus)
      if ( associated(this%mus(i)%dq) ) then
        this%mus(i)%dq(:) = 0._dp
      end if
    end do
      
  end subroutine ts_dq_t_initialize_dq

  subroutine ts_dq_t_delete(this)
    class(ts_dq_t), intent(inout) :: this

    integer :: i

    ! Quick return if already cleaned
    if ( .not. associated(this%mus) ) return

    do i = 1, size(this%mus)
      if ( associated(this%mus(i)%dq) ) then
        deallocate(this%mus(i)%dq)
        deallocate(this%mus(i)%eta)
        deallocate(this%mus(i)%idxE)
      end if
    end do

    deallocate(this%mus)
    nullify(this%mus)

  end subroutine ts_dq_t_delete

  subroutine ts_dq_t_calculate_dEf(this, Qtot, dit, nspin, n_nzs, DM, S, dEf)

    use parallel, only: Node
    use units, only : eV
    use class_OrbitalDistribution
    use alloc, only: de_alloc
    use m_interpolate
#ifdef MPI
    use mpi_siesta
#endif

    !< The dq-type which contain charges calculated at individual chemical potentials
    class(ts_dq_t), intent(inout) :: this
    !< The target total charge
    real(dp), intent(in) :: Qtot
    type(OrbitalDistribution), intent(inout) :: dit
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrices and overlap
    real(dp), intent(in) :: DM(n_nzs,nspin), S(n_nzs)
    real(dp), intent(out) :: dEf

! ******************* Local arrays *******************
    real(dp) :: dQ, dQ_Ef, tQ
    integer :: imu, ind, N_dq
#ifdef MPI
    real(dp), allocatable :: tmp(:)
    integer :: comm
    integer :: MPIerror
#endif

    ! Return if not runned
    if ( .not. this%run ) return

    ! Calculate charge for current density matrix
    dQ = 0._dp
!$OMP parallel do default(shared), private(ind), reduction(+:dQ)
    do ind = 1, n_nzs
      dQ = dQ + sum(DM(ind,:)) * S(ind)
    end do
!$OMP end parallel do

#ifdef MPI
    ! Reductions etc.
    comm = dist_comm(dit)

    ! Now reduce all entries
    do imu = 1, size(this%mus)
      N_dq = size(this%mus(imu)%dq)
      allocate(tmp(N_dq))
      call MPI_Reduce(this%mus(imu)%dq(1), tmp(1), N_dq, &
          MPI_Double_Precision, MPI_SUM, 0, comm, MPIerror)
      if ( Node == 0 ) then
        this%mus(imu)%dq(:) = tmp(:)
      end if
      deallocate(tmp)
    end do
    
    call MPI_Reduce(dQ, dQ_Ef, 1, MPI_Double_Precision, &
        MPI_SUM, 0, comm, MPIerror)
    dQ = dQ_Ef
#endif

    if ( Node == 0 ) then
      ! The additional charge
      dQ = dQ - Qtot
      write(*,'(a,es11.4)') 'ts-dq: dq = ', dQ

      dQ_Ef = 0._dp
      
      do imu = 1, size(this%mus)
        N_dq = size(this%mus(imu)%dq)
        ! Interpolate the charge at the given eta value
        call interp_spline(N_dq, this%mus(imu)%eta, this%mus(imu)%dq, TS_DQ_FERMI_ETA, tQ)
        dQ_Ef = dQ_Ef + tQ
      end do

      ! Average contribution from each chemical potential
      dQ_Ef = dQ_Ef / size(this%mus)

      ! Now we have the difference in Q
      ! Correct the Fermi level so that a change dE DM would 
      ! account for the missing/excess charge.
      ! dQ + dE * DM@(Ef) = 0 => dE = -dQ / DM@(Ef)
      dEf = - (dQ / dQ_Ef) * TS_DQ_FACTOR

      ! If Ef lies in the middle of bands we will have no DOS
      ! right at the Fermi level.
      ! If this is the case we truncate the change in Fermi-level
      ! to the maximum allowed shift...
      call ts_dq_truncate(0._dp, TS_DQ_FERMI_MAX, dEf)

    end if

#ifdef MPI
    ! B-cast to other nodes
    call MPI_Bcast(dEf, 1, MPI_Double_Precision, 0, &
        MPI_Comm_World, MPIerror)
#endif
    
  end subroutine ts_dq_t_calculate_dEf


  function ts_dq_t_get_index(this, imu, idxE) result(index_dq)
    class(ts_dq_t), intent(inout) :: this
    integer, intent(in) :: imu, idxE
    integer :: index_dq

    integer :: i

    if ( this%run ) then
      do i = 1, size(this%mus(imu)%idxE)
        if ( this%mus(imu)%idxE(i) == idxE ) then
          index_dq = i
          return
        end if
      end do
    end if

    index_dq = 0

  end function ts_dq_t_get_index
  
  subroutine ts_dq_Fermi_file(Ef, dEf)

    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    use parallel, only : Node
    use units, only : eV
    use m_interpolate
#ifdef MPI
    use mpi_siesta
#endif
    !< Ef is the -dq / dq@Ef corrected Fermi-level
    real(dp), intent(inout) :: Ef
    !< The shift in the Fermi-level for this routine
    real(dp), intent(inout) :: dEf
    
    integer :: iu, ioerr, i
    real(dp) :: cur(2)

    real(dp) :: Ef_new
    real(dp), allocatable :: Q_Ef(:,:), first_Q(:)
    
    integer :: N, max_itt, tmp
    integer :: i_start
    character(len=2) :: char2
    logical :: is_pos

#ifdef MPI
    integer :: MPIerror
#endif

    dEf = 0._dp

    ! First we need to read in the entire 
    ! file and figure out the highest number of iterations
    ! performed.
    if ( Node == 0 ) then

      max_itt = 0
      N = 0
      ! open file
      call io_assign(iu)
      open(unit=iu, file='TS_FERMI', form='formatted', status='old', &
          iostat=ioerr)
      if ( ioerr /= 0 ) then
        call die('The file has not been created, this should not happen')
      end if
      rewind(iu)
       
      ! Figure out number of fermi iterations
      do while ( ioerr /= -1 ) 
        N = N + 1
        read(iu,*) ! # TSiscf line
        read(iu,'(a2,i15)') char2, tmp
        max_itt = max(max_itt,tmp)
        do i = 1 , tmp
          read(iu,*) ! data line
        end do
        read(iu,*,iostat=ioerr) ! empty line
      end do

      ! If only one data point is present,
      ! we cannot do anything...
      if ( N > 1 ) then

        ! First calculated charge for every iteration
        allocate(first_Q(N))
        ! Q_Ef(iteration, [Ef, Q])
        allocate(Q_Ef(N,2))

        ! Rewind and read data
        rewind(iu)
        ioerr = 0

        ! Read in the data and move it to Q_Ef
        N = 0
        do while ( ioerr /= -1 )
          N = N + 1
          read(iu,*) ! # TSiscf line
          read(iu,'(a2,i15)') char2, tmp

          read(iu,'(2(tr1,e20.10))') cur(:)
          first_Q(N) = cur(2)
          ! Convert to Ry
          cur(1) = cur(1) * eV
          Q_Ef(N,:) = cur(:)
          do i = 2 , tmp
            read(iu,*) !
          end do
          read(iu,*,iostat=ioerr) ! empty line

        end do

        ! Figure out where we should use it
        is_pos = is_positive(first_Q(N))
        do i_start = N , 2, -1
          if ( is_pos .neqv. is_positive(first_Q(i_start-1)) ) exit
        end do
        i_start = max(i_start, 1)

        if ( N - i_start > 1 ) then ! len(i_start:N) >= 2
          ! Interpolate the new fermi level by using first entry
          call interp_spline(N-i_start+1,Q_Ef(i_start:N,2),Q_Ef(i_start:N,1),0._dp,Ef_new)
          if ( IEEE_IS_NAN(Ef_new) ) Ef_new = Ef
          dEf = TS_DQ_FACTOR * (Ef_new - Ef)
        end if

        deallocate(Q_Ef, first_Q)

        ! Truncate to the maximum allowed difference
        call ts_dq_truncate(0._dp, TS_DQ_FERMI_MAX, dEf)

      end if ! N > 1

      call io_close(iu)

    end if

#ifdef MPI
    call MPI_Bcast(dEf, 1, MPI_Double_Precision, 0, &
        MPI_Comm_World, MPIerror)
#endif

  contains

    pure function is_positive(V) result(is)
      real(dp), intent(in) :: V
      logical :: is
      is = V > 0._dp
    end function is_positive

  end subroutine ts_dq_Fermi_file

  subroutine ts_dq_truncate(Ef, max_diff, Ef_new, truncated)
    real(dp), intent(in) :: Ef, max_diff
    real(dp), intent(inout) :: Ef_new
    logical, intent(out), optional :: truncated

    if ( abs(Ef_new - Ef) > max_diff ) then
      if ( present(truncated) ) truncated = .true.
      if ( Ef_new - Ef > 0._dp ) then
        Ef_new =   max_diff + Ef
      else
        Ef_new = - max_diff + Ef
      end if
    else if ( present(truncated) ) then
      truncated = .false.
    end if

  end subroutine ts_dq_truncate

  subroutine ts_dq_scale_EDM_elec(N_Elec, Elecs, dit, sp, nspin, n_nzs, DM, EDM, &
      dEf_elec, dEf_device)

    use m_ts_method
    use parallel, only : Node
#ifdef MPI
    use mpi_siesta
#endif
    use class_OrbitalDistribution
    use class_Sparsity
    use m_ts_electype

! **********************
! * INPUT variables    *
! **********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    type(OrbitalDistribution), intent(inout) :: dit
    ! SIESTA local sparse pattern (not changed)
    type(Sparsity), intent(inout) :: sp
    ! Number of non-zero elements
    integer, intent(in) :: nspin, n_nzs
    ! The density matrix
    real(dp), intent(in) :: DM(n_nzs,nspin)
    ! The energy density matrix
    real(dp), intent(inout) :: EDM(n_nzs,nspin)
    !< Fermi-level changes in the electrode region
    real(dp), intent(in) :: dEf_elec, dEf_device

! **********************
! * LOCAL variables    *
! **********************
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: no_l, no_u, lio, io, ind
    integer :: ir_dm, jr_dm
    integer :: ir, jr

    ! Retrieve information about the sparsity pattern
    call attach(sp, &
        n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
        nrows=no_l,nrows_g=no_u)

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ind,ir,ir_dm,jr,jr_dm)
    do lio = 1 , no_l

      ! obtain the global index of the orbital.
      io = index_local_to_global(dit,lio,Node)
      ir = orb_type(io)
      if ( ir == TYP_BUFFER ) then
        ir_dm = 0
      else if ( ir == TYP_DEVICE ) then
        ir_dm = 2
      else
        ir_dm = Elecs(ir)%DM_update
      end if

      ! Loop number of entries in the row... (index frame)
      do ind = l_ptr(lio) + 1 , l_ptr(lio) + l_ncol(lio)

        ! as the local sparsity pattern is a super-cell pattern,
        ! we need to check the unit-cell orbital
        ! The unit-cell column index
        jr = orb_type(l_col(ind))
        if ( jr == TYP_BUFFER ) then
          jr_dm = 0
        else if ( jr == TYP_DEVICE ) then
          jr_dm = 2
        else
          jr_dm = Elecs(jr)%DM_update
        end if

        if ( ir_dm + jr_dm <= 2 ) then
          ! In a non-updated region
          EDM(ind,:) = EDM(ind,:) + DM(ind,:) * dEf_elec
        else
          ! In an updated region
          EDM(ind,:) = EDM(ind,:) + DM(ind,:) * dEf_device
        end if

      end do
    end do
!$OMP end parallel do

  end subroutine ts_dq_scale_EDM_elec

end module ts_dq_m
