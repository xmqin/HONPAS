!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

module m_transiesta

  use precision, only : dp

  use files, only : slabel

  use m_ts_sparse, only : ts_sparse_init
  use m_ts_method, only : ts_method
  use m_ts_method, only : ts_A_method, TS_BTD_A_COLUMN

  use m_ts_method, only : TS_FULL, TS_BTD
#ifdef SIESTA__MUMPS
  use m_ts_method, only : TS_MUMPS
#endif

  use m_ts_tri_init, only : ts_tri_init

  use m_ts_fullg
  use m_ts_fullk

  use m_ts_trig
  use m_ts_trik

#ifdef SIESTA__MUMPS
  use m_ts_mumpsg
  use m_ts_mumpsk
#endif

  implicit none

  public :: transiesta
  private

contains

  subroutine transiesta(TSiscf,nspin, &
       sp_dist, sparse_pattern, &
       no_aux_cell, ucell, nsc, isc_off, no_u, na_u, lasto, xa, n_nzs, &
       H, S, DM, EDM, Ef, &
       Qtot)

    use units, only : eV
    use alloc, only : re_alloc, de_alloc

    use parallel, only : IONode
    use m_os, only: file_delete

    use class_OrbitalDistribution
    use class_Sparsity

    use m_ts_kpoints, only : ts_nkpnt, ts_Gamma

    use m_ts_electype

    use m_ts_options, only : N_Elec, Elecs
    use m_ts_options, only : IsVolt, Calc_Forces

    use m_ts_options, only : BTD_method

    use m_ts_contour_eq , only : N_Eq_E
    use m_ts_contour_neq, only : N_nEq_E

    use ts_charge_m
    use ts_dq_m
    
    use m_ts_gf, only : read_Green
    use m_interpolate

    use m_energies, only: NEGF_DE

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in)  :: TSiscf
    integer, intent(in)  :: nspin
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    logical, intent(in)  :: no_aux_cell
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in)  :: nsc(3), no_u, na_u
    integer, intent(in) :: isc_off(3,product(nsc))
    integer, intent(in)  :: lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: H(n_nzs,nspin), S(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Qtot
    real(dp), intent(inout) :: Ef

! ******************** IO descriptors ************************
    integer, allocatable :: uGF(:)
! ************************************************************

! ****************** Electrode variables *********************
    integer, allocatable :: nq(:)
! ************************************************************

    ! * local variables
    integer :: iEl, NEn, no_used, no_used2
    logical :: converged, dq_run, truncated
    ! In case of Fermi-correction, we save the previous steps
    ! and do a spline interpolation... :)
    integer :: N_F, i_F, ioerr
    real(dp) :: dEf
    real(dp), pointer :: Q_Ef(:,:) => null()

    if ( nspin > 2 ) then
      call die('transiesta does not work for non-collinear or spin-orbit')
    end if

    ! Open GF files...
    ! Read-in header of Green functions
    ! Prepare for the calculation
    ! We read in the k-points that the electrode was generated with.
    ! Furthermore we read in the expansion q-points
    ! They are communicated in the routine

    ! Initialize the DE_NEGF energy
    ! Note that this will *only* be updated in case V /= 0.
    ! This energy corresponds to the non-equilibrium energies:
    !   e \sum_i N_i * \mu_i
    NEGF_DE = 0._dp

    if ( TSiscf == 1 ) then
      ! We need to initialize TRANSIESTA
      
      call timer('TS_init',1)

      if ( IONode ) then
        if ( file_delete('TS_FERMI') ) then
          write(*,'(/a/)') 'ts: Deleted TS_FERMI'
        end if
      end if
      
      call ts_sparse_init(slabel, IsVolt, N_Elec, Elecs, &
          ucell, nsc, na_u, xa, lasto, sp_dist, sparse_pattern, no_aux_cell, &
            isc_off)
      
      if ( ts_method == TS_BTD ) then
        ! initialize the tri-diagonal partition
        call ts_tri_init( sp_dist, sparse_pattern , N_Elec, &
            Elecs, IsVolt, ucell, na_u, xa, lasto ,nsc, isc_off, &
            BTD_method )
      end if

      ! print out estimated memory usage...
      call ts_print_memory(ts_Gamma)

      call ts_charge_print(N_Elec,Elecs, Qtot, sp_dist, sparse_pattern, &
          nspin, n_nzs, DM, S, TS_Q_INFO_FULL)
      
      if ( .not. Calc_Forces ) then
        if ( IONode ) then
          write(*,'(a)') 'transiesta: *** The forces are NOT updated ***'
          write(*,'(a)') 'transiesta: ***  Will set the energy-density matrix to 0!  ***'
        end if
        EDM(:,:) = 0._dp
      end if

      call timer('TS_init',2)
      
    end if

    
    call timer('TS',1)

    
    ! Total number of energy-points...
    NEn = N_Eq_E() + N_nEq_E()

    ! in case the file-descriptor is negative it basically 
    ! means "out-of-core" calculation.
    allocate(uGF(N_Elec),nq(N_Elec))
    uGF(:) = -1
    do iEl = 1 , N_Elec

       ! Calculate number of Bloch expansion k-points
       nq(iEl) = Elecs(iEl)%Bloch%size()

       ! Allocate the electrode quantities
       nullify(Elecs(iEl)%HA,Elecs(iEl)%SA,Elecs(iEl)%Gamma)

       ! We allocate for once as much space as needed,

       ! Allocate the non-repeated hamiltonian and overlaps...
       no_used = Elecs(iEl)%no_used
       if ( Elecs(iEl)%pre_expand > 1 ) then ! > 1 also expand H, S before writing
          no_used = TotUsedOrbs(Elecs(iEl))
          nq(iEl) = 1
       end if

       if ( IsVolt .or. .not. Elecs(iEl)%Bulk ) then
          ! If we using bulk electrodes, we need not the Hamiltonian, 
          ! nor the overlap...
          call re_alloc(Elecs(iEl)%HA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
          call re_alloc(Elecs(iEl)%SA,1,no_used,1,no_used,1,nq(iEl),routine='transiesta')
       end if

       no_used = TotUsedOrbs(Elecs(iEl))
       if ( IsVolt ) then
          ! We need Gamma's with voltages (now they are both GAA and GammaT)
          no_used2 = no_used
       else 
          ! This is only for having space for GA
          if ( Elecs(iEl)%pre_expand > 0 ) then
             no_used2 = no_used
          else
             no_used2 = Elecs(iEl)%no_used
          end if
       end if
       call re_alloc(Elecs(iEl)%Gamma,1,no_used*no_used2,routine='transiesta')

       ! This seems stupid, however, we never use the expansion array and
       ! GammaT at the same time. Hence it will be safe
       ! to have them point to the same array.
       ! When the UC_expansion_Sigma_GammaT is called:
       ! first the GAA is "emptied" of information and then
       ! Gamma is filled.
       if ( Elecs(iEl)%pre_expand == 0 ) no_used2 = Elecs(iEl)%no_used
       Elecs(iEl)%GA => Elecs(iEl)%Gamma(1:no_used*no_used2)

    end do

    ! start calculation
    converged = .false.
    i_F = -1
    dq_run = ts_dq%run
    if ( dq_run ) then

      ! we will utilize the old Fermi-level to correct the 
      ! EDM matrix (just in case the 
      ! electrode region elements are not taken care of)
      
      ! Allocate for interpolation, this is typically
      ! converging in 3-5 steps
      N_F = 10
      i_F = 0
      call re_alloc(Q_Ef,1,N_F,1,2)
    end if

    do while ( .not. converged )
      
      call open_GF(N_Elec,Elecs,uGF,NEn)

      select case ( ts_method )
      case ( TS_FULL )
        if ( ts_Gamma ) then
          call ts_fullg(N_Elec,Elecs, &
              nq, uGF, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        else
          call ts_fullk(N_Elec,Elecs, &
              nq, uGF, &
              ucell, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        end if
      case ( TS_BTD )
        if ( ts_Gamma ) then
          call ts_trig(N_Elec,Elecs, &
              nq, uGF, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        else
          call ts_trik(N_Elec,Elecs, &
              nq, uGF, &
              ucell, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        end if
#ifdef SIESTA__MUMPS
      case ( TS_MUMPS )
        if ( ts_Gamma ) then
          call ts_mumpsg(N_Elec,Elecs, &
              nq, uGF, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        else
          call ts_mumpsk(N_Elec,Elecs, &
              nq, uGF, &
              ucell, nspin, na_u, lasto, &
              sp_dist, sparse_pattern, &
              no_u, n_nzs, &
              H, S, DM, EDM, Ef, NEGF_DE)
        end if
#endif
      case default
        call die('Error in code')
      end select
      
      ! Close files
      do iEl = 1 , N_Elec
        if ( IONode .and. Elecs(iEl)%out_of_core ) then
          call io_close(uGF(iEl))
        end if
      end do
      
      if ( dq_run ) then

        ! This is the 1. part of the Fermi correction
        !
        ! In this section we estimate the new Fermi level by
        ! calculating the charge Q(E_F) and then correct
        ! E_F such that dQ = Q(TS) - Qtot -> 0.
        i_F = i_F + 1
        if ( N_F < i_F ) then
          N_F = N_F + 10
          call re_alloc(Q_Ef, 1, N_F, 1, 2, copy=.true.)
        end if

        ! Save current fermi level and charge
        call ts_charge_get(N_Elec, sp_dist, sparse_pattern, &
            nspin, n_nzs, DM, S, Qtot = Q_Ef(i_F,1) )
        Q_Ef(i_F,2) = Ef

        ! Determine whether the charge *has* converged
        converged = abs(Q_Ef(i_F,1) - Qtot) < TS_DQ_FERMI_TOLERANCE

        if ( i_F < 2 ) then

          call ts_dq%calculate_dEf(Qtot, &
              sp_dist, nspin, n_nzs, DM, S, dEf)
          ! In case the charge *is* converged
          ! Then we further reduce dEf since
          ! it causes SCF fluctuations
          if ( converged ) then
            Ef = Ef + dEf * min(abs(Q_Ef(i_F,1) - Qtot), 0.05_dp)
          else
            Ef = Ef + dEf
          end if

        else

          ! Instead of doing Q(E_F) for every change
          ! we will perform spline interpolation ones we have 2 estimated
          ! Q(E_F).
          ! This tends to drastically speed up the convergence of the dQ -> 0.

          ! In case we have accumulated 2 or more points
          dEf = Ef
          call interp_spline(i_F,Q_Ef(1:i_F,1),Q_Ef(1:i_F,2),Qtot,Ef)
          dEf = Ef - dEf

          ! Truncate to the maximum allowed change in Fermi-level
          call ts_dq_truncate(Q_Ef(i_F,2), &
              TS_DQ_FERMI_MAX, Ef, truncated=truncated)
          converged = converged .and. (.not. truncated)

          if ( IONode ) then
            write(*,'(a,es11.4)') 'ts-dq-iscf: cubic spline dq = ', &
                Q_Ef(i_F,1) - qtot
          end if

        end if

        ! Signal only doing this for the first run
        ! for subsequent dq-scf we interpolate
        ts_dq%run = .false.

      else

        ! If no Fermi-correction, we are converged
        converged = .true.

      end if

    end do

    ! We will only do major extrapolation when we are
    ! far from the charge neutrality point.
    ! Otherwise we will iterate forever...
    if ( dq_run ) then

      ! Store the actual number of charge iterations runned
      N_F = i_F

      if ( N_F > 1 ) then

        if ( IONode ) then
          ! After converge we write out the convergence
          call io_assign(iEl)
          inquire(file='TS_FERMI', exist=converged)
          if ( converged ) then
            open(unit=iEl,file='TS_FERMI',position='append',form='formatted', &
                status='old',iostat=ioerr)
            write(iEl,'(/,a,i0)') '# TSiscf = ', TSiscf
          else
            open(unit=iEl,file='TS_FERMI',form='formatted', &
                status='new')
            write(iEl,'(a,i0)') '# TSiscf = ', TSiscf
          end if
          write(iEl,'(a,i0)')'# ', N_F ! Number of iterations
          do i_F = 1 , N_F
            write(iEl,'(2(tr1,e20.10))') Q_Ef(i_F,2)/eV, Q_Ef(i_F,1) - Qtot
          end do

          call io_close(iEl)

        end if

        ! When N_F > 1 we *know* that the initial run did
        ! not converge the charges

        ! This is the 2nd step of dEF correction.
        ! At this point we have corrected E_f for the current
        ! iteration. But generally the Hartree potential will counter
        ! the change in E_F. So to speed up convergence
        ! we do a spline interpolation of the dE_F by doing a spline
        ! interpolation of the ISCF corrections.
        ! Say TS corrects EF at iterations 50 and 80
        ! which means TS_FERMI may look like this:
        !#####
        ! # TSiscf = 50
        ! # 3
        !   -0.204782E+01    0.172948E-01
        !   -0.205297E+01    0.834017E-03
        !   -0.205323E+01   -0.250016E-05
        !
        ! # TSiscf = 80
        !# 3
        !   -0.207423E+01    0.967930E-02
        !   -0.207710E+01    0.490200E-03
        !   -0.207726E+01   -0.817259E-06
        !#####
        !

        ! Guess-stimate the actual Fermi-shift
        ! typically will the above be "too" little
        ! So we interpolate between all previous
        ! estimations for this geometry...
        call ts_dq_Fermi_file(Ef, dEf)
        Ef = Ef + dEf

      end if

      ! We have now calculated the new Ef
      if ( IONode ) then
        write(*,'(a,es11.4,a)') 'ts-dq: dEf = ', (Ef - Q_Ef(1,2))/eV, ' eV'
      end if

      ! We shift it EDM to the correct level
      ! Only correct energy density matrix in the respective
      ! places. For instance the above
      ! Only correct for the latest one, which have not
      ! been applied during the charge-convergence ISCF
      ! Although the DM/EDM for the device region
      ! have been calculated with the
      ! fermi level stored in Q_Ef(N_F,2) it seems the forces
      ! in the electrode regions are better if one uses
      ! the converged Ef as the baseline
      ! I.e. if one wants EDM in the electrode regions
      ! to be commensurate with the EDM in the device region
      ! one should use Q_Ef(N_F,2) - Q_Ef(1,2)
      call ts_dq_scale_EDM_elec(N_elec, Elecs, &
          sp_dist, sparse_pattern, nspin, n_nzs, DM, EDM, &
          Ef - Q_Ef(1,2), Ef - Q_Ef(N_F,2))

      call de_alloc(Q_Ef)

    end if

    !***********************
    !       Clean up
    !***********************
    do iEl = 1 , N_Elec
      if ( .not. Elecs(iEl)%out_of_core ) then
        call delete(Elecs(iEl))
      end if
    end do

    !***********************
    !  Clean up electrodes
    !***********************
    do iEl = 1 , N_Elec
      if ( associated(Elecs(iEl)%HA) ) then
        call de_alloc(Elecs(iEl)%HA,routine='transiesta')
        call de_alloc(Elecs(iEl)%SA,routine='transiesta')
      end if
      call de_alloc(Elecs(iEl)%Gamma,routine='transiesta')
    end do

    deallocate(uGF,nq)

    ! We do the charge correction of the transiesta
    ! computation here (notice that the routine will automatically
    ! return if no charge-correction is requested)
    if ( TS_DQ_METHOD == TS_DQ_METHOD_BUFFER ) then
      call ts_dq_buffer(N_Elec, sp_dist, &
          sparse_pattern, nspin, n_nzs, DM, EDM, S, Qtot)
    end if

    call ts_charge_print(N_Elec, Elecs, Qtot, sp_dist, sparse_pattern, &
        nspin, n_nzs, DM, S, &
        method = TS_Q_INFO_SCF)

    call timer('TS',2)

  contains

    subroutine init_Electrode_HS(El)
      use class_Sparsity
      use class_dSpData1D
      use class_dSpData2D
      use alloc, only : re_alloc
      type(Elec), intent(inout) :: El
      
      ! If already initialized, return immediately
      if ( initialized(El%sp) ) return

      ! Read-in and create the corresponding transfer-matrices
      call delete(El) ! ensure clean electrode
      call read_Elec(El,Bcast=.true., IO = .false.)
      
      if ( .not. associated(El%isc_off) ) then
        call die('An electrode file needs to be a non-Gamma calculation. &
            &Ensure at least two k-points in the T-direction.')
      end if
      
      call create_sp2sp01(El, IO = .false.)

      ! Clean-up, we will not need these!
      ! we should not be very memory hungry now, but just in case...
      call delete(El%H)
      call delete(El%S)
      
      ! We do not accept onlyS files
      if ( .not. initialized(El%H00) ) then
        call die('An electrode file must contain the Hamiltonian')
      end if

      call delete(El%sp)

    end subroutine init_Electrode_HS

    subroutine open_GF(N_Elec,Elecs,uGF,NEn)
      integer, intent(in) :: N_Elec
      type(Elec), intent(inout) :: Elecs(N_Elec)
      integer, intent(out) :: uGF(N_Elec)
      integer, intent(in) :: NEn

      ! Local variables
      integer :: iEl
      
      do iEl = 1 , N_Elec

        ! Initialize k-points (never seen k-point)
        Elecs(iEl)%bkpt_cur(:) = 2352345._dp

        if ( Elecs(iEl)%out_of_core ) then

          if ( IONode ) then
            call io_assign(uGF(iEl))
            open(file=Elecs(iEl)%GFfile,unit=uGF(iEl),form='unformatted')
          end if

        else

          ! prepare the electrode to create the surface self-energy
          call init_Electrode_HS(Elecs(iEl))

        end if

        if ( Elecs(iEl)%out_of_core ) then
          call read_Green(uGF(iEl),Elecs(iEl), ts_nkpnt, NEn )
        end if
         
      end do
      
    end subroutine open_GF

  end subroutine transiesta

  subroutine ts_print_memory(ts_Gamma)
    
    use parallel, only : IONode
    use precision, only : i8b

#ifdef MPI
    use mpi_siesta, only : MPI_Comm_World
    use mpi_siesta, only : MPI_Max
    use mpi_siesta, only : MPI_Double_Precision
#endif 

    use class_Sparsity
    use m_ts_options, only : IsVolt, Calc_Forces
    use m_ts_options, only : N_mu, N_Elec, Elecs
    use m_ts_contour_neq, only : N_nEq_id
    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc, ltsup_sp_sc
    use m_ts_electype

    use m_ts_tri_init, only : c_Tri
    use m_ts_tri_common, only : GFGGF_needed_worksize
    use m_ts_tri_common, only : nnzs_tri_i8b
    use m_ts_method, only : no_Buf

    logical, intent(in) :: ts_Gamma ! transiesta Gamma
    integer :: i, no_E, no_used
    integer(i8b) :: nel
    integer :: padding, worksize
    real(dp) :: mem, dmem, zmem
#ifdef MPI
    integer :: MPIerror
#endif

    ! Estimate electrode sizes
    zmem = 0._dp
    do i = 1 , N_Elec

      no_used = Elecs(i)%no_used
      no_E = TotUsedOrbs(Elecs(i))

      if ( IsVolt .or. .not. Elecs(i)%Bulk ) then
        ! Hamiltonian and overlap
        if ( Elecs(i)%pre_expand > 1 ) then
          zmem = zmem + no_E ** 2 * 2
        else
          zmem = zmem + no_E * no_used * 2
        end if
      end if

      if ( IsVolt ) then
        zmem = zmem + no_E ** 2 ! GS/Gamma
      else
        if ( Elecs(i)%pre_expand > 0 ) then
          zmem = zmem + no_E ** 2
        else
          zmem = zmem + no_E * no_used
        end if
      end if
       
    end do
    zmem = zmem * 16._dp / 1024._dp ** 2
    if ( IONode ) then
      write(*,'(/,a,t55,f10.2,a)') &
          'transiesta: mem of electrodes (static): ', zmem,'MB'
    end if
    mem = zmem

    ! Global arrays
    dmem = 0._dp
    zmem = 0._dp
    nel = nnzs(ts_sp_uc) * 2
    if ( ts_Gamma ) then
      dmem = dmem + nel
    else
      zmem = zmem + nel
    end if

    ! global sparsity update
    nel = nnzs(tsup_sp_uc)
    if ( Calc_Forces ) then
      i = max(N_mu,N_nEq_id) + N_mu
    else
      i = max(N_mu,N_nEq_id)
    end if
    if ( ts_Gamma ) then
      dmem = dmem + nel * i
    else
      zmem = zmem + nel * i
    end if
    ! Convert to MB
    dmem = dmem * 8._dp / 1024._dp ** 2
    zmem = zmem * 16._dp / 1024._dp ** 2
    mem = mem + dmem + zmem

    if ( IONode ) then
      write(*,'(a,t55,f10.2,a)') &
          'transiesta: mem of global update arrays (static): ', dmem+zmem,'MB'
    end if

    ! Local sparsity update
    if ( IsVolt ) then
      nel = nnzs(ltsup_sp_sc)
      if ( Calc_Forces ) then
        dmem = nel * ( 2 * N_mu + N_nEq_id )
      else
        dmem = nel * ( N_mu + N_nEq_id )
      end if
      ! Bias local sparsity pattern is always
      ! in double precision
      dmem = dmem * 8._dp / 1024._dp ** 2
      if ( IONode ) then
        write(*,'(a,t55,f10.2,a)') &
            'transiesta: mem of master node sparse arrays: ', dmem,'MB'
      end if
      mem = mem + dmem
    end if

    if ( ts_method == TS_BTD ) then

      ! initialize padding and work-size query
      padding = 0
      worksize = 0

      if ( ts_A_method == TS_BTD_A_COLUMN .and. IsVolt ) then
        ! Calculate size of the tri-diagonal matrix
        call GFGGF_needed_worksize(c_Tri%n,c_Tri%r, &
            N_Elec, Elecs, padding, worksize)
      end if

      nel = nnzs_tri_i8b(c_Tri%n,c_Tri%r)
      if ( nel > huge(1) ) then
        call die('transiesta: Memory consumption is too large!')
      end if
      zmem = (nel * 2._dp + padding + worksize ) * 16._dp / 1024._dp ** 2
      if ( IONode ) &
          write(*,'(a,t55,f10.2,a)') &
          'transiesta: mem of tri-diagonal matrices: ', zmem,'MB'
      mem = mem + zmem

    else if ( ts_method == TS_FULL ) then
      ! Calculate size of the full matrices
      ! Here we calculate number of electrodes not needed to update the cross-terms
      no_E = sum(TotUsedOrbs(Elecs),Elecs(:)%DM_update==0)
      i = nrows_g(ts_sp_uc) - no_Buf
      ! LHS
      zmem = i ** 2
      ! RHS
      if ( IsVolt ) then
        zmem = zmem + i * max(i-no_E,sum(TotUsedOrbs(Elecs)))
      else
        zmem = zmem + i * (i-no_E)
      end if
      zmem = zmem * 16._dp / 1024._dp ** 2
      if ( IONode ) &
          write(*,'(a,t55,f10.2,a)') &
          'transiesta: mem of full matrices: ', zmem,'MB'
      mem = mem + zmem
#ifdef SIESTA__MUMPS
    else if ( ts_method == TS_MUMPS ) then
      if ( IONode ) then
        write(*,'(a)')'transiesta: mem is determined by MUMPS.'
        write(*,'(a)')'transiesta: Search in TS_MUMPS_<Node>.dat for: ### Minimum memory.'
      end if
#endif
    end if

#ifdef MPI
    call MPI_Reduce(mem,zmem,1,MPI_Double_Precision, &
        MPI_MAX, 0, MPI_Comm_World, MPIerror)
#else
    zmem = mem
#endif

    if ( IONode ) then
      write(*,'(a,t55,f10.2,a)') &
          'transiesta: Total memory usage: ', zmem,'MB'
    end if

  end subroutine ts_print_memory

end module m_transiesta

