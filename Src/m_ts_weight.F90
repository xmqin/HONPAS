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

! A module for doing different weighting schemes for the Transiesta voltage calculations
! *Notice* that for the Gamma-point they are all equivalent!

! In the following we describe the weights by:
!  w_L = the weight for the left correction term
!  w_R = the weight for the right correction term
!  w   = w_L / ( w_L + w_R )

! In order to accomodate several methods and precisions we introduce the following 3 methods
! 1)
!   do a single weight for the full correction term, i.e.:
!     w_L = (\sum_k \Delta_L(k))^2
!     w_R = (\sum_k \Delta_R(k))^2
!   this will assume that the k-point correction terms are correlated in some way

! 2)
!   do a weight for each k
!     w_L^k = \Delta_L^2(k)
!     w_R^k = \Delta_R^2(k)
!   this will actually save us more memory, as we don't need the full sum of the k-points.

module m_ts_weight

  use precision, only: dp

  implicit none

  private :: dp

  ! Generic methods for the k-weighting method
  ! Method: 1)
  integer, parameter :: TS_W_K_CORRELATED = 1
  ! Method: 2)
  integer, parameter :: TS_W_K_UNCORRELATED = 2

  ! General weighting method
  ! 1) orb-orb weighting
  integer, parameter :: TS_W_ORB_ORB = 1
  ! 2) atom-atom weighting using the trace of the nEq DM
  integer, parameter :: TS_W_TR_ATOM_ATOM = 2
  ! 3) atom-atom weighting using the sum of the nEq DM
  integer, parameter :: TS_W_SUM_ATOM_ATOM = 3
  ! 4) atom-orb weighting using the trace of the nEq DM
  integer, parameter :: TS_W_TR_ATOM_ORB = 4
  ! 5) atom-orb weighting using the sum of the nEq DM
  integer, parameter :: TS_W_SUM_ATOM_ORB = 5
  ! 6) Simple mean of contributions
  integer, parameter :: TS_W_MEAN = 6

  ! The general weighting can be:
  !   UNCORRELATED or CORRELATED
  ! The correlated is the default
  integer, parameter :: TS_W_CORRELATED = 100

  integer, save :: TS_W_METHOD = TS_W_ORB_ORB

  ! we default weight of uncorrelated as that should be 
  ! the most accurate
  integer, save :: TS_W_K_METHOD = TS_W_K_UNCORRELATED

contains

  subroutine read_ts_weight( )

    use fdf, only : fdf_get, leqi
    character(len=200) :: chars
    integer :: i

    ! Update the weight function
    chars = fdf_get('TS.Weight.k.Method','correlated')
    if ( leqi(chars,'correlated') ) then
      TS_W_K_METHOD = TS_W_K_CORRELATED
    else if ( leqi(chars,'uncorrelated') ) then
      TS_W_K_METHOD = TS_W_K_UNCORRELATED
    else
      call die('Could not determine flag TS.Weight.k.Method, &
          &please see manual.')
    end if

    ! The default weighting method is correlated if
    ! atom-atom is utilised
    TS_W_METHOD = TS_W_CORRELATED
    chars = fdf_get('TS.Weight.Method','orb-orb')
    ! first check whether we have correlated weighting
    i = index(chars,'+')
    if ( i > 0 ) then
      ! we do have something else
      if ( leqi(chars(1:i-1),'correlated') .or. &
          leqi(chars(1:i-1),'corr') ) then
        TS_W_METHOD = TS_W_CORRELATED
      else if ( leqi(chars(1:i-1),'uncorrelated') .or. &
          leqi(chars(1:i-1),'uncorr') ) then
        TS_W_METHOD = 0 ! non-correlated
      else
        call die('Unrecognized second option for TS.Weight.Method &
            &must be [[un]correlated+][orb-orb|tr-atom-atom|sum-atom-atom|mean]')
      end if
      chars = chars(i+1:)
    end if
    if ( leqi(chars,'orb-orb') ) then
      TS_W_METHOD = TS_W_ORB_ORB
      ! this does not make sense to make correlated, hence always assign
    else if ( leqi(chars,'tr-atom-atom') ) then
      TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ATOM 
    else if ( leqi(chars,'tr-atom-orb') ) then
      TS_W_METHOD = TS_W_METHOD + TS_W_TR_ATOM_ORB
    else if ( leqi(chars,'sum-atom-atom') ) then
      TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ATOM 
    else if ( leqi(chars,'sum-atom-orb') ) then
      TS_W_METHOD = TS_W_METHOD + TS_W_SUM_ATOM_ORB
    else if ( leqi(chars,'mean') ) then
      TS_W_METHOD = TS_W_MEAN
    else
      call die('Unrecognized option for TS.Weight.Method &
          &must be [[un]correlated+|][orb-orb|tr-atom-[atom|orb]|sum-atom-[atom|orb]|mean]')
    end if

  end subroutine read_ts_weight

  subroutine weight_DM(N_Elec,Elecs,N_mu, mus, na_u, lasto, &
      sp_dist, sparse_pattern, S, &
      spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)

    use units, only: eV
#ifdef MPI
    use mpi_siesta
#endif
    use parallel,  only: IONode, Node, Nodes
    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D

    use m_ts_electype
    use m_ts_chem_pot, only : ts_mu
    use m_ts_contour_neq, only : N_nEq_ID, ID2mu
#ifdef TRANSIESTA_DEBUG_NELEC
    use m_ts_contour_neq, only : nEq_ID
#endif
    use geom_helper, only : iaorb, ucorb

    use intrinsic_missing, only: SFIND

    implicit none

    ! *********************
    ! * OUTPUT variables  *
    ! *********************
    integer,            intent(in) :: N_Elec
    type(Elec),         intent(in) :: Elecs(N_Elec)
    integer,            intent(in) :: N_mu
    type(ts_mu),        intent(in) :: mus(N_mu)
    ! The last-orbital of each atom
    integer, intent(in) :: na_u, lasto(0:na_u)
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    real(dp), intent(in) :: S(:)
    ! Contour part of DM integration
    type(dSpData2D), intent(inout) :: spDM
    ! Real-axis part of DM integration
    type(dSpData2D), intent(inout) :: spDMneq
    ! Estimates of EDM
    type(dSpData2D), intent(inout) :: spEDM
    ! Number of supercells
    integer, intent(in) :: n_s
    ! the offsets
    real(dp), intent(in) :: sc_off(3,0:n_s-1)
    ! The correction to the total energy
    real(dp), intent(inout) :: DE_NEGF

    ! *********************
    ! * LOCAL variables   *
    ! *********************
    real(dp) :: w(N_mu), neq(N_mu)
    ! The total charge per chemical potential.
    real(dp) :: q(N_mu)
    real(dp), parameter :: EPS = 1.e-4_dp

    ! arrays for looping in the sparsity pattern
    integer,  pointer :: s_ncol(:), s_ptr(:), s_col(:)
    integer :: snr, sj, sind
    type(Sparsity), pointer :: sp
    type(OrbitalDistribution), pointer :: dit
    real(dp), pointer :: DM(:,:), DMneq(:,:), EDM(:,:)
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer :: nr, n_nzs
    integer :: io, jo, ind, j, is
    integer :: mu_i, mu_j
    integer, allocatable :: ID_mu(:)
    ! For error estimation
    integer  :: eM_i, eM_j
    real(dp) :: eM, DMe, ee, tmp, m_err, e_f, ew
    integer  :: EeM_i, EeM_j
    real(dp) :: EDMe, Em_err, Eee, EeM, Eew, Ee_f
    logical :: hasEDM, is_correlated, is_trace
    ! collecting the error contribution for each atom
    real(dp), allocatable :: atom_w(:,:), atom_neq(:,:), cor(:)
    integer :: ng, lio, ia1, ia2, ia, TS_W

#ifdef MPI
    integer :: MPIerror
    integer :: status(MPI_STATUS_SIZE)
#endif

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE weightDM' )
#endif

    ! TODO Enforce that sparsity is the same
    ! (however, we know that they are the same)
    sp => spar(spDM)
    call attach(sp,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
        nrows=nr,nrows_g=ng,nnzs=n_nzs)


    ! Obtain the values in the arrays...
    DM     => val(spDM)
    DMneq  => val(spDMneq)
    hasEDM = initialized(spEDM)
    if ( hasEDM ) EDM => val(spEDM)

    ! point to the orbital-distribution
    dit => dist(spDM)

    ! Retrieve information about the overlap matrix sparse pattern.
    call attach(sparse_pattern,n_col=s_ncol,list_ptr=s_ptr,list_col=s_col, &
        nrows=snr)

    if ( nr /= snr ) then
      call die('Error in weight_DM, non-equivalent sparse patterns.')
    end if
    if ( .not. same(sp_dist, dit) ) then
      call die('Currently weight_DM requires the same S and spDM distribution')
    end if

    allocate(cor(N_nEq_ID))
    allocate(ID_mu(N_nEq_ID))
    do io = 1 , N_nEq_ID
      ID_mu(io) = ID2mu(io)
    end do

#ifdef TRANSIESTA_DEBUG_NELEC
    allocate(atom_w(N_mu,1))
    allocate(atom_neq(N_nEq_ID,1))
    do io = 1 , N_nEq_ID
      atom_neq(io,1) = io
    end do

    ! calculate fake contribution
    call calc_neq_weight(N_mu,N_nEq_ID,ID_mu,atom_neq,neq,w)
    call calc_neq(N_mu,N_nEq_ID,ID_mu,atom_neq,atom_w)

    do io = 1 , N_mu
      ia1 = 0
      do jo = 1 , N_mu
        if ( mus(jo)%ID == io ) then
          ia1 = jo
          exit
        end if
      end do
      if ( ia1 == 0 ) call die('Error')
      if ( ia1 /= io ) call die('Error')

      ! Get correct mu
      print '(a)','Electrode Eq: '//trim(mus(ia1)%name)
      print '(a,2(tr1,e10.5))','  nEq, w: ',neq(ia1), w(ia1)
      print '(a,2(tr1,e10.5))','  nEq: ',atom_w(ia1,1)

      ! Get all corrections
      do jo = 1 , N_nEq_ID
        if ( ia1 == ID_mu(jo) ) then
          print '(a,2(tr1,i0))','  nEq-ID:',jo, nEq_ID(jo)%ID
          print '(a,a)','   mu   : ',trim(nEq_ID(jo)%mu%name)
          print '(a,a)','   Gamma: ',trim(nEq_ID(jo)%El%name)
        end if
      end do

    end do

    deallocate(atom_w,atom_neq)
#endif


    ! initialize the errors
    eM  = 0._dp
    ew  = 0._dp
    m_err = 0._dp
    ! energy density matrix
    EeM = 0._dp
    Eew = 0._dp
    Em_err = 0._dp

    ! Is the data correlated
    is_correlated = TS_W_METHOD >= TS_W_CORRELATED
    TS_W = TS_W_METHOD
    if ( is_correlated ) TS_W = TS_W - TS_W_CORRELATED

    if ( TS_W /= TS_W_ORB_ORB .and. TS_W /= TS_W_MEAN ) then
      ! we are doing weighting per trace/sum of each atom

      ! this will not be that large an array... :)
      allocate(atom_neq(N_nEq_ID,na_u))
      atom_neq(:,:) = 0._dp

      ! determine whether it is trace or sum
      select case ( TS_W )
      case ( TS_W_SUM_ATOM_ATOM , TS_W_SUM_ATOM_ORB )
        is_trace = .false.
      case ( TS_W_TR_ATOM_ATOM , TS_W_TR_ATOM_ORB )
        is_trace = .true.
      case default
        call die('Error in weights: determine trace')
      end select

!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ia,j,ind,ia2,is), &
!$OMP&reduction(+:atom_neq)
      do lio = 1 , nr

        ! We are in a buffer region...
        if ( l_ncol(lio) /= 0 ) then
          io = index_local_to_global(dit,lio)

          ia = iaorb(io,lasto) ! atom-index

          lio_connect: do j = 1 , l_ncol(lio)

            ind = l_ptr(lio) + j

            if ( .not. is_Trace ) then ! TS_W == TS_W_SUM_ATOM_?

              ia2 = iaorb(l_col(ind),lasto)
              ! Only allow the same atom to contribute
              if ( ia2 /= ia ) cycle lio_connect

            else ! TS_W == TS_W_TR_ATOM_?

              ! This is a SC sparsity pattern

              ! Only allow the diagonal entry of
              ! the density matrix
              is = (l_col(ind)-1) / ng
              ! Check the unit-cell offset
              if ( sum(abs(sc_off(:,is))) > EPS ) cycle lio_connect

            end if

            if ( is_correlated ) then
              atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:)
            else
              atom_neq(:,ia) = atom_neq(:,ia) + DMneq(ind,:) ** 2
            end if

          end do lio_connect
        end if
      end do
!$OMP end parallel do

#ifdef MPI
      allocate(atom_w(N_nEq_ID,na_u))
      ! We need to reduce the things
      call MPI_AllReduce(atom_neq(1,1),atom_w(1,1),N_nEq_ID*na_u, &
          MPI_Double_Precision, MPI_Sum, MPI_Comm_World, MPIerror)
      if ( is_correlated ) then
        atom_neq(:,:) = atom_w(:,:) ** 2
      else
        atom_neq(:,:) = atom_w(:,:)
      end if
      deallocate(atom_w)
#else
      if ( is_correlated ) then
        atom_neq(:,:) = atom_neq(:,:) ** 2
      end if
#endif

      ! in case of Bulk or DM_update /= update all we
      ! can set the equivalent atom_w to 2 * maximum value
      ! This will force the nearest electrode to contribute the
      ! most. :)
      tmp = maxval(atom_neq)
      allocate(atom_w(N_mu,na_u))
      l_atom: do ia = 1 , na_u
        do io = 1 , N_Elec
          ! If we DO NOT use bulk electrodes we
          ! do have access to the diagonal correction
          ! contribution. Hence we only overwrite the electrode
          ! weight if Elec%Bulk
          if ( (.not. Elecs(io)%Bulk) .and. Elecs(io)%DM_update /= 2 ) cycle
          ! if we are not in the electrode we do not correct weight
          if ( .not. AtomInElec(Elecs(io),ia) ) cycle
          ! in case of sum with the off-diagonal terms we have to sum (otherwise it should be (:,ia) = tmp ; (mu%ID,ia) = 0._dp) 
          atom_w(:,ia) = 0._dp
          atom_w(Elecs(io)%mu%ID,ia) = 1._dp

          !if (node==0) &
          !     write(*,'(a,i2,2(tr1,g10.5))')'W: ', ia,atom_w(:,ia)

          cycle l_atom

        end do

        ! Calculate weights for this atom
        call calc_weight(N_mu,N_nEq_ID,ID_mu, &
            atom_neq(:,ia), atom_w(:,ia) )
        !if (node==0) &
        !     write(*,'(a,i2,4(tr1,g10.5))')'W: ', ia,atom_w(:,ia),atom_neq(:,ia)

      end do l_atom

      ! clean up
      deallocate(atom_neq)

    end if

    if ( TS_W == TS_W_MEAN ) then
      ! The weight will always be divided
      w(:) = 1._dp / real(N_mu,dp)
    end if

    ! Initialize the charges
    q = 0._dp

    ! Now we need to loop the overlap matrix to perform the calculation
    ! of the charges originating from each of the chemical potentials
    ! It should be noted that even though the formalism specifies a
    ! contribution from each of the electrodes it is superfluous due
    ! to two electrodes having the same chemical potential will have:
    !    (N_1 + N_2) * \mu_{1,2}
    ! which is easily seen to be valid for the full write out.

    ! DM is accessed individually, so we will never have a data race
!$OMP parallel do default(shared), &
!$OMP&private(lio,io,ia1,ia2,j,ind,jo,neq,cor), &
!$OMP&private(sj,sind), &
!$OMP&private(ee,e_f,mu_i,mu_j,tmp), &
!$OMP&private(Eee,Ee_f), &
!$OMP&firstprivate(w), reduction(+:m_err,Em_err,q)
    do lio = 1 , nr
      ! We are in a buffer region...
      if ( l_ncol(lio) /= 0 ) then

        ! The global orbital
        io = index_local_to_global(dit,lio)

        ! Update the weight of the row-atom
        ia1 = iaorb(io,lasto)

        ! Here we will loop the overlap matrix
        ! because the spDM sparse patterns are sorted
        ! we can more easily search them.
        do sj = 1, s_ncol(lio)

          sind = s_ptr(lio) + sj

          ! Find the equivalent orbital, or skip
          ind = l_ptr(lio) + &
              SFIND(l_col(l_ptr(lio)+1:l_ptr(lio)+l_ncol(lio)), s_col(sind))
          if ( ind <= l_ptr(lio) ) cycle ! The element does not exist

          ! Retrieve the connecting orbital
          jo = l_col(ind)
          cor(:) = DMneq(ind,:)

          if ( TS_W == TS_W_ORB_ORB ) then

            ! Get the non-equilibrium contribution and the weight associated
            call calc_neq_weight(N_mu,N_nEq_ID,ID_mu, &
                cor,neq,w)

          else if ( TS_W == TS_W_MEAN ) then

            ! "w" already set
            ! Get the non-equilibrium contribution
            call calc_neq(N_mu,N_nEq_ID,ID_mu,cor,neq)

          else ! we have weight per atom "somewhere"

            ! To compare the weights... For DEBUGging purposes...
            !call get_neq_weight(N_mu,N_nEq_ID,ID_mu, &
            !     cor,neq,w)
            !write(*,'(a,i2,tr1,i2,4(tr1,g10.5))')'Wi: ', ia1,ia2,w

            ! Re-calculate the weight for special weighting...
            ia2 = iaorb(jo,lasto)

            ! [[ this is the geometric mean method...
            ! This should probably be used as the geometric mean
            ! retains the tendency of the data, however I am not 
            ! quite sure of its arguments...
            ! For test:
            ! A1:    [ 0.95  0.05 ]   # weight 1
            ! A2:    [ 0.5   0.5  ]   # weight 2
            ! GM-N:  [ 0.813 0.187]   # geometric mean of weights
            ! AM-N:  [ 0.725 0.275]   # arithmetic mean of weights

            w = sqrt(atom_w(:,ia1) * atom_w(:,ia2))

            ! geometric mean does not retain normalization
            w = w / sum(w) 
            ! ]]

            ! [[ arithmetic mean
            ! the mean value between the atomic weights
            ! will be used as the actual weight
            ! w = (atom_w(:,ia1) + atom_w(:,ia2)) * 0.5_dp
            ! ]] 

            select case ( TS_W )                   
            case ( TS_W_TR_ATOM_ATOM , TS_W_SUM_ATOM_ATOM )

              ! do nothing...

            case ( TS_W_TR_ATOM_ORB , TS_W_SUM_ATOM_ORB )

              ! this ensures that the atomic weights are
              ! both taken into account, as well as the orb-orb
              call calc_weight(N_mu,N_nEq_ID,ID_mu, &
                  cor ** 2, neq )

              ! see above for arguments
              w = sqrt(w(:) * neq(:))
              w = w / sum(w) 

            case default
              call die('Error in weights...')
            end select

            ! Get the non-equilibrium contribution
            call calc_neq(N_mu,N_nEq_ID,ID_mu,cor,neq)

            !write(*,'(a,i2,tr1,i2,4(tr1,g10.5))')'Wt: ', ia1,ia2,w,sum(w)
          end if

#ifdef TRANSIESTA_WEIGHT_DEBUG
          if ( io == ucorb(jo,ng) .and. io == 28 ) then
            print '(2(a7,3(tr1,f10.5)))','Left',DM(ind,1),neq(1),w(1), &
                'Right',DM(ind,2),neq(2),w(2)
          end if
#endif

          ! Calculate each contribution
          do mu_i = 1 , N_mu
            DM(ind,mu_i) = DM(ind,mu_i) + neq(mu_i)
          end do

          ! Do error estimation (capture before update)
          ee = 0._dp
          Eee = 0._dp
          if ( hasEDM ) then

            do mu_i = 1 , N_mu - 1
              do mu_j = mu_i + 1 , N_mu
                ! Density matrix
                tmp = DM(ind,mu_i) - DM(ind,mu_j)
                ! Calculate sum of all errors
                m_err = m_err + tmp
                if ( abs(tmp) > abs(ee) ) ee = tmp

                ! Energy density matrix
                tmp = EDM(ind,mu_i) - EDM(ind,mu_j)
                Em_err = Em_err + tmp
                if ( abs(tmp) > abs(Eee) ) Eee = tmp

              end do
            end do

            ! Store for later estimation of the "final" error
            e_f  = DM(ind,1)
            Ee_f = EDM(ind,1)

            ! Also calculate the charge after weighting
            q(1) = q(1) + w(1) * DM(ind,1) * S(sind)
            DM(ind,1) = w(1) * DM(ind,1)
            EDM(ind,1) = w(1) * EDM(ind,1)
            do mu_i = 2 , N_mu
              q(mu_i) = q(mu_i) + w(mu_i) * DM(ind,mu_i) * S(sind)
              DM(ind,1) = DM(ind,1) + w(mu_i) * DM(ind,mu_i)
              EDM(ind,1) = EDM(ind,1) + w(mu_i) * EDM(ind,mu_i)
            end do

            ! Calculate error from estimated density
            e_f = e_f - DM(ind,1)
            Ee_f = Ee_f - EDM(ind,1)
            do mu_i = 2 , N_mu
              tmp = DM(ind,mu_i) - DM(ind,1)
              if ( abs(tmp) > abs(e_f) ) e_f = tmp
              tmp = EDM(ind,mu_i) - EDM(ind,1)
              if ( abs(tmp) > abs(Ee_f) ) Ee_f = tmp
            end do

            if ( abs(Eee) > abs(EeM) ) then
!$OMP critical
              if ( abs(Eee) > abs(EeM) ) then
                ! Energy density matrix
                EeM_i = io
                EeM_j = jo
                EeM   = Eee
                EDMe  = EDM(ind,1)
                Eew   = Ee_f
              end if
!$OMP end critical
            end if

          else ! no EDM calculation

            do mu_i = 1 , N_mu - 1
              do mu_j = mu_i + 1 , N_mu
                ! Density matrix
                tmp = DM(ind,mu_i) - DM(ind,mu_j)
                ! Calculate sum of all errors
                m_err = m_err + tmp
                if ( abs(tmp) > abs(ee) ) ee = tmp

              end do
            end do

            ! Store for later estimation of the "final" error
            e_f = DM(ind,1)

            ! Also calculate the charge after weighting
            q(1) = q(1) + w(1) * DM(ind,1) * S(sind)
            DM(ind,1) = w(1) * DM(ind,1)
            do mu_i = 2 , N_mu
              q(mu_i) = q(mu_i) + w(mu_i) * DM(ind,mu_i) * S(sind)
              DM(ind,1) = DM(ind,1) + w(mu_i) * DM(ind,mu_i)
            end do

            ! Calculate error from estimated density
            e_f = e_f - DM(ind,1)
            do mu_i = 2 , N_mu
              tmp = DM(ind,mu_i) - DM(ind,1)
              if ( abs(tmp) > abs(e_f) ) e_f = tmp
            end do

          end if

          if ( abs(ee) > abs(eM) ) then
!$OMP critical
            ! This double if takes care of OpenMP 
            ! critical region
            if ( abs(ee) > abs(eM) ) then
              ! Density matrix
              eM_i = io
              eM_j = jo
              eM   = ee
              DMe  = DM(ind,1)
              ew   = e_f
            end if
!$OMP end critical
          end if

        end do
      end if
    end do
!$OMP end parallel do

    ! Calculate mean of mean difference
    ! First calculate number of differences used
    io = 0
    do mu_i = 1, N_mu - 1
      do mu_j = mu_i + 1 , N_mu
        io = io + 1
      end do
    end do
    m_err  = m_err / real(io,dp)
    Em_err = Em_err / real(io,dp)

    if ( TS_W /= TS_W_ORB_ORB .and. TS_W /= TS_W_MEAN ) then
      deallocate(atom_w)
    end if

    deallocate(ID_mu)
    deallocate(cor)

#ifdef MPI
    if ( Nodes > 1 ) then
      ! nullify pointer
      nullify(DM)

      ! Gather all maximum errors on the IO-node
      allocate(DM(2,Nodes))
      w(1) = eM
      w(2) = EeM

      ! First figure out which process it is
      call MPI_Gather(w(1),2,MPI_Double_Precision, &
          DM(1,1),2,MPI_Double_Precision, &
          0, MPI_Comm_World, MPIerror)
      if ( Node == 0 ) then
        ! Figure out which node contains the largest error
        ia1 = 0
        neq(1) = DM(1,1) ! density matrix
        ia2 = 0
        neq(2) = DM(2,1) ! energy density matrix
        do io = 2 , Nodes
          if ( abs(DM(1,io)) > abs(neq(1)) ) then
            neq(1) = DM(1,io)
            ia1 = io - 1
          end if
          if ( abs(DM(2,io)) > abs(neq(2)) ) then
            neq(2) = DM(2,io)
            ia2 = io - 1
          end if
        end do
      end if

      ! Clean-up
      deallocate(DM)

      ! Reduce number of updated elements
      allocate(DM(3,2))
      DM(1,1) = real(n_nzs,dp)
      DM(2,1) = m_err   ! DM mean-error
      DM(3,1) = Em_err  ! EDM mean-error
      call MPI_Reduce(DM(1,1),DM(1,2),3,MPI_Double_Precision, &
          MPI_SUM, 0, MPI_Comm_World, MPIerror)
      ! Get total number of updated elements
      n_nzs  = int(DM(1,2))
      ! total error
      m_err  = DM(2,2)
      Em_err = DM(3,2)

      ! Clean-up
      deallocate(DM)

      ! B-cast nodes that has the highest error
      allocate(ID_mu(2))
      ID_mu(1) = ia1
      ID_mu(2) = ia2
      call MPI_Bcast(ID_mu,2,MPI_Integer,0, MPI_Comm_World, MPIerror)
      ia1 = ID_mu(1)
      ia2 = ID_mu(2)
      deallocate(ID_mu)

      ! Allocate for sending of data
      allocate(DM(5,1))

      ! Density matrix
      DM(1,1) = real(eM_i,dp)
      DM(2,1) = real(eM_j,dp)
      DM(3,1) = eM ! maximum difference between mu_i
      DM(4,1) = DMe ! final density at max-diff
      DM(5,1) = ew ! maximum difference from final density
      if ( Node == 0 .and. Node == ia1 ) then
        ! do nothing, everything is updated
      else if ( Node == 0 ) then
        ! recv from main node
        call MPI_Recv(DM(1,1),5, MPI_Double_Precision, &
            ia1, 0, MPI_Comm_World, status, MPIerror)
        eM_i = int(DM(1,1))
        eM_j = int(DM(2,1))
        eM   = DM(3,1)
        DMe  = DM(4,1)
        ew   = DM(5,1)
      else if ( Node == ia1 ) then
        ! send to main node
        call MPI_Send(DM(1,1),5, MPI_Double_Precision, &
            0, 0, MPI_Comm_World, MPIerror)
      end if

      ! Energy density matrix
      DM(1,1) = real(EeM_i,dp)
      DM(2,1) = real(EeM_j,dp)
      DM(3,1) = EeM
      DM(4,1) = EDMe
      DM(5,1) = Eew
      if ( Node == 0 .and. Node == ia2 ) then
        ! do nothing, everything is updated
      else if ( Node == 0 ) then
        ! recv from main node
        call MPI_Recv(DM(1,1),5, MPI_Double_Precision, &
            ia2, 0, MPI_Comm_World, status, MPIerror)
        EeM_i = int(DM(1,1))
        EeM_j = int(DM(2,1))
        EeM   = DM(3,1)
        EDMe  = DM(4,1)
        Eew   = DM(5,1)
      else if ( Node == ia2 ) then
        ! send to main node
        call MPI_Send(DM(1,1),5, MPI_Double_Precision, &
            0, 0, MPI_Comm_World, MPIerror)
      end if

      deallocate(DM)

    end if

    ! AllReduce the charges for calculating the DE_NEGF
    call MPI_AllReduce(q,w,N_mu,MPI_Double_Precision, &
        MPI_SUM, MPI_Comm_World, MPIerror)
    q = w
#endif

    ! Calculate mean
    m_err  =  m_err / real(n_nzs,dp)
    Em_err = Em_err / real(n_nzs,dp)

    call print_error_estimate(IONode,'ts-err-D:', &
        eM,ew,eM_i,eM_j,DMe,m_err)
    if ( hasEDM ) then
      EeM    = EeM / eV
      Eew    = Eew / eV
      EDMe   = EDMe / eV
      Em_err = Em_err / eV
      call print_error_estimate(IONode,'ts-err-E:', &
          EeM,Eew,EeM_i,EeM_j,EDMe,Em_err)
    end if

    ! Calculate the energy contribution to the total energy due
    ! to the electrons from the baths
    !   Etot = Etot - e \sum_i N_i \mu_i
    call print_charge_weight(IONode, 'ts-w-q:', q)
    do mu_i = 1, N_mu
      DE_NEGF = DE_NEGF + q(mu_i) * mus(mu_i)%mu
    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS weightDM' )
#endif

  end subroutine weight_DM

  ! ************* Commonly used modules ****************

  ! Write out the charge contributions from each of the chemical potentials
  subroutine print_charge_weight(IONode,a,q)
    logical, intent(in)  :: IONode
    character(len=*), intent(in) :: a
    real(dp), intent(in) :: q(:)
    integer :: i

    if ( IONode ) then
      write(*,'(a,tr8)',advance='no') trim(a)
      do i = 1 , size(q)
        write(*,'(tr1,a8,i0)',advance='no') 'qP',i
      end do
      ! Ensure there is a new-line
      write(*,'(/a,tr8,1000(tr1,f9.3))') trim(a), q
    end if

  end subroutine print_charge_weight

  ! Write out the error-estimate for the current iteration (or, k-point)
  subroutine print_error_estimate(IONode,a,eM,ew,eM_i,eM_j,DM,m_err)
    logical, intent(in)  :: IONode
    character(len=*), intent(in) :: a
    real(dp), intent(in) :: eM, ew, DM, m_err
    integer, intent(in)  :: eM_i,eM_j

    if ( IONode ) then
      write(*,'(a,tr1,a,i5,'','',i6,3(a,g10.5e1)&
          &,a,g11.5)') trim(a), &
          'ij(',eM_i,eM_j, '), M = ',DM, &
          ', ew = ',ew, ', em = ',eM, &
          '. avg_em = ',m_err
    end if

  end subroutine print_error_estimate


  ! Calculate the theta values
  pure subroutine calc_weight(N_mu,N_id,ID_mu,w_ID,w)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: w_ID(N_id)
    real(dp), intent(out) :: w(N_mu)
    real(dp) :: theta(N_mu), tmp
    integer :: i, j

    call calc_neq(N_mu,N_id,ID_mu,w_ID,theta)
    ! We need to do it this complicated since
    ! one of theta's may be 0.
    tmp = 0._dp
    do i = 1 , N_mu
      w(i) = 1._dp
      do j = 1, N_mu
        if ( i /= j ) then
          w(i) = w(i) * theta(j)
        end if
      end do
      tmp = tmp + w(i)
    end do

    ! The denominator
    if ( tmp == 0._dp ) then
      w = 1._dp / real(N_mu,dp)
    else
      w = w / tmp
    end if

  end subroutine calc_weight

  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  pure subroutine calc_neq_weight(N_mu,N_id,ID_mu,neq_ID,neq,w)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    real(dp), intent(out) :: w(N_mu)
    integer :: i, j
    real(dp) :: tmp

    ! TODO check that this is correct for several electrodes
    call calc_neq(N_mu,N_id,ID_mu,neq_ID**2,neq)
    tmp = 0._dp
    do i = 1 , N_mu
      ! Calculate weight for chemical potential i
      w(i) = 1._dp
      do j = 1, N_mu
        if ( i /= j ) then
          w(i) = w(i) * neq(j)
        end if
      end do
      tmp = tmp + w(i)
    end do

    ! The denominator
    if ( tmp == 0._dp ) then
      w = 1._dp / real(N_mu,dp)
    else
      w = w / tmp
    end if

    call calc_neq(N_mu,N_id,ID_mu,neq_ID,neq)

  end subroutine calc_neq_weight

  ! Calculate both the non-equilibrium contribution
  ! and the weight associated with those points
  ! Note that the non-equilibrium contribution calculation
  ! is the same as the theta calculation.
  pure subroutine calc_neq(N_mu,N_id,ID_mu,neq_ID,neq)
    integer,  intent(in)  :: N_mu, N_id, ID_mu(N_id)
    real(dp), intent(in)  :: neq_ID(N_id)
    real(dp), intent(out) :: neq(N_mu)
    integer :: ID

    ! TODO check that this is correct for several electrodes
    neq(:) = 0._dp
    do ID = 1 , N_id
      neq(ID_mu(ID)) = neq(ID_mu(ID)) + neq_ID(ID)
    end do

  end subroutine calc_neq

end module m_ts_weight
