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
! * It has been heavily inspired by the original authors of the 
!   Transiesta code (hence the references here are still remaining) *

module m_ts_mumpsg

  use precision, only : dp

#ifdef SIESTA__MUMPS  
  use m_ts_sparse_helper

  use m_ts_dm_update, only : init_DM
  use m_ts_dm_update, only : update_DM
  use m_ts_dm_update, only : add_Gamma_DM
  
  use m_ts_weight, only : weight_DM

  use m_ts_method, only: orb_offset, no_Buf

  use m_ts_mumps_init
  
  implicit none
  
  public :: ts_mumpsg
  
  private
  
contains
  
  subroutine ts_mumpsg(N_Elec,Elecs, &
       nq, uGF, nspin, na_u, lasto, &
       sp_dist, sparse_pattern, &
       no_u, n_nzs, &
       Hs, Ss, DM, EDM, Ef, DE_NEGF)

    use units, only : eV, Pi
    use parallel, only : Node, Nodes
#ifdef MPI
    use mpi_siesta
#endif

    use alloc, only : re_alloc, de_alloc

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D

    use m_ts_electype
    ! Self-energy read
    use m_ts_gf
    ! Self-energy expansion
    use m_ts_elec_se

    use m_ts_options, only : Calc_Forces
    use m_ts_options, only : N_mu, mus

    use m_ts_options, only : IsVolt

    use m_ts_sparse, only : ts_sp_uc, tsup_sp_uc
    use m_ts_sparse, only : ltsup_sp_sc, ltsup_sc_pnt
    use m_ts_sparse, only : sc_off

    use m_ts_cctype
    use m_ts_contour_eq,  only : Eq_E, ID2idx, c2weight_eq
    use m_ts_contour_neq, only : nEq_E, has_cE_nEq
    use m_ts_contour_neq, only : N_nEq_ID, c2weight_neq
    
    use m_iterator

    use m_ts_mumps_scat

    use ts_dq_m, only: ts_dq

#ifdef TRANSIESTA_DEBUG
    use m_ts_debug
#endif

    include 'zmumps_struc.h'

! ********************
! * INPUT variables  *
! ********************
    integer, intent(in) :: N_Elec
    type(Elec), intent(inout) :: Elecs(N_Elec)
    integer, intent(in) :: nq(N_Elec), uGF(N_Elec)
    integer, intent(in) :: nspin, na_u, lasto(0:na_u)
    type(OrbitalDistribution), intent(inout) :: sp_dist
    type(Sparsity), intent(inout) :: sparse_pattern
    integer, intent(in)  :: no_u
    integer, intent(in)  :: n_nzs
    real(dp), intent(in) :: Hs(n_nzs,nspin), Ss(n_nzs)
    real(dp), intent(inout) :: DM(n_nzs,nspin), EDM(n_nzs,nspin)
    real(dp), intent(in) :: Ef
    real(dp), intent(inout) :: DE_NEGF

! ******************* Computational arrays *******************
    integer :: ndwork, nzwork, n_s
    real(dp), pointer :: dwork(:,:)
    ! The solution arrays
    type(zMUMPS_STRUC) :: mum
    complex(dp), pointer :: zwork(:), Gf(:)

    ! A local orbital distribution class (this is "fake")
    type(OrbitalDistribution) :: fdist
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D) :: spH, spS
    ! local sparsity pattern in local SC pattern
    type(dSpData2D) :: spDM, spDMneq
    type(dSpData2D) :: spEDM ! only used if calc_forces
    ! The different sparse matrices that will surmount to the integral
    ! These two lines are in global update sparsity pattern (UC)
    type(dSpData2D) ::  spuDM
    type(dSpData2D) :: spuEDM ! only used if calc_forces
! ************************************************************

! ******************* Computational variables ****************
    type(ts_c_idx) :: cE
    integer :: index_dq !< Index for the current charge calculation @ E == mu
    real(dp) :: kw, dq_mu
    complex(dp) :: W, ZW
! ************************************************************

! ******************** Loop variables ************************
    type(itt1) :: Sp
    integer, pointer :: ispin
    integer :: iEl, iID, up_nzs
    integer :: iE, imu, idx
! ************************************************************

! ******************* Miscalleneous variables ****************
    integer :: no_u_TS, off, no
    real(dp), parameter :: bkpt(3) = (/0._dp,0._dp,0._dp/)
! ************************************************************

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'PRE transiesta mem' )
#endif

    ! Number of supercells (even though its gamma we
    ! can have different schemes...)
    n_s = size(sc_off,dim=2)

    ! Number of orbitals in TranSIESTA
    no_u_TS = no_u - no_Buf

    ! Number of elements that are transiesta updated
    up_nzs = nnzs(tsup_sp_uc)

    nullify(Gf)

    ! We do need the full GF AND a single work array to handle the
    ! left-hand side of the inversion...
    ! We will provide all work arrays as single dimension arrays.
    ! This will make interfaces more stringent and allow for
    ! re-use in several other places.
    ! However, this comes at the cost of programmer book-keeping.
    ! It should be expected that the work arrays return GARBAGE
    ! on ALL routines, i.e. they are not used for anything other
    ! than, yes, work.

#ifdef TRANSIESTA_TIMING
    call timer('TS_MUMPS_INIT',1)
#endif

    ! initialize MUMPS
    call init_MUMPS(mum,Node)
    
    ! Initialize for the correct size
    mum%N    = no_u_TS ! order of matrix
    mum%NRHS = no_u_TS ! RHS-vectors

    ! prepare the LHS of MUMPS-solver.
    ! This is the same no matter the contour type
    call prep_LHS(IsVolt,mum,N_Elec,Elecs)
    zwork => mum%A(:)
    nzwork = size(zwork)

    ! analyzation step
    call analyze_MUMPS(mum)

    if ( .not. IsVolt ) then
       ! We can allocate the equilibrium now
       call prep_RHS_Eq(mum,no_u_TS,up_nzs,N_Elec,Elecs,Gf)
    end if

#ifdef TRANSIESTA_TIMING
    call timer('TS_MUMPS_INIT',2)
#endif

    ! Create the Fake distribution
    ! The Block-size is the number of orbitals, i.e. all on the first processor
    ! Notice that we DO need it to be the SIESTA size.
#ifdef MPI
    call newDistribution(no_u,MPI_COMM_WORLD,fdist,name='TS-fake dist')
#else
    call newDistribution(no_u,-1            ,fdist,name='TS-fake dist')
#endif

    ! The Hamiltonian and overlap matrices (in Gamma calculations
    ! we will not have any phases, hence, it makes no sense to
    ! have the arrays in complex)
    call newdSpData1D(ts_sp_uc,fdist,spH,name='TS spH')
    call newdSpData1D(ts_sp_uc,fdist,spS,name='TS spS')

    ! If we have a bias calculation we need additional arrays.
    ! If not bias we don't need the update arrays (we already have
    ! all information in tsup_sp_uc (spDMu))

    ! Allocate space for global sparsity arrays
    no = max(N_mu,N_nEq_id)
    call newdSpData2D(tsup_sp_uc,no,fdist, spuDM, name='TS spuDM')
    ! assign dwork, this will problably come at the expence of
    ! two full reductions, however, we save some memory...
    ndwork = nnzs(tsup_sp_uc)
    dwork => val(spuDM)
    if ( Calc_Forces ) then
       call newdSpData2D(tsup_sp_uc,N_mu,fdist, spuEDM, name='TS spuEDM')
    end if
    
    if ( IsVolt ) then
       ! Allocate space for update arrays, local sparsity arrays
       call newdSpData2D(ltsup_sp_sc,N_mu,    sp_dist,spDM   ,name='TS spDM')
       call newdSpData2D(ltsup_sp_sc,N_nEq_id,sp_dist,spDMneq,name='TS spDM-neq')
       if ( nnzs(ltsup_sp_sc) > ndwork ) then
          ! only update if this array is larger (should only happen in 
          ! few processor setups
          ndwork = nnzs(ltsup_sp_sc)
          dwork => val(spDMneq)
       end if
       if ( Calc_Forces ) then
          call newdSpData2D(ltsup_sp_sc,N_mu, sp_dist,spEDM  ,name='TS spEDM')
       end if
    end if

    ! Initialize the charge correction scheme (will return if not used)
    call ts_dq%initialize_dq()

    ! start the itterators
    call itt_init  (Sp,end=nspin)
    ! point to the index iterators
    call itt_attach(Sp,cur=ispin)

    do while ( .not. itt_step(Sp) )

       write(mum%ICNTL(1),'(/,a,i0,/)') '### Solving for spin: ',ispin

       call init_DM(sp_dist, sparse_pattern, &
            n_nzs, DM(:,ispin), EDM(:,ispin), &
            tsup_sp_uc, Calc_Forces)

       ! Include spin factor and 1/\pi
       kw = 2._dp / (Pi * nspin)

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',1)
#endif

       ! Work-arrays are for MPI distribution...
       call create_HS(sp_dist,sparse_pattern, &
            Ef, &
            N_Elec, Elecs, no_u, & ! electrodes, SIESTA size
            n_nzs, Hs(:,ispin), Ss, &
            spH, spS, &
            ndwork, dwork(:,1)) ! annoyingly we can't pass the full array!!!!!

#ifdef TRANSIESTA_TIMING
       call timer('TS_HS',2)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',1)
#endif

       if ( IsVolt ) &
            call prep_RHS_Eq(mum,no_u_TS,up_nzs,N_Elec,Elecs,Gf)

       ! ***************
       ! * EQUILIBRIUM *
       ! ***************
       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       no = no_u_TS
       do iEl = 1 , N_Elec
          if ( Elecs(iEl)%DM_update == 0 ) then
             no = no - TotUsedOrbs(Elecs(iEl))
          end if
       end do
       iE = Nodes - Node
       cE = Eq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(ispin, 1, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .false. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
#ifdef TRANSIESTA_TIMING
          call timer('TS_PREP',1)
#endif
          call prepare_invGF(cE, mum, &
               N_Elec, Elecs, &
               spH=spH , spS=spS)
#ifdef TRANSIESTA_TIMING
          call timer('TS_PREP',2)
#endif

          ! *******************
          ! * calc GF         *
          ! *******************
          if ( .not. cE%fake ) then
#ifdef TRANSIESTA_TIMING
             call timer('TS_MUMPS_SOLVE',1)
#endif
             write(mum%ICNTL(1),'(/,a,i0,2(a,i0),/)') &
                  '### Solving Eq Node/iC: ',Node,'/',cE%idx(2),',',cE%idx(3)
             mum%JOB = 5
             call zMUMPS(mum)
             call mum_err(mum, &
                  'MUMPS failed the Eq. inversion, check the output log')
#ifdef TRANSIESTA_TIMING
             call timer('TS_MUMPS_SOLVE',2)
#endif   
             
             ! ** At this point we have calculated the Green function
             
             ! ****************
             ! * save GF      *
             ! ****************
             do imu = 1 , N_mu
               call ID2idx(cE,mus(imu)%ID,idx)
               if ( idx < 1 ) cycle

               call c2weight_eq(cE,idx, kw, W ,ZW)

               ! Figure out if this point is a charge-correction energy
               index_dq = ts_dq%get_index(imu, iE)
               if ( index_dq > 0 ) then
                 ! Accummulate charge at the electrodes chemical potentials
                 ! Note that this dq_mu does NOT have the prefactor 1/Pi
                 call add_DM( spuDM, W, spuEDM, ZW, &
                     mum, &
                     N_Elec, Elecs, &
                     DMidx=mus(imu)%ID, &
                     spS=spS, q=dq_mu)
                 ts_dq%mus(imu)%dq(index_dq) = ts_dq%mus(imu)%dq(index_dq) + dq_mu * kw
               else
                 call add_DM( spuDM, W, spuEDM, ZW, &
                     mum, &
                     N_Elec, Elecs, &
                     DMidx=mus(imu)%ID)
               end if
             end do
          end if
             
          ! step energy-point
          iE = iE + Nodes
          cE = Eq_E(iE,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_EQ',2)
#endif

#ifdef MPI
       ! We need to reduce all the arrays
       call MPI_Barrier(MPI_Comm_World,iE)
       call timer('TS_comm',1)
       call my_full_G_reduce(spuDM,nzwork*2,zwork,N_mu)
       if ( Calc_Forces ) then
          call my_full_G_reduce(spuEDM,nzwork*2,zwork,N_mu)
       end if
       call timer('TS_comm',2)
#endif

       if ( .not. IsVolt ) then
          call update_DM(sp_dist,sparse_pattern, n_nzs, &
               DM(:,ispin), spuDM, Ef=Ef, &
               EDM=EDM(:,ispin), spEDM=spuEDM, &
               UpSpGlobal = .true.)

          ! The remaining code segment only deals with 
          ! bias integration... So we skip instantly

          cycle

       end if

       ! *****************
       ! * only things with non-Equilibrium contour...
       ! *****************

       ! initialize to zero
       ! local sparsity update patterns
       call init_val(spDM)
       call init_val(spDMneq)
       if ( Calc_Forces ) call init_val(spEDM)

       ! transfer data to local sparsity arrays
       call add_Gamma_DM(spDM,   spuDM, D_dim2=N_mu, &
            spEDM=spEDM, spuEDM=spuEDM, E_dim2=N_mu)

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',1)
#endif

       call prep_RHS_nEq(mum,no_u_TS,N_Elec,Elecs,Gf)

       ! *******************
       ! * NON-EQUILIBRIUM *
       ! *******************

       ! We have the definition of: Gamma = i(\Sigma - \Sigma^\dagger)
       ! (not with one half)
       ! Hence we need to half the contribution for the non-equilibrium
       kw = 0.5_dp * kw

       call init_val(spuDM)
       if ( Calc_Forces ) call init_val(spuEDM)
       iE = Nodes - Node
       cE = nEq_E(iE,step=Nodes) ! we read them backwards
       do while ( cE%exist )

          ! *******************
          ! * prep Sigma      *
          ! *******************
          call read_next_GS(ispin, 1, bkpt, &
               cE, N_Elec, uGF, Elecs, &
               nzwork, zwork, .false., forward = .false. )
          do iEl = 1 , N_Elec
             call UC_expansion(cE, Elecs(iEl), nzwork, zwork, &
                  non_Eq = .true. )
          end do

          ! *******************
          ! * prep GF^-1      *
          ! *******************
          call prepare_invGF(cE, mum, &
               N_Elec, Elecs, &
               spH =spH , spS =spS)
          
          ! *******************
          ! * calc GF         *
          ! *******************
          if ( .not. cE%fake ) then
             
             write(mum%ICNTL(1),'(/,a,i0,2(a,i0),/)') &
                  '### Solving nEq Node/iC: ',Node,'/',cE%idx(2),',',cE%idx(3)
             mum%JOB = 5
             call zMUMPS(mum)
             call mum_err(mum, &
                  'MUMPS failed the nEq. inversion, check the output log')
             
             ! ** At this point we have calculated the Green function
             
             ! ****************
             ! * save GF      *
             ! ****************
             off = 0
             do iEl = 1 , N_Elec
                
                ! The mumps solver is initialized to always
                ! solve for all electrode columns... (not very sparse :( )
                
                ! offset and number of orbitals
                no = TotUsedOrbs(Elecs(iEl))
                
                ! step to the next electrode position
                off = off + no
                
                ! *notice* we correct the Gf index for the column
                call GF_Gamma_GF(Elecs(iEl), mum, no_u_TS, no, &
                     Gf(no_u_TS*(off-no)+1:))
                
                do iID = 1 , N_nEq_ID
                   
                   if ( .not. has_cE_nEq(cE,iEl,iID) ) cycle
                   
                   call c2weight_neq(cE,iID,kw,W,imu,ZW)
                   
                   call add_DM( spuDM, W, spuEDM, ZW, &
                        mum, &
                        N_Elec, Elecs, &
                        DMidx=iID, EDMidx=imu, is_eq = .false. )
                end do
             end do
          end if
          
          ! step energy-point
          iE = iE + Nodes
          cE = nEq_E(iE,step=Nodes) ! we read them backwards
       end do

#ifdef TRANSIESTA_TIMING
       call timer('TS_NEQ',2)
#endif

#ifdef MPI
       ! We need to reduce all the arrays
       call MPI_Barrier(MPI_Comm_World,iE)
       call timer('TS_comm',1)
       call my_full_G_reduce(spuDM, nzwork*2, zwork, N_nEq_id)
       if ( Calc_Forces ) then
          call my_full_G_reduce(spuEDM, nzwork*2, zwork, N_mu)
       end if
       call timer('TS_comm',2)
#endif

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',1)
#endif

       ! 1. move from global UC to local SC
       ! 2. calculate the correct contribution by applying the weight
       ! 3. add the density to the real arrays
       call add_Gamma_DM(spDMneq, spuDM, D_dim2=N_nEq_id, &
            spEDM=spEDM,  spuEDM=spuEDM, E_dim2=N_mu)
       
       call weight_DM( N_Elec, Elecs, N_mu, mus, na_u, lasto, &
            sp_dist, sparse_pattern, Ss, &
            spDM, spDMneq, spEDM, n_s, sc_off, DE_NEGF)
       
       call update_DM(sp_dist,sparse_pattern, n_nzs, &
            DM(:,ispin), spDM, Ef=Ef, &
            EDM=EDM(:,ispin), spEDM=spEDM, ipnt=ltsup_sc_pnt)

#ifdef TRANSIESTA_TIMING
       call timer('TS_weight',2)
#endif
       
       ! We don't need to do anything here..

    end do ! spin

    call itt_destroy(Sp)

#ifdef TRANSIESTA_DEBUG
    write(*,*) 'Completed TRANSIESTA SCF'
#endif

!***********************
! CLEAN UP
!***********************

    call delete(spH)
    call delete(spS)

    call delete(spuDM)
    call delete(spuEDM)

    call delete(spDM)
    call delete(spDMneq)
    call delete(spEDM)

    ! We can safely delete the orbital distribution, it is local
    call delete(fdist)

    ! Deallocate user data
    call memory('D','I',size(mum%IRN)*2, 'prep_LHS')
    deallocate( mum%IRN, mum%JCN )
    nullify( mum%IRN, mum%JCN )
    call memory('D', 'Z', size(mum%A), 'prep_LHS')
    deallocate( mum%A ) ; nullify( mum%A )
    call memory('D', 'I', size(mum%IRHS_PTR), 'allocate_num')
    deallocate( mum%IRHS_PTR ) ; nullify( mum%IRHS_PTR )
    call memory('D', 'I', size(mum%IRHS_SPARSE), 'allocate_num')
    deallocate( mum%IRHS_SPARSE ) ; nullify( mum%IRHS_SPARSE )
    call memory('D', 'Z', size(Gf), 'allocate_num')
    deallocate( Gf ) ! mum%RHS_SPARSE is => GF
    nullify( mum%RHS_SPARSE )
    ! retain IO, killing makes a last print-out
    no = mum%ICNTL(1)

    ! Destroy the instance (deallocate internal data structures)
    mum%JOB = -2
    call zMUMPS(mum)
    call mum_err(mum, 'Cleaning of MUMPS failed.')

    ! close the file
    call io_close(no)

    ! Nullify external pointers
    do iEl = 1, N_Elec
      nullify(Elecs(iEl)%Sigma)
    end do

#ifdef TRANSIESTA_DEBUG
    call write_debug( 'POS transiesta mem' )
#endif

  end subroutine ts_mumpsg


  ! Update DM
  ! These routines are supplied for easy update of the update region
  ! sparsity patterns
  ! Note that these routines implement the usual rho(Z) \propto - GF
  subroutine add_DM(DM, DMfact,EDM, EDMfact, &
       mum, &
       N_Elec,Elecs, &
       DMidx, EDMidx, &
       spS, q, &
       is_eq)

    use intrinsic_missing, only: SFIND

    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    use m_ts_electype
    use m_ts_method, only : ts2s_orb

    include 'zmumps_struc.h'

    ! The DM and EDM equivalent matrices
    type(dSpData2D), intent(inout) :: DM
    complex(dp), intent(in) :: DMfact
    type(dSpData2D), intent(inout) :: EDM
    complex(dp), intent(in) :: EDMfact
    ! the mumps structure
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! the index of the partition
    integer, intent(in) :: DMidx
    integer, intent(in), optional :: EDMidx
    !< Overlap matrix setup for a k-point is needed for calculating q
    type(dSpData1D), intent(in), optional :: spS
    !< Charge calculated at this energy-point
    !!
    !! This does not contain the additional factor 1/Pi
    real(dp), intent(inout), optional :: q
    logical, intent(in), optional :: is_eq

    ! Arrays needed for looping the sparsity
    type(Sparsity), pointer :: s
    integer,  pointer :: l_ncol(:), l_ptr(:), l_col(:)
    integer, pointer :: s_ncol(:), s_ptr(:), s_col(:)
    real(dp), pointer :: D(:,:), E(:,:), Sg(:)
    complex(dp), pointer :: GF(:)
    integer :: io, jo, ind, ir, nr, Hn, ind_H
    integer :: s_ptr_begin, s_ptr_end, sin
    integer :: i1, i2
    logical :: hasEDM, lis_eq, calc_q

    lis_eq = .true.
    if ( present(is_eq) ) lis_eq = is_eq

    calc_q = present(q) .and. present(spS)
    
    ! Remember that this sparsity pattern HAS to be in Global UC
    s => spar(DM)
    call attach(s,n_col=l_ncol,list_ptr=l_ptr,list_col=l_col, &
         nrows=nr)
    D => val(DM)
    hasEDM = initialized(EDM)
    if ( hasEDM ) E => val(EDM)

    i1 = DMidx
    i2 = i1
    if ( present(EDMidx) ) i2 = EDMidx

    if ( lis_eq ) then

      GF => mum%RHS_SPARSE(:)

      if ( calc_q ) then
        q = 0._dp
        s => spar(spS)
        Sg => val(spS)
        call attach(s, n_col=s_ncol, list_ptr=s_ptr, list_col=s_col)
      end if

      if ( calc_q .and. hasEDM ) then

!$OMP parallel do default(shared), private(ir,jo,ind,io,Hn,ind_H,s_ptr_begin,s_ptr_end,sin)
        do ir = 1 , mum%NRHS
          
          ! this is column index
          jo = ts2s_orb(ir)

          ! The update region equivalent GF part
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1

            io = ts2s_orb(mum%IRHS_SPARSE(ind))

            Hn    = l_ncol(io)
            ind_H = l_ptr(io)
            ! Requires that l_col is sorted
            ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
            if ( ind_H > l_ptr(io) ) then

              s_ptr_begin = s_ptr(io) + 1
              s_ptr_end = s_ptr(io) + s_ncol(io)

              ! Search for overlap index
              ! spS is transposed, so we have to conjugate the
              ! S value, then we may take the imaginary part.
              sin = s_ptr_begin - 1 + SFIND(s_col(s_ptr_begin:s_ptr_end), jo)
              if ( sin >= s_ptr_begin ) q = q - aimag(GF(ind) * Sg(sin))

              D(ind_H,i1) = D(ind_H,i1) - dimag( GF(ind) * DMfact  )
              E(ind_H,i2) = E(ind_H,i2) - dimag( GF(ind) * EDMfact )
            end if

          end do
        end do
!$OMP end parallel do

      else if ( hasEDM ) then

!$OMP parallel do default(shared), private(ir,jo,ind,io,Hn,ind_H)
        do ir = 1 , mum%NRHS
          jo = ts2s_orb(ir)
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1
            io = ts2s_orb(mum%IRHS_SPARSE(ind))
            Hn    = l_ncol(io)
            ind_H = l_ptr(io)
            ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
            if ( ind_H > l_ptr(io) ) then
              D(ind_H,i1) = D(ind_H,i1) - dimag( GF(ind) * DMfact  )
              E(ind_H,i2) = E(ind_H,i2) - dimag( GF(ind) * EDMfact )
            end if
          end do
        end do
!$OMP end parallel do

      else if ( calc_q ) then

!$OMP parallel do default(shared), private(ir,jo,ind,io,Hn,ind_H,s_ptr_begin,s_ptr_end,sin)
        do ir = 1 , mum%NRHS
          jo = ts2s_orb(ir)
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1
            io = ts2s_orb(mum%IRHS_SPARSE(ind))
            Hn    = l_ncol(io)
            ind_H = l_ptr(io)
            ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
            if ( ind_H > l_ptr(io) ) then
              s_ptr_begin = s_ptr(io) + 1
              s_ptr_end = s_ptr(io) + s_ncol(io)
              sin = s_ptr_begin - 1 + SFIND(s_col(s_ptr_begin:s_ptr_end), jo)
              if ( sin >= s_ptr_begin ) q = q - aimag(GF(ind) * Sg(sin))
              D(ind_H,i1) = D(ind_H,i1) - dimag( GF(ind) * DMfact  )
            end if
          end do
        end do
!$OMP end parallel do

      else

!$OMP parallel do default(shared), private(ir,jo,ind,io,Hn,ind_H)
        do ir = 1 , mum%NRHS
          jo = ts2s_orb(ir)
          do ind = mum%IRHS_PTR(ir) , mum%IRHS_PTR(ir+1)-1
            io = ts2s_orb(mum%IRHS_SPARSE(ind))
            Hn    = l_ncol(io)
            ind_H = l_ptr(io)
            ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
            if ( ind_H > l_ptr(io) ) then
              D(ind_H,i1) = D(ind_H,i1) - dimag( GF(ind) * DMfact  )
            end if
          end do
        end do
!$OMP end parallel do

      end if

    else ! lis_eq

      GF => mum%A(:)

      if ( hasEDM ) then

!$OMP parallel do default(shared), private(ind,io,jo,Hn,ind_H)
        do ind = 1 , mum%NZ ! looping A
          
          ! collect the two indices
          io = ts2s_orb(mum%JCN(ind)) ! this is row index for Gf.G.Gf^\dagger
          jo = ts2s_orb(mum%IRN(ind)) ! this is column index for Gf.G.Gf^\dagger

          Hn    = l_ncol(io)
          ind_H = l_ptr(io)
          ! Check that the entry exists
          ! Requires that l_col is sorted
          ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
                    
          if ( ind_H > l_ptr(io) ) then ! this occurs as mum%A contains
                                        ! the electrode as well
          
            D(ind_H,i1) = D(ind_H,i1) + real( GF(ind) * DMfact  ,dp)
            E(ind_H,i2) = E(ind_H,i2) + real( GF(ind) * EDMfact ,dp)

          end if
             
        end do
!$OMP end parallel do
       
      else

!$OMP parallel do default(shared), private(ind,io,jo,Hn,ind_H)
        do ind = 1 , mum%NZ
          io = ts2s_orb(mum%JCN(ind))
          jo = ts2s_orb(mum%IRN(ind))
          Hn    = l_ncol(io)
          ind_H = l_ptr(io)
          ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)
          if ( ind_H > l_ptr(io) ) then
            D(ind_H,i1) = D(ind_H,i1) + real( GF(ind) * DMfact ,dp)
          end if
        end do
!$OMP end parallel do

      end if

    end if

    ! For ts_dq we should not multiply by 2 since we don't do G + G^\dagger for Gamma-only
    ! this is because G is ensured symmetric for Gamma-point and thus it is not needed.
    ! So here the weights are not scaled

  end subroutine add_DM


  ! creation of the GF^{-1}.
  ! this routine will insert the zS-H and \Sigma_{LR} terms in the GF 
  subroutine prepare_invGF(cE, mum, N_Elec, Elecs, spH, spS)

    use intrinsic_missing, only : SFIND
    use class_dSpData1D
    use class_Sparsity
    use m_ts_electype
    use m_ts_cctype, only : ts_c_idx
    use m_ts_method, only : ts2s_orb
    include 'zmumps_struc.h'

    ! the current energy point
    type(ts_c_idx), intent(in) :: cE
    type(zMUMPS_STRUC), intent(inout) :: mum
    integer, intent(in) :: N_Elec
    type(Elec), intent(in) :: Elecs(N_Elec)
    ! The Hamiltonian and overlap sparse matrices
    type(dSpData1D), intent(inout) :: spH,  spS

    ! Local variables
    complex(dp) :: Z
    type(Sparsity), pointer :: sp
    integer, pointer :: l_ncol(:), l_ptr(:), l_col(:)
    real(dp), pointer :: H(:), S(:)
    complex(dp), pointer :: iG(:)
    integer :: io, jo, ind, Hn, ind_H

    if ( cE%fake ) return

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',1)
#endif

    Z = cE%e
    
    sp => spar(spH)
    H => val(spH)
    S => val(spS)

    l_ncol => n_col(sp)
    l_ptr => list_ptr(sp)
    l_col => list_col(sp)

    ! Initialize
    iG => mum%A(:)

!$OMP parallel default(shared), private(ind,io,jo,Hn,ind_H)

!$OMP do
    do ind = 1, mum%NZ

       iG(ind) = 0._dp

       io = ts2s_orb(mum%JCN(ind))
       jo = ts2s_orb(mum%IRN(ind))

       Hn = l_ncol(io)
       if ( Hn /= 0 ) then

       ind_H = l_ptr(io)
       ! Requires that l_col is sorted
       ind_H = ind_H + SFIND(l_col(ind_H+1:ind_H+Hn),jo)

       if ( ind_H > l_ptr(io) ) then
       
          ! Notice that we transpose S and H back here
          ! See symmetrize_HS_Gamma (H is hermitian)
          iG(ind) = Z * S(ind_H) - H(ind_H)

       end if
       end if

    end do
!$OMP end do

    do io = 1 , N_Elec
       call insert_Self_Energies(mum, Elecs(io))
    end do

!$OMP end parallel

#ifdef TRANSIESTA_TIMING
    call timer('TS-prep',2)
#endif

  end subroutine prepare_invGF
   
#endif
end module m_ts_mumpsg
