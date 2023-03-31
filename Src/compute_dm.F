! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_compute_dm
      private
      public :: compute_dm

      logical, public, save :: PreviousCallDiagon=.false.

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      use units, only: eV
      USE siesta_options
      use class_dSpData1D, only : val
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, no_u, no_s, indxuo, 
     &                    qtot, Qtots, no_l
      use sys, only: die, bye
      use Kpoint_grid
      use m_energies,   only: Ebs, Ecorrec, Entropy
      use m_energies,   only: Ef, Efs
      use m_rmaxh
      use m_eo
      use m_spin,   only: spin
      use m_diagon, only: diagon
      use parallel, only: IONode
      use parallel, only: SIESTA_worker
      use m_compute_ebs_shift, only: compute_ebs_shift
#ifdef SIESTA__PEXSI
      use m_pexsi_solver,  only: pexsi_solver
#endif
      use m_hsx, only      : write_hs_formatted
#ifdef MPI
      use mpi_siesta
#endif
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
      use m_dminim,        only : dminim
      use m_zminim,        only : zminim
      use m_ordern,        only : ordern
      use m_steps,         only : istp
      use m_normalize_dm,  only : normalize_dm
      use m_ts_global_vars, only : TSmode, TSinit, TSrun
      use m_transiesta,     only : transiesta
      
      implicit none
      
!     Input variables
      integer, intent(in) :: iscf

      real(dp)            :: delta_Ebs, delta_Ef
      logical             :: CallDiagon
      integer             :: nnz
      real(dp), pointer :: H_kin(:)

      ! e1>e2 to signal that we do not want DOS weights
      real(dp), parameter :: e1 = 1.0_dp, e2 = -1.0_dp
      real(dp)            :: buffer1
      integer             :: mpierr

!-------------------------------------------------------------------- BEGIN

      if (SIESTA_worker) call timer( 'compute_dm', 1 )

#ifdef MPI
      call MPI_Bcast(isolve,1,MPI_integer,0,true_MPI_Comm_World,mpierr)
#endif     

      if (SIESTA_worker) then
        ! Save present density matrix
        call dcopy(size(Dscf), Dscf(1,1), 1, Dold(1,1), 1)
        if ( converge_EDM ) then
          call dcopy(size(Escf), Escf(1,1), 1, Eold(1,1), 1)
        end if
      end if

      ! Compute shift in Tr(H*DM) for fermi-level bracketting
      ! Use the current H, the previous iteration H, and the
      ! previous iteration DM

      if (SIESTA_worker) then
         if (iscf > 1) then
            call compute_Ebs_shift(Dscf,H,Hold,delta_Ebs)
            delta_Ef = delta_Ebs / qtot
            if (ionode.and.isolve.eq.SOLVE_PEXSI) then
               write(6,"(a,f16.5)")
     $              "Estimated change in band-structure energy:",
     $              delta_Ebs/eV, "Estimated shift in E_fermi: ",
     $              delta_Ef/eV
            endif
               
         else
            delta_Ebs = 0.0_dp
            delta_Ef  = 0.0_dp
         endif
      endif

#ifdef SIESTA__PEXSI
      if (isolve .eq. SOLVE_PEXSI) then
        ! This test done in node 0 since NonCol and SpOrb
        ! are not set for PEXSI-solver-only processes
         if (ionode) then
            if (spin%NCol .or. spin%SO) call die(
     $               "The PEXSI solver does not implement "//
     $               "non-coll spins or Spin-orbit yet")
         endif
         if (ionode) then
            ! This should never happen for large-scale calculations...
            if (no_s /= no_u) call die(
     $               "The PEXSI solver cannot work " //
     $               "with an auxiliary supercell")
         endif
         call pexsi_solver(iscf, no_u, no_l, spin%spinor,
     $              maxnh, numh, listhptr, listh,
     $              H, S, qtot, Dscf, Escf,
     $              ef, Entropy, temp, delta_Ef)
      endif
      if (.not. SIESTA_worker) RETURN
#endif
      ! Here we decide if we want to calculate one or more SCF steps by
      ! diagonalization before proceeding with the OMM routine
      CallDiagon=.false.
      if (isolve .eq. SOLVE_MINIM) then
        if (istp .eq. 1) then
          if ((iscf .le. call_diagon_first_step) .or.
     &        (call_diagon_first_step<0)) CallDiagon=.true.
        else
          if ((iscf .le. call_diagon_default) .or.
     &        (call_diagon_default<0)) CallDiagon=.true.
        endif
      endif

      if (isolve .eq. MATRIX_WRITE) then
!             write(indexstr,'(I15)') iscf
!             write(filename,fnameform) 'H_', trim(adjustl(indexstr)), 
!      &                                '.matrix'
!             call write_global_matrix( no_s, no_l, maxnh, numh, listh,
!      &           H(1:maxnh,1), filename )
!  
!             write(filename,fnameform) 'S_', trim(adjustl(indexstr)), 
!      &                                '.matrix'

!        Note: only one-shot for now
         call write_hs_formatted(no_u, spin%H,
     $        maxnh, numh, listhptr, listh, H, S)
         call bye("End of run after writing H.matrix and S.matrix")

c$$$        call write_global_matrix_singlenodewrite( 
c$$$     &           no_u, no_s, maxnh, numh, listhptr, listh, 
c$$$     &           H(:,1), 'H.matrix')
c$$$
c$$$        call write_global_matrix_singlenodewrite( 
c$$$     &           no_u, no_s, maxnh, numh, listhptr, listh, 
c$$$     &           S, 'S.matrix')

      elseif ((isolve .eq. SOLVE_DIAGON) .or. (CallDiagon)) then
        call diagon(no_s, spin%spinor, 
     &              no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              H, S, qtot, fixspin, qtots, temp, e1, e2,
     $              xijo, indxuo, gamma_SCF, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted,
     &              dealloc_psi=.not. CallDiagon)
        Ecorrec = 0.0_dp
        PreviousCallDiagon=.true.
      elseif (isolve .eq. SOLVE_ORDERN) then
        if ( .not. gamma_SCF ) call die("Cannot do O(N) with k-points.")
        if ( spin%NCol .or. spin%SO )
     .      call die("Cannot do O(N) with non-coll spins or Spin-orbit")
        call ordern(usesavelwf, ioptlwf, na_u, no_u, no_l, lasto,
     &               isa, qa, rcoor, rmaxh, ucell, xa, iscf,
     &               istp, ncgmax, etol, eta, qtot, maxnh, numh,
     &               listhptr, listh, H, S, chebef, noeta, rcoorcp,
     &               beta, pmax, Dscf, Escf, Ecorrec, spin%H, qtots )
        Entropy = 0.0_dp
      elseif (isolve .eq. SOLVE_MINIM) then
        if ( spin%NCol .or. spin%SO ) 
     &      call die('ERROR: Non-collinear spin calculations
     &                       not yet implemented with OMM!')
        H_kin => val(H_kin_1D)

        ! Decide which version of OMM to use.
        ! Test based on use of auxiliary supercell
        ! It might still be possible to avoid the complex version
        
        if ( no_u == no_s ) then    ! Not using an auxiliary supercell
          call dminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 spin%H, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qtots, H, S, H_kin)
        else
          ! When using an auxiliary supercell
          ! (even for gamma point; not optimized yet)
          call zminim(.false., PreviousCallDiagon, iscf, istp, no_l,
     &                 spin%H, no_u, maxnh, numh, listhptr, listh, Dscf,
     &                 eta, qtots, no_s, xijo, indxuo, nkpnt, kpoint,
     &                 kweight, H, S, H_kin)
        end if
        Ecorrec = 0.0_dp
        Entropy = 0.0_dp
        PreviousCallDiagon=.false.
      elseif (TSmode .and. TSinit) then
        call diagon(no_s, spin%spinor, 
     &              no_l, maxnh, maxnh, no_u,
     &              numh, listhptr, listh, numh, listhptr, listh,
     &              H, S, qtot, fixspin, qtots, temp, e1, e2,
     $              xijo, indxuo, gamma_SCF, nkpnt, kpoint, kweight,
     &              eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &              occtol, iscf, neigwanted)

        Ecorrec = 0._dp

      else if (TSrun) then

         call transiesta(iscf,spin%H, block_dist, sparse_pattern,
     &        no_u == no_s, ucell, nsc, isc_off, no_u, na_u,
     &        lasto, xa, maxnh,
     &        H, S, Dscf, Escf, Ef, Qtot)

         Ecorrec = 0._dp
         Entropy = 0.0_dp
         
      else
        !call die('siesta: ERROR: wrong solution method')
      endif

#ifdef CDF
      if ( writedmhs_cdf_history) then
        call write_dmh_netcdf( no_l, maxnh, spin%H, Dold, H, Dscf )
      else if (writedmhs_cdf) then
        call write_dmh_netcdf( no_l, maxnh, spin%H, Dold, H, Dscf,
     &                         overwrite=.true. )
      endif
#endif

!     Normalize density matrix to exact charge
!     Placed here for now to avoid disturbing EHarris
      if ( .not. TSrun ) then
         call normalize_dm( first= .false. )
      end if

      call timer( 'compute_dm', 2 )
#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after compute_DM")
#endif

!-----------------------------------------------------------------------END
      END subroutine compute_dm
      END MODULE m_compute_dm
