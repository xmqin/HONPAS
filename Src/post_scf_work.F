! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_post_scf_work
      private
      public :: post_scf_work

      CONTAINS

      subroutine post_scf_work( istep, iscf, SCFconverged )

      USE m_final_H_f_stress, only: final_H_f_stress
      USE siesta_options
      use sparse_matrices, only: Dscf, Escf, maxnh
      use sparse_matrices, only: H, S
      use sparse_matrices, only: xijo, listh, listhptr, numh
      use siesta_geom

      use sparse_matrices, only: DM_2D, DM_history
      use class_Geometry,  only: newGeometry, Geometry, delete
      use class_Pair_Geometry_dSpData2D, only: Pair_Geometry_dSpData2D
      use class_Pair_Geometry_dSpData2D, only: new, delete
      use class_Fstack_Pair_Geometry_dSpData2D, only: new, max_size
      use class_Fstack_Pair_Geometry_dSpData2D, only: push

      use parallel, only: ionode
      use atomlist, only: lasto, rmaxo, datm, indxuo, no_s, no_u,
     &                    iphorb, no_l, qtot, qtots
      use m_energies
      use neighbour,   only : maxna=>maxnna   ! For plcharge...
      use m_spin,         only : nspin, h_spin_dim
      use m_spin,         only : spinor_dim 
      use m_diagon,       only : diagon 
      use m_dminim,       only : dminim
      use m_zminim,       only : zminim
      use m_steps,        only : istp, inicoor
      use m_compute_dm,   only : PreviousCallDiagon
      use m_eo
      use Kpoint_grid
      implicit none

      ! MD-step, SCF-step
      integer, intent(in) :: istep , iscf
      logical, intent(in) :: SCFconverged

      character(len=20) :: msg
      type(Pair_Geometry_dSpData2D)  :: pair
      type(Geometry)                :: geom

      call timer( 'PostSCF', 1 )

!     If converged, make one last iteration to find forces and stress

!     If we use the minimization routine, the energy-density
!     matrix is not calculated inside the SCF step, and
!     so this must be done now

!     Which Hamiltonian is used in this call? It must be the H
!     that generated the DM, for consistency

      if ((isolve .eq. SOLVE_MINIM) .and.
     &    (.not. PreviousCallDiagon)) then
        if (minim_calc_eigenvalues) then
           if (MixH .and. mix_after_convergence) then
              ! probably better to have a logical 'h_is_pure'
              if (ionode) then
                 write(6,"(a)") ":!: SOLVE_MINIM: " //
     $                "Generating Escf with a modified H..."
              endif
           endif
          call diagon(no_s, spinor_dim, 
     &                no_l, maxnh, maxnh, no_u,
     &                numh, listhptr, listh, numh, listhptr, listh, 
     &                H, S, qtot, fixspin, qtots, temp, 1.0_dp, -1.0_dp,
     &                xijo, indxuo, gamma_SCF, nkpnt, kpoint, kweight,
     &                eo, qo, Dscf, Escf, ef, efs, Entropy, no_u,
     &                occtol, iscf, neigwanted)
          Ecorrec = 0.0_dp
        else
            if ( no_u == no_s ) then ! Not using an auxiliary supercell
            call dminim(.true., .false., iscf, istp, no_l, nspin, no_u,
     &                  maxnh, numh, listhptr, listh, Escf, eta, qtots)
          else
            ! When using an auxiliary supercell
            ! (even for gamma point; not optimized yet)
            call zminim(.true., .false., iscf, istp, no_l, nspin, no_u,
     &                  maxnh, numh, listhptr, listh, Escf, eta, qtots,
     &                  no_s, xijo, indxuo, nkpnt, kpoint, kweight)
          end if
        endif
      endif

!     Find final energy, forces and stress.

!     Note that here we are using the Dscf coming out of the
!     scf cycle, which could be DM_out or DM_mixed, depending
!     on whether we are mixing H or DM, and on whether we
!     mix one final time after convergence

      if (ionode) then
         if (mix_after_convergence .and. (.not. MixH)
     $                             .and. (.not. mix_charge)) then
            write(6,"(/,a)")
     $      ":!: Using DM_mixed to compute the final energy and forces"
            write(6,"(a)")
     $      "Consider using 'SCF.MixAfterConvergence F', mixing H" //
     $           " or mixing the charge density"
         else
            write(6,"(/,a)")
     $      "Using DM_out to compute the final energy and forces"
         endif
      endif

      call final_H_f_stress( istep , iscf , SCFconverged )

      ! These energies are computed with a different DM from those in the
      ! last SCF step, except if mixing the Hamiltonian, in which case they
      ! are the same.

      call update_DEna()
      call update_Etot()

!     Since this Etot is computed with DM_mixed (the DM predicted by
!     'mixer' for a hypothetical next SCF iteration) or with DM_out
!     (depending on whether a final mixing is done at the end of the
!     cycle) this value for Eharrs is suspect. This value was never
!     actually printed, since printing routines used 'Eharrs1', which is
!     the last value of Eharrs computed in the SCF cycle (set in
!     'scfconvergence_test'). To reproduce that behavior, the following
!     statement has been removed, and Eharris1 replaced by Eharrs
!     everywhere.

      !! Eharrs = Etot + DEharr
      !
      ! Update the Free electronic energy with the new Etot
      ! and the latest entropy. This is not completely 
      ! consistent if using a mixed DM to compute Etot, but
      ! is a reasonable approximation
      call update_FreeE( Temp )
      
      
      if (dumpcharge) then
        call plcharge( no_s, na_s, no_u, maxnh, maxna, h_spin_dim,
     &       isa, iphorb, indxuo, lasto,
     &       scell, nsc, xa, rmaxo, datm )
      endif

      ! Store the density matrix to the coordinate stack
      if ( max_size(DM_history) > 0 ) then

      ! In case the run is an FC run we do not 
      ! save any additional things in the stack, thus
      ! the extraction will be the "zero" structure
      if ( (idyn==6 .or. idyn==7).and. istep > inicoor ) then
         ! do nothing, we already have the first DM
      else            
         ! Save (geom,DM) pair to history stack
         write(msg,"(a,i0)") "step: ",istep
         call newGeometry(geom,na_u,ucell,xa,isa,
     $                    name="Geometry for " // trim(msg))
         call new(pair,geom,DM_2D,name="Geom, DM for "//trim(msg))
         call push(DM_history,pair)
         call delete(pair)
         call delete(geom)
      end if
      
      end if

      call timer( 'PostSCF', 2 )

#ifdef DEBUG
      call write_debug( '    POS post_scf_work' )
#endif

      END subroutine post_scf_work

      END MODULE m_post_scf_work
