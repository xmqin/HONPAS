! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_siesta_forces
      implicit none
      private

      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
#ifdef MPI      
      use mpi_siesta
#endif      
      use precision, only: dp
      use files, only: slabel
      use siesta_cml
#ifdef SIESTA__FLOOK
      use flook_siesta, only : slua_call
      use flook_siesta, only : LUA_INIT_MD, LUA_SCF_LOOP
      use siesta_dicts, only : dict_variable_add
      use m_ts_options, only : ts_scf_mixs
      use variable, only: cunpack
#ifndef NCDF_4
      use dictionary, only: assign
#endif      
      use m_mixing, only: mixers_history_init
#endif
      use m_state_init
      use m_setup_hamiltonian
      use m_setup_H0
      use m_compute_dm
      use m_compute_max_diff
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer, only: mixer
      use m_mixing_scf, only: mixing_scf_converged
      use m_mixing_scf, only: mixers_scf_history_init
      use m_mixing_scf, only: scf_mixs, scf_mix

      use m_rhog,                only: mix_rhog, compute_charge_diff
      use siesta_options
      use parallel,     only : IOnode, SIESTA_worker
      use m_state_analysis
      use m_steps
      use sys, only : bye
      use sparse_matrices, only: DM_2D, S_1D
      use sparse_matrices, only: H, Hold, Dold, Dscf, Eold, Escf
      use m_spin, only: spin
      use atomlist, only: Qtot
      use m_dm_charge, only: dm_charge
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
      use siesta_geom,   only: na_u           ! Number of atoms in unit cell
      use m_energies,    only: Etot           ! Total energy
      use m_forces,      only: fa, cfa        ! Forces and constrained forces
      use m_stress,      only: cstress        ! Constrained stress tensor
      use siesta_master, only: forcesToMaster ! Send forces to master prog
      use siesta_master, only: siesta_server  ! Is siesta a server?
      use m_save_density_matrix, only: save_density_matrix
      use m_iodm_old,    only: write_spmatrix
!
      use units,         only: eV, Ang
      use m_pexsi_solver,        only: prevDmax
      use m_forces,              only: fa
      use units,                 only: eV, Ang
      use write_subs,            only: siesta_write_forces
      use write_subs,            only: transiesta_write_forces
      use write_subs,            only: siesta_write_stress_pressure
#ifdef NCDF_4
      use m_ncdf_siesta, only : cdf_save_settings
#endif
      use m_compute_energies, only: compute_energies
      use m_mpi_utils, only: broadcast
      use fdf
      use m_check_walltime
      use m_mpi_utils, only: barrier
      
#ifdef SIESTA__PEXSI      
      use m_pexsi, only: pexsi_finalize_scfloop
#endif     

      use m_ts_options, only : N_Elec
      use m_ts_method
      use m_ts_global_vars,      only: TSmode, TSinit, TSrun
      use siesta_geom,           only: nsc, na_u, xa, ucell, isc_off
      use atomlist,              only: no_u, lasto
      use sparse_matrices,       only: sparse_pattern, block_dist
      use sparse_matrices,       only: Escf, S, maxnh
      use m_ts_options,          only: N_mu, mus
      use ts_dq_m,               only: TS_DQ_METHOD
      use ts_dq_m,               only: TS_DQ_METHOD_FERMI
      use ts_dq_m,               only: TS_DQ_FERMI_TOLERANCE
      use ts_dq_m,               only: TS_DQ_FERMI_SCALE
      use ts_dq_m,               only: ts_dq
      use m_transiesta,          only: transiesta
      use kpoint_grid, only : gamma_scf
      use m_energies, only : Ef


#ifdef HAVE_LIBINT
      use atomlist,      only: iaorb, indxua, iphorb, no_l, no_s
      use siesta_geom,   only: isa, na_s
      use sparse_matrices, only: listh, listhptr, numh
      use nao2gto_types, only: l_max
      use nao2gto_data, only: hfx_gtos, hfx_libint
      use nao2gto_data, only: hfx_options, hfx_system
      use nao2gto_data, only: hfx_libderiv
      use nao2gto_parallel
      use nao2gto_prescreen, only: init_prescreen_eri
      use nao2gto_wrappers, only: nao2gto_system_init
      use nao2gto_wrappers, only: nao2gto_system_free
      use nao2gto_wrappers, only: nao2gto_libint_init
      use nao2gto_wrappers, only: nao2gto_libderiv_init
#endif

      implicit none

      integer, intent(in)  :: istep
      integer :: iscf
      logical :: first, last_step
      logical :: SCFconverged
      real(dp) :: dDmax     ! Max. change in DM elements
      real(dp) :: dHmax     ! Max. change in H elements
      real(dp) :: dEmax     ! Max. change in EDM elements
      real(dp) :: drhog     ! Max. change in rho(G) (experimental)
      real(dp) :: dQ        ! Change charge
      real(dp), target :: G2max      ! the actual mesh-cutoff used
      type(converger_t)  ::  conv_harris, conv_freeE
!---------------------------------------------------------------------- BEGIN
      logical               :: time_is_up
      character(len=40)     :: tmp_str
      
      external :: die, message

      real(dp) :: Qcur
#ifdef SIESTA__FLOOK
      ! len=24 from m_mixing.F90
      character(len=1), target :: next_mixer(24)
      character(len=24) :: nnext_mixer
      integer :: imix
#endif
#ifdef MPI
      integer :: MPIerror
#endif

#ifdef DEBUG
      call write_debug( '    PRE siesta_forces' )
#endif

#ifdef SIESTA__PEXSI      
      ! Broadcast relevant things for program logic
      ! These were set in read_options, called only by "SIESTA_workers".
      call broadcast(nscf,comm=true_MPI_Comm_World)
#endif
      
      if (SIESTA_worker)  then
      !     Initialization tasks for a given geometry
         call state_init( istep )
#ifdef HAVE_LIBINT
      if ( hfx_wanted ) then
         call nao2gto_system_init(hfx_system, maxnh, na_s, no_s, 
     &          spin%Grid, na_u, no_l, no_u, iaorb, indxua, iphorb, 
     &          isa, listh, listhptr, nsc, numh, ucell, xa)
#ifdef MPI
         call MPI_Barrier(MPI_Comm_World, MPIerror)
         call Broadcast_NAO2GTO()
#endif
         call nao2gto_libint_init(hfx_libint, l_max)
         call nao2gto_libderiv_init(hfx_libderiv, l_max)
         call init_prescreen_eri(hfx_libint, hfx_options, hfx_system)
       endif
#endif
      endif

#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after state_init")
#endif

      if (fdf_get("Sonly",.false.)) then
         if (SIESTA_worker) then
            call timer( 'all', 2 )
            call timer( 'all', 3 )
         endif
         call bye("S only")
      endif

      if ( TS_DQ_METHOD == TS_DQ_METHOD_FERMI ) then
        ! Initialize the charge correction
        call ts_dq%initialize(N_mu, mus)
      end if

      Qcur = Qtot
      
#ifdef SIESTA__FLOOK
      ! Add the iscf constant to the list of variables
      ! that are available only in this part of the
      ! routine.
      call dict_variable_add('SCF.iteration',iscf)
      call dict_variable_add('SCF.converged',SCFConverged)
      call dict_variable_add('SCF.charge',Qcur)
      call dict_variable_add('SCF.dD',dDmax)
      call dict_variable_add('SCF.dH',dHmax)
      call dict_variable_add('SCF.dE',dEmax)
      call dict_variable_add('SCF.drhoG',drhog)
      ! We have to set the meshcutoff here
      ! because the asked and required ones are not
      ! necessarily the same
      call dict_variable_add('Mesh.Cutoff.Minimum',G2cut)
      call dict_variable_add('Mesh.Cutoff.Used',G2max)

      if ( mix_charge ) then
        call dict_variable_add('SCF.Mixer.Weight',wmix)
      else
        call dict_variable_add('SCF.Mixer.Weight',scf_mix%w)
        call dict_variable_add('SCF.Mixer.Restart',scf_mix%restart)
        call dict_variable_add('SCF.Mixer.Iterations',scf_mix%n_itt)
        ! Just to populate the table in the dictionary
        call dict_variable_add('SCF.Mixer.Switch',next_mixer)
      end if

      ! Initialize to no switch
      next_mixer = ' '
#endif

!     The current structure of the loop tries to reproduce the
!     historical Siesta usage. It should be made more clear.
!     Two changes: 
!
!       -- The number of scf iterations performed is exactly
!          equal to the number specified (i.e., the "forces"
!          phase is not counted as a final scf step)
!
!       -- At the change to a TranSiesta GF run the variable "first"
!          is implicitly reset to "true".

!     Start of SCF loop
      iscf = 0
      
      ! This call computes the non-scf part of H and initializes the
      ! real-space grid structures.  It might be better to split the two,
      ! putting the grid initialization into state_init and moving the
      ! calculation of H_0 to the body of the loop, done if first=.true.  This
      ! would suit "analysis" runs in which nscf = 0
      if (SIESTA_worker) call setup_H0(G2max)

#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after setup_H0")
#endif

#ifdef SIESTA__FLOOK
      ! Communicate with lua, just before entering the SCF loop
      ! This is mainly to be able to communicate
      ! mesh-related quantities (g2max)
      call slua_call(LUA, LUA_INIT_MD)
#endif

#ifdef NCDF_4
      ! Initialize the NC file
      if ( write_cdf ) then

!       Save the settings (important to do here since mesh-cutoff may
!       have changed).
        call cdf_save_settings(trim(slabel)//'.nc')
#ifdef MPI
        call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

      end if
#endif

      ! The dHmax variable only has meaning for Hamiltonian
      ! mixing, or when requiring the Hamiltonian to be converged.
      dDmax = -1._dp
      dHmax = -1._dp
      dEmax = -1._dp
      dQ = 0._dp

      DO 
        ! Conditions of exit:
        !  -- At the top, to catch a non-positive nscf
        !  -- At the bottom, based on convergence or # of iterations

        iscf = iscf + 1
        if ( iscf > nscf ) EXIT

        ! Note implications for TranSiesta when mixing H
        ! Now H will be recomputed instead of simply being
        ! inherited, however, this is required as the 
        ! if we have bias calculations as the electric
        ! field across the junction needs to be present.
        first = (iscf == 1)

        if (SIESTA_worker) then
           if (first) then
           SCFConverged = .false.
              if ( converge_Eharr ) then
                 call reset(conv_harris)
                 call set_tolerance(conv_harris,tolerance_Eharr)
              endif
              if ( converge_FreeE ) then
                 call reset(conv_FreeE)
                 call set_tolerance(conv_FreeE,tolerance_FreeE)
              endif
           endif

           ! Check whether we are short of time to continue
           call check_walltime(time_is_up)
           if (time_is_up) then
           ! Save DM/H if we were not saving it...
           !     Do any other bookeeping not done by "die"
           call timer('all',2)
           call timer('all',3)
             if (.not. SCFConverged) then
              call message('WARNING',
     $             'SCF_NOT_CONV: SCF did not converge'//
     $             ' before wall time exhaustion')
              write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, prevDmax
              call message(' (info)',"Geom step, scf iteration, dmax:"
     $                            //trim(tmp_str))
             endif
           call barrier() ! A non-root node might get first to the 'die' call
             call die("OUT_OF_TIME: Time is up.")
           endif
        
          call timer( 'IterSCF', 1 )
          if (cml_p)
     &         call cmlStartStep( xf=mainXML, type='SCF', index=iscf )

       endif ! SIESTA_worker

       if (mixH) then
          if ( first .and. SIESTA_worker) then
             if (fdf_get("Read-H-from-file",.false.)) then
                call get_H_from_file()
             else
                call setup_hamiltonian( iscf )
             endif
          endif
          
          call compute_dm( iscf )
          if (SIESTA_worker) then
             ! Maybe set Dold to zero if reading charge or H...
             call compute_max_diff(Dold,Dscf,dDmax)
             if ( converge_EDM ) 
     &            call compute_max_diff(Eold,Escf,dEmax)
             call setup_hamiltonian( iscf )
             call compute_max_diff(Hold,H,dHmax)
          end if
       else
          if (SIESTA_worker) then
             call setup_hamiltonian( iscf )
             call compute_max_diff(Hold,H,dHmax)
          end if
          call compute_dm( iscf )
          if (SIESTA_worker) then
             call compute_max_diff(Dold,Dscf,dDmax)
             if ( converge_EDM ) 
     &            call compute_max_diff(Eold,Escf,dEmax)
          end if
       endif


       if (SIESTA_worker) then

         ! Calculate current charge based on the density matrix
         call dm_charge(spin, DM_2D, S_1D, Qcur)
         dQ = Qcur - Qtot

         call compute_energies(iscf)
         if (mix_charge) then
           call compute_charge_diff(drhog)
         endif
         ! Note: For DM and H convergence checks. At this point:
         ! If mixing the DM:
         !        Dscf=DM_out, Dold=DM_in(mixed), H=H_in, Hold=H_in(prev step)
         !        dDmax=maxdiff(DM_out,DM_in)
         !        dHmax=maxdiff(H_in - H_in(prev step))
         ! If mixing the Hamiltonian:
         !        Dscf=DM_out, Dold=DM_in, H=H_(DM_out), Hold=H_in(mixed)
         !        dDmax=maxdiff(DM_out,DM_in)
         !        dHmax=maxdiff(H(DM_out),H_in)
         call scfconvergence_test( first, iscf, 
     &                             dDmax, dHmax, dEmax, dQ,
     &                             conv_harris, conv_freeE,
     &                             SCFconverged )
        !
        ! ** Check this heuristic
        if ( MixH ) then
           prevDmax = dHmax
        else
           prevDmax = dDmax
        endif
        
        if (SCFconverged .and. (iscf < min_nscf)) then
           SCFconverged = .false.
           if ( IONode ) then
              write(6,"(a,i0)")
     $         "SCF cycle continued for minimum number of iterations: ",
     $            min_nscf
           end if
        end if
        
        ! In case the user has requested a Fermi-level correction
        ! Then we start by correcting the fermi-level
        if ( TSrun .and. TS_DQ_METHOD == TS_DQ_METHOD_FERMI ) then
          ! Signal for next SCF
          ts_dq%run = .true.
          if ( converge_DM )
     &        ts_dq%run = ts_dq%run .and.
     &        dDtol * TS_DQ_FERMI_SCALE > dDmax
          if ( converge_H )
     &        ts_dq%run = ts_dq%run .and.
     &        dHtol * TS_DQ_FERMI_SCALE > dHmax
          if ( converge_EDM )
     &        ts_dq%run = ts_dq%run .and.
     &        tolerance_EDM * TS_DQ_FERMI_SCALE > dEmax
          if ( abs(dQ) > TS_DQ_FERMI_TOLERANCE ) then
            if ( IONode .and. SCFconverged ) then
              write(6,"(2a)") "SCF cycle continued due ",
     &            "to TranSiesta charge deviation"
            end if
            SCFconverged = .false.
          end if
        end if

        ! Check whether we should step to the next mixer
        call mixing_scf_converged( SCFconverged )

        if (monitor_forces_in_scf) call compute_forces()
        
        ! Mix_after_convergence preserves the old behavior of
        ! the program.
        if ((.not. SCFconverged) .or. mix_after_convergence) then
           ! Mix for next step
           if ( mix_charge ) then
              call mix_rhog(iscf)
           else
              call mixer( iscf )
           endif
           ! Save for possible restarts
           if (MixH) then
             call write_spmatrix(H,file="H_MIXED",when=writeH)
             call save_density_matrix(SCFconverged,
     &           file="DM_OUT",when=writedm)
           else
             call save_density_matrix(SCFconverged,
     &           file="DM_MIXED",when=writedm)
             call write_spmatrix(H,file="H_DMGEN",when=writeH)
           endif
        endif

        call timer( 'IterSCF', 2 )
        call print_timings( first, istep == inicoor )
        if (cml_p) call cmlEndStep(mainXML)

#ifdef SIESTA__FLOOK
        ! Communicate with lua
        call slua_call(LUA, LUA_SCF_LOOP)

        ! Retrieve an easy character string
        nnext_mixer = cunpack(next_mixer)
        if ( len_trim(nnext_mixer) > 0 .and. .not. mix_charge ) then
           if ( TSrun ) then
              do imix = 1 , size(ts_scf_mixs)
                 if ( ts_scf_mixs(imix)%name == nnext_mixer ) then
                    call mixers_history_init(ts_scf_mixs)
                    scf_mix => ts_scf_mixs(imix)
                    exit
                 end if
              end do
           else 
              do imix = 1 , size(scf_mixs)
                 if ( scf_mixs(imix)%name == nnext_mixer ) then
                    call mixers_history_init(scf_mixs)
                    scf_mix => scf_mixs(imix)
                    exit
                 end if
              end do
           end if

           ! Check that we indeed have changed the mixer
           if ( IONode .and. scf_mix%name /= nnext_mixer ) then
              write(*,'(2a)') 'siesta-lua: WARNING: trying to change ',
     &             'to a non-existing mixer! Not changing anything!'
              
           else if ( IONode ) then
              write(*,'(2a)') 'siesta-lua: Switching mixer method to: ',
     &             trim(nnext_mixer)

           end if
           ! Reset for next loop
           next_mixer = ' '

           ! Update the references
           call dict_variable_add('SCF.Mixer.Weight',scf_mix%w)
           call dict_variable_add('SCF.Mixer.Restart',scf_mix%restart)
           call dict_variable_add('SCF.Mixer.Iterations',scf_mix%n_itt)
           
        end if
#endif
        
        ! This should be the last step if we have converged...
        last_step = SCFConverged

        ! ... except that we might continue for TranSiesta
        if ( last_step ) then
           call transiesta_switch() ! might reset 'last_step'
                                    ! and iscf
        end if

        if ( iscf == nscf ) then
           last_step = .true.
        end if

       endif ! SIESTA_worker

#ifdef SIESTA__PEXSI       
       call broadcast(last_step,comm=true_MPI_Comm_World)
       call broadcast(iscf,comm=true_MPI_Comm_World)
#endif
       if ( last_step ) exit

      end do

#ifdef SIESTA__PEXSI      
      if (isolve == SOLVE_PEXSI) then
         call pexsi_finalize_scfloop()
      endif
#endif            

      ! Clean up the charge correction object
      call ts_dq%delete()
      
      if (.not. SIESTA_worker) RETURN

      call end_of_cycle_save_operations(SCFconverged)
      
      if (.not.SCFconverged) then
         if (SCFMustConverge) then
            call message('FATAL','SCF_NOT_CONV: SCF did not converge' //
     $           ' in maximum number of steps (required).')
            write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, prevDmax
            call message(' (info)',"Geom step, scf iteration, dmax:"
     $          //trim(tmp_str))
            if ( TSrun ) then
              write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, dQ
              call message(' (info)',"Geom step, scf iteration, dq:"
     $            //trim(tmp_str))
            end if
            call timer( 'all', 2 ) ! New call to close the tree
            call timer( 'all', 3 )
            call barrier()
            call die('ABNORMAL_TERMINATION')
         else
            call message('WARNING',
     $           'SCF_NOT_CONV: SCF did not converge'//
     $           ' in maximum number of steps.')
            write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, prevDmax
            call message(' (info)',"Geom step, scf iteration, dmax:"
     $                            //trim(tmp_str))
            if ( TSrun ) then
              write(tmp_str,"(i5,1x,i5,f12.6)") istep, iscf, dQ
              call message(' (info)',"Geom step, scf iteration, dq:"
     $            //trim(tmp_str))
            end if
         endif
      end if

      if ( TSmode.and.TSinit.and.(.not. SCFConverged) ) then
         call die('SCF did not converge before proceeding to transiesta
     &calculation')
      end if

      ! Clean-up here to limit memory usage
      call mixers_scf_history_init( )

      
!     End of standard SCF loop.
!     Do one more pass to compute forces and stresses

!     Note that this call will no longer overwrite H while computing the
!     final energies, forces and stresses...

      if (fdf_get("compute-forces",.true.)) then
         call post_scf_work( istep, iscf , SCFconverged )
#ifdef SIESTA__PEXSI
         if (ionode) call memory_snapshot("after post_scf_work")
#endif
      endif

!     ... so H at this point is the latest generator of the DM, except
!     if mixing H beyond self-consistency or terminating the scf loop
!     without convergence while mixing H

      call state_analysis( istep )
#ifdef SIESTA__PEXSI
      if (ionode) call memory_snapshot("after state_analysis")
#endif

      ! If siesta is running as a subroutine, send forces to master program
      if (siesta_server) 
     &  call forcesToMaster( na_u, Etot, cfa, cstress )

#ifdef HAVE_LIBINT
      if ( hfx_wanted ) then
         call nao2gto_system_free(hfx_system)
      end if
#endif

#ifdef DEBUG
      call write_debug( '    POS siesta_forces' )
#endif

      CONTAINS

!-----------------------------------------------------
      subroutine get_H_from_file()
      use sparse_matrices, only: maxnh, numh, listh, listhptr
      use atomlist,        only: no_l
      use m_spin,          only: spin
      use m_iodm_old,      only: read_spmatrix

      logical :: found

      call read_spmatrix(maxnh, no_l, spin%H, numh,
     .     listhptr, listh, H, found, userfile="H_IN")
      if (.not. found) call die("Could not find H_IN")
      end subroutine get_H_from_file

!-----------------------------------------------------
      subroutine compute_forces()
      ! Computes forces and stresses with the current DM_out

      use siesta_options, only: recompute_H_after_scf
      use m_final_H_f_stress, only: final_H_f_stress
      use write_subs

      real(dp), allocatable  :: fa_old(:,:), Hsave(:,:)

      allocate(fa_old(size(fa,dim=1),size(fa,dim=2)))
      fa_old(:,:) = fa(:,:)
      if ( recompute_H_after_scf ) then
         allocate(Hsave(size(H,dim=1),size(H,dim=2)))
         Hsave(:,:) = H(:,:)
      end if
      call final_H_f_stress( istep , iscf , .false. )
      if ( recompute_H_after_scf ) then
         H(:,:) = Hsave(:,:)
         deallocate(Hsave)
      end if
      if (ionode) then
        write(6,'(a,f11.6)') "Max diff in force (eV/Ang): ",
     $      maxval(abs(fa-fa_old))*Ang/eV
        call siesta_write_forces(-1)
        if ( TSrun ) then
          call transiesta_write_forces()
        end if
        call siesta_write_stress_pressure()
      endif
      deallocate(fa_old)

      end subroutine compute_forces

      subroutine print_timings(first, first_md)
      use timer_options, only: use_tree_timer
      use m_ts_global_vars, only : TSrun
      logical, intent(in) :: first, first_md
      character(len=20) :: routine

      ! If this is not the first iteration,
      ! we immediately return.
      if ( .not. first ) return
      if ( .not. first_md ) return

      routine = 'IterSCF'

      if ( TSrun ) then
         ! with Green function generation
         ! The tree-timer requires direct
         ! children of the routine to be 
         ! queried.
         ! This is not obeyed in the TS case... :(
         if ( .not. use_tree_timer ) then
            routine = 'TS'
         end if
      endif
      call timer( routine, 3 ) 

      end subroutine 
      
      !
      ! Depending on various conditions, save the DMin
      ! or the DMout, and possibly keep a copy of H

      ! NOTE: Only if the scf cycle converged before exit it
      ! is guaranteed that the DM is "pure out" and that
      ! we can recover the right H if mixing H.
      !
      subroutine end_of_cycle_save_operations(SCFconverged)
      logical, intent(in) :: SCFconverged

      if (mix_after_convergence) then
            ! If we have been saving them, there is no point in doing
            ! it one more time
         if (MixH) then
            call save_density_matrix(SCFconverged, file="DM_OUT",
     $          when=((.not. writedm) .and. write_dm_at_end_of_cycle))
            call write_spmatrix(H,file="H_MIXED", 
     $          when=((.not. writeH) .and. write_H_at_end_of_cycle))
         else
            call save_density_matrix(SCFconverged, file="DM_MIXED",
     $          when=((.not. writedm) .and. write_dm_at_end_of_cycle))
            call write_spmatrix(H,file="H_DMGEN", 
     $          when=((.not. writeH) .and. write_H_at_end_of_cycle))
         endif
      else
         call save_density_matrix(SCFconverged, file="DM_OUT",
     $                            when=write_dm_at_end_of_cycle)
         if (MixH) then
             ! Make sure that we keep the H actually used 
             ! to generate the last DM, if needed.
             H = Hold
         endif
         call write_spmatrix(H,file="H_DMGEN", 
     $                       when=write_H_at_end_of_cycle)
      endif

      end subroutine end_of_cycle_save_operations

!-----------------------------------------------------
      subroutine transiesta_switch()

      use precision,             only: dp
      use parallel,              only: IONode
      use class_dSpData2D
      use class_Fstack_dData1D
      use densematrix, only: resetDenseMatrix

      use siesta_options,        only: fire_mix, broyden_maxit
      use siesta_options,        only: dDtol, dHtol

      use sparse_matrices, only : DM_2D, EDM_2D
      use atomlist, only: lasto
      use siesta_geom, only: nsc, isc_off, na_u, xa, ucell
      use m_energies, only : Ef
      use m_mixing, only: mixers_history_init
      use m_mixing_scf, only: scf_mix, scf_mixs
      use m_rhog, only: resetRhoG

      use m_ts_global_vars,      only: TSinit, TSrun
      use m_ts_global_vars,      only: ts_print_transiesta
      use m_ts_method
      use m_ts_options,          only: N_Elec, Elecs
      use m_ts_options,          only: val_swap
      use m_ts_options,          only: ts_Dtol, ts_Htol
      use m_ts_options,          only: ts_hist_keep
      use m_ts_options,          only: ts_siesta_stop
      use m_ts_options,          only: ts_scf_mixs
      use m_ts_electype

      integer :: iEl, na_a
      integer, allocatable :: allowed_a(:)
      real(dp), pointer :: DM(:,:), EDM(:,:)

      ! We are done with the initial diagon run
      ! Now we start the TRANSIESTA (Green functions) run
      if ( .not. TSmode ) return
      if ( .not. TSinit ) return

      ! whether we are in siesta initialization step
      TSinit = .false.
      ! whether transiesta is running
      TSrun  = .true.

      ! If transiesta should stop immediately
      if ( ts_siesta_stop ) then
         
        if ( IONode ) then
          write(*,'(a)') 'ts: Stopping transiesta (user option)!'
        end if
        
        return
         
      end if

      ! Reduce memory requirements
      call resetDenseMatrix()

      ! Signal to continue...
      last_step = .false.
      iscf   = 0

      ! DANGER (when/if going back to the DIAGON run, we should
      ! re-instantiate the original mixing value)
      call val_swap(dDtol,ts_Dtol)
      call val_swap(dHtol,ts_Htol)

      ! Clean up mixing history.
      if ( mix_charge ) then
        call resetRhoG(.true.)
      else
        if ( associated(ts_scf_mixs, target=scf_mixs) ) then
          do iel = 1 , size(scf_mix%stack)
            call reset(scf_mix%stack(iel), -ts_hist_keep)
            ! Reset iteration count as certain 
            ! mixing schemes require this for consistency
            scf_mix%cur_itt = n_items(scf_mix%stack(iel))
          end do
        else
          call mixers_history_init(scf_mixs)
        end if
      end if
      
      ! Transfer scf_mixing to the transiesta mixing routine
      scf_mix => ts_scf_mixs(1)
#ifdef SIESTA__FLOOK
      if ( .not. mix_charge ) then
        call dict_variable_add('SCF.Mixer.Weight',scf_mix%w)
        call dict_variable_add('SCF.Mixer.Restart',scf_mix%restart)
        call dict_variable_add('SCF.Mixer.Iterations',scf_mix%n_itt)
      end if
#endif
      
      call ts_print_transiesta()

      ! In case of transiesta and DM_bulk.
      ! In case we ask for initialization of the DM in bulk
      ! we read in the DM files from the electrodes and 
      ! initialize the bulk to those values
      if ( any(Elecs(:)%DM_init > 0) ) then
         
        if ( IONode ) then
          write(*,'(/,2a)') 'transiesta: ',
     &        'Initializing bulk DM in electrodes.'
        end if
        
        ! The electrode EDM is aligned at Ef == 0
        ! We need to align the energy matrix to Ef == 0, then we switch
        ! it back later.
        DM => val(DM_2D)
        EDM => val(EDM_2D)
        iEl = size(DM)
        call daxpy(iEl,-Ef,DM(1,1),1,EDM(1,1),1)
        
        na_a = 0
        do iEl = 1 , na_u
          if ( a_isBuffer(iEl) ) then
            na_a = na_a + 1
          else if ( a_isDev(iEl) ) then
            ! do nothing, not allowed overwriting
          else if ( Elecs(atom_type(iEl))%DM_init > 0 ) then
            na_a = na_a + 1
          end if
        end do
        allocate(allowed_a(na_a))
        na_a = 0 
        do iEl = 1 , na_u
          if ( a_isBuffer(iEl) ) then
            na_a = na_a + 1
            allowed_a(na_a) = iEl
          else if ( a_isDev(iEl) ) then
            ! do nothing, not allowed overwriting
          else if ( Elecs(atom_type(iEl))%DM_init > 0 ) then
            na_a = na_a + 1
            allowed_a(na_a) = iEl
          end if
        end do
        
        do iEl = 1 , N_Elec
          if ( Elecs(iEl)%DM_init == 0 ) cycle
          
          if ( IONode ) then
            write(*,'(/,2a)') 'transiesta: ',
     &          'Reading in electrode DM for '//
     &          trim(Elecs(iEl)%Name)
          end if
          
          ! Copy over the DM in the lead
          ! Notice that the EDM matrix that is copied over
          ! will be equivalent at Ef == 0
          call copy_DM(Elecs(iEl),na_u,xa,lasto,nsc,isc_off, 
     &        ucell, DM_2D, EDM_2D, na_a, allowed_a)
           
        end do
        
        ! Clean-up
        deallocate(allowed_a)
        
        if ( IONode ) then
          write(*,*) ! new-line
        end if
        
        ! The electrode EDM is aligned at Ef == 0
        ! We need to align the energy matrix
        iEl = size(DM)
        call daxpy(iEl,Ef,DM(1,1),1,EDM(1,1),1)
        
      end if

      end subroutine transiesta_switch

      END subroutine siesta_forces

      END module m_siesta_forces
