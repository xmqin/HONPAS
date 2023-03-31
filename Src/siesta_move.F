! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_move

      private
      public :: siesta_move

      CONTAINS

      subroutine siesta_move( istep, relaxd )

      use zmatrix,         only: lUseZmatrix, iofaZmat, write_Zmatrix,
     .                           CartesianForce_to_ZmatForce
      use atomlist,        only: iza, amass
      use m_ioxv,          only: ioxv
      USE siesta_options
      use units,           only: Kelvin
      use parallel,        only: IOnode
      use siesta_cml,      only: cml_p, cmlEndStep, mainXML
      use sys,             only: die

      use m_broyden_optim, only: broyden_optimizer
      use m_fire_optim, only: fire_optimizer
      use m_zm_broyden_optim,      only: zm_broyden_optimizer
      use m_zm_fire_optim,      only: zm_fire_optimizer
      use m_cell_broyden_optim,      only: cell_broyden_optimizer
      use m_cell_fire_optim,      only: cell_fire_optimizer

      use m_dynamics,      only: nose, verlet2, npr, anneal, pr

      use m_energies,      only: Ekinion
      use write_subs
      use m_steps
      use m_kinetic
      use siesta_geom,     only: na_u       ! Number of atoms in unit cell
      use siesta_geom,     only: na_s       ! Number of atoms in super-cell
      use siesta_geom,     only: isa        ! Atomic species
      use siesta_geom,     only: xa         ! Atomic positions
      use siesta_geom,     only: xa_last    ! Previous atomic positions
      use siesta_geom,     only: va         ! Atomic velocities
      use siesta_geom,     only: ucell      ! Unit cell vectors
      use siesta_geom,     only: ucell_last ! Previous unit cell vectors
      use siesta_geom,     only: scell, nsc ! Super-cell vectors and multip.
      use siesta_geom,     only: scell_last ! Previous super-cell vectors
      use siesta_geom,     only: vcell      ! Time derivative of unit cell vecs
      use m_forces,        only: cfa        ! Constrained forces
      use m_forces,        only: ntcon      ! Number of geometry constraints
      use m_stress,        only: cstress    ! Constrained stress tensor
      use siesta_master,   only: coordsFromMaster ! Get coords from master prog
      use m_target_stress, only: subtract_target_stress

      use m_check_walltime
      use m_mpi_utils,     only: barrier
      
#ifdef NCDF_4
      use m_exp_coord, only : exp_coord_next
#endif
#ifdef SIESTA__FLOOK
      use siesta_dicts, only : dict_variable_add
      use flook_siesta, only : slua_call, LUA_MOVE
#endif

      implicit none

      integer, intent(in)   :: istep
      logical, intent(out)  :: relaxd


      integer               :: ix, iadispl, ixdispl
      integer, parameter    :: iunit = 2   ! Physical-units option for MD:
                                       ! 1=>(eV,Ang), 2=>(Ry,Bohr)
      real(dp)              :: Pint    ! Instantaneous pressure 
      real(dp)              :: eff_stress(3,3)  
      real(dp)              :: tp_pr
                
      logical               :: foundxv  ! dummy for call to ioxv
      logical               :: foundzm  ! dummy for call to iozm

      logical               :: time_is_up
      character(len=40)     :: tmp_str
      
      real(dp), external :: volcel

#ifdef DEBUG
      call write_debug( '    PRE siesta_move' )
#endif

!------------------------------------------------------------------ BEGIN

      call timer( 'siesta_move', 1 )

      ! Save the last geometry for which the density matrix 
      ! (electronic structure) has been calculated

      xa_last(1:3,1:na_s) = xa(1:3,1:na_s)
      ucell_last(1:3,1:3) = ucell(1:3,1:3)
      scell_last(1:3,1:3) = scell(1:3,1:3)

! Move atoms ..........................................................
!
!  ** Regarding the output of coordinates:
!     
!     -- 'state_init' should record the actual coordinates of a step (CG, MD, whatever)
!     -- this routine (siesta_move) should only record checkpointing information, such as
!        the predicted coordinates for a possible next geometry step.
!        Historically, the checkpointing has been done with the (highly confusing)
!        XV file  (routine ioxv), but there is the newer option (useful only for
!        relaxations) of writing a '.STRUCT_NEXT_ITER' file (routine write_positions with move="true")
!     -- 'siesta_analysis' should only write final coordinates.

!      There is a pletora of output routines: (and more in 'siesta_analysis')
!
!      -- ioxv
!      -- iozm
!      -- write_positions  (which calls write_struct and zm_canonical...)
!      -- write_md_record
!
!      They should be rationalized.
!
!      In the following, it is better to put all the logic in the individual blocks, avoiding
!      branching and re-checks of the 'idyn' variable (wich, by the way, should have symbolic
!      values instead of hardwired numerical ones)

      select case(idyn)

      case(0)   
         ! The original "relaxation" case, but note that it also covers
         ! pure single-point calculations (with MD.NumCGsteps = 0)


        ! Pure single-point calculations will not enter this block
        if (nmove .ne. 0) then    ! That is, if requesting "CG" steps

         ! Here we want checkpointing, except if the structure is already relaxed

         ! Note that these routines do not update the coordinates if
         ! the force/stress criterion for relaxation is satisfied

          if (RelaxCellOnly) then
             if (broyden_optim) then
                call cell_broyden_optimizer( na_u, xa, ucell,
     $            cstress, tp, strtol,
     $            varcel, relaxd)
             elseif (fire_optim) then
                call cell_fire_optimizer( na_u, xa, ucell,
     $            cstress, tp, strtol,
     $            varcel, relaxd)
             else
                call die("Cell-only optim needs Broyden or FIRE")
             endif

          else   ! Coordinate relaxation (and maybe cell)

           if (lUseZmatrix) then
             if (broyden_optim) then
                call zm_broyden_optimizer( na_u, xa, ucell, cstress, tp,
     &                                     strtol, varcel, relaxd )
             elseif (fire_optim) then
                call zm_fire_optimizer( na_u, xa, cfa, ucell,
     $               cstress, dxmax, tp, ftol, strtol,
     $               varcel, relaxd)
             else
                call cgvc_zmatrix( na_u, xa, cfa, ucell, cstress,
     $               dxmax, tp, ftol, strtol, varcel,
     $               relaxd, usesavecg )
             endif
           else
             if (broyden_optim) then
                call broyden_optimizer( na_u, xa, cfa, ucell,
     $               cstress, tp, ftol, strtol, varcel, relaxd )
             elseif (fire_optim) then
                call fire_optimizer( na_u, xa, cfa, ucell,
     $               cstress, dxmax, tp, ftol, strtol,
     $               varcel, relaxd )
             else
                call cgvc( na_u, xa, cfa, ucell, cstress,
     $               dxmax, tp, ftol, strtol, varcel,
     $               relaxd, usesavecg )
             endif
           endif
          endif ! RelaxCellOnly

          if (relaxd) then

            ! Will not call ioxv et al in the block below, so that:
            !  - The XV file will contain the xa coords
            !    written in the previous step (if any, otherwise
            !    those written by state_init)
            !  - The .STRUCT_NEXT_ITER file would be that
            !    written in the previous step (if any), and it
            !    will contain the same coords as the .STRUCT_OUT
            !    file produced in state_init in this step.
            !  - The Zmatrix files...
          else
             call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
             if (lUseZmatrix) call iozm('write',ucell,vcell,xa,foundzm)
             ! This writes the "next_iter" STRUCT and canonical zmatrix files
             call siesta_write_positions(moved=.true.)
             ! Save atomic positions and velocities accumulatively and
             ! accumulate coor in Xmol file for animation 
             call write_md_record( istep )
             call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )
          endif

        endif

      case(1)  ! Micro-canonical MD with fixed cell (Verlet)

         relaxd = .false.    ! It is intent(out)
        ! Check convergence for quenching runs (which are
        ! really relaxations)
        ! Avoid calling the routine if relaxed, so that the coordinates
        ! are not changed.
         if ( iquench .ne. 0 ) then
            relaxd = all(abs(cfa(1:3,1:na_u)) < ftol)
            ! We might want to fall back on a different relaxation
            ! scheme when we reach the slow-moving part
            ! (check temp_ion?;  check cfa behavior?)
         endif
         if (.not. relaxd) then
            call verlet2(istp, iunit, iquench, na_u, cfa, dt,
     .       amass, ntcon, va, xa, Ekinion, tempion)
            if (IOnode) then
               write(6,'(/,a,f12.3,a)')
     $              'siesta: Temp_ion =', tempion, ' K'
            endif
            call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
            call siesta_write_positions(moved=.true.)
            call write_md_record( istep )
            call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )
         endif


      case (2) ! Canonical MD with fixed cell (Nose-Hoover)

         relaxd = .false.    ! It is intent(out)
         call nose( istp, iunit, na_u, cfa, tt, dt, amass, mn,
     .              ntcon, va, xa, Ekinion, kn, vn, tempion )
         if (IOnode) then
            write(6,'(/,a,f12.3,a)') 'siesta: Temp_ion =', tempion, ' K'
         endif
         call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
         call siesta_write_positions(moved=.true.)
         call write_md_record( istep )
         call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )


      case (3) ! Micro-canonical variable-cell MD (Parrinello-Rahman)

        ! Check convergence for quenching runs (which are
        ! really relaxations)
        ! Avoid calling the routine if relaxed, so that the coordinates
        ! are not changed.
         if ( iquench .ne. 0 ) then
            ! Allow a general target stress in this relaxation mode
            call subtract_target_stress(cstress,eff_stress)
            relaxd = all(abs(cfa(1:3,1:na_u)) < ftol)
            relaxd = relaxd .AND. (all(abs(eff_stress) < strtol))
            cstress = eff_stress
            tp_pr = 0.0_dp  ! As we are already passing the modified stress
         else
            ! Standard MD mode
            relaxd = .false.    
            tp_pr = tp      ! Keep the original target pressure
         endif
         if (.not. relaxd) then
            call pr(istp, iunit, iquench, na_u, cfa, cstress, tp_pr, dt,
     .           amass, mpr, ntcon, va, xa, vcell, ucell, Ekinion, 
     .           kpr, vpr, tempion, Pint)
            if (IOnode) then
               write(6,'(/,a,f12.3,a)')
     .              'siesta: E_kin PR =', kpr/Kelvin, ' K'
               write(6,'(/,a,f12.3,a)')
     $              'siesta: Temp_ion =', tempion, ' K'
            endif
            call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &             foundxv)
            call siesta_write_positions(moved=.true.)
            call write_md_record( istep )
            call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )
         endif

      case (4)  ! Canonical variable-cell MD (Nose-Parrinello-Rahman)

         relaxd = .false.    
         call npr(istp, iunit, na_u, cfa, cstress, tp, tt, dt,
     .            amass, mn, mpr, ntcon, va, xa, vcell, ucell, 
     .            Ekinion, kn, kpr, vn, vpr, tempion, Pint)
         if (IOnode) then
            write(6,'(/,a,f12.3,a)')
     $           'siesta: Temp_ion =', tempion, ' K'
         endif
         call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &        foundxv)
         call siesta_write_positions(moved=.true.)
         call write_md_record( istep )
         call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )

      case (5)  ! Annealings

         relaxd = .false.    ! It is intent(out)
         call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .       na_u, cfa, cstress, tp, tt, dt, amass, ntcon,
     .       va, xa, ucell, Ekinion, tempion, Pint)
         if (IOnode) then
            write(6,'(/,a,f12.3,a)')
     $           'siesta: Temp_ion =', tempion, ' K'
         endif
         call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &        foundxv)
         call siesta_write_positions(moved=.true.)
         call write_md_record( istep )
         call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )

      case (6) ! Force-constant-matrix calculation

         relaxd = .false.    ! It is intent(out)

        ! Output of coordinates is meaningless here in general,
        ! except maybe to checkpoint a calculation to be restarted
        ! later, by reading the XV file and choosing inicoor and
        ! fincoor appropriately

        ! Save current (displaced) atomic positions
         call ioxv( 'write', ucell, vcell, na_u, isa, iza,  
     &        xa, va, foundxv)

        ! Undo the last atom displacement
        if ( istep > inicoor ) then ! inicoor == 0
          iadispl = (istep-mod(istep-1,6))/6+ia1
          ix = mod(istep-1,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) - dx
        end if

        ! Displace atom by dx
        ! NOTE: MOVE is before istep increment, hence
        ! fincoor == istep will be the last step and 
        ! the initial geometry should be retained.
        if ( istep < fincoor ) then ! fincoor = (ia2-ia1+1)*6
          iadispl = ((istep+1)-mod(istep,6))/6+ia1
          ix = mod(istep,6)+1
          ixdispl = (ix - mod(ix-1,2) +1)/2
          dx = -dx
          xa(ixdispl,iadispl) = xa(ixdispl,iadispl) + dx
        end if

      case (7)  ! PHONON interface --- removed

         relaxd = .true.

      case (8)  ! Server mode

         relaxd = .false.

         ! It can be argued that in server operation the
         ! responsibility of writing coordinates falls on the
         ! client.

          ! Save atomic positions and velocities accumulatively and
          ! accumulate coor in Xmol file for animation 
          call write_md_record( istep )

          ! Get coordinates from driver program 
          call coordsFromMaster( na_u, xa, ucell )
          if (volcel(ucell) < 1.0e-8_dp) then
            call automatic_cell(ucell,scell,na_u,xa,isa,charnet)
          endif

#ifdef NCDF_4
       case (9)                 ! explicit coordinate mode

          relaxd = .false.

          ! Save the md-record
          call write_md_record( istep )
          if ( istep < fincoor ) then
             call exp_coord_next( istep + 1 , na_u, xa)
          else
             relaxd = .true.
          end if
#endif
          
#ifdef SIESTA__FLOOK
       case (10) ! LUA hook coordinates

          relaxd = .false.

          ! Add the quest for relaxed list of variables
          call dict_variable_add('MD.Relaxed',relaxd)

          ! Communicate with lua
          call slua_call(LUA, LUA_MOVE)

          if ( volcel(ucell) < 1.0e-8_dp ) then
             ! In case the user requests a recalculation of the
             ! unit-cell size (for molecules)
             call automatic_cell(ucell,scell,na_u,xa,isa,charnet)
          end if

          if ( relaxd ) then
             ! see comments for idyn == 1
          else
             call ioxv( 'write', ucell, vcell, na_u, isa, iza, xa, va, 
     &            foundxv)
             
             call siesta_write_positions( moved=.true. )
             call write_md_record( istep )
             call superx_if_compat( ucell, nsc, na_u, na_s, xa, scell )
          end if
#endif

       end select

! Output memory use at the end of this geometry step
      if (cml_p) then
         call cmlEndStep(mainXML)
      endif

      call timer( 'siesta_move', 2 )
      call timer( 'IterGeom', 2 )
! End of one MD step - flush stdout
      if (ionode) call pxfflush(6)

        ! Check whether we are short of time to continue
        call check_walltime(time_is_up)
        if (time_is_up) then
           ! Do any other bookeeping not done by "die"
           call timer('all',2)
           call timer('all',3)
           if (.not. relaxd) then
              call message('WARNING','GEOM_NOT_CONV: '//
     $             'Geometry relaxation not converged')
              write(tmp_str,"(a,i5)") 'Geometry step: ', istep
              call message(' (info)', trim(tmp_str))
           endif
           call barrier()
           call die("OUT_OF_TIME: Time is up at end of geometry step. ")
        endif

!--------------------------------------------------------------------------- END

#ifdef DEBUG
      call write_debug( '    POS siesta_move' )
#endif

      contains

      ! For backwards compatibility only

      subroutine superx_if_compat( ucell, nsc, na, maxa, xa, scell )
      use siesta_options, only : compat_pre_v4_dynamics
      use atomlist, only : superx

      integer :: maxa, na, nsc(3)
      real(dp) :: ucell(3,3), xa(3,MAXA), scell(3,3)
      ! Only xa and scell are overwritten upon exit
      ! but superx does not have intent (so we ca not either)

      if (compat_pre_v4_dynamics) then
         call superx( UCELL, NSC, NA, MAXA, XA, SCELL )
      end if
      end subroutine superx_if_compat
      
      end subroutine siesta_move

      end module m_siesta_move
