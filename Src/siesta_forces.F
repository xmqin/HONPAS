      module m_siesta_forces
      private
      
      public :: siesta_forces

      CONTAINS

      subroutine siesta_forces(istep)
      use siesta_cml
      use m_state_init
      use m_setup_hamiltonian
      use m_compute_dm
      use m_scfconvergence_test
      use m_post_scf_work
      use m_mixer
      USE siesta_options
      use parallel, only: IOnode
      use m_state_analysis
      use m_steps
      use sys, only : die
      use m_convergence, only: converger_t
      use m_convergence, only: reset, set_tolerance
#ifdef TRANSIESTA
      use m_ts_options, only : mixH
      use m_ts_global_vars, only : ts_istep
#endif /* TRANSIESTA */


      implicit none

      integer, intent(in)  :: istep

      integer  :: iscf
      logical  :: first, last, SCFconverged
      real(dp) :: dDmax     ! Max. change in density matrix elements
      type(converger_t)  ::  conv_harris, conv_etot

#ifdef TRANSIESTA
! Used in force constant runs
      ts_istep = istep
#endif /* TRANSIESTA */
!------------------------------------------------------------------------- BEGIN

        ! Initialization tasks for a given geometry
        call state_init( istep )

        ! The current structure of the loop tries to reproduce the
        ! historical Siesta usage. It should be made more clear.

        ! Start of SCF loop

        iscf  = 1
        first = .true.
        last  = .false.
        SCFconverged = .false.
        if (require_harris_convergence) then
           call reset(conv_harris)
           call set_tolerance(conv_harris,harris_tolerance)
        endif
        if (require_energy_convergence) then
           call reset(conv_etot)
           call set_tolerance(conv_etot,energy_tolerance)
        endif

        ! The SCF loop is executed at least once, 
        ! except if nscf <= 0.
        if (nscf <= 0) then
           last = .true.
        endif

        DO WHILE (.not. last)
          call timer( 'IterSCF', 1 )

#ifndef TRANSIESTA
          if (cml_p)
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)

          call setup_hamiltonian( first, last, iscf )

          call compute_dm( iscf )

          call mixer( iscf, dDmax )
#else /* TRANSIESTA */

          if (cml_p)
     $         call cmlStartStep(xf=mainXML, type='SCF', index=iscf)
! Adding the possibility of Mixinh H.
          if (mixH) then

            if ( first ) then
             call setup_hamiltonian( first, last, iscf )
            end if

            call compute_dm( iscf )

            call setup_hamiltonian( first, last, iscf )

          else ! Mix DM

            call setup_hamiltonian( first, last, iscf )

            call compute_dm( iscf )

          end if ! mixH

          call mixer( iscf, dDmax )
#endif /* TRANSIESTA */

          call scfconvergence_test( first, iscf, last, dDmax,
     $                              conv_harris, conv_etot )
          ! At this point last is true if SCF has converged. But this
          ! changed below
          SCFconverged = last

          call timer( 'IterSCF', 2 )
          if (istep.eq.inicoor .and. first) call timer( 'IterSCF', 3 )

          iscf = iscf + 1
          first = .false.

          ! The last scheduled iteration is used to compute forces
          ! and stresses, hence the exit if iscf==nscf.
          ! The ">" clause catches the nscf==1 case. For this,

          if (iscf >= nscf) then
            last = .true.
          endif

          if (cml_p) call cmlEndStep(mainXML)

        ENDDO

        ! If we have not converged, but get to here, we are out
        ! of steps. We may want to error out.
        if ((SCFMustConverge).and.(.not.SCFconverged)) then
          call die('SCF did not converge in maximum number of steps.')
        endif

        ! End of standard SCF loop.
        ! Do one more pass to compute forces and stresses

        if (last) then
          call post_scf_work( first, last, iscf )
        else
          STOP 'Bad logic in SCF loop'
        endif

        call state_analysis( istep )

!--------------------------------------------------------------------------- END
        END subroutine siesta_forces
      END module m_siesta_forces
