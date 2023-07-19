      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( first, iscf, last, dDmax,
     $                                conv_harris, conv_etot)
      USE siesta_options
      use siesta_cml
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
      use write_subs
      use m_energies
      use m_convergence, only: converger_t
      use m_convergence, only: add_value, is_converged
#ifdef TRANSIESTA
      use m_ts_global_vars, only: TSinit, TSrun, TSiscf
      use m_pulay, only: init_pulay_arrays
#endif /* TRANSIESTA */

      implicit none

      integer :: iscf
      logical :: first
      logical :: last

      real(dp), intent(in) :: dDmax     ! Max. change in density matrix elements
      type(converger_t), intent(inout)  :: conv_harris, conv_etot
!------------------------------------------------------------------------- BEGIN
          ! convergence test
          DEna = Enascf - Enaatm
          Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec+Emad+Emm+
     &           Emeta
          Eharrs = Etot + DEharr

          call add_value(conv_harris, eharrs)
          call add_value(conv_etot, etot)

          FreeE  = Etot - Temp * Entropy

          ! Recalculating the energy in the last iter (for gridcellsampling)
          ! but preserving the value of Eharrs1
          Eharrs1 = Eharrs

          Elast = Etot
          if (require_energy_convergence) then
             if (   dDmax.lt.dDtol
     $            .and. is_converged(conv_etot)) then
                last = .true.
                if (IOnode) then
                   write(6,"(a)") "SCF Convergence by Etot+DM criteria"
                endif
             endif
          else if (require_harris_convergence) then
             if (is_converged(conv_harris)) then
                last = .true.
                if (IOnode) then
                   write(6,"(a)") "SCF Convergence by Harris criterion"
                endif
             endif
          else
            if (dDmax.lt.dDtol) last = .true.
          endif


          ! Print energies
          if (IOnode) then
            call siesta_write_energies( first, iscf, dDmax )

            if (harrisfun) then
              write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ',Eharrs/eV
              if (cml_p) then
                call cmlStartPropertyList(mainXML, title='SCF Cycle')
                call cmlAddProperty(xf=mainXML, value=Eharrs/eV,
     .           units="siestaUnits:eV", dictRef="siesta:Eharrs", 
     .           fmt="r7")
                call cmlEndPropertyList(mainXML)
              endif
            endif
           
            ! flush stdout
            call pxfflush(6)
            call wallclock("-------------- end of scf step")
          endif
#ifdef TRANSIESTA
! TSS Begin
        if ((isolve. eq. 2) .and. last .and. TSinit) then
           last = .false.       ! we're only done with the initial diagon run
           TSinit = .false.      ! now we start the TRANSIESTA run
           TSrun = .true.
           !
           ! From now on, mixing is done on potentials,
           ! so we have to reset the Pulay record keeping
           !
           call init_pulay_arrays()  
           !
           if(IONode) then
              write(6,'(/a)')
     .  '                          ************************       '
              write(6,'(a)')
     .  '                          *   TRANSIESTA BEGIN   *       '
          write(6,'(a)')
     .  '                          ************************       '

          end if
        end if
        if (TSiscf .eq. 1) call timer( 'IterSCF', 3 ) !TSS
        if (TSiscf .eq. 1) call timer( 'TS', 3 ) !TSS
        if (last) TSiscf = 0

! TSS End
#endif /* TRANSIESTA */

!--------------------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
