! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_scfconvergence_test
      private
      public :: scfconvergence_test
      CONTAINS

      subroutine scfconvergence_test( first, iscf,
     &                                dDmax, dHmax, dEmax, dQ,
     $                                conv_harris, conv_freeE,
     $                                converged)
      use precision, only: dp
      USE siesta_options
      use siesta_cml
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
      use write_subs
      use m_energies
      use units, only: eV
      use m_convergence, only: converger_t, tolerance
      use m_convergence, only: add_value, is_converged
      use dftu_specs,    only: switch_dftu, dTol_dftupop
      use m_dftu,        only: DFTU_dPop

      use m_ts_global_vars, only: TSrun
      use m_ts_options,  only: ts_dQtol, ts_converge_dQ

      use sparse_matrices, only:  Dscf, S, numh, listhptr, listh, maxnh
      use siesta_geom, only: na_u, isa
      use atomlist, only:    no_u, lasto, iaorb, iphorb, indxuo
      use m_spin,           only: h_spin_dim, SpOrb

      implicit none

      integer :: iscf
      logical :: first

      real(dp), intent(in) :: dDmax ! Max. change in dens. matrix
      real(dp), intent(in) :: dHmax ! Max. change in H
      real(dp), intent(in) :: dEmax ! Max. change in EDM
      real(dp), intent(in) :: dQ ! Max. change in charge (mainly for TS)

      type(converger_t), intent(inout)  :: conv_harris, conv_freeE
      logical, intent(out)              :: converged

      logical :: is_DFTU_conv
      real(dp) :: dummy_qspin(8)
      character(len=64) :: conv_text

!------------------------------------------------------------------- BEGIN

      ! convergence test

      call add_value(conv_harris, eharrs)
      call add_value(conv_freeE, freeE)

      ! Print energies
      if (IOnode) then
        call siesta_write_energies( iscf, dDmax, dHmax )

        if (harrisfun) then
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ',Eharrs/eV
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='SCF Cycle')
            call cmlAddProperty(xf=mainXML, value=Eharrs/eV,
     .       units="siestaUnits:eV", dictRef="siesta:Eharrs",
     .       fmt="r7")
            call cmlEndPropertyList(mainXML)
          endif
        endif

        ! flush stdout
        call pxfflush(6)
        call wallclock("-------------- end of scf step")
      endif

!     Print spin polarization at each SCF step if requested before mixing
      if (spndeb) then
         call print_spin(dummy_qspin)
      endif
!     Print populations at each SCF step if requested before mixing ......
      if (muldeb) then
        if (ionode) write (6,"(/a)") 'Using DM_out for analysis:'
        if ( SpOrb .and. orbmoms) then
           call moments( 1, na_u, no_u, maxnh, numh, listhptr,
     .           listh, S, Dscf, isa, lasto, iaorb, iphorb,
     .           indxuo )
        endif
        ! Call this unconditionally
        call mulliken( mullipop, na_u, no_u, maxnh,
     &                 numh, listhptr, listh, S, Dscf, isa,
     &                 lasto, iaorb, iphorb )
      endif

      is_DFTU_conv = switch_dftu .and.
     &     DFTU_dPop < dTol_dftupop

      ! Assume initial convergence
      ! Since 4.1 this logic has changed to be fully determined
      ! by the user.
      ! The defaults are however as in pre 4.1 versions so there
      ! should be no surprises.
      converged = .true.
      conv_text = ' '
      if ( converge_Eharr ) then
        conv_text = trim(conv_text)//'+Harris'
        converged = converged .and. is_converged(conv_harris)
      end if
      if ( converge_FreeE ) then
        conv_text = trim(conv_text)//'+FreeE'
        converged = converged .and. is_converged(conv_FreeE)
      end if
      if ( converge_EDM ) then
        conv_text = trim(conv_text)//'+EDM'
        converged = converged .and. dEmax < tolerance_EDM
      end if
      if ( converge_DM ) then
        conv_text = trim(conv_text)//'+DM'
        converged = converged .and. dDmax < dDtol
      end if
      if ( converge_H ) then
        conv_text = trim(conv_text)//'+H'
        converged = converged .and. dHmax < dHtol
      end if
      if ( TSrun .and. ts_converge_dQ ) then
        conv_text = trim(conv_text)//'+dQ'
        converged = converged .and. abs(dQ) < ts_dQtol
      end if

      if (converged .and. IOnode) then

         ! Remove the initial '+'
         conv_text = conv_text(2:)

         write(6,"(/,3a)") "SCF Convergence by ",trim(conv_text),
     &        " criterion"

         if ( converge_Eharr ) then
            write(6,"(a,f14.8)") "|EH_i-EH_(i-1)| (eV) < ",
     &           tolerance(conv_harris)/eV
         end if

         if ( converge_FreeE ) then
            write(6,"(a,f14.8)") "|FreeE_i-FreeE_(i-1)| (eV) < ",
     &           tolerance(conv_freeE)/eV
         end if

         ! No matter what we print out the differences of DM and H
         if ( mix_charge ) then
            write(6,"(a,f16.10)") "max |DM_i - DM_(i-1)|        : ",
     &           dDmax
            if ( dEmax >= 0._dp )
     &      write(6,"(a,f16.10)") "max |EDM_i - EDM_(i-1)| (eV) : ",
     &              dEmax/eV
         else
            write(6,"(a,f16.10)") "max |DM_out - DM_in|         : ",
     &           dDmax
            if ( dEmax >= 0._dp )
     &      write(6,"(a,f16.10)") "max |EDM_out - EDM_in|  (eV) : ",
     &              dEmax/eV
         endif
         if ( dHmax >= 0._dp ) then
            write(6,"(a,f16.10)") "max |H_out - H_in|      (eV) : ",
     &           dHmax/eV
         end if
         if ( TSrun ) then
           write(6,"(a,f16.10)") "     Q(TS) - Q(dev)          : ",
     &         dQ
         end if

	 if ( is_DFTU_conv ) then
           ! Inform of LDA+U population convergence
           write(6,"(a)")
     &         "LDA+U population converged by dTol_pop criteria:"
           write(6,"(a,f16.10)")
     &         "max |pop(LDA+U)_i - pop(LDA+U)_j|: ", DFTU_dPop
         end if

         write(6,"(a,i0,a)") "SCF cycle converged after ",
     &           iscf," iterations"

         call pxfflush(6)

      endif

      if (harrisfun) then
         converged = .true.
         if (IOnode) then
            write(6,"(a)") "Harris-functional calculation " //
     $                     "considered 'converged'"
            call pxfflush(6)
         endif
      endif

!---------------------------------------------------------------- END
      END subroutine scfconvergence_test
      end module m_scfconvergence_test
