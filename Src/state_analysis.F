! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_state_analysis
      use write_subs

      private
      public :: state_analysis

      CONTAINS

      subroutine state_analysis( istep )
      use siesta_cml
      use m_born_charge, only : born_charge
      use parallel,      only : IOnode
      use m_wallclock,   only : wallclock
      use zmatrix,       only : lUseZmatrix, iofaZmat,
     &                          CartesianForce_to_ZmatForce
      use atomlist,      only : iaorb, iphorb, amass, no_u, lasto
      use atomlist,      only : indxuo, Qtot
      use m_spin,        only : nspin, SpOrb
      use m_fixed,       only : fixed
      use sparse_matrices
      use siesta_geom

      USE siesta_options
      use units, only: amu, eV
      use m_stress
      use m_energies, only: Etot, FreeE, Eharrs, FreeEHarris, Entropy
      use m_energies, only: Ebs, Ef
      use m_ntm
      use m_forces
      use m_energies, only: update_FreeE, update_FreeEHarris
      use m_intramol_pressure, only: remove_intramol_pressure
      use m_ts_global_vars, only: TSrun
      use m_ts_options, only: N_elec, Elecs
      use ts_charge_m, only: ts_charge_print, TS_Q_INFO_FULL

#ifdef SIESTA__FLOOK
      use flook_siesta, only : slua_call, LUA_FORCES
#endif

      implicit none

      integer  :: istep
      integer  :: ia, jx, ix
      real(dp) :: volume
      logical  :: eggbox_block=.true. ! Read eggbox info from data file?
      real(dp) :: qspin
      
      external :: eggbox, mulliken, moments
      real(dp), external :: volcel

!------------------------------------------------------------------------- BEGIN
      call timer( 'state_analysis', 1 )
#ifdef DEBUG
      call write_debug( '  PRE state_analysis' )
#endif

      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='SCF Finalization')   
      endif

!     Write final Kohn-Sham and Free Energy
      
      FreeE       = Etot - Temp * Entropy
      FreeEHarris = Eharrs - Temp * Entropy

      if (cml_p) call cmlStartPropertyList(mainXML,
     &                        title='Energies and spin')
      if (IOnode) then
        if ( .not. harrisfun)
     &      write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
        if (cml_p) then
           call cmlAddProperty(xf=mainXML, value=Etot/eV,
     &       dictref='siesta:E_KS', units='siestaUnits:eV', 
     .       fmt='r6')
           call cmlAddProperty(xf=mainXML, value=FreeE/eV,
     &       dictref='siesta:FreeE', units='siestaUnits:eV', 
     .       fmt='r6')
           call cmlAddProperty(xf=mainXML, value=Ebs/eV,
     &       dictref='siesta:Ebs', units='siestaUnits:eV', 
     .       fmt='r6')
           call cmlAddProperty(xf=mainXML, value=Ef/eV,
     &       dictref='siesta:E_Fermi', units='siestaUnits:eV', 
     .       fmt='r6')
        endif
      endif

!     Substract egg box effect from energy
      if (eggbox_block) then
        call eggbox( 'energy', ucell, na_u, isa, ntm, xa, fa, Etot,
     &               eggbox_block )
        FreeE  = Etot - Temp * Entropy
        if (IOnode)
     &    write(6,"(/a,f14.4)") 'siesta: E_KS - E_eggbox = ',Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, value=Etot/eV,
     &         dictref='siesta:E_KS_egg', units='siestaUnits:eV', 
     .         fmt='r6')
      endif

      call update_FreeE( Temp )
      call update_FreeEHarris( Temp )
      call print_spin(qspin)
      
      if (cml_p) call cmlEndPropertyList( mainXML )

!     Substract egg box effect from the forces 
      if (eggbox_block) then
        call eggbox('forces',ucell,na_u,isa,ntm,xa,fa,Etot,eggbox_block)
      endif

      if (IOnode) call write_raw_efs(stress,na_u,fa,FreeE)

!     Compute stress without internal molecular pressure
      call remove_intramol_pressure(ucell,stress,na_u,xa,fa,mstress)

!     Impose constraints to atomic movements by changing forces ...........
      if (RemoveIntraMolecularPressure) then
!        Consider intramolecular pressure-removal as another
!        kind of constraint
         call fixed( ucell, mstress, na_u, isa, amass, xa, fa,
     &               cstress, cfa, ntcon , 
     &               magnitude_usage = idyn==0 )
      else
         call fixed( ucell, stress, na_u, isa, amass, xa, fa,
     &               cstress, cfa, ntcon ,
     &               magnitude_usage = idyn==0 )
      endif

#ifdef SIESTA__FLOOK
      ! We call it right after using the
      ! geometry constraints.
      ! In that way we can use both methods on top
      ! of each other!
      ! The easy, already implemented methods in fixed,
      ! and custom ones in Lua :)
      call slua_call(LUA, LUA_FORCES)
#endif

!     Calculate and output Zmatrix forces
      if (lUseZmatrix .and. (idyn.eq.0)) then
        call CartesianForce_to_ZmatForce(na_u,xa,fa)
        if (IOnode) call iofaZmat()
      endif

!     Compute kinetic contribution to stress
      kin_stress(1:3,1:3) = 0.0_dp
      volume = volcel(ucell)
      do ia = 1,na_u
        do jx = 1,3
          do ix = 1,3
            kin_stress(ix,jx) = kin_stress(ix,jx) -
     &             amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
          enddo
        enddo
      enddo
!     Add kinetic term to stress tensor
      tstress = stress + kin_stress

!     Force output 
      if (IOnode) then
        call siesta_write_forces(istep)
        if ( TSrun ) then
          call transiesta_write_forces()
        end if
        call siesta_write_stress_pressure()
        call wallclock('--- end of geometry step')
      endif

!     Population and moment analysis 
      if ( SpOrb .and. orbmoms) then
         call moments( 1, na_u, no_u, maxnh, numh, listhptr,
     .           listh, S, Dscf, isa, lasto, iaorb, iphorb,
     .           indxuo )
      endif
      ! Call this unconditionally
      call mulliken( mullipop, na_u, no_u, maxnh,
     &               numh, listhptr, listh, S, Dscf, isa,
     &               lasto, iaorb, iphorb )

      ! Also write TS charges
      if ( TSrun ) then
        call ts_charge_print(N_Elec,Elecs, Qtot,
     &      block_dist, sparse_pattern,
     &      nspin, maxnh, Dscf, S, TS_Q_INFO_FULL)
      end if

!
!     Call the born effective charge routine only in those steps (even) 
!     in which the dx  is positive.
      if (bornz .and. (mod(istep,2) .eq. 0)) then
         call born_charge()
      endif

!     End the xml module corresponding to the analysis
      if (cml_p) then
         call cmlEndModule(mainXML)         
      endif 
      call timer( 'state_analysis', 2 )

!--------------------------------------------------------------------------- END
      END subroutine state_analysis

      END MODULE m_state_analysis
