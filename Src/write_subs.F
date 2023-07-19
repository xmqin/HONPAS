!
!     Some of the functionality of the former write_subs module 
!     has been split into  three sub-modules to work around a compiler bug.
!     The submodules appear first in this file, and module
!     write_subs is at the end.
!
      Module write_subs_positions

      use m_steps

      private
      public :: siesta_write_positions
      CONTAINS
      
      subroutine siesta_write_positions()
      use siesta_geom
      use atomlist, only: elem, iza
      use siesta_cml
      use units, only: Ang
      use zmatrix, only: lUseZmatrix, write_canonical_ucell_and_Zmatrix
      use m_iostruct, only: write_struct
      use siesta_options, only: idyn

      implicit none

!
!      Write out structural information in "crystallography" format,
!      and in "canonical" (see Zmatrix module) zmatrix format.
!      By now the structure might have changed in response to forces
!      and stresses.

        call write_struct( ucell, na_u, isa, iza, xa, moved=.true.)
        if (lUseZmatrix .and. (idyn .eq. 0)) then
           call write_canonical_ucell_and_Zmatrix
     $                         (filename="NEXT_ITER.UCELL.ZMATRIX")
        endif

        if (cml_p) then
          call cmlAddMolecule(xf=mainXML, natoms=na_u, elements=elem,
     .           atomRefs=cisa, coords=xa/Ang)
          call cmlAddLattice(xf=mainXML, cell=ucell/Ang, 
     .           units='siestaUnits:Ang', dictref='siesta:ucell')
        endif
      end subroutine siesta_write_positions

      END MODULE write_subs_positions
      Module write_subs_pressure

      use m_steps

      private
      public :: siesta_write_stress_pressure

      CONTAINS
      
      subroutine siesta_write_stress_pressure()

      use parallel, only: IOnode
      use precision
      USE siesta_options
      use siesta_geom
      use atomlist, only: iza
      use m_iostruct,   only: write_struct
      use siesta_cml
      use units
      use m_energies, only: FreeE
      use m_stress, only: stress, kin_stress, mstress, cstress, tstress
      use m_forces,   only: fa   ! This does not really belong here  

      implicit none

      integer :: ia, jx, ix

      real(dp):: Pmol       ! "Molecular pressure" (discounting Virial term)
      real(dp):: Psol       ! Pressure of "solid"
      real(dp):: Press      ! Pressure
      real(dp):: ps(3,3)    ! Auxiliary array

! Stress tensor and pressure:
      
      if (.not.final) then
!
!           Write Voigt components of total stress tensor 
!
            ps = stress + kin_stress
            write(6,'(/,a,6f12.2)')
     .           'Stress-tensor-Voigt (kbar):',
     .           (ps(jx,jx)/kbar,jx=1,3),
     $            ps(1,2)/kbar,
     $            ps(2,3)/kbar,
     $            ps(1,3)/kbar
            Press = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
            write(6,"(a,f14.4)") "(Free)E + p*V (eV/cell)",
     $           (FreeE + Press*volume_of_some_cell)/eV

            if (RemoveIntraMolecularPressure) then
               ps = mstress + kin_stress
               write(6,'(/,a,6f12.2)')
     .           'Inter-Molecular-Stress-Voigt (kbar):',
     .           (ps(jx,jx)/kbar,jx=1,3),
     $            ps(1,2)/kbar,
     $            ps(2,3)/kbar,
     $            ps(1,3)/kbar
               Press = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
               write(6,"(a,f14.4)")
     $              "(Free)E + p_inter_molec * V  (eV/cell)",
     $              (FreeE + Press*volume_of_some_cell)/eV
            endif
!
!      This use of the volume is OK, as it is called from state_analysis,
!      before possibly changing the cell.
!      Write "target enthalpy" (E + pV, where p is the *target* pressure)
            write(6,"(a,f14.4)") "Target enthalpy (eV/cell)",
     $           (FreeE + tp*volume_of_some_cell)/eV

      ! Output depends on dynamics option
        select case (idyn)
        case(0:5,8)
          if (idyn==0 .and. (.not.varcel)) then
            continue
          else
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (static) (eV/Ang**3):',
     .           ('     ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((stress(1,1) + stress(2,2) + stress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (static):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=stress*Ang**3, 
     .             dictref='siesta:stress', title='Stress',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol, 
     .             dictref='siesta:psol', title='Pressure (Static)',
     .             units='siestaUnits:kBar')
            endif !cml_p
!
            write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (total) (eV/Ang**3):',
     .           ('     ',(tstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
            Psol = - ((tstress(1,1)+tstress(2,2) +tstress(3,3))/3.0_dp)
            write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (total):', Psol/kBar, '  kBar'
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=tstress*Ang**3, 
     .             dictref='siesta:tstress', title='Total Stress',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol,
     .             dictref='siesta:tpsol', title='Pressure (Total)',
     .             units='siestaUnits:kBar')
            endif !cml_p

            if (RemoveIntraMolecularPressure) then
             ps = mstress
             write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (nonmol) (eV/Ang**3):',
     .           ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
             Psol = - ((ps(1,1) + ps(2,2) + ps(3,3))/3.0_dp)
             write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (nonmol):', Psol/kBar, '  kBar'
             if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=ps*Ang**3, 
     .          dictref='siesta:mstress', 
     .          title='Stress tensor (normal)',
     .          units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol, 
     .          dictref='siesta:pmol', title='Pressure (Nonmol)',
     .          units='siestaUnits:kBar')
             endif                !cml_p
!
             ps = mstress + kin_stress
             write(6,'(/,a,3(/,a,3f12.6))')
     .           'siesta: Stress tensor (nonmol+kin) (eV/Ang**3):',
     .           ('     ',(ps(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
             Psol = - ((ps(1,1)+ps(2,2) +ps(3,3))/3.0_dp)
             write(6,'(/,a,f20.8,a)')
     .           'siesta: Pressure (nonmol+kin):', Psol/kBar, '  kBar'
             if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=ps*Ang**3, 
     .             dictref='siesta:tmstress', 
     .             title='Stress tensor (nonmol+kin)',
     .             units='siestaUnits:evpa3')
              call cmlAddProperty(xf=mainXML, value=Psol,
     .         dictref='siesta:tpmol', title='Pressure (Nonmol+Kin)',
     .         units='siestaUnits:kBar')
             endif                !cml_p
            endif                 ! Remove intramolecular pressure

          endif                  !varcel

        ! Write Force Constant matrix if FC calculation ...
        case(6)
          call ofc(fa,dx,na_u)
        case(7)
!          call phonon_write_forces(fa,na_u,ucell,istep)
           if (IOnode) write(*,*) 'phonon suport deactivated'
        end select !idyn

      else !final

           ! AG: Possible BUG 
           ! The volume here refers to the "old" cell, and
           ! might be out of date if the cell has changed.

! Print stress tensor unconditionally
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(stress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, value=stress*Ang**3/eV, 
     .         dictref='siesta:stress', units='siestaUnits:eV_Ang__3')
        endif !cml_p

! Print constrained stress tensor if different from unconstrained
        if (Any(cstress /= stress )) then
             write(6,'(/,a,3(/,a,3f12.6))')
     .       'siesta: Constrained stress tensor (static) (eV/Ang**3):',
     .       ('siesta: ',(cstress(jx,ix)*Ang**3/eV,jx=1,3),ix=1,3)
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, value=cstress*Ang**3/eV, 
     .           dictref='siesta:cstress', 
     .           units='siestaUnits:eV_Ang__3')
          endif !cml_p
        endif

! Find pressure

        Psol = - (( stress(1,1) + stress(2,2) + stress(3,3) )/3.0_dp)
        Pmol = - (( mstress(1,1) + mstress(2,2) + mstress(3,3) )/3.0_dp)

        write(6,'(/,a,f18.6,a)')
     .       'siesta: Cell volume =',
     $       volume_of_some_cell/Ang**3, ' Ang**3'
        write(6,'(/,a,/,a,2a20,a,3(/,a,2f20.8,a))')
     .       'siesta: Pressure (static):',
     .       'siesta: ','Solid',        'Molecule',      '  Units',
     .       'siesta: ', Psol,           Pmol,           '  Ry/Bohr**3',
     .       'siesta: ', Psol*Ang**3/eV, Pmol*Ang**3/eV, '  eV/Ang**3',
     .       'siesta: ', Psol/kBar,      Pmol/kBar,      '  kBar'
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Final Pressure')
          call cmlAddProperty(xf=mainXML,
     $         value=volume_of_some_cell/Ang**3, 
     .         title='cell volume', dictref='siesta:cellvol', 
     .         units='siestaUnits:Ang__3')
          call cmlAddProperty(xf=mainXML, value=Psol/kBar, 
     .         title='Pressure of Solid', dictref='siesta:pressSol', 
     .         units='siestaUnits:kbar')
          call cmlAddProperty(xf=mainXML, value=Pmol/kBar,       
     .         title='Pressure of Molecule', dictref='siesta:pressMol', 
     .         units='siestaUnits:kbar')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
        
      endif !final for stress & pressure

      end subroutine siesta_write_stress_pressure
      
      
      End Module write_subs_pressure

      Module write_subs_energies

      use m_steps

      private
      public :: siesta_write_energies
      CONTAINS
      
      subroutine siesta_write_energies( first, iscf, dDmax )
      USE siesta_options
      use siesta_cml
      use units
      use m_energies 
      use m_spin
#ifdef TRANSIESTA
      use m_ts_global_vars, only: TSinit, TSrun, TSiscf 
#endif /* TRANSIESTA */
      implicit none

      logical :: first
      integer :: iscf
      real(dp), intent(in) :: dDmax     ! Max. change in density matrix elements
#ifdef TRANSIESTA
      character(11) :: tmpname
#endif /* TRANSIESTA */

      integer :: i 


      ! Only print out full decomposition at very beginning and end.
      if ((istp==1.and.first).or.final) then
        write(6,'(/,a,/,(a,f17.6))')
     .     'siesta: Program''s energy decomposition (eV):',
     .     'siesta: Ebs     =', Ebs/eV, 
     .     'siesta: Eions   =', Eions/eV,
     .     'siesta: Ena     =', Ena/eV,
     .     'siesta: Ekin    =', Ekin/eV,
     .     'siesta: Enl     =', Enl/eV,
     .     'siesta: DEna    =', DEna/eV,
     .     'siesta: DUscf   =', DUscf/eV,
     .     'siesta: DUext   =', DUext/eV,
     .     'siesta: Exc     =', Exc/eV,
     .     'siesta: eta*DQ  =', Ecorrec/eV,
     .     'siesta: Emadel  =', Emad/eV,
     .     'siesta: Emeta   =', Emeta/eV,
     .     'siesta: Emolmec =', Emm/eV,
     .     'siesta: Ekinion =', Ekinion/eV,
     .     'siesta: Eharris =', (Eharrs1+Ekinion)/eV,
     .     'siesta: Etot    =', (Etot+Ekinion)/eV,
     .     'siesta: FreeEng =', (FreeE+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(mainXML,
     .         title='Energy Decomposition')
          call cmlAddProperty(xf=mainXML, value=Ebs/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ebs', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Eions/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Eions', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Ena/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ena', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Ekin/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Enl/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Enl', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DEna/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DEna', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DUscf/eV, 
     .         units='siestaUnits:eV',
     .         dictref='siesta:DUscf', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Exc/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='r6')
          call cmlAddProperty(xf=mainXML,value=Ecorrec/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ecorrec', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emad/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emad', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emeta/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emeta', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=Emm/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Emm', fmt='r6')
          call cmlAddProperty(xf=mainXML,value=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(Eharrs1+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EharrsK', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(Etot+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:EtotK', fmt='r6')
          call cmlAddProperty(xf=mainXML, value=(FreeE+Ekinion)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:FreeEK', fmt='r6')
          call cmlEndPropertyList(mainXML)
        endif
      endif
      ! On all SCF steps, print out the current energy (format depending on type of run)
      if (.not.final) then
        ! Print total energy and density matrix error .........................
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='SCF Cycle')
          ! Eharrs is always output
          call cmlAddProperty(xf=mainXML, value=Eharrs/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:Eharrs", fmt="r7")
        endif
        ! This chain of if statements determines which properties are output.
        if (harrisfun) then 
          write(6,"(/a,f14.6,/)") 'siesta: Eharris(eV) = ', Eharrs/eV
            ! No need for further cml output
        elseif (isolve==0) then
          if (cml_p) 
     .         call cmlAddProperty(xf=mainXML, value=FreeE/eV, 
     .         units="siestaUnits:eV", 
     .         dictRef="siesta:FreeE", fmt="r7")
          if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=Etot/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Etot", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=FreeE/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:FreeE", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=dDmax/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:dDmax",  fmt="r7")
            endif
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,a8,a7,a11)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax, 
     .           (Efs(i)/eV,i=1,2)
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=Efs(1)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=Efs(2)/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Efs", fmt="r7")
            endif
          else !fixspin
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,2a8)')
     .           'siesta: iscf', '   Eharris(eV)', 
     .           '      E_KS(eV)', '   FreeEng(eV)', 
     .           '   dDmax', '  Ef(eV)'
            write(6,'(a8,i4,3f14.4,2f8.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, 
     .           dDmax, Ef/eV
            if (cml_p) then
              call cmlAddProperty(xf=mainXML,value=Ef/eV, 
     .             units="siestaUnits:eV", 
     .             dictRef="siesta:Ef", fmt="r7")
            endif !cml_p
          endif !fixspin
        elseif (isolve==1) then
          write(6,'(/,a15,i4)') 'siesta: iscf = ',iscf
          write(6,'(a14,f15.4,a13,f15.4,a10,f7.4/)') 
     .         'Eharris(eV) = ',Eharrs/eV,
     .         '  E_KS(eV) = ',Etot/eV,'  dDmax = ',dDmax
          if (cml_p) then
            call cmlAddProperty(xf=mainXML, value=Etot/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:Etot", fmt="r7")
            call cmlAddProperty(xf=mainXML, value=dDmax/eV, 
     .           units="siestaUnits:eV", 
     .           dictRef="siesta:dDmax", fmt="r7")
          endif
        endif !harrisfun/isolve

#ifdef TRANSIESTA
! TSS Begin
        if (TSinit) then     !TSS
         if (cml_p)
     .         call cmlAddProperty(xf=mainXML, value=FreeE/eV,
     .         units="siestaUnits:eV",
     .         dictRef="siesta:FreeE",  fmt="r7")
         if (fixspin) then
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=Etot/eV,
     .             units="siestaUnits:eV",
     .             dictRef="siesta:Etot",   fmt="r7")
              call cmlAddProperty(xf=mainXML, value=FreeE/eV,
     .             units="siestaUnits:eV",
     .             dictRef="siesta:FreeE", fmt="r7")
              call cmlAddProperty(xf=mainXML, value=dDmax/eV,
     .             units="siestaUnits:eV",
     .             dictRef="siesta:dDmax",  fmt="r7")
            endif
            if ((iscf .eq. 1).or.muldeb)
     .           write(6,'(/,a12,3a14,a8,a7,a11)')
     .           'siesta: iscf', '   Eharris(eV)',
     .           '      E_KS(eV)', '   FreeEng(eV)',
     .           '   dDmax', '  Ef_up', '  Ef_dn(eV)'
            write(6,'(a8,i4,3f14.4,f8.4,2f9.4)')
     .           'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV, dDmax,
     .           (Efs(i)/eV,i=1,2)
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=Efs(1)/eV,
     .             units="siestaUnits:eV",
     .             dictRef="siesta:Efs",    fmt="r7")
              call cmlAddProperty(xf=mainXML, value=Efs(2)/eV,
     .             units="siestaUnits:eV",
     .             dictRef="siesta:Efs",    fmt="r7")
            endif
         else
          if (iscf .eq. 1) write(6,'(/,a12,3a14,2a8)')
     .             'siesta: iscf', '   Eharris(eV)',
     .             '      E_KS(eV)', '   FreeEng(eV)',
     .             '   dDmax', '  Ef(eV)'
          write(6,'(a8,i4,3f14.4,2f8.4)')
     .             'siesta: ',iscf, Eharrs/eV, Etot/eV, FreeE/eV,
     .             dDmax, Ef/eV
         end if !fixspin
        endif                !TSinit
        if (TSrun) then      !TSS
              tmpname = 'transiesta:'
              if (fixspin) then
                 if (TSiscf .eq. 1)
     .                write(6,'(/,a11,a7,3a14,a8,a7,a11)')
     .                tmpname,' TSiscf', '   Eharris(eV)',
     .                '      E_KS(eV)', '   FreeEng(eV)',
     .                '   dDmax', '  Ef_up', '  Ef_dn(eV)'
                 write(6,'(a11,i4,3f14.4,f8.4,2f9.4)')
     .                tmpname,TSiscf, Eharrs/eV, Etot/eV, FreeE/eV,
     .                dDmax, (Efs(i)/eV,i=1,2)
              else
                 if (TSiscf .eq. 1) write(6,'(/,a11,a7,3a14,2a8)')
     .                tmpname,' TSiscf', '   Eharris(eV)',
     .                '      E_KS(eV)', '   FreeEng(eV)',
     .                '   dDmax', '  Ef(eV)'
                 write(6,'(a11,i4,3f14.4,2f8.4)')
     .                tmpname,TSiscf, Eharrs/eV, Etot/eV,
     .                FreeE/eV, dDmax, Ef/eV
              endif             !fixspin
        endif                !TSrun
! TSS End
#endif /* TRANSIESTA */
        
        if (cml_p) then
          call cmlEndPropertyList(mainXML)
        endif
        
      else !final
      ! Print out additional information in finalization.

        write(6,'(/,a)') 'siesta: Final energy (eV):'
        write(6,'(a,a15,f15.6)')
     .    'siesta: ', 'Band Struct. =', Ebs/eV,
     .    'siesta: ',      'Kinetic =', Ekin/eV,
     .    'siesta: ',      'Hartree =', Uscf/eV,
     .    'siesta: ',   'Ext. field =', DUext/eV,
     .    'siesta: ',  'Exch.-corr. =', Exc/eV,
     .    'siesta: ', 'Ion-electron =', (Enascf+Enl+DUscf-Uscf-Uatm)/eV,
     .    'siesta: ',      'Ion-ion =', (Ena+Uatm-Enaatm-Eions)/eV,
     .    'siesta: ',      'Ekinion =', Ekinion/eV,
     .    'siesta: ',        'Total =', (Etot+Ekinion)/eV
        if (cml_p) then
          call cmlStartPropertyList(xf=mainXML, title='Final Energy')
          call cmlAddProperty(xf=mainXML, value=Ebs/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ebs', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Ekin/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekin', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Uscf/eV,  
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Uscf', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=DUext/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:DUext', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Exc/eV,   
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Exc', fmt='r7')
          call cmlAddProperty(xf=mainXML, 
     .         value=(Enascf+Enl+DUscf-Uscf-Uatm)/eV, 
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-e', fmt='r7')
          call cmlAddProperty(xf=mainXML, 
     .         value=(Ena+Uatm-Enaatm-Eions)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:I-I', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=Ekinion/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Ekinion', fmt='r7')
          call cmlAddProperty(xf=mainXML, value=(Etot+Ekinion)/eV,
     .         units='siestaUnits:eV', 
     .         dictref='siesta:Etot', fmt='r7')
          call cmlEndPropertyList(mainXML)
        endif !cml_p
      endif !final

      end subroutine siesta_write_energies

      END MODULE write_subs_energies
!
!     Main module
!
      Module write_subs

      use m_steps
      use write_subs_pressure,  only : siesta_write_stress_pressure
      use write_subs_energies,  only : siesta_write_energies
      use write_subs_positions, only : siesta_write_positions

      private
      public :: siesta_write_forces, siesta_write_stress_pressure,
     &          siesta_write_energies, siesta_write_positions
      CONTAINS
      
      subroutine siesta_write_forces()
      USE siesta_options
      use siesta_geom
      use siesta_cml
      use units
      use m_forces

      implicit none

      integer :: ia, ix
      real(dp) :: cfmax, fmax, fres
      real(dp) :: ftot(3), cftot(3)

        fmax = maxval(abs(fa))
        ftot = sum(fa, dim=2)
        cfmax = maxval(abs(cfa))
        cftot = sum(cfa, dim=2)
        fres = sqrt( sum(fa**2) / (na_u*3.0_dp))

      ! Almost the same forces output whether during simulation
      ! or at the end. Unfortunately not quite, therefore slightly
      ! tortuous logic below. If we are content to change format
      ! of output file slightly, this can be simplified.
      if (.not.final) then
        ! print forces to xml every step. 
        ! output forces to stdout depending on writef
        if (cml_p) then
          call cmlStartPropertyList(mainXML, title='Forces')
          call cmlAddProperty(xf=mainXML, value=fa*Ang/eV,
     .         dictref='siesta:forces', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, value=ftot*Ang/eV,
     .         dictref='siesta:ftot', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, value=fmax*Ang/eV, 
     .         dictref='siesta:fmax', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, value=fres*Ang/eV,
     .         dictref='siesta:fres', units='siestaUnits:evpa')
          call cmlAddProperty(xf=mainXML, value=cfmax*Ang/eV, 
     .         dictref='siesta:cfmax', units='siestaUnits:evpa')
          call cmlEndPropertyList(mainXML)
        endif

        write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
        if (writef) then
          write(6,'(i6,3f12.6)')(ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
        endif
        ! Always write forces to .FA file
        call iofa( na_u, fa )

        write(6,'(40("-"),/,a6,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
        write(6,'(40("-"),/,a6, f12.6)') 'Max',fmax*Ang/eV
        write(6,'(a6,f12.6,a)')'Res',fres*Ang/eV,
     .       '    sqrt( Sum f_i^2 / 3N )'
        write(6,'(40("-"),/,a6, f12.6,a)') 'Max',cfmax*Ang/eV, 
     .       '    constrained'
      else !not final
! In finalization, only print forces if sufficiently large.
        if (fmax .gt. ftol) then
          write(6,'(/,a)') 'siesta: Atomic forces (eV/Ang):'
          write(6,'(a,i6,3f12.6)')
     .         ('siesta: ', ia,(fa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
          write(6,'(a,40("-"),/,a,a6,3f12.6)')
     .         'siesta: ','siesta: ','Tot',(ftot(ix)*Ang/eV,ix=1,3)
          if (cml_p) then
            call cmlStartPropertyList(mainXML, title='Force Summary')
            call cmlAddProperty(xf=mainXML, value=fa*Ang/eV,
     .           dictref='siesta:forces', units='siestaUnits:evpa')
            call cmlAddProperty(xf=mainXML, value=ftot*Ang/eV, 
     .           dictref='siesta:ftot', units='siestaUnits:evpa')
            call cmlEndPropertyList(mainXML)
          endif !cml_p
        endif
        if (Any(cfa /= fa)) then
          if (cfmax .gt. ftol) then
            write(6,'(/,a)') 'siesta: Constrained forces (eV/Ang):'
            write(6,'(a,i6,3f12.6)')
     .           ('siesta: ',ia,(cfa(ix,ia)*Ang/eV,ix=1,3),ia=1,na_u)
            write(6,'(a,40("-"),/,a,a4,3f12.6)')
     .           'siesta: ','siesta: ','Tot',(cftot(ix)*Ang/eV,ix=1,3)
            if (cml_p) then
              call cmlStartPropertyList(mainXML,
     .             title='Constrained Force Summary')
              call cmlAddProperty(xf=mainXML, value=cfa*Ang/eV, 
     .             dictref='siesta:cforces', units='siestaUnits:evpa')
              call cmlAddProperty(xf=mainXML, value=cftot*Ang/eV, 
     .             dictref='siesta:cftot', units='siestaUnits:evpa')
              call cmlEndPropertyList(mainXML)
            endif !cml_p
          endif
        endif
      endif !final for forces

      end subroutine siesta_write_forces

      End Module write_subs
