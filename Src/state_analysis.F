      MODULE m_state_analysis
      use write_subs

      private
      public :: state_analysis

      CONTAINS

      subroutine state_analysis( istep  )
      use siesta_cml
      use m_born_charge, only: born_charge
      use parallel, only: IOnode
      use m_wallclock, only : wallclock
      use zmatrix,      only: lUseZmatrix, iofaZmat, write_Zmatrix,
     .                        CartesianForce_to_ZmatForce,
     .                        write_canonical_ucell_and_Zmatrix
      use m_iostruct,   only: write_struct
      use sparse_matrices
      use atomlist, only: iaorb, iphorb, iza, amass, no_u, lasto, superx
      use siesta_geom
      USE siesta_options
      use m_stress
      use m_energies, only: Etot, FreeE, Eharrs, FreeEHarris, Entropy
      use m_ntm
      use m_spin,     only: nspin
      use m_forces
      use m_intramol_pressure, only: remove_intramol_pressure

      implicit none

      integer  :: istep
      integer  :: ia, jx, ix
      real(dp) :: volume
      logical  :: eggbox_block=.true. ! Read eggbox info from data file?

      external :: eggbox, mulliken, fixed
      real(dp), external :: volcel

!------------------------------------------------------------------------- BEGIN

!     Save structural information in crystallographic format,
!     and canonical Zmatrix (if applicable) before moving the atoms
!

      if (cml_p) then
         call cmlStartModule(xf=mainXML, title='SCF Finalization')   
      endif

      if (IOnode) then
         call write_struct( ucell, na_u, isa, iza, xa )
         if (lUseZmatrix .and. (idyn.eq.0)) then
            call write_canonical_ucell_and_Zmatrix()
         endif
      endif
!

! Write final Kohn-Sham Energy 
      if (cml_p) call cmlStartPropertyList(mainXML,
     .                                   title='Final KS Energy')
      if (IOnode) then
        if ( .not. harrisfun) 
     .    write(6,"(/a,f14.4)")  'siesta: E_KS(eV) =        ', Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, value=Etot/eV,
     .       dictref='siesta:E_KS', units='siestaUnits:eV', 
     .       fmt='r6')
      endif

! Substract egg box effect form energy 
      if (eggbox_block) then
        call eggbox('energy',ucell,na_u,isa,ntm,xa,fa,Etot,
     .               eggbox_block)
        if (IOnode)
     .    write(6,"(/a,f14.4)") 'siesta: E_KS - E_eggbox = ',Etot/eV
        if (cml_p) call cmlAddProperty(xf=mainXML, value=Etot/eV,
     .         dictref='siesta:E_KS_egg', units='siestaUnits:eV', 
     .         fmt='r6')
      endif
      if (cml_p) call cmlEndPropertyList(mainXML)

      FreeE = Etot - Temp * Entropy
      FreeEHarris = Eharrs - Temp * Entropy

! Substract egg box effect from the forces 
      if (eggbox_block) then
        call eggbox('forces',ucell,na_u,isa,ntm,xa,fa,Etot,eggbox_block)
      endif

      if (IOnode) call write_raw_efs(stress,na_u,fa,FreeE)

      !  Compute stress without internal molecular pressure
      call remove_intramol_pressure(ucell,stress,na_u,xa,fa,mstress)

C Impose constraints to atomic movements by changing forces ...........
      if (RemoveIntraMolecularPressure) then
         ! Consider intramolecular pressure-removal as another
         ! kind of constraint
         call fixed(ucell,mstress,na_u,isa, amass, xa, fa,
     $        cstress, cfa, ntcon )
      else
         call fixed(ucell,stress,na_u,isa, amass, xa, fa,
     $        cstress, cfa, ntcon )
      endif
C ...................

! Calculate and output Zmatrix forces
      if (lUseZmatrix .and. (idyn.eq.0)) then
        call CartesianForce_to_ZmatForce(na_u,xa,fa)
        if (IOnode) call iofaZmat()
      endif

      !  Compute kinetic contribution to stress
      kin_stress(1:3,1:3) = 0.0_dp
      volume = volcel(ucell)
      do ia = 1,na_u
        do jx = 1,3
          do ix = 1,3
            kin_stress(ix,jx) = kin_stress(ix,jx) -
     .             amu * amass(ia) * va(ix,ia) * va(jx,ia) / volume
          enddo
        enddo
      enddo
      !  Add kinetic term to stress tensor
      tstress = stress + kin_stress

! Force output 
      if (IOnode) then
        call siesta_write_forces()
        call siesta_write_stress_pressure()
        call wallclock('--- end of geometry step')
      endif

! Mulliken population analysis 
      call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .               numh, listhptr, listh, S, Dscf, isa, 
     .               lasto, iaorb, iphorb )
!

!     Call the born effective charge routine only in those steps (even) 
!     in which the dx  is positive.
      if (bornz .and. (mod(istep,2) .eq. 0)) then
         call born_charge()
      endif

      !End the xml module corresponding to the analysis
      if (cml_p) then
         call cmlEndModule(mainXML)         
      endif 
!--------------------------------------------------------------------------- END
      END subroutine state_analysis

      END MODULE m_state_analysis
