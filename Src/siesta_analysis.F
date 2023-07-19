      module m_siesta_analysis
      use write_subs

      private
      public :: siesta_analysis

      CONTAINS

      subroutine siesta_analysis( relaxd )
      USE band, only: nbk, bk, maxbk, bands
      USE writewave, only: nwk, wfk, wwave, maxwk
      USE writewave, only: setup_wfs_list, wfs_filename
      USE m_ksvinit, only: nkpol, kpol, wgthpol
      use m_ksv
      USE m_projected_DOS, only: projected_DOS
      USE m_local_DOS, only: local_DOS
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use m_fdf_global, only: fdf_global_get
      use atomlist, only: indxuo, iaorb, lastkb, lasto, datm, no_l,
     .                    iphkb, no_u, no_s, iza, iphorb, rmaxo, indxua
      use fdf, only : fdf_boolean
      use writewave,    only: wwave
      use siesta_cml
      use files,        only : slabel
      use zmatrix,      only: lUseZmatrix, write_zmatrix
      use Kpoint_grid
      use parallel, only: IOnode
      use files,       only : label_length
      use parse
      use m_energies
      use m_steps
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use m_eo
      use m_forces,          only: fa
      use m_gamma 
      use alloc,             only: re_alloc, de_alloc
      use basis_enthalpy,    only: write_basis_enthalpy
      use m_partial_charges, only: want_partial_charges
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none

      character(len=label_length+5), external :: paste

      logical :: relaxd, getPSI

      character(len=label_length+5) :: fildrh
      character(len=label_length+5) :: filepsch
      character(len=label_length+5) :: filetoch
      character(len=label_length+5) :: filevh
      character(len=label_length+5) :: filevt
      character(len=label_length+5) :: filrho

      real(dp), dimension(3,3)      :: dummy_str, dummy_strl
      real(dp)                      :: qspin(4)         ! Local

      real(dp)                      :: polxyz(3,nspin)  ! Automatic, small
      real(dp)                      :: polR(3,nspin)    ! Automatic, small
      real(dp)                      :: qaux             ! Auxiliary to call spnvec (not used later)
      real(dp)                      :: stot             ! Total spin polarization
      real(dp)                      :: svec(3)          ! Total spin vector in non-collinear-spin calculations
      real(dp), pointer             :: ebk(:,:,:)       ! Band energies

      integer                       :: j, ix, ind, ik, io, ispin
      integer                       :: wfs_band_min, wfs_band_max

      real(dp) :: g2max
#ifdef MPI
      real(dp):: qtmp(4)   ! Temporary for globalized spin charge
#endif

!------------------------------------------------------------------------- BEGIN
      final = .true.

      if (cml_p) then
        call cmlStartModule(xf=mainXML, title='Finalization')
      endif
!
!     NOTE that the geometry output by the following sections
!     is that "predicted" for the next MD or relaxation step.
!     This can be considered a feature, or a bug.
!
      if (IOnode) then
! Print atomic coordinates (and also unit cell for ParrRah.)
        if (nmove .ne. 0) then
          if (relaxd) 
     .      call outcoor(ucell, xa, na_u, 'Relaxed', .true. )
          if (.not.relaxd) 
     .      call outcoor(ucell, xa, na_u,
     .                 'Final (unrelaxed)', .true. )
        endif
        call siesta_write_positions()
        if (lUseZmatrix) call write_Zmatrix
        if ( varcel .or. (idyn.eq.8)) call outcell(ucell)

! Print coordinates in xmol format in a separate file

        if (fdf_boolean('WriteCoorXmol',.false.)) 
     .     call coxmol(iza, xa, na_u )

! Print coordinates in cerius format in a separate file

        if (fdf_boolean('WriteCoorCerius',.false.))
     .     call coceri(iza, xa, ucell, na_u, sname )

! Find interatomic distances (output in file BONDS_FINAL)

        call bonds( ucell, na_u, isa, xa,
     $       rmax_bonds, trim(slabel) // ".BONDS_FINAL" )

       endif ! IONode
!--- end output of geometry information
!
!
! NOTE: In the following sections, wavefunction generation, computation
!       of band structure, etc, are carried out using the last Hamiltonian
!       generated in the SCF run for the last geometry considered.

! Find and print wavefunctions at selected k-points
      if (nwk.gt.0) then
        wfs_filename = trim(slabel)//".selected.WFSX"
        if (IONode) print "(a)",
     $              "Writing WFSX for selected k-points in "
     $                           // trim(wfs_filename)
        call wwave( no_s, nspin, nspin, no_u, no_l, maxnh, maxwk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              nwk, wfk, no_u, gamma, occtol )
      endif


      if (write_coop) then
        ! Output the wavefunctions for the kpoints in the SCF set
        ! Note that we allow both a band number and an energy range
        ! The user is responsible for using appropriate values.
        call fdf_global_get(wfs_band_min,"WFS.BandMin",1)
        call fdf_global_get(wfs_band_max,"WFS.BandMax",no_u)
        call setup_wfs_list(nkpnt,no_u,wfs_band_min,wfs_band_max,
     $                      use_scf_weights=.true.,
     $                      use_energy_window=.true.)
         wfs_filename = trim(slabel)//".fullBZ.WFSX"
         if (IONode) print "(a)", "Writing WFSX for COOP/COHP in " 
     $                           // trim(wfs_filename)
         call wwave( no_s, nspin, nspin, no_u, no_l, maxnh, maxk,
     .        numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .        nkpnt, kpoint, no_u, gamma, occtol)
      endif

      ! Compute bands

      nullify( ebk )
      call re_alloc( ebk, 1, no_u, 1, nspin, 1, maxbk,
     &               name    = 'ebk',
     &               routine = 'siesta_analysis' )

      if (nbk.gt.0) then
        if (IONode) print "(a)", "Computing bands..."
        call fdf_global_get(getPSI,'WFS.Write.For.Bands', .false.)
        if (getPSI) then
           wfs_filename = trim(slabel)//".bands.WFSX"
         if (IONode) print "(a)", "Writing WFSX for bands in "
     $                           // trim(wfs_filename)
           call fdf_global_get(wfs_band_min,"WFS.BandMin",1)
           call fdf_global_get(wfs_band_max,"WFS.BandMax",no_u)
           call setup_wfs_list(nbk,no_u,wfs_band_min,wfs_band_max,
     $          use_scf_weights=.false.,
     $          use_energy_window=.false.)
        endif
        call bands( no_s, nspin, nspin, no_u, no_l, maxnh, maxbk,
     .              numh, listhptr, listh, H, S, Ef, xijo, indxuo,
     .              .true., nbk, bk, ebk, occtol, getPSI )
        if (IOnode) then
          if ( writbk ) then
            write(6,'(/,a,/,a4,a12)')
     .       'siesta: Band k vectors (Bohr**-1):', 'ik', 'k'
            do ik = 1,nbk
              write(6,'(i4,3f12.6)') ik, (bk(ix,ik),ix=1,3)
            enddo
          endif
        
          if ( writb ) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Band energies (eV):', 'ik', 'is', 'eps'
            do ispin = 1,min(nspin,2)
              do ik = 1,nbk
                write(6,'(i4,i3,10f7.2)')
     .            ik, ispin, (ebk(io,ispin,ik)/eV,io=1,min(10,no_u))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .              (ebk(io,ispin,ik)/eV,io=11,no_u)
              enddo
            enddo
          endif
        endif
      endif

! Print eigenvalues
      if (IOnode .and. writeig) then
        if (isolve.eq.0 .and. no_l.lt.1000) then
          if (nspin .le. 2) then
            write(6,'(/,a,/,a4,a3,a7)')
     .       'siesta: Eigenvalues (eV):', 'ik', 'is', 'eps'
            do ik = 1,nkpnt
              do ispin = 1,nspin
                write(6,'(i4,i3,10f7.2)')
     .            ik,ispin,(eo(io,ispin,ik)/eV,io=1,min(10,neigwanted))
                if (no_u.gt.10) write(6,'(7x,10f7.2)')
     .            (eo(io,ispin,ik)/eV,io=11,neigwanted)
              enddo
            enddo
          else
            write(6,'(/,a)') 'siesta: Eigenvalues (eV):'
            do ik = 1,nkpnt
              write(6,'(a,i6)') 'ik =', ik
              write(6,'(10f7.2)')
     .          ((eo(io,ispin,ik)/eV,io=1,neigwanted),ispin=1,2)
            enddo
          endif
          write(6,'(a,f15.6,a)') 'siesta: Fermi energy =', ef/eV, ' eV'
        endif
      endif

      if (isolve.eq.0.and.IOnode)
     .     call ioeig(eo,ef,neigwanted,nspin,nkpnt,no_u,nspin,maxk,
     .                kpoint, kweight)


      call projected_DOS()

! Print program's energy decomposition and final forces
      if (IOnode) then
        call siesta_write_energies( first=.false., iscf=0, dDmax=0.0_dp)
        call siesta_write_forces()
        call siesta_write_stress_pressure()
        call write_basis_enthalpy(FreeE,FreeEHarris)
      endif

! NOTE: Here, the spin polarization is computed using Dscf, which is
!       a density matrix obtained after mixing the "in" and "out"
!       DMs of the SCF run for the last geometry considered.
!       This can be considered a feature or a bug.

! Print spin polarization
      if (nspin .ge. 2) then
        do ispin = 1,nspin
          qspin(ispin) = 0.0_dp
          do io = 1,no_l
            do j = 1,numh(io)
              ind = listhptr(io)+j
              qspin(ispin) = qspin(ispin) + Dscf(ind,ispin) * S(ind)
            enddo
          enddo
        enddo

#ifdef MPI
! Global reduction of spin components
        call globalize_sum(qspin(1:nspin),qtmp(1:nspin))
        qspin(1:nspin) = qtmp(1:nspin)
#endif
        if (nspin .eq. 2) then
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', 
     .       qspin(1) - qspin(2)
          endif
          if (cml_p) call cmlAddProperty(xf=mainXML,
     .         value=qspin(1)-qspin(2), dictref='siesta:qspin',
     .         units='siestaUnits:spin')
        elseif (nspin .eq. 4) then
          call spnvec( nspin, qspin, qaux, stot, svec )
          if (IOnode) then
            write(6,'(/,a,f12.6)')
     .       'siesta: Total spin polarization (Qup-Qdown) =', stot
            write(6,'(a,3f12.6)') 'siesta: Spin vector =', svec
            if (cml_p) then
              call cmlAddProperty(xf=mainXML, value=stot,
     .             dictref='siesta:stot',
     .             units='siestaUnits:spin')
              call cmlAddProperty(xf=mainXML, value=svec,
     .             dictref='siesta:svec',
     .             units='siestaUnits:spin')
            endif !cml_p
          endif
        endif
      endif

! Print electric dipole
      if (shape .ne. 'bulk') then
        if (IOnode) then
          write(6,'(/,a,3f12.6)')
     .      'siesta: Electric dipole (a.u.)  =', dipol
          write(6,'(a,3f12.6)')
     .      'siesta: Electric dipole (Debye) =', 
     .      (dipol(ix)/Debye,ix=1,3)
        endif
        if (cml_p) then
          call cmlAddProperty(xf=mainXML, value=dipol/Debye,
     .         title='Electric dipole', dictref='siesta:dipol',
     .         units='siestaUnits:Debye')
        endif !cml_p
      endif

! NOTE: The use of *_last geometries in the following sections
!       guarantees that the analysis of the electronic structure
!       is done for the geometry for which it was computed.

! Calculation of the bulk polarization using the Berry phase
! formulas by King-Smith and Vanderbilt
! Attention H0 is used as an auxiliary array
      if (nkpol.gt.0 .and. .not.bornz) then
        call KSV_pol(na_u, na_s, xa_last, rmaxo, scell_last, ucell_last,
     .               no_u, no_l, no_s, nspin, qspin, 
     .               maxnh, nkpol, numh, listhptr, listh, 
     .               H, S, H0, xijo, indxuo, isa, iphorb, 
     .               iaorb, lasto, shape,
     .               nkpol,kpol,wgthpol, polR, polxyz ) 
      endif

! Calculation of the optical conductivity
! Attention H0, Eold, Dold are used as auxiliary arrays
      call optical(na_u, na_s, xa_last, scell_last, ucell_last,
     .             no_u, no_l, no_s, nspin, qspin,
     .             maxnh, numh, listhptr, listh, H, S, H0,
     .             Eold(1,1), Dold(1,1),
     .             xijo, indxuo, indxua, ebk, ef, temp,
     .             isa, iphorb, iphKB, iaorb, lasto, lastkb,
     .             shape )

!...................................

!
!  NOTE: Dscf here is the mixed one (see above).
!
      want_partial_charges = (hirshpop .or. voropop) .AND.
     $                   (.not. partial_charges_at_every_geometry)
!
! Save electron density and potential
      if (savrho .or. savdrh .or. savevh .or. savevt .or.
     .    savepsch .or. savetoch
     $     .or. want_partial_charges) then
        filrho = ' '
        fildrh = ' '
        filevh = ' '
        filevt = ' '
        filepsch = ' '
        filetoch = ' '
        if (savrho) filrho = paste( slabel, '.RHO' )
        if (savdrh) fildrh = paste( slabel, '.DRHO' )
        if (savevh) filevh = paste( slabel, '.VH'  )
        if (savevt) filevt = paste( slabel, '.VT'  )
        if (savepsch) filepsch = paste( slabel, '.IOCH'  )
        if (savetoch) filetoch = paste( slabel, '.TOCH'  )
        g2max = g2cut
        call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa_last, indxua, 
     .              ucell_last, mscell, g2max, ntm,
     .              0, 0, 0, filrho, fildrh, filevh, filevt,
     .              filepsch, filetoch, 
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dipol, fa, dummy_str, fa, dummy_strl )
                    ! next to last argument is dummy here,
                    ! as no forces are calculated
                    ! todo: make all these optional
      endif

! Find local density of states
      call local_DOS()

      call de_alloc( ebk, name="ebk")

!--------------------------------------------------------------------------- END
      END subroutine siesta_analysis
      END module m_siesta_analysis
