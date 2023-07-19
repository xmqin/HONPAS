      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( first, last, iscf )
      USE m_normalize_dm, only: normalize_dm
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_overfsm,    only: overfsm
      use m_stress
      use m_energies
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: nspin
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use m_gamma
      use m_hsx, only: write_hsx
      use sys, only: die
      use m_partial_charges, only: want_partial_charges
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
#ifdef TRANSIESTA
      use files, only: slabel, label_length
      use m_ts_options, only: savetshs
      use m_ts_global_vars, only: TSrun, ts_istep
      use m_ts_kpoints, only: ts_gamma_scf, ts_kscell, ts_kdispl
      use m_ts_io, only: ts_iohs
#endif /* TRANSIESTA */


      implicit none

      logical, intent(in) :: first, last
      integer, intent(in) :: iscf

      real(dp)              :: stressl(3,3)
      real(dp), pointer     :: fal(:,:)   ! Local-node part of atomic forces

#ifdef MPI
      real(dp)              :: stresstmp(3,3)
      real(dp), pointer     :: fatmp(:,:)
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif

      integer :: ia, io, is, ispin

      integer:: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer:: istr    ! Calculate stress? 0=>no, 1=>yes
      integer:: ihmat   ! Calculate hamiltonian? 0=>no, 1=>yes

      real(dp) :: g2max
#ifdef TRANSIESTA
! TSS Begin
      character(len=label_length+13) :: fname
      character(8) :: numstring
      character(8) :: iotask
! TSS End
#endif /* TRANSIESTA */
#ifdef TRANSIESTA
! TSS Begin
! SIESTA External
      character(len=label_length+5) PASTE
! TSS End
#endif /* TRANSIESTA */



!------------------------------------------------------------------------- BEGIN

#ifndef TRANSIESTA      
      ! Normalize density matrix to exact charge, except if
      ! we are just using the DM to generate the initial rho and exit
      !
      if (.not. save_initial_charge_density) then
         call normalize_dm(first)
      endif
#else
      if (.not.TSrun ) call normalize_dm(first)
#endif /* TRANSIESTA */

      ! Initialize Hamiltonian ........................................
      H = 0.0_dp

! Initialize forces and stress ...................

      nullify(fal)
      call re_alloc(fal,1,3,1,na_u,name="fal",
     $              routine="setup_hamiltonian")

        if (first.or.last) then
          fa(1:3,1:na_u) = 0.0_dp
          stress(1:3,1:3) = 0.0_dp
          stressl(1:3,1:3) = 0.0_dp  ! fal already initialied
        endif
! ..................

! Self-energy of isolated ions ........................................
        if (first) then
          Eions = 0.0_dp
          do ia = 1,na_u
            is = isa(ia)
            Eions = Eions + uion(is)
          enddo
        endif
! ..................

! Neutral-atom: energy, forces and stress ............................
! First time for energy, last time for forces
        if (first.or.last) then
          call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .               isa, Ena, fa, stress)
          call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     .                isa, DEna, fa, stress) 
          Ena = Ena + DEna
        endif
! ..................

! Kinetic: energy, forces, stress and matrix elements .................
        if (first.or.last) then
          call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo, no_l,
     .                 maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Dscf, Ekin, fal, stressl, H ) 
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Ekin,buffer1)
          Ekin = buffer1
#endif
        endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
        if (first.or.last) then
          call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     .                maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                min(nspin,2), Dscf, Enl, fal, stressl, H)
        
#ifdef MPI
! Global reduction of energy terms
          call globalize_sum(Enl,buffer1)
          Enl = buffer1
#endif
        endif
! ..................

! Save or get partial Hamiltonian (non-SCF part) ......................
        if (first.or.last) then
          do io = 1,nh
            H0(io) = H(io,1)
          enddo
        else
          do ispin = 1,nspin
            if (ispin .le. 2) then
              do io = 1,nh
                H(io,ispin) = H0(io)
              enddo
            else
              do io = 1,nh
                H(io,ispin) = 0.0_dp
              enddo
            endif
          enddo          
        endif
! ..................

! Non-SCF part of total energy .......................................
        if (first.or.last) then
          E0 = -Eions + Ena + Ekin + Enl
        else
          E0 = 0.0_dp
          do ispin = 1,min(nspin,2)
            do io = 1,nh
              E0 = E0 + H0(io) * Dscf(io,ispin)
            enddo
          enddo
#ifdef MPI
! Global reduction of E0
          call globalize_sum(E0,buffer1)
          E0 = buffer1
#endif
          E0 = E0 - Eions + Ena
        endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
! Add SCF contribution to energy and matrix elements ..................
        g2max = g2cut
        if (last) then
! Last call to dhscf and grid-cell sampling if requested
          ifa  = 1
          istr = 1
          call grdsam( nspin, no_s, iaorb, iphorb, 
     .                 no_l, no_u, na_u, na_s, isa, xa, indxua,
     .                 ucell, mscell, g2max, ntm, ifa, istr, maxnh,
     .                 maxnh, numh, listhptr, listh, Dscf, Datm, H,
     .                 Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                 Exc, Dxc, dipol, fa, stress, fal, stressl)
        else
          ifa  = 0
          istr = 0
          ihmat = 1
          if ((hirshpop .or. voropop)
     $         .and. partial_charges_at_every_scf_step) then
             want_partial_charges = .true.
          endif
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .                no_u, na_u, na_s, isa, xa, indxua, 
     .                ucell, mscell, g2max, ntm,
     .                ifa, istr, ihmat, ' ', ' ', ' ', ' ', ' ', ' ',
     .                maxnh, numh, listhptr, listh, Dscf, Datm,
     .                maxnh, numh, listhptr, listh, H,
     .                Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .                Exc, Dxc, dipol, fa, stress, fal, stressl)
          want_partial_charges = .false.
        endif
            
! Output memory use after first call to dhscf
        if (istp.eq.1 .and. iscf.eq.1) call printmemory( 6, 0 )

*       if (istp.eq.1 .and. iscf.eq.1) write(6,'(/,a,f10.3,a)')
*    .    'siesta: dhscf mesh cutoff =', g2max, ' Ry'

! ..................

! Orthonormalization forces ...........................................
        if (last) then
           call overfsm(na_u, na_s, no_s, scell, xa, indxua,rmaxo,no_l,
     .                 maxnh, maxnh, lasto, iphorb, isa, 
     .                 numh, listhptr, listh, numh, listhptr, listh, 
     .                 min(nspin,2), Escf, fal, stressl, S ) 
        endif
! ..................

        ! Metadynamics forces
        if (lMetaForce.and.(first.or.last)) then
          call meta(xa,na_u,ucell,Emeta,fa,stress,last,last)
        endif

        ! Add on force field contribution if required
        if (first .or. last) then
           call twobody(na_u,xa,isa,ucell,Emm,
     $                  ifa=1,fa=fa,istr=1,stress=stress)
        endif

        if (first .or. last) then    ! Could it really be "last" only?
#ifdef MPI
           !Global reduction of forces and stresses
           nullify(fatmp)
           call re_alloc(fatmp,1,3,1,na_u,name="fatmp",
     $              routine="setup_hamiltonian")
           call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
           call globalize_sum(fal(1:3,1:na_u),fatmp(1:3,1:na_u))
           stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
           fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
           call de_alloc(fatmp,routine="setup_hamiltonian")
#else
           stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
           fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif
        endif

        call de_alloc(fal,routine="setup_hamiltonian")

! Find entropy ........................................................
! Entropy is returned from the call to diagon. To add to the energy
! the entropy computed from the input charge, here it is assigned to the one
! of the former SCF step

        Entropy = 0.0_dp
        if (isolve .eq. 0) then
          if (istp.gt.1 .or. iscf.gt.1) then
            Entropy = Entrop
          endif
        endif

! Save present density matrix ........................................
        do is = 1,nspin
          do io = 1,nh
            Dold(io,is) = Dscf(io,is)
            Eold(io,is) = Escf(io,is)
          enddo
        enddo

! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file

        if (write_coop .or. savehs) then
          call write_hsx( gamma, no_u, no_s, nspin, indxuo,
     $               maxnh, numh, listhptr, listh, H, S, qtot, temp,
     $               xijo)
        endif

#ifdef TRANSIESTA
        if (last.and.savetshs) then
           if (idyn.eq.6) then ! Force constants run
              write(numstring,'(2i4)') ia1,ts_istep
              fname = paste(slabel, numstring)
              fname = paste(fname, '.TSHS')
              iotask='write'
              call ts_iohs(iotask, gamma, .false., no_u, no_s, nspin,
     .             indxuo, maxnh, numh, listhptr, listh, H, S, qtot, 
     .             temp, xijo,
     .             label_length+5, fname, na_u, lasto, isa, ef, ucell,
     .             ts_kscell, ts_kdispl, ts_gamma_scf, xa, 
     .             ts_istep, ia1)
           else
              fname = paste( slabel, '.TSHS' )
              iotask='write'
              call ts_iohs(iotask, gamma, .false., no_u, no_s, nspin,
     .             indxuo, maxnh, numh, listhptr, listh, H, S, qtot, 
     .             temp, xijo,
     .             label_length+5, fname, na_u, lasto, isa, ef, ucell,
     .             ts_kscell, ts_kdispl, ts_gamma_scf, xa, 0, 0)
           endif
        endif
#endif /* TRANSIESTA */


!--------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
