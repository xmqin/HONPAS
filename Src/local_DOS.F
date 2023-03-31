! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      MODULE m_local_DOS
      private
      public :: local_DOS

      CONTAINS

      subroutine local_DOS( )
      use m_energies

      use sparse_matrices
      USE siesta_options
      use siesta_geom
      use atomlist,       only: indxuo, indxua           
      use atomlist,       only: qtot, qtots, no_u, no_l
      use atomlist,       only: iphorb                   
      use atomlist,       only: datm, no_s, iaorb        
      use fdf
      use sys,            only: die                      
      use files,          only: slabel     ! system label
      use files,          only: filesOut_t ! derived type for output file names
      use Kpoint_grid
      use parallel,       only: IOnode                   
      use files,          only : label_length            
      use m_ntm
      use m_forces,       only: fa
      use m_energies,     only: Ef, Efs
      use m_eo
      use m_spin,         only: nspin
      use m_spin,         only: spinor_dim 
      use m_diagon,       only: diagon
      use m_dhscf,        only: dhscf
      use m_partial_charges, only: want_partial_charges
      implicit none

      integer :: dummy_iscf = 1

      real(dp):: e1  ! Lower bound of energy range
      real(dp):: e2  ! Upper bound of energy range

      real(dp)  :: dummy_str(3,3), dummy_strl(3,3)  ! for dhscf call
      real(dp)  :: dummy_dipol(3)

      real(dp)  :: factor, g2max, dummy_Entrop
      logical :: genlogic ! Generate local density of states?
      logical :: want_partial_charges_save ! to skip Voronoi/Hirshfeld in LDOS calculation

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline
      type(filesOut_t)           :: filesOut  ! blank output file names

#ifdef DEBUG
      call write_debug( '  PRE local_DOS' )
#endif
! Find local density of states
      genlogic = fdf_block('LocalDensityOfStates',bfdf)

      if ( genlogic ) then

        want_partial_charges_save = want_partial_charges
        want_partial_charges = .false.

! Find the desired energy range
        if (.not. fdf_bline(bfdf,pline))
     .    call die('local_DOS: ERROR in LocalDensityOfStates block')
        if (.not. fdf_bmatch(pline,'vvn'))
     .    call die("Wrong format in LocalDensityofStates")
        factor = fdf_convfac( fdf_bnames(pline,1), 'Ry' )
        e1 = fdf_bvalues(pline,1)*factor
        e2 = fdf_bvalues(pline,2)*factor
        call fdf_bclose(bfdf)

        !Find the density matrix for states between e1 and e2
        if ((isolve .eq. SOLVE_DIAGON) .or.
     .      ((isolve .eq. SOLVE_MINIM) .and.
     .       minim_calc_eigenvalues)) then
          call diagon(no_s, spinor_dim, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qtots, temp, e1, e2,
     .                xijo, indxuo, gamma_SCF, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, dummy_Entrop, no_u,
     .                occtol, dummy_iscf, neigwanted)

          !Find the LDOS in the real space mesh
          filesOut%rho = trim(slabel) // '.LDOS'
          g2max = g2cut
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa_last, indxua, 
     .              ntm, 0, 0, 0, filesOut,
     .              maxnh, numh, listhptr, listh, Dscf, Datm, maxnh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dummy_dipol, dummy_str, fa, dummy_strl )
                    ! next to last argument is dummy here,
                    ! as no forces are calculated
                    ! todo: make all these optional
        else
          if (IOnode)  write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
        endif

        want_partial_charges = want_partial_charges_save

      endif ! genlogic


#ifdef DEBUG
      call write_debug( '  POS local_DOS' )
#endif
      END subroutine local_DOS

      END module m_local_DOS
