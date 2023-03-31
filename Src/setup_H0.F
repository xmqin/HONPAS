! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_H0
      private
      public :: setup_H0
      CONTAINS

      subroutine setup_H0(G2max)
      
C     Computes non-self-consistent part of the Hamiltonian
C     and initializes data structures on the grid.
      
      USE siesta_options,  only: g2cut
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_so_2D
      use sparse_matrices, only: Dscf

      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_dhscf,      only: dhscf_init
      use m_energies,   only: Eions, Ena, DEna, Emm, Emeta, Eso
      use m_ntm
      use m_spin,       only: spin
      use spinorbit,    only: spinorb
      use alloc, only: re_alloc, de_alloc
      use class_dSpData1D, only: val
      use class_dSpData2D, only: val

#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none
      real(dp), intent(inout) :: g2max
      
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1), dummy_dm(1,1)
      real(dp) :: dummy_E
      integer  :: ia, is

      real(dp), pointer :: H_val(:), H_so(:,:)

#ifdef DEBUG
      call write_debug( '    PRE setup_H0' )
#endif

!----------------------------------------------------------------------BEGIN
      call timer('Setup_H0',1)

C     Self-energy of isolated ions
      Eions = 0.0_dp
      do ia = 1,na_u
        is = isa(ia)
        Eions = Eions + uion(is)
      enddo

!     In these routines, add a flag to tell them NOT to compute
!     forces and stresses in this first pass, only energies.

!     Neutral-atom: energy

      call naefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &           isa, Ena, dummy_fa, dummy_stress,
     &           forces_and_stress=.false.)
      call dnaefs(na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, dummy_fa, dummy_stress,
     &            forces_and_stress=.false.) 
      Ena = Ena + DEna
 
C     Metadynamics energy
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,dummy_fa,dummy_stress,
     $           .false.,.false.)
      endif

C     Add on force field contribution to energy
      call twobody( na_u,xa,isa,ucell,Emm,
     &              ifa=0,fa=dummy_fa,istr=0,stress=dummy_stress)

!
!     Now we compute matrix elements of the Kinetic and Non-local
!     parts of H

!     Kinetic: matrix elements only
      H_val => val(H_kin_1D)
      H_val(:) = 0.0_dp
      call kinefsm(na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &             maxnh, maxnh, lasto, iphorb, isa, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             1,
     &             dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &             H_val,
     &             matrix_elements_only=.true.) 

!     Non-local-pseudop:  matrix elements only
      H_val => val(H_vkb_1D)
      H_val(:) = 0.0_dp
      call nlefsm(scell, na_u, na_s, isa, xa, indxua, 
     &            maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &            numh, listhptr, listh, numh, listhptr, listh, 
     &            1,
     &            dummy_dm, dummy_E, dummy_fa, dummy_stress,
     &            H_val,
     &            matrix_elements_only=.true.) 


! ..................

! If in the future the spin-orbit routine is able to compute
! forces and stresses, then "last" will be needed. If we are not
! computing forces and stresses, calling it in the first iteration
! should be enough
!
      if ( spin%SO ) then
         H_so => val(H_so_2D)
         H_so = 0._dp
         call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     &        maxnh,numh,listhptr,listh,Dscf,H_so,Eso)
      else
         Eso = 0._dp
      end if
      

C     This will take care of possible changes to the mesh and atomic-related
C     mesh structures for geometry changes
      g2max = g2cut
      call dhscf_init( spin%Grid, no_s, iaorb, iphorb,
     &                 no_l, no_u, na_u, na_s,
     &                 isa, xa, indxua, ucell,
     &                 mscell, G2max, ntm,
     &                 maxnh, numh, listhptr, listh, datm,
     &                 dummy_fa, dummy_stress)
  
      call timer('Setup_H0',2)

#ifdef DEBUG
      call write_debug( '    POS setup_H0' )
#endif

!---------------------------------------------------------------------- END
      END subroutine setup_H0
      END module m_setup_H0
