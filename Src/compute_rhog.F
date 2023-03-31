! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_compute_rhog
      private
      public :: compute_rhog
      CONTAINS

      subroutine compute_rhog(DM)

!      USE siesta_options
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H
      use siesta_geom
      use siesta_options, only: g2cut
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l

      use m_dhscf,      only: dhscf
      use m_stress
      use m_energies
      use m_ntm
      use m_spin,         only: nspin, h_spin_dim
      use m_dipol
      use alloc, only: re_alloc, de_alloc
      use files, only : filesOut_t    ! derived type for output file names

      implicit none

      real(dp), intent(in):: DM(maxnh,h_spin_dim)

      real(dp)            :: stressl(3,3)
      real(dp), pointer   :: fal(:,:)   ! Local-node part of atomic F
      integer             :: io, is, ispin
      integer             :: ifa     ! Calc. forces?      0=>no, 1=>yes
      integer             :: istr    ! Calc. stress?      0=>no, 1=>yes
      integer             :: ihmat   ! Calc. hamiltonian? 0=>no, 1=>yes
      real(dp)            :: g2max
      type(filesOut_t)    :: filesOut  ! blank output file names

!------------------------------------------------------------------------- BEGIN

      call timer('compute_rhog',1)

! Add SCF contribution to energy and matrix elements ..................
      g2max = g2cut

      nullify(fal)
      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

      ifa  = 0
      istr = 0
      ihmat = 0

      ! Will compute rhog if mix_charge is .true., and
      ! then return
      call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, filesOut,
     .            maxnh, numh, listhptr, listh, DM, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl,
     .            charge_density_only=.true.)

      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 

      call timer('compute_rhog',2)


!------------------------------------------------------------------------- END
      END subroutine compute_rhog
      END module m_compute_rhog
