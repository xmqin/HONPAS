! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_setup_hamiltonian
      private
      public :: setup_hamiltonian
      CONTAINS

      subroutine setup_hamiltonian( iscf )

      USE siesta_options
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: H_dftu_2D, H_so_2D
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H, S, Hold
      use sparse_matrices, only: Dscf, Escf, xijo
      use class_dSpData1D,  only: val
      use class_dSpData2D,  only: val

      use siesta_geom
      use atmfuncs, only: uion
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     .                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     .                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use dftu_specs,   only: switch_dftu     ! This variable determines whether
                                              !   the subroutine to compute the
                                              !   Hubbard terms should be called
                                              !   or not
      use m_dftu,       only: hubbard_term    ! Subroutine that compute the
                                              !   Hubbard terms
      use m_dhscf,      only: dhscf
      use m_stress
      use m_energies
      use parallel, only: Node
      use m_steps, only: istp
      use m_ntm
      use m_spin,         only: spin
      use m_dipol
      use alloc, only: re_alloc, de_alloc
      use m_hsx, only: write_hsx
      use sys, only: die, bye
      use m_partial_charges, only: want_partial_charges
      use files, only : filesOut_t    ! derived type for output file names
      use m_rhog, only: rhog_in, rhog
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
#ifdef HAVE_LIBINT
      use nao2gto_data, only: hfx_libint, hfx_options, hfx_system
      use nao2gto_hfx,  only: setup_hfx
#endif

      implicit none
      integer, intent(in) :: iscf
      real(dp)            :: stressl(3,3)
      real(dp), pointer   :: fal(:,:)   ! Local-node part of atomic F
#ifdef MPI
      real(dp)            :: buffer1
#endif
      integer             :: io, is, ispin
      integer             :: ifa     ! Calc. forces?      0=>no, 1=>yes
      integer             :: istr    ! Calc. stress?      0=>no, 1=>yes
      integer             :: ihmat   ! Calc. hamiltonian? 0=>no, 1=>yes
      real(dp)            :: g2max
      type(filesOut_t)    :: filesOut  ! blank output file names
      logical             :: use_rhog_in

      real(dp), pointer   :: H_vkb(:), H_kin(:), H_dftu(:,:), H_so(:,:)

!------------------------------------------------------------------------- BEGIN

      call timer('setup_H',1)

      ! Nullify pointers
      nullify(fal)

!$OMP parallel default(shared), private(ispin,io)

!     Save present H matrix
      do ispin = 1, spin%H
!$OMP do
         do io = 1,maxnh
            Hold(io,ispin) = H(io,ispin)
         enddo
!$OMP end do nowait
      enddo

!$OMP single
      H_kin => val(H_kin_1D)
      H_vkb => val(H_vkb_1D)
      if ( spin%SO ) then
         ! Sadly some compilers (g95), does
         ! not allow bounds for pointer assignments :(
         H_so => val(H_so_2D)
      end if
!$OMP end single ! keep wait

      ! We do not need to set the non-spinor components
      ! For non-colinear they are set down below,
      ! while for spin-orbit they are set to the H_so initial
      ! spin-orbit.
      
      do ispin = 1, spin%spinor
!$OMP do
          do io = 1,maxnh
            H(io,ispin) = H_kin(io) + H_vkb(io)
          end do
!$OMP end do nowait
      end do

      if ( spin%NCol ) then

         do ispin = 3 , spin%H
!$OMP do
            do io = 1, maxnh
               H(io,ispin) = 0._dp
            end do
!$OMP end do nowait
         end do
         
      else if ( spin%SO ) then

         do ispin = 3 , spin%H
!$OMP do
            do io = 1, maxnh
               H(io,ispin) = H_so(io,ispin-2)
            end do
!$OMP end do nowait
         end do

      end if
         
! ..................

! Non-SCF part of total energy .......................................
! Note that these will be "impure" for a mixed Dscf

! If mixing the charge, Dscf is the previous step's DM_out. Since
! the "scf" components of the energy are computed with the (mixed)
! charge, this introduces an inconsistency. In this case the energies
! coming out of this routine need to be corrected.
! 
!$OMP single
      Ekin = 0.0_dp
      Enl = 0.0_dp
!$OMP end single ! keep wait
      
!$OMP do reduction(+:Ekin,Enl)
      do io = 1,maxnh
        do ispin = 1, spin%spinor
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
        end do
      end do
!$OMP end do nowait

!$OMP single
      Eso = 0._dp
!$OMP end single
      if ( spin%SO ) then
!$OMP do reduction(+:Eso)
         do io = 1, maxnh
            Eso = Eso + H_so(io,1)*Dscf(io,7) + H_so(io,2)*Dscf(io,8)
     .           + H_so(io,5)*Dscf(io,3) + H_so(io,6)*Dscf(io,4)
     .           - H_so(io,3)*Dscf(io,5) - H_so(io,4)*Dscf(io,6)
         end do
!$OMP end do nowait
      end if

!$OMP end parallel
      
#ifdef MPI
      ! Global reduction of Ekin, Enl
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
      if ( spin%SO ) then
         ! Global reduction of Eso 
         call globalize_sum(Eso,buffer1)
         Eso = buffer1
      end if
#endif

!     Non-SCF part of total energy
      call update_E0()

! Hubbard term for LDA+U: energy, forces, stress and matrix elements ....
      if( switch_dftu ) then
        if ( spin%NCol ) then
         call die('LDA+U cannot be used with non-colinear spin.')
        end if
        if ( spin%SO ) then
         call die('LDA+U cannot be used with spin-orbit coupling.')
        end if

        call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

        H_dftu => val(H_dftu_2D)
        call hubbard_term(scell, na_u, na_s, isa, xa, indxua,
     .                    maxnh, maxnh, lasto, iphorb, no_u, no_l, 
     .                    numh, listhptr, listh, numh, listhptr, listh,
     .                    spin%spinor, Dscf, Edftu, DEdftu, H_dftu, 
     .                    fal, stressl, H, iscf,
     .                    matrix_elements_only=.true.)

#ifdef MPI
        ! Global reduction of energy terms
        call globalize_sum(Edftu,buffer1)
        Edftu = buffer1
        ! DEdftu should not be globalized
        ! as it is based on globalized occupations
#endif
        Edftu = Edftu + DEdftu

        call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 
      endif
! ..................


! Add SCF contribution to energy and matrix elements ..................
      g2max = g2cut

      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'setup_hamiltonian' )

      ifa  = 0
      istr = 0
      ihmat = 1
      if ((hirshpop .or. voropop)
     $     .and. partial_charges_at_every_scf_step) then
         want_partial_charges = .true.
      endif
      use_rhog_in =  (mix_charge .and. iscf > 1)
  
      call dhscf( spin%Grid, no_s, iaorb, iphorb, no_l,
     .            no_u, na_u, na_s, isa, xa, indxua, 
     .            ntm, ifa, istr, ihmat, filesOut,
     .            maxnh, numh, listhptr, listh, Dscf, Datm,
     .            maxnh, H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     .            Exc, Dxc, dipol, stress, fal, stressl,
     .            use_rhog_in)
#ifdef HAVE_LIBINT
      ! Add Hartree-Fock exchange potential to the Hamiltonian
      if ( hfx_wanted ) then
        call setup_hfx(hfx_libint, hfx_options, hfx_system, H)
      endif
#endif
      ! This statement will apply to iscf = 1, for example, when
      ! we do not use rhog_in. Rhog here is always the charge used to
      ! build H, that is, rhog_in.
      if (mix_charge) rhog_in = rhog

      want_partial_charges = .false.
      call de_alloc( fal, 'fal', 'setup_hamiltonian' ) 
!  It is wasteful to write over and over H and S, as there are
!  no different files.
! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file

      if (savehs .or. write_coop) then
         call write_hsx( no_u, no_s, spin%H, indxuo,
     &        maxnh, numh, listhptr, listh, H, S, qtot,
     &        temp, xijo)
      endif

      call timer('setup_H',2)
#ifdef SIESTA__PEXSI
      if (node==0) call memory_snapshot("after setup_H")
#endif

      if ( h_setup_only ) then
         call timer( 'all', 2 )  ! New call to close the tree
         call timer( 'all', 3 )
         call bye("H-Setup-Only requested")
         STOP
      endif

!------------------------------------------------------------------------- END
      END subroutine setup_hamiltonian
      END module m_setup_hamiltonian
