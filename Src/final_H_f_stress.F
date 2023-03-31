! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_final_H_f_stress
      private
      public :: final_H_f_stress
      CONTAINS

      subroutine final_H_f_stress( istep, iscf , SCFconverged )
      use units, only: eV

      use files,        only : slabel
#ifdef NCDF_4
      use siesta_options, only: write_cdf
#endif
      use siesta_options, only: fixspin
      use siesta_options, only: g2cut, savehs, temp, nmove
      use siesta_options, only: recompute_H_after_scf
      use sparse_matrices, only: numh, listh, listhptr
      use sparse_matrices, only: H, S, Dscf, Escf, maxnh, xijo
      use sparse_matrices, only: H_dftu_2D
      use class_dSpData1D, only: val
      use class_dSpData2D, only: val

      use siesta_geom
      use atomlist, only: no_u, iaorb, iphkb, qtot, indxuo, datm, 
     &                    lastkb, no_s, rmaxv, indxua, iphorb, lasto,
     &                    rmaxo, no_l
      use metaforce, only: lMetaForce, meta
      use molecularmechanics, only : twobody
      use m_nlefsm,     only: nlefsm
      use m_overfsm,    only: overfsm
      use m_kinefsm,    only: kinefsm
      use m_naefs,      only: naefs
      use m_dnaefs,     only: dnaefs
      use m_grdsam,     only: grdsam
      use dftu_specs,   only: switch_dftu     ! This variable determines whether
                                              !   the subroutine to compute the
                                              !   Hubbard terms should be called
                                              !   or not
      use m_dftu,       only: hubbard_term    ! Subroutine that compute the
                                              !   Hubbard terms
      use m_stress
      use m_energies
      use m_ntm
      use m_spin,         only: spin
      use spinorbit,      only: spinorb
      use m_dipol
      use m_forces,          only: fa
      use alloc, only: re_alloc, de_alloc
      use m_hsx, only: write_hsx
      use sys, only: die
      use fdf
#ifdef HAVE_LIBINT
      use siesta_options, only: hfx_wanted
      use nao2gto_data, only: hfx_libint, hfx_libderiv, 
     &                        hfx_options, hfx_system
      use nao2gto_hfx,  only: setup_hfx
      use nao2gto_hfx_force, only: setup_hfx_force
#endif
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
      use parallel, only: IOnode
#ifdef CDF
#ifdef NCDF_4
      use dictionary
      use m_ncdf_siesta
#endif
#endif
      use siesta_options, only : idyn, ia1, write_tshs_history
      use sparse_matrices, only: H_2D, S_1D
      use files,        only : label_length
      use m_ts_options, only : TS_HS_save
      use m_ts_kpoints, only : ts_Gamma, ts_kscell, ts_kdispl
      use m_ts_io,      only : ts_write_TSHS,fname_TSHS, FC_index

#ifdef FINAL_CHECK_HS
      use m_compute_max_diff, only: compute_max_diff
#endif
      
      implicit none

      ! MD-step, SCF-step
      integer,   intent(in) :: istep, iscf
      logical,   intent(in) :: SCFconverged
      real(dp)              :: stressl(3,3)
      real(dp), pointer     :: fal(:,:)   ! Local-node part of atomic F
      integer               :: ifa     ! Calculate forces? 0=>no, 1=>yes
      integer               :: istr    ! Calculate stress? 0=>no, 1=>yes
      real(dp)              :: g2max
      ! To avoid overwriting the current Hamiltonian and S
      real(dp), pointer :: H_tmp(:,:) => null()
      real(dp), pointer :: S_tmp(:) => null()
      real(dp), pointer :: H_dftu(:,:)
#ifdef FINAL_CHECK_HS
      real(dp) :: diff_H
#endif
#ifdef MPI
      real(dp)              :: stresstmp(3,3)
      real(dp), pointer     :: fatmp(:,:)
      real(dp)              :: buffer1
#endif
      character(len=label_length+13) :: fname
      integer :: io_istep, io_ia1
      logical :: not_using_auxcell
#ifdef CDF
#ifdef NCDF_4
      type(dictionary_t) :: dic_save
#endif
#endif
!------------------------------------------------------------------------- BEGIN

      not_using_auxcell = (no_s == no_u)
      
!     Initialize Hamiltonian ........................................
      call re_alloc(H_tmp, 1,maxnh, 1, spin%H,
     $             'H_tmp','final_H_f_stress')

      H_tmp(:,:) = 0.0_dp

!     Initialize forces and stress ...................
      nullify(fal)
      call re_alloc( fal, 1, 3, 1, na_u, 'fal', 'final_H_f_stress' )
      
      fa(1:3,1:na_u)   = 0.0_dp
      fal(1:3,1:na_u)  = 0.0_dp
      stress(1:3,1:3)  = 0.0_dp
      stressl(1:3,1:3) = 0.0_dp
      

! Neutral-atom: energy, forces and stress ............................
      call naefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, Ena, fa, stress, forces_and_stress=.true. )
      call dnaefs( na_u, na_s, scell, xa, indxua, rmaxv,
     &            isa, DEna, fa, stress, forces_and_stress=.true. )
      Ena = Ena + DEna

!     Kinetic: energy, forces, stress and matrix elements .................
      call kinefsm( na_u, na_s, no_s, scell, xa, indxua, rmaxo,
     &              maxnh, maxnh, lasto, iphorb, isa, 
     &              numh, listhptr, listh, numh, listhptr, listh, 
     &              spin%spinor, Dscf, Ekin, fal, stressl, H_tmp,
     &              matrix_elements_only=.false. ) 
#ifdef MPI
! Global reduction of energy terms
      call globalize_sum(Ekin,buffer1)
      Ekin = buffer1
#endif
! ..................

! Non-local-pseudop: energy, forces, stress and matrix elements .......
      call nlefsm( scell, na_u, na_s, isa, xa, indxua, 
     &             maxnh, maxnh, lasto, lastkb, iphorb, iphKB, 
     &             numh, listhptr, listh, numh, listhptr, listh, 
     &             spin%spinor, Dscf, Enl, fal, stressl, H_tmp,
     &             matrix_elements_only=.false. ) 

#ifdef MPI
! Global reduction of energy terms
      call globalize_sum(Enl,buffer1)
      Enl = buffer1
#endif

! If in the future the spin-orbit routine is able to compute
! forces and stresses, then "last" will be needed. If we are not
! computing forces and stresses, calling it in the first iteration
! should be enough
!
      if ( spin%SO ) then
         call spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     &        maxnh,numh,listhptr,listh,Dscf,H_tmp(:,3:),Eso)
      else
         Eso = 0._dp
      endif
      
!     Non-SCF part of total energy
      call update_E0()

! Hubbard term for LDA+U: energy, forces, stress and matrix elements ....
      if( switch_dftu ) then
         H_dftu => val(H_dftu_2D)
         call hubbard_term(scell, na_u, na_s, isa, xa, indxua,
     .        maxnh, maxnh, lasto, iphorb, no_u, no_l, 
     .        numh, listhptr, listh, numh, listhptr, listh,
     .        spin%spinor, Dscf, Edftu, DEdftu, H_dftu, 
     .        fal, stressl, H_tmp, iscf,
     .        matrix_elements_only=.false.)
#ifdef MPI
         ! Global reduction of energy terms
         call globalize_sum(Edftu,buffer1)
         Edftu = buffer1
         ! DEdftu should not be globalized as it 
         ! is based on globalized occupations
#endif
         Edftu = Edftu + DEdftu
      endif
! ..................


!     Non-local-pseudop: energy, forces, stress and matrix elements
!     Add SCF contribution to energy and matrix elements 
      g2max = g2cut
!     Last call to dhscf and grid-cell sampling if requested
      ifa  = 1
      istr = 1

      ! This will call dhscf with the final DM coming out of
      ! the scf cycle

      call grdsam( spin%Grid, no_s, iaorb, iphorb, 
     &             no_l, no_u, na_u, na_s, isa, xa, indxua,
     &             ucell, ntm, ifa, istr, maxnh,
     &             maxnh, numh, listhptr, listh, Dscf, Datm, H_tmp,
     &             Enaatm, Enascf, Uatm, Uscf, DUscf, DUext,
     &             Exc, Dxc, dipol, stress, fal, stressl )

#ifdef HAVE_LIBINT
      ! Update the Hartree-Fock hamiltonian and forces
      if ( hfx_wanted ) then
        call setup_hfx(hfx_libint, hfx_options, hfx_system, H_tmp)
        call setup_hfx_force(hfx_libint, hfx_libderiv, hfx_options,
     &    hfx_system, fal)
      endif
#endif

      ! Orthonormalization forces
      call overfsm( na_u, na_s, no_l, no_s, scell, xa, indxua, rmaxo,
     &              lasto, iphorb, isa, maxnh, numh, listhptr, listh,
     &              spin, Escf, fal, stressl ) 

!     Metadynamics forces
      if (lMetaForce) then
        call meta(xa,na_u,ucell,Emeta,fa,stress,.true.,.true.)
      endif

!     Add on force field contribution if required
      call twobody( na_u, xa, isa, ucell, Emm,
     &              ifa=1, fa=fa, istr=1, stress=stress )

#ifdef MPI
!     Global reduction of forces and stresses
      nullify(fatmp)
      call re_alloc( fatmp, 1, 3, 1, na_u, 'fatmp', 'final_H_f_stress' )
      call globalize_sum(stressl(1:3,1:3),stresstmp(1:3,1:3))
      call globalize_sum( fal(1:3,1:na_u), fatmp(1:3,1:na_u) )
      stress(1:3,1:3) = stress(1:3,1:3) + stresstmp(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fatmp(1:3,1:na_u)
      call de_alloc( fatmp, 'fatmp', 'final_H_f_stress' )
#else
      stress(1:3,1:3) = stress(1:3,1:3) + stressl(1:3,1:3)
      fa(1:3,1:na_u) = fa(1:3,1:na_u) + fal(1:3,1:na_u)
#endif

      ! If backward compatibility is requested, recover
      ! the old behavior with respect to H

#ifdef FINAL_CHECK_HS
      ! We also print-out the max-difference between the
      ! final H after SCF and the final computed H used for
      ! the forces.
      ! This value should correspond to dHmax
      call compute_max_diff(H, H_tmp, diff_H)
      if ( IONode ) then
         ! Print out the final differences
         write(*,"(a,f10.6)") ":!: Final H_scf - H_force (eV) : ",
     &        diff_H/eV
      end if
#endif
      
      if (recompute_H_after_scf) then
         if (ionode) then
            write(6,"(a)") ":!: Updating H after scf cycle" //
     $                     " as requested by compat. option"
         endif
         H = H_tmp
      endif

      call de_alloc( fal, 'fal', 'final_H_f_stress' )
      call de_alloc( H_tmp, 'H_tmp', 'final_H_f_stress' )

      ! Determine whether anything should be saved
      ! If not, return immediately
      if ( .not. SCFconverged ) return

! Save Hamiltonian and overlap matrices ............................
! Only in HSX format now.  Use Util/HSX/hsx2hs to generate an HS file
!
! Note that we save the Hamiltonian coming out of the scf cycle,
! not the one computed (from DM_out) in this routine.
! This call could be moved to a more appropriate place
!
      if (savehs) then
         call write_hsx( no_u, no_s, spin%H, indxuo,
     &        maxnh, numh, listhptr, listh, H, S, Qtot,
     &        Temp, xijo)
      endif
#ifdef CDF
#ifdef NCDF_4
      if ( write_cdf ) then
        if ( idyn == 6 ) then
          call cdf_save_fc(trim(slabel)//'.nc', istep)
          if ( istep == 0 ) then
            dic_save = ('fa'.kv.1)//('stress'.kv.1)//('Ef'.kv.1)
            if ( savehs .or. TS_HS_save )
     &          dic_save = dic_save//('H'.kv.1)
            call cdf_save_state(trim(slabel)//'.nc',dic_save)
            call delete(dic_save)
          end if
        else if ( idyn /=0 .or. nmove /= 0 ) then
          !call cdf_save_md(trim(slabel)//'.nc')
        else
          dic_save = ('fa'.kv.1)//('stress'.kv.1)//('Ef'.kv.1)
          if ( savehs .or. TS_HS_save ) dic_save = dic_save //('H'.kv.1)
          call cdf_save_state(trim(slabel)//'.nc',dic_save)
          call delete(dic_save)
        end if
      end if
#endif
#endif

      if ( fixspin .and. (write_tshs_history .or. TS_HS_save) ) then
        ! For fixed spin calculations we shift the Hamiltonian according
        ! to the first spin such that we have a "common" E_F, after writing
        ! we shift back again.
        H_tmp => val(H_2D)
        S_tmp => val(S_1D)
        call daxpy(size(S_tmp),Efs(1)-Efs(2),S_tmp(1),1,H_tmp(1,2),1)
        Ef = Efs(1)
      end if
      
      if ( write_tshs_history ) then
         ! This is "pure" MD and we only write consecutive numbers
         ! Together with this you cannot also save FC
         fname = fname_TSHS(slabel, istep = istep )
         call ts_write_tshs(fname, .false., not_using_auxcell, ts_Gamma,
     &        ucell, nsc, isc_off, na_u, no_s, spin%H, 
     &        ts_kscell, ts_kdispl, 
     &        xa, lasto, 
     &        H_2D, S_1D, indxuo,
     &        Ef, Qtot, Temp, istep, 0)
      else if ( TS_HS_save ) then
         if ( idyn == 6 ) then
            ! Correct the FC indices such that 
            ! the TSHS format contains, atom displaced, index of 
            ! displacement:
            !  1 = -x
            !  2 = +x
            !  3 = -y
            !  4 = +y
            !  5 = -z
            !  6 = +z
            call FC_index(istep,ia1,io_istep,io_ia1)
            fname = fname_TSHS(slabel,istep = io_istep, ia1 = io_ia1)
            call ts_write_tshs(fname, .false.,
     &           not_using_auxcell, ts_Gamma,
     &           ucell, nsc, isc_off, na_u, no_s, spin%H, 
     &           ts_kscell, ts_kdispl, 
     &           xa, lasto, 
     &           H_2D, S_1D, indxuo,
     &           Ef, Qtot, Temp, io_istep, io_ia1)
         else
            fname = fname_TSHS(slabel)
            call ts_write_tshs(fname, .false.,
     &           not_using_auxcell, ts_Gamma,
     &           ucell, nsc, isc_off, na_u, no_s, spin%H,
     &           ts_kscell, ts_kdispl, 
     &           xa, lasto,
     &           H_2D, S_1D, indxuo,
     &           Ef, Qtot, Temp, 0, 0)
         end if
      endif

      if ( fixspin .and. (write_tshs_history .or. TS_HS_save) ) then
        ! Shift back
        call daxpy(size(S_tmp),Efs(2)-Efs(1),S_tmp(1),1,H_tmp(1,2),1)
      end if

!----------------------------------------------------------------------- END
      END subroutine final_H_f_stress
      END module m_final_H_f_stress
