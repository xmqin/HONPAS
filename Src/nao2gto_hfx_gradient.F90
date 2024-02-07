! This file is part of the HONPAS package.
!

! Module coded by xmqin, October 2013
!
! calculate HFX(u,v) 
!
! HFX(u,v) = (um|vn) * Dm (m,n)
!
      module hfx_gradient
      use kinds,           only : dp, int_8
      use atm_types,       only : l_max, nspecies, &
                                  maxn_contract, maxn_orbnl
      use hfx_types
      use alloc,           only :  re_alloc, de_alloc
      use parallel,        only :  Node, Nodes
      use atmfuncs,        only : lofio, mofio
      use listsc_module,   only :  listsc
      use extended_index,  only :  extended_index_init 
      use prescreen
      use nao2gto_gradient,    only :  evaluate_gradient
      use atomlist,        only :  indxuo
      use gamma_F,         only :  init_md_ftable
      use contract_eri,    only :  calc_contract_eri
      use libint_wrapper_types,   only : lib_int, lib_deriv
      use libint_wrapper,         only : initialize_libint, &
                                         initialize_libderiv

      implicit none
      private

      public build_hfx_gradient !, on_the_fly

      contains

!------------------------------------------------------------------------------------

      subroutine build_hfx_gradient( nspin, norb, iaorb, iphorb, nuo, nuotot,nua, na, isa, &
                                      xa, indxua, cell, nsc, maxnh, numh, listhptr,     &
                                      listh, Dscf, Dscf_max, Fal )

!---------------------------- INPUT  VARIABLES -------------------------------------
      integer, intent(in)     :: &
        maxnh, nua, na, norb, nspin, nuo, nuotot, &
        iaorb(norb), indxua(na), iphorb(norb), isa(na),  &
        listh(maxnh), listhptr(nuo), &
        numh(nuo), nsc(3)
      real(dp),  intent(in)   :: xa(3,na), cell(3,3)
      real(dp),  intent(in)   :: Dscf(maxnh,nspin)
      real(dp),  intent(in)   :: Dscf_max(maxnh)

!--------------------------- INOUT VARIABLES ---------------------------------------
      real(dp), intent(inout) :: Fal(3,nua)

!------------------------- TEMPOS, INTERNAL VARIABLES ------------------------------
      integer ncells, i, ind, j, io, iu, ju, iuo, juo, jo, jjo, ioa, &
              is, l, m, nshells, ispin, ia,js, ishell,jshell,ipgf, jpgf,&
              l_j,m_j,joa

      real(dp) tmax, gint, time_start, time_end
!      logical, pointer, save    ::  um_cut(:,:)

      real(dp)                  ::  scell(3,3), rscell(3,3)
      real(dp), save            ::  max_eri
!      real(dp), pointer, save   ::  eri_prescreen(:)  !,S_tmp
      type(hfx_input_parameter), save :: hfx_parameter
      type(lib_int)             ::  lib
      type(lib_deriv)           ::  deriv
!      logical, save             ::  frstme = .true.
!      type(hfx_screen_coeff_type), &
!      dimension(:, :, :, :, :, :), pointer, save  :: &
!                                   pair_dist_radii_pgf, &
!                                   sfc_pgf
!      type(hfx_screen_coeff_type), &
!      dimension(:, :, :, :), pointer, save  :: &
!                                   sfc_shell
!      type(hfx_screen_coeff_type), &
!      dimension(:, :), pointer, save  :: &
!                                   sfc_kind
      real(dp), save :: coeffs_kind_max0, log10_eps_schwarz



      external memory, timer
!
! Find supercell and ncells

      ncells=nsc(1)*nsc(2)*nsc(3)

      do  i = 1,3
        do j = 1,3
          scell(j,i) = cell(j,i) * nsc(i)
        enddo
      enddo

      call reclat(scell, rscell, 0)

! If this routine is called first time, read and print hfx input information
        call read_hfx_info(hfx_parameter)
!       call print_hfx_info(hfx_parameter)

! supercell orbital initialisation
!        call extended_index_init(nsc, nuotot)

! gamma function table and Libint initialisation 

        call init_md_ftable(4*l_max)
        call initialize_libint(lib, l_max)
        call initialize_libderiv(deriv,l_max)        

        log10_eps_schwarz = LOG10(hfx_parameter%eps_schwarz)



        if(node.eq.0) then
!          write(6,*) ncells*nuotot*(nuotot+1)/2, ncells**2*nuotot*(nuotot+1)/2
          write(6,'(a,2x,I12)')  'list_ij:',list_ij%nelement
          write(6,'(a,2x,I12)')  'list_kl:',list_kl%nelement
!          write(6,'(a,f12.9)') 'max_eri', max_eri
        endif

!      endif
      
      call timer('HFX_gradient',1)    
!      call cpu_time(time_start)
      call  evaluate_gradient( deriv, nspin, norb, iaorb, iphorb, nuotot, nua, na, isa,    &
                               maxnh, numh, listhptr, listh, &
                               scell, rscell, hfx_parameter, Dscf, Dscf_max, &
                               !pair_dist_radii_pgf, sfc_pgf, sfc_shell,   &
                               log10_eps_schwarz, Fal )


!      call cpu_time(time_end)
!      if(node.eq.0)write(6,'(a, f12.6)') "ERIs time = ", time_end-time_start
!                                   "(Sec), on Node =", node

      call timer('HFX_gradient',2)

      end subroutine build_hfx_gradient

      end module hfx_gradient

