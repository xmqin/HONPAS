! This file is part of the HONPAS package.
!

! Module coded by xmqin, October 2013
!
! calculate HFX(u,v) 
!
! HFX(u,v) = (um|vn) * Dm (m,n)
!
      module hfx_potential
      use kinds,           only : dp, int_8
      use atm_types,       only : l_max, nspecies, &
                                  maxn_contract, maxn_orbnl
      use hfx_types
      use alloc,           only :  re_alloc, de_alloc
      use parallel,        only :  Node, Nodes, BSize
      use parallelsubs,    only :  set_bsize_NAO2GTO
      use atmfuncs,        only :  lofio, mofio
      use listsc_module,   only :  listsc
      use extended_index,  only :  extended_index_init 
      use prescreen
      use atomlist,        only :  indxuo
      use nao2gto_eri ,    only :  evaluate_eri
      use gamma_F,         only :  init_md_ftable
      use contract_eri,    only :  calc_contract_eri
      use libint_wrapper_types,   only : lib_int 
      use libint_wrapper,         only : initialize_libint
      use fdf

#ifdef MPI
      use mpi_siesta
#endif
!      use atm_
      implicit none
      private

      public build_hfx_potential !, on_the_fly

      contains

!------------------------------------------------------------------------------------

      subroutine build_hfx_potential( nspin, norb, iaorb, iphorb, nuo, nuotot, na, isa, &
                                      xa, indxua, cell, nsc, maxnh, numh, listhptr,     &
                                      listh, samexa, Dscf, Dscf_max, HFX )

!---------------------------- INPUT  VARIABLES -------------------------------------
      integer, intent(in)     :: &
        maxnh, na, norb, nspin, nuo, nuotot, &
        iaorb(norb), indxua(na), iphorb(norb), isa(na),  &
        listh(maxnh), listhptr(nuo), &
        numh(nuo), nsc(3)
      real(dp),  intent(in)   :: xa(3,na), cell(3,3)
      logical,   intent(in)   :: samexa
      real(dp),  intent(in)   :: Dscf(maxnh,nspin)
      real(dp),  intent(in)   :: Dscf_max(maxnh)

!--------------------------- INOUT VARIABLES ---------------------------------------
      real(dp), intent(inout) :: HFX(maxnh, nspin)

!------------------------- TEMPOS, INTERNAL VARIABLES ------------------------------
      integer ncells, nua, i, ind, j, io, iu, ju, iuo, juo, jo, jjo, ioa, &
              is, l, m, nshells, ispin, ia, js, ishell,jshell,ipgf, jpgf, &
              l_j,m_j,joa

      real(dp) tmax, gint, time_start, time_end
!      logical, pointer, save    ::  um_cut(:,:)  

      real(dp)                  ::  scell(3,3), rscell(3,3)
      real(dp), save            ::  max_eri
!      real(dp), pointer, save   ::  eri_prescreen(:)  !,S_tmp
      type(hfx_input_parameter), save :: hfx_parameter
      type(lib_int)             ::  lib
      logical, save             ::  frstme = .true.
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

#ifdef MPI
      integer       :: bsizedefault
      integer       :: MPIerror
#endif

      external memory, timer
!
! Find supercell and ncells

      ncells=nsc(1)*nsc(2)*nsc(3)
      nua = na/ncells

      do  i = 1,3
        do j = 1,3
          scell(j,i) = cell(j,i) * nsc(i)
        enddo
      enddo

      call reclat(scell, rscell, 0)

! If this routine is called first time, read and print hfx input information
      if(frstme) then
        call read_hfx_info(hfx_parameter)
        call print_hfx_info(hfx_parameter)

! supercell orbital initialisation
        call extended_index_init(nsc, nuotot)

! gamma function table and Libint initialisation 

        call init_md_ftable(4*l_max)
        call initialize_libint(lib, l_max)

!  orbital shell initialisation
        nullify(subshell)
        call re_alloc(subshell,1,norb)
        subshell(1:norb) = 0
! 2l+1 orbitals form a shell
        i = 0
        do io = 1,norb
           is = isa(iaorb(io))
           ioa = iphorb(io)
           l = lofio(is,ioa)
           m = mofio(is,ioa)
           if(m.eq.-l) then
              i = i + 1
              subshell(io) = i
           else
              subshell(io) = i
           endif
        enddo
        nshells = subshell(norb)
!        if(node.eq.0) &
!           write(6,'(A, I8, A, I8)') &
!           'Supercell orbitals = ', norb, ' shells = ', nshells
      endif

        
      if(frstme .or. (.not.samexa)) then
        nullify(eri_prescreen)
        call re_alloc( eri_prescreen, 1, maxnh,  &
               name='eri_prescreen', routine='build_hfx_potential' )  
        eri_prescreen(maxnh)=0.0d0

!  To calculate ERI prescreen matrix

        call calc_prescreen_eri( lib, norb, iaorb, iphorb, nuotot, nuotot, &
                                 na, isa, xa, indxua, scell, rscell, &
                                 maxnh, numh, listhptr, listh, &
                                 hfx_parameter, max_eri)


! --- xmqin add screenfunc metond ! -----------------------------
        nullify(pair_dist_radii_pgf)
        nullify(sfc_pgf)
        nullify(sfc_shell)
        nullify(sfc_kind)

        allocate(pair_dist_radii_pgf(maxn_contract, maxn_contract, &
                       maxn_orbnl, maxn_orbnl, nspecies, nspecies))
        allocate(sfc_pgf(maxn_contract, maxn_contract, &
                       maxn_orbnl, maxn_orbnl, nspecies, nspecies))
        allocate(sfc_shell(maxn_orbnl, maxn_orbnl, &
                                                 nspecies, nspecies))
        allocate(sfc_kind(nspecies, nspecies))

      DO is = 1,nspecies
        DO js = 1,nspecies
           sfc_kind(js,is)%x(:) = 0.0_dp
          DO ishell=1,maxn_orbnl
            DO jshell=1,maxn_orbnl
               sfc_shell(jshell,ishell,js,is)%x(:) = 0.0_dp
              DO ipgf=1,maxn_contract
                DO jpgf=1,maxn_contract
                   pair_dist_radii_pgf(jpgf,ipgf,jshell,&
                                       ishell,js,is)%x(:) = 0.0_dp
                   sfc_pgf(jpgf,ipgf,jshell, &
                                            ishell,js,is)%x(:) = 0.0_dp
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

      call calc_pair_dist_radii(hfx_parameter, pair_dist_radii_pgf)

      call calc_screening_functions(lib, hfx_parameter, scell, rscell, &
                                    sfc_pgf, sfc_shell,&
                                    pair_dist_radii_pgf)

      coeffs_kind_max0=MAXVAL(sfc_shell(:,:,:,:)%x(2))
      if(node.eq.0) write(6,*) "max0 = ", coeffs_kind_max0

!! build left and right orbital pair_list
!        max_element1 = ncells*nuotot*(nuotot+1)/2
!        max_element2 = ncells**2*nuotot*(nuotot+1)/2
       if (allocated(list_ij%element)) then
         deallocate(list_ij%element)
       endif
       if (allocated(list_kl%element)) then
         deallocate(list_kl%element)
       endif
       allocate(list_ij%element(ncells*nuotot*(nuotot+1)/2))
       allocate(list_kl%element(ncells**2*nuotot*(nuotot+1)/2))

        call build_pair_list( nua, na, nua, isa, xa, indxua, nuotot, norb, iaorb, &
                              iphorb, maxnh, numh, listhptr, listh, scell, rscell, hfx_parameter,   &
                              max_eri, list_ij)

        call build_pair_list( na, na, nua, isa, xa, indxua, nuotot, norb, iaorb,      &
                              iphorb, maxnh, numh, listhptr, listh, scell, rscell, hfx_parameter, &
                              max_eri, list_kl)

        log10_eps_schwarz = LOG10(hfx_parameter%eps_schwarz)

        if(node.eq.0) then
!          write(6,*) ncells*nuotot*(nuotot+1)/2, ncells**2*nuotot*(nuotot+1)/2
          write(6,'(a,2x,I12)')  'list_ij:',list_ij%nelement
          write(6,'(a,2x,I12)')  'list_kl:',list_kl%nelement
          write(6,'(a,f12.9)') 'max_eri', max_eri
          write(6,*) hfx_parameter%eps_schwarz, log10_eps_schwarz

        endif

!#ifdef MPI
!        if (Node.eq.0) then
!          call set_bsize_NAO2GTO(Nodes,list_kl%nelement,bsizedefault)
!          BSize = fdf_integer('blocksize_NAO2GTO',bsizedefault)
!         write(6,*) "Blocksize for ERIs distribution :", BSize
!       endif
!        call MPI_Bcast(BSize,1,MPI_integer,0,MPI_Comm_World,MPIerror)
!#endif

      endif

      
      frstme = .false.

!! calculate ERIs and store them on RAM:

      call timer('HFX',1)
!      call cpu_time(time_start) 
      call evaluate_eri( lib, nspin, norb, iaorb, iphorb, nuotot, na, isa,   &
                          maxnh, numh, listhptr, listh, &
                          scell, rscell, hfx_parameter, Dscf, Dscf_max, &
          !                pair_dist_radii_pgf, sfc_pgf, sfc_shell,            &
                          log10_eps_schwarz, HFX )

!      call cpu_time(time_end)
!      write(6,'(a, f12.6, a, I6)') "ERIs time = ", time_end-time_start, &
!                                   "(Sec), on Node =", node
!
      call timer('HFX',2)

      end subroutine build_hfx_potential

      end module hfx_potential

