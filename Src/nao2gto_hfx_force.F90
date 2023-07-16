! This file is part of the HONPAS package.
!
!
! Module coded by xmqin, Oct.12, 2013
! Edited by  xmqin, Nov. 27, 2016
!
module nao2gto_hfx_force
  use precision, only: dp
  use atm_types, only: maxn_orbnl, nspecies
  use atm_types, only: species, species_info
  use atmfuncs, only: lofio, mofio
  use atomlist, only: indxuo
  use listsc_module, only :  listsc
  use parallel, only: Node, Nodes
  use parallelsubs, only: GetNodeOrbs, GlobalToLocalOrb, &
                          LocalToGlobalOrb, WhichNodeOrb

!  use nao2gto_common, only: dp, l_max, maxn_contract
  use nao2gto_contract, only: calc_contract_eri, calc_contract_deriv_eri
  use nao2gto_data, only: nco, nso
  use nao2gto_data, only: D2Sindx
  use nao2gto_index
  use nao2gto_prescreen
  use nao2gto_types
  use nao2gto_utils
#ifdef MPI
    use mpi_siesta
#endif

  implicit none

  private

  public :: setup_hfx_force

  real(dp), save :: max_eri_force = 0.0_dp

contains

  subroutine setup_hfx_force(libint_data, deriv_data, hfx_optdata, &
               hfx_sysdata, Fal)

    use alloc, only: re_alloc, de_alloc
    use parallel, only: Node, Nodes
    use parallelsubs, only: GetNodeOrbs, GlobalToLocalOrb, &
                            LocalToGlobalOrb, WhichNodeOrb
    use sparse_matrices, only: Dscf

    use nao2gto_dm
    use nao2gto_libint
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    type(Libderiv_t), intent(inout) :: deriv_data
    type(hfx_options_type), intent(in) :: hfx_optdata
    type(hfx_system_type), intent(in) :: hfx_sysdata
    real(dp), intent(inout) :: fal(3,hfx_sysdata%nua)

    ! Local variables
    integer :: i, ia, io, iio, is, ispin, j, n
    integer :: iu, ju, jo, ind, num_u ,num_v, num_m, num_n
    real(dp), pointer :: DM_tmp(:,:,:) => null(), P_max(:,:) => null()

#ifdef MPI
    ! Global buffers for the storage of sparse matrix
    integer :: BNode, MPIerror, maxnhg, maxndg, maxnumh, nuog
    integer, pointer :: listhg(:) => null(), listhptrg(:) => null(), &
      numhg(:) => null()
    real(dp), pointer :: Dscfg(:,:) => null(), Dscf_max(:) => null()
    integer,  dimension(:),     pointer :: listdptrg => null()
    integer,  dimension(:),     pointer :: numdg => null()
    integer,  dimension(:),     pointer :: listdg => null()
#endif

    external :: timer

    ! -------------------------------------------------------------------------

    call timer('HFX_FORCE', 1)

#ifdef MPI
    call re_alloc(numhg, 1, hfx_sysdata%nuotot, name='numhg', &
      routine='setup_hfx_force')
    call re_alloc(listhptrg, 1, hfx_sysdata%nuotot, name='listhptrg', &
      routine='setup_hfx_force')
    call re_alloc(numdg,     1, hfx_sysdata%nuotot, name='numdg',     &
      routine='setup_hfx_force')
    call re_alloc(listdptrg, 1, hfx_sysdata%nuotot, name='listdptrg', &
      routine='setup_hfx_force')


    ! Globalise numh
    do io = 1, hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        numhg(io) = hfx_sysdata%numh(iio)
        numdg(io) = hfx_sysdata%numh(iio)
      endif
      call MPI_Bcast(numhg(io), 1, MPI_integer, BNode, &
        MPI_Comm_World, MPIerror)
      call MPI_Bcast(numdg(io), 1, MPI_integer, BNode, &
        MPI_Comm_World, MPIerror)
    enddo

!   Build global listhptr
    listhptrg(1) = 0
    listdptrg(1) = 0
    do io = 2, hfx_sysdata%nuotot
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      listdptrg(io) = listdptrg(io-1) + numdg(io-1)
    enddo

    ! Globalse listh
    maxnhg = listhptrg(hfx_sysdata%nuotot) + numhg(hfx_sysdata%nuotot)
    maxndg = listdptrg(hfx_sysdata%nuotot) + numdg(hfx_sysdata%nuotot)
    call re_alloc(listhg, 1, maxnhg, name='listhg', routine='setup_hfx_force')
    call re_alloc(listdg, 1, maxndg, name='listdg', routine='setup_hfx_force')

    do io = 1, hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        do jo = 1, numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
            hfx_sysdata%listh(hfx_sysdata%listhptr(iio)+1: &
              hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
          listdg(listdptrg(io)+1:listdptrg(io)+numdg(io)) = hfx_sysdata%listh( &
            hfx_sysdata%listhptr(iio)+1:hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
        enddo
      endif

      call MPI_Bcast(listhg(listhptrg(io)+1), numhg(io), MPI_integer, &
        BNode, MPI_Comm_World, MPIerror)
      call MPI_Bcast(listdg(listdptrg(io)+1), numdg(io), MPI_integer,&
        BNode, MPI_Comm_World, MPIerror)
    enddo

    call re_alloc(Dscfg, 1, maxndg, 1, hfx_sysdata%nspin, name='Dscfg', &
      routine='setup_hfx_force')
    Dscfg(:,:) = 0.0_dp

    do io = 1, hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        do ispin = 1, hfx_sysdata%nspin
          do jo = 1, hfx_sysdata%numh(iio)
            Dscfg(listdptrg(io)+jo,ispin) = &
              Dscf(hfx_sysdata%listhptr(iio)+jo,ispin)
          enddo
        enddo
      endif
      do ispin = 1, hfx_sysdata%nspin
        call MPI_Bcast(Dscfg(listdptrg(io)+1,ispin), numdg(io), &
          MPI_double_precision, BNode, MPI_Comm_World, MPIerror)
      enddo
    enddo
#endif

    !--------------------end for global Dscfg and Hmatg--------------------


#ifdef MPI
    ! Transform sparse Dm to full matrix
    call re_alloc(Dscf_max, 1, maxndg, name='Dscf_max', routine='setup_hfx_force')
    Dscf_max(:) = 0.0_dp

    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, maxndg, numdg, listdptrg, &
      listdg, D2Sindx, Dscfg, Dscf_max)
    call timer('ERI_deriv',1)

    call build_hfx_gradient(libint_data, deriv_data, hfx_optdata, &
      hfx_sysdata, maxnhg, numhg, listhptrg, listhg, Dscfg, Dscf_max, Fal)

    call timer('ERI_deriv',2)
#else
    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%maxnh, &
      hfx_sysdata%numh, hfx_sysdata%listhptr, hfx_sysdata%listh, &
      D2Sindx, Dscf, Dscf_max )

    call timer('ERI_deriv',1)

    call build_hfx_gradient(libint_data, deriv_data, hfx_optdata, &
      hfx_sysdata, hfx_sysdata%maxnh, hfx_sysdata%numh, &
      hfx_sysdata%listhptr, hfx_sysdata%listh, Dscf, Dscf_max, Fal)
 
    call timer('ERI_deriv',2)
#endif

    call de_alloc(Dscf_max, name='P_max', routine='setup_hfx_force')

#ifdef MPI
    call de_alloc(listhg, name='listhg', routine='setup_hfx_force')
    call de_alloc(listhptrg, name='listhptrg', routine='setup_hfx_force')
    call de_alloc(numhg, name='numhg', routine='setup_hfx_force')
    call de_alloc(Dscfg, name='Dscfg', routine='setup_hfx_force')
#endif

!!   For debugging
!    do io = 1, hfx_sysdata%nua
!      write(6,'(a,i5,3f12.5)')             &
! &      'setup_hfx_force: io, fal = ',     &
! &                        io, fal(:,io)
!    enddo 
!!   Ecd debugging

    call timer('HFX_FORCE', 2)

  end subroutine setup_hfx_force

  ! ***** HFX_GRADIENT MODULE BEGIN *****
  ! This file is part of the HONPAS package.
  !
  ! Module coded by xmqin, October 2013
  !
  ! calculate HFX(u,v)
  !
  ! HFX(u,v) = (um|vn) * Dm (m,n)
  !
  ! ***** HFX_GRADIENT MODULE END *****
  ! ***** GRADIENT MODULE BEGIN *****
  ! Written by xmqin, October 2013
  ! ***** GRADIENT MODULE END *****
  subroutine build_hfx_gradient(libint_data, deriv_data, hfx_opts, hfx_sys, &
               maxnh, numh, listhptr, listh, Dscf, Dscf_max, Fal)

    use precision, only: i8b 
    use alloc, only: de_alloc, re_alloc
    use parallel, only: IOnode
!    use nao2gto_common
    use nao2gto_data, only: D2Sindx, eri_prescreen, hfx_gtos, list_ij, list_kl, &
                            pair_dist_radii_pgf, sfc_pgf, sfc_shell
    use nao2gto_libint
!    use nao2gto_prescreen
    use nao2gto_wrappers, only: nao2gto_libderiv_init

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    type(Libderiv_t), intent(inout) :: deriv_data
    type(hfx_options_type), intent(in) :: hfx_opts
    type(hfx_system_type), intent(in) :: hfx_sys
    integer, intent(in) :: maxnh
    integer, intent(in) :: numh(hfx_sys%nuotot)
    integer, intent(in) :: listhptr(hfx_sys%nuotot)
    integer, intent(in) :: listh(maxnh)
    real(dp), intent(in) :: Dscf(maxnh,hfx_sys%nspin)
    real(dp), intent(in) :: Dscf_max(maxnh)
    real(dp), intent(inout) :: Fal(3,hfx_sys%nua)

    ! Local variables
    real(dp) :: ri(3), rj(3), rk(3), rl(3), rij2,  rkl2
    real(dp) :: eps_temp, DM_max,symm_factor
    real(dp) :: max_contraction_val, max_val, max_val1, &
                max_val2, pmax_entry
    real(dp), pointer :: eri_deriv(:,:,:,:,:) => null()
    type(hfx_screen_coeff_type), dimension(:,:), pointer :: &
      tmp_R_1 => null(), tmp_R_2 => null(), &
      tmp_screen_pgf1 => null(), tmp_screen_pgf2 => null()
    real(dp) :: max_val2_set, log10_pmax
    real(dp) :: nao_eri_deriv(12)
    integer :: io, jo, ko, lo, is, js, ks, ls, ioa, joa, koa, loa,   &
               l_i, l_j, l_k, l_l, m_i, m_j, m_k, m_l,               &
               ncoi, ncoj, ncok, ncol, npgfi, npgfj, npgfk, npgfl,   &
               nsoi, nsoj, nsok, nsol, num_a, num_b, num_c, num_d,   &
               i_list_ij, i_list_kl, i_list_kl_local, list_kl_local, &
               index_ij, index_kl, ia, ja, ka, la , ishell, jshell,  &
               kshell, lshell
    integer :: ispin
    type(species_info), pointer :: ispp => null(), jspp => null(), &
      kspp => null(), lspp => null()
    type(gto_info_type), pointer :: igto => null(), jgto => null(), &
      kgto => null(), lgto => null()
#ifdef MPI
      integer :: MPIerror, Request, Status(MPI_Status_Size)
      integer :: num_loc
#endif
    integer(i8b) :: number_task
    logical :: flag
    integer :: i, buf, tag, task_frag
    integer :: i_list_ijkl, temp_i_list_kl, count, temp_count
    integer :: list_pos(2), temp_list_pos(2)

    ! -------------------------------------------------------------------------

    ! Calculate ERI prescreen matrix
    ! FIXME: Do we need to recalculate eri_prescreen now? (YP)
!    call calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
!      maxnh, numh, listhptr, listh, max_eri_force)

!    if ( IOnode ) write(*,'("max_eri_force = ",E18.6)') max_eri_force

      if(hfx_opts%Dynamic_parallel) then

      buf  =  11
      request = 1000
      tag = 111
      task_frag = hfx_opts%frag_size

      if( Node .eq. 0 ) then
        number_task = list_ij%nelement*list_kl%nelement
!        write(6,*) node, list_ij%nelement, list_kl%nelement, number_task

!        call cpu_time(time_start3)

        call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, MPI_ANY_TAG,&
                       MPI_COMM_WORLD, request, MPIerror)
        list_pos(1) = 1
        list_pos(2) = 1

        do while(number_task > 0)
            CALL MPI_TEST(request, flag, status, MPIerror)
!            CALL MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status, MPIerror)
            if(flag) THEN
                if(number_task >= task_frag ) then
                   call MPI_SEND(task_frag, 1, MPI_integer, status(MPI_SOURCE), &
                               status(MPI_TAG), MPI_COMM_WORLD, MPIerror)
                   call MPI_SEND(list_pos(1), 2, MPI_integer, status(MPI_SOURCE), &
                               status(MPI_TAG), MPI_COMM_WORLD, MPIerror)

                    number_task = number_task - task_frag
                    list_pos(1) = list_pos(1) + task_frag
                    do while(list_pos(1) > list_ij%nelement)
                        list_pos(1) = list_pos(1) - list_ij%nelement
                        list_pos(2) = list_pos(2) + 1
                    enddo
                    call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                                   MPI_COMM_WORLD, request, MPIerror)
                else
                    count = number_task
                    call MPI_SEND(count , 1, MPI_integer, status(MPI_SOURCE), &
                                status(MPI_TAG), MPI_Comm_World, MPIerror)
                    call MPI_SEND(list_pos(1), 2, MPI_integer, status(MPI_SOURCE), &
                                status(MPI_TAG), MPI_COMM_WORLD, MPIerror)

                   number_task = 0
                    call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                          MPI_COMM_WORLD, request, MPIerror)
                endif
            endif

        enddo

        call MPI_WAIT(request, status, MPIerror)
        call MPI_SEND(0, 1, MPI_integer, status(MPI_SOURCE), &
                      status(MPI_TAG), MPI_COMM_WORLD, MPIerror)

        do i = 1 , Nodes-2
            call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                        MPI_COMM_WORLD, request, MPIerror)
            call MPI_WAIT(request, status, MPIerror)
            call MPI_SEND(0, 1, MPI_integer, status(MPI_SOURCE), &
                          status(MPI_TAG), MPI_Comm_World, MPIerror)
        enddo

      else
        do while(.true.)
            call MPI_ISEND(buf, 1, MPI_integer, 0, tag, MPI_COMM_WORLD, request, MPIerror)
            CALL MPI_WAIT(request, status, MPIerror)
            CALL MPI_RECV(count, 1, MPI_integer, 0, tag, MPI_COMM_WORLD, status, MPIerror)
            if(count .eq. 0) then
                exit
            endif
            CALL MPI_RECV(list_pos(1), 2, MPI_integer, 0, tag, MPI_COMM_WORLD, status, MPIerror)
         i_list_kl = 0
         do i_list_ijkl = 1, count
            i_list_ij = list_pos(1) + i_list_ijkl - 1
            temp_i_list_kl = list_pos(2)

            do while(i_list_ij > list_ij%nelement)
                i_list_ij = i_list_ij - list_ij%nelement
                temp_i_list_kl = temp_i_list_kl + 1
            enddo

          io = list_ij%element(i_list_ij)%pair(1)
          jo = list_ij%element(i_list_ij)%pair(2)
          ri = list_ij%element(i_list_ij)%r1
          rj = list_ij%element(i_list_ij)%r2
          rij2 = list_ij%element(i_list_ij)%dist2

          ia = hfx_sys%iaorb(io)
          ja = hfx_sys%iaorb(jo)

          is = hfx_sys%isa(ia)
          ispp => species(is)
          igto => hfx_gtos(is)
          ioa = hfx_sys%iphorb(io)
          l_i = lofio(is, ioa)
          npgfi = igto%orbnl_contract(ispp%orb_index(ioa))
          ncoi = nco(l_i)*npgfi

          js = hfx_sys%isa(ja)
          jspp => species(js)
          jgto => hfx_gtos(js)
          joa = hfx_sys%iphorb(jo)
          l_j = lofio(js, joa)
          npgfj = jgto%orbnl_contract(jspp%orb_index(joa))
          ncoj = nco(l_j)*npgfj

          ishell = ispp%orb_index(ioa)
          jshell = jspp%orb_index(joa)
          index_ij = list_ij%element(i_list_ij)%nl_index
          max_val1 = sfc_shell(jshell,ishell,js,is)%x(1)*rij2 + &
                     sfc_shell(jshell,ishell,js,is)%x(2)

         if(temp_i_list_kl .ne. i_list_kl) then
            i_list_kl = temp_i_list_kl
            ko = list_kl%element(i_list_kl)%pair(1)
            lo = list_kl%element(i_list_kl)%pair(2)
            rk = list_kl%element(i_list_kl)%r1
            rl = list_kl%element(i_list_kl)%r2
            rkl2 =  list_kl%element(i_list_kl)%dist2

            ka = hfx_sys%iaorb(ko)
            la = hfx_sys%iaorb(lo)

            ks = hfx_sys%isa(ka)
            kspp => species(ks)
            kgto => hfx_gtos(ks)
            koa = hfx_sys%iphorb(ko)
            l_k = lofio(ks, koa)
            npgfk = kgto%orbnl_contract(kspp%orb_index(koa))
            ncok =nco(l_k)*npgfk

            ls = hfx_sys%isa(la)
            lspp => species(ls)
            lgto => hfx_gtos(ls)
            loa = hfx_sys%iphorb(lo)
            l_l = lofio(ls, loa)
            npgfl = lgto%orbnl_contract(lspp%orb_index(loa))
            ncol = nco(l_l)*npgfl

            kshell = kspp%orb_index(koa)
            lshell = lspp%orb_index(loa)

            index_kl = list_kl%element(i_list_kl)%nl_index
         endif

         if(ia.eq.ja.and.ia.eq.ka.and.ka.eq.la) cycle
         if ( index_kl .le. index_ij ) then

            if((D2Sindx(io,ko).eq.0).and.(D2Sindx(io,lo).eq.0).and.(D2Sindx(jo,ko).eq.0) &
              .and.(D2Sindx(jo,lo).eq.0) )  cycle

          max_val2_set = (sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
                          sfc_shell(lshell,kshell,ls,ks)%x(2) )

          max_val = max_val1 + max_val2_set

          eps_temp = eri_prescreen(D2Sindx(io,jo))*eri_prescreen(D2Sindx(ko,lo))
          eps_temp = dsqrt(eps_temp)

          if ( hfx_opts%DM_trunc ) then
             DM_max = 2.0_dp*max( Dscf_max(D2Sindx(io, ko))*Dscf_max(D2Sindx(jo,lo)), &
                        Dscf_max(D2Sindx(io, lo))*Dscf_max(D2Sindx(jo, ko)) )

            IF ( DM_max <= 0.0_dp ) THEN
              log10_pmax = log_zero
            ELSE
              log10_pmax = LOG10(DM_max)
            END IF
            eps_temp = DM_max*eps_temp
          else
            DM_max = 1.0d0
          end if

          if ( eps_temp .gt. hfx_opts%eps_schwarz ) then

            tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
            tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
            tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
            tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

            max_contraction_val = igto%orbnl_contraction_coeff(ishell) * &
                                  jgto%orbnl_contraction_coeff(jshell)* &
                                  kgto%orbnl_contraction_coeff(kshell)* &
                                  lgto%orbnl_contraction_coeff(lshell)* &
                                  DM_max

            call re_alloc(eri_deriv, 1, nso(l_i), 1, nso(l_j), 1, nso(l_k), &
              1, nso(l_l), 1, 12, name="eri_deriv", &
              routine="evaluate_gradient")
            eri_deriv(:,:,:,:,:) = 0.0_dp

            call calc_contract_deriv_eri(deriv_data, &
              hfx_sys%cell, hfx_sys%cell_r, &
              ri, rj, rk, rl, npgfi, npgfj, npgfk, npgfl, &
              l_i, l_j, l_k, l_l, ncoi, ncoj, ncok, ncol, &
              igto%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
              jgto%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
              kgto%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
              lgto%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
              igto%sphi(1:ncoi, ioa:ioa+nso(l_i)-1),        &
              jgto%sphi(1:ncoj, joa:joa+nso(l_j)-1),        &
              kgto%sphi(1:ncok, koa:koa+nso(l_k)-1),        &
              lgto%sphi(1:ncol, loa:loa+nso(l_l)-1),        &
              hfx_opts, eri_deriv, max_contraction_val, &
              max_val2_set, log10_pmax, tmp_R_1, tmp_R_2, &
              tmp_screen_pgf1, tmp_screen_pgf2)

            do nsoi = 1, nso(l_i)
              do nsoj = 1, nso(l_j)
                do nsok = 1, nso(l_k)
                  do nsol = 1, nso(l_l)

                    if ( DM_max*maxval(abs( &
                           eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2)) &
                             .gt. hfx_opts%eps_stored ) then
                      num_a = io + nsoi-1
                      num_b = jo + nsoj-1
                      num_c = ko + nsok-1
                      num_d = lo + nsol-1
                      nao_eri_deriv(1:12) =  &
                        eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2

                      call hfx_gradient_matrix(hfx_sys, num_a, num_b, &
                        num_c, num_d, maxnh, numh, listhptr, listh, &
                        nao_eri_deriv, Dscf, Fal)
                    endif

                  enddo
                enddo
              enddo
            enddo

            call de_alloc(eri_deriv, name="eri_deriv", &
              routine="evaluate_gradient")

          endif ! Schwarz inequality
        endif
      enddo   ! i_list_kl
    enddo   ! i_list_ij

   endif

  else


#ifdef MPI
    call GetNodeOrbs(list_kl%nelement, Node, Nodes, list_kl_local)
#endif

    do i_list_ij = 1, list_ij%nelement

      io = list_ij%element(i_list_ij)%pair(1)
      jo = list_ij%element(i_list_ij)%pair(2)
      ri = list_ij%element(i_list_ij)%r1
      rj = list_ij%element(i_list_ij)%r2
      rij2 = list_ij%element(i_list_ij)%dist2

      ia = hfx_sys%iaorb(io)
      ja = hfx_sys%iaorb(jo)

      is = hfx_sys%isa(ia)
      ispp => species(is)
      igto => hfx_gtos(is)
      ioa = hfx_sys%iphorb(io)
      l_i = lofio(is, ioa)
      npgfi = igto%orbnl_contract(ispp%orb_index(ioa))
      ncoi = nco(l_i)*npgfi

      js = hfx_sys%isa(ja)
      jspp => species(js)
      jgto => hfx_gtos(js)
      joa = hfx_sys%iphorb(jo)
      l_j = lofio(js, joa)
      npgfj = jgto%orbnl_contract(jspp%orb_index(joa))
      ncoj = nco(l_j)*npgfj

      ishell = ispp%orb_index(ioa)
      jshell = jspp%orb_index(joa)
      index_ij = list_ij%element(i_list_ij)%nl_index
      max_val1 = sfc_shell(jshell,ishell,js,is)%x(1)*rij2 + &
                 sfc_shell(jshell,ishell,js,is)%x(2)


#ifdef MPI
      do i_list_kl_local = 1, list_kl_local

        call LocalToGlobalOrb(i_list_kl_local, Node, Nodes, i_list_kl)
#else
      do i_list_kl = 1, list_kl%nelement
#endif
        ko = list_kl%element(i_list_kl)%pair(1)
        lo = list_kl%element(i_list_kl)%pair(2)
        rk = list_kl%element(i_list_kl)%r1
        rl = list_kl%element(i_list_kl)%r2
        rkl2 =  list_kl%element(i_list_kl)%dist2

        ka = hfx_sys%iaorb(ko)
        la = hfx_sys%iaorb(lo)

        ks = hfx_sys%isa(ka)
        kspp => species(ks)
        kgto => hfx_gtos(ks)
        koa = hfx_sys%iphorb(ko)
        l_k = lofio(ks, koa)
        npgfk = kgto%orbnl_contract(kspp%orb_index(koa))
        ncok = nco(l_k)*npgfk

        ls = hfx_sys%isa(la)
        lspp => species(ls)
        lgto => hfx_gtos(ls)
        loa = hfx_sys%iphorb(lo)
        l_l = lofio(ls, loa)
        npgfl = lgto%orbnl_contract(lspp%orb_index(loa))
        ncol = nco(l_l)*npgfl

        kshell = kspp%orb_index(koa)
        lshell = lspp%orb_index(loa)

        index_kl = list_kl%element(i_list_kl)%nl_index

        if ( index_kl .le. index_ij ) then

            if((D2Sindx(io,ko).eq.0).and.(D2Sindx(io,lo).eq.0).and.(D2Sindx(jo,ko).eq.0) &
              .and.(D2Sindx(jo,lo).eq.0) )  cycle

          ! Four centers along to the same atom
          if (ia .eq. ja .and. ia .eq. ka .and. ka .eq. la ) cycle

          max_val2_set = (sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
                          sfc_shell(lshell,kshell,ls,ks)%x(2) )

          max_val = max_val1 + max_val2_set

          eps_temp = eri_prescreen(D2Sindx(io,jo))*eri_prescreen(D2Sindx(ko,lo))
          eps_temp = dsqrt(eps_temp)

          if ( hfx_opts%DM_trunc ) then
             DM_max = 2.0_dp*max( Dscf_max(D2Sindx(io, ko))*Dscf_max(D2Sindx(jo,lo)), &
                        Dscf_max(D2Sindx(io, lo))*Dscf_max(D2Sindx(jo, ko)) )

            IF ( DM_max <= 0.0_dp ) THEN
              log10_pmax = log_zero
            ELSE
              log10_pmax = LOG10(DM_max)
            END IF
            eps_temp = DM_max*eps_temp
          else
            DM_max = 1.0d0
          end if

          if ( eps_temp .gt. hfx_opts%eps_schwarz ) then

            tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
            tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
            tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
            tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

            max_contraction_val = igto%orbnl_contraction_coeff(ishell) * &
                                  jgto%orbnl_contraction_coeff(jshell)* &
                                  kgto%orbnl_contraction_coeff(kshell)* &
                                  lgto%orbnl_contraction_coeff(lshell)* &
                                  DM_max

            call re_alloc(eri_deriv, 1, nso(l_i), 1, nso(l_j), 1, nso(l_k), &
              1, nso(l_l), 1, 12, name="eri_deriv", &
              routine="evaluate_gradient")
            eri_deriv(:,:,:,:,:) = 0.0_dp

            call calc_contract_deriv_eri(deriv_data, &
              hfx_sys%cell, hfx_sys%cell_r, &
              ri, rj, rk, rl, npgfi, npgfj, npgfk, npgfl, &
              l_i, l_j, l_k, l_l, ncoi, ncoj, ncok, ncol, &
              igto%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
              jgto%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
              kgto%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
              lgto%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
              igto%sphi(1:ncoi, ioa:ioa+nso(l_i)-1),        &
              jgto%sphi(1:ncoj, joa:joa+nso(l_j)-1),        &
              kgto%sphi(1:ncok, koa:koa+nso(l_k)-1),        &
              lgto%sphi(1:ncol, loa:loa+nso(l_l)-1),        &
              hfx_opts, eri_deriv, max_contraction_val, &
              max_val2_set, log10_pmax, tmp_R_1, tmp_R_2, &
              tmp_screen_pgf1, tmp_screen_pgf2)

            do nsoi = 1, nso(l_i)
              do nsoj = 1, nso(l_j)
                do nsok = 1, nso(l_k)
                  do nsol = 1, nso(l_l)

                    if ( DM_max*maxval(abs( &
                           eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2)) &
                             .gt. hfx_opts%eps_stored ) then
                      num_a = io + nsoi-1
                      num_b = jo + nsoj-1
                      num_c = ko + nsok-1
                      num_d = lo + nsol-1
                      nao_eri_deriv(1:12) =  &
                        eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2

                      call hfx_gradient_matrix(hfx_sys, num_a, num_b, &
                        num_c, num_d, maxnh, numh, listhptr, listh, &
                        nao_eri_deriv, Dscf, Fal)
                    endif

                  enddo
                enddo
              enddo
            enddo

            call de_alloc(eri_deriv, name="eri_deriv", &
              routine="evaluate_gradient")

          endif ! Schwarz inequality
        endif
      enddo   ! i_list_kl
    enddo   ! i_list_ij

endif

  end subroutine build_hfx_gradient

  ! ***************************************************************************
  !> \brief To be defined
  ! ***************************************************************************
  subroutine hfx_gradient_matrix(hfx_sys, io, jo, ko, lo, &
                maxnh, numh, listhptr, listh, nao_eri_deriv, Dscf, Fal)

    use nao2gto_data, only: subshell, D2Sindx

    implicit none

    ! Arguments
    type(hfx_system_type), intent(in) :: hfx_sys
    integer, intent(in) :: io, jo, ko, lo
    integer, intent(in) :: maxnh, listh(maxnh), listhptr(hfx_sys%nuotot), numh(hfx_sys%nuotot)
    real(dp), intent(in) :: nao_eri_deriv(12)
    real(dp), intent(in) :: Dscf(maxnh,hfx_sys%nspin)
    real(dp), intent(inout) :: Fal(1:3,hfx_sys%nua)

    ! Local variables
    integer :: ispin, ncells
    integer :: iuo, juo, kuo, luo, llo, ishell, jshell, kshell, lshell, &
      iushell, jushell, kushell, lushell, index_ij, index_kl, &
      io_trans, jo_trans, ko_trans, lo_trans, ia, ja, ka, la
    real(dp) :: gint_deriv(1:12), spin_factor

    ! -------------------------------------------------------------------------

    if(nspin.eq.1) then
       spin_factor = 1.0d0
    else
       spin_factor = 2.0d0
    endif

    ncells = hfx_sys%norb / hfx_sys%nuotot

    gint_deriv = nao_eri_deriv
    iuo = io ! u is always u0
    juo = indxuo(jo)
    kuo = indxuo(ko)
    luo = indxuo(lo)
    llo = indexsc(ko, kuo, lo)
    ! lo have to trans to play with m0 to get kl and
    ! compared to ij, so there is llo

    ia = hfx_sys%iaorb(iuo)
    ja = hfx_sys%iaorb(juo)
    ka = hfx_sys%iaorb(kuo)
    la = hfx_sys%iaorb(luo)

    iushell = subshell(iuo)
    jushell = subshell(juo)
    kushell = subshell(kuo)
    lushell = subshell(luo)
    jshell  = subshell(jo)
    lshell  = subshell(llo)

    index_ij = ncells*iushell*(iushell-1)/2 + &
               ((jshell-1)/subshell(hfx_sys%nuotot))*iushell + jushell

    index_kl = ncells*kushell*(kushell-1)/2 + &
               ((lshell-1)/subshell(hfx_sys%nuotot))*kushell + lushell

    if ( iushell .eq. jushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
    if ( kushell .eq. lushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
    if ( index_ij .eq. index_kl ) gint_deriv(1:12) = gint_deriv(1:12)*0.5d0

    ! HFX
    !
    !         (u0vR|mR'nR")        =       (u0vR|nR"mR')
    ! = (v0u[-R]|m[R'-R]n[R"-R])   = (v0u[-R]|n[R"-R]m[R'-R])
    ! = (m0n[R"-R']|u[-R']v[R-R']) = (m0n[R"-R']|v[R-R']u[-R'])
    ! = (n0m[R'-R"]|u[-R"]v[R-R"]) = (n0m[R'-R"]|v[R-R"]u[-R"])
    do ispin=1,hfx_sys%nspin

      ! 1.VEE(1[0]  2[H] | 3[G]  4[N])  (u0v[R]|m[R']n[R"])
      !   VEE(1[0]  2[H] | 4[N]  3[G])  (u0v[R])|n[R"]m[R'])
      if((D2Sindx(jo,lo).ne.0) .and. (D2Sindx(io,ko).ne.0)) then
      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(jo,lo),ispin)*Dscf(D2Sindx(io,ko),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(jo,lo),ispin)*Dscf(D2Sindx(io,ko),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(jo,lo),ispin)*Dscf(D2Sindx(io,ko),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(jo,lo),ispin)*Dscf(D2Sindx(io,ko),ispin)
      endif

      if((D2Sindx(jo,ko).ne.0) .and. (D2Sindx(io,lo).ne.0)) then
      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(jo,ko),ispin)*Dscf(D2Sindx(io,lo),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(jo,ko),ispin)*Dscf(D2Sindx(io,lo),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(jo,ko),ispin)*Dscf(D2Sindx(io,lo),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(jo,ko),ispin)*Dscf(D2Sindx(io,lo),ispin)
      endif


      ! 2.VEE(2[0] 1[-H]| 3[G-H] 4[N-H])  (v0u[-R]|m[R'-R]n[R"-R])
      !   VEE(2[0] 1[-H]| 4[N-H] 3[G-H])  (v0u[-R]|n[R"-R]m[R'-R])

         io_trans = indexsc( jo, juo, io )
         ko_trans = indexsc( jo, juo, ko )
         lo_trans = indexsc( jo, juo, lo )

      if((D2Sindx(io_trans,lo_trans).ne.0) .and. (D2Sindx(juo,ko_trans).ne.0)) then
      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(io_trans,lo_trans),ispin)*Dscf(D2Sindx(juo,ko_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(io_trans,lo_trans),ispin)*Dscf(D2Sindx(juo,ko_trans),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(io_trans,lo_trans),ispin)*Dscf(D2Sindx(juo,ko_trans),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(io_trans,lo_trans),ispin)*Dscf(D2Sindx(juo,ko_trans),ispin)
      endif

      if((D2Sindx(io_trans,ko_trans).ne.0) .and.(D2Sindx(juo,lo_trans).ne.0)) then
      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(io_trans,ko_trans),ispin)*Dscf(D2Sindx(juo,lo_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(io_trans,ko_trans),ispin)*Dscf(D2Sindx(juo,lo_trans),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(io_trans,ko_trans),ispin)*Dscf(D2Sindx(juo,lo_trans),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(io_trans,ko_trans),ispin)*Dscf(D2Sindx(juo,lo_trans),ispin)
      endif

      ! 3.VEE(3[0]  4[N-G] | 1[-G] 2[H-G]) (m0n[R"-R']|u[-R']v[R-R'])
      !   VEE(3[0]  4[N-G] |2[H-G] 1[-G] ) (m0n[R"-R']|v[R-R']u[-R'])

         io_trans = indexsc( ko, kuo, io )
         jo_trans = indexsc( ko, kuo, jo )
         lo_trans = indexsc( ko, kuo, lo )

      if((D2Sindx(lo_trans,jo_trans).ne.0) .and. (D2Sindx(kuo,io_trans).ne.0))then
      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin)
      endif

      if((D2Sindx(lo_trans,io_trans).ne.0) .and. (D2Sindx(kuo,jo_trans).ne.0))then
      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(lo_trans,io_trans),ispin)*Dscf(D2Sindx(kuo,jo_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(lo_trans,io_trans),ispin)*Dscf(D2Sindx(kuo,jo_trans),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(lo_trans,io_trans),ispin)*Dscf(D2Sindx(kuo,jo_trans),ispin)

      Fal(1:3,la) = Fal(1:3,la) &
        +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(lo_trans,io_trans),ispin)*Dscf(D2Sindx(kuo,jo_trans),ispin)
      endif

! 4.VEE(4[0]  3[G-N] | 1[-N] 2[H-N])  (n0m[R'-R"]|u[-R"]v[R-R"])
!   VEE(4[0]  3[G-N] | 2[H-N] 1[-N])  (n0m[R'-R"]|v[R-R"]u[-R"])

         io_trans = indexsc( lo, luo, io )
         jo_trans = indexsc( lo, luo, jo )
         ko_trans = indexsc( lo, luo, ko )

      if((D2Sindx(ko_trans,jo_trans).ne.0) .and. (D2Sindx(luo,io_trans).ne.0)) then
      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(ko_trans,jo_trans),ispin)*Dscf(D2Sindx(luo,io_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(ko_trans,jo_trans),ispin)*Dscf(D2Sindx(luo,io_trans),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(ko_trans,jo_trans),ispin)*Dscf(D2Sindx(luo,io_trans),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(ko_trans,jo_trans),ispin)*Dscf(D2Sindx(luo,io_trans),ispin)

      endif

      if((D2Sindx(ko_trans,io_trans).ne.0) .and. (D2Sindx(luo,jo_trans).ne.0)) then
      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(ko_trans,io_trans),ispin)*Dscf(D2Sindx(luo,jo_trans),ispin)

      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(ko_trans,io_trans),ispin)*Dscf(D2Sindx(luo,jo_trans),ispin)

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(ko_trans,io_trans),ispin)*Dscf(D2Sindx(luo,jo_trans),ispin)

      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(ko_trans,io_trans),ispin)*Dscf(D2Sindx(luo,jo_trans),ispin)
      endif

    enddo

  end subroutine hfx_gradient_matrix
  ! ***** GRADIENT MODULE END *****

end module nao2gto_hfx_force
