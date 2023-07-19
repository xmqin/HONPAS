! Writen by xmqin, October 2013
      module nao2gto_gradient
      use kinds,          only : dp, int_8
      use parallel,       only : Node, Nodes
      use parallelsubs,   only : GetNodeOrbs, GlobalToLocalOrb, &
                                 LocalToGlobalOrb,WhichNodeOrb, &
                                 GetNodeOrbs_NAO2GTO, &
                                 LocalToGlobalOrb_NAO2GTO
#ifdef MPI
      use mpi_siesta
#endif
      use atm_types,      only : species, species_info, l_max, nco, nso
      use hfx_types,      only : ERI_link, head, tail, ptr
      use hfx_types,      only : hfx_input_parameter, subshell
      use hfx_types,      only : pair_dist_radii_pgf, sfc_pgf, sfc_shell
      use atomlist,       only : indxuo
      use hfx_types,      only : list_ij, list_kl
      use hfx_types,      only : hfx_screen_coeff_type
      use hfx_types,      only : log_zero, D2Sindx
      use hfx_types,      only : eri_prescreen
      use atmfuncs,       only : lofio, mofio
      use listsc_module,  only : listsc
      use extended_index, only : indexsc 
      use alloc,          only : re_alloc, de_alloc
      use gamma_F,       only : init_md_ftable
      use contract_eri,  only : calc_contract_deriv_eri
      use libint_wrapper_types,  only : lib_int, lib_deriv
!      use get_DM,        only : dense2sparse

      implicit none 
      private

      public evaluate_gradient,  hfx_gradient_matrix

      contains

      subroutine evaluate_gradient( deriv, nspin, norb, iaorb, iphorb, nuotot,  &
                              nua, na, isa, maxnh, numh, listhptr, listh, &
                              cell, rcell,  &
                              hfx_parameter, Dscf, Dscf_max,  &
!                              radii_pgf, screen_coeffs_pgf, screen_coeffs_shell,  &
                              log10_eps_schwarz, Fal )


! ----------------------------INPUT--------------------------------
      !type(lib_int)        :: lib
      type(lib_deriv)      :: deriv
      integer,  intent(in) :: &
        maxnh, nspin, nua, na, norb, nuotot,  &
        iaorb(norb), iphorb(norb), isa(na), &
        listh(maxnh), listhptr(nuotot),  numh(nuotot)


      type(hfx_input_parameter) :: hfx_parameter

      real(dp), intent(in) :: Dscf(maxnh,nspin)
      real(dp), intent(in) :: Dscf_max(maxnh)
!      logical,  intent(in) :: um_cut(norb,norb)

      real(dp), intent(in) :: cell(3,3), rcell(3,3)
      real(dp), intent(inout) :: Fal(3, nua)


      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:,:,:,:,:), POINTER        :: screen_coeffs_pgf, radii_pgf
      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:,:,:), POINTER            :: screen_coeffs_shell

      real(dp), intent(in)   :: log10_eps_schwarz

!------------------------------------------------------------------
!------------------------------------------------------------------
      type(species_info), pointer :: ispp, jspp, kspp, lspp
      real(dp)  ri(3), rj(3), rk(3), rl(3), rij2,  rkl2

      real(dp), dimension(:,:,:,:,:), allocatable :: eri_deriv
      real(dp) eps_temp, DM_max,symm_factor

      real(dp) max_contraction_val, max_val, max_val1, &
               max_val2, pmax_entry

      TYPE(hfx_screen_coeff_type), &
        DIMENSION(:,:), POINTER                :: tmp_R_1, tmp_R_2, &
                                     tmp_screen_pgf1, tmp_screen_pgf2
      REAL(dp)       :: max_val2_set, log10_pmax

      
      real(dp) :: nao_eri_deriv(12)
      integer io, jo, ko, lo, is, js, ks, ls, ioa, joa, koa, loa,   &
              l_i, l_j, l_k, l_l, m_i, m_j, m_k, m_l,               &
              ncoi, ncoj, ncok, ncol, npgfi, npgfj, npgfk, npgfl,   &
              nsoi, nsoj, nsok, nsol, num_a, num_b, num_c, num_d,   &
              i_list_ij, i_list_kl, i_list_kl_local, list_kl_local, &
              index_ij, index_kl, ia, ja, ka, la , ishell, jshell,  &
              kshell, lshell

              
      integer ispin, ncells

      real(dp) time_start, time_end
!               time_start1, time_end1

      integer :: ind
!      logical :: find_in_Dscf
      logical :: um_cut_1, um_cut_2, um_cut_3, um_cut_4
      real(dp) :: eri_prescreen_value_1, eri_prescreen_value_2

      integer(int_8) shell_eri, number_task

#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   num_loc
#endif
      logical  :: flag
      integer  :: list_pos(2), temp_list_pos(2)
      integer  :: buf, i, i_list_ijkl, temp_i_list_kl, tag, &
                  task_frag, &
                  count, temp_count

!      tot_time = 0.0_dp
!      eri_time = 0.0_dp
!      hfx_time = 0.0_dp
     
       call cpu_time(time_start)

!      shell_eri = 0_int_8
      ncells = norb / nuotot

      if(hfx_parameter%parallel) then

      buf  =  11
      request = 1000
      tag = 111
      task_frag = hfx_parameter%fragsize

      if( Node .eq. 0 ) then
        number_task = list_ij%nelement*list_kl%nelement

        write(6,*) list_ij%nelement, list_kl%nelement, number_task
        call MPI_IRECV(buf, 1, MPI_integer, MPI_ANY_SOURCE, MPI_ANY_TAG,&
                       MPI_COMM_WORLD, request, MPIerror)
        list_pos(1) = 1
        list_pos(2) = 1
        do while(number_task > 0)
            CALL MPI_TEST(request, flag, status, MPIerror)
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
            CALL MPI_RECV(count, 1, MPI_integer, 0, tag, MPI_COMM_WORLD,status, MPIerror)
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

         ia = iaorb(io)
         ja = iaorb(jo)

         is = isa(iaorb(io))
         ispp => species(is)
         ioa = iphorb(io)
         l_i = lofio(is, ioa)
         m_i = mofio(is, ioa)
         npgfi = ispp%orbnl_contract(ispp%orb_index(ioa))
         ncoi = nco(l_i)*npgfi

         js = isa(iaorb(jo))
         jspp => species(js)
         joa = iphorb(jo)
         l_j = lofio(js, joa)
         m_j = mofio(js, joa)
         npgfj = jspp%orbnl_contract(jspp%orb_index(joa))
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
            ka = iaorb(ko) 
            la = iaorb(lo)

            ks = isa(iaorb(ko))
            kspp => species(ks)
            koa = iphorb(ko)
            l_k = lofio(ks, koa)
            m_k = mofio(ks, koa)
            npgfk = kspp%orbnl_contract(kspp%orb_index(koa))
            ncok =nco(l_k)*npgfk

            ls = isa(iaorb(lo))
            lspp => species(ls)
            loa = iphorb(lo)
            l_l = lofio(ls, loa)
            m_l = mofio(ls, loa)
            npgfl = lspp%orbnl_contract(lspp%orb_index(loa))
            ncol = nco(l_l)*npgfl

            kshell = kspp%orb_index(koa)
            lshell = lspp%orb_index(loa)
            
            index_kl = list_kl%element(i_list_kl)%nl_index
           endif

           if(ia.eq.ja.and.ia.eq.ka.and.ka.eq.la) cycle             
            if( index_kl .le. index_ij ) then

            if((D2Sindx(io,ko).eq.0).and.(D2Sindx(io,lo).eq.0).and.(D2Sindx(jo,ko).eq.0) &
              .and.(D2Sindx(jo,lo).eq.0) )  cycle

            max_val2_set = (sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
                        sfc_shell(lshell,kshell,ls,ks)%x(2) )

            max_val = max_val1 + max_val2_set

            eps_temp = eri_prescreen(D2Sindx(io, jo))*eri_prescreen(D2Sindx(ko,lo))
            eps_temp = dsqrt(eps_temp)
!            if(node.eq.0) write(60,*) io,jo,ko,lo,eps_temp
            if(hfx_parameter%DM_trunc) then
               DM_max = 2.0_dp*max( Dscf_max(D2Sindx(io, ko))*Dscf_max(D2Sindx(jo, lo)), &
                        Dscf_max(D2Sindx(io, lo))*Dscf_max(D2Sindx(jo, ko)) )
               eps_temp = DM_max*eps_temp
             else
               DM_max = 1.0_dp
            endif


              IF(DM_max <= 0.0_dp) THEN
                  log10_pmax = log_zero
              ELSE
                  log10_pmax = LOG10(DM_max)
              END IF

               if( eps_temp .ge.hfx_parameter%eps_schwarz) then  ! hfx_parameter%eps_schwarz) then
    
                  tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
                  tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
                  tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
                  tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

                  max_contraction_val = ispp%orbnl_contraction_coeff(ishell) * &
                                        jspp%orbnl_contraction_coeff(jshell)* &
                                        kspp%orbnl_contraction_coeff(kshell)* &
                                        lspp%orbnl_contraction_coeff(lshell)*DM_max

                  allocate( eri_deriv(nso(l_i),nso(l_j),nso(l_k),nso(l_l),12) )
                  eri_deriv = 0.0_dp

!                  call cpu_time(time_start)

                  call calc_contract_deriv_eri( deriv, cell, rcell, ri, rj, rk, rl,   &
                                            npgfi, npgfj, npgfk, npgfl,  &
                                            l_i, l_j, l_k, l_l,          &
                                            ncoi, ncoj, ncok, ncol,      &
                           ispp%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
                           jspp%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
                           kspp%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
                           lspp%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
                           ispp%sphi(1:ncoi, ioa:ioa+nso(l_i)-1),        &
                           jspp%sphi(1:ncoj, joa:joa+nso(l_j)-1),        &
                           kspp%sphi(1:ncok, koa:koa+nso(l_k)-1),        &
                           lspp%sphi(1:ncol, loa:loa+nso(l_l)-1),        &
                           hfx_parameter, eri_deriv,max_contraction_val, &
                           max_val2_set,log10_eps_schwarz, log10_pmax,   &
                           tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2 )
!                   call cpu_time(time_end)
!                   eri_time = eri_time + time_end-time_start

                do nsoi = 1, nso(l_i)
                   do nsoj = 1, nso(l_j)
                      do nsok = 1, nso(l_k)
                         do nsol = 1, nso(l_l)                   

                          if(DM_max*maxval(dabs(eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2.0_dp)) &
                              .ge.hfx_parameter%eps_stored ) then 
                               num_a = io + nsoi-1
                               num_b = jo + nsoj-1
                               num_c = ko + nsok-1
                               num_d = lo + nsol-1
                               nao_eri_deriv(1:12) =  &
                                  eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2.0_dp
                               
!                               call cpu_time(time_start)
                               call  hfx_gradient_matrix( nspin, ncells, nuotot, norb, &
                                               iaorb, nua, num_a, num_b, num_c, num_d, &
                                               maxnh, numh, listhptr, listh, &
                                               nao_eri_deriv, Dscf, Fal )
!                               call cpu_time(time_end)
!                               hfx_time = hfx_time + time_end-time_start
!                               shell_eri = shell_eri + 1
                             endif

                         enddo
                      enddo
                   enddo
                enddo

                deallocate(eri_deriv)

                endif ! Schwarz inequality
            endif
         enddo !mn
       enddo   !uv
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

         ia = iaorb(io)
         ja = iaorb(jo)

         is = isa(iaorb(io))
         ispp => species(is)
         ioa = iphorb(io)
         l_i = lofio(is, ioa)
         m_i = mofio(is, ioa)
         npgfi = ispp%orbnl_contract(ispp%orb_index(ioa))
         ncoi = nco(l_i)*npgfi

         js = isa(iaorb(jo))
         jspp => species(js)
         joa = iphorb(jo)
         l_j = lofio(js, joa)
         m_j = mofio(js, joa)
         npgfj = jspp%orbnl_contract(jspp%orb_index(joa))
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
            ka = iaorb(ko) 
            la = iaorb(lo)

            ks = isa(iaorb(ko))
            kspp => species(ks)
            koa = iphorb(ko)
            l_k = lofio(ks, koa)
            m_k = mofio(ks, koa)
            npgfk = kspp%orbnl_contract(kspp%orb_index(koa))
            ncok =nco(l_k)*npgfk

            ls = isa(iaorb(lo))
            lspp => species(ls)
            loa = iphorb(lo)
            l_l = lofio(ls, loa)
            m_l = mofio(ls, loa)
            npgfl = lspp%orbnl_contract(lspp%orb_index(loa))
            ncol = nco(l_l)*npgfl

            kshell = kspp%orb_index(koa)
            lshell = lspp%orb_index(loa)
            
            index_kl = list_kl%element(i_list_kl)%nl_index

            if(ia.eq.ja.and.ia.eq.ka.and.ka.eq.la) cycle 
            if( index_kl .le. index_ij ) then

            if((D2Sindx(io,ko).eq.0).and.(D2Sindx(io,lo).eq.0).and.(D2Sindx(jo,ko).eq.0) &
              .and.(D2Sindx(jo,lo).eq.0) )  cycle

            max_val2_set = (sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
                        sfc_shell(lshell,kshell,ls,ks)%x(2) )

            max_val = max_val1 + max_val2_set

            eps_temp = eri_prescreen(D2Sindx(io, jo))*eri_prescreen(D2Sindx(ko,lo))
            eps_temp = dsqrt(eps_temp)

            if(hfx_parameter%DM_trunc) then
               DM_max = 2.0_dp*max( Dscf_max(D2Sindx(io, ko))*Dscf_max(D2Sindx(jo, lo)), &
                        Dscf_max(D2Sindx(io, lo))*Dscf_max(D2Sindx(jo, ko)) )

               eps_temp = DM_max*eps_temp

             else
               DM_max = 1.0_dp

            endif

              IF(DM_max <= 0.0_dp) THEN
                  log10_pmax = log_zero
              ELSE
                  log10_pmax = LOG10(DM_max)
              END IF

               if( eps_temp .ge.hfx_parameter%eps_schwarz) then  ! hfx_parameter%eps_schwarz) then
    
                  tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
                  tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
                  tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
                  tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

                  max_contraction_val = ispp%orbnl_contraction_coeff(ishell) * &
                                        jspp%orbnl_contraction_coeff(jshell)* &
                                        kspp%orbnl_contraction_coeff(kshell)* &
                                        lspp%orbnl_contraction_coeff(lshell)* &
                                        DM_max


                   allocate( eri_deriv(nso(l_i),nso(l_j),nso(l_k),nso(l_l),12) )

                   eri_deriv = 0.0_dp
!                   call cpu_time(time_start)
                   call calc_contract_deriv_eri( deriv, cell, rcell, ri, rj, rk, rl,   &
                                            npgfi, npgfj, npgfk, npgfl,  &
                                            l_i, l_j, l_k, l_l,          &
                                            ncoi, ncoj, ncok, ncol,      &
                           ispp%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
                           jspp%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
                           kspp%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
                           lspp%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
                           ispp%sphi(1:ncoi, ioa:ioa+nso(l_i)-1),        &
                           jspp%sphi(1:ncoj, joa:joa+nso(l_j)-1),        &
                           kspp%sphi(1:ncok, koa:koa+nso(l_k)-1),        &
                           lspp%sphi(1:ncol, loa:loa+nso(l_l)-1),        &
                           hfx_parameter, eri_deriv,max_contraction_val, &
                           max_val2_set,log10_eps_schwarz, log10_pmax,   &
                           tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2 )
!                   call cpu_time(time_end)
!                   eri_time = eri_time + time_end-time_start
                    
                do nsoi = 1, nso(l_i)
                   do nsoj = 1, nso(l_j)
                      do nsok = 1, nso(l_k)
                         do nsol = 1, nso(l_l)                   

                          if(DM_max*maxval(dabs(eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2)) &
                              .ge.hfx_parameter%eps_stored ) then 
                               num_a = io + nsoi-1
                               num_b = jo + nsoj-1
                               num_c = ko + nsok-1
                               num_d = lo + nsol-1
                               nao_eri_deriv(1:12) =  &
                                  eri_deriv(nsoi,nsoj,nsok,nsol,1:12)*2 
                               
!                               call cpu_time(time_start)
                               call  hfx_gradient_matrix( nspin, ncells, nuotot, norb, &
                                               iaorb, nua, num_a, num_b, num_c, num_d, &
                                               maxnh, numh, listhptr, listh, &
                                               nao_eri_deriv, Dscf, Fal )
!                               call cpu_time(time_end)
!                               hfx_time = hfx_time + time_end-time_start
!                               shell_eri = shell_eri + 1
                             endif

                         enddo
                      enddo
                   enddo
                enddo

                deallocate(eri_deriv)

                endif ! Schwarz inequality
            endif
         enddo !mn
       enddo   !uv

      endif

       call cpu_time(time_end)
      if(node .eq.0) then
         write(6,*) "HFX Force time : ", time_end-time_start, " s "
      endif

!       tot_time = time_end1-time_start1
!             call cpu_time(time_start)

!       write(6,'(a)')"--------------PrintHFX_Time_Information-----------------"
!       write(6,'(a21,2x,f12.6,2x,a3)') " ERI calculate time =", eri_time, "[s]"
!       write(6,'(a21,2x,f12.6,2x,a3)') " HFX Force time     =", hfx_time, "[s]"
!       write(6,'(a21,2x,f12.6,2x,a3)') " Total Force time   =",tot_time, "[s]"
!       write(6,'(a)')"----------------------------------------------------"
!
      end subroutine evaluate_gradient


      subroutine hfx_gradient_matrix( nspin, ncells, nuotot, norb, &
                                      iaorb, nua, io, jo, ko, lo,  &
                                      maxnh, numh, listhptr, listh, &
                                      nao_eri_deriv, Dscf, Fal )

!---------------------------- INPUT  VARIABLES -------------------------------------
      integer, intent(in)     :: nspin, ncells, nuotot, norb, iaorb(norb), &
                                 nua, io, jo, ko, lo, &
                                 maxnh, listh(maxnh), listhptr(nuotot), numh(nuotot)

      real(dp), intent(in)    :: nao_eri_deriv(12)
      real(dp), intent(in)    :: Dscf(maxnh,nspin)
!--------------------------- INOUT VARIABLES ---------------------------------------
      real(dp), intent(inout) :: Fal(1:3,nua)

!------------------------- TEMPOS, INTERNAL VARIABLES ------------------------------
      real(dp) gint_deriv(1:12), spin_factor
      integer iuo, juo, kuo, luo, llo,  &
              ishell, jshell, kshell, lshell, iushell,  &
              jushell, kushell, lushell, index_ij, index_kl, &
              io_trans, jo_trans, ko_trans, lo_trans, &
              ia, ja, ka, la

      integer ispin
      integer ind_1, ind_2, ind_3, ind_4
      logical :: find_in_1, find_in_2, find_in_3, find_in_4

      if(nspin.eq.1) then
         spin_factor = 1.0d0
      else
         spin_factor = 2.0d0
      endif

            gint_deriv = nao_eri_deriv
            iuo = io ! u is always u0
            juo = indxuo(jo)
            kuo = indxuo(ko)
            luo = indxuo(lo)
            llo = indexsc(ko, kuo, lo)
            ! lo have to trans to play with m0 to get kl and
            ! campared to ij, so there is llo

            ia=iaorb(iuo)
            ja=iaorb(juo)
            ka=iaorb(kuo)
            la=iaorb(luo)
                                                 
            iushell = subshell(iuo)
            jushell = subshell(juo)
            kushell = subshell(kuo)
            lushell = subshell(luo)
            jshell  = subshell(jo)
            lshell  = subshell(llo)

            index_ij = ncells*iushell*(iushell-1)/2 + &
                       ((jshell-1)/subshell(nuotot))*iushell + jushell 

            index_kl = ncells*kushell*(kushell-1)/2 + &
                       ((lshell-1)/subshell(nuotot))*kushell + lushell

            if( iushell .eq. jushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
            if( kushell .eq. lushell )   gint_deriv(1:12) = gint_deriv(1:12)*0.5d0
            if( index_ij .eq. index_kl ) gint_deriv(1:12) = gint_deriv(1:12)*0.5d0

!! HFX
        do ispin=1,nspin
!
!         (u0vR|mR'nR")        =       (u0vR|nR"mR')
! = (v0u[-R]|m[R'-R]n[R"-R])   = (v0u[-R]|n[R"-R]m[R'-R])
! = (m0n[R"-R']|u[-R']v[R-R']) = (m0n[R"-R']|v[R-R']u[-R'])
! = (n0m[R'-R"]|u[-R"]v[R-R"]) = (n0m[R'-R"]|v[R-R"]u[-R"])
!
       ! 1.VEE(1[0]  2[H] | 3[G]  4[N])  (u0v[R]|m[R']n[R"])
       !   VEE(1[0]  2[H] | 4[N]  3[G])  (u0v[R])|n[R"]m[R'])

!      write(60,*) io, jo, ko, lo, Dscf(D2Sindx(jo,lo),ispin), Dscf(D2Sindx(io,ko),ispin), &
!                   Dscf(D2Sindx(jo,ko),ispin), Dscf(D2Sindx(io,lo),ispin)

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

      if((D2Sindx(lo_trans,jo_trans).ne.0) .and. (D2Sindx(kuo,io_trans).ne.0)) then
      Fal(1:3,ka) = Fal(1:3,ka) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(7:9)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin) 

      Fal(1:3,ia) = Fal(1:3,ia) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(1:3)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin) 

      Fal(1:3,ja) = Fal(1:3,ja) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(4:6)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin) 

      Fal(1:3,la) = Fal(1:3,la) &
       +0.25d0*0.25d0*spin_factor*gint_deriv(10:12)*Dscf(D2Sindx(lo_trans,jo_trans),ispin)*Dscf(D2Sindx(kuo,io_trans),ispin) 
      endif 

      if((D2Sindx(lo_trans,io_trans).ne.0) .and. (D2Sindx(kuo,jo_trans).ne.0)) then
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

      end module nao2gto_gradient
