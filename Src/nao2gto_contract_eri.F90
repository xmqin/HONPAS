! This file is part of HONPAS package!
! Coded by cp2k group
! Edited by shanghui, xmqin  
module contract_eri

    use kinds,                           only : dp, int_8
    use primitive_eri,                   only : calc_primitive_eri
    use primitive_eri,                   only : calc_primitive_eri2
    use primitive_eri,                   only : calc_primitive_deriv_eri
    use libint_wrapper_types,            only : lib_int, lib_deriv
    use atm_types,                       only : l_max, nso, nco
    use cell_pbc
    use hfx_types,                       only : hfx_input_parameter, &
                                                hfx_screen_coeff_type


    implicit none
    private
    public calc_contract_eri, calc_contract_eri2, calc_contract_deriv_eri 

  contains

    subroutine calc_contract_eri( lib, cell, rcell, ra, rb, rc, rd,  &
                                  npgfa, npgfb, npgfc, npgfd, &
                                  la, lb, lc, ld,             & 
                                  ncoa, ncob, ncoc, ncod,     & 
                                  zeta, zetb, zetc, zetd,     &
                                  sphia, sphib, sphic, sphid, &
                                  hfx_parameter, eri, neris_tmp, &
                                  max_contraction, max_val2_set,&
                                  log10_eps_schwarz, log10_pmax,R1_pgf, R2_pgf, pgf1, pgf2 ) 

    !input and output variables 
    type(lib_int)                            :: lib 
    integer npgfa,npgfb,npgfc,npgfd,la,lb,lc,ld,ncoa,ncob,ncoc,ncod
    type(hfx_input_parameter)  :: hfx_parameter    
    real(dp) ::  cell(3,3), rcell(3,3),                            &
                 ra(3),rb(3),rc(3),rd(3), rab(3), rcd(3),          &
                 zeta(npgfa),zetb(npgfb),zetc(npgfc),zetd(npgfd),  & 
                 sphia(ncoa,nso(la)),sphib(ncob,nso(lb)), &
                 sphic(ncoc,nso(lc)),sphid(ncod,nso(ld))           
    real(dp) ::  eri( nso(la),nso(lb),nso(lc),nso(ld) )

    !internal variables
    real(dp),dimension(:),allocatable  ::  primitive_integrals
    real(dp),dimension(:),allocatable  ::  T1
 
    integer ipgf,jpgf,kpgf,lpgf,offset_a,offset_b,offset_c,offset_d,  &
            am,bm,cm,dm,ieri,i,j,k,l,index_primitive_integrals

    integer(int_8)                     ::  neris_tmp
    REAL(dp), INTENT(IN)      ::  max_contraction
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: R1_pgf, R2_pgf, pgf1, pgf2

    REAL(dp)                                 :: max_val2_set, &
                                                log10_eps_schwarz, log10_pmax


    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, &
      RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
      Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv,rab2,rcd2,      &
      r_temp(3),r_pbc_temp(3) ,   &
      alpha_P,alpha_Q,R_P,R_Q,Kab,Kcd,theta_w, far_eri, &
      rc_trans(3),rd_trans(3) , R1, R2, pgf_max_1, pgf_max_2, cart_estimate

    real(dp), parameter  ::  pi=3.14159265358979323846264338_dp 
    real(dp), parameter  ::  erfc_epsilon = 0.99998871620832942100_dp      

    allocate( primitive_integrals(ncoa*ncob*ncoc*ncod) )    
    allocate( T1(ncoa*ncob*ncoc*ncod) ) 
    primitive_integrals(1:ncoa*ncob*ncoc*ncod)=0.0d0
    T1(1:ncoa*ncob*ncoc*ncod)=0.0d0

    neris_tmp = 0
    cart_estimate = 0.0_dp

    rab(1:3)=ra(1:3)-rb(1:3)
    rcd(1:3)=rc(1:3)-rd(1:3)
    rab2=rab(1)**2+rab(2)**2+rab(3)**2
    rcd2=rcd(1)**2+rcd(2)**2+rcd(3)**2

    DO ipgf = 1,npgfa
       offset_a = (ipgf-1)*nco(la)
       Zeta_A = zeta(ipgf)

      DO jpgf = 1,npgfb
         pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
!         write(111,*) pgf_max_1,max_val2_set,log10_pmax
         IF( pgf_max_1 + max_val2_set + log10_pmax < log10_eps_schwarz) CYCLE
!         R1 = MAX(0.0_dp, R1_pgf(jpgf,ipgf)%x(1)*rab2 + R1_pgf(jpgf,ipgf)%x(2))

         Zeta_B = zetb(jpgf)
         Zeta1 = Zeta_A + Zeta_B
         ZetaInv = 1.0d0/Zeta1
         S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
         P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv
         offset_b = (jpgf-1)*nco(lb)

         DO kpgf = 1,npgfc
            offset_c = (kpgf-1)*nco(lc)
            Zeta_C = zetc(kpgf)

            DO lpgf = 1,npgfd
               pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
               IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz ) CYCLE
!               R2 = MAX(0.0_dp, R2_pgf(lpgf,kpgf)%x(1)*rcd2 + R2_pgf(lpgf,kpgf)%x(2))

               Zeta_D = zetd(lpgf)
               Eta  =  Zeta_C + Zeta_D
               EtaInv = 1.0d0/Eta
               ZetapEtaInv = Zeta1+Eta
               ZetapEtaInv = 1.0d0/ZetapEtaInv
               Rho = Zeta1*Eta*ZetapEtaInv
               RhoInv = 1.0d0/Rho
               S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
               Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
               r_temp(1:3) = Q(1:3)-P(1:3)
               call trans_pbc(r_temp,cell,rcell,r_pbc_temp)

               Q(1:3)=P(1:3)+r_pbc_temp(1:3)
               rc_trans(1:3)=rc(1:3)+r_pbc_temp(1:3)-r_temp(1:3)
               rd_trans(1:3)=rd(1:3)+r_pbc_temp(1:3)-r_temp(1:3)            

               rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
               W(1:3) = (Zeta1*P(1:3)+Eta*Q(1:3))*ZetapEtaInv
               offset_d = (lpgf-1)*nco(ld)

               tmp_max = 0.0_dp

               if((hfx_parameter%potential_type.eq.1) .or. (.not.hfx_parameter%far_field)) then

                  call calc_primitive_eri(lib, ra, rb, rc_trans, rd_trans,&
                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                           la, lb, lc ,ld,&
                           ncoa,ncob,ncoc,ncod,&
                           offset_a,offset_b,offset_c,offset_d, &
                           primitive_integrals, hfx_parameter,  &
                           max_contraction, tmp_max, neris_tmp,     &
                           Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                           P,Q,W,rpq2,rab,rcd,R1,R2)

               else ! farfield screening for short-range primitive ERIs

                  alpha_P=Zeta_A*Zeta_B*ZetaInv
                  R_P=INT(1.0_dp/(SQRT(2.0_dp*alpha_P)*erfc_epsilon))+1
                  alpha_Q=Zeta_C*Zeta_D*EtaInv
                  R_Q=INT(1.0_dp/(SQRT(2.0_dp*alpha_Q)*erfc_epsilon))+1

                  if(dsqrt(rpq2).gt.R_P+R_Q) then 
      
                    Kab=dsqrt(2.0d0)*pi**1.25d0*ZetaInv*dexp(-alpha_P*rab2)
                    Kcd=dsqrt(2.0d0)*pi**1.25d0*EtaInv*dexp(-alpha_Q*rcd2)
                    theta_w=1.0d0/(1.0d0/alpha_P+1.0d0/alpha_Q+82.644628099173553719)


                    far_eri = max_contraction*(Kab*Kcd*derfc(theta_w**0.5d0*dsqrt(rpq2))/dsqrt(rpq2))
                    if(far_eri.gt.hfx_parameter%eps_far) then

                     call calc_primitive_eri(lib, ra, rb, rc_trans, rd_trans,&
                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                           la, lb, lc ,ld,&
                           ncoa,ncob,ncoc,ncod,&
                           offset_a,offset_b,offset_c,offset_d, &
                           primitive_integrals, hfx_parameter,  &
                           max_contraction, tmp_max, neris_tmp,     &
                           Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                           P,Q,W,rpq2,rab,rcd,R1,R2)

                    endif
       
                  else  !near field
                     call calc_primitive_eri(lib, ra, rb, rc_trans, rd_trans,&
                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                           la, lb, lc ,ld,&
                           ncoa,ncob,ncoc,ncod,&
                           offset_a,offset_b,offset_c,offset_d, &
                           primitive_integrals, hfx_parameter,  &
                           max_contraction,tmp_max, neris_tmp,     &
                           Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                           P,Q,W,rpq2,rab,rcd,R1,R2)            

                  endif 
      
               endif

               cart_estimate = MAX(tmp_max,cart_estimate)


             enddo
          enddo
       enddo
    enddo

    IF(cart_estimate .ge. hfx_parameter%eps_schwarz) THEN

    CALL dgemm("T","N",ncob*ncoc*ncod,nso(la),ncoa,&
               1.0_dp, primitive_integrals(1),ncoa,&
               sphia(1,1), ncoa ,&
               0.0_dp, T1(1),ncob*ncoc*ncod)
  
    CALL dgemm("T","N",nso(la)*ncoc*ncod,nso(lb),ncob,&
               1.0_dp, T1(1),ncob,&
               sphib(1,1), ncob ,&
               0.0_dp, primitive_integrals(1),nso(la)*ncoc*ncod)

    CALL dgemm("T","N",nso(la)*nso(lb)*ncod,nso(lc),ncoc,&
               1.0_dp, primitive_integrals(1),ncoc,&
               sphic(1,1), ncoc ,&
               0.0_dp, T1(1),nso(la)*nso(lb)*ncod)

    CALL dgemm("T","N",nso(la)*nso(lb)*nso(lc),nso(ld),ncod,&
               1.0_dp, T1(1),ncod,&
               sphid(1,1), ncod ,&
               0.0_dp, primitive_integrals(1),nso(la)*nso(lb)*nso(lc))
    ENDIF

    do dm=1,nso(ld)
       do cm=1,nso(lc)
          do bm=1,nso(lb)
             do am=1,nso(la)
                index_primitive_integrals=(dm-1)*nso(lc)*nso(lb)*nso(la) +     &
                                          (cm-1)*nso(lb)*nso(la)+              &
                                          (bm-1)*nso(la)+                      &
                                           am   
                eri(am,bm,cm,dm)=primitive_integrals(index_primitive_integrals)   
             enddo
          enddo
       enddo
    enddo

    deallocate(primitive_integrals)
    deallocate(T1) 
              
    end subroutine 


    subroutine calc_contract_deriv_eri( deriv, cell, rcell, ra, rb, rc, rd,  &
                                        npgfa, npgfb, npgfc, npgfd,    &
                                        la, lb, lc, ld,                & 
                                        ncoa, ncob, ncoc, ncod,        & 
                                        zeta, zetb, zetc, zetd,        & 
                                        sphia, sphib, sphic, sphid,    &
                                        hfx_parameter, eri_force,      &
                                        max_contraction, max_val2_set, &
                                        log10_eps_schwarz, log10_pmax,R1_pgf, R2_pgf, pgf1, pgf2 )
 


    !input and output variables 
    implicit none
    type(lib_deriv)  :: deriv
    integer npgfa,npgfb,npgfc,npgfd,la,lb,lc,ld,ncoa,ncob,ncoc,ncod
    real(dp) ::  cell(3,3), rcell(3,3), ra(3),rb(3),rc(3),rd(3),   &
                 zeta(npgfa),zetb(npgfb),zetc(npgfc),zetd(npgfd),  & 
                 sphia(ncoa,nso(la)),sphib(ncob,nso(lb)), &
                 sphic(ncoc,nso(lc)),sphid(ncod,nso(ld))           
    type(hfx_input_parameter)  :: hfx_parameter
    real(dp) ::  eri_force( nso(la),nso(lb),nso(lc),nso(ld),12 )

    REAL(dp), INTENT(IN)      ::  max_contraction
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: R1_pgf, R2_pgf, pgf1, pgf2

    REAL(dp)                                 :: max_val2_set, &
                                                log10_eps_schwarz, log10_pmax

    !internal variables
    real(dp),dimension(:),allocatable,TARGET :: primitive_force
    real(dp),dimension(:),allocatable        :: work_forces
    real(dp),dimension(:),POINTER            :: T1, T2
 
    integer ipgf,jpgf,kpgf,lpgf,offset_a,offset_b,offset_c,offset_d,  &
            am,bm,cm,dm,ieri,jeri,keri,leri,i,j,k,l,index_primitive_force, &
            coord
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, &
      RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
      Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv,rab(3),rcd(3),rab2,rcd2,      &
      r_temp(3),r_pbc_temp(3) ,   &
      alpha_P,alpha_Q,R_P,R_Q,Kab,Kcd,theta_w, far_eri, &
      rc_trans(3),rd_trans(3),R1, R2, pgf_max_1, pgf_max_2, cart_estimate


    real(dp),parameter :: pi=3.14159265358979323846264338_dp 
    REAL(dp), PARAMETER      :: erfc_epsilon = 0.99998871620832942100_dp      
!    real(dp) max_contraction, max_contraction_a, max_contraction_b, &
!             max_contraction_c, max_contraction_d, far_eri

    cart_estimate = 0.0_dp
 
    allocate( primitive_force(ncoa*ncob*ncoc*ncod*12) )    
    allocate( work_forces(nco(la)*nco(lb)*nco(lc)*nco(ld)*12) )
    allocate( T1(ncoa*ncob*ncoc*ncod) )
    primitive_force(1:ncoa*ncob*ncoc*ncod*12)=0.0d0
    T1(1:ncoa*ncob*ncoc*ncod)=0.0d0
 
    rab(1:3)=ra(1:3)-rb(1:3)
    rcd(1:3)=rc(1:3)-rd(1:3)
    rab2=rab(1)**2+rab(2)**2+rab(3)**2
    rcd2=rcd(1)**2+rcd(2)**2+rcd(3)**2

    DO ipgf = 1,npgfa
       offset_a = (ipgf-1)*nco(la)
       Zeta_A = zeta(ipgf)

       DO jpgf = 1,npgfb
          pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
          IF( pgf_max_1 + max_val2_set + log10_pmax < log10_eps_schwarz) CYCLE
          Zeta_B = zetb(jpgf)
          Zeta1 = Zeta_A + Zeta_B
          ZetaInv = 1.0d0/Zeta1
          S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
          P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv
          offset_b = (jpgf-1)*nco(lb)
 
          DO kpgf = 1,npgfc
             offset_c = (kpgf-1)*nco(lc)
             Zeta_C = zetc(kpgf)
           
             DO lpgf = 1,npgfd
                pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
                IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz ) CYCLE
                Zeta_D = zetd(lpgf)
                Eta  =  Zeta_C + Zeta_D
                EtaInv = 1.0d0/Eta
                ZetapEtaInv = Zeta1+Eta
                ZetapEtaInv = 1.0d0/ZetapEtaInv
                Rho = Zeta1*Eta*ZetapEtaInv
                RhoInv = 1.0d0/Rho
                S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
                Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
                r_temp(1:3)=Q(1:3)-P(1:3)
                r_pbc_temp = pbc(r_temp,cell,rcell)

                Q(1:3)=P(1:3)+r_pbc_temp(1:3)
                rc_trans(1:3)=rc(1:3)+r_pbc_temp(1:3)-r_temp(1:3)
                rd_trans(1:3)=rd(1:3)+r_pbc_temp(1:3)-r_temp(1:3)            

                rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
                W(1:3) = (Zeta1*P(1:3)+Eta*Q(1:3))*ZetapEtaInv
                offset_d = (lpgf-1)*nco(ld)
               
                tmp_max = 0.0_dp

                if(.not.hfx_parameter%far_field) then

                    CALL calc_primitive_deriv_eri(deriv, ra, rb, rc_trans, rd_trans,&
                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                           la, lb, lc ,ld,work_forces,&
                           ncoa,ncob,ncoc,ncod,&
                           primitive_force, hfx_parameter, &
                           max_contraction, tmp_max,       &
                           offset_a,offset_b,offset_c,offset_d, &
                           Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                           P,Q,W,rpq2,rab,rcd)
 
                else

                   alpha_P=Zeta_A*Zeta_B*ZetaInv
                   R_P=INT(1.0_dp/(SQRT(2.0_dp*alpha_P)*erfc_epsilon))+1
                   alpha_Q=Zeta_C*Zeta_D*EtaInv
                   R_Q=INT(1.0_dp/(SQRT(2.0_dp*alpha_Q)*erfc_epsilon))+1

                   if(dsqrt(rpq2).gt.R_P+R_Q) then ! far field

                      Kab=dsqrt(2.0d0)*pi**1.25d0*ZetaInv*dexp(-alpha_P*rab2)
                      Kcd=dsqrt(2.0d0)*pi**1.25d0*EtaInv*dexp(-alpha_Q*rcd2)
                      theta_w=1.0d0/(1.0d0/alpha_P+1.0d0/alpha_Q+82.644628099173553719)  !1.0d0/(0.11*0.11)

                      far_eri = max_contraction*(Kab*Kcd*derfc(theta_w**0.5d0*dsqrt(rpq2))/dsqrt(rpq2))
          
                      if(far_eri.gt.hfx_parameter%eps_far) then

                         CALL calc_primitive_deriv_eri(deriv, ra, rb, rc_trans, rd_trans,&
                                   zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                   la, lb, lc ,ld,work_forces,&
                                   ncoa,ncob,ncoc,ncod,&
                                   primitive_force, hfx_parameter, &
                                   max_contraction, tmp_max,       &
                                   offset_a,offset_b,offset_c,offset_d, &
                                   Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                                   P,Q,W,rpq2,rab,rcd)
   
                      endif

                   else !near field
             
                     CALL calc_primitive_deriv_eri(deriv, ra, rb, rc_trans, rd_trans,&
                               zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                               la, lb, lc ,ld,work_forces,&
                               ncoa,ncob,ncoc,ncod,&
                               primitive_force, hfx_parameter, &
                               max_contraction, tmp_max,       &
                               offset_a,offset_b,offset_c,offset_d, &
                               Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                               P,Q,W,rpq2,rab,rcd)
                   endif

                endif

                cart_estimate = MAX(tmp_max,cart_estimate)


             enddo
          enddo
       enddo
    enddo

     IF(cart_estimate .ge. hfx_parameter%eps_schwarz) THEN

     do  coord = 1,12
       T2=>primitive_force((coord-1)*ncoa*ncob*ncoc*ncod+1:coord*ncoa*ncob*ncoc*ncod)  

            CALL dgemm("T","N",ncob*ncoc*ncod,nso(la),ncoa,&
                       1.0_dp, T2(1),ncoa,&
                       sphia(1,1), ncoa ,&
                       0.0_dp, T1(1),ncob*ncoc*ncod)

            CALL dgemm("T","N",nso(la)*ncoc*ncod,nso(lb),ncob,&
                       1.0_dp, T1(1),ncob,&
                       sphib(1,1), ncob ,&
                       0.0_dp, T2(1),nso(la)*ncoc*ncod)

            CALL dgemm("T","N",nso(la)*nso(lb)*ncod,nso(lc),ncoc,&
                       1.0_dp, T2(1),ncoc,&
                       sphic(1,1), ncoc ,&
                       0.0_dp, T1(1),nso(la)*nso(lb)*ncod)

            CALL dgemm("T","N",nso(la)*nso(lb)*nso(lc),nso(ld),ncod,&
                       1.0_dp, T1(1),ncod,&
                       sphid(1,1), ncod ,&
                       0.0_dp, T2(1),nso(la)*nso(lb)*nso(lc))

     do dm=1,nso(ld)
         do cm=1,nso(lc)
            do bm=1,nso(lb)
               do am=1,nso(la)
               
               index_primitive_force = (dm-1)*nso(lc)*nso(lb)*nso(la) +     &
                                       (cm-1)*nso(lb)*nso(la)+              &
                                       (bm-1)*nso(la)+                      &
                                        am   
               eri_force(am,bm,cm,dm,coord)=T2(index_primitive_force)   
              enddo
            enddo
         enddo
      enddo
      enddo !coord
    ENDIF  
  
    deallocate(primitive_force)
    deallocate(work_forces)
    deallocate(T1) 
              
    end subroutine 

    subroutine calc_contract_eri2( lib, cell, rcell, ra, rb, rc, rd,  &
                                  npgfa, npgfb, npgfc, npgfd, &
                                  la, lb, lc, ld,             & 
                                  ncoa, ncob, ncoc, ncod,     & 
                                  zeta, zetb, zetc, zetd,     &
                                  sphia, sphib, sphic, sphid, &
                                  hfx_parameter, eri ) 

    !input and output variables 
    type(lib_int)                            :: lib 
    integer npgfa,npgfb,npgfc,npgfd,la,lb,lc,ld,ncoa,ncob,ncoc,ncod
    type(hfx_input_parameter)  :: hfx_parameter    
    real(dp) ::  cell(3,3), rcell(3,3),                            &
                 ra(3),rb(3),rc(3),rd(3), rab(3), rcd(3),          &
                 zeta(npgfa),zetb(npgfb),zetc(npgfc),zetd(npgfd),  & 
                 sphia(ncoa,nso(la)),sphib(ncob,nso(lb)), &
                 sphic(ncoc,nso(lc)),sphid(ncod,nso(ld))           
    real(dp) ::  eri( nso(la),nso(lb),nso(lc),nso(ld) )

    !internal variables
    real(dp),dimension(:),allocatable  ::  primitive_integrals
    real(dp),dimension(:),allocatable  ::  T1
 
    integer ipgf,jpgf,kpgf,lpgf,offset_a,offset_b,offset_c,offset_d,  &
            am,bm,cm,dm,ieri,i,j,k,l,index_primitive_integrals
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, &
      RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
      Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv,rab2,rcd2,      &
      r_temp(3),r_pbc_temp(3) ,   &
      alpha_P,alpha_Q,R_P,R_Q,Kab,Kcd,theta_w,&
      rc_trans(3),rd_trans(3)
    real(dp), parameter  ::  pi=3.14159265358979323846264338_dp 
    real(dp), parameter  ::  erfc_epsilon = 0.99998871620832942100_dp      

    allocate( primitive_integrals(ncoa*ncob*ncoc*ncod) )    
    allocate( T1(ncoa*ncob*ncoc*ncod) ) 
    primitive_integrals(1:ncoa*ncob*ncoc*ncod)=0.0d0
    T1(1:ncoa*ncob*ncoc*ncod)=0.0d0

    rab(1:3)=ra(1:3)-rb(1:3)
    rcd(1:3)=rc(1:3)-rd(1:3)
    rab2=rab(1)**2+rab(2)**2+rab(3)**2
    rcd2=rcd(1)**2+rcd(2)**2+rcd(3)**2

    DO ipgf = 1,npgfa
       offset_a = (ipgf-1)*nco(la)
       Zeta_A = zeta(ipgf)

      DO jpgf = 1,npgfb
         Zeta_B = zetb(jpgf)
         Zeta1 = Zeta_A + Zeta_B
         ZetaInv = 1.0d0/Zeta1
         S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
         P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv
         offset_b = (jpgf-1)*nco(lb)

         DO kpgf = 1,npgfc
            offset_c = (kpgf-1)*nco(lc)
            Zeta_C = zetc(kpgf)

            DO lpgf = 1,npgfd
               Zeta_D = zetd(lpgf)
               Eta  =  Zeta_C + Zeta_D
               EtaInv = 1.0d0/Eta
               ZetapEtaInv = Zeta1+Eta
               ZetapEtaInv = 1.0d0/ZetapEtaInv
               Rho = Zeta1*Eta*ZetapEtaInv
               RhoInv = 1.0d0/Rho
               S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
               Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
               r_temp(1:3)=Q(1:3)-P(1:3)
               r_pbc_temp = pbc(r_temp,cell,rcell)

               Q(1:3)=P(1:3)+r_pbc_temp(1:3)
               rc_trans(1:3)=rc(1:3)+r_pbc_temp(1:3)-r_temp(1:3)
               rd_trans(1:3)=rd(1:3)+r_pbc_temp(1:3)-r_temp(1:3)            

               rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
               W(1:3) = (Zeta1*P(1:3)+Eta*Q(1:3))*ZetapEtaInv
               offset_d = (lpgf-1)*nco(ld)

               call calc_primitive_eri2(lib, ra, rb, rc_trans, rd_trans,&
                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                           la, lb, lc ,ld,&
                           ncoa,ncob,ncoc,ncod,&
                           offset_a,offset_b,offset_c,offset_d, &
                           primitive_integrals, hfx_parameter, &
                           Zeta1,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                           P,Q,W,rpq2,rab,rcd)

             enddo
          enddo
       enddo
    enddo


    CALL dgemm("T","N",ncob*ncoc*ncod,nso(la),ncoa,&
               1.0_dp, primitive_integrals(1),ncoa,&
               sphia(1,1), ncoa ,&
               0.0_dp, T1(1),ncob*ncoc*ncod)
  
    CALL dgemm("T","N",nso(la)*ncoc*ncod,nso(lb),ncob,&
               1.0_dp, T1(1),ncob,&
               sphib(1,1), ncob ,&
               0.0_dp, primitive_integrals(1),nso(la)*ncoc*ncod)

    CALL dgemm("T","N",nso(la)*nso(lb)*ncod,nso(lc),ncoc,&
               1.0_dp, primitive_integrals(1),ncoc,&
               sphic(1,1), ncoc ,&
               0.0_dp, T1(1),nso(la)*nso(lb)*ncod)

    CALL dgemm("T","N",nso(la)*nso(lb)*nso(lc),nso(ld),ncod,&
               1.0_dp, T1(1),ncod,&
               sphid(1,1), ncod ,&
               0.0_dp, primitive_integrals(1),nso(la)*nso(lb)*nso(lc))

    do dm=1,nso(ld)
       do cm=1,nso(lc)
          do bm=1,nso(lb)
             do am=1,nso(la)
                index_primitive_integrals=(dm-1)*nso(lc)*nso(lb)*nso(la) +     &
                                          (cm-1)*nso(lb)*nso(la)+              &
                                          (bm-1)*nso(la)+                      &
                                           am   
                eri(am,bm,cm,dm)=primitive_integrals(index_primitive_integrals)   
             enddo
          enddo
       enddo
    enddo

    deallocate(primitive_integrals)
    deallocate(T1) 
              
    end subroutine 
end module contract_eri
