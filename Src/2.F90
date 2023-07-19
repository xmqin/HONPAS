!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the Libint-Library 
!> \note
!>      IMPORTANT NOTE : this file currently is for a libint configured for
!>                       LIBINT_MAX_AM 5 and LIBINT_MAX_AM1 4
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE primitive_eri

  !USE f77_blas
  USE gamma_F,                         ONLY: fgamma
  USE libint_wrapper_types,            ONLY: lib_int, lib_deriv, &
                                             prim_data
  USE libint_wrapper,                  ONLY: get_eris,get_derivs
  USE kinds,                           ONLY: dp 
  USE mathconstants
  USE atm_types,                       ONLY: nco

  IMPLICIT NONE
  PRIVATE
  PUBLIC calc_primitive_eri,& 
         calc_primitive_deriv_eri
         !set_eps_cutoff
  INTEGER, DIMENSION(12), PARAMETER :: full_perm1 = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm2 = (/4,5,6,1,2,3,7,8,9,10,11,12/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm3 = (/1,2,3,4,5,6,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm4 = (/4,5,6,1,2,3,10,11,12,7,8,9/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm5 = (/7,8,9,10,11,12,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm6 = (/7,8,9,10,11,12,4,5,6,1,2,3/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm7 = (/10,11,12,7,8,9,1,2,3,4,5,6/)
  INTEGER, DIMENSION(12), PARAMETER :: full_perm8 = (/10,11,12,7,8,9,4,5,6,1,2,3/)

  !CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_interface'
  !REAL(KIND=dp), SAVE                :: eps_cutoff

!***

  CONTAINS

! *****************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_primitive_eri(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,&
                          ncoa,ncob,ncoc,ncod,&
                          offset_a,offset_b,offset_c,offset_d,&
                          primitives,         &
                          Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                          P,Q,W,rpq2,AB,CD)
    
    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, offset_a, &
                                                offset_b, offset_c, offset_d
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitives
    REAL(dp), INTENT(IN)                     :: Zeta, ZetaInv, Eta, EtaInv, &
                                                ZetapEtaInv, Rho, RhoInv, &
                                                S1234, P(3), Q(3), W(3), &
                                                rpq2, AB(3), CD(3)

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, p1, p2, p3, &
                                                p4, perm_case, q1, q2, q3, q4
    LOGICAL                                  :: do_it
    REAL(dp), DIMENSION(:), POINTER          :: p_work
    TYPE(prim_data), TARGET                  :: prim
  
    m_max = n_a+n_b+n_c+n_d

    !write(6,*) '-------------------------------------------------'
    !write(6,*) 'n_a~n_d:',n_a,n_b,n_c,n_d
    !write(6,*) 'nco(n_a)~nco(n_d):',nco(n_a),nco(n_b),nco(n_c),nco(n_d)
    !write(6,*) 'ncoa~ncod:',ncoa, ncob, ncoc, ncod    

    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize

    do_it = .TRUE.
    IF(m_max/=0) THEN
      perm_case = 1
      IF(n_a<n_b) THEN
        perm_case = perm_case + 1
      END IF
      IF(n_c<n_d) THEN
        perm_case = perm_case + 2
      END IF
      IF(n_a+n_b > n_c+n_d) THEN
        perm_case = perm_case + 4
      END IF

     

      SELECT CASE(perm_case)
        CASE(1)
          CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                             prim, do_it, n_a+n_b, n_c+n_d,  &
                            Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=AB!A-B
          lib%CD=CD!C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)

          DO i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            q1 = offset_a + i
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              q2 = offset_b + j
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c + k
                DO l = 1,nco(n_d)
                  q4 = offset_d + l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(2)
          CALL build_quartet_data(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                             prim, do_it, n_b+n_a, n_c+n_d, &
                            Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-AB!B-A
          lib%CD=CD!C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b+j
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              q1 = offset_a+i
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c+k
                DO l = 1,nco(n_d)
                  q4 = offset_d+l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(3)
          CALL build_quartet_data(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            prim, do_it, n_a+n_b, n_d+n_c, &
                            Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=AB!A-B
          lib%CD=-CD!D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)

          DO i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            q1 = offset_a + i
            DO j = 1,nco(n_b)
              q2 = offset_b + j
              p2 = (p1 + j-1)*nco(n_d)
              DO l = 1,nco(n_d)
                q4 = offset_d + l
                p3 = (p2 + l-1) * nco(n_c)
                DO k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3+k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(4)
          CALL build_quartet_data(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            prim, do_it, n_b+n_a, n_d+n_c, &
                            Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-AB!B-A
          lib%CD=-CD!D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b + j
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              q1 = offset_a + i
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                q4 = offset_d + l
                DO k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3 + k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(5)
          CALL build_quartet_data(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            prim, do_it, n_c+n_d, n_a+n_b ,&
                            Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=CD!C-D
          lib%CD=AB!A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)

          DO k = 1,nco(n_c)
            q3 = offset_c + k
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_a)
              DO i = 1,nco(n_a)
                q1 = offset_a + i
                p3 = (p2 + i-1)*nco(n_b)
                DO j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(6)
          CALL build_quartet_data(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            prim, do_it, n_c+n_d, n_b+n_a ,&
                            Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=CD!C-D
          lib%CD=-AB!B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            q3 = offset_c + k
            DO l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_b)
              DO j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1)*nco(n_a)
                DO i = 1,nco(n_a)
                  p4 = p3+i
                  q1 = offset_a + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(7)
          CALL build_quartet_data(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            prim, do_it, n_d+n_c, n_a+n_b ,&
                            Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-CD!D-C
          lib%CD=AB!A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            q4 = offset_d + l
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              q3 = offset_c + k
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                q1 = offset_a + i
                DO j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
        CASE(8)
          CALL build_quartet_data(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            prim, do_it, n_d+n_c, n_b+n_a ,&
                            Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            Q,P,W,rpq2)
          IF( .NOT. do_it ) RETURN
          lib%AB=-CD!D-C
          lib%CD=-AB!B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)

          DO l = 1,nco(n_d)
            q4 = offset_d + l
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              q3 = offset_c + k
              p2 = (p1 + k-1) * nco(n_b)
              DO j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1) * nco(n_a)
                DO i = 1,nco(n_a)
                  q1 = offset_a + i
                  p4 = p3 + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                END DO
              END DO
            END DO
          END DO
      END SELECT
    ELSE
      !write(6,*) 'before build_quartet_data'
      CALL build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            prim, do_it, 0, 0,         &
                            Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
      !write(6,*) 'after build_quartet_data_1',prim%F(1)
      !write(6,*) 'after build_quartet_data_2',offset_a,offset_b,offset_c,offset_d
      IF( .NOT. do_it ) RETURN
      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = prim%F(1)
    END IF
    
   !write(6,*) primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) 


  END SUBROUTINE calc_primitive_eri



! *****************************************************************************
!> \brief Fill data structure used in libint
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                prim, do_it, np, nq ,&
                                Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                                P,Q,W,PQ2)
     
    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(prim_data)                          :: prim
    LOGICAL, INTENT(INOUT)                   :: do_it
    INTEGER                                  :: np, nq
    REAL(dp)                                 ::  Zeta, ZetaInv, Eta, &
                                                EtaInv, ZetapEtaInv, Rho, &
                                                RhoInv, S1234, P(3), Q(3), &
                                                W(3), PQ2, omega, omega2, &
                                                omega_corr, omega_corr2

    INTEGER                                  :: i
    LOGICAL                                  :: use_gamma
    REAL(KIND=dp)                            :: factor,T, tmp
    REAL(KIND=dp), DIMENSION(17)    :: Fm

    omega=0.11d0  !0.11

    T = Rho*PQ2

    do_it = .TRUE.

      !CASE(do_hfx_potential_coulomb)
       ! CALL fgamma(m_max,T,prim%F)
       ! factor = 2.0_dp*Pi*RhoInv

      !CASE(do_hfx_potential_short)
       CALL fgamma(m_max,T,prim%F)
        omega2 = omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = - omega_corr
        DO i=1,m_max+1
          prim%F(i)=prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        END DO
        factor = 2.0_dp*Pi*RhoInv



    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(1:3,1) = P-A
!    prim%U(1:3,2) = P-B                 !Not used in libint
    prim%U(1:3,3) = Q-C
!    prim%U(1:3,4) = Q-D                 !Not used in libint
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
!    prim%twozeta_a = 0.0_dp                      !Not used in libint
!    prim%twozeta_b = 0.0_dp                      !Not used in libint
!    prim%twozeta_c = 0.0_dp                      !Not used in libint 
!    prim%twozeta_d = 0.0_dp                      !Not used in libint
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
!    prim%ss_r12_ss = 0.0_dp                      !Not used in libint,libderiv
  END SUBROUTINE build_quartet_data

! *****************************************************************************
!> \brief Fill data structure used in libderiv
!> \par History
!>      03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE build_deriv_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                               prim, do_it, np, nq,&
                              Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              P,Q,W,PQ2 )
     
    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    !TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(prim_data)                          :: prim
    LOGICAL                                  :: do_it
    INTEGER                                  :: np, nq
    REAL(dp)                                 ::  Zeta, ZetaInv, Eta, &
                                                EtaInv, ZetapEtaInv, Rho, &
                                                RhoInv, S1234, P(3), Q(3), &
                                                W(3), PQ2

    INTEGER                                  :: i
    LOGICAL                                  :: use_gamma
    REAL(KIND=dp)                            :: factor,omega, omega2, omega_corr, &
                                                omega_corr2, T, tmp
    REAL(KIND=dp), DIMENSION(17)             :: Fm


    omega=0.11d0  !0.11 
    T = Rho*PQ2
    do_it = .TRUE.

      !CASE(do_hfx_potential_short)
        CALL fgamma(m_max,T,prim%F)
        omega2 = omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = - omega_corr
        DO i=1,m_max+1
          prim%F(i)=prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2 
        END DO
        factor = 2.0_dp*Pi*RhoInv
    tmp = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    DO i=1,m_max+1
       prim%F(i)=prim%F(i)*factor
    ENDDO
    prim%U(:,1) = P-A
    prim%U(:,2) = P-B
    prim%U(:,3) = Q-C
    prim%U(:,4) = Q-D
    prim%U(:,5) = W-P
    prim%U(:,6) = W-Q
    prim%twozeta_a = 2.0_dp*Zeta_A
    prim%twozeta_b = 2.0_dp*Zeta_B
    prim%twozeta_c = 2.0_dp*Zeta_C
    prim%twozeta_d = 2.0_dp*Zeta_D
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv
!    prim%ss_r12_ss = 0.0_dp                      !Not used in libint,libderiv

!    C = C1
!    D = D1
  END SUBROUTINE build_deriv_data

! *****************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive quartet
!> \par History
!>      03.2007 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_primitive_deriv_eri(deriv,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                n_a,n_b,n_c,n_d,work_forces,&
                                ncoa, ncob, ncoc, ncod,&
                                primitive_forces,&
                                offset_a, offset_b, offset_c, offset_d,&
                                Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                                P,Q,W,rpq2,AB,CD)
    
    TYPE(lib_deriv)                          :: deriv
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*&
      nco(n_c)*nco(n_d), 12)                 :: work_forces
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    INTEGER, INTENT(IN)                      :: offset_a, offset_b, offset_c, &
                                                offset_d
    REAL(dp), INTENT(IN)                     :: Zeta, ZetaInv, Eta, EtaInv, &
                                                ZetapEtaInv, Rho, RhoInv, &
                                                S1234, P(3), Q(3), W(3), &
                                                rpq2, AB(3), CD(3)

    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, n, p1, p2, p3, &
                                                perm_case
    LOGICAL                                  :: do_it
    REAL(dp)                                 :: tmp_max
    TYPE(prim_data), TARGET                  :: prim

!permutation of configuration

    perm_case = 1
    IF(n_a<n_b) THEN
      perm_case = perm_case + 1
    END IF

    IF(n_c<n_d) THEN
      perm_case = perm_case + 2
    END IF
 
    IF( n_a+n_b > n_c+n_d) THEN
      perm_case = perm_case + 4
    END IF

    do_it = .TRUE.
    m_max = n_a+n_b+n_c+n_d
    m_max = m_max + 1
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize=mysize
    SELECT CASE(perm_case)
      CASE(1)
        CALL build_deriv_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                              prim, do_it, n_a+n_b, n_c+n_d,  &
                              Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              P,Q,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=AB!A-B
        deriv%CD=CD!C-D
        CALL get_derivs(n_d, n_c, n_b, n_a, deriv, prim, work_forces, a_mysize)
        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                                  work_forces(i,k+3) + &
                                                  work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm1(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm1(n))+&
                                   work_forces(p3+l,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(2)
        CALL build_deriv_data(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                              prim, do_it, n_b+n_a, n_c+n_d,&
                              Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              P,Q,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=-AB!B-A
        deriv%CD=CD!C-D
        CALL get_derivs(n_d, n_c, n_a, n_b, deriv, prim, work_forces, a_mysize)

        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                                  work_forces(i,k+3) + &
                                                  work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              DO k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                DO l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm2(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm2(n))+&
                                   work_forces(p3+l,n)
                END DO
              END DO
            END DO
          END DO
        END DO
     CASE(3)
        CALL build_deriv_data(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                              prim, do_it, n_a+n_b, n_d+n_c,&
                              Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              P,Q,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=AB!A-B
        deriv%CD=-CD!D-C
        CALL get_derivs(n_c, n_d, n_b, n_a, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
 !         tmp_max = 0.0_dp
 !         DO i=1,mysize
 !           tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
 !         END DO
 !         tmp_max = tmp_max*max_contraction
 !         tmp_max_all = MAX(tmp_max_all, tmp_max)
 !         IF(tmp_max<eps_schwarz) CYCLE

          DO i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            DO j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1) * nco(n_c)
                DO k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm3(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm3(n))+&
                                   work_forces(p3+k,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(4)
        CALL build_deriv_data(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                              prim, do_it, n_b+n_a, n_d+n_c,&
                              Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              P,Q,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=-AB!B-A
        deriv%CD=-CD!D-C
        CALL get_derivs(n_c, n_d, n_a, n_b, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12          
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            DO i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              DO l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                DO k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm4(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm4(n))+&
                                   work_forces(p3+k,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(5)
        CALL build_deriv_data(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                              prim, do_it, n_c+n_d, n_a+n_b,&
                               Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                               Q,P,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=CD!C-D
        deriv%CD=AB!A-B
        CALL get_derivs(n_b, n_a, n_d, n_c, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1)*nco(n_b)
                DO j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm5(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm5(n))+&
                                   work_forces(p3+j,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(6)
        CALL build_deriv_data(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                              prim, do_it, n_c+n_d, n_b+n_a,&
                              Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              Q,P,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=CD!C-D
        deriv%CD=-AB!B-A
        CALL get_derivs(n_a, n_b, n_d, n_c, deriv, prim, work_forces, a_mysize)

        DO k=4,6
          DO i=1,mysize
             work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
          ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            DO l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1)*nco(n_a)
                DO i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm6(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm6(n))+&
                                   work_forces(p3+i,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(7)
        CALL build_deriv_data(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                              prim, do_it, n_d+n_c, n_a+n_b,&
                              Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              Q,P,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=-CD!D-C
        deriv%CD=AB!A-B
        CALL get_derivs(n_b, n_a, n_c, n_d, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              DO i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                DO j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm7(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm7(n))+&
                                   work_forces(p3+j,n)
                END DO
              END DO
            END DO
          END DO
        END DO
      CASE(8)
        CALL build_deriv_data(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                              prim, do_it, n_d+n_c,n_b+n_a,&
                              Eta,EtaInv,Zeta,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                              Q,P,W,rpq2)
        IF( .NOT. do_it ) RETURN
        deriv%AB=-CD!D-C
        deriv%CD=-AB!B-A
        CALL get_derivs(n_a, n_b, n_c, n_d, deriv, prim, work_forces, a_mysize)

        DO k=4,6
           DO i=1,mysize
              work_forces(i,k) = - 1.0_dp* (work_forces(i,k-3) + &
                                               work_forces(i,k+3) + &
                                               work_forces(i,k+6) )
           ENDDO
        END DO
        DO n=1,12
!          tmp_max = 0.0_dp
!          DO i=1,mysize
!            tmp_max = MAX(tmp_max, ABS(work_forces(i,n)))
!          END DO
!          tmp_max = tmp_max*max_contraction
!          tmp_max_all = MAX(tmp_max_all, tmp_max)
!          IF(tmp_max<eps_schwarz) CYCLE

          DO l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            DO k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_b)
              DO j = 1,nco(n_b)
                p3 = (p2 + j-1) * nco(n_a)
                DO i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm8(n))=&
!                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l,full_perm8(n))+&
                                   work_forces(p3+i,n)
                END DO
              END DO
            END DO
          END DO
        END DO
    END SELECT

!shanghui_force
     ! write(6,*) primitive_forces(offset_a+1,offset_b+1,offset_c+1,offset_d+1,1:12)

    ! stop
  END SUBROUTINE calc_primitive_deriv_eri




! *****************************************************************************
!> \brief sets the threshold for truncated calculations
!> \par History
!>      12.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
!  SUBROUTINE set_eps_cutoff(eps_schwarz)
!    REAL(dp), INTENT(IN)                     :: eps_schwarz
!
!    eps_cutoff = SQRT(-LOG(eps_schwarz))
!  END SUBROUTINE set_eps_cutoff



END MODULE primitive_eri
