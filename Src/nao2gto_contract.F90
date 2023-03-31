! *** Module: nao2gto_contract ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief NAO2GTO contraction routines
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.2010 Edited [Xinming Qin]
!!      - 01.2016 Edited [Honghui Shang]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
! *****************************************************************************
module nao2gto_contract

  use precision

  implicit none

  private

  real(dp), parameter :: erfc_epsilon = 0.99998871620832942100_dp

  public :: calc_contract_eri, calc_contract_deriv_eri, calc_contract_eri2

contains

! *****************************************************************************
!> \brief Computes contractions of spherical Gaussians
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] cell: supercell in real space 
!! \param[in] rcell: supercell in reciprocal space
!! \param[in] ra: ...
!! \param[in] rb: ...
!! \param[in] rc: ...
!! \param[in] rd: ...
!! \param[in] npgfa: ...
!! \param[in] npgfb: ...
!! \param[in] npgfc: ...
!! \param[in] npgfd: ...
!! \param[in] la: ...
!! \param[in] lb: ...
!! \param[in] lc: ...
!! \param[in] ld: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] zeta: ...
!! \param[in] zetb: ...
!! \param[in] zetc: ...
!! \param[in] zetd: ...
!! \param[in] sphia: ...
!! \param[in] sphib: ...
!! \param[in] sphic: ...
!! \param[in] sphid: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri: ...
! *****************************************************************************
  subroutine calc_contract_eri(libint_data, cell, rcell, ra, rb, rc, rd,  &
&                npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                ncoa, ncob, ncoc, ncod, zeta, zetb, zetc, zetd, &
&                sphia, sphib, sphic, sphid, neris_tmp, max_contraction, &
&                max_val2_set, log10_pmax, R1_pgf, R2_pgf, pgf1, &
&                pgf2, hfx_opts, eri)

    use units,              only: pi
    use alloc,              only: de_alloc, re_alloc
    use nao2gto_data
    use nao2gto_libint,     only: Libint_t
    use nao2gto_pbc,        only: trans_pbc
    use nao2gto_primitive,  only: calc_primitive_eri
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libint_t), intent(inout) :: libint_data
    integer, intent(in)   :: npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                            ncoa, ncob, ncoc, ncod
    integer(int_8), intent(out) :: neris_tmp
    real(dp), intent(in)  :: cell(3,3), rcell(3,3), &
&     ra(3), rb(3), rc(3), rd(3), &
&     zeta(npgfa), zetb(npgfb), zetc(npgfc), zetd(npgfd), &
&     sphia(ncoa,nso(la)), sphib(ncob,nso(lb)), &
&     sphic(ncoc,nso(lc)), sphid(ncod,nso(ld))
    real(dp), intent(in) :: max_contraction, max_val2_set, log10_pmax
    real(dp), intent(out) :: eri(nso(la),nso(lb),nso(lc),nso(ld))
    type(hfx_screen_coeff_type), dimension(:,:), intent(in) :: &
&     R1_pgf, R2_pgf, pgf1, pgf2
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: ipgf     ! Counter for a loop over the first  primitive Gaussians
    integer :: jpgf     ! Counter for a loop over the second primitive Gaussians
    integer :: kpgf     ! Counter for a loop over the third  primitive Gaussians
    integer :: lpgf     ! Counter for a loop over the fourth primitive Gaussians
    integer :: offset_a, offset_b, &
&              offset_c, offset_d, am, bm, cm, dm, i, &
&              index_primitive_integrals
    integer :: ieri, j, k, l
    real(dp) :: Zeta_A  ! Exponent of the first  primitive Gaussian
    real(dp) :: Zeta_B  ! Exponent of the second primitive Gaussian
    real(dp) :: Zeta_C  ! Exponent of the third  primitive Gaussian
    real(dp) :: Zeta_D  ! Exponent of the fourth primite Gaussian
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, RhoInv, rpq2, S1234, S1234a, &
&               tmp_max, W(3), Zeta1,  &
&               ZetaInv, ZetapEtaInv, rab2, rcd2, r_temp(3), r_pbc_temp(3), &
&               alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, theta_w, &
&               rc_trans(3), rd_trans(3), rab(3), rcd(3), R1, R2
    real(dp) :: far_eri, pgf_max_1, pgf_max_2, cart_estimate

    real(dp), dimension(:), pointer :: primitive_integrals => null()
    real(dp), dimension(:), pointer :: T1 => null()

    external :: dgemm

    ! -------------------------------------------------------------------------

    call re_alloc(primitive_integrals, 1, ncoa*ncob*ncoc*ncod, &
&     name='primitive_integrals', routine='calc_contract_eri')
    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='T1', routine='calc_contract_eri')
    primitive_integrals(:) = 0.0_dp
    T1(:) = 0.0_dp

    neris_tmp = 0
    cart_estimate = 0.0_dp
    eri(:,:,:,:) = 0.0_dp

    rab(1:3) = ra(1:3) - rb(1:3)
    rcd(1:3) = rc(1:3) - rd(1:3)
    rab2 = (rab(1)**2) + (rab(2)**2) + (rab(3)**2)
    rcd2 = (rcd(1)**2) + (rcd(2)**2) + (rcd(3)**2)

    do ipgf = 1, npgfa
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

      do jpgf = 1,npgfb
        offset_b = (jpgf-1)*nco(lb)
        pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
        if ( pgf_max_1+max_val2_set+log10_pmax < log10_eps_schwarz ) cycle

        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0_dp/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv

        do kpgf = 1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)

          do lpgf = 1,npgfd
            offset_d = (lpgf - 1)*nco(ld)
            pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
            if ( pgf_max_1+pgf_max_2+log10_pmax < log10_eps_schwarz ) cycle

            Zeta_D = zetd(lpgf)
            Eta = Zeta_C + Zeta_D
            EtaInv = 1.0_dp/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
            r_temp(1:3) = Q(1:3) - P(1:3)
            call trans_pbc(r_temp, cell,rcell, r_pbc_temp)

            Q(1:3) = P(1:3) + r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)

            rpq2 = ((P(1)-Q(1))**2) + ((P(2)-Q(2))**2) + ((P(3)-Q(3))**2)
            W(1:3) = (Zeta1*P(1:3) + Eta*Q(1:3))*ZetapEtaInv

            tmp_max = 0.0_dp

            if ( (hfx_opts%potential_type .ne. 2).or.(.not. hfx_opts%farfield )) then

              call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&               la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&               offset_a, offset_b, offset_c, offset_d, &
&               primitive_integrals, hfx_opts, &
&               max_contraction, tmp_max, neris_tmp, &
&               ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&               P, Q, W, rpq2, rab, rcd, R1, R2)

            else ! farfield screening for GTO ERIs

              alpha_P = Zeta_A*Zeta_B*ZetaInv
              R_P = aint(1.0_dp/(sqrt(2.0_dp*alpha_P)*erfc_epsilon)) + 1.0_dp
              alpha_Q = Zeta_C*Zeta_D*EtaInv
              R_Q = aint(1.0_dp/(sqrt(2.0_dp*alpha_Q)*erfc_epsilon)) + 1.0_dp

              if ( sqrt(rpq2) .gt. (R_P + R_Q) ) then

                Kab = sqrt(2.0_dp)*pi**1.25_dp*ZetaInv*exp(-alpha_P*rab2)
                Kcd = sqrt(2.0_dp)*pi**1.25_dp*EtaInv*exp(-alpha_Q*rcd2)
                theta_w = 1.0_dp/ &
&                 (1.0_dp/alpha_P+1.0_dp/alpha_Q+82.644628099173553719)

                far_eri = max_contraction * &
&                 (Kab*Kcd*erfc((theta_w**0.5_dp)*sqrt(rpq2))/sqrt(rpq2))

                if ( far_eri .gt. hfx_opts%eps_farfield ) then

                  call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&                   la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&                   offset_a, offset_b, offset_c, offset_d, &
&                   primitive_integrals, hfx_opts, &
&                   max_contraction, tmp_max, neris_tmp, &
&                   ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                   P, Q, W, rpq2, rab, rcd, R1, R2)

                endif

              else ! near field

                call calc_primitive_eri(libint_data, ra, rb, rc_trans, rd_trans, &
&                 la, lb, lc ,ld, ncoa, ncob, ncoc, ncod, &
&                 offset_a, offset_b, offset_c, offset_d, &
&                 primitive_integrals, hfx_opts, &
&                 max_contraction, tmp_max, neris_tmp, &
&                 ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                 P, Q, W, rpq2, rab, rcd, R1, R2)


              endif ! sqrt(rpq2) .gt. (R_P + R_Q)

            endif ! .not. hfx_opts%farfield

            cart_estimate = max(tmp_max, cart_estimate)

          enddo ! lpgf
        enddo ! kpgf
      enddo ! jpgf
    enddo ! ipgf

    if ( cart_estimate >= hfx_opts%eps_schwarz ) then

      call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
&       1.0_dp, primitive_integrals(1), ncoa, &
&       sphia(1,1), ncoa, 0.0_dp, T1(1), ncob*ncoc*ncod)

      call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
&       1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
&       0.0_dp, primitive_integrals(1), nso(la)*ncoc*ncod)

      call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
&       1.0_dp, primitive_integrals(1), ncoc, &
&       sphic(1,1), ncoc, 0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

      call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
&       1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
&       0.0_dp, primitive_integrals(1), nso(la)*nso(lb)*nso(lc))

    end if

    do dm=1,nso(ld)
      do cm=1,nso(lc)
        do bm=1,nso(lb)
          do am=1,nso(la)
            index_primitive_integrals = &
&             (dm-1)*nso(lc)*nso(lb)*nso(la) + (cm-1)*nso(lb)*nso(la) + &
&             (bm-1)*nso(la) + am
            eri(am,bm,cm,dm) = primitive_integrals(index_primitive_integrals)
          enddo
        enddo
      enddo
    enddo

    call de_alloc(primitive_integrals, &
&     name='primitive_integrals', routine='calc_contract_eri')
    call de_alloc(T1, &
&     name='primitive_integrals', routine='calc_contract_eri')

  end subroutine calc_contract_eri

! *****************************************************************************
!> \brief Computes contracted derivatives of spherical Gaussians
!!
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] deriv_data: Libderiv data structure
!! \param[in] cell: supercell in real space
!! \param[in] rcell: supercell in reciprocal space
!! \param[in] ra: ...
!! \param[in] rb: ...
!! \param[in] rc: ...
!! \param[in] rd: ...
!! \param[in] npgfa: ...
!! \param[in] npgfb: ...
!! \param[in] npgfc: ...
!! \param[in] npgfd: ...
!! \param[in] la: ...
!! \param[in] lb: ...
!! \param[in] lc: ...
!! \param[in] ld: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[in] zeta: ...
!! \param[in] zetb: ...
!! \param[in] zetc: ...
!! \param[in] zetd: ...
!! \param[in] sphia: ...
!! \param[in] sphib: ...
!! \param[in] sphic: ...
!! \param[in] sphid: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri_force: ...
! *****************************************************************************
  subroutine calc_contract_deriv_eri(deriv_data, cell, rcell, ra, rb, rc, rd, &
&                npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                ncoa, ncob, ncoc, ncod, zeta, zetb, zetc, zetd, &
&                sphia, sphib, sphic, sphid, hfx_opts, eri_force, &
&                max_contraction, max_val2_set, log10_pmax, &
&                R1_pgf, R2_pgf, pgf1, pgf2)

    use units,             only: pi
    use alloc,             only: de_alloc, re_alloc
    use nao2gto_data
    use nao2gto_libint,    only: Libderiv_t
    use nao2gto_pbc,       only: cell_pbc
    use nao2gto_primitive, only: calc_primitive_deriv_eri
    use nao2gto_types

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout) :: deriv_data
    integer, intent(in)  :: npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, &
&                           ncoa, ncob, ncoc, ncod
    real(dp), intent(in)  :: cell(3,3), rcell(3,3), ra(3), rb(3), &
&     rc(3), rd(3), zeta(npgfa), zetb(npgfb), zetc(npgfc), zetd(npgfd), &
&     sphia(ncoa,nso(la)), sphib(ncob,nso(lb)), &
&     sphic(ncoc,nso(lc)), sphid(ncod,nso(ld))
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(out) :: eri_force(nso(la),nso(lb),nso(lc),nso(ld),12)
    real(dp), intent(in) :: max_contraction, max_val2_set, log10_pmax
    type(hfx_screen_coeff_type), dimension(:, :), pointer :: &
&     R1_pgf, R2_pgf, pgf1, pgf2

    ! Local variables
    integer  :: ipgf, jpgf, kpgf, lpgf, offset_a, offset_b, &
&     offset_c,offset_d, am, bm, cm, dm, i, index_primitive_force, coord
    real(dp) :: Eta, EtaInv, P(3), Q(3), Rho, &
&     RhoInv, rpq2, S1234, S1234a, tmp_max, W(3), Zeta1, Zeta_A, Zeta_B, &
&     Zeta_C, Zeta_D, ZetaInv, ZetapEtaInv, rab(3), rcd(3), rab2, rcd2, &
&     r_temp(3), r_pbc_temp(3), alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, &
&     theta_w, rc_trans(3), rd_trans(3)
    real(dp) :: far_eri, R1, R2, pgf_max_1, pgf_max_2, cart_estimate

    real(dp), dimension(:), pointer :: primitive_force => null()
    real(dp), dimension(:,:), pointer :: work_forces => null()
    real(dp), dimension(:), pointer :: T1 => null(), T2 => null()

    external :: dgemm

    ! -------------------------------------------------------------------------

    call re_alloc(primitive_force, 1, ncoa*ncob*ncoc*ncod*12, &
&     name='primitive_force', routine='calc_contract_deriv_eri')
    primitive_force(:) = 0.0_dp
    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='T1', routine='calc_contract_deriv_eri')
    T1(:) = 0.0_dp
    call re_alloc(work_forces, 1, nco(la)*nco(lb)*nco(lc)*nco(ld), 1, 12, &
&     name='work_forces', routine='calc_contract_deriv_eri')

    cart_estimate = 0.0_dp
    rab(1:3) = ra(1:3) - rb(1:3)
    rcd(1:3) = rc(1:3) - rd(1:3)
    rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
    rcd2 = rcd(1)**2 + rcd(2)**2 + rcd(3)**2

    do ipgf = 1,npgfa
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

      do jpgf = 1,npgfb
        offset_b = (jpgf-1)*nco(lb)
        pgf_max_1 = pgf1(jpgf,ipgf)%x(1)*rab2+pgf1(jpgf,ipgf)%x(2)
        if ( pgf_max_1+max_val2_set+log10_pmax < log10_eps_schwarz ) cycle

        Zeta_B = zetb(jpgf)
        Zeta1 = Zeta_A + Zeta_B
        ZetaInv = 1.0_dp/Zeta1
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv

        do kpgf = 1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)

          do lpgf = 1,npgfd
            offset_d = (lpgf-1)*nco(ld)
            pgf_max_2 = pgf2(lpgf,kpgf)%x(1)*rcd2+pgf2(lpgf,kpgf)%x(2)
            if (pgf_max_1+pgf_max_2+log10_pmax < log10_eps_schwarz ) cycle

            Zeta_D = zetd(lpgf)
            Eta = Zeta_C + Zeta_D
            EtaInv = 1.0_dp/Eta
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv
            Rho = Zeta1*Eta*ZetapEtaInv
            RhoInv = 1.0_dp/Rho
            S1234 = EXP(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
            r_temp(1:3) = Q(1:3)-P(1:3)
            r_pbc_temp = cell_pbc(r_temp,cell,rcell)

            Q(1:3) = P(1:3) + r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)

            rpq2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2
            W(1:3) = (Zeta1*P(1:3)+Eta*Q(1:3))*ZetapEtaInv

            tmp_max = 0.0_dp

            if ( .not. hfx_opts%farfield ) then

              call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&               rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&               zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&               ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&               max_contraction, tmp_max, &
&               offset_a, offset_b, offset_c, offset_d, ZetaInv, EtaInv, &
&               ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, rpq2, rab, rcd)

            else

              alpha_P = Zeta_A*Zeta_B*ZetaInv
              R_P = aint(1.0_dp/(sqrt(2.0_dp*alpha_P)*erfc_epsilon)) + 1.0_dp
              alpha_Q = Zeta_C*Zeta_D*EtaInv
              R_Q = aint(1.0_dp/(sqrt(2.0_dp*alpha_Q)*erfc_epsilon)) + 1.0_dp

              if ( sqrt(rpq2) > (R_P + R_Q) ) then ! far field

                Kab = sqrt(2.0_dp)*pi**1.25_dp*ZetaInv*exp(-alpha_P*rab2)
                Kcd = sqrt(2.0_dp)*pi**1.25_dp*EtaInv*exp(-alpha_Q*rcd2)

                ! Note: 1.0_dp/(0.11*0.11)
                theta_w = 1.0_dp/ &
&                 (1.0_dp/alpha_P+1.0_dp/alpha_Q+82.644628099173553719)

                far_eri = max_contraction * &
&                 (Kab*Kcd*erfc((theta_w**0.5_dp)*sqrt(rpq2))/sqrt(rpq2))

                if ( far_eri > hfx_opts%eps_farfield ) then

                  call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&                   rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&                   zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&                   ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&                   max_contraction, tmp_max, &
&                   offset_a, offset_b, offset_c, offset_d, ZetaInv, &
&                   EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                   rpq2, rab, rcd)

                endif

              else !near field

                call calc_primitive_deriv_eri(deriv_data, ra, rb, &
&                 rc_trans, rd_trans, zeta(ipgf), zetb(jpgf), &
&                 zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, work_forces, &
&                 ncoa, ncob, ncoc, ncod, primitive_force, hfx_opts, &
&                 max_contraction, tmp_max, &
&                 offset_a, offset_b, offset_c, offset_d, ZetaInv, &
&                 EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                 rpq2, rab, rcd)

              endif ! sqrt(rpq2) .gt. (R_P + R_Q)

            endif ! .not. hfx_opts%farfield

            cart_estimate = max(tmp_max, cart_estimate)

          enddo ! lpgf
        enddo ! kpgf
      enddo ! jpgf
    enddo ! ipgf

    if ( cart_estimate > hfx_opts%eps_schwarz ) then

      do coord = 1,12

        T2 => primitive_force( &
&         (coord-1)*ncoa*ncob*ncoc*ncod+1:coord*ncoa*ncob*ncoc*ncod)

        call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
&         1.0_dp, T2(1), ncoa, sphia(1,1), ncoa, &
&         0.0_dp, T1(1), ncob*ncoc*ncod)

        call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
&         1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
&         0.0_dp, T2(1), nso(la)*ncoc*ncod)

        call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
&         1.0_dp, T2(1), ncoc, sphic(1,1), ncoc, &
&         0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

        call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
&         1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
&         0.0_dp, T2(1), nso(la)*nso(lb)*nso(lc))

        do dm=1,nso(ld)
          do cm=1,nso(lc)
            do bm=1,nso(lb)
              do am=1,nso(la)
                index_primitive_force = (dm-1)*nso(lc)*nso(lb)*nso(la) + &
&                                       (cm-1)*nso(lb)*nso(la) + &
&                                       (bm-1)*nso(la) + am
                eri_force(am,bm,cm,dm,coord) = T2(index_primitive_force)
              enddo ! am
            enddo ! bm
          enddo ! cm
        enddo ! dm

      enddo ! coord

    end if

    nullify(T2)
    call de_alloc(primitive_force, &
&     name='primitive_force', routine='calc_contract_deriv_eri')
    call de_alloc(work_forces, &
&     name='work_forces', routine='calc_contract_deriv_eri')
    call de_alloc(T1, &
&     name='T1', routine='calc_contract_deriv_eri')

  end subroutine calc_contract_deriv_eri

! *****************************************************************************
!> \brief Computes contractions of spherical Gaussians
!!
!! We shall make use of many expressions for Molecular Integrals Evaluation
!! taken from Ref. \cite Fermann:2020
!! \par History
!!      - 01.2010 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] cell: supercell in real space (units in Bohrs) 
!!                  first index:  cartesian index
!!                  second index: vector index
!! \param[in] rcell: supercell in reciprocal space (units in Bohrs^-1)
!!                  first index: cartesian index
!!                  second index: vector index
!! \param[in] ra: Cartesian position of the the center of the 
!!                first  atomic orbital (in Bohrs)
!! \param[in] rb: Cartesian position of the the center of the 
!!                second atomic orbital (in Bohrs)
!! \param[in] rc: Cartesian position of the the center of the 
!!                third atomic orbital (in Bohrs)
!! \param[in] rd: Cartesian position of the the center of the 
!!                fourth atomic orbital (in Bohrs)
!! \param[in] npgfa: Number of Gaussians in the expansion of the radial part
!!                   of the first numerical atomic orbital
!! \param[in] npgfb: Number of Gaussians in the expansion of the radial part
!!                   of the second numerical atomic orbital
!! \param[in] npgfc: Number of Gaussians in the expansion of the radial part
!!                   of the third numerical atomic orbital
!! \param[in] npgfd: Number of Gaussians in the expansion of the radial part
!!                   of the fourth numerical atomic orbital
!! \param[in] la: Angular momentum of the first atomic orbital
!! \param[in] lb: Angular momentum of the second atomic orbital
!! \param[in] lc: Angular momentum of the third atomic orbital
!! \param[in] ld: Angular momentum of the fourth atomic orbital
!! \param[in] ncoa:  Total number of Cartesian Gaussians 
!!                   required to expand a NAO of a
!!                   given angular momentum shell.
!!                   For a s-shell = npgfa.
!!                   For a p-shell = npgfa * 3, 
!!                   since to expand the s-shell we 
!!                   need three cartesian functions.
!!                   For a d-shell = npgfa * 6
!!                   since to expand the d-shell we
!!                   need six cartesian functions.
!!                   For a f-shell = npgfa * 10
!!                   since to expand the f-shell we
!!                   need six cartesian functions.
!! \param[in] ncob: Same as ncoa for the second atomic orbital
!! \param[in] ncoc: Same as ncoa for the third atomic orbital
!! \param[in] ncod: Same as ncoa for the fourth atomic orbital
!! \param[in] zeta: Exponents of the Gaussians required
!!                  to expand the radial part of the
!!                  first numerical atomic orbital
!! \param[in] zetb: Same as zeta for the second atomic orbital
!! \param[in] zetc: Same as zeta for the third atomic orbital
!! \param[in] zetd: Same as zeta for the fourth atomic orbital
!! \param[in] sphia: For the first NAO, with angular momentum la, 
!!            and magnetic quantum number ranging between 1 and nso(la) 
!!            [second index of the array]
!!            sphia stores the coefficients of all the Cartesian Gaussian 
!!            functions required to expand that NAO
!! \param[in] sphib: Same as sphia but for the second NAO
!! \param[in] sphic: Same as sphia but for the third NAO
!! \param[in] sphid: Same as sphia but for the fourth NAO
!! \param[in] hfx_opts: data structure containing Hartree-Fok exchange
!!                           parameters
!! \param[out] eri: ...
! *****************************************************************************
  subroutine calc_contract_eri2(lib, cell, rcell, ra, rb, rc, rd, &
      npgfa, npgfb, npgfc, npgfd, la, lb, lc, ld, ncoa, ncob, ncoc, ncod, &
      zeta, zetb, zetc, zetd, sphia, sphib, sphic, sphid, hfx_opts, eri)

    use alloc
    use nao2gto_data, only: nso, nco
    use nao2gto_libint
    use nao2gto_pbc
    use nao2gto_primitive
    use nao2gto_types

!   For debugging
#ifdef MPI
    use mpi_siesta
#endif
 
!   End debugging

    implicit none

!   For debugging
#ifdef MPI
    integer     :: MPIerror
#endif
!   End debugging


!   Arguments
    type(Libint_t), intent(inout) :: lib
                                        ! Libint data structure
    integer,  intent(in) :: npgfa       ! Number of Gaussians in the expansion 
                                        !    of the radial part of the first 
                                        !    numerical atomic orbital
    integer,  intent(in) :: npgfb       ! Number of Gaussians in the expansion 
                                        !    of the radial part of the second
                                        !    numerical atomic orbital
    integer,  intent(in) :: npgfc       ! Number of Gaussians in the expansion 
                                        !    of the radial part of the third
                                        !    numerical atomic orbital
    integer,  intent(in) :: npgfd       ! Number of Gaussians in the expansion 
                                        !    of the radial part of the fourth
                                        !    numerical atomic orbital
    integer,  intent(in) :: la          ! Angular momentum of the first 
                                        !    atomic orbital
    integer,  intent(in) :: lb          ! Angular momentum of the second
                                        !    atomic orbital
    integer,  intent(in) :: lc          ! Angular momentum of the third
                                        !    atomic orbital
    integer,  intent(in) :: ld          ! Angular momentum of the fourth
                                        !    atomic orbital
    integer,  intent(in) :: ncoa        ! Total number of Cartesian Gaussians 
                                        !    required to expand a NAO of a 
                                        !    given angular momentum shell
                                        !    For a s-shell = npgfa
                                        !    For a p-shell = npgfa * 3, 
                                        !    since to expand the s-shell we 
                                        !    need three cartesian functions
                                        !    For a d-shell = npgfa * 6
                                        !    since to expand the d-shell we 
                                        !    need six cartesian functions
                                        !    For a f-shell = npgfa * 10
                                        !    since to expand the f-shell we 
                                        !    need six cartesian functions
    integer,  intent(in) :: ncob        ! Same as ncoa for the second orbital
    integer,  intent(in) :: ncoc        ! Same as ncoa for the third orbital
    integer,  intent(in) :: ncod        ! Same as ncoa for the fourth orbital
    real(dp), intent(in) :: cell(3,3)   ! Supercell lattice vectors in 
                                        !    real space (in Bohrs)
                                        !    First index:  cartesian index
                                        !    Second index: vector index
    real(dp), intent(in) :: rcell(3,3)  ! Reciprocal lattice vectors of the 
                                        !    supercell (in Bohrs^-1)
                                        !    The factor of 2*pi is not included
                                        !    First index:  cartesian index
                                        !    Second index: vector index
    real(dp), intent(in) :: ra(3)       ! Cartesian position of the center of 
                                        !    the first  atomic orbital
    real(dp), intent(in) :: rb(3)       ! Cartesian position of the center of 
                                        !   the  second atomic orbital
    real(dp), intent(in) :: rc(3)       ! Cartesian position of the center of 
                                        !    the third  atomic orbital
    real(dp), intent(in) :: rd(3)       ! Cartesian position of the center of
                                        !    the fourth atomic orbital
                                        !    All these cartesian positions 
                                        !    in Bohrs
    real(dp), intent(in) :: zeta(npgfa) ! Exponents of the Gaussians required
                                        !    to expand the radial part of the
                                        !    first numerical atomic orbital
    real(dp), intent(in) :: zetb(npgfb) ! Exponents of the Gaussians required
                                        !    to expand the radial part of the
                                        !    second numerical atomic orbital
    real(dp), intent(in) :: zetc(npgfc) ! Exponents of the Gaussians required
                                        !    to expand the radial part of the
                                        !    third numerical atomic orbital
    real(dp), intent(in) :: zetd(npgfd) ! Exponents of the Gaussians required
                                        !    to expand the radial part of the
                                        !    fourth numerical atomic orbital
    real(dp), intent(in) :: sphia(ncoa,nso(la)) 
                                        ! For the first NAO, 
                                        !    with angular momentum la, 
                                        !    and magnetic quantum number ranging
                                        !    between 1 and nso(la) 
                                        !    [second index of the array]
                                        !    sphia stores the coefficients
                                        !    of all the Cartesian Gaussian 
                                        !    functions required to expand 
                                        !    that NAO
    real(dp), intent(in) :: sphib(ncob,nso(lb)) 
                                        ! Same as sphia, but for the second
                                        !    NAO
    real(dp), intent(in) :: sphic(ncoc,nso(lc)) 
                                        ! Same as sphia, but for the third
                                        !    NAO
                                        ! 
    real(dp), intent(in) :: sphid(ncod,nso(ld)) 
                                        ! Same as sphia, but for the fourth
                                        !    NAO
    type(hfx_options_type),intent(in) :: hfx_opts  
                                        ! Data structure storing Hartree-Fock 
                                        !    exchange options
    real(dp), intent(out) :: eri(nso(la),nso(lb),nso(lc),nso(ld))

!   Local variables
    integer :: ipgf       ! Counter for a loop over all the 
                          !    cartesian Gaussian functions required
                          !    to expand the first NAO
    integer :: jpgf       ! Same as jpgf for the second NAO
    integer :: kpgf       ! Same as kpgf for the second NAO
    integer :: lpgf       ! Same as lpgf for the second NAO
    integer :: offset_a   ! Index to order the primitive integrals 
                          !    between the primitive Cartesian Gaussian 
                          !    functions for the first atomic orbital.
                          !    If we assume the same number of Gaussians
                          !    in the expansion of all the radial NAO,
                          !    then
                          !    For a s-shell 
                          !      (=> One Cartessian Gaussian function)
                          !      offset_a = 0 for the first exponential
                          !      offset_a = 1 for the second exponential
                          !      offset_a = NG-1 for the NG exponential
                          !    For a p-shell 
                          !      (=> Three Cartessian Gaussian function)
                          !      offset_a=0 for the px of the first exponential
                          !      offset_a=1 for the py of the first exponential
                          !      offset_a=2 for the pz of the first exponential
                          !      offset_a=3 for the px of the second exponential
                          !      offset_a=4 for the py of the second exponential
                          !      offset_a=5 for the pz of the second exponential
                          !      ...
                          !      offset_a=3*NG-3 for the px of the NG exponenti 
                          !      offset_a=3*NG-2 for the py of the NG exponenti 
                          !      offset_a=3*NG-1 for the pz of the NG exponenti 
                          !    For a d-shell 
                          !      (=> Six Cartessian Gaussian function)
                          !      offset_a=0 for the dxx of the first exponential
                          !      offset_a=1 for the dxy of the first exponential
                          !      offset_a=2 for the dxz of the first exponential
                          !      offset_a=3 for the dyy of the first exponential
                          !      offset_a=4 for the dyz of the second exponentia
                          !      offset_a=5 for the dzz of the second exponentia
                          !      ...
                          !      offset_a=6*NG-6 for the dxx of the NG exponenti
                          !      offset_a=6*NG-5 for the dxy of the NG exponenti
                          !      offset_a=6*NG-4 for the dxz of the NG exponenti
                          !      offset_a=6*NG-3 for the dyy of the NG exponenti
                          !      offset_a=6*NG-2 for the dyz of the NG exponenti
                          !      offset_a=6*NG-1 for the dzz of the NG exponenti
    integer :: offset_b   ! Same as offset_a for the second NAO
    integer :: offset_c   ! Same as offset_a for the third NAO
    integer :: offset_d   ! Same as offset_a for the fourth NAO

    integer ::  am, bm, cm, dm, ieri, i, j, k, l, index_primitive_integrals
    real(dp) :: rab(3)    ! Relative position between the centers of the
                          !   first and the second NAO
    real(dp) :: rcd(3)    ! Relative position between the centers of the
                          !   third and the fourth NAO
    real(dp) :: rab2      ! Square of the distance between the centers of the 
                          !   first and the second NAO
    real(dp) :: rcd2      ! Square of the distance between the centers of the 
                          !   third and the fourth NAO
    real(dp) :: Zeta_A    ! Exponent of the first  primitive Cartesian Gaussian
    real(dp) :: Zeta_B    ! Exponent of the second primitive Cartesian Gaussian
    real(dp) :: Zeta_C    ! Exponent of the third  primitive Cartesian Gaussian
    real(dp) :: Zeta_D    ! Exponent of the fourth primitive Cartesian Gaussian
    real(dp) :: Zeta1     ! Zeta_A + Zeta_B 
                          !    In the LIBINT Manual, this variable is called 
                          !    zeta.
                          !    In the notes by Fermann et al, Eq. (2.33),
                          !    it is also referred to as \alpha_p
    real(dp) :: ZetaInv   ! Inverse of Zeta
    real(dp) :: S1234a    ! Eq. (15) of the LIBINT Manual:
                          !   \f$ \left[ \left( \frac{\pi}{\zeta} \right)^{3/2}
                          !   \exp \left( - \frac{\zeta_{a} \zeta_{b}}{\zeta}
                          !   \vert \vec{AB} \vert^{2} \right) \right]
    real(dp) :: P(3)      ! Center of the Gaussian that results from the 
                          !    product of the first and second primitive
                          !    Cartesian Gaussian functions
    real(dp) :: Q(3)      ! Center of the Gaussian that results from the 
                          !    product of the third and fourth primitive 
                          !    Cartesian Gaussian functions
    real(dp) :: ZetapEtaInv 
                          ! 1.0/(Zeta + Eta) following the notation of 
                          !    LIBINT Manual.
                          ! 1.0/(gamma_p + \gamma_q) in the notation of 
                          !    Eq. (6.4) of Ref. Fermann:2020
    real(dp) :: Rho       ! [(\zeta * \eta)/(\zeta + \eta)] in the notation of
                          !    LIBINT Manual.
                          ! Or [(\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)]
                          !    that appears in Eq. (6.4) of Ref. Fermann:2020
    real(dp) :: RhoInv    ! Inverse of Rho
    real(dp) :: rpq2      ! The variable rpq2 corresponds with 
                          !   $\overline{\vec{P}\vec{Q}}^{2}$
                          !   in Eq. (6.4) of Ref. Fermann:2020
    real(dp) :: Eta       ! Zeta_C + Zeta_D 
                          !    In the LIBINT Manual, this variable is called 
                          !    eta.
                          !    In the notes by Fermann et al, Eq. (2.33),
                          !    it is also referred to as \alpha_q
    real(dp) :: EtaInv    ! Inverse of Eta
    real(dp) :: S1234     ! Product of the Eq. (15) and (16) of LIBINT Manual.
                          !    \f$\left[\exp\left(-\frac{\zeta_{a}\zeta_{b}}
                          !     {\zeta} \vert \vec{AB} \vert^{2} \right) \right]
                          !     \times \left[ 
                          !     \exp \left( - \frac{\zeta_{c} \zeta_{d}}{\eta}
                          !     \vert \vec{CD} \vert^{2} \right) \right]    \f$
    real(dp) :: W(3)      ! Eq. (12) of LIBINT Manual.
                          !     \f$ \vec{W}=\frac{\zeta \vec{P} + \eta \vec{Q}}
                          !     {\zeta + \eta} \f$
    real(dp) :: tmp_max, r_temp(3), r_pbc_temp(3), &
&     alpha_P, alpha_Q, R_P, R_Q, Kab, Kcd, theta_w, rc_trans(3), rd_trans(3)

    real(dp), pointer :: primitive_integrals(:) => null()
    real(dp), pointer :: T1(:) => null()

! -------------------------------------------------------------------------

    call re_alloc(primitive_integrals, 1, ncoa*ncob*ncoc*ncod, &
&     name='primitive_integrals', routine='calc_contract_eri2')
    primitive_integrals(:) = 0.0_dp

    call re_alloc(T1, 1, ncoa*ncob*ncoc*ncod, &
&     name='T1', routine='calc_contract_eri2')
    T1(:) = 0.0_dp

!!   For debugging
!    do ipgf = 1, 3
!      write(6,'(a,i5,3f12.5)')     &
! &     'calc_contract_eri2: cell,  ivector = ', ipgf, cell(ipgf,1:3)
!    enddo 
!    do ipgf = 1, 3
!      write(6,'(a,i5,3f12.5)')     &
! &     'calc_contract_eri2: rcell, ivector = ', ipgf, rcell(ipgf,1:3)
!    enddo 
!    write(6,'(a,3f12.5)')     &
! &   'calc_contract_eri2: ra = ', ra(:)
!    write(6,'(a,3f12.5)')     &
! &   'calc_contract_eri2: rb = ', rb(:)
!    write(6,'(a,3f12.5)')     &
! &   'calc_contract_eri2: rc = ', rc(:)
!    write(6,'(a,3f12.5)')     &
! &   'calc_contract_eri2: rd = ', rd(:)
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: npgfa = ', npgfa
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: npgfb = ', npgfb
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: npgfc = ', npgfc
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: npgfd = ', npgfd
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: la = ', la
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: lb = ', lb
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: lc = ', lc
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ld = ', ld
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ncoa = ', ncoa
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ncob = ', ncob
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ncoc = ', ncoc
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ncod = ', ncod
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: nso(la) = ', nso(la)
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: nso(lb) = ', nso(lb)
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: nso(lc) = ', nso(lc)
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: nso(ld) = ', nso(ld)
!    write(6,'(a,i5)')     &
! &   'calc_contract_eri2: ncoa*ncob*ncoc*ncod = ', ncoa*ncob*ncoc*ncod
!!   End debugging

!   Compute the relative position between the centers of the 
!   first and the second NAO
    rab(1:3) = ra(1:3) - rb(1:3)
!   Compute the relative position between the centers of the 
!   third and the fourth NAO
    rcd(1:3) = rc(1:3) - rd(1:3)
!   Square of the distance between the centers of the first and the second NAO
    rab2 = rab(1)**2 + rab(2)**2 + rab(3)**2
!   Square of the distance between the centers of the third and the fourth NAO
    rcd2 = rcd(1)**2 + rcd(2)**2 + rcd(3)**2

!   Loop over all the cartesian Gaussian functions required to expand the 
!   first NAO
    do ipgf = 1, npgfa
!     This primitive Gaussian is used in the expansion of the radial part of 
!     a NAO. This radial part depends on the angular momentum, but not on
!     the magnetic quantum number.
!     In other words, all the numerical atomic orbitals with the same l, but
!     different m will share the same radial part.
!     Here we leave space to accomodate all the Cartesian Gaussian functions
!     where this primitive Gaussian enters
      offset_a = (ipgf-1)*nco(la)
      Zeta_A = zeta(ipgf)

!!     For debugging
!      write(6,'(a,2i5,f12.5)')                              &
! &      'calc_contract_eri2: ipgf, offset_a, zeta_a = ',    &
! &      ipgf, offset_a, zeta_a
!!     End debugging

!     Loop over all the cartesian Gaussian functions required to expand the 
!     second NAO
      do jpgf = 1, npgfb
        offset_b = (jpgf-1)*nco(lb)
        Zeta_B = zetb(jpgf)
!!       For debugging
!        write(6,'(a,2i5,f12.5)')                              &
! &        'calc_contract_eri2: jpgf, offset_b, zeta_b = ',    &
! &        jpgf, offset_b, zeta_b
!!       End debugging

!       Compute the sum of exponents of the first and second Gaussian functions,
!       equivalent to \zeta in the previous doxygen documentation
!       also called (\gamma) in Eq. (2.33) of the notes by Fermann, 
!       Ref. \cite Fermann:2020,
!       on the Gaussian Product Theorem
        Zeta1 = Zeta_A + Zeta_B
!       Compute the center of the Gaussian that results from the product of
!       the first and second Gaussians,
!       equivalent to \vec{P} in the previous doxygen documentation
        ZetaInv = 1.0_dp/Zeta1
        P(1:3) = (Zeta_A*ra(1:3) + Zeta_B*rb(1:3))*ZetaInv
!       Compute the argument that appears in the first exponential 
!       in the expression of the Gaussian Product Theorem
!       See doxygen documentation
        S1234a = (-Zeta_A*Zeta_B*ZetaInv*rab2)

!       Loop over all the cartesian Gaussian functions required to expand the 
!       third NAO
        do kpgf=1,npgfc
          offset_c = (kpgf-1)*nco(lc)
          Zeta_C = zetc(kpgf)
!!         For debugging
!          write(6,'(a,2i5,f12.5)')                              &
! &          'calc_contract_eri2: kpgf, offset_c, zeta_c = ',    &
! &          kpgf, offset_c, zeta_c
!!         End debugging

!         Loop over all the cartesian Gaussian functions required to expand the 
!         fourth NAO
          do lpgf=1,npgfd
            offset_d = (lpgf-1)*nco(ld)
            Zeta_D = zetd(lpgf)
!!           For debugging
!            write(6,'(a,2i5,f12.5)')                              &
! &            'calc_contract_eri2: lpgf, offset_d, zeta_d = ',    &
! &            lpgf, offset_d, zeta_d
!!           End debugging

!           Compute the sum of exponents of the 
!           third and fourth Gaussian functions,
!           equivalent to \gamma in the previous doxygen documentation
!           on the Gaussian Product Theorem
            Eta = Zeta_C + Zeta_D
!           Compute the center of the Gaussian that results from the product of
!           the first and second Gaussians,
!           equivalent to \vec{P} in the previous doxygen documentation
            EtaInv = 1.0_dp/Eta
            Q(1:3) = (Zeta_C*rc(1:3) + Zeta_D*rd(1:3))*EtaInv
!           Compute the argument that appears in the first exponential 
!           in the expression of the Gaussian Product Theorem
!           See doxygen documentation
            S1234 = exp(S1234a-Zeta_C*Zeta_D*EtaInv*rcd2)

!           We compute the term (\gamma_p + \gamma_q) that appears in 
!           Eq. (6.4) of Ref. Fermann:2020
            ZetapEtaInv = Zeta1+Eta
            ZetapEtaInv = 1.0_dp/ZetapEtaInv

!           We compute the term [(\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)]
!           that appears in Eq. (6.4) of Ref. Fermann:2020
            Rho = Zeta1 * Eta * ZetapEtaInv
!           We invert Rho
            RhoInv = 1.0_dp / Rho

!           Compute the square of the distances between the centers
!           of the two Gaussians that come from:
!           The product of the first and second primitive Gaussian functions,
!           (\vec{P})
!           the product of the third and fourth primitive Gaussian functions,
!           (\vec{Q})
!           The variable rpq2 corresponds with $\overline{\vec{P}\vec{Q}}^{2}$
!           in Eq. (6.4) of Ref. Fermann:2020
            r_temp(1:3) = Q(1:3)-P(1:3)
            call trans_pbc(r_temp, cell, rcell, r_pbc_temp)

            Q(1:3) = P(1:3)+r_pbc_temp(1:3)
            rc_trans(1:3) = rc(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rd_trans(1:3) = rd(1:3) + r_pbc_temp(1:3) - r_temp(1:3)
            rpq2 = (P(1)-Q(1))**2 + (P(2)-Q(2))**2 + (P(3)-Q(3))**2

            W(1:3) = (Zeta1*P(1:3) + Eta*Q(1:3))*ZetapEtaInv

            call calc_primitive_eri2(lib, ra, rb, rc_trans, rd_trans, &
&             zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf), la, lb, lc ,ld, &
&             ncoa, ncob, ncoc, ncod, offset_a, offset_b, offset_c, offset_d, &
&             primitive_integrals, hfx_opts, Zeta1, ZetaInv, Eta, EtaInv, &
&             ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, rpq2, rab, rcd)

          enddo
        enddo
      enddo
    enddo


    call dgemm("T", "N", ncob*ncoc*ncod, nso(la), ncoa, &
      1.0_dp, primitive_integrals(1), ncoa, sphia(1,1), ncoa, &
      0.0_dp, T1(1), ncob*ncoc*ncod)

    call dgemm("T", "N", nso(la)*ncoc*ncod, nso(lb), ncob, &
      1.0_dp, T1(1), ncob, sphib(1,1), ncob, &
      0.0_dp, primitive_integrals(1), nso(la)*ncoc*ncod)

    call dgemm("T", "N", nso(la)*nso(lb)*ncod, nso(lc), ncoc, &
      1.0_dp, primitive_integrals(1), ncoc, sphic(1,1), ncoc, &
      0.0_dp, T1(1), nso(la)*nso(lb)*ncod)

    call dgemm("T", "N", nso(la)*nso(lb)*nso(lc), nso(ld), ncod, &
      1.0_dp, T1(1), ncod, sphid(1,1), ncod, &
      0.0_dp, primitive_integrals(1), nso(la)*nso(lb)*nso(lc))

    do dm=1,nso(ld)
      do cm=1,nso(lc)
        do bm=1,nso(lb)
          do am=1,nso(la)
            index_primitive_integrals = (dm-1)*nso(lc)*nso(lb)*nso(la) + &
              (cm-1)*nso(lb)*nso(la) + (bm-1)*nso(la) + am
            eri(am,bm,cm,dm) = primitive_integrals(index_primitive_integrals)
!!           For debugging
!            write(6,'(a,5i5,f20.12)') &
! &            'calc_contract_eri2: dm, cm, bm, am, index_primitive_integrals, eri = ' , &
! &            dm, cm, bm, am, index_primitive_integrals, eri(am,bm,cm,dm)
!!           End debugging
          enddo
        enddo
      enddo
    enddo

    call de_alloc(primitive_integrals, name='primitive_integrals', &
      routine='calc_contract_eri2')
    call de_alloc(T1, name='T1', routine='calc_contract_eri2')

!!   For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!    call die()
!!   End debugging


  end subroutine calc_contract_eri2

end module nao2gto_contract
