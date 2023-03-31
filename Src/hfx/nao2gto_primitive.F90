! *** Module: nao2gto_primitive ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Interface to the Libint and Libderiv libraries
!!
!!  This module was originally extracted from CP2K: A general program to
!!  perform molecular dynamics simulations. It was then adapted for use within
!!  SIESTA.
!!
!! \note
!!      This file currently works with a version of Libint configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2000-2009 CP2K Developers Group
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2006 Created [Manuel Guidon]
!!      - 01.2018 Refactored for SIESTA [Yann Pouillon]
! *****************************************************************************
module nao2gto_primitive

  implicit none

  private

  integer, parameter :: full_perm1(12) = [1,2,3,4,5,6,7,8,9,10,11,12]
  integer, parameter :: full_perm2(12) = [4,5,6,1,2,3,7,8,9,10,11,12]
  integer, parameter :: full_perm3(12) = [1,2,3,4,5,6,10,11,12,7,8,9]
  integer, parameter :: full_perm4(12) = [4,5,6,1,2,3,10,11,12,7,8,9]
  integer, parameter :: full_perm5(12) = [7,8,9,10,11,12,1,2,3,4,5,6]
  integer, parameter :: full_perm6(12) = [7,8,9,10,11,12,4,5,6,1,2,3]
  integer, parameter :: full_perm7(12) = [10,11,12,7,8,9,1,2,3,4,5,6]
  integer, parameter :: full_perm8(12) = [10,11,12,7,8,9,4,5,6,1,2,3]

  public :: &
&   calc_primitive_deriv_eri, &
&   calc_primitive_eri, &
&   calc_primitive_eri2, &
&   calc_primitive_screen

contains

! *****************************************************************************
!> \brief Evaluates derivatives of electron repulsion integrals for a primitive
!!        quartet.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 08.2007 Refactured permutation part [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] deriv_data: Libderiv data structure
!! \param[in] A: ...
!! \param[in] B: ...
!! \param[in] C: ...
!! \param[in] D: ...
!! \param[in] ZetaA: ...
!! \param[in] ZetaB: ...
!! \param[in] ZetaC: ...
!! \param[in] ZetaD: ...
!! \param[in] n_a: ...
!! \param[in] n_b: ...
!! \param[in] n_c: ...
!! \param[in] n_d: ...
!! \param[in,out] work_forces: ...
!! \param[in] ncoa: ...
!! \param[in] ncob: ...
!! \param[in] ncoc: ...
!! \param[in] ncod: ...
!! \param[out] primitive_forces: ...
!! \param[in] hfx_opts: ...
!! \param[in] offset_a: ...
!! \param[in] offset_b: ...
!! \param[in] offset_c: ...
!! \param[in] offset_d: ...
!! \param[in] ZetaInv: ...
!! \param[in] EtaInv: ...
!! \param[in] ZetapEtaInv: ...
!! \param[in] Rho: ...
!! \param[in] RhoInv: ...
!! \param[in] S1234: ...
!! \param[in] P: ...
!! \param[in] Q: ...
!! \param[in] W: ...
!! \param[in] rpq2: ...
!! \param[in] AB: ...
!! \param[in] CD: ...
! *****************************************************************************
  subroutine calc_primitive_deriv_eri( &
&                deriv_data, A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&                n_a, n_b, n_c, n_d, work_forces, ncoa, ncob, ncoc, ncod, &
&                primitive_forces, hfx_opts, max_contraction, tmp_max_all, &
&                offset_a, offset_b, offset_c,offset_d, &
&                ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
&                P, Q, W, rpq2, AB, CD)

    use sys, only: die
    use precision, only: dp
    use nao2gto_data, only: nco
    use nao2gto_libint
    use nao2gto_types, only: hfx_options_type
    use nao2gto_wrappers, only: get_derivs

    implicit none

    ! Arguments
    type(Libderiv_t), intent(inout)       :: deriv_data
    real(dp), intent(in)                  :: A(3), B(3), C(3), D(3), &
&                                            Zeta_A, Zeta_B, Zeta_C, Zeta_D
    integer, intent(in)                   :: n_a, n_b, n_c, n_d
    real(dp),&
&     dimension(nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d), 12), &
&     intent(inout)                       :: work_forces
    integer, intent(in)                   :: ncoa, ncob, ncoc, ncod
    real(dp), &
&     dimension(ncoa, ncob, ncoc, ncod, 12), & 
&     intent(out)                         :: primitive_forces
    real(dp), intent(in) :: max_contraction
    real(dp), intent(out) :: tmp_max_all
    integer, intent(in)                   :: offset_a, offset_b, offset_c, &
&                                            offset_d
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in)                  :: ZetaInv, EtaInv, ZetapEtaInv, &
&                                            Rho, RhoInv, S1234, P(3), Q(3),&
&                                            W(3), rpq2, AB(3), CD(3)

    ! Local variables
    character(len=120)      :: msg
    integer                 :: a_mysize(1), i, j, k, l, m_max, mysize, n, &
&                              p1, p2, p3, perm_case
    logical                 :: do_it
    real(dp) :: tmp_max
    type(prim_data), target :: prim

    ! -------------------------------------------------------------------------

    ! Permutation of configurations
    perm_case = 1
    if(n_a<n_b) then
      perm_case = perm_case + 1
    end if
    if(n_c<n_d) then
      perm_case = perm_case + 2
    end if
    if( n_a+n_b > n_c+n_d) then
      perm_case = perm_case + 4
    end if

    do_it = .true.
    m_max = n_a+n_b+n_c+n_d
!   Since we are going to compute first derivatives of the ERIS,
!   we have to increase m_max by 1.
!   See page 7 of LIBINT Manual.
    m_max = m_max + 1
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize=mysize

    select case(perm_case)

      case(1)
        call build_deriv_data(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = AB ! A-B
        deriv_data%CD = CD ! C-D
        call get_derivs(n_d, n_c, n_b, n_a, deriv_data, prim, work_forces, a_mysize)
        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                do l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm1(n)) = work_forces(p3+l,n)
                end do
              end do
            end do
          end do
        end do

      case(2)
        call build_deriv_data(B, A, C, D, Zeta_B, Zeta_A, Zeta_C, Zeta_D, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = -AB ! B-A
        deriv_data%CD =  CD ! C-D
        call get_derivs(n_d, n_c, n_a, n_b, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                do l = 1,nco(n_d)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm2(n)) = work_forces(p3+l,n)
                end do
              end do
            end do
          end do
        end do

     case(3)
        call build_deriv_data(A, B, D, C, Zeta_A, Zeta_B, Zeta_D, Zeta_C, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB =  AB ! A-B
        deriv_data%CD = -CD ! D-C
        call get_derivs(n_c, n_d, n_b, n_a, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp* (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_d)
              do l = 1,nco(n_d)
                p3 = (p2 + l-1) * nco(n_c)
                do k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm3(n)) = work_forces(p3+k,n)
                end do
              end do
            end do
          end do
        end do

      case(4)
        call build_deriv_data(B, A, D, C, Zeta_B, Zeta_A, Zeta_D, Zeta_C, &
&         m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, P, Q, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = -AB ! B-A
        deriv_data%CD = -CD ! D-C
        call get_derivs(n_c, n_d, n_a, n_b, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              do l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                do k = 1,nco(n_c)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm4(n)) = work_forces(p3+k,n)
                end do
              end do
            end do
          end do
        end do

      case(5)
        call build_deriv_data(C, D, A, B, Zeta_C, Zeta_D, Zeta_A, Zeta_B, &
&          m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&          Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB = CD ! C-D
        deriv_data%CD = AB ! A-B
        call get_derivs(n_b, n_a, n_d, n_c, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_a)
              do i = 1,nco(n_a)
                p3 = (p2 + i-1)*nco(n_b)
                do j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm5(n)) = work_forces(p3+j,n)
                end do
              end do
            end do
          end do
        end do

      case(6)
        call build_deriv_data(C, D, B, A, Zeta_C, Zeta_D, Zeta_B, Zeta_A, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB =  CD ! C-D
        deriv_data%CD = -AB ! B-A
        call get_derivs(n_a, n_b, n_d, n_c, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
          do i=1,mysize
             work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                          work_forces(i,k+3) + &
&                                          work_forces(i,k+6))
          end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              p2 = (p1 + l-1)*nco(n_b)
              do j = 1,nco(n_b)
                p3 = (p2 + j-1)*nco(n_a)
                do i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm6(n)) = work_forces(p3+i,n)
                end do
              end do
            end do
          end do
        end do

      case(7)
        call build_deriv_data(D, C, A, B, Zeta_D, Zeta_C, Zeta_A, Zeta_B, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB=-CD ! D-C
        deriv_data%CD=AB ! A-B
        call get_derivs(n_b, n_a, n_c, n_d, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              do i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                do j = 1,nco(n_b)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm7(n)) = work_forces(p3+j,n)
                end do
              end do
            end do
          end do
        end do

      case(8)
        call build_deriv_data(D, C, B, A, Zeta_D, Zeta_C, Zeta_B, Zeta_A, &
&         m_max, hfx_opts, prim, do_it, EtaInv, ZetaInv, ZetapEtaInv, &
&         Rho, RhoInv, S1234, Q, P, W, rpq2)

        if( .not. do_it ) return

        deriv_data%AB=-CD ! D-C
        deriv_data%CD=-AB ! B-A
        call get_derivs(n_a, n_b, n_c, n_d, deriv_data, prim, work_forces, a_mysize)

        do k=4,6
           do i=1,mysize
              work_forces(i,k) = -1.0_dp * (work_forces(i,k-3) + &
&                                           work_forces(i,k+3) + &
&                                           work_forces(i,k+6))
           end do
        end do
        do n=1,12
          tmp_max = 0.0_dp
          do i=1,mysize
            tmp_max = max(tmp_max, abs(work_forces(i,n)))
          end do
          tmp_max = tmp_max*max_contraction
          tmp_max_all = max(tmp_max_all, tmp_max)
          if ( tmp_max < hfx_opts%eps_schwarz ) cycle

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_b)
              do j = 1,nco(n_b)
                p3 = (p2 + j-1) * nco(n_a)
                do i = 1,nco(n_a)
                  primitive_forces(offset_a+i,offset_b+j,offset_c+k,offset_d+l, &
&                   full_perm8(n)) = work_forces(p3+i,n)
                end do
              end do
            end do
          end do
        end do

      case default
        write(msg, fmt='(A,": ",I2)') "Invalid HFX permutation case", &
          perm_case
        call die(msg)

    end select

  end subroutine calc_primitive_deriv_eri

! *****************************************************************************
!> \brief Evaluates electron repulsion integrals for a primitive quartet.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 11.2006 Created [Manuel Guidon]
!!      - 08.2007 Refactured permutation part [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] lib: Libint data structure
!! \param[in] A: Center of the first  Cartesian Gaussian function
!! \param[in] B: Center of the second Cartesian Gaussian function
!! \param[in] C: Center of the third  Cartesian Gaussian function 
!! \param[in] D: Center of the fourth Cartesian Gaussian function 
!! \param[in] n_a: Angular momentum of the first  primitive Cartesian Gaussian
!! \param[in] n_b: Angular momentum of the second primitive Cartesian Gaussian
!! \param[in] n_c: Angular momentum of the third  primitive Cartesian Gaussian
!! \param[in] n_d: Angular momentum of the fourth primitive Cartesian Gaussian
!! \param[in] ncoa: Total number of Cartesian Gaussians
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
!! \param[in] offset_a: ...
!! \param[in] offset_b: ...
!! \param[in] offset_c: ...
!! \param[in] offset_d: ...
!! \param[out] primitives: ...
!! \param[in] hfx_opts: data structure containing Hartree-Fock exchange
!!                           parameters
!! \param[in] ZetaInv: \f$ 1.0/ ( \zeta_a + \zeta_b ) \f$
!!            where \f$ \zeta_{a} \f$ and \f$ \zeta_{b} \f$ are the exponents
!!            of the first and second primitive Cartesian Gaussian functions
!! \param[in] EtaInv: \f$ 1.0/ ( \zeta_c + \zeta_d ) \f$
!!            where \f$ \zeta_{c} \f$ and \f$ \zeta_{d} \f$ are the exponents
!!            of the third and fourth primitive Cartesian Gaussian functions
!! \param[in] ZetapEtaInv: \f$ 1.0/(\zeta + \eta) \f$ in the notation of
!!                   LIBINT Manual, where \f$ \zeta = \zeta_{a} + \zeta_{b} \f$,
!!                   and \f$ \eta = \zeta_{c} + \zeta_{b} \f$.
!!                   Or \f$ 1.0/(\gamma_p + \gamma_q) \f$ in the notation of
!!                   Eq. (6.4) of
!!                   Ref. \cite Fermann:2020
!! \param[in] Rho: \f$ (\zeta \times \eta) / (\zeta+\eta) \f$ in the notation of
!!                   LIBINT Manual.
!!                   Or \f$ (\gamma_p \times \gamma_q) / (\gamma_p+\gamma_q) \f$
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] RhoInv: \f$ (\zeta+\eta)/(\zeta \times \eta) \f$ in the notation
!!                   of LIBINT Manual.
!!                   Or \f$ (\gamma_p+\gamma_q)/(\gamma_p \times \gamma_q) \f$
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] S1234: Product of \f$ K_{1}\times K_{2} \f$ in Eq. (6.4) of
!!                   Ref. \cite Fermann:2020, also specified in
!!                   Eq. (15) and Eq. (16) of
!!                   LIBINT Manual (here called \f$S_{12}\f$ and \f$S_{34}\f$).
!!                   \f$ \left[ \exp \left( - \frac{\zeta_{a} \zeta_{b}}{\zeta}
!!                        \vert \vec{AB} \vert^{2} \right) \right]
!!                   \times \left[
!!                   \exp \left( - \frac{\zeta_{c} \zeta_{d}}{\eta}
!!                        \vert \vec{CD} \vert^{2} \right) \right]    \f$
!! \param[in] P: \f$ \vec{P} = \frac{\zeta_{a} \vec{A} + \zeta_{b} \vec{B}}
!!                             {\zeta} \f$. See Eq. (10) of LIBINT Manual.
!! \param[in] Q: \f$ \vec{Q} = \frac{\zeta_{c} \vec{C} + \zeta_{d} \vec{D}}
!!                             {\eta} \f$. See Eq. (11) of LIBINT Manual.
!! \param[in] W: \f$ \vec{W} = \frac{\zeta \vec{P} + \eta \vec{Q}}
!!                             {\zeta + \eta} \f$. See Eq. (12) of
!!                             LIBINT Manual.
!! \param[in] rpq2: The variable PQ2 corresponds with
!!                 \f$ \overline{\vec{PQ}}^{2} \f$
!!                 in Eq. (6.4) of Ref. \cite Fermann:2020, 
!!                 or inside the argument of the incomplete Gamma function
!!                 in Eq. (14) of LIBINT Manual
!! \param[in] AB: \f$ \vec{A} - \vec{B} \f$
!! \param[in] CD: \f$ \vec{C} - \vec{D} \f$
! *****************************************************************************
  subroutine calc_primitive_eri(lib, A, B, C, D, n_a, n_b, n_c, n_d, ncoa, &
&                ncob, ncoc, ncod, offset_a, offset_b, offset_c, offset_d, &
&                primitives, hfx_opts, max_contraction, tmp_max, neris, &
&                ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, &
&                rpq2, AB, CD, R1, R2)

    use alloc, only: de_alloc, re_alloc
    use sys, only: die
    use precision, only: dp, int_8
    use nao2gto_data, only: nco              ! Number of Cartesian Gaussian 
                                             !   functions per angular momentum 
                                             !   shell
    use nao2gto_libint
    use nao2gto_types, only: hfx_options_type
    use nao2gto_wrappers, only: get_eris

    implicit none

!   Arguments
    type(Libint_t), intent(inout)         :: lib
    real(dp), intent(in)                  :: A(3)
                                             ! Center of the first primitive 
                                             !    Cartesian Gaussian function
    real(dp), intent(in)                  :: B(3)
                                             ! Center of the second primitive 
                                             !    Cartesian Gaussian function
    real(dp), intent(in)                  :: C(3)
                                             ! Center of the third primitive 
                                             !    Cartesian Gaussian function
    real(dp), intent(in)                  :: D(3)
                                             ! Center of the fourth primitive 
                                             !    Cartesian Gaussian function
    integer, intent(in)                   :: n_a 
                                             ! Angular momentum of the first
                                             !   primitive Cartesian Gaussian
                                             !   function
    integer, intent(in)                   :: n_b 
                                             ! Angular momentum of the second
                                             !   primitive Cartesian Gaussian
                                             !   function
    integer, intent(in)                   :: n_c 
                                             ! Angular momentum of the third
                                             !   primitive Cartesian Gaussian
                                             !   function
    integer, intent(in)                   :: n_d 
                                             ! Angular momentum of the fourth
                                             !   primitive Cartesian Gaussian
                                             !   function
    integer,  intent(in)                  :: ncoa 
                                             ! Total number of Cartesian 
                                             !   Gaussians
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
    integer,  intent(in) :: ncob             ! Same as ncoa for the 2nd orbital
    integer,  intent(in) :: ncoc             ! Same as ncoa for the 3rd orbital
    integer,  intent(in) :: ncod             ! Same as ncoa for the 4th orbital

    integer, intent(in)                   :: offset_a, &
&                                            offset_b, offset_c, offset_d
    real(dp), &
&     dimension(ncoa, ncob, ncoc, ncod), &
&     intent(out)                         :: primitives
    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in) :: max_contraction
    real(dp), intent(out) :: tmp_max
    integer(int_8), intent(inout) :: neris
    real(dp), intent(in)                  :: AB(3) 
                           ! \vec{A} - \vec{B}
    real(dp), intent(in)                  :: CD(3) 
                           ! \vec{C} - \vec{D}
    real(dp), intent(in)                 :: Rho
                           !   [(\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)]
                           !    that appears in Eq. (6.4) of Ref. Fermann:2020
   real(dp), intent(in)                   :: RhoInv
                           !   RhoInv = (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q)
    real(dp), intent(in)                  :: ZetapEtaInv
                           !   ZetapEtaInv = 1.0/(Zeta + Eta),
                           !       where Zeta = Zeta_A + Zeta_B,
                           !       and Eta = Zeta_C + Zeta_D,
                           !   ZetapEtaInv = 1.0/(gamma_p + \gamma_q)
    real(dp), intent(in)                  :: S1234
                           !   Product of K_{1}*K_{2} in Eq. (6.4) of
                           !     Ref. \cite{Fermann:2020}
    real(dp), intent(in)                 :: P(3)
                           !   \vec{P}=\frac{\zeta_{a}\vec{A}+\zeta_{b}\vec{B}}
                           !                {\zeta}
    real(dp), intent(in)                 :: Q(3)
                           !   \vec{P}=\frac{\zeta_{c}\vec{C}+\zeta_{d}\vec{D}}
                           !                {\eta}
    real(dp), intent(in)                 :: W(3)
                           !   \vec{W}=\frac{\zeta\vec{P}+\eta\vec{Q}}
                           !                {\zeta + \eta}
    real(dp), intent(in)                 :: rpq2
                           !   The variable PQ2 corresponds with
                           !   $\overline{\vec{P}\vec{Q}}^{2}$
                           !   in Eq. (6.4) of Ref. Fermann:2020



    real(dp), intent(in)                  :: ZetaInv, EtaInv, &
&                                             R1, R2

!   Local variables
    character(len=120)      :: msg
    integer                         :: mysize 
                                             ! Number of total primitive ERI 
                                             !   integrals per quartet
    integer                         :: m_max
                                             ! m_max is the maximum index to 
                                             !   compute the values of the 
                                             !   auxiliary primitive
                                             !   integrals (00|00)^{(m)}. 
                                             !   When computing ERIS, 
                                             !   it runs from 
                                             !   0 \le m \le \lambda_{a} + 
                                             !               \lambda_{b} + 
                                             !               \lambda_{c} + 
                                             !               \lambda_{d}   
                                             !   where the \lambda are the 
                                             !   angular momentums of the four 
                                             !   primitive Cartesian Gaussian 
                                             !   functions
    integer                         :: perm_case
                                             ! Permutational case to compute
                                             !   the ERIs  
    integer                         :: a_mysize(1), i, j, k, l, &
&                                      p1, p2, p3, &
&                                      p4, q1, q2, q3, q4
    logical                         :: do_it
    real(dp), dimension(:), pointer :: p_work => null()
    type(prim_data), target         :: prim
                                             ! Primitive quartet data
                                             !   The structure is defined in
                                             !   pages 7-8 of LIBINT Manual

! ------------------------------------------------------------------------------

!   m_max is the maximum index to compute the values of the auxiliary primitive
!   integrals (00|00)^{(m)}. 
!   When computing eris, it runs from 
!   0 \le m \le \lambda_{a} + \lambda_{b} + \lambda{c} + \lambda{d},
!   where the \lambda are the angular momentums of the four primitive
!   Cartesian Gaussian functions
    m_max = n_a + n_b + n_c + n_d

!   Number of total primitive ERI integrals per quartet
!   See page 4 of the LIBINT Manual
!   For instance, a (ps|sd) class or quartet, consists of 
!   3 x 1 x 1 x 6 = 18 integrals
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize(:) = mysize

    do_it = .true.

    if( m_max /= 0 ) then     ! At least one of the Cartesian Gaussian function
                              !   is not of s-type
!     Note that currently LIBINT has a very important restriction on the 
!     angular momentum ordering of the functions in shell quartets that it 
!     can handle. 
!     LIBINT can evaluate a shell quartet (ab|cd) if 
!     \lambda(a) \ge \lambda(b),
!     \lambda(c) \ge \lambda(d),
!     and \lambda(c) + \lambda(d) \ge \lambda(a) + \lambda(b). 
!     If one needs to compute a quartet that doesn’t conform the rule, e.g. 
!     of type (pf|sd), permutational symmetry of integrals can be utilized 
!     to compute such quartet
!     (pq|rs)=(pq|sr)=(qp|rs)=(qp|sr)=(rs|pq)=(rs|qp)=(sr|pq)=(sr|qp)
!     In the case of (pf|sd) shell quartet, one computes quartet (ds|fp) 
!     instead, and then permutes function indices back to obtain the 
!     desired (pf|sd).

      perm_case = 1
      if( n_a < n_b ) then
        perm_case = perm_case + 1
      end if
      if( n_c < n_d ) then
        perm_case = perm_case + 2
      end if
      if( n_a+n_b > n_c+n_d ) then
        perm_case = perm_case + 4
      end if

      select case(perm_case)

        case(1)
          call build_quartet_data(A, B, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=AB ! A-B
          lib%CD=CD ! C-D
!         Here we compute the electron repulsion integrals 
!         written in Eq. (6.2) of Ref. Fermann:2020,
!         \f$ I =  (ij \vert kl) = \int 
!         \phi_{i}(\vec{r}_{1}) 
!         \phi_{j}(\vec{r}_{1}) 
!         \frac{1}{r_{12}} 
!         \phi_{k}(\vec{r}_{2}) 
!         \phi_{l}(\vec{r}_{2}) 
!         d \vec{r}_{1} d \vec{r}_{2} \f$
!         by a recursive application of a simple relation.
          call get_eris( n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize )

          neris = neris + 1
          do i = 1, mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do i = 1,nco(n_a)
            p1 = (i-1) *nco(n_b)
            q1 = offset_a + i
            do j = 1,nco(n_b)
              p2 = (p1 + j-1)*nco(n_c)
              q2 = offset_b + j
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c + k
                do l = 1,nco(n_d)
                  q4 = offset_d + l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(2)
          call build_quartet_data(B, A, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-AB ! B-A
          lib%CD=CD ! C-D
          call get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b+j
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_c)
              q1 = offset_a+i
              do k = 1,nco(n_c)
                p3 = (p2 + k-1)*nco(n_d)
                q3 = offset_c+k
                do l = 1,nco(n_d)
                  q4 = offset_d+l
                  p4 = p3 + l
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(3)
          call build_quartet_data(A, B, D, C, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=AB ! A-B
          lib%CD=-CD ! D-C
          call get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do i = 1,nco(n_a)
            p1 = (i-1)*nco(n_b)
            q1 = offset_a + i
            do j = 1,nco(n_b)
              q2 = offset_b + j
              p2 = (p1 + j-1)*nco(n_d)
              do l = 1,nco(n_d)
                q4 = offset_d + l
                p3 = (p2 + l-1) * nco(n_c)
                do k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3+k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(4)
          call build_quartet_data(B, A, D, C, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-AB ! B-A
          lib%CD=-CD ! D-C
          call get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do j = 1,nco(n_b)
            p1 = (j-1)*nco(n_a)
            q2 = offset_b + j
            do i = 1,nco(n_a)
              p2 = (p1 + i-1)*nco(n_d)
              q1 = offset_a + i
              do l = 1,nco(n_d)
                p3 = (p2 + l-1)*nco(n_c)
                q4 = offset_d + l
                do k = 1,nco(n_c)
                  q3 = offset_c + k
                  p4 = p3 + k
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(5)
          call build_quartet_data(C, D, A, B, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=CD ! C-D
          lib%CD=AB ! A-B
          call get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do k = 1,nco(n_c)
            q3 = offset_c + k
            p1 = (k-1)*nco(n_d)
            do l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_a)
              do i = 1,nco(n_a)
                q1 = offset_a + i
                p3 = (p2 + i-1)*nco(n_b)
                do j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(6)
          call build_quartet_data(C, D, B, A, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=CD ! C-D
          lib%CD=-AB ! B-A
          call get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do k = 1,nco(n_c)
            p1 = (k-1)*nco(n_d)
            q3 = offset_c + k
            do l = 1,nco(n_d)
              q4 = offset_d + l
              p2 = (p1 + l-1)*nco(n_b)
              do j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1)*nco(n_a)
                do i = 1,nco(n_a)
                  p4 = p3+i
                  q1 = offset_a + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(7)
          call build_quartet_data(D, C, A, B, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-CD ! D-C
          lib%CD=AB ! A-B
          call get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do l = 1,nco(n_d)
            p1 = (l-1)*nco(n_c)
            q4 = offset_d + l
            do k = 1,nco(n_c)
              p2 = (p1 + k-1) * nco(n_a)
              q3 = offset_c + k
              do i = 1,nco(n_a)
                p3 = (p2 + i-1) *nco(n_b)
                q1 = offset_a + i
                do j = 1,nco(n_b)
                  q2 = offset_b + j
                  p4 = p3+j
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

        case(8)
          call build_quartet_data(D, C, B, A, m_max, &
                            hfx_opts, prim, do_it, &
                            EtaInv, ZetaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            Q, P, W, rpq2)
          if( .not. do_it ) return
          lib%AB=-CD ! D-C
          lib%CD=-AB ! B-A
          call get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)

          neris = neris + 1
          do i=1,mysize
            tmp_max = max(tmp_max, abs(p_work(i)))
          end do
          tmp_max = tmp_max*max_contraction
          if ( tmp_max < hfx_opts%eps_schwarz ) return

          do l = 1,nco(n_d)
            q4 = offset_d + l
            p1 = (l-1)*nco(n_c)
            do k = 1,nco(n_c)
              q3 = offset_c + k
              p2 = (p1 + k-1) * nco(n_b)
              do j = 1,nco(n_b)
                q2 = offset_b + j
                p3 = (p2 + j-1) * nco(n_a)
                do i = 1,nco(n_a)
                  q1 = offset_a + i
                  p4 = p3 + i
                  primitives(q1,&
                             q2,&
                             q3,&
                             q4)=p_work(p4)
                end do
              end do
            end do
          end do

      case default
        write(msg, fmt='(A,": ",I2)') "Invalid HFX permutation case", &
          perm_case
        call die(msg)

      end select

    else   ! if m_max = 0 (all the Cartesian Gaussian functions are s-type

      call build_quartet_data(A, B, C, D, m_max, &
                            hfx_opts, prim, do_it, &
                            ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, S1234, &
                            P, Q, W, rpq2)
      neris = neris + 1

      if( .not. do_it ) return

      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = prim%F(1)
      tmp_max = max_contraction*abs(prim%F(1))

    end if

  end subroutine calc_primitive_eri

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

! *****************************************************************************
  SUBROUTINE calc_primitive_eri2(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,&
                          ncoa,ncob,ncoc,ncod,&
                          offset_a,offset_b,offset_c,offset_d,&
                          primitives, hfx_opts,    &
                          Zeta,ZetaInv,Eta,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                          P,Q,W,rpq2,AB,CD)

    use sys, only: die
    use precision, only: dp
    use nao2gto_data, only: nco
    use nao2gto_libint
    use nao2gto_types
    use nao2gto_wrappers, only: get_eris
!   For debugging
    use parallel,    only : Node
#ifdef MPI
    use mpi_siesta
#endif
!   End debugging

    implicit none
    
    TYPE(Libint_t)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d, ncoa, &
                                                ncob, ncoc, ncod, offset_a, &
                                                offset_b, offset_c, offset_d
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitives
    TYPE(hfx_options_type)                :: hfx_opts

    REAL(dp), INTENT(IN)                     :: Zeta, ZetaInv, Eta, EtaInv, &
                                                ZetapEtaInv, Rho, RhoInv, &
                                                S1234, P(3), Q(3), W(3), &
                                                rpq2, AB(3), CD(3)

    character(len=120)      :: msg
    INTEGER                                  :: a_mysize(1), i, j, k, l, &
                                                m_max, mysize, p1, p2, p3, &
                                                p4, perm_case, q1, q2, q3, q4
    LOGICAL                                  :: do_it
    REAL(dp), DIMENSION(:), POINTER          :: p_work => null()
    TYPE(prim_data), TARGET                  :: prim

!   For debugging
    integer     :: ia, ib, ic, id
#ifdef MPI
    integer     :: MPIerror
#endif
!   End debugging

    m_max = n_a+n_b+n_c+n_d
    mysize = nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)
    a_mysize = mysize
    do_it = .TRUE.
 
!!   For debugging
!    if( Node .eq. 0 ) then
!    write(6,'(a)')                                                       &
! &    'calc_primitive_eri2: Data of the first primitive function: '   
!    write(6,'(a,i5,4f12.5,3i7)')                                         &
! &    'calc_primitive_eri2: Node, A,    ZetaA, la, nco, offset ',        &
! &                          Node, A(:), Zeta_A, n_a, nco(n_a), offset_a
!    write(6,'(a)')                                                       &
! &    'calc_primitive_eri2: Data of the second primitive function: '  
!    write(6,'(a,i5,4f12.5,3i7)')                                         &
! &    'calc_primitive_eri2: Node, B,    ZetaB, lb, nco, offset ',        &
! &                          Node, B(:), Zeta_B, n_b, nco(n_b), offset_b
!    write(6,'(a)')                                                       &
! &    'calc_primitive_eri2: Data of the third primitive function: '   
!    write(6,'(a,i5,4f12.5,3i7)')                                         &
! &    'calc_primitive_eri2: Node, C,    ZetaC, lc, nco, offset ',        &
! &                          Node, C(:), Zeta_C, n_c, nco(n_c), offset_c
!    write(6,'(a)')                                                       &
! &    'calc_primitive_eri2: Data of the fourth primitive function: '  
!    write(6,'(a,i5,4f12.5,3i7)')                                         &
! &    'calc_primitive_eri2: Node, D,    ZetaD, ld, nco, offset ',        &
! &                          Node, D(:), Zeta_D, n_d, nco(n_d), offset_d
!    write(6,'(a,3i7)')                                                   &
! &    'calc_primitive_eri2: Node, m_max, mysize = ',                     &
! &                          Node, m_max, mysize
!    write(6,'(a)') &
! &    'calc_primitive_eri2: ---------------------------------------'
!    endif
!!   End debugging


    IF ( m_max == 0 ) THEN
      CALL build_quartet_data(A,B,C,D, m_max,&
                            hfx_opts, prim, do_it, &
                            ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
                            P,Q,W,rpq2)
      IF( .NOT. do_it ) RETURN
      primitives(offset_a+1,offset_b+1,offset_c+1,offset_d+1) = prim%F(1)
      return
    END IF

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
        CALL build_quartet_data(A,B,C,D, m_max,&
                          hfx_opts, prim, do_it,  &
                          ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(B,A,C,D, m_max,&
                          hfx_opts, prim, do_it, &
                          ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(A,B,D,C, m_max,&
                          hfx_opts, prim, do_it, &
                          ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(B,A,D,C, m_max,&
                          hfx_opts, prim, do_it, &
                          ZetaInv,EtaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(C,D,A,B, m_max,&
                          hfx_opts, prim, do_it ,&
                          EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(C,D,B,A, m_max,&
                          hfx_opts, prim, do_it ,&
                          EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(D,C,A,B, m_max,&
                          hfx_opts, prim, do_it ,&
                          EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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
        CALL build_quartet_data(D,C,B,A, m_max,&
                          hfx_opts, prim, do_it ,&
                          EtaInv,ZetaInv,ZetapEtaInv,Rho,RhoInv,S1234,&
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

      case default
        write(msg, fmt='(A,": ",I2)') "Invalid HFX permutation case", &
          perm_case
        call die(msg)

    END SELECT

!!   For debugging
!    if( Node .eq. 0 ) then
!    write(6,'(a)')                                                       &
! &    'calc_primitive_eri2: Primitive integrals: '                      
!    do ia = 1, ncoa
!      do ib = 1, ncob
!        do ic = 1, ncoc
!          do id = 1, ncod
!            if( primitives(ia,ib,ic,id) .gt. 1.d-12) then 
!            write(6,'(a,5i5,f20.12)')                                     &
! &            'calc_primitive_eri2: Node, ia, ib, ic, id, primitives = ', &
! &            Node, ia, ib, ic, id, primitives(ia,ib,ic,id)
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!    endif
!
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!    call die()
!   End debugging


  END SUBROUTINE calc_primitive_eri2

! *****************************************************************************
!> \brief Evaluate electron repulsion integrals for a primitive quartet
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      08.2007 refactured permutation part [Manuel Guidon]
!> \author Manuel Guidon

! *****************************************************************************
  SUBROUTINE calc_primitive_screen(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                   n_a,n_b,n_c,n_d,&
                                   max_val, hfx_opts)

    use sys, only: die
    use precision, only: dp
    use nao2gto_data, only: nco
    use nao2gto_libint
    use nao2gto_types
    use nao2gto_wrappers, only: get_eris

    implicit none

    TYPE(Libint_t)                            :: lib
    REAL(dp), INTENT(IN)                     :: A(3), B(3), C(3), D(3), &
                                                Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: n_a, n_b, n_c, n_d
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_options_type)                :: hfx_opts

    character(len=120)      :: msg
    INTEGER                                  :: a_mysize(1), i, m_max, &
                                                mysize, perm_case
    LOGICAL                                  :: do_it
    REAL(dp), DIMENSION(:), POINTER          :: p_work => null()
    TYPE(prim_data), TARGET                  :: prim

!permutation of configuration

    m_max = n_a+n_b+n_c+n_d
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
          CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                            hfx_opts, prim, do_it, n_a+n_b, n_c+n_d)
          lib%AB=A-B
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(2)
          CALL build_quartet_data_screen(B,A,C,D,Zeta_B, Zeta_A, Zeta_C, Zeta_D, m_max,&
                            hfx_opts, prim, do_it, n_b+n_a, n_c+n_d)
          lib%AB=B-A
          lib%CD=C-D
          CALL get_eris(n_d, n_c, n_a, n_b, lib, prim, p_work, a_mysize)

          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(3)
          CALL build_quartet_data_screen(A,B,D,C,Zeta_A, Zeta_B, Zeta_D, Zeta_C, m_max,&
                            hfx_opts, prim, do_it, n_a+n_b, n_d+n_c)
          lib%AB=A-B
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_b, n_a, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(4)
          CALL build_quartet_data_screen(B,A,D,C,Zeta_B, Zeta_A, Zeta_D, Zeta_C, m_max,&
                            hfx_opts, prim, do_it, n_b+n_a, n_d+n_c)
          lib%AB=B-A
          lib%CD=D-C
          CALL get_eris(n_c, n_d, n_a, n_b, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(5)
          CALL build_quartet_data_screen(C,D,A,B,Zeta_C, Zeta_D, Zeta_A, Zeta_B, m_max,&
                            hfx_opts, prim, do_it, n_c+n_d, n_a+n_b)
          lib%AB=C-D
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(6)
          CALL build_quartet_data_screen(C,D,B,A,Zeta_C, Zeta_D, Zeta_B, Zeta_A, m_max,&
                            hfx_opts, prim, do_it, n_c+n_d, n_b+n_a)
          lib%AB=C-D
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_d, n_c, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(7)
          CALL build_quartet_data_screen(D,C,A,B,Zeta_D, Zeta_C, Zeta_A, Zeta_B, m_max,&
                            hfx_opts, prim, do_it, n_d+n_c, n_a+n_b)
          lib%AB=D-C
          lib%CD=A-B
          CALL get_eris(n_b, n_a, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO
        CASE(8)
          CALL build_quartet_data_screen(D,C,B,A,Zeta_D, Zeta_C, Zeta_B, Zeta_A, m_max,&
                            hfx_opts, prim, do_it, n_d+n_c, n_b+n_a)
          lib%AB=D-C
          lib%CD=B-A
          CALL get_eris(n_a, n_b, n_c, n_d, lib, prim, p_work, a_mysize)
          DO i=1,mysize
            max_val = MAX(max_val, ABS(p_work(i)))
          END DO

      case default
        write(msg, fmt='(A,": ",I2)') "Invalid HFX permutation case", &
          perm_case
        call die(msg)

      END SELECT

    ELSE
      CALL build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max,&
                              hfx_opts, prim, do_it, 0, 0)
      max_val = ABS(prim%F(1))
    END IF
  END SUBROUTINE calc_primitive_screen

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! *****************************************************************************
!> \brief Fills the data structure used in Libderiv.
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in,out] lib: Libint data structure
!! \param[in] A: Center of the first primitive Cartesian Gaussian function
!! \param[in] B: Center of the second primitive Cartesian Gaussian function
!! \param[in] C: Center of the third primitive Cartesian Gaussian function
!! \param[in] D: Center of the fourth primitive Cartesian Gaussian function
!! \param[in] ZetaA: Exponent of the first primitive Cartesian Gaussian function
!! \param[in] ZetaB: Exponent of the second primitive Cartesian Gaussian functio
!! \param[in] ZetaC: Exponent of the third primitive Cartesian Gaussian function
!! \param[in] ZetaD: Exponent of the fourth primitive Cartesian Gaussian functio
!! \param[in] m_max: Maximum \f$ m \f$ value of the incomplete Gamma function
!!            \f$F_{m}(T)\f$, defined in Eq. (6.4) of Ref. \cite Fermann:2020 as
!!            \f$ F_{m} (T) = \int_{0}^{1} dt \:\: t^{2m} \exp(-Tt^{2}) \f$.
!!            In page 7 of the LIBINT Manual, it is specified how:
!!             - When computing electron repulsion integrals, 
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) \f$,
!!            where \f$ l \f$ stands for the angular momentum of the 
!!            corresponding Cartesian Gaussian function.
!!             - When computing first derivatives of the 
!!            electron repulsion integrals, 
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) + 1 \f$.
!!             - When computing second derivatives of the electron repulsion 
!!            integrals, and integrals for linear R12 methods, 
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) + 2 \f$.
!! \param[in,out] prim: Internal coefficients for Libint and Libderiv
!! \param[in,out] do_it: ...
!! \param[in] ZetaInv: \f$ 1.0/ ( \zeta_a + \zeta_b ) \f$ 
!!            where \f$ \zeta_{a} \f$ and \f$ \zeta_{b} \f$ are the exponents
!!            of the first and second primitive Cartesian Gaussian functions
!! \param[in] EtaInv: \f$ 1.0/ ( \zeta_c + \zeta_d ) \f$ 
!!            where \f$ \zeta_{c} \f$ and \f$ \zeta_{d} \f$ are the exponents
!!            of the third and fourth primitive Cartesian Gaussian functions
!! \param[in] ZetapEtaInv: \f$ 1.0/(\zeta + \eta) \f$ in the notation of 
!!                   LIBINT Manual, where \f$ \zeta = \zeta_{a} + \zeta_{b} \f$,
!!                   and \f$ \eta = \zeta_{c} + \zeta_{b} \f$.
!!                   Or \f$ 1.0/(\gamma_p + \gamma_q) \f$ in the notation of 
!!                   Eq. (6.4) of
!!                   Ref. \cite Fermann:2020
!! \param[in] Rho: \f$ (\zeta \times \eta) / (\zeta+\eta) \f$ in the notation of
!!                   LIBINT Manual. 
!!                   Or \f$ (\gamma_p \times \gamma_q) / (\gamma_p+\gamma_q) \f$
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] RhoInv: \f$ (\zeta+\eta)/(\zeta \times \eta) \f$ in the notation 
!!                   of LIBINT Manual.
!!                   Or \f$ (\gamma_p+\gamma_q)/(\gamma_p \times \gamma_q) \f$ 
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] S1234: Product of \f$ K_{1}\times K_{2} \f$ in Eq. (6.4) of 
!!                   Ref. \cite Fermann:2020, and Eq. (15) and Eq. (16) of 
!!                   LIBINT Manual (here called \f$S_{12}\f$ and \f$S_{34}\f$).
!!                   \f$ \left[ \exp \left( - \frac{\zeta_{a} \zeta_{b}}{\zeta} 
!!                        \vert \vec{AB} \vert^{2} \right) \right] 
!!                   \times \left[ 
!!                   \exp \left( - \frac{\zeta_{c} \zeta_{d}}{\eta}  
!!                        \vert \vec{CD} \vert^{2} \right) \right]    \f$
!! \param[in] P: \f$ \vec{P} = \frac{\zeta_{a} \vec{A} + \zeta_{b} \vec{B}}
!!                             {\zeta} \f$
!! \param[in] Q: \f$ \vec{Q} = \frac{\zeta_{c} \vec{C} + \zeta_{d} \vec{D}}
!!                             {\eta} \f$
!! \param[in] W: \f$ \vec{W} = \frac{\zeta \vec{P} + \eta \vec{Q}}
!!                             {\zeta + \eta} \f$
!! \param[in] PQ2: The variable PQ2 corresponds with
!!                 \f$ \overline{\vec{PQ}}^{2} \f$
!!                 in Eq. (6.4) of Ref. \cite Fermann:2020
! *****************************************************************************
  subroutine build_deriv_data(A, B, C, D, Zeta_A, Zeta_B, Zeta_C, Zeta_D, &
&                m_max, hfx_opts, prim, do_it, ZetaInv, EtaInv, &
&                ZetapEtaInv, Rho, RhoInv, S1234, P, Q, W, PQ2)

    use sys, only: die
    use units,          only: pi
    use precision,  only: dp
    use nao2gto_libint
    use nao2gto_utils,  only: fgamma => fgamma_0
                           ! Calculation of the incomplete Gamma function 
                           !    F(T) for multicenter integrals over 
                           !    Gaussian functions. 
                           !    f returns a vector with all \f$F_m(T)\f$ 
                           !    values for 0 <= m <= m_max.
    use nao2gto_types,  only: hfx_options_type
    use nao2gto_types,  only: do_hfx_potential_coulomb, &
   &                          do_hfx_potential_short,   &
   &                          do_hfx_potential_truncated

    implicit none

    ! Arguments
    real(dp), intent(in)               :: A(3)
                           ! Center of the first primitive Gaussian
    real(dp), intent(in)               :: B(3)
                           ! Center of the second primitive Gaussian
    real(dp), intent(in)               :: C(3)
                           ! Center of the third primitive Gaussian
    real(dp), intent(in)               :: D(3)
                           ! Center of the fourth primitive Gaussian
    real(dp), intent(in)               :: Zeta_A
                           ! Exponent of the first  primitive Gaussian
    real(dp), intent(in)               :: Zeta_B
                           ! Exponent of the second primitive Gaussian
    real(dp), intent(in)               :: Zeta_C
                           ! Exponent of the third  primitive Gaussian
    real(dp), intent(in)               :: Zeta_D
                           ! Exponent of the fourth primitive Gaussian
    integer,  intent(in)               :: m_max
                           ! Maximum \f$ m \f$ value of the 
                           !    incomplete Gamma function, defined in 
                           !    Eq. (6.4) of Ref. \cite{Fermann:2020}.
                           !    In page 34 of this Reference, 
                           !    m is said to run from
                           !    0 \le m \le (la + lb + lc + ld) + 1
                           !    where (la + lb + lc + ld) are the angular
                           !    momenta of the four primitive 
                           !    Cartesian Gaussian functions 
                           !    (see page 28 of that Reference)
                           !    The +1 comes from the fact that we are
                           !    computing the derivatives.
                           !    See page 7 of LIBINT Manual
    type(hfx_options_type), intent(in) :: hfx_opts
    type(prim_data), intent(out)       :: prim
                           !    Internal coefficients for Libint and Libderiv
    logical, intent(inout)             :: do_it
    real(dp), intent(in)               :: Rho
                           !   [(\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)]
                           !    that appears in Eq. (6.4) of Ref. Fermann:2020
    real(dp), intent(in)               :: PQ2
                           !   The variable PQ2 corresponds with
                           !   $\overline{\vec{P}\vec{Q}}^{2}$
                           !   in Eq. (6.4) of Ref. Fermann:2020
    real(dp), intent(in)               :: RhoInv 
                           !   RhoInv = (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q)
    real(dp), intent(in)               :: ZetapEtaInv
                           !   ZetapEtaInv = 1.0/(Zeta + Eta),
                           !       where Zeta = Zeta_A + Zeta_B,
                           !       and Eta = Zeta_C + Zeta_D,
                           !   ZetapEtaInv = 1.0/(gamma_p + \gamma_q)
    real(dp), intent(in)               :: S1234
                           !   Product of K_{1}*K_{2} in Eq. (6.4) of 
                           !     Ref. \cite{Fermann:2020}
    real(dp), intent(in)               :: ZetaInv, EtaInv,  &
&                                         P(3), Q(3), &
&                                         W(3)

!   Local variables
    character(len=120)      :: msg ! Error message
    integer                 :: i
    real(dp)                :: T  ! T = \overline{\vec{P}\vec{Q}}^{2} *
                                  !  (\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)
                                  !  This is the argument of the function 
                                  !  F_{m}(T), defined in page 26
                                  !  of Ref. Fermann:2020
    real(dp)                :: factor
                                  ! Prefactor of the incomplete Gamma factor
    real(dp)                :: tmp   
                                  ! Temporary variable to compute the prefactor
                                  !   of the incomplete Gamma function
    real(dp)                :: omega2, omega_corr, omega_corr2
    real(dp), dimension(17) :: Fm
                                  ! Incomplete Gamma function
                                  ! F_m (T) = 
                                  !   \int_{0}^{1} dt \: t^{2m} \exp(-Tt^{2}} 

! -------------------------------------------------------------------------
!   Define the variable of the function F_m, defined in page 26 of 
!   Ref. Fermann:2020
!   T = \overline{\vec{P}\vec{Q}}^{2} *
!       (\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)
!   or in the notation of the LIBINT Manual
!   T = \overline{\vec{P}\vec{Q}}^{2} * (\zeta * \eta)/(\zeta + \eta)
    T = Rho*PQ2

    do_it = .true.

    select case (hfx_opts%potential_type)

      case(do_hfx_potential_coulomb)
        call fgamma(m_max, T, prim%F(1))
        factor = 2.0_dp*Pi*RhoInv

      case(do_hfx_potential_short)
        call fgamma(m_max, T, prim%F)
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = dsqrt(omega_corr2)
        T = T*omega_corr2
        call fgamma(m_max,T,Fm)
        tmp = - omega_corr
        do i=1,m_max+1
          prim%F(i) = prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        end do
        factor = 2.0_dp*Pi*RhoInv

      case default
        write(msg, fmt='(A,": ",I2)') "Invalid Hartree-Fock potential type", &
          hfx_opts%potential_type
        call die(msg)

    end select

!   Factor is the term that appears right after the equal sign in Eq. (6.4)
!   of Ref. \cite Fermann:2020
!   $\frac {2 \pi^{5/2} K_{1} K_{2} } 
!          {\gamma_p \gamma_q (\gamma_p + \gamma_q)^{1/2}} $
!   Up to this point, factor = 2 \pi * RhoInv =
!                            = 2 * \pi * (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q)
!   Now, we multiply it by tmp = (Pi*ZetapEtaInv)**3
!                            = \pi**{3} * 1.0/(gamma_p + \gamma_q)**{3/2}
!   The result is 2 * \pi^{5/2} * 
!                 1.0/((\gamma_p*\gamma_q)*(gamma_p + \gamma_q)**{1/2})
!   Finally, we multiply this by K_1* K_{2},
!   where K_{1} = exp(-\alpha_{1} * \alpha_{2} * {\overline{\vec{PQ}}}^{2} 
!                    / gamma_p}
!   and K_{2} = exp(-\alpha_{3} * \alpha_{4} * {\overline{\vec{PQ}}}^{2} 
!                    / gamma_q}
!   This product of the two K is stored in the variable S1234
    tmp = (Pi*ZetapEtaInv)**3
    factor = factor * S1234 * dsqrt(tmp)

!   Populate the entries of each component of prim_data, defined
!   as in page 7 of the LIBINT Programmer's Manual
!   prim%F(i) is an auxiliary integral over primitive s-functions,
!   referred to as (0 0 | 0 0)^{(m)} in LIBINT Programmer's Manual
    do i = 1, m_max + 1
      prim%F(i) = prim%F(i) * factor
    end do

    prim%U(:,1)    = P-A
    prim%U(:,2)    = P-B
    prim%U(:,3)    = Q-C
    prim%U(:,4)    = Q-D
    prim%U(:,5)    = W-P
    prim%U(:,6)    = W-Q
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

  end subroutine build_deriv_data

! *****************************************************************************
!> \brief Fills the data structure used in Libint (see page 7 of LIBINT Manual).
!! In particular, we compute the \f$(00|00)^{m} \f$ integrals defined
!! in Eq. (14) of LIBINT Manual, and stored here in 
!! prim%F(i)
!!
!! \author Manuel Guidon
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2007 Created [Manuel Guidon]
!!      - 01.2018 Changed interface for SIESTA [Yann Pouillon]
!!
!! \param[in] A: Center of the first primitive Cartesian Gaussian function
!! \param[in] B: Center of the second primitive Cartesian Gaussian function
!! \param[in] C: Center of the third primitive Cartesian Gaussian function
!! \param[in] D: Center of the fourth primitive Cartesian Gaussian function
!! \param[in] m_max: Maximum \f$ m \f$ value of the incomplete Gamma function
!!            \f$F_{m}(T)\f$, defined in Eq. (6.4) of Ref. \cite Fermann:2020 as
!!            \f$ F_{m} (T) = \int_{0}^{1} dt \:\: t^{2m} \exp(-Tt^{2}) \f$.
!!            In page 7 of the LIBINT Manual, it is specified how:
!!             - When computing electron repulsion integrals,
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) \f$,
!!            where \f$ l \f$ stands for the angular momentum of the
!!            corresponding Cartesian Gaussian function.
!!             - When computing first derivatives of the
!!            electron repulsion integrals,
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) + 1 \f$.
!!             - When computing second derivatives of the electron repulsion
!!            integrals, and integrals for linear R12 methods,
!!            \f$ 0 \le m \le  l(a)+l(b)+l(c)+l(d) + 2 \f$.
!! \param[in] hfx_opts: data structure containing Hartree-Fock exchange
!!                           parameters
!! \param[in,out] prim: Internal coefficients for Libint and Libderiv.
!!                      Contains all the information of a primitive quartet
!! \param[in,out] do_it: ...
!! \param[in] ZetaInv: \f$ 1.0/ ( \zeta_a + \zeta_b ) \f$ 
!!            where \f$ \zeta_{a} \f$ and \f$ \zeta_{b} \f$ are the exponents
!!            of the first and second primitive Cartesian Gaussian functions
!! \param[in] EtaInv: \f$ 1.0/ ( \zeta_c + \zeta_d ) \f$ 
!!            where \f$ \zeta_{c} \f$ and \f$ \zeta_{d} \f$ are the exponents
!!            of the third and fourth primitive Cartesian Gaussian functions
!! \param[in] ZetapEtaInv: \f$ 1.0/(\zeta + \eta) \f$ in the notation of 
!!                   LIBINT Manual, where \f$ \zeta = \zeta_{a} + \zeta_{b} \f$,
!!                   and \f$ \eta = \zeta_{c} + \zeta_{b} \f$.
!!                   Or \f$ 1.0/(\gamma_p + \gamma_q) \f$ in the notation of 
!!                   Eq. (6.4) of
!!                   Ref. \cite Fermann:2020
!! \param[in] Rho: \f$ (\zeta \times \eta) / (\zeta+\eta) \f$ in the notation of
!!                   LIBINT Manual. 
!!                   Or \f$ (\gamma_p \times \gamma_q) / (\gamma_p+\gamma_q) \f$
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] RhoInv: \f$ (\zeta+\eta)/(\zeta \times \eta) \f$ in the notation 
!!                   of LIBINT Manual.
!!                   Or \f$ (\gamma_p+\gamma_q)/(\gamma_p \times \gamma_q) \f$ 
!!                   in Eq. (6.4) of Ref. \cite Fermann:2020
!! \param[in] S1234: Product of \f$ K_{1}\times K_{2} \f$ in Eq. (6.4) of
!!                   Ref. \cite Fermann:2020, also specified in 
!!                   Eq. (15) and Eq. (16) of 
!!                   LIBINT Manual (here called \f$S_{12}\f$ and \f$S_{34}\f$).
!!                   \f$ \left[ \exp \left( - \frac{\zeta_{a} \zeta_{b}}{\zeta}
!!                        \vert \vec{AB} \vert^{2} \right) \right]
!!                   \times \left[
!!                   \exp \left( - \frac{\zeta_{c} \zeta_{d}}{\eta}
!!                        \vert \vec{CD} \vert^{2} \right) \right]    \f$
!! \param[in] P: \f$ \vec{P} = \frac{\zeta_{a} \vec{A} + \zeta_{b} \vec{B}}
!!                             {\zeta} \f$. See Eq. (10) of LIBINT Manual.
!! \param[in] Q: \f$ \vec{Q} = \frac{\zeta_{c} \vec{C} + \zeta_{d} \vec{D}}
!!                             {\eta} \f$. See Eq. (11) of LIBINT Manual.
!! \param[in] W: \f$ \vec{W} = \frac{\zeta \vec{P} + \eta \vec{Q}}
!!                             {\zeta + \eta} \f$. See Eq. (12) of 
!!                             LIBINT Manual.
!! \param[in] PQ2: The variable PQ2 corresponds with
!!                 \f$ \overline{\vec{PQ}}^{2} \f$
!!                 in Eq. (6.4) of Ref. \cite Fermann:2020
!!                 or inside the argument of the incomplete Gamma function
!!                 in Eq. (14) of LIBINT Manual
! *****************************************************************************
  subroutine build_quartet_data(A, B, C, D, m_max, hfx_opts, &
&                prim, do_it, ZetaInv, EtaInv, ZetapEtaInv, Rho, RhoInv, &
&                S1234, P, Q, W, PQ2)

    use units,          only: pi
    use sys,            only: die
    use precision
    use nao2gto_libint
    use nao2gto_utils,  only: fgamma => fgamma_0
    use nao2gto_types,  only: hfx_options_type
    use nao2gto_types,  only: do_hfx_potential_coulomb, &
   &                          do_hfx_potential_short,   &
   &                          do_hfx_potential_truncated

    implicit none

!   Arguments
    real(dp), intent(in)               :: A(3)
                           ! Center of the first primitive Gaussian
    real(dp), intent(in)               :: B(3)
                           ! Center of the second primitive Gaussian
    real(dp), intent(in)               :: C(3)
                           ! Center of the third primitive Gaussian
    real(dp), intent(in)               :: D(3)
                           ! Center of the fourth primitive Gaussian
    integer,  intent(in)               :: m_max
                           ! Maximum \f$ m \f$ value of the
                           !    incomplete Gamma function, defined in
                           !    Eq. (6.4) of Ref. \cite{Fermann:2020}.
                           !    In page 34 of this Reference,
                           !    m is said to run from
                           !    0 \le m \le (la + lb + lc + ld)
                           !    where (la + lb + lc + ld) are the angular
                           !    momenta of the four primitive
                           !    Cartesian Gaussian functions
                           !    (see page 28 of that Reference)
                           !    and page 7 of LIBINT Manual
    type(prim_data), intent(out)       :: prim
                           !    Internal coefficients for Libint and Libderiv
    logical, intent(inout)             :: do_it
    real(dp), intent(in)               :: Rho
                           !   [(\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)]
                           !    that appears in Eq. (6.4) of Ref. Fermann:2020
    real(dp), intent(in)               :: PQ2
                           !   The variable PQ2 corresponds with
                           !   $\overline{\vec{P}\vec{Q}}^{2}$
                           !   in Eq. (6.4) of Ref. Fermann:2020
    real(dp), intent(in)               :: RhoInv 
                           !   RhoInv = (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q)
    real(dp), intent(in)               :: ZetapEtaInv
                           !   ZetapEtaInv = 1.0/(Zeta + Eta),
                           !       where Zeta = Zeta_A + Zeta_B,
                           !       and Eta = Zeta_C + Zeta_D,
                           !   ZetapEtaInv = 1.0/(gamma_p + \gamma_q)
    real(dp), intent(in)               :: S1234
                           !   Product of K_{1}*K_{2} in Eq. (6.4) of 
                           !     Ref. \cite{Fermann:2020}
    real(dp), intent(in)               :: P(3)
                           !   \vec{P}=\frac{\zeta_{a}\vec{A}+\zeta_{b}\vec{B}}
                           !                {\zeta} 
    real(dp), intent(in)               :: Q(3)
                           !   \vec{P}=\frac{\zeta_{c}\vec{C}+\zeta_{d}\vec{D}}
                           !                {\eta} 
    real(dp), intent(in)               :: W(3)
                           !   \vec{W}=\frac{\zeta\vec{P}+\eta\vec{Q}}
                           !                {\zeta + \eta} 


    type(hfx_options_type), intent(in) :: hfx_opts
    real(dp), intent(in)               :: ZetaInv, EtaInv

!   Local variables
    character(len=120)      :: msg
    integer                 :: i
    real(dp)                :: T  ! T = \overline{\vec{P}\vec{Q}}^{2} *
                                  !  (\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)
                                  !  This is the argument of the function
                                  !  F_{m}(T), defined in page 26
                                  !  of Ref. Fermann:2020
    real(dp)                :: factor
                                  ! Prefactor of the incomplete Gamma factor
    real(dp)                :: tmp
                                  ! Temporary variable to compute the prefactor
                                  !   of the incomplete Gamma function
    real(dp)                :: omega2, omega_corr, omega_corr2
   real(dp), dimension(17) :: Fm
                                  ! Incomplete Gamma function
                                  ! F_m (T) =
                                  !   \int_{0}^{1} dt \: t^{2m} \exp(-Tt^{2}}


!-------------------------------------------------------------------------------

!   Factor is the term that appears right after the equal sign in Eq. (6.4)
!   of Ref. \cite Fermann:2020
!   $\frac {2 \pi^{5/2} K_{1} K_{2} }
!          {\gamma_p \gamma_q (\gamma_p + \gamma_q)^{1/2}} $
!   Up to this point, factor = 2 \pi * RhoInv =
!                            = 2 * \pi * (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q
    factor = 2.0_dp*Pi*RhoInv

!   Define the variable of the function F_m, defined in page 26 of
!   Ref. Fermann:2020 or in Eq.(14) of LIBINT Manual
!   (see the argument of F_{m} there).
!   T = \overline{\vec{P}\vec{Q}}^{2} *
!       (\gamma_p * \gamma_q)/(\gamma_p + \gamma_q)
!   or in the notation of the LIBINT Manual
!   T = \overline{\vec{P}\vec{Q}}^{2} * (\zeta * \eta)/(\zeta + \eta)
    T = Rho*PQ2

    do_it = .true.

    select case (hfx_opts%potential_type)

      case(do_hfx_potential_coulomb)
        call fgamma(m_max, T, prim%F(1))

      case(do_hfx_potential_short)
        call fgamma(m_max, T, prim%F(1))
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        call fgamma(m_max,T,Fm)
        tmp = - omega_corr
        do i=1, m_max+1
          prim%F(i) = prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        end do

      case default
        write(msg, fmt='(A,": ",I2)') "HFX potential type not implemented", &
          hfx_opts%potential_type
        call die(msg)

    end select

!   Factor is the term that appears right after the equal sign in Eq. (6.4)
!   of Ref. \cite Fermann:2020
!   $\frac {2 \pi^{5/2} K_{1} K_{2} }
!          {\gamma_p \gamma_q (\gamma_p + \gamma_q)^{1/2}} $
!   Up to this point, factor = 2 \pi * RhoInv =
!                            = 2 * \pi * (\gamma_p+\gamma_q)/(\gamma_p*\gamma_q
!   Now, we multiply it by tmp = (Pi*ZetapEtaInv)**3
!                            = \pi**{3} * 1.0/(gamma_p + \gamma_q)**{3/2}
!   The result is 2 * \pi^{5/2} *
!                 1.0/((\gamma_p*\gamma_q)*(gamma_p + \gamma_q)**{1/2})
!   Finally, we multiply this by K_1* K_{2},
!   where K_{1} = exp(-\alpha_{1} * \alpha_{2} * {\overline{\vec{PQ}}}^{2}
!                    / gamma_p}
!   and K_{2} = exp(-\alpha_{3} * \alpha_{4} * {\overline{\vec{PQ}}}^{2}
!                    / gamma_q}
!   This product of the two K is stored in the variable S1234
    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor * S1234 * dsqrt(tmp)

!   Populate the entries of each component of prim_data, defined
!   as in page 7 of the LIBINT Programmer's Manual
!   prim%F(i) is an auxiliary integral over primitive s-functions,
!   referred to as (0 0 | 0 0)^{(m)} in LIBINT Programmer's Manual
    do i = 1, m_max + 1
      prim%F(i) = prim%F(i)*factor
    end do
    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z     = 0.5_dp*ZetaInv
    prim%oo2n     = 0.5_dp*EtaInv
    prim%oo2zn    = 0.5_dp*ZetapEtaInv
    prim%poz      = Rho*ZetaInv
    prim%pon      = Rho*EtaInv
    prim%oo2p     = 0.5_dp*RhoInv

  end subroutine build_quartet_data

! *****************************************************************************
  SUBROUTINE build_quartet_data_screen(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, &
                                       hfx_opts, prim, do_it, np, nq)

    use sys, only: die
    use units, only: pi
    use precision
    use nao2gto_libint
    use nao2gto_types
    use nao2gto_utils, only: fgamma => fgamma_0

    implicit none

    REAL(KIND=dp)                            :: A(3), B(3), C(3), D(3)
    REAL(KIND=dp), INTENT(IN)                :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                      :: m_max
    TYPE(hfx_options_type)                :: hfx_opts
    TYPE(prim_data)                          :: prim
    LOGICAL, INTENT(INOUT)                   :: do_it
    INTEGER                                  :: np, nq

    character(len=120) :: msg
    LOGICAL :: use_gamma
    INTEGER :: i
    REAL(KIND=dp) :: AB(3), AB2, CD(3), CD2, Eta, EtaInv, factor, omega2, &
      omega_corr, omega_corr2, P(3), PQ(3), PQ2, Q(3), R, R1, R2, Rho, &
      RhoInv, S1234, T, tmp, W(3), Zeta, ZetaInv, ZetapEtaInv
    REAL(KIND=dp), DIMENSION(17)             :: Fm
    Zeta = Zeta_A + Zeta_B
    ZetaInv = 1.0_dp/Zeta
    Eta  = Zeta_C + Zeta_D
    EtaInv = 1.0_dp/Eta
    ZetapEtaInv = Zeta+Eta
    ZetapEtaInv = 1.0_dp/ZetapEtaInv
    Rho  = Zeta*Eta*ZetapEtaInv
    RhoInv = 1.0_dp/Rho

    do i=1,3
      P(i) = (Zeta_A*A(i) + Zeta_B*B(i))*ZetaInv
      Q(i) = (Zeta_C*C(i) + Zeta_D*D(i))*EtaInv
      AB(i) = A(i)-B(i)
      CD(i) = C(i)-D(i)
      PQ(i) = P(i)-Q(i)
      W(i) = (Zeta*P(i) + Eta*Q(i))*ZetapEtaInv
    end do

    AB2 = DOT_PRODUCT(AB,AB)
    CD2 = DOT_PRODUCT(CD,CD)
    PQ2 = DOT_PRODUCT(PQ,PQ)

    factor = 2.0_dp*Pi*RhoInv
    S1234= EXP((-Zeta_A*Zeta_B*ZetaInv*AB2)+(-Zeta_C*Zeta_D*EtaInv*CD2))
    T = Rho*PQ2

    do_it = .TRUE.

    select case(hfx_opts%potential_type)

      case (do_hfx_potential_coulomb)
        CALL fgamma(m_max,T,prim%F)

      case (do_hfx_potential_short)
        CALL fgamma(m_max,T,prim%F)
        omega2 = hfx_opts%omega**2
        omega_corr2 = omega2/(omega2+Rho)
        omega_corr = SQRT(omega_corr2)
        T = T*omega_corr2
        CALL fgamma(m_max,T,Fm)
        tmp = - omega_corr
        do i=1,m_max+1
          prim%F(i)=prim%F(i) + Fm(i)*tmp
          tmp = tmp * omega_corr2
        end do

      case default
        write(msg, fmt='(A,": ",I2)') "HFX potential type not implemented", &
          hfx_opts%potential_type
        call die(msg)

    end select

    tmp    = (Pi*ZetapEtaInv)**3
    factor = factor*S1234*SQRT(tmp)

    do i=1,m_max+1
       prim%F(i) = prim%F(i)*factor
    end do
    prim%U(1:3,1) = P-A
    prim%U(1:3,3) = Q-C
    prim%U(1:3,5) = W-P
    prim%U(1:3,6) = W-Q
    prim%oo2z      = 0.5_dp*ZetaInv
    prim%oo2n      = 0.5_dp*EtaInv
    prim%oo2zn     = 0.5_dp*ZetapEtaInv
    prim%poz       = Rho*ZetaInv
    prim%pon       = Rho*EtaInv
    prim%oo2p      = 0.5_dp*RhoInv

  END SUBROUTINE build_quartet_data_screen

end module nao2gto_primitive
