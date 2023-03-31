! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!< Bloch expansion module implementing a way to expand a Hermitian matrix obeying Bloch's theorem.
!!
!! Bloch's theorem states that only a constant phase-factor is accummulated upon a translation
!! vector `T`.
!!
!! This module implements a type `bloch_unfold_t` implementing the required algorithms
!! for unfolding a matrix.
module bloch_unfold_m

  implicit none
  private

  integer, parameter :: dp = selected_real_kind(p=15)

  type bloch_unfold_t

    !< Number of times a matrix is unfolded
    integer, public :: B(3) = 1

    !< Local variable for quick query of product(this%B)
    integer, private :: prod_B = 1

  contains

    procedure, pass :: initialize => bloch_unfold_init
    procedure, pass :: size => bloch_unfold_size
#ifdef __BLOCH_UNFOLD_M_TEST
    procedure, pass :: unfold_M_original => bloch_unfold_M_original
    procedure, pass :: unfold_M_do => bloch_unfold_M_do
    procedure, pass :: unfold_M_manual => bloch_unfold_M_manual
    procedure, pass :: unfold_M_manual_simd => bloch_unfold_M_manual_simd
    procedure, pass :: unfold_M_task => bloch_unfold_M_task
    procedure, pass :: unfold_M_task_bad => bloch_unfold_M_task_bad
    procedure, pass :: unfold_M_workshare => bloch_unfold_M_workshare
#else
# ifdef _OPENMP
    procedure, pass :: unfold_M => bloch_unfold_M_manual
# else
    procedure, pass :: unfold_M => bloch_unfold_M_do
# endif
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
    procedure, pass :: unfold_HS_G_original => bloch_unfold_HS_G_original
    procedure, pass :: unfold_HS_G_do => bloch_unfold_HS_G_do
    procedure, pass :: unfold_HS_G_manual => bloch_unfold_HS_G_manual
#else
# ifdef _OPENMP
    procedure, pass :: unfold_HS_G => bloch_unfold_HS_G_manual
# else
    procedure, pass :: unfold_HS_G => bloch_unfold_HS_G_do
# endif
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
    procedure, pass :: unfold_HS_original => bloch_unfold_HS_original
    procedure, pass :: unfold_HS_do => bloch_unfold_HS_do
    procedure, pass :: unfold_HS_manual => bloch_unfold_HS_manual
#else
# ifdef _OPENMP
    procedure, pass :: unfold_HS => bloch_unfold_HS_manual
# else
    procedure, pass :: unfold_HS => bloch_unfold_HS_do
# endif
#endif
    procedure, pass :: get_k => bloch_unfold_get_k
    procedure, pass :: unravel_index => bloch_unfold_unravel_index

  end type bloch_unfold_t
  public :: bloch_unfold_t

#ifdef __BLOCH_UNFOLD_M_TEST
  public :: unfold_k2pi
#endif

contains

  subroutine bloch_unfold_init(this, B)
    class(bloch_unfold_t), intent(inout) :: this
    integer, intent(in) :: B(3)

    this%B(:) = B(:)

    if ( any(this%B < 1) ) then
      call die('bloch_unfold requires unfolding to be >= 1!')
    end if

    this%prod_B = product(this%B)

  end subroutine bloch_unfold_init

  pure function bloch_unfold_size(this) result(size)
    class(bloch_unfold_t), intent(in) :: this
    integer :: size

    size = this%prod_B

  end function bloch_unfold_size


#ifdef __BLOCH_UNFOLD_M_TEST
  subroutine bloch_unfold_M_original(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%B(1),this%B(2),this%B(3))
    complex(dp), intent(inout) :: uM(N,this%B(1),this%B(2),this%B(3),N,this%B(1),this%B(2),this%B(3))

    integer :: B(3), i1, i2, i3
    integer :: ia1, ia2, ia3,iuo
    integer :: ja1, ja2, ja3,juo
    real(dp) :: wq, rPi(3)
    complex(dp) :: qPi, p(3)

!$OMP parallel default(shared) private(wq,rPi,qPi,p) &
!$OMP&  private(B,i1,i2,i3) &
!$OMP&  private(ia1,ia2,ia3,iuo) &
!$OMP&  private(ja1,ja2,ja3,juo)

!$OMP workshare
    uM(:,:,:,:,:,:,:,:) = 0._dp
!$OMP end workshare

  ! Save some multiplications
    B(:) = this%B(:)
    wq = log(1._dp / real(product(B),dp))

    do i3 = 1, B(3)
      rPi(3) = unfold_k2pi(bk(3), B(3), i3)
      do i2 = 1, B(2)
        rPi(2) = unfold_k2pi(bk(2), B(2), i2)
        do i1 = 1, B(1)
          rPi(1) = unfold_k2pi(bk(1), B(1), i1)

          qPi = exp(cmplx(0._dp,rPi(1), dp))

!$OMP do collapse(3) schedule(static)
          do ia3 = 1 , B(3)
            do ia2 = 1 , B(2)
              do ia1 = 1 , B(1)

                p(3) = exp(cmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3), dp))
                do iuo = 1, N
                  do ja3 = 1 , B(3)
                    p(2) = p(3)*exp(cmplx(0._dp,ja3*rPi(3), dp))
                    do ja2 = 1 , B(2)
                      p(1) = p(2)*exp(cmplx(0._dp,ja2*rPi(2), dp))
                      do ja1 = 1 , B(1)
                        p(1) = p(1)*qPi
                        do juo = 1, N
                          uM(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) = &
                              uM(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) + p(1) * M(juo,iuo,i1,i2,i3)
                        end do !juo
                      end do !ja1
                    end do !ja2
                  end do !ja3
                end do !iuo
              end do !ia1
            end do !ia2
          end do !ia3
!$OMP end do

        end do !i1
      end do !i2
    end do !i3

!$OMP end parallel

  end subroutine bloch_unfold_M_original
#endif

#if defined (__BLOCH_UNFOLD_M_TEST) || ! defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_task_z` but uses
  !! the do loop construct for OpenMP parallelization.
  subroutine bloch_unfold_M_do(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_do_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_do_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_do_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_do_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_do_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_do_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_do

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_do_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,cph)

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0, init="M_do")
#endif

    ! Initialize un-folded matrix
!$OMP do schedule(static)
    do j = 1, N
      do iMA = 1, NA
        do i = 1, N
          uM(i,iMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( j == 1 ) call print_matrix(NA,iMA,1, c='0')
#endif
      end do
    end do
!$OMP end do nowait
!$OMP do schedule(static)
    do j = 1, N
      do iMA = 2, NA
        do i = 1, N
          uM(i,1,j,iMA) = cmplx(0._dp, 0._dp, dp)
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( j == 1 ) call print_matrix(NA,1,iMA, c='0')
#endif
      end do
    end do
!$OMP end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0, out=.true.)
#endif

    ! Full weight
    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
      if ( TMA == 1 ) call print_matrix(NA,1,1)
#endif

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        cph = conjg(ph)
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * cph
          end do
        end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( TMA == 1 ) call print_matrix(NA,iMA,1)
        if ( TMA == 1 ) call print_matrix(NA,1,iMA)
#endif
        ph = ph * phA_step
      end do

!$OMP barrier

    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0, out=.true.)
#endif

    ! At this point the following has been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do TMA = 2, NA
      do iMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uM(i,iMA,j,TMA) = uM(i,iMA-1,j,TMA-1)
          end do
        end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
        call print_matrix(NA,iMA,TMA, c='-')
#endif
      end do
!$OMP barrier
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0, out=.true.)
#endif

!$OMP end parallel

  end subroutine bloch_unfold_M_do_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_do_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, cph

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,cph)

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, init="M_do")
#endif

    ! Initialize un-folded matrix
!$OMP do schedule(static)
    do j = 1, N
      do i = 1, N
        uM(i,1,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,1,1,NB,1,1, c='0')
#endif

    do iMA = 2, NA
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uM(i,iMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uM(i,1,1,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
      call print_matrix(NA,iMA,1,NB,1,1, c='0')
      call print_matrix(NA,1,iMA,NB,1,1, c='0')
#endif
    end do

    do iMB = 2, NB
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uM(i,1,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uM(i,1,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
      call print_matrix(NA,1,1,NB,iMB,1, c='0')
      call print_matrix(NA,1,1,NB,1,iMB, c='0')
#endif

      do iMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uM(i,iMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,iMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,1,iMB,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,1,1,j,iMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
        call print_matrix(NA,iMA,1,NB,iMB,1, c='0')
        call print_matrix(NA,iMA,1,NB,1,iMB, c='0')
        call print_matrix(NA,1,iMA,NB,iMB,1, c='0')
        call print_matrix(NA,1,iMA,NB,1,iMB, c='0')
#endif
      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

!$OMP barrier

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,1)
#endif

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA
          ph = phA
          cph = conjg(phA)
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,1)
#endif
          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          cph = conjg(phB)
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,iMB,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,iMB)
#endif

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA
            ph = phB * phA
            cph = phB * conjg(phA)
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,iMB,1)
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,iMB,1)
#endif

            ph = conjg(phB) * phA
            cph = conjg(phB * phA)
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,iMB)
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,iMB)
#endif

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

        ! Since all above do loops have nowait, we really have
        ! to be sure we don't begin the next iteration, so we have to have
        ! an explicit barrier. Otherwise we may end up with a data-race
!$OMP barrier

      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

    ! Now we have for all NB calculated the
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! components. Now fill each NA block

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do TMA = 2, NA
      do iMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uM(i,iMA,1,j,TMA,1) = uM(i,iMA-1,1,j,TMA-1,1)
          end do
        end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
        call print_matrix(NA,iMA,TMA,NB,1,1, c='-')
#endif
      end do
      ! we need a wait since it copies from the previous iMA
!$OMP barrier
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do iMA = 2, NA
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uM(i,iMA,TMB,j,TMA,1) = uM(i,iMA-1,TMB,j,TMA-1,1)
            end do
          end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
          call print_matrix(NA,iMA,TMA,NB,TMB,1, c='-')
#endif
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uM(i,iMA,1,j,TMA,TMB) = uM(i,iMA-1,1,j,TMA-1,TMB)
            end do
          end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
          call print_matrix(NA,iMA,TMA,NB,1,TMB, c='-')
#endif
        end do
      end do
    end do

!$OMP barrier

    ! Now we have filled all
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! for NB blocks
    do TMB = 2, NB
      do TMA = 1, NA
        do iMB = 2, NB
          do iMA = 1, NA
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uM(i,iMA,iMB,j,TMA,TMB) = uM(i,iMA,iMB-1,j,TMA,TMB-1)
              end do
            end do
!$OMP end do nowait
#ifdef __BLOCH_UNFOLD_M_DEBUG
            call print_matrix(NA,iMA,TMA,NB,iMB,TMB, c='-')
#endif
          end do
        end do
      end do
!$OMP barrier
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

!$OMP end parallel

  end subroutine bloch_unfold_M_do_2
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_do` but uses
  !! the task construct for OpenMP parallelization.
  !! Currently, this method is slower.
  subroutine bloch_unfold_M_task(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_task_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_task_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_task_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_task_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_task_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_task_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_task

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_task_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) private(w)

    ! Full weight
    w = 1._dp / real(NA, dp)

!$OMP single private(TMA,iMA,k_A,phA_step,ph)

!$OMP task shared(uM) private(i,j) &
!$OMP&  depend(out:uM(1,1,1,1))
    do j = 1, N
      do i = 1, N
        uM(i,1,j,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end task

    do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(1,iMA,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(1,1,1,iMA))
      do j = 1, N
        do i = 1, N
          uM(i,1,j,iMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
    end do

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,w) &
!$OMP&  depend(inout:uM(1,1,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do
!$OMP end task

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j) firstprivate(iMA,TMA,ph) &
!$OMP&  depend(inout:uM(1,iMA,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
        end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,cph) firstprivate(iMA,TMA,ph) &
!$OMP&  depend(inout:uM(1,1,1,iMA))
        cph = conjg(ph)
        do j = 1, N
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * cph
          end do
        end do
!$OMP end task

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA

!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(in:uM(1,1,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,j,iMA) = uM(i,1,j,1)
        end do
      end do
!$OMP end task

      do TMA = 1, NA - iMA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(1,iMA,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,TMA+iMA,j,TMA+1) = uM(i,iMA,j,1)
          end do
        end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(1,1,1,iMA))
        do j = 1, N
          do i = 1, N
            uM(i,TMA+1,j,TMA+iMA) = uM(i,1,j,iMA)
          end do
        end do
!$OMP end task
      end do

    end do

!$OMP end single nowait

!$OMP end parallel

  end subroutine bloch_unfold_M_task_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_task_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph

!$OMP parallel default(shared) private(w)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

!$OMP single private(TMA,iMA,k_A,phA_step,phA,TMB,iMB,k_B,phB_step,phB,ph)

!$OMP task shared(uM) private(i,j) &
!$OMP&  depend(out:uM(1,1,1,1,1,1))
    do j = 1, N
      do i = 1, N
        uM(i,1,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end task

    do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(1,iMA,1,1,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(1,1,1,1,iMA,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,1,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
    end do

    do iMB = 2, NB
!$OMP task shared(uM) private(i,j) firstprivate(iMB) &
!$OMP&  depend(out:uM(1,1,iMB,1,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMB) &
!$OMP&  depend(out:uM(1,1,1,1,1,iMB))
      do j = 1, N
        do i = 1, N
          uM(i,1,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task

      do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(1,iMA,iMB,1,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(1,iMA,1,1,1,iMB))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(1,1,iMB,1,iMA,1))
        do j = 1, N
          do i = 1, N
            uM(i,1,iMB,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(1,1,1,1,iMA,iMB))
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,iMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
      end do
    end do

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,TMB,w) &
!$OMP&  depend(inout:uM(1,1,1,1,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end task

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j) firstprivate(iMA,TMA,TMB,phA) &
!$OMP&  depend(inout:uM(1,iMA,1,1,1,1))
          do j = 1, N
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * phA
            end do
          end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,TMB,phA) &
!$OMP&  depend(inout:uM(1,1,1,1,iMA,1))
          ph = conjg(phA)
          do j = 1, N
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * ph
            end do
          end do
!$OMP end task

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,iMB,TMB,phB) &
!$OMP&  depend(inout:uM(1,1,iMB,1,1,1))
          do j = 1, N
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * phB
            end do
          end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(TMA,iMB,TMB,phB) &
!$OMP&  depend(inout:uM(1,1,1,1,1,iMB))
          ph = conjg(phB)
          do j = 1, N
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
            end do
          end do
!$OMP end task

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(1,iMA,iMB,1,1,1))
            ph = phB * phA
            do j = 1, N
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(1,1,iMB,1,iMA,1))
            ph = phB * conjg(phA)
            do j = 1, N
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(1,iMA,1,1,1,iMB))
            ph = conjg(phB) * phA
            do j = 1, N
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(1,1,1,1,iMA,iMB))
            ph = conjg(phB * phA)
            do j = 1, N
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do TMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(1,TMA-1,1,1,iMA-1,1)) &
!$OMP&  depend(out:uM(1,TMA,1,1,iMA,1))
        do j = 1, N
          do i = 1, N
            uM(i,TMA,1,j,iMA,1) = uM(i,TMA-1,1,j,iMA-1,1)
          end do
        end do
!$OMP end task
      end do
    end do

    do iMB = 2, NB
      do iMA = 2, NA
        do TMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB) &
!$OMP&  depend(in:uM(1,TMA-1,iMB,1,iMA-1,1)) &
!$OMP&  depend(out:uM(1,TMA,iMB,1,iMA,1))
          do j = 1, N
            do i = 1, N
              uM(i,TMA,iMB,j,iMA,1) = uM(i,TMA-1,iMB,j,iMA-1,1)
            end do
          end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB) &
!$OMP&  depend(in:uM(1,TMA-1,1,1,iMA-1,iMB)) &
!$OMP&  depend(out:uM(1,TMA,1,1,iMA,iMB))
          do j = 1, N
            do i = 1, N
              uM(i,TMA,1,j,iMA,iMB) = uM(i,TMA-1,1,j,iMA-1,iMB)
            end do
          end do
!$OMP end task
        end do
      end do
    end do

    do iMB = 2, NB
      do iMA = 1, NA
        do TMB = 2, NB
          do TMA = 1, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB,TMB) &
!$OMP&  depend(in:uM(1,TMA,TMB-1,1,iMA,iMB-1)) &
!$OMP&  depend(out:uM(1,TMA,TMB,1,iMA,iMB))
            do j = 1, N
              do i = 1, N
                uM(i,TMA,TMB,j,iMA,iMB) = uM(i,TMA,TMB-1,j,iMA,iMB-1)
              end do
            end do
!$OMP end task
          end do
        end do
      end do
    end do

!$OMP end single nowait

!$OMP end parallel

  end subroutine bloch_unfold_M_task_2

  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_task` but uses
  !! bad dependency lists.
  subroutine bloch_unfold_M_task_bad(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_task_bad_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_task_bad_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_task_bad_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_task_bad_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_task_bad_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_task_bad_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_task_bad

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_task_bad_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) private(w)

    ! Full weight
    w = 1._dp / real(NA, dp)

!$OMP single private(TMA,iMA,k_A,phA_step,ph)

!$OMP task shared(uM) private(i,j) &
!$OMP&  depend(out:uM(:,1,:,1))
    do j = 1, N
      do i = 1, N
        uM(i,1,j,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end task

    do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(:,iMA,:,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(:,1,:,iMA))
      do j = 1, N
        do i = 1, N
          uM(i,1,j,iMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
    end do

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,w) &
!$OMP&  depend(inout:uM(:,1,:,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do
!$OMP end task

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j) firstprivate(iMA,TMA,ph) &
!$OMP&  depend(inout:uM(:,iMA,:,1))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
        end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,cph) firstprivate(iMA,TMA,ph) &
!$OMP&  depend(inout:uM(:,1,:,iMA))
        cph = conjg(ph)
        do j = 1, N
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * cph
          end do
        end do
!$OMP end task

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA

!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(in:uM(:,1,:,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,j,iMA) = uM(i,1,j,1)
        end do
      end do
!$OMP end task

      do TMA = 1, NA - iMA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(:,iMA,:,1))
        do j = 1, N
          do i = 1, N
            uM(i,TMA+iMA,j,TMA+1) = uM(i,iMA,j,1)
          end do
        end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(:,1,:,iMA))
        do j = 1, N
          do i = 1, N
            uM(i,TMA+1,j,TMA+iMA) = uM(i,1,j,iMA)
          end do
        end do
!$OMP end task
      end do

    end do

!$OMP end single nowait

!$OMP end parallel

  end subroutine bloch_unfold_M_task_bad_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_task_bad_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph

!$OMP parallel default(shared) private(w)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

!$OMP single private(TMA,iMA,k_A,phA_step,phA,TMB,iMB,k_B,phB_step,phB,ph)

!$OMP task shared(uM) private(i,j) &
!$OMP&  depend(out:uM(:,1,1,:,1,1))
    do j = 1, N
      do i = 1, N
        uM(i,1,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end task

    do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(:,iMA,1,:,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,iMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA) &
!$OMP&  depend(out:uM(:,1,1,:,iMA,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,1,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
    end do

    do iMB = 2, NB
!$OMP task shared(uM) private(i,j) firstprivate(iMB) &
!$OMP&  depend(out:uM(:,1,iMB,:,1,1))
      do j = 1, N
        do i = 1, N
          uM(i,1,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMB) &
!$OMP&  depend(out:uM(:,1,1,:,1,iMB))
      do j = 1, N
        do i = 1, N
          uM(i,1,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end task

      do iMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(:,iMA,iMB,:,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(:,iMA,1,:,1,iMB))
        do j = 1, N
          do i = 1, N
            uM(i,iMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(:,1,iMB,:,iMA,1))
        do j = 1, N
          do i = 1, N
            uM(i,1,iMB,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
!$OMP task shared(uM) private(i,j) firstprivate(iMA,iMB) &
!$OMP&  depend(out:uM(:,1,1,:,iMA,iMB))
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,iMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end task
      end do
    end do

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,TMB,w) &
!$OMP&  depend(inout:uM(:,1,1,:,1,1))
        do j = 1, N
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end task

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j) firstprivate(iMA,TMA,TMB,phA) &
!$OMP&  depend(inout:uM(:,iMA,1,:,1,1))
          do j = 1, N
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * phA
            end do
          end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,TMB,phA) &
!$OMP&  depend(inout:uM(:,1,1,:,iMA,1))
          ph = conjg(phA)
          do j = 1, N
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * ph
            end do
          end do
!$OMP end task

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
!$OMP task shared(uM,M) private(i,j) firstprivate(TMA,iMB,TMB,phB) &
!$OMP&  depend(inout:uM(:,1,iMB,:,1,1))
          do j = 1, N
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * phB
            end do
          end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(TMA,iMB,TMB,phB) &
!$OMP&  depend(inout:uM(:,1,1,:,1,iMB))
          ph = conjg(phB)
          do j = 1, N
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
            end do
          end do
!$OMP end task

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(:,iMA,iMB,:,1,1))
            ph = phB * phA
            do j = 1, N
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(:,1,iMB,:,iMA,1))
            ph = phB * conjg(phA)
            do j = 1, N
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(:,iMA,1,:,1,iMB))
            ph = conjg(phB) * phA
            do j = 1, N
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

!$OMP task shared(uM,M) private(i,j,ph) firstprivate(iMA,TMA,iMB,TMB,phA,phB) &
!$OMP&  depend(inout:uM(:,1,1,:,iMA,iMB))
            ph = conjg(phB * phA)
            do j = 1, N
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
!$OMP end task

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do TMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA) &
!$OMP&  depend(in:uM(:,TMA-1,1,:,iMA-1,1)) &
!$OMP&  depend(out:uM(:,TMA,1,:,iMA,1))
        do j = 1, N
          do i = 1, N
            uM(i,TMA,1,j,iMA,1) = uM(i,TMA-1,1,j,iMA-1,1)
          end do
        end do
!$OMP end task
      end do
    end do

    do iMB = 2, NB
      do iMA = 2, NA
        do TMA = 2, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB) &
!$OMP&  depend(in:uM(:,TMA-1,iMB,:,iMA-1,1)) &
!$OMP&  depend(out:uM(:,TMA,iMB,:,iMA,1))
          do j = 1, N
            do i = 1, N
              uM(i,TMA,iMB,j,iMA,1) = uM(i,TMA-1,iMB,j,iMA-1,1)
            end do
          end do
!$OMP end task

!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB) &
!$OMP&  depend(in:uM(:,TMA-1,1,:,iMA-1,iMB)) &
!$OMP&  depend(out:uM(:,TMA,1,:,iMA,iMB))
          do j = 1, N
            do i = 1, N
              uM(i,TMA,1,j,iMA,iMB) = uM(i,TMA-1,1,j,iMA-1,iMB)
            end do
          end do
!$OMP end task
        end do
      end do
    end do

    do iMB = 2, NB
      do iMA = 1, NA
        do TMB = 2, NB
          do TMA = 1, NA
!$OMP task shared(uM) private(i,j) firstprivate(iMA,TMA,iMB,TMB) &
!$OMP&  depend(in:uM(:,TMA,TMB-1,:,iMA,iMB-1)) &
!$OMP&  depend(out:uM(:,TMA,TMB,:,iMA,iMB))
            do j = 1, N
              do i = 1, N
                uM(i,TMA,TMB,j,iMA,iMB) = uM(i,TMA,TMB-1,j,iMA,iMB-1)
              end do
            end do
!$OMP end task
          end do
        end do
      end do
    end do

!$OMP end single nowait

!$OMP end parallel

  end subroutine bloch_unfold_M_task_bad_2
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_do` but uses
  !! the workshare construct for OpenMP parallelization.
  !! Currently, this method is slower.
  subroutine bloch_unfold_M_workshare(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_workshare_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_workshare_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_workshare_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_workshare_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_workshare_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_workshare_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_workshare

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_workshare_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) private(w,k_A,ph,cph,phA_step,TMA,iMA)

!$OMP workshare
    uM(:,:,:,1) = cmplx(0._dp, 0._dp, dp)
    uM(:,1,:,2:) = cmplx(0._dp, 0._dp, dp)
!$OMP end workshare

    ! Full weight
    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP workshare
      uM(:,1,:,1) = uM(:,1,:,1) + M(:,:,TMA) * w
!$OMP end workshare nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

        cph = conjg(ph)
!$OMP workshare
        uM(:,iMA,:,1) = uM(:,iMA,:,1) + M(:,:,TMA) * ph
        uM(:,1,:,iMA) = uM(:,1,:,iMA) + M(:,:,TMA) * cph
!$OMP end workshare nowait

        ph = ph * phA_step
      end do

!$OMP barrier
    end do

    ! At this point the following have been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA

!$OMP workshare
      uM(:,iMA,:,iMA) = uM(:,1,:,1)
!$OMP end workshare nowait

      do TMA = 1, NA - iMA
!$OMP workshare
        uM(:,TMA+iMA,:,TMA+1) = uM(:,iMA,:,1)
        uM(:,TMA+1,:,TMA+iMA) = uM(:,1,:,iMA)
!$OMP end workshare nowait
      end do

    end do

!$OMP end parallel

  end subroutine bloch_unfold_M_workshare_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_workshare_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, ph1, ph2, ph3

!$OMP parallel default(shared) private(w) &
!$OMP&   private(TMA,iMA,k_A,phA_step,phA) &
!$OMP&   private(TMB,iMB,k_B,phB_step,phB) &
!$OMP&   private(ph,ph1,ph2,ph3)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

!$OMP workshare
    uM(:,:,:,:,1,1) = cmplx(0._dp, 0._dp, dp)
    uM(:,:,1,:,1,2:) = cmplx(0._dp, 0._dp, dp)
    uM(:,1,:,:,2:,1) = cmplx(0._dp, 0._dp, dp)
    uM(:,1,1,:,2:,2:) = cmplx(0._dp, 0._dp, dp)
!$OMP end workshare

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

!$OMP workshare
        uM(:,1,1,:,1,1) = uM(:,1,1,:,1,1) + M(:,:,TMA,TMB) * w
!$OMP end workshare nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

          ph = conjg(phA)
!$OMP workshare
          uM(:,iMA,1,:,1,1) = uM(:,iMA,1,:,1,1) + M(:,:,TMA,TMB) * phA
          uM(:,1,1,:,iMA,1) = uM(:,1,1,:,iMA,1) + M(:,:,TMA,TMB) * ph
!$OMP end workshare nowait

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = conjg(phB)
!$OMP workshare
          uM(:,1,iMB,:,1,1) = uM(:,1,iMB,:,1,1) + M(:,:,TMA,TMB) * phB
          uM(:,1,1,:,1,iMB) = uM(:,1,1,:,1,iMB) + M(:,:,TMA,TMB) * ph
!$OMP end workshare nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            ph1 = phB * conjg(phA)
            ph2 = conjg(phB) * phA
            ph3 = conjg(phB * phA)

!$OMP workshare
            uM(:,iMA,iMB,:,1,1) = uM(:,iMA,iMB,:,1,1) + M(:,:,TMA,TMB) * ph
            uM(:,1,iMB,:,iMA,1) = uM(:,1,iMB,:,iMA,1) + M(:,:,TMA,TMB) * ph1
            uM(:,iMA,1,:,1,iMB) = uM(:,iMA,1,:,1,iMB) + M(:,:,TMA,TMB) * ph2
            uM(:,1,1,:,iMA,iMB) = uM(:,1,1,:,iMA,iMB) + M(:,:,TMA,TMB) * ph3
!$OMP end workshare nowait

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do TMA = 2, NA
      do iMA = 2, NA
!$OMP workshare
        uM(:,iMA,1,:,TMA,1) = uM(:,iMA-1,1,:,TMA-1,1)
!$OMP end workshare nowait
      end do
!$OMP barrier
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do iMA = 2, NA
!$OMP workshare
          uM(:,iMA,TMB,:,TMA,1) = uM(:,iMA-1,TMB,:,TMA-1,1)
          uM(:,iMA,1,:,TMA,TMB) = uM(:,iMA-1,1,:,TMA-1,TMB)
!$OMP end workshare nowait
        end do
      end do
    end do

!$OMP barrier

    do TMB = 2, NB
      do iMB = 2, NB
!$OMP workshare
        uM(:,:,iMB,:,:,TMB) = uM(:,:,iMB-1,:,:,TMB-1)
!$OMP end workshare nowait
      end do
!$OMP barrier
    end do

!$OMP end parallel

  end subroutine bloch_unfold_M_workshare_2
#endif

#if defined(__BLOCH_UNFOLD_M_TEST) || defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_do` but uses
  !! the manual parallelization constructs for OpenMP parallelization.
  !! Currently, this method is faster.
  subroutine bloch_unfold_M_manual(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_manual_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_manual_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_manual_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_manual_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_manual_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_manual_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_manual

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_manual_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,k_A,phA_step,ph,cph,TMA,iMA) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uM(i,TMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uM(i,1,j,TMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      do j = jmin, jmax
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

        cph = conjg(ph)
        do j = jmin, jmax
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * cph
          end do
        end do

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do j = jmin, jmax
        do TMA = 2, NA
          do i = 1, N
            uM(i,TMA,j,iMA) = uM(i,TMA-1,j,iMA-1)
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_M_manual_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_manual_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph) &
!$OMP&   private(k_A,phA_step,phA,TMA,iMA) &
!$OMP&   private(k_B,phB_step,phB,TMB,iMB) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, init="M_manual")
#endif

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uM(i,TMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( j == jmin ) call print_matrix(NA,TMA,1,NB,1,1, c='0')
#endif
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uM(i,1,1,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
      call print_matrix(NA,1,TMA,NB,1,1, c='0')
#endif
    end do

    do iMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do i = 1, N
            uM(i,TMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,1,1,j,TMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,TMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uM(i,1,iMB,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        call print_matrix(NA,TMA,1,NB,iMB,1, c='0')
        call print_matrix(NA,1,TMA,NB,1,iMB, c='0')
        call print_matrix(NA,TMA,1,NB,1,iMB, c='0')
        call print_matrix(NA,1,TMA,NB,iMB,1, c='0')
#endif
      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        do j = jmin, jmax
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,1)
#endif

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

          ph = conjg(phA)

          do j = jmin, jmax
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * phA
            end do
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * ph
            end do
          end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,1)
#endif

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ph = conjg(phB)
          ! Diagonal B (A has zero phase)
          do j = jmin, jmax
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * phB
            end do
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
            end do
          end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,iMB,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,iMB)
#endif

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            do j = jmin, jmax
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,iMB,1)
#endif

            ph = phB * conjg(phA)
            do j = jmin, jmax
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,iMB,1)
#endif

            ph = conjg(phB) * phA
            do j = jmin, jmax
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,iMB)
#endif

            ph = conjg(phB * phA)
            do j = jmin, jmax
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,iMB)
#endif

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do


    do TMA = 2, NA
      do j = jmin, jmax
        do iMA = 2, NA
          do i = 1, N
            uM(i,iMA,1,j,TMA,1) = uM(i,iMA-1,1,j,TMA-1,1)
          end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,1,1, c='-')
#endif
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do j = jmin, jmax
          do iMA = 2, NA
            do i = 1, N
              uM(i,iMA,TMB,j,TMA,1) = uM(i,iMA-1,TMB,j,TMA-1,1)
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,TMB,1, c='-')
#endif

            do i = 1, N
              uM(i,iMA,1,j,TMA,TMB) = uM(i,iMA-1,1,j,TMA-1,TMB)
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,1,TMB, c='-')
#endif
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do iMB = 2, NB
            do iMA = 1, NA
              do i = 1, N
                uM(i,iMA,iMB,j,TMA,TMB) = uM(i,iMA,iMB-1,j,TMA,TMB-1)
              end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
              if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,iMB,TMB, c='-')
#endif
            end do
          end do
        end do
      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

!$OMP end parallel

  end subroutine bloch_unfold_M_manual_2
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  !!
  !! This method is equivalent to `bloch_unfold_M_do` but uses
  !! the manual parallelization constructs for OpenMP parallelization.
  !! Currently, this method is faster.
  subroutine bloch_unfold_M_manual_simd(this, bk, N, M, uM)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), intent(in) :: M(N,N,this%prod_B)
    complex(dp), intent(inout) :: uM(N,this%prod_B,N,this%prod_B)

    if ( this%prod_B == 1 ) then

      call zcopy(N*N, M(1,1,1), 1, uM(1,1,1,1), 1)
      return

    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_M_manual_simd_1(N, M, bk(3), this%B(3), uM)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_M_manual_simd_1(N, M, bk(2), this%B(2), uM)
      else
        call bloch_unfold_M_manual_simd_2(N, M, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_M_manual_simd_1(N, M, bk(1), this%B(1), uM)
      else
        call bloch_unfold_M_manual_simd_2(N, M, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uM)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_M_manual_simd_2(N, M, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uM)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_M_manual_simd

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_manual_simd_1(N, M, kA, NA, uM)
    integer, intent(in) :: N, NA
    complex(dp), intent(in) :: M(N,N,NA)
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uM(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phA_step

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,k_A,phA_step,ph,cph,TMA,iMA) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
!$OMP simd collapse(2)
      do TMA = 1, NA
        do i = 1, N
          uM(i,TMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
!$OMP simd
        do i = 1, N
          uM(i,1,j,TMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      do j = jmin, jmax
!$OMP simd
        do i = 1, N
          uM(i,1,j,1) = uM(i,1,j,1) + M(i,j,TMA) * w
        end do
      end do

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

        cph = conjg(ph)
        do j = jmin, jmax
!$OMP simd
          do i = 1, N
            uM(i,iMA,j,1) = uM(i,iMA,j,1) + M(i,j,TMA) * ph
          end do
!$OMP simd
          do i = 1, N
            uM(i,1,j,iMA) = uM(i,1,j,iMA) + M(i,j,TMA) * cph
          end do
        end do

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do j = jmin, jmax
!$OMP simd collapse(2)
        do TMA = 2, NA
          do i = 1, N
            uM(i,TMA,j,iMA) = uM(i,TMA-1,j,iMA-1)
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_M_manual_simd_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_M_manual_simd_2(N, M, kA, NA, kB, NB, uM)
    integer, intent(in) :: N, NA, NB
    complex(dp), intent(in) :: M(N,N,NA,NB)
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uM(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph) &
!$OMP&   private(k_A,phA_step,phA,TMA,iMA) &
!$OMP&   private(k_B,phB_step,phB,TMB,iMB) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, init="M_manual")
#endif

    ! Initialize
    do j = jmin, jmax
!$OMP simd collapse(2)
      do TMA = 1, NA
        do i = 1, N
          uM(i,TMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( j == jmin ) call print_matrix(NA,TMA,1,NB,1,1, c='0')
#endif
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
!$OMP simd
        do i = 1, N
          uM(i,1,1,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
      call print_matrix(NA,1,TMA,NB,1,1, c='0')
#endif
    end do

    do iMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
!$OMP simd
          do i = 1, N
            uM(i,TMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
!$OMP simd
          do i = 1, N
            uM(i,1,1,j,TMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
!$OMP simd
          do i = 1, N
            uM(i,TMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
!$OMP simd
          do i = 1, N
            uM(i,1,iMB,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        call print_matrix(NA,TMA,1,NB,iMB,1, c='0')
        call print_matrix(NA,1,TMA,NB,1,iMB, c='0')
        call print_matrix(NA,TMA,1,NB,1,iMB, c='0')
        call print_matrix(NA,1,TMA,NB,iMB,1, c='0')
#endif
      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        do j = jmin, jmax
!$OMP simd
          do i = 1, N
            uM(i,1,1,j,1,1) = uM(i,1,1,j,1,1) + M(i,j,TMA,TMB) * w
          end do
        end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
        if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,1)
#endif

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

          ph = conjg(phA)

          do j = jmin, jmax
!$OMP simd
            do i = 1, N
              uM(i,iMA,1,j,1,1) = uM(i,iMA,1,j,1,1) + M(i,j,TMA,TMB) * phA
            end do
!$OMP simd
            do i = 1, N
              uM(i,1,1,j,iMA,1) = uM(i,1,1,j,iMA,1) + M(i,j,TMA,TMB) * ph
            end do
          end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,1)
#endif

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ph = conjg(phB)
          ! Diagonal B (A has zero phase)
          do j = jmin, jmax
!$OMP simd
            do i = 1, N
              uM(i,1,iMB,j,1,1) = uM(i,1,iMB,j,1,1) + M(i,j,TMA,TMB) * phB
            end do
!$OMP simd
            do i = 1, N
              uM(i,1,1,j,1,iMB) = uM(i,1,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
            end do
          end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,iMB,1)
          if ( TMA*TMB == 1 ) call print_matrix(NA,1,1,NB,1,iMB)
#endif

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            do j = jmin, jmax
!$OMP simd
              do i = 1, N
                uM(i,iMA,iMB,j,1,1) = uM(i,iMA,iMB,j,1,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,iMB,1)
#endif

            ph = phB * conjg(phA)
            do j = jmin, jmax
!$OMP simd
              do i = 1, N
                uM(i,1,iMB,j,iMA,1) = uM(i,1,iMB,j,iMA,1) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,iMB,1)
#endif

            ph = conjg(phB) * phA
            do j = jmin, jmax
!$OMP simd
              do i = 1, N
                uM(i,iMA,1,j,1,iMB) = uM(i,iMA,1,j,1,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,iMA,1,NB,1,iMB)
#endif

            ph = conjg(phB * phA)
            do j = jmin, jmax
!$OMP simd
              do i = 1, N
                uM(i,1,1,j,iMA,iMB) = uM(i,1,1,j,iMA,iMB) + M(i,j,TMA,TMB) * ph
              end do
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( TMA*TMB == 1 ) call print_matrix(NA,1,iMA,NB,1,iMB)
#endif

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    do TMA = 2, NA
      do j = jmin, jmax
!$OMP simd collapse(2)
        do iMA = 2, NA
          do i = 1, N
            uM(i,iMA,1,j,TMA,1) = uM(i,iMA-1,1,j,TMA-1,1)
          end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
          if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,1,1, c='-')
#endif
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do j = jmin, jmax
          do iMA = 2, NA
!$OMP simd
            do i = 1, N
              uM(i,iMA,TMB,j,TMA,1) = uM(i,iMA-1,TMB,j,TMA-1,1)
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,TMB,1, c='-')
#endif

!$OMP simd
            do i = 1, N
              uM(i,iMA,1,j,TMA,TMB) = uM(i,iMA-1,1,j,TMA-1,TMB)
            end do
#ifdef __BLOCH_UNFOLD_M_DEBUG
            if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,1,TMB, c='-')
#endif
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
!$OMP simd collapse(3)
          do iMB = 2, NB
            do iMA = 1, NA
              do i = 1, N
                uM(i,iMA,iMB,j,TMA,TMB) = uM(i,iMA,iMB-1,j,TMA,TMB-1)
              end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
              if ( j == jmin ) call print_matrix(NA,iMA,TMA,NB,iMB,TMB, c='-')
#endif
            end do
          end do
        end do
      end do
    end do

#ifdef __BLOCH_UNFOLD_M_DEBUG
    call print_matrix(NA,0,0,NB,0,0, out=.true.)
#endif

!$OMP end parallel

  end subroutine bloch_unfold_M_manual_simd_2
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
  subroutine bloch_unfold_HS_G_original(this, bk, N, H, S, G, Z, uSZmH, uG)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%B(1),this%B(2),this%B(3)), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%B(1),this%B(2),this%B(3),N,this%B(1),this%B(2),this%B(3)), intent(inout) :: uSZmH, uG

    integer :: B(3), i1, i2, i3
    integer :: ia1, ia2, ia3,iuo
    integer :: ja1, ja2, ja3,juo
    real(dp) :: wq, rPi(3)
    complex(dp) :: qPi, p(3), pZ

!$OMP parallel default(shared) private(wq,rPi,qPi,p,pZ) &
!$OMP&  private(B,i1,i2,i3) &
!$OMP&  private(ia1,ia2,ia3,iuo) &
!$OMP&  private(ja1,ja2,ja3,juo)

!$OMP workshare
    uSZmH(:,:,:,:,:,:,:,:) = 0._dp
    uG(:,:,:,:,:,:,:,:) = 0._dp
!$OMP end workshare

  ! Save some multiplications
    B(:) = this%B(:)
    wq = log(1._dp / real(product(B),dp))

    do i3 = 1, B(3)
      rPi(3) = unfold_k2pi(bk(3), B(3), i3)
      do i2 = 1, B(2)
        rPi(2) = unfold_k2pi(bk(2), B(2), i2)
        do i1 = 1, B(1)
          rPi(1) = unfold_k2pi(bk(1), B(1), i1)

          qPi = exp(cmplx(0._dp,rPi(1), dp))

!$OMP do collapse(3) schedule(static)
          do ia3 = 1 , B(3)
            do ia2 = 1 , B(2)
              do ia1 = 1 , B(1)

                p(3) = exp(cmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3), dp))
                do iuo = 1, N
                  do ja3 = 1 , B(3)
                    p(2) = p(3)*exp(cmplx(0._dp,ja3*rPi(3), dp))
                    do ja2 = 1 , B(2)
                      p(1) = p(2)*exp(cmplx(0._dp,ja2*rPi(2), dp))
                      do ja1 = 1 , B(1)
                        p(1) = p(1) * qPi
                        pZ = p(1) * Z
                        do juo = 1, N
                          uSZmH(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) = uSZmH(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) &
                              + pZ * S(juo,iuo,i1,i2,i3) - p(1) * H(juo,iuo,i1,i2,i3)
                          uG(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) = uG(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) &
                              + p(1) * G(juo,iuo,i1,i2,i3)
                        end do !juo
                      end do !ja1
                    end do !ja2
                  end do !ja3
                end do !iuo
              end do !ia1
            end do !ia2
          end do !ia3
!$OMP end do

        end do !i1
      end do !i2
    end do !i3

!$OMP end parallel

  end subroutine bloch_unfold_HS_G_original
#endif

#if defined(__BLOCH_UNFOLD_M_TEST) || ! defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_HS_G_do(this, bk, N, H, S, G, Z, uSZmH, uG)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH, uG

    if ( this%prod_B == 1 ) then

      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)
      call zcopy(N*N, G(1,1,1), 1, uG(1,1,1,1), 1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_HS_G_do_1(N, H, S, G, Z, bk(3), this%B(3), uSZmH, uG)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_G_do_1(N, H, S, G, Z, bk(2), this%B(2), uSZmH, uG)
      else
        call bloch_unfold_HS_G_do_2(N, H, S, G, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_G_do_1(N, H, S, G, Z, bk(1), this%B(1), uSZmH, uG)
      else
        call bloch_unfold_HS_G_do_2(N, H, S, G, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_HS_G_do_2(N, H, S, G, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH, uG)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_HS_G_do

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_G_do_1(N, H, S, G, Z, kA, NA, uSZmH, uG)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), dimension(N,NA,N,NA), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: phA_step
    complex(dp) :: ph, cph, phZ, cphZ

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,cph,phZ,cphZ) firstprivate(Z)

    ! Initialize un-folded matrix
!$OMP do schedule(static)
    do j = 1, N
      uSZmH(:,:,j,1) = cmplx(0._dp, 0._dp, dp)
      uG(:,:,j,1) = cmplx(0._dp, 0._dp, dp)
      uSZmH(:,1,j,2:) = cmplx(0._dp, 0._dp, dp)
      uG(:,1,j,2:) = cmplx(0._dp, 0._dp, dp)
    end do
!$OMP end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
      phZ = Z * w
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * w
        end do
        do i = 1, N
          uG(i,1,j,1) = uG(i,1,j,1) + G(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        cph = conjg(ph)
        phZ = ph * Z
        cphZ = cph * Z
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + S(i,j,TMA) * cphZ - H(i,j,TMA) * cph
          end do
          do i = 1, N
            uG(i,iMA,j,1) = uG(i,iMA,j,1) + G(i,j,TMA) * ph
          end do
          do i = 1, N
            uG(i,1,j,iMA) = uG(i,1,j,iMA) + G(i,j,TMA) * cph
          end do
        end do
!$OMP end do nowait
        ph = ph * phA_step
      end do

!$OMP barrier

    end do

    do iMA = 2, NA
      do TMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,TMA,j,iMA) = uSZmH(i,TMA-1,j,iMA-1)
          end do
          do i = 1, N
            uG(i,TMA,j,iMA) = uG(i,TMA-1,j,iMA-1)
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_G_do_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_G_do_2(N, H, S, G, Z, kA, NA, kB, NB, uSZmH, uG)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), dimension(N,NA,NB,N,NA,NB), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, cph, phZ, cphZ

    uSZmH(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)
    uG(:,:,:,:,:,:) = cmplx(0._dp, 0._dp, dp)

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,cph,phZ,cphZ) firstprivate(Z)

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
        phZ = Z * w
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
          do i = 1, N
            uG(i,1,1,j,1,1) = uG(i,1,1,j,1,1) + G(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA
          ph = phA
          cph = conjg(phA)
          phZ = ph * Z
          cphZ = cph * Z
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do

            do i = 1, N
              uG(i,iMA,1,j,1,1) = uG(i,iMA,1,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,iMA,1) = uG(i,1,1,j,iMA,1) + G(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait
          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          cph = conjg(phB)
          phZ = ph * Z
          cphZ = cph * Z
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do

            do i = 1, N
              uG(i,1,iMB,j,1,1) = uG(i,1,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,1,iMB) = uG(i,1,1,j,1,iMB) + G(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA
            ph = phB * phA
            cph = phB * conjg(phA)
            phZ = ph * Z
            cphZ = cph * Z
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do

              do i = 1, N
                uG(i,iMA,iMB,j,1,1) = uG(i,iMA,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,iMB,j,iMA,1) = uG(i,1,iMB,j,iMA,1) + G(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait

            ph = conjg(phB) * phA
            cph = conjg(phB) * conjg(phA)
            phZ = ph * Z
            cphZ = cph * Z
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do

              do i = 1, N
                uG(i,iMA,1,j,1,iMB) = uG(i,iMA,1,j,1,iMB) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,1,j,iMA,iMB) = uG(i,1,1,j,iMA,iMB) + G(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

!$OMP barrier

      end do
    end do

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do TMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,TMA,1,j,iMA,1) = uSZmH(i,TMA-1,1,j,iMA-1,1)
          end do
          do i = 1, N
            uG(i,TMA,1,j,iMA,1) = uG(i,TMA-1,1,j,iMA-1,1)
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier
    end do

    do iMB = 2, NB
      do iMA = 2, NA
        do TMA = 2, NA
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,TMA,iMB,j,iMA,1) = uSZmH(i,TMA-1,iMB,j,iMA-1,1)
            end do
            do i = 1, N
              uG(i,TMA,iMB,j,iMA,1) = uG(i,TMA-1,iMB,j,iMA-1,1)
            end do
          end do
!$OMP end do nowait
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,TMA,1,j,iMA,iMB) = uSZmH(i,TMA-1,1,j,iMA-1,iMB)
            end do
            do i = 1, N
              uG(i,TMA,1,j,iMA,iMB) = uG(i,TMA-1,1,j,iMA-1,iMB)
            end do
          end do
!$OMP end do nowait
        end do
      end do
!$OMP barrier
    end do

    ! Now we have filled all
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! for NB blocks
    do iMB = 2, NB
      do iMA = 1, NA
        do TMB = 2, NB
          do TMA = 1, NA
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,TMA,TMB,j,iMA,iMB) = uSZmH(i,TMA,TMB-1,j,iMA,iMB-1)
              end do
              do i = 1, N
                uG(i,TMA,TMB,j,iMA,iMB) = uG(i,TMA,TMB-1,j,iMA,iMB-1)
              end do
            end do
!$OMP end do nowait
          end do
        end do
      end do
!$OMP barrier
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_G_do_2
#endif


#if defined(__BLOCH_UNFOLD_M_TEST) || defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_HS_G_manual(this, bk, N, H, S, G, Z, uSZmH, uG)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH, uG

    if ( this%prod_B == 1 ) then

      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)
      call zcopy(N*N, G(1,1,1), 1, uG(1,1,1,1), 1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_HS_G_manual_1(N, H, S, G, Z, bk(3), this%B(3), uSZmH, uG)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_G_manual_1(N, H, S, G, Z, bk(2), this%B(2), uSZmH, uG)
      else
        call bloch_unfold_HS_G_manual_2(N, H, S, G, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_G_manual_1(N, H, S, G, Z, bk(1), this%B(1), uSZmH, uG)
      else
        call bloch_unfold_HS_G_manual_2(N, H, S, G, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH, uG)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_HS_G_manual_2(N, H, S, G, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH, uG)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_HS_G_manual

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_G_manual_1(N, H, S, G, Z, kA, NA, uSZmH, uG)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), dimension(N,NA,N,NA), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phZ, cphZ, phA_step

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph,cph,phZ,cphZ) &
!$OMP&   private(k_A,phA_step,TMA,iMA) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uSZmH(i,TMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uG(i,TMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,j,TMA) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uG(i,1,j,TMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      phZ = Z * w
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + phZ * S(i,j,TMA) - H(i,j,TMA) * w
        end do
        do i = 1, N
          uG(i,1,j,1) = uG(i,1,j,1) + G(i,j,TMA) * w
        end do
      end do

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

        cph = conjg(ph)
        phZ = Z * ph
        cphZ = Z * cph
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + phZ * S(i,j,TMA) - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + cphZ * S(i,j,TMA) - H(i,j,TMA) * cph
          end do
        end do
        do j = jmin, jmax
          do i = 1, N
            uG(i,iMA,j,1) = uG(i,iMA,j,1) + ph * G(i,j,TMA)
          end do
          do i = 1, N
            uG(i,1,j,iMA) = uG(i,1,j,iMA) + cph * G(i,j,TMA)
          end do
        end do

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uSZmH(:,:,:,1)
    !   uSZmH(:,1,:,:)
    ! Due to uSZmH being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do j = jmin, jmax
        do TMA = 2, NA
          do i = 1, N
            uSZmH(i,TMA,j,iMA) = uSZmH(i,TMA-1,j,iMA-1)
          end do
          do i = 1, N
            uG(i,TMA,j,iMA) = uG(i,TMA-1,j,iMA-1)
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_G_manual_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_G_manual_2(N, H, S, G, Z, kA, NA, kB, NB, uSZmH, uG)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S, G
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), dimension(N,NA,NB,N,NA,NB), intent(inout) :: uSZmH, uG

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, cph, phZ, cphZ

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph,cph,phZ,cphZ) &
!$OMP&   private(k_A,phA_step,phA,TMA,iMA) &
!$OMP&   private(k_B,phB_step,phB,TMB,iMB) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uSZmH(i,TMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uG(i,TMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,1,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uG(i,1,1,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    do iMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,TMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,1,j,TMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uG(i,TMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uG(i,1,1,j,TMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
      end do

      do TMA = 1, NA
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,TMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,iMB,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uG(i,TMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uG(i,1,iMB,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
      end do
    end do

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        phZ = Z * w
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
          do i = 1, N
            uG(i,1,1,j,1,1) = uG(i,1,1,j,1,1) + G(i,j,TMA,TMB) * w
          end do
        end do

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

          ph = phA
          cph = conjg(ph)
          phZ = Z * ph
          cphZ = Z * cph

          do j = jmin, jmax
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
            do i = 1, N
              uG(i,iMA,1,j,1,1) = uG(i,iMA,1,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,iMA,1) = uG(i,1,1,j,iMA,1) + G(i,j,TMA,TMB) * cph
            end do
          end do

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          cph = conjg(phB)
          phZ = Z * ph
          cphZ = Z * cph

          do j = jmin, jmax
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
            do i = 1, N
              uG(i,1,iMB,j,1,1) = uG(i,1,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uG(i,1,1,j,1,iMB) = uG(i,1,1,j,1,iMB) + G(i,j,TMA,TMB) * cph
            end do
          end do

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            cph = phB * conjg(phA)
            phZ = Z * ph
            cphZ = Z * cph
            do j = jmin, jmax
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
              do i = 1, N
                uG(i,iMA,iMB,j,1,1) = uG(i,iMA,iMB,j,1,1) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,iMB,j,iMA,1) = uG(i,1,iMB,j,iMA,1) + G(i,j,TMA,TMB) * cph
              end do
            end do

            ph = conjg(phB) * phA
            cph = conjg(phB * phA)
            phZ = Z * ph
            cphZ = Z * cph
            do j = jmin, jmax
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
              do i = 1, N
                uG(i,iMA,1,j,1,iMB) = uG(i,iMA,1,j,1,iMB) + G(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uG(i,1,1,j,iMA,iMB) = uG(i,1,1,j,iMA,iMB) + G(i,j,TMA,TMB) * cph
              end do
            end do

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    do TMA = 2, NA
      do j = jmin, jmax
        do iMA = 2, NA
          do i = 1, N
            uSZmH(i,iMA,1,j,TMA,1) = uSZmH(i,iMA-1,1,j,TMA-1,1)
          end do
          do i = 1, N
            uG(i,iMA,1,j,TMA,1) = uG(i,iMA-1,1,j,TMA-1,1)
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do j = jmin, jmax
          do iMA = 2, NA
            do i = 1, N
              uSZmH(i,iMA,TMB,j,TMA,1) = uSZmH(i,iMA-1,TMB,j,TMA-1,1)
            end do
            do i = 1, N
              uSZmH(i,iMA,1,j,TMA,TMB) = uSZmH(i,iMA-1,1,j,TMA-1,TMB)
            end do
            do i = 1, N
              uG(i,iMA,TMB,j,TMA,1) = uG(i,iMA-1,TMB,j,TMA-1,1)
            end do
            do i = 1, N
              uG(i,iMA,1,j,TMA,TMB) = uG(i,iMA-1,1,j,TMA-1,TMB)
            end do
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do iMB = 2, NB
            do iMA = 1, NA
              do i = 1, N
                uSZmH(i,iMA,iMB,j,TMA,TMB) = uSZmH(i,iMA,iMB-1,j,TMA,TMB-1)
              end do
              do i = 1, N
                uG(i,iMA,iMB,j,TMA,TMB) = uG(i,iMA,iMB-1,j,TMA,TMB-1)
              end do
            end do
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_G_manual_2
#endif

#ifdef __BLOCH_UNFOLD_M_TEST
  subroutine bloch_unfold_HS_original(this, bk, N, H, S, Z, uSZmH)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%B(1),this%B(2),this%B(3)), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%B(1),this%B(2),this%B(3),N,this%B(1),this%B(2),this%B(3)), intent(inout) :: uSZmH

    integer :: B(3), i1, i2, i3
    integer :: ia1, ia2, ia3,iuo
    integer :: ja1, ja2, ja3,juo
    real(dp) :: wq, rPi(3)
    complex(dp) :: qPi, p(3), pZ

!$OMP parallel default(shared) private(wq,rPi,qPi,p,pZ) &
!$OMP&  private(B,i1,i2,i3) &
!$OMP&  private(ia1,ia2,ia3,iuo) &
!$OMP&  private(ja1,ja2,ja3,juo)

!$OMP workshare
    uSZmH(:,:,:,:,:,:,:,:) = 0._dp
!$OMP end workshare

  ! Save some multiplications
    B(:) = this%B(:)
    wq = log(1._dp / real(product(B),dp))

    do i3 = 1, B(3)
      rPi(3) = unfold_k2pi(bk(3), B(3), i3)
      do i2 = 1, B(2)
        rPi(2) = unfold_k2pi(bk(2), B(2), i2)
        do i1 = 1, B(1)
          rPi(1) = unfold_k2pi(bk(1), B(1), i1)

          qPi = exp(cmplx(0._dp,rPi(1), dp))

!$OMP do collapse(3) schedule(static)
          do ia3 = 1 , B(3)
            do ia2 = 1 , B(2)
              do ia1 = 1 , B(1)

                p(3) = exp(cmplx(wq,-ia1*rPi(1)-ia2*rPi(2)-ia3*rPi(3), dp))
                do iuo = 1, N
                  do ja3 = 1 , B(3)
                    p(2) = p(3)*exp(cmplx(0._dp,ja3*rPi(3), dp))
                    do ja2 = 1 , B(2)
                      p(1) = p(2)*exp(cmplx(0._dp,ja2*rPi(2), dp))
                      do ja1 = 1 , B(1)
                        p(1) = p(1)*qPi
                        pZ = p(1) * Z
                        do juo = 1, N
                          uSZmH(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) = uSZmH(juo,ja1,ja2,ja3,iuo,ia1,ia2,ia3) + &
                              pZ * S(juo,iuo,i1,i2,i3) - p(1) * H(juo,iuo,i1,i2,i3)
                        end do !juo
                      end do !ja1
                    end do !ja2
                  end do !ja3
                end do !iuo
              end do !ia1
            end do !ia2
          end do !ia3
!$OMP end do

        end do !i1
      end do !i2
    end do !i3

!$OMP end parallel

  end subroutine bloch_unfold_HS_original
#endif

#if defined(__BLOCH_UNFOLD_M_TEST) || ! defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_HS_do(this, bk, N, H, S, Z, uSZmH)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH

    if ( this%prod_B == 1 ) then

      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_HS_do_1(N, H, S, Z, bk(3), this%B(3), uSZmH)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_do_1(N, H, S, Z, bk(2), this%B(2), uSZmH)
      else
        call bloch_unfold_HS_do_2(N, H, S, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_do_1(N, H, S, Z, bk(1), this%B(1), uSZmH)
      else
        call bloch_unfold_HS_do_2(N, H, S, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_HS_do_2(N, H, S, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_HS_do

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_do_1(N, H, S, Z, kA, NA, uSZmH)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uSZmH(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    real(dp) :: w, k_A
    complex(dp) :: phA_step
    complex(dp) :: ph, cph, phZ, cphZ

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step) &
!$OMP& private(w,ph,cph,phZ,cphZ) firstprivate(Z)

    ! Initialize un-folded matrix
!$OMP do schedule(static)
    do j = 1, N
      do iMA = 1, NA
        do i = 1, N
          uSZmH(i,iMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
!$OMP end do nowait
!$OMP do schedule(static)
    do j = 1, N
      do iMA = 2, NA
        do i = 1, N
          uSZmH(i,1,j,iMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
!$OMP end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      ! We don't collapse and hope the compiler will use vector
      ! instructions for the inner loop!
      phZ = Z * w
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * w
        end do
      end do
!$OMP end do nowait

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA
        cph = conjg(ph)
        phZ = Z * ph
        cphZ = Z * cph
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + S(i,j,TMA) * phZ - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + S(i,j,TMA) * cphZ - H(i,j,TMA) * cph
          end do
        end do
!$OMP end do nowait
        ph = ph * phA_step
      end do

!$OMP barrier

    end do

    do iMA = 2, NA
      do TMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,TMA,j,iMA) = uSZmH(i,TMA-1,j,iMA-1)
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_do_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_do_2(N, H, S, Z, kA, NA, kB, NB, uSZmH)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uSZmH(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, cph, phZ, cphZ

!$OMP parallel default(shared) private(i,j) &
!$OMP& private(TMA,iMA,k_A,phA_step,phA) &
!$OMP& private(TMB,iMB,k_B,phB_step,phB) &
!$OMP& private(w,ph,cph,phZ,cphZ) firstprivate(Z)

    ! Initialize un-folded matrix
!$OMP do schedule(static)
    do j = 1, N
      do i = 1, N
        uSZmH(i,1,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
      end do
    end do
!$OMP end do nowait

    do iMA = 2, NA
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uSZmH(i,iMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uSZmH(i,1,1,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end do nowait
    end do

    do iMB = 2, NB
!$OMP do schedule(static)
      do j = 1, N
        do i = 1, N
          uSZmH(i,1,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
        do i = 1, N
          uSZmH(i,1,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
!$OMP end do nowait

      do iMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,iMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,iMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,iMB,j,iMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,1,j,iMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
!$OMP end do nowait
      end do
    end do

!$OMP barrier

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        ! We don't collapse and hope the compiler will use vector
        ! instructions for the inner loop!
        phZ = Z * w
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
        end do
!$OMP end do nowait

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA
          ph = phA
          cph = conjg(phA)
          phZ = Z * ph
          cphZ = Z * cph
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait
          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          cph = conjg(phB)
          phZ = Z * ph
          cphZ = Z * cph
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
          end do
!$OMP end do nowait

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            cph = phB * conjg(phA)
            phZ = Z * ph
            cphZ = Z * cph
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait

            ph = conjg(phB) * phA
            cph = conjg(phB * phA)
            phZ = Z * ph
            cphZ = Z * cph
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
            end do
!$OMP end do nowait

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

!$OMP barrier

      end do
    end do

    ! Due to uM being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do TMA = 2, NA
!$OMP do schedule(static)
        do j = 1, N
          do i = 1, N
            uSZmH(i,TMA,1,j,iMA,1) = uSZmH(i,TMA-1,1,j,iMA-1,1)
          end do
        end do
!$OMP end do nowait
      end do
!$OMP barrier
    end do

    do iMB = 2, NB
      do iMA = 2, NA
        do TMA = 2, NA
!$OMP do schedule(static)
          do j = 1, N
            do i = 1, N
              uSZmH(i,TMA,iMB,j,iMA,1) = uSZmH(i,TMA-1,iMB,j,iMA-1,1)
            end do
            do i = 1, N
              uSZmH(i,TMA,1,j,iMA,iMB) = uSZmH(i,TMA-1,1,j,iMA-1,iMB)
            end do
          end do
!$OMP end do nowait
        end do
      end do
    end do

!$OMP barrier

    ! Now we have filled all
    !   uM(:,:,:,1)
    !   uM(:,1,:,:)
    ! for NB blocks
    do iMB = 2, NB
      do iMA = 1, NA
        do TMB = 2, NB
          do TMA = 1, NA
!$OMP do schedule(static)
            do j = 1, N
              do i = 1, N
                uSZmH(i,TMA,TMB,j,iMA,iMB) = uSZmH(i,TMA,TMB-1,j,iMA,iMB-1)
              end do
            end do
!$OMP end do nowait
          end do
        end do
      end do
!$OMP barrier
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_do_2
#endif


#if defined(__BLOCH_UNFOLD_M_TEST) || defined(_OPENMP)
  !< Unfold a specific index for the matrix unfold machinery from M -> uM
  subroutine bloch_unfold_HS_manual(this, bk, N, H, S, Z, uSZmH)
    class(bloch_unfold_t), intent(in) :: this
    real(dp), intent(in) :: bk(3)
    integer, intent(in) :: N
    complex(dp), dimension(N,N,this%prod_B), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    complex(dp), dimension(N,this%prod_B,N,this%prod_B), intent(inout) :: uSZmH

    if ( this%prod_B == 1 ) then

      uSZmH(:,1,:,1) = Z * S(:,:,1) - H(:,:,1)

      return
    end if

    ! Now do the actual expansion
    if ( this%B(1) == 1 ) then
      if ( this%B(2) == 1 ) then
        call bloch_unfold_HS_manual_1(N, H, S, Z, bk(3), this%B(3), uSZmH)
      else if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_manual_1(N, H, S, Z, bk(2), this%B(2), uSZmH)
      else
        call bloch_unfold_HS_manual_2(N, H, S, Z, &
            bk(2), this%B(2), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(2) == 1 ) then
      if ( this%B(3) == 1 ) then
        call bloch_unfold_HS_manual_1(N, H, S, Z, bk(1), this%B(1), uSZmH)
      else
        call bloch_unfold_HS_manual_2(N, H, S, Z, &
            bk(1), this%B(1), &
            bk(3), this%B(3), uSZmH)
      end if
    else if ( this%B(3) == 1 ) then
      call bloch_unfold_HS_manual_2(N, H, S, Z, &
          bk(1), this%B(1), &
          bk(2), this%B(2), uSZmH)
    else
      call die('currently not implemented')
    end if

  end subroutine bloch_unfold_HS_manual

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_manual_1(N, H, S, Z, kA, NA, uSZmH)
    integer, intent(in) :: N, NA
    complex(dp), dimension(N,N,NA), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA
    complex(dp), intent(inout) :: uSZmH(N,NA,N,NA)

    integer :: TMA, iMA
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A
    complex(dp) :: ph, cph, phZ, cphZ, phA_step

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph,cph,phZ,cphZ) &
!$OMP&   private(k_A,phA_step,TMA,iMA) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uSZmH(i,TMA,j,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,j,TMA) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    w = 1._dp / real(NA, dp)

    ! TMA loop is over different matrices
    do TMA = 1, NA

      ! The first diagonal part is *always* phase-less
      ! So there is no reason to multiply with *extra* stuff
      ! (1, :, 1, :) = sum(M(1, 1, :))

      k_A = unfold_k2pi(kA, NA, TMA)

      ! exp(2\pi k) where k == K/N
      phA_step = cmplx(cos(k_A), sin(k_A), dp)

      phZ = Z * w
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,j,1) = uSZmH(i,1,j,1) + phZ * S(i,j,TMA) - H(i,j,TMA) * w
        end do
      end do

      ! perform average in this phase-factor
      ph = w * phA_step
      do iMA = 2, NA

        cph = conjg(ph)
        phZ = Z * ph
        cphZ = Z * cph
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,iMA,j,1) = uSZmH(i,iMA,j,1) + phZ * S(i,j,TMA) - H(i,j,TMA) * ph
          end do
          do i = 1, N
            uSZmH(i,1,j,iMA) = uSZmH(i,1,j,iMA) + cphZ * S(i,j,TMA) - H(i,j,TMA) * cph
          end do
        end do

        ph = ph * phA_step
      end do

    end do

    ! At this point the following have been calculated:
    !   uSZmH(:,:,:,1)
    !   uSZmH(:,1,:,:)
    ! Due to uSZmH being a Toeplitz matrix, we simply copy around the data!
    do iMA = 2, NA
      do j = jmin, jmax
        do TMA = 2, NA
          do i = 1, N
            uSZmH(i,TMA,j,iMA) = uSZmH(i,TMA-1,j,iMA-1)
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_manual_1

  !< Unfold a given matrix for only one un-fold point
  subroutine bloch_unfold_HS_manual_2(N, H, S, Z, kA, NA, kB, NB, uSZmH)
    integer, intent(in) :: N, NA, NB
    complex(dp), dimension(N,N,NA,NB), intent(in) :: H, S
    complex(dp), intent(in) :: Z
    real(dp), intent(in) :: kA, kB
    complex(dp), intent(inout) :: uSZmH(N,NA,NB,N,NA,NB)

    integer :: TMA, iMA, TMB, iMB
    integer :: i, j
    integer :: NT, IT
    integer :: jmin, jmax
    real(dp) :: w, k_A, k_B
    complex(dp) :: phA, phA_step
    complex(dp) :: phB, phB_step
    complex(dp) :: ph, cph, phZ, cphZ

!$OMP parallel default(shared) &
!$OMP&   private(i,j,w,ph,cph,phZ,cphZ) &
!$OMP&   private(k_A,phA_step,phA,TMA,iMA) &
!$OMP&   private(k_B,phB_step,phB,TMB,iMB) &
!$OMP&   private(NT,IT,jmin,jmax)

    call bloch_unfold_openmp_init(nt, it)
    call bloch_unfold_openmp_loop_split(nt, it, N, jmin, jmax)

    ! Initialize
    do j = jmin, jmax
      do TMA = 1, NA
        do i = 1, N
          uSZmH(i,TMA,1,j,1,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do
    do TMA = 2 , NA
      do j = jmin, jmax
        do i = 1, N
          uSZmH(i,1,1,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
        end do
      end do
    end do

    do iMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,TMA,iMB,j,1,1) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,1,j,TMA,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
      end do

      do TMA = 1, NA
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,TMA,1,j,1,iMB) = cmplx(0._dp, 0._dp, dp)
          end do
          do i = 1, N
            uSZmH(i,1,iMB,j,TMA,1) = cmplx(0._dp, 0._dp, dp)
          end do
        end do
      end do
    end do

    ! Full weight
    w = 1._dp / real(NA*NB, dp)

    ! TMB/A loop is over different matrices
    do TMB = 1, NB

      k_B = unfold_k2pi(kB, NB, TMB)
      phB_step = cmplx(cos(k_B), sin(k_B), dp)

      do TMA = 1, NA

        ! The first diagonal part is *always* phase-less
        ! So there is no reason to multiply with *extra* stuff
        ! (1, :, TMB, 1, :, TMB) = sum(M(1, 1, :, TMB))

        k_A = unfold_k2pi(kA, NA, TMA)
        phA_step = cmplx(cos(k_A), sin(k_A), dp)

        phZ = Z * w
        do j = jmin, jmax
          do i = 1, N
            uSZmH(i,1,1,j,1,1) = uSZmH(i,1,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * w
          end do
        end do

        ! The Toeplitz nature of a double Bloch expanded matrix
        ! makes this a bit more challenging
        ! Within each sub-block where NA is expanded we find
        ! the same structure as in the single Bloch expanded matrix

        ! for iMB == 1 and iMA in [1...NA] we need phB == w
        phA = w * phA_step

        ! This is the first iMB == 1 sub-block
        do iMA = 2, NA

          ph = phA
          cph = conjg(ph)
          phZ = Z * ph
          cphZ = Z * cph

          do j = jmin, jmax
            do i = 1, N
              uSZmH(i,iMA,1,j,1,1) = uSZmH(i,iMA,1,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,iMA,1) = uSZmH(i,1,1,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
          end do

          phA = phA * phA_step
        end do

        ! Now handle all iMB > 1
        phB = w * phB_step
        do iMB = 2, NB

          ! Diagonal B (A has zero phase)
          ph = phB
          cph = conjg(phB)
          phZ = Z * ph
          cphZ = Z * cph

          do j = jmin, jmax
            do i = 1, N
              uSZmH(i,1,iMB,j,1,1) = uSZmH(i,1,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
            end do
            do i = 1, N
              uSZmH(i,1,1,j,1,iMB) = uSZmH(i,1,1,j,1,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
            end do
          end do

          ! Now do all A off-diagonals (phB has weight)
          phA = phA_step
          do iMA = 2, NA

            ph = phB * phA
            cph = phB * conjg(phA)
            phZ = Z * ph
            cphZ = Z * cph
            do j = jmin, jmax
              do i = 1, N
                uSZmH(i,iMA,iMB,j,1,1) = uSZmH(i,iMA,iMB,j,1,1) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,iMB,j,iMA,1) = uSZmH(i,1,iMB,j,iMA,1) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
            end do

            ph = conjg(phB) * phA
            cph = conjg(phB * phA)
            phZ = Z * ph
            cphZ = Z * cph
            do j = jmin, jmax
              do i = 1, N
                uSZmH(i,iMA,1,j,1,iMB) = uSZmH(i,iMA,1,j,1,iMB) + S(i,j,TMA,TMB) * phZ - H(i,j,TMA,TMB) * ph
              end do
              do i = 1, N
                uSZmH(i,1,1,j,iMA,iMB) = uSZmH(i,1,1,j,iMA,iMB) + S(i,j,TMA,TMB) * cphZ - H(i,j,TMA,TMB) * cph
              end do
            end do

            phA = phA * phA_step
          end do

          phB = phB * phB_step
        end do

      end do
    end do

    do TMA = 2, NA
      do j = jmin, jmax
        do iMA = 2, NA
          do i = 1, N
            uSZmH(i,iMA,1,j,TMA,1) = uSZmH(i,iMA-1,1,j,TMA-1,1)
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 2, NA
        do j = jmin, jmax
          do iMA = 2, NA
            do i = 1, N
              uSZmH(i,iMA,TMB,j,TMA,1) = uSZmH(i,iMA-1,TMB,j,TMA-1,1)
            end do
            do i = 1, N
              uSZmH(i,iMA,1,j,TMA,TMB) = uSZmH(i,iMA-1,1,j,TMA-1,TMB)
            end do
          end do
        end do
      end do
    end do

    do TMB = 2, NB
      do TMA = 1, NA
        do j = jmin, jmax
          do iMB = 2, NB
            do iMA = 1, NA
              do i = 1, N
                uSZmH(i,iMA,iMB,j,TMA,TMB) = uSZmH(i,iMA,iMB-1,j,TMA,TMB-1)
              end do
            end do
          end do
        end do
      end do
    end do

!$OMP end parallel

  end subroutine bloch_unfold_HS_manual_2
#endif

  pure function bloch_unfold_get_k(this, iA, iB, iC) result(bk)
    class(bloch_unfold_t), intent(in) :: this
    integer, intent(in) :: iA, iB, iC
    real(dp) :: bk(3)

    ! TODO, the current implementation assumes k-symmetry
    ! of the electrode electronic structure.
    ! Using Bloch expansion with non-symmetry will, likely, produce
    ! wrong results.
    ! Luckily this is not a problem currently.
    ! Perhaps one should consider this in tbtrans

    bk(1) = real(iA - 1, dp) / real(this%B(1), dp)
    bk(2) = real(iB - 1, dp) / real(this%B(2), dp)
    bk(3) = real(iC - 1, dp) / real(this%B(3), dp)

  end function bloch_unfold_get_k

  pure function unfold_k2pi(bk, Nk, ik) result(k)
#ifdef __BLOCH_UNFOLD_M_TEST
    real(dp), parameter :: pi2 = 6.28318530717958623199592693708837032318115234375_dp
#else
    use units, only: Pi
    real(dp), parameter :: pi2 = Pi * 2._dp
#endif
    real(dp), intent(in) :: bk
    integer, intent(in) :: Nk, ik
    real(dp) :: k
    ! bk should already be in *local* coordinates
    k = Pi2 * (bk + real(ik - 1, dp) / real(Nk, dp))
  end function unfold_k2pi

  subroutine bloch_unfold_openmp_init(nt, it)
!$  use omp_lib
    !< Number of threads
    integer, intent(out) :: nt
    !< ID of thread
    integer, intent(out) :: it

#ifdef _OPENMP
    nt = omp_get_num_threads()
    it = omp_get_thread_num()
#else
    nt = 1
    it = 0
#endif

  end subroutine bloch_unfold_openmp_init

  pure subroutine bloch_unfold_openmp_loop_split(nt, it, N, Nmin, Nmax)
    !< Number of threads in the construct
    integer, intent(in) :: nt
    !< Thread ID in the construct
    integer, intent(in) :: it
    !< Number of elements to split
    integer, intent(in) :: N
    !< Minimum/max element to use
    integer, intent(out) :: Nmin, Nmax

#ifdef _OPENMP
    integer :: i

    Nmin = N / nt
    Nmax = (it + 1) * Nmin
    Nmin = Nmax - Nmin + 1
    i = mod(N, nt)
    if ( i > 0 ) then
      if ( it < i ) then
        Nmin = Nmin + it
        Nmax = Nmax + it + 1
      else
        Nmin = Nmin + i
        Nmax = Nmax + i
      end if
    end if
    if ( N < nt ) then
      Nmin = it + 1
      Nmax = it + 1
      if ( it >= N ) Nmax = N
    end if
#else
    Nmin = 1
    Nmax = N
#endif

  end subroutine bloch_unfold_openmp_loop_split

  !< Unravel a linear bloch-expansion coefficient to its Bloch-indices
  pure subroutine bloch_unfold_unravel_index(this, lin, iB)
    class(bloch_unfold_t), intent(in) :: this
    integer, intent(in) :: lin
    integer, intent(out) :: iB(3)

    integer :: i, j, k

    i = this%B(1)
    j = i * this%B(2)
    k = j * this%B(3)
    if ( lin <= i ) then
      iB(3) = 1
      iB(2) = 1
      iB(1) = lin
    else if ( lin <= j ) then
      j = lin / i
      if ( mod(lin, i) /= 0 ) j = j + 1
      iB(3) = 1
      iB(2) = j
      iB(1) = lin - (j-1) * i
    else if ( lin <= k ) then
      k = lin / j
      if ( mod(lin, j) /= 0 ) k = k + 1
      iB(3) = k
      k = lin - (k-1) * j
      j = k / i
      if ( mod(k, i) /= 0 ) j = j + 1
      iB(2) = j
      iB(1) = k - (j-1) * i
    else
      ! signal an error
      iB(:) = 0
    end if

  end subroutine bloch_unfold_unravel_index

#ifdef __BLOCH_UNFOLD_M_DEBUG
  subroutine print_matrix(B1,i1,j1,B2,i2,j2, init, out, column, c)
!$  use omp_lib
    integer, intent(in) :: B1, i1, j1
    integer, intent(in), optional :: B2, i2, j2
    character(len=*), optional :: init
    logical, intent(in), optional :: out
    integer, intent(in), optional :: column
    character(len=1), intent(in), optional :: c

    character(len=1) :: lc
    character(len=64) :: fmt_v, fmt

    logical :: linit, lout
    integer :: lcolumn

    ! Quick return if not printing thread
!$ if ( omp_get_thread_num() > 0 ) return

    lc = '*'
    if ( present(c) ) lc = c
    linit = present(init)
    lout = .false.
    if ( present(out) ) lout = out

    lcolumn = 1
    if ( present(column) ) lcolumn = column
    write(fmt, '("(t",i0,",a)")') lcolumn

    ! Clean
    if ( present(B2) ) then
      call print_2D()
    else
      call print_1D()
    end if

  contains

    subroutine print_1d()
      character(len=1), save, allocatable :: stars(:,:)
      integer :: i

      ! A 1D structure
      if ( linit ) then
        write(*,'(/3a,tr1,i0)') "Debugging ", init, "_1 = ", B1
        if ( allocated(stars) ) deallocate(stars)
        allocate(stars(B1,B1))
        stars = ' '
      else if ( lout ) then
        ! Print out matrix
        write(*,fmt) repeat('-', B1 + 2)
        do i = 1, B1
          write(*,'(1000a)') '|', stars(i,:), '|'
        end do
        write(*,fmt) repeat('-', B1 + 2)
      else
        if ( stars(i1,j1) /= ' ' .and. &
            (stars(i1,j1) /= '0' .and. lc /= '0') ) then
          write(*,'(2a," -> ",a,4(tr1,i0))') "DOUBLE WRITING ",stars(i1,j1),lc,i1,j1
        end if
        stars(i1,j1) = lc
      end if

    end subroutine print_1d

    subroutine print_2d()
      character(len=1), save, allocatable :: stars(:,:,:,:)
      integer :: i, j

      ! A 2D structure
      if ( linit ) then
        write(*,'(/3a,2(tr1,i0))') "Debugging ", init, "_2 = ", B1, B2
        if ( allocated(stars) ) deallocate(stars)
        allocate(stars(B1,B2,B1,B2))
        ! Ensure everything has nothing
        stars = ' '

      else if ( lout ) then

        ! Clarify format
        !"|",B2("|",B1A,"|"),"|"
        write(fmt_v,'("(t",i0,",''|'',",i0,"(''|'',",i0,"a),''|'',''|'')")') lcolumn, B2, B1
        !"|",B2*B1("|",B1A,"|"),"|"
        !write(fmt_h,'("('''',",i0,"(''-'',",i0,"a),''|'',''|'')")') B2, B1

        ! Print out matrix
        write(*,fmt) repeat('-', B2 * (B1 + 1) + 3)
        do j = 1, B2
          do i = 1, B1
            write(*,fmt_v) stars(i,j,:,:)
          end do
          write(*,fmt) repeat('-', B2 * (B1 + 1) + 3)
        end do

      else
        if ( stars(i1,i2,j1,j2) /= ' ' .and. &
            (stars(i1,i2,j1,j2) /= '0' .and. lc /= '0') ) then
          write(*,'(2a," -> ",a,4(tr1,i0))') "DOUBLE WRITING ",stars(i1,i2,j1,j2),lc,i1,j1
        end if
        stars(i1,i2,j1,j2) = lc

      end if

    end subroutine print_2d


  end subroutine print_matrix
#endif

end module bloch_unfold_m


#ifdef __BLOCH_UNFOLD_M_TEST
program bloch_unfold
  use bloch_unfold_m
!$  use omp_lib

  implicit none

  integer, parameter :: dp = selected_real_kind(p=15)

  type(bloch_unfold_t) :: bu

  integer :: N_min = 8
  integer :: N_max = 128
  integer :: N_step = 8
  integer :: N_itt

  ! arguments
  integer :: single_B = 0
  logical :: run_B(3) = .true.
  logical :: benchmark(3) = .false.
  logical :: in_check = .false.

  ! Different methods and their timing
  integer, parameter :: N_methods = 6
  logical :: run_methods(N_METHODS)
  real(dp) :: t(0:N_METHODS)

  real(dp) :: t0, t1, k(3), trun
  complex(dp), dimension(:,:,:), allocatable :: zH, zS, zG
  complex(dp), dimension(:,:,:,:), allocatable :: zHS1, zG1, zHS2, zG2
  complex(dp), parameter :: Z = cmplx(0.1_dp, 0.2_dp, dp)

  integer :: iu
  integer :: N, i
  integer :: B1, B2, B3

#ifdef _OPENMP
  write(*,'(a,i0)') 'OpenMP version ', _OPENMP
!$OMP parallel
!$OMP single
  N = omp_get_num_threads()
!$OMP end single
!$OMP end parallel
  write(*,'(a,i0)') 'OpenMP threads ', N
  call omp_get_schedule(i,N)
  select case ( i )
  case ( OMP_SCHED_STATIC )
    write(*,'(a,i0)') '* OpenMP runtime schedule STATIC, chunks ',N
  case ( OMP_SCHED_DYNAMIC )
    write(*,'(a,i0)') '* OpenMP runtime schedule DYNAMIC, chunks ',N
    if ( N == 1 ) then
      ! this is the default scheduling, probably the user
      ! have not set the value, predefine it to 32
      N = 32
      write(*,'(a,i0)')'** OpenMP runtime schedule DYNAMIC, chunks ',N
    end if
  case ( OMP_SCHED_GUIDED )
    write(*,'(a,i0)') '* OpenMP runtime schedule GUIDED, chunks ',N
  case ( OMP_SCHED_AUTO )
    write(*,'(a,i0)') '* OpenMP runtime schedule AUTO, chunks ',N
  case default
    write(*,'(a,i0)') '* OpenMP runtime schedule UNKNOWN, chunks ',N
  end select
  call omp_set_schedule(i, N)
#endif

  ! Parse arguments to program
  ! Run only for ~ 0.001 seconds per method
  trun = 0.01_dp
  call parse_args()

  iu = 13524
  k = [0.25_dp, 0.25_dp, 0.25_dp]

  t0 = 0._dp
  t1 = 0._dp

  ! initialize BLAS library (in case there are initialization routines)
  call dcopy(1, t0, 1, t1, 1)
  call zcopy(1, k, 1, t, 1)

  run_methods(1) = .true. ! do
  run_methods(2) = .true. ! manual
  run_methods(3) = .true. ! task
  run_methods(4) = .true. ! workshare
  run_methods(5) = .true. ! task-bad
  run_methods(6) = .true. ! manual-simd

  if ( any(benchmark) ) then
    write(*,'(a,tr1,"[",2(i0,":"),i0,"]")') "Running N in", N_min, N_max, N_step
  end if

  in_check = .false.
  if ( benchmark(1) ) then
    if ( run_B(1) ) call benchmark_M_b1()
    if ( run_B(2) ) call benchmark_M_b2()
  end if
  if ( benchmark(2) ) then
    if ( run_B(1) ) call benchmark_HS_b1()
    if ( run_B(2) ) call benchmark_HS_b2()
  end if
  if ( benchmark(3) ) then
    if ( run_B(1) ) call benchmark_HS_G_b1()
    if ( run_B(2) ) call benchmark_HS_G_b2()
  end if

  if ( .not. any(benchmark) ) then
    in_check = .true.
    if ( run_B(1) ) call check_b1()
    if ( run_B(2) ) call check_b2()
  end if

contains

  subroutine parse_args()
    integer :: Narg, iarg
    character(len=64) :: arg

    Narg = command_argument_count()
    iarg = 1
    do while ( iarg <= Narg )

      call get_command_argument(iarg, arg)
      iarg = iarg + 1

      select case ( trim(arg) )
      case ( '-skip-b1', '-skip-B1' )
        run_B(1) = .false.
      case ( '-skip-b2', '-skip-B2' )
        run_B(2) = .false.
      case ( '-M' )
        benchmark(1) = .true.
      case ( '-HS' )
        benchmark(2) = .true.
      case ( '-HSG', '-HS_G', '-HS-G' )
        benchmark(3) = .true.
      case ( '-Nmin' )
        call get_command_argument(iarg, arg)
        read(arg, *) N_min
        iarg = iarg + 1
      case ( '-Nmax' )
        call get_command_argument(iarg, arg)
        read(arg, *) N_max
        iarg = iarg + 1
      case ( '-Nstep' )
        call get_command_argument(iarg, arg)
        read(arg, *) N_step
        iarg = iarg + 1
      case ( '-t' )
        call get_command_argument(iarg, arg)
        read(arg, *) trun
        iarg = iarg + 1
      case default
        call get_command_argument(iarg, arg)
        write(*,*) "Unknown option: ", arg
        stop
      end select

    end do

  end subroutine parse_args

  subroutine benchmark_M_b1()
    B2 = 1
    B3 = 1
    write(*,*) 'Benchmarking rank(B) == 1 [M] examples...'

    t(:) = 0._dp
    open(iu, file="M.B1.timings", status="unknown", action="write")
    write(iu, '(a,a7,tr1,a7,10(tr1,a12))') '#', 'N', 'B1', 'simple [s]', 'do [s]', 'manual [s]', &
        'task [s]', 'work [s]', 'task-bad [s]', 'man-simd [s]'

    B1 = 1
    do while ( B1 < 20 )
      write(*,'(tr2,i0)', advance='no') B1
      do N = N_min, N_max, N_step
        call bu%initialize([B1,1,1])
        allocate(zG(N, N, bu%size()))
        allocate(zG1(N, bu%size(), N, bu%size()))
        allocate(zG2(N, bu%size(), N, bu%size()))

        do i = 1, bu%size()
          call initialize(N, zG(:,:,i))
        end do
        zG1(:, :, :, :) = 0._dp
        zG2(:, :, :, :) = 0._dp

        N_itt = calc_iterations(N, B1, B2, B3)

        t0 = my_time()
        do i = 1, N_itt
          call bu%unfold_M_original(k, N, zG, zG1)
        end do
        t1 = my_time()
        t(0) = (t1 - t0) / N_itt
        !print *, t(0), t(0) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

        if ( run_methods(1) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_do(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "do")
          t(1) = (t1 - t0) / N_itt
          !  print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
        end if

        if ( run_methods(2) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_manual(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "manual")
          t(2) = (t1 - t0) / N_itt
        end if

        if ( run_methods(3) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_task(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "task")
          t(3) = (t1 - t0) / N_itt
        end if

        if ( run_methods(4) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_workshare(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "workshare")
          t(4) = (t1 - t0) / N_itt
        end if

        if ( run_methods(5) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_task_bad(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "taskbad")
          t(5) = (t1 - t0) / N_itt
        end if

        if ( run_methods(6) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_manual_simd(k, N, zG, zG2)
          end do
          t1 = my_time()
          call cmp(zG1, zG2, "man-simd")
          t(6) = (t1 - t0) / N_itt
        end if

        ! Write out data
        write(iu, '(2(tr1,i7),10(tr1,e12.7))') N, B1, t

        deallocate(zG, zG1, zG2)
      end do
      write(iu, *) ! newline
      call step_B(B1)
    end do
    close(iu)
    write(*,*) ! newline

  end subroutine benchmark_M_b1

  subroutine benchmark_HS_b1()
    B2 = 1
    B3 = 1
    write(*,*) 'Benchmarking rank(B) == 1 [HS] examples...'

    t(:) = 0._dp
    open(iu, file="HS.B1.timings", status="unknown", action="write")
    write(iu, '(a,a7,tr1,a7,10(tr1,a12))') '#', 'N', 'B1', 'simple [s]', 'do [s]', 'manual [s]'

    B1 = 1
    do while ( B1 < 20 )
      write(*,'(tr2,i0)', advance='no') B1
      do N = N_min, N_max, N_step
        call bu%initialize([B1,1,1])
        allocate(zH(N, N, bu%size()))
        allocate(zS(N, N, bu%size()))
        allocate(zHS1(N, bu%size(), N, bu%size()))
        allocate(zHS2(N, bu%size(), N, bu%size()))

        do i = 1, bu%size()
          call initialize(N, zH(:,:,i))
          call initialize(N, zS(:,:,i))
        end do
        zHS1(:, :, :, :) = 0._dp
        zHS2(:, :, :, :) = 0._dp

        N_itt = calc_iterations(N, B1, B2, B3)

        t0 = my_time()
        do i = 1, N_itt
          call bu%unfold_HS_original(k, N, zH, zS, Z, zHS1)
        end do
        t1 = my_time()
        t(0) = (t1 - t0) / N_itt
        !print *, t(0), t(0) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

        if ( run_methods(1) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_do(k, N, zH, zS, Z, zHS2)
          end do
          t1 = my_time()
          call cmp(zHS1, zHS2, "do")
          t(1) = (t1 - t0) / N_itt
          !  print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
        end if

        if ( run_methods(2) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_manual(k, N, zH, zS, Z, zHS2)
          end do
          t1 = my_time()
          call cmp(zHS1, zHS2, "manual")
          t(2) = (t1 - t0) / N_itt
        end if

        ! Write out data
        write(iu, '(2(tr1,i7),10(tr1,e12.7))') N, B1, t(:2)

        deallocate(zH, zS, zHS1, zHS2)
      end do
      write(iu, *) ! newline
      call step_B(B1)
    end do
    close(iu)
    write(*,*) ! newline

  end subroutine benchmark_HS_b1

  subroutine benchmark_HS_G_b1()
    B2 = 1
    B3 = 1
    write(*,*) 'Benchmarking rank(B) == 1 [HSG] examples...'

    t(:) = 0._dp
    open(iu, file="HSG.B1.timings", status="unknown", action="write")
    write(iu, '(a,a7,tr1,a7,10(tr1,a12))') '#', 'N', 'B1', 'simple [s]', 'do [s]', 'manual [s]'

    B1 = 1
    do while ( B1 < 20 )
      write(*,'(tr2,i0)', advance='no') B1
      do N = N_min, N_max, N_step
        call bu%initialize([B1,1,1])
        allocate(zH(N, N, bu%size()))
        allocate(zS(N, N, bu%size()))
        allocate(zG(N, N, bu%size()))
        allocate(zHS1(N, bu%size(), N, bu%size()))
        allocate(zG1(N, bu%size(), N, bu%size()))
        allocate(zHS2(N, bu%size(), N, bu%size()))
        allocate(zG2(N, bu%size(), N, bu%size()))

        do i = 1, bu%size()
          call initialize(N, zH(:,:,i))
          call initialize(N, zS(:,:,i))
          call initialize(N, zG(:,:,i))
        end do
        zHS1(:, :, :, :) = 0._dp
        zHS2(:, :, :, :) = 0._dp
        zG1(:, :, :, :) = 0._dp
        zG2(:, :, :, :) = 0._dp

        N_itt = calc_iterations(N, B1, B2, B3)

        t0 = my_time()
        do i = 1, N_itt
          call bu%unfold_HS_G_original(k, N, zH, zS, zG, Z, zHS1, zG1)
        end do
        t1 = my_time()
        t(0) = (t1 - t0) / N_itt
        !print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

        if ( run_methods(1) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_G_do(k, N, zH, zS, zG, Z, zHS2, ZG2)
          end do
          t1 = my_time()
          call cmp(zHS1, zHS2, "do")
          call cmp(zG1, zG2, "do")
          t(1) = (t1 - t0) / N_itt
          !  print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
        end if

        if ( run_methods(2) ) then
          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_G_manual(k, N, zH, zS, zG, Z, zHS2, zG2)
          end do
          t1 = my_time()
          call cmp(zHS1, zHS2, "manual")
          call cmp(zG1, zG2, "manual")
          t(2) = (t1 - t0) / N_itt
        end if

        ! Write out data
        write(iu, '(2(tr1,i7),10(tr1,e12.7))') N, B1, t(:2)

        deallocate(zH, zS, zG, zHS1, zHS2, zG1, zG2)
      end do
      write(iu, *) ! newline
      call step_B(B1)
    end do
    close(iu)
    write(*,*) ! newline

  end subroutine benchmark_HS_G_b1

  subroutine check_b1()
    integer, parameter :: max_N = 300
    integer :: maxB(3), startB(3)
    write(*,'(a)') 'Checking rank(B) == 1 examples...'
    startB(:) = 1
#ifdef __BLOCH_UNFOLD_M_DEBUG
    maxB(:) = 3
#else
    maxB(:) = 20
#endif

    B3 = startB(3)
    do while ( B3 <= maxB(3) )
      B2 = startB(2)
      do while ( B2 <= maxB(2) )
        B1 = startB(1)
        do while ( B1 <= maxB(1) )
          write(*,'(tr2,i0)', advance='no') B1*B2*B3
          ! Just have an N different from B
          ! Just have an N different from B
          N = max_N / B1
          N = N / B2
          N = N / B3
          call bu%initialize([B1,B2,B3])
          allocate(zH(N, N, bu%size()))
          allocate(zS(N, N, bu%size()))
          allocate(zG(N, N, bu%size()))
          allocate(zG1(N, bu%size(), N, bu%size()))
          allocate(zG2(N, bu%size(), N, bu%size()))
          allocate(zHS1(N, bu%size(), N, bu%size()))
          allocate(zHS2(N, bu%size(), N, bu%size()))

          do i = 1, bu%size()
            call initialize(N, zH(:,:,i))
            call initialize(N, zS(:,:,i))
            call initialize(N, zG(:,:,i))
          end do

          call bu%unfold_M_original(k, N, zG, zG1)
          call bu%unfold_M_do(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[do]")
          call bu%unfold_M_manual(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[manual]")
          call bu%unfold_M_task(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[task]")
          call bu%unfold_M_workshare(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[work]")
          call bu%unfold_M_task_bad(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[taskbad]")
          call bu%unfold_M_manual_simd(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[man-simd]")


          ! HS_G
          call bu%unfold_HS_G_original(k, N, zH, zS, zG, Z, zHS1, zG1)
          call bu%unfold_HS_G_do(k, N, zH, zS, zG, Z, zHS2, zG2)
          call cmp(zG1, zG2, "G[do]")
          call cmp(zHS1, zHS2, "HSG[do]")
          call bu%unfold_HS_G_manual(k, N, zH, zS, zG, Z, zHS2, zG2)
          call cmp(zG1, zG2, "G[manual]")
          call cmp(zHS1, zHS2, "HSG[manual]")


          ! HS
          call bu%unfold_HS_original(k, N, zH, zS, Z, zHS1)
          call bu%unfold_HS_do(k, N, zH, zS, Z, zHS2)
          call cmp(zHS1, zHS2, "HS[do]")
          call bu%unfold_HS_do(k, N, zH, zS, Z, zHS2)
          call cmp(zHS1, zHS2, "HS[manual]")

          deallocate(zH, zS, zG, zHS1, zG1, zHS2, zG2)

          if ( B3 * B2 > 1 ) then
            B1 = huge(1)
          else
            call step_B(B1)
          end if
        end do
        if ( B3 > 1 ) then
          B2 = huge(1)
        else
          call step_B(B2)
        end if
      end do
      call step_B(B3)
    end do
    write(*,*) ! newline
    write(*,*) ! newline

  end subroutine check_b1

  subroutine benchmark_M_b2()
    B3 = 1
    t(:) = 0._dp

    write(*,*) 'Benchmarking rank(B) == 2 [M] examples...'

    open(iu, file="M.B2.timings", status="unknown", action="write")
    write(iu, '(a,a7,2(tr1,a7),10(tr1,a12))') '#', 'N', 'B1', 'B2', 'simple [s]', 'do [s]', 'manual [s]', &
        'task [s]', 'work [s]', 'taskbad [s]', 'man-simd [s]'

    B2 = 2
    do while ( B2 <= 9 )
      write(*,'(tr2,i0,":")') B2
      B1 = 2
      write(*,'(tr2)', advance='no')
      do while ( B1 <= 9 )
        write(*,'(tr2,i0)', advance='no') B1
        do N = N_min, N_max, N_step
          call bu%initialize([B1,B2,1])
          allocate(zG(N, N, bu%size()))
          allocate(zG1(N, bu%size(), N, bu%size()))
          allocate(zG2(N, bu%size(), N, bu%size()))

          do i = 1, bu%size()
            call initialize(N, zG(:,:,i))
          end do
          zG1(:, :, :, :) = 0._dp
          zG2(:, :, :, :) = 0._dp

          N_itt = calc_iterations(N, B1, B2, B3)

          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_M_original(k, N, zG, zG1)
          end do
          t1 = my_time()
          t(0) = (t1 - t0) / N_itt
!          print *, t(0), t(0) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

          if ( run_methods(1) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_do(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "do")
            t(1) = (t1 - t0) / N_itt
!            print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
          end if

          if ( run_methods(2) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_manual(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "manual")
            t(2) = (t1 - t0) / N_itt
          end if

          if ( run_methods(3) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_task(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "task")
            t(3) = (t1 - t0) / N_itt
          end if

          if ( run_methods(4) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_workshare(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "workshare")
            t(4) = (t1 - t0) / N_itt
          end if

          if ( run_methods(5) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_task_bad(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "taskbad")
            t(5) = (t1 - t0) / N_itt
          end if

          if ( run_methods(6) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_M_manual_simd(k, N, zG, zG2)
            end do
            t1 = my_time()
            call cmp(zG1, zG2, "man-simd")
            t(6) = (t1 - t0) / N_itt
          end if

          ! Write out data
          write(iu, '(3(tr1,i7),10(tr1,e12.7))') N, B1, B2, t

          deallocate(zG, zG1, zG2)
        end do
        write(iu, *) ! newline
        call step_B(B1)
      end do
      call step_B(B2)
      write(*,*) ! newline
    end do
    close(iu)

  end subroutine benchmark_M_b2

  subroutine benchmark_HS_b2()
    B3 = 1
    t(:) = 0._dp

    write(*,*) 'Benchmarking rank(B) == 2 [HS] examples...'

    open(iu, file="HS.B2.timings", status="unknown", action="write")
    write(iu, '(a,a7,2(tr1,a7),10(tr1,a12))') '#', 'N', 'B1', 'B2', 'simple [s]', 'do [s]', 'manual [s]'

    B2 = 2
    do while ( B2 <= 9 )
      write(*,'(tr2,i0,":")') B2
      B1 = 2
      write(*,'(tr2)', advance='no')
      do while ( B1 <= 9 )
        write(*,'(tr2,i0)', advance='no') B1
        do N = N_min, N_max, N_step
          call bu%initialize([B1,B2,1])
          allocate(zH(N, N, bu%size()))
          allocate(zS(N, N, bu%size()))
          allocate(zHS1(N, bu%size(), N, bu%size()))
          allocate(zHS2(N, bu%size(), N, bu%size()))

          do i = 1, bu%size()
            call initialize(N, zH(:,:,i))
            call initialize(N, zS(:,:,i))
          end do
          zHS1(:, :, :, :) = 0._dp
          zHS2(:, :, :, :) = 0._dp

          N_itt = calc_iterations(N, B1, B2, B3)

          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_original(k, N, zH, zS, Z, zHS1)
          end do
          t1 = my_time()
          t(0) = (t1 - t0) / N_itt
!          print *, t(0), t(0) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

          if ( run_methods(1) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_HS_do(k, N, zH, zS, Z, zHS2)
            end do
            t1 = my_time()
            call cmp(zHS1, zHS2, "do")
            t(1) = (t1 - t0) / N_itt
!            print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
          end if

          if ( run_methods(2) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_HS_manual(k, N, zH, zS, Z, zHS2)
            end do
            t1 = my_time()
            call cmp(zHS1, zHS2, "manual")
            t(2) = (t1 - t0) / N_itt
          end if

          ! Write out data
          write(iu, '(3(tr1,i7),10(tr1,e12.7))') N, B1, B2, t(:2)

          deallocate(zH, zS, zHS1, zHS2)
        end do
        write(iu, *) ! newline
        call step_B(B1)
      end do
      call step_B(B2)
      write(*,*) ! newline
    end do
    close(iu)

  end subroutine benchmark_HS_b2

  subroutine benchmark_HS_G_b2()
    B3 = 1
    t(:) = 0._dp

    write(*,*) 'Benchmarking rank(B) == 2 [HS_G] examples...'

    open(iu, file="HSG.B2.timings", status="unknown", action="write")
    write(iu, '(a,a7,2(tr1,a7),10(tr1,a12))') '#', 'N', 'B1', 'B2', 'simple [s]', 'do [s]', 'manual [s]'

    B2 = 2
    do while ( B2 <= 9 )
      write(*,'(tr2,i0,":")') B2
      B1 = 2
      write(*,'(tr2)', advance='no')
      do while ( B1 <= 9 )
        write(*,'(tr2,i0)', advance='no') B1
        do N = N_min, N_max, N_step
          call bu%initialize([B1,B2,1])
          allocate(zH(N, N, bu%size()))
          allocate(zS(N, N, bu%size()))
          allocate(zG(N, N, bu%size()))
          allocate(zHS1(N, bu%size(), N, bu%size()))
          allocate(zG1(N, bu%size(), N, bu%size()))
          allocate(zHS2(N, bu%size(), N, bu%size()))
          allocate(zG2(N, bu%size(), N, bu%size()))

          do i = 1, bu%size()
            call initialize(N, zH(:,:,i))
            call initialize(N, zS(:,:,i))
            call initialize(N, zG(:,:,i))
          end do
          zHS1(:, :, :, :) = 0._dp
          zHS2(:, :, :, :) = 0._dp
          zG1(:, :, :, :) = 0._dp
          zG2(:, :, :, :) = 0._dp

          N_itt = calc_iterations(N, B1, B2, B3)

          t0 = my_time()
          do i = 1, N_itt
            call bu%unfold_HS_G_original(k, N, zH, zS, zG, Z, zHS1, zG1)
          end do
          t1 = my_time()
          t(0) = (t1 - t0) / N_itt
!          print *, t(0), t(0) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0

          if ( run_methods(1) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_HS_G_do(k, N, zH, zS, zG, Z, zHS2, zG2)
            end do
            t1 = my_time()
            call cmp(zHS1, zHS2, "do")
            call cmp(zG1, zG2, "do")
            t(1) = (t1 - t0) / N_itt
!            print *, t(1), t(1) / calc_operations(N, B1, B2, B3), N_itt, t1 - t0
          end if

          if ( run_methods(2) ) then
            t0 = my_time()
            do i = 1, N_itt
              call bu%unfold_HS_G_manual(k, N, zH, zS, zG, Z, zHS2, zG2)
            end do
            t1 = my_time()
            call cmp(zHS1, zHS2, "manual")
            call cmp(zG1, zG2, "manual")
            t(2) = (t1 - t0) / N_itt
          end if

          ! Write out data
          write(iu, '(3(tr1,i7),10(tr1,e12.7))') N, B1, B2, t(:2)

          deallocate(zH, zS, zG, zHS1, zHS2, zG1, zG2)
        end do
        write(iu, *) ! newline
        call step_B(B1)
      end do
      call step_B(B2)
      write(*,*) ! newline
    end do
    close(iu)

  end subroutine benchmark_HS_G_b2

  subroutine check_b2()
    integer, parameter :: max_N = 1000
    integer :: maxB(3), startB(3)
    write(*,'(a)') 'Checking rank(B) == 2 examples...'
    startB(:) = 1
#ifdef __BLOCH_UNFOLD_M_DEBUG
    maxB(:) = 3
#else
    maxB(:) = 9
#endif

    B3 = startB(3)
    do while ( B3 <= maxB(3) )
      B2 = startB(2)
      do while ( B2 <= maxB(2) )
        if ( B2 == 1 .or. B3 == 1 ) then
          B1 = 2
        else
          B1 = startB(1)
        end if
        do while ( B1 <= maxB(1) .and. (B2 > 1 .or. B3 > 1) )
          write(*,'(tr2,"[",2(i0,","),i0,"]")', advance='no') B1,B2,B3

          ! Just have an N different from B
          N = max_N / B1
          N = N / B2
          N = N / B3
          call bu%initialize([B1,B2,B3])
          allocate(zH(N, N, bu%size()))
          allocate(zS(N, N, bu%size()))
          allocate(zG(N, N, bu%size()))
          allocate(zG1(N, bu%size(), N, bu%size()))
          allocate(zG2(N, bu%size(), N, bu%size()))
          allocate(zHS1(N, bu%size(), N, bu%size()))
          allocate(zHS2(N, bu%size(), N, bu%size()))

          do i = 1, bu%size()
            call initialize(N, zH(:,:,i))
            call initialize(N, zS(:,:,i))
            call initialize(N, zG(:,:,i))
          end do

          call bu%unfold_M_original(k, N, zG, zG1)
          call bu%unfold_M_do(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[do]")
          call bu%unfold_M_manual(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[manual]")
          call bu%unfold_M_task(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[task]")
          call bu%unfold_M_workshare(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[workshare]")
          call bu%unfold_M_task_bad(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[taskbad]")
          call bu%unfold_M_manual_simd(k, N, zG, zG2)
          call cmp(zG1, zG2, "M[man-simd]")

          ! HS_G
          call bu%unfold_HS_G_original(k, N, zH, zS, zG, Z, zHS1, zG1)
          call bu%unfold_HS_G_do(k, N, zH, zS, zG, Z, zHS2, zG2)
          call cmp(zG1, zG2, "G[do]")
          call cmp(zHS1, zHS2, "HSG[do]")
          call bu%unfold_HS_G_manual(k, N, zH, zS, zG, Z, zHS2, zG2)
          call cmp(zG1, zG2, "G[manual]")
          call cmp(zHS1, zHS2, "HSG[manual]")


          ! HS
          call bu%unfold_HS_original(k, N, zH, zS, Z, zHS1)
          call bu%unfold_HS_do(k, N, zH, zS, Z, zHS2)
          call cmp(zHS1, zHS2, "HS[do]")
          call bu%unfold_HS_manual(k, N, zH, zS, Z, zHS2)
          call cmp(zHS1, zHS2, "HS[manual]")

          deallocate(zH, zS, zG, zHS1, zG1, zHS2, zG2)

          if ( B3 > 1 .and. B2 > 1 ) then
            B1 = huge(1)
          else
            call step_B(B1)
          end if
        end do
        call step_B(B2)
      end do
      call step_B(B3)
      write(*,*) ! newline
    end do
    write(*,*) ! newline

  end subroutine check_b2

  subroutine step_B(B)
    integer, intent(inout) :: B
    ! step exponentially (this will take all from 1 - 4)
    B = B + max(1,B / 2)
  end subroutine step_B

  function calc_iterations(N,B1,B2,B3) result(nitt)
    integer, intent(in) :: N, B1, B2, B3
    real(dp) :: ops
    integer :: nitt
    ops = calc_operations(N,B1,B2,B3)
    ! We want each method to run for ~ trun seconds
    nitt = max(nint(trun / (ops * 1.e-9_dp)), 1)
  end function calc_iterations

  function calc_operations(N,B1,B2,B3) result(ops)
    integer, intent(in) :: N, B1, B2, B3
    integer :: mult_ops, copy_ops
    real(dp) :: ops

    ! Number of blocks where we calculate stuff
    mult_ops = (2 * B3 - 1) * (2 * B2 - 1) * (2 * B1 - 1)
    if ( mult_ops == 1 ) mult_ops = 0

    copy_ops = (B1 * B2 * B3) ** 2 - mult_ops

    ! Then we have all the matrices we need to do the
    ! operations on
    mult_ops = mult_ops * B1 * B2 * B3

    ! This factor scheme makes the timing per element roughly the same for
    ! B2 == B3 == 1
    ! Need to check for only B3 == 1
    ops = (real(mult_ops * 4, dp) + real(copy_ops, dp)) * N ** 2

  end function calc_operations

  subroutine cmp(A1, A2, msg)
    complex(dp), intent(inout) :: A1(:,:,:,:), A2(:,:,:,:)
    character(len=*), intent(in) :: msg

    integer :: lo(4), loB(3,2)
    real(dp) :: diff, abs1, abs2

    diff = maxval(abs(A1 - A2))
    if ( diff > 1.e-14_dp ) then
      lo = maxloc(abs(A1 - A2))
      call bu%unravel_index(lo(2), loB(:,1))
      call bu%unravel_index(lo(4), loB(:,2))
      abs1 = abs(A1(lo(1),lo(2),lo(3),lo(4)))
      abs2 = abs(A2(lo(1),lo(2),lo(3),lo(4)))
      if ( diff > 1.e-13_dp ) then
        write(*,'(/a,tr1,i0,tr2,"[",2(i0,","),i0,"]",tr3,i0,",",i0,2(tr1,"[",2(i0,","),i0,"]"),tr2,3(tr1,e15.5)/)') &
            'ERROR '//trim(msg), N,B1, B2, B3, lo(1),lo(3),loB(:,1),loB(:,2), abs1, abs2, diff
      end if
      if ( abs1 > 1.e-14_dp .and. abs2 > 1.e-14_dp .and. diff / abs1 > 1.e-10_dp ) stop
    end if
    if ( in_check ) then
      call initialize(size(A2, 1) * size(A2, 2), A2)
    else
      A2(:,:,:,:) = -2._dp
    end if

  end subroutine cmp

  subroutine initialize(N, zA)
    integer, intent(in) :: N
    complex(dp), intent(inout) :: zA(N,N)
    integer :: i, j
    ! Seed for LAPACK *larnv
    integer, save :: seed(4) = [124469, 54327, 46457, 2903583]

    ! real and imag in [-1,1]
    call zlarnv(2, seed, N*N, zA)

  end subroutine initialize

  function my_time() result(t)
    real(dp) :: t
#ifdef _OPENMP
!$  t = omp_get_wtime()
#else
    call cpu_time(t)
#endif
  end function my_time

end program bloch_unfold

subroutine die(str)
  character(len=*) :: str
  print *, str
  stop
end subroutine die
#endif
