! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! Module to calculate the total charge.

module m_dm_charge

  implicit none

  private
  public :: dm_charge

contains

  subroutine dm_charge(spin, DM_2D, S_1D, Q)
    use precision,       only: dp

    use t_spin, only: tSpin
    use class_Sparsity
    use class_dSpData1D
    use class_dSpData2D
    use class_OrbitalDistribution

#ifdef MPI
    use m_mpi_utils,     only: globalize_sum
#endif

    ! Containers for data
    type(tSpin), intent(in) :: spin
    type(dSpData2D), intent(inout) :: DM_2D
    type(dSpData1D), intent(inout) :: S_1D
    real(dp), intent(out) :: Q

    !  New interface data
    type(Sparsity), pointer :: sp

    integer :: is, ir, ind
    integer :: lnr, nr
    integer, pointer :: lptr(:), ncol(:), lcol(:)
    real(dp), pointer :: DM(:,:), S(:)
#ifdef MPI
    real(dp) :: Qr
#endif

    ! get the distribution
    sp => spar(DM_2D)
    call attach(sp ,n_col=ncol,list_ptr=lptr,nrows=lnr,nrows_g=nr)
    
    S => val(S_1D)
    DM => val(DM_2D)

    ! Initialize the total charge
    Q = 0._dp

!$OMP parallel do default(shared), &
!$OMP& private(is,ir,ind), reduction(+:Q)
    do is = 1, spin%spinor
       do ir = 1, lnr
          do ind = lptr(ir) + 1, lptr(ir) + ncol(ir)
             Q = Q + DM(ind,is) * S(ind)
          end do
       end do
    end do
!$OMP end parallel do

#ifdef MPI
    call globalize_sum(Q, Qr)
    Q = Qr
#endif
    
  end subroutine dm_charge
  
end module m_dm_charge
