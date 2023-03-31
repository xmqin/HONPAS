module m_compute_ebs_shift
implicit none
public :: compute_ebs_shift
CONTAINS
  subroutine compute_ebs_shift(Dscf,H,Hprev,delta_Ebs)
    use precision,   only: dp
    use m_mpi_utils, only: globalize_sum
    use parallel,    only: Siesta_comm

    real(dp), intent(in) :: Dscf(:,:)
    real(dp), intent(in) :: H(:,:)
    real(dp), intent(in) :: Hprev(:,:)
    real(dp), intent(out) :: delta_Ebs

    integer  :: ispin, nspin, io, nnzs
    real(dp) :: const
#ifdef MPI
    real(dp) :: buffer1
#endif

    nnzs = size(Dscf,dim=1)
    nspin = size(Dscf,dim=2)

    delta_Ebs = 0.0_dp
    do ispin = 1,nspin
!       const factor takes into account that there are two nondiagonal
!       elements in non-collinear spin density matrix, stored as
!       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
!       ispin=4 => Imag(D12)
       const = 1._dp
       if (ispin .gt. 2) const = 2._dp
       do io = 1,nnzs
          delta_Ebs = delta_Ebs + (H(io,ispin)-Hprev(io,ispin)) * const * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
      call globalize_sum(delta_Ebs,buffer1,comm=Siesta_comm)
      delta_Ebs = buffer1
#endif
    end subroutine compute_ebs_shift
  end module m_compute_ebs_shift
