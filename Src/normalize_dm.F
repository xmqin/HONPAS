      MODULE m_normalize_dm
      private
      public :: normalize_dm

      CONTAINS

      subroutine normalize_dm( first )
      use precision, only: dp
      use sparse_matrices, only: Dscf, Escf, nh, S
      use atomlist, only: qtot
      use parallel, only: IOnode
      use m_spin,   only: nspin
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      implicit none

      logical :: first
      integer :: io, ispin
      real(dp):: qsol       ! Total unnormalized electron charge

#ifdef MPI
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif

      ! Normalize density matrix to exact charge
      qsol = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,nh
          qsol = qsol + Dscf(io,ispin) * s(io)
        enddo
      enddo
#ifdef MPI
        call globalize_sum(qsol,buffer1)
        qsol = buffer1
#endif
      if (IOnode) then
        if (.not.first .and.
     .     abs(qsol/qtot-1._dp).gt.1.d-2) write(6,'(a,2f15.6)')
     .    'siesta: WARNING: Qtot, Tr[D*S] =', qtot, qsol
      endif
      do ispin = 1,nspin
        do io = 1,nh
          Dscf(io,ispin) = Dscf(io,ispin) * qtot/qsol
          Escf(io,ispin) = Escf(io,ispin) * qtot/qsol
        enddo
      enddo

      END subroutine normalize_dm

      END module m_normalize_dm




