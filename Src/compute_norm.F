! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_norm
      CONTAINS
      function compute_norm(DM,nnz,nspin) result(norm)
      use precision, only: dp
      use sparse_matrices, only: S
      use m_spin, only: spinor_dim
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif

      integer, intent(in)   ::  nnz, nspin
      real(dp), intent(in)  ::  DM(:,:)
      real(dp)  ::       norm

#ifdef MPI
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif
      integer   ::  io, ispin

      norm = 0.0_dp
      do ispin = 1,spinor_dim
        do io = 1,nnz
          norm = norm + DM(io,ispin) * s(io)
        enddo
      enddo
#ifdef MPI
        call globalize_sum(norm,buffer1)
        norm = buffer1
#endif
      end function compute_norm
      end module m_norm
