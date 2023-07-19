      MODULE m_post_scf_work
      private
      public :: post_scf_work

      CONTAINS

      subroutine post_scf_work( first, last, iscf )
      USE m_setup_hamiltonian, only: setup_hamiltonian
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: lasto, rmaxo, datm, indxuo, no_s, no_u, iphorb
      use m_energies
      use neighbour,   only: maxna=>maxnna   ! For plcharge...
      use m_spin,         only: nspin
      implicit none

      logical :: first, last
      integer :: iscf

      external :: plcharge

!------------------------------------------------------------------------- BEGIN
          ! If converged, make one last iteration to find forces and stress

      call setup_hamiltonian( first, last, iscf )

! Print energies
      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad+Emm+
     .     Emeta
      Eharrs = Etot + DEharr
      FreeE  = Etot - Temp * Entropy

!
!     WHY IS THIS?? To make Dscf (which might be re-used in the next
!     geometry iteration) equal to the last DM used to construct the
!     hamiltonian, instead of the output DM.
!
      Dscf(1:nh,1:nspin) = Dold(1:nh,1:nspin)
      Escf(1:nh,1:nspin) = Eold(1:nh,1:nspin)

      if (dumpcharge) then
         call plcharge( no_s, na_s, no_u, maxnh, maxna, nspin,
     .        isa, iphorb, indxuo, lasto,
     .        scell, nsc, xa, rmaxo, datm )
      endif
      
!--------------------------------------------------------------------------- END
      END subroutine post_scf_work

      END MODULE m_post_scf_work
