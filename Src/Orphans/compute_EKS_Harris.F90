module m_EKS_Harris
  implicit none
  public :: compute_EKS_Harris
CONTAINS
  subroutine  compute_EKS_Harris (E_Harris, E_KS_Good)

    ! This routine computes the Harris energy from E_KS(DM_in):
    !
    !    E_Harris = E_KS(DM_in) + Tr[H_in*(DM_out-DM_in)]
    !
    ! and E_KS(DM_out) as
    !
    !    E_KS(DM_out) = Tr[H_in*DM_out] + E_HXC(DM_out)
    !
    ! Note that E_KS(DM_in), as computed in setup_hamiltonian, is not variational, 
    ! as the kinetic energy term is using a DM (DM_in) which is not derived from wave-functions.
    !
    ! E_KS(DM_out) computed in setup_hamiltonian with DM_out is variational. Rather than 
    ! calling again setup_hamiltonian, it is enough to compute the "scf" part of the energy
    ! by calling dhscf with DM_out.

      use precision,       only: dp
      use fdf,             only: fdf_get
      use siesta_options,  only: g2cut
      use sparse_matrices, only: H_kin_1D, H_vkb_1D
      use sparse_matrices, only: listh, listhptr, numh, maxnh
      use sparse_matrices, only: H
      use sparse_matrices, only: Dscf, Dold
      use class_dSpData1D,  only: val
      use m_dhscf,         only: dhscf
      use m_energies
      use atomlist,        only: no_u, iaorb, iphkb, qtot, indxuo, datm,   &
                                 lastkb, no_s, rmaxv, indxua, iphorb, lasto, &
                                 rmaxo, no_l
      use m_ntm,           only: ntm
      use m_spin,          only: nspin
      use m_dipol,         only: dipol
      use siesta_geom,     only: na_u, na_s, xa, isa
#ifdef MPI
      use m_mpi_utils,     only: globalize_sum
#endif


      real(dp), intent(out) :: E_Harris
      real(dp), intent(out) :: E_KS_Good

      integer               :: ihmat, ifa, istr, ispin, io
      real(dp), pointer     :: H_vkb(:), H_kin(:)
#ifdef MPI
      real(dp) :: buffer1
#endif

      real(dp) :: const, Escf_out
      real(dp) :: dummy_stress(3,3), dummy_fa(1,1)
      real(dp) :: dummy_E, g2max, dummy_H(1,1)

!     Compute the band-structure energy and the correction for E_Harris

      DEharr = 0.0_dp
      Ebs = 0.0_dp
      do ispin = 1,nspin
!       const factor takes into account that there are two nondiagonal
!       elements in non-collinear spin density matrix, stored as
!       ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
!       ispin=4 => Imag(D12)
        const = 1._dp
        if (ispin .gt. 2) const = 2._dp
        do io = 1,maxnh
          DEharr = DEharr + H(io,ispin) * const * ( Dscf(io,ispin) - Dold(io,ispin) )
          Ebs = Ebs + H(io,ispin) * const * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
!     Global reduction of DEharr
      call globalize_sum( DEharr, buffer1 )
      DEharr = buffer1
      call globalize_sum(Ebs,buffer1)
      Ebs = buffer1
#endif

      ! These energies were calculated in the latest call to
      ! setup_hamiltonian, using as ingredient D_in 

      ! Ecorrec comes from O(N)...
      ! DUext is the energy of the charge from DM_in in a field.
      ! Emad, Emm, Emeta are extra terms that are added for
      ! consistency of the total energy.

      DEna = Enascf - Enaatm
      Etot = E0 + DEna + DUscf + DUext + Exc + Ecorrec + Emad + Emm + Emeta

      ! This is correct
      E_Harris = Etot + DEharr

      if (fdf_get("DoNotCorrectEKS",.false.)) then
         E_KS_good = Etot
         return
      endif

      ! Now for E_KS(DM_out)

      g2max = g2cut
      ifa  = 0
      istr = 0
      ihmat = 0

      ! Pass DM_out to compute E_HXC(out)

      ! Remove unwanted arguments...

      call dhscf( nspin, no_s, iaorb, iphorb, no_l,                         &
                  no_u, na_u, na_s, isa, xa, indxua,                        &
                  ntm, ifa, istr, ihmat, ' ', ' ', ' ', ' ', ' ', ' ',      &
                  maxnh, numh, listhptr, listh, Dscf, Datm,                 &
                  maxnh, dummy_H, Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, &
                  Exc, Dxc, dipol, dummy_stress, dummy_fa, dummy_stress)

      DEna     = Enascf - Enaatm

      ! Clarify: 
      !          DUext (external electric field) -- should it be in or out?

      Escf_out = DEna + DUscf + DUext + Exc 

!     Compute Tr[H_0*DM_out] = Ekin + Enl with DM_out

      H_kin => val(H_kin_1D)
      H_vkb => val(H_vkb_1D)
      Ekin = 0.0_dp
      Enl  = 0.0_dp
      do ispin = 1,min(nspin,2)
        do io = 1,maxnh
          Ekin = Ekin + H_kin(io) * Dscf(io,ispin)
          Enl  = Enl  + H_vkb(io) * Dscf(io,ispin)
        enddo
      enddo
#ifdef MPI
!     Global reduction of Ekin, Enl
      call globalize_sum( Ekin, buffer1 )
      Ekin = buffer1
      call globalize_sum( Enl, buffer1 )
      Enl = buffer1
#endif

      ! E0 = Ekin + Enl - Eions + Ena

      ! Clarify: Ecorrec (from O(N))
      !          
      E_KS_good = Ekin + Enl - Eions + Ena + Escf_out + Ecorrec + Emad + Emm + Emeta

    end subroutine compute_EKS_Harris


end module m_EKS_Harris
