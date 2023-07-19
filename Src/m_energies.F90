module m_energies
  use precision, only: dp
	implicit none

  public

  real(dp):: DEharr
  real(dp):: DEna       ! Neutral-atom energy term, calculated  in dnaefs
  real(dp):: DUext      ! Interaction energy with external  electric field,
                        ! calculated in dhscf

  real(dp):: DUscf      ! Electrostatic energy of (rhoscf- rhoatm), calc. in dhscf
  real(dp):: Dxc        ! Integral((epsxc-Vxc)*Rho), calculated  in dhscf (not used)
  real(dp):: Ecorrec    ! Energy term eta*DeltaQ, calculated in  ordern
  real(dp):: ef         ! Fermi energy
  real(dp):: Eharrs     ! Harris-functional total energy
  real(dp):: Eharrs1    ! Same as Eharrs, but preserved in grid  cell sampling
  real(dp):: Eions      ! Self-energy of isolated ions
  real(dp):: Ekin       ! Kinetic energy of electrons,  calculated in kinefsm
  real(dp):: Ekinion    ! Kinetic energy of ions
  real(dp):: Elast      ! Total energy in the previous SCF  iteration
  real(dp):: Emad       ! Madelung energy term, calculated in  madelung
  real(dp):: Ena        ! Neutral-atom term in the total energy,  calculated in naefs
  real(dp):: Enaatm     ! Integral of Vna * rhoatm, calculated  in dhscf
  real(dp):: Enascf     ! Integral of Vna * rhoscf, calculated  in dhscf
  real(dp):: Enl        ! Non-local pseudopot. energy term,  calculated in nlefsm
  real(dp):: Emeta      ! Metadynamics energy contribution  calculated in meta
  real(dp):: Entrop     ! Temporary to call diagon
  real(dp):: Entropy    ! Entropy due to electron state  occupations, calc. in diagon
  real(dp):: Etot       ! Total electronic energy
  real(dp):: Exc        ! Exchange-correlation energy,  calculated in dhscf
  real(dp):: E0         ! Non-SCF part of total energy
  real(dp):: Emm        ! Classical two-body term, calculated in  twobody
  real(dp):: FreeE      ! Free energy
  real(dp):: FreeEharris! Free energy computed with Harris total energy
  real(dp):: Uatm       ! Harris hartree electron energy,  calculated in dhscf
  real(dp):: Uscf       ! SCF hartree electron energy,  calculated in dhscf
  real(dp):: Ebs        ! Band-structure energy, Tr(DM*H), calculated in compute_dm

end module m_energies



