! *** Module: nao2gto_prescreen ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Prescreen ERIs using the Schwarz inequality to get list_ij,
!!        and list_kl
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!      - 09.2021 Map dense to sparse matrix for reduing memory [Xinming Qin]       -
! *****************************************************************************
module nao2gto_prescreen

  use precision,  only: dp
  use parallel,   only: Node, Nodes

  implicit none

  private

  public ::                    & 
&   build_pair_list,           &
&   calc_prescreen_eri,        &
&   init_prescreen_eri,        &
&   calc_pair_dist_radii,      &
&   calc_screening_functions

contains

! *****************************************************************************
!> \brief This routine uses the Schwarz inequality to compute the requested values.
!!
!! Following Sec. III A 1 of Ref. \cite Shang:2011, we are going to compute
!! \f$ \left( \mu \nu \vert \mu \nu \right)_{\rm SR} \f$, where according
!! to the heading of the Section, and Eq. (14) of Ref. \cite Shang:2011,
!! corresponds to 
!!
!! \f{eqnarray*}{
!!   \left( \mu \nu \vert \mu \nu \right)_{\rm SR} \equiv 
!!   \left(\phi_{\mu}^{\vec{0}} \phi_{\nu}^{\vec{N}} \vert 
!!         \phi_{\mu} ^{\vec{0}} \phi_{\nu}^{\vec{N}}\right)_{\rm SR}
!!   = \int \int \frac{\phi_{\mu}(\vec{r}-\vec{R}_{\mu})
!!                     \phi_{\nu}(\vec{r}-\vec{N}-\vec{R}_{\nu}) \:\:
!!                     {\rm erfc}(\omega \vert \vec{r} - \vec{r}^{\prime}\vert)
!!                     \:\: \phi_{\mu}(\vec{r}^{\prime}-\vec{R}_{\mu})
!!                     \phi_{\nu}(\vec{r}^{\prime}-\vec{N}-\vec{R}_{\nu})
!!   }
!!     {\vert \vec{r} - \vec{r}^{\prime} \vert} d^{3}r d^{3}r^{\prime},
!! \f}
!!
!! where we use the compact index notation 
!! \f$ \mu \equiv \lbrace{I l m \zeta \rbrace} \f$, 
!! the orbital \f$\mu\f$ is centered
!! at the position \f$ \vec{R}_{\mu} \f$ within the home unit cell,
!! and the orbital \f$ \nu \f$ is centered on the position 
!! \f$ \vec{R}_{\nu} \f$ at a periodic replica of the unit cell within the 
!! auxiliary supercell displaced by the unit cell lattice vector 
!! \f$ \vec{N} \f$.
!!
!! \f{eqnarray*}{
!! \phi_{\mu} (\vec{r}) & = \phi_{Ilm\zeta} (\vec{r}) 
!!   \nonumber \\
!!   & = R_{I l \zeta} (r_{I}) Y_{lm}(\hat{r}_{I})
!!   \nonumber \\
!!   & = \frac{ R_{I l \zeta} (r_{I})}{r_{I}^{l}} 
!!       \left( r_{I} ^{l} Y_{lm}(\hat{r}_{I}) \right)
!!   \nonumber \\
!!   & = \phi_{I l \zeta} (r_{I}) \left( r_{I} ^{l} Y_{lm}(\hat{r}_{I}) \right),
!! \f}
!!
!! where \f$ \vec{r}_{I} = \vec{r} - \vec{R}_{I} \f$,
!! and we have defined \f$ \phi_{I l \zeta} (r_{I}) = 
!! \frac{ R_{I l \zeta} (r_{I})}{r_{I}^{l}} \f$.
!! This radial part is expanded as a linear combination of \f$ N_{\rm G} \f$ 
!! primitive Gaussians. 
!!
!! \f{eqnarray*}{
!! \phi_{I l \zeta} (r_{I}) = \sum_{a=1}^{N_{\rm G}} D_{Il\zeta a} e^{-\alpha_{I l \zeta a} r_{I}^{2}}.
!! \f}
!!
!! It is important to mention here how this expansion
!! is independent of the angular quantum number \f$ m \f$, i.e. 
!! the coefficients and the exponents of the Gaussians that enter in the
!! expansion are the same for the 2\f$p_{x}\f$, and for a 2\f$p_{y}\f$
!! atomic orbital (this will appear below in the subroutine when we 
!! cycle with expressions like "if ( m_i .ne. -l_i )".
!!
!! Therefore
!! \f{eqnarray*}{
!! \phi_{\mu} (\vec{r}) & = \phi_{Ilm\zeta} (\vec{r}) 
!!   \nonumber \\
!!   & = \left(\sum_{a=1}^{N_{\rm G}} D_{Il\zeta a} e^{-\alpha_{I l \zeta a} 
!!        r_{I}^{2}} \right) 
!!       \times \left( r_{I} ^{l} Y_{lm}(\hat{r}_{I}) \right) = 
!!       \sum_{a=1}^{N_{\rm G}} D_{Il\zeta a} 
!!       \left( r_{I} ^{l} Y_{lm}(\hat{r}_{I}) e^{-\alpha_{I l \zeta a} 
!!       r_{I}^{2}} \right)
!! \f}
!! 
!! However, efficient algorithms to calculate integrals over Gaussians, 
!! like the ones implemented in LIBINT, expects the cartesian Gaussian-type 
!! orbitals (CGTO) introduced by Boys, and defined as 
!! \f{eqnarray*}{
!!   G_{I l \zeta a l_x l_y l_z} (\vec{r}) = (x-R_{Ix})^{l_x} (y-R_{Iy})^{l_y} 
!!     (z-R_{Iz})^{l_z} e^{-\alpha_{I l \zeta a} r_{I}^{2}},
!! \f}
!! subject to the condition \f$ l_{x} + l_{y} + l_{z} = l \f$.
!!
!! The pure spherical harmonics \f$ r_{I} ^{l} Y_{lm}(\hat{r}_{I}) 
!!   e^{-\alpha_{I l \zeta a} r_{I}^{2}} \f$ can be constructed from the
!! the appropriate Cartesian Gaussians \cite Schlegel-95
!! \f{eqnarray*}{
!!    r_{I} ^{l} Y_{lm}(\hat{r}_{I}) 
!!   e^{-\alpha_{I l \zeta a} r_{I}^{2}} = 
!!   \sum_{l_{x},l_{y},l_{z} / l_{x} + l_{y} + l_{z} = l}
!!   F_{lml_{x}l_{y}l_{z}} G_{I l \zeta a l_x l_y l_z} (\vec{r})
!! \f}
!!
!! All in all, we end up with the following expansion of a numerical atomic
!! orbital into Cartesian Gaussian Type Orbitals
!!
!! \f{eqnarray*}{
!! \phi_{\mu} (\vec{r}) = \phi_{Ilm\zeta} (\vec{r}) = 
!!    \sum_{a=1}^{N_{\rm G}} 
!!    \sum_{l_{x},l_{y},l_{z} / l_{x} + l_{y} + l_{z} = l}  D_{Il\zeta a} 
!!    F_{lml_{x}l_{y}l_{z}} G_{I l \zeta a l_x l_y l_z} (\vec{r})
!! \f}
!!
!! If we require \f$ N_{\rm G} = six \f$ Gaussians to expand the radial part
!! of an atomic orbital, 
!! if that orbital corresponds to a \f$ p_{x} \f$ (\f$ l = 1, m = +1 \f$), 
!! we will need 18 terms in the expansion (six Gaussians \f$ \times \f$ 
!! three CGTO to expand the real spherical harmonic).
!! If the orbital corresponds to a \f$ d_{xy} \f$ (\f$ l = 2, m = +1 \f$),
!! then we have 36 terms in the expansion (6 Gaussians \f$ \times \f$ 
!! six CGTO to expand the real spherical harmonic)
!! 
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] hfx_sys: data structure storing system information
!! \param[in] maxnh: maximum number of orbitals interacting
!! \param[in] numh: number of nonzero element of each row of the hamiltonian 
!! \param[in] listhptr: pointer to start of each row of the hamiltonian matrix
!! \param[in] listh: nonzero hamiltonian-matrix elements
!! \param[out] max_eri: maximum value of eri_prescreen. It corresponds
!!             to the variable \f$ M_{e} \f$ defined after Eq. (25) in
!!             Ref. \cite Shang:2011
!!
!! NOTE: In parallel runs, and as is implemented now, 
!! all the Nodes know everything about the neighbour list.
!! This can be highly improved for efficienty
! *****************************************************************************
  subroutine calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&   maxnh, numh, listhptr, listh, max_eri)

    use atm_types,      only: species_info  ! Derived type with all the info
                                            !   about the radial functions
                                            !   (PAOs, KB projectors,
                                            !   LDA+U proj,
                                            !   VNA potentials, etc)
                                            !   for a given atomic specie
    use atm_types,      only: species       ! Actual array where the
                                            !   previous information is
                                            !   stored
    use parallel,       only: Node          ! Working Node
    use parallel,       only: Nodes         ! Total number of nodes
    use atmfuncs,       only: lofio         ! Returns the angular momentum
                                            !   quantum number of a given
                                            !   atomic basis orbital
    use atmfuncs,       only: mofio         ! Returns the magnetic quantum
                                            !   number of a given
                                            !   atomic basis orbital
    use atomlist,       only: indxuo        ! Index of equivalent orbital in
                                            !   the unit cell
    use listsc_module,  only: listsc        ! Orbital, equivalent to JUO,
                                            !   which is related to IO
                                            !   exactly as JUO is related to IUO
    use alloc,          only: re_alloc      ! Reallocation routines
    use alloc,          only: de_alloc      ! Deallocation routines
    use nao2gto_contract
    use nao2gto_data,   only: hfx_gtos      ! Information about GTOs
    use nao2gto_data,   only: nco           ! Number of Cartesian Gaussians type
                                            !   orbitals for a given angular 
                                            !   momentum shell:
                                            ! 1 for the s
                                            ! 3 for the p
                                            ! 6 FOR THE d (be careful with this)
                                            !10 FOR THE f (be careful with this)
                                            !    ...       
                                            !   The total number amounts to 
                                            !   (l+1)(l+2)/2
    use nao2gto_data,   only: nso           ! Number of Spherical Gaussians type
                                            !   orbitals for a given angular 
                                            !   momentum shell:
                                            !   1 for the s
                                            !   3 for the p
                                            !   5 for the d 
                                            !   7 for the f 
                                            !    ...       
    use nao2gto_data,   only: D2Sindx
    use nao2gto_data,   only: eri_prescreen ! (phi_mu phi_nu|phi_mu phi_nu). 
                                            !   It appears in Eq. (25) of 
                                            !   Ref. \cite Shang:2011.
                                            !   This integral is computed as 
                                            !   in Eq. (14) of 
                                            !   Ref. \cite Shang:2011.
    use nao2gto_libint
    use nao2gto_pbc,    only: trans_pbc
    use nao2gto_types,  only: hfx_options_type
                                            ! Data type to store Hartree-Fock 
                                            !    exchange options that can be
                                            !    read from SIESTA input files
    use nao2gto_types,  only: hfx_system_type
                                            ! Data type to point 
                                            !    NAO2GTO routines to 
                                            !    relevant SIESTA data
    use nao2gto_types,  only: gto_info_type
                                            ! Data type with the information
                                            !    of the Gaussians required to 
                                            !    expand the radial part of 
                                            !    the Numerical Atomic Orbitals

!   For debugging
#ifdef MPI
    use mpi_siesta
#endif
!   End debugging

    ! Arguments
    type(Libint_t),         intent(inout) :: libint_data
    type(hfx_options_type), intent(in)    :: hfx_opts
                                            ! Data structure storing 
                                            !  Hartree-Fock exchange option
    type(hfx_system_type),  intent(in)    :: hfx_sys
                                            ! Data structure storing 
                                            !   system information
    integer,                intent(in)    :: maxnh
                                            ! Maximum number of orbitals
                                            !   interacting
    integer,                intent(in)    :: numh(hfx_sys%nuotot)
                                            ! Number of nonzero element of
                                            !   each row of the
                                            !   hamiltonian matrix
    integer,                intent(in)    :: listhptr(hfx_sys%nuotot)
                                            ! Pointer to start of each row
                                            !   of the hamiltonian matrix
    integer,                intent(in)    :: listh(maxnh)
                                            ! Nonzero hamiltonian-matrix
                                            !   elements
    real(dp),               intent(out)   :: max_eri
                                            ! Maximum value of eri_prescreen. 
                                            !   It corresponds to the variable 
                                            !   \f$ M_{e} \f$ defined after 
                                            !   Eq.(25) in Ref. \cite Shang:2011

!   Local variables
    integer  :: io        ! Counter for loop on atomic orbitals
    integer  :: j         ! Counter for loop on neighbours
    integer  :: i1        ! Counter for loop 
    integer  :: i2        ! Counter for loop 
    integer  :: ia        ! Atomic index
    integer  :: is        ! Atomic species
    integer  :: ioa       ! Orbital index of each orbital in its atom
    integer  :: joa       ! Orbital index of each neighbour orbital in its atom
    integer  :: l_i       ! Angular momentum of the
                          !   orbital whose neighbours
                          !   are sought
    integer  :: m_i       ! Magnetic quantum number of
                          !   the orbital whose
                          !   neighbours are sought
    integer  :: l_j       ! Angular momentum of the
                          !   neighbour orbital
    integer  :: m_j       ! Magnetic quantum number of
                          !   the neighbour orbital
    integer  :: npgfi     ! Number of Gaussians in the expansion of the radial
                          !    part of a given orbital
                          !    This corresponds to N_G in the 
                          !    doxygen documentation
    integer  :: npgfj     ! Number of Gaussians in the expansion of the radial
                          !    part of a given neighbour orbital
                          !    This corresponds to N_G in the 
                          !    doxygen documentation
    integer  :: ncoi      ! Number of Cartesian Gaussian type orbitals in the 
                          !    expansion of a givel shell of orbitals
                          !    This corresponds to the numbers "18" or "36"
                          !    as explained in the doxygen documentation
    integer  :: ncoj      ! Number of Cartesian Gaussian type orbitals in the 
                          !    expansion of a givel shell of orbitals in the
                          !    neighbonur atom
                          !    This corresponds to the numbers "18" or "36"
                          !    as explained in the doxygen documentation
    integer  :: ind       ! Index of the neighbour
    integer  :: jo        ! Neighbour orbital
    integer  :: ja        ! Atom where to which neighbour orbital belongs to
    integer  :: js        ! Species of the atom where the neighbour orbital 
                          !   is centered

    integer  :: iu, ju
    real(dp) :: ri(3)     ! Position of the atom whose neighbours are sought
    real(dp) :: rj(3)     ! Position of the neighbour orbital
    real(dp) :: r_temp(3) ! Relative position between a given atom and 
                          !   its neighbours, computed from the atomic
                          !   positions in the array XA.
                          !   Remember that, due to the PBC, a folding 
                          !   might have to be applied:
                          !   the actual neighbour might be a periodic replica 
                          !   of the atom included in the neighbour list,
                          !   located in other supercell
    real(dp) :: r_pbc(3)  ! Actual relative position  between a given atom and
                          !   its neighbours, after performing the folding
    real(dp), dimension(:,:,:,:), pointer :: eri => null()
    type(species_info),  pointer :: ispp => null()   ! Pointer to species type
                                                     !   for the atom in the
                                                     !   unit cell
    type(species_info),  pointer :: jspp => null()   ! Pointer to species type
                                                     !   for the neighbour atom
    type(gto_info_type), pointer :: igto => null()   ! Pointer to the Gaussian
                                                     !   expansion of the radial
                                                     !   part of the orbital
                                                     !   in the unit cell
    type(gto_info_type), pointer :: jgto => null()   ! Pointer to the Gaussian
                                                     !   expansion of the radial
                                                     !   part of the neighbour  
                                                     !   orbital

!   For debugging
#ifdef MPI
    integer     :: MPIerror
#endif
!   End debugging

! ------------------------------------------------------------------------------
!!   For debugging
!    write(6,'(a,3i5)') &
! &    'calc_prescreen_eri: Node, Nodes, hfx_sys_nuotot = ', &
! &    Node, Nodes, hfx_sys%nuotot
!!   End debugging

!   Loop over all the orbitals in the unit cell 
!   (mu in the equation of the doxygen documentation)
    do io = 1, hfx_sys%nuotot

!     Identify the atomic index of the orbital
      ia   =  hfx_sys%iaorb(io)

!     Identify the atomic species
      is   =  hfx_sys%isa(ia)

!     Pointer to species type
      ispp => species(is)

!     Pointer to the information on the Gaussian Type Orbitals
      igto => hfx_gtos(is)

!     Identify the orbital index of each orbital in its atom
      ioa  =  hfx_sys%iphorb(io)

!     Identify the orbital angular momentum
      l_i  =  lofio(is,ioa)
!     Identify the magnetic angular momentum
      m_i  =  mofio(is,ioa)

!!     For debugging
!      if( Node .eq. 0 ) then
!        do j = 1, numh(io)
!          ind = listhptr(io) + j
!          jo  = listh(ind)
!          write(6,'(a,4i7)') &
! &        'calc_prescreen_eri: Node, io, j, ind, jo = ',    &
! &                 io, j, ind, jo 
!        enddo
!      endif
!!     End debugging

!     We follow only the first time that a given angular momentum
!     of a given atom appears
!     In other words the loop will continue only for different angular
!     momentum shells
      if ( m_i .ne. -l_i ) cycle

!     Identify the center of the atomic orbital
!     This is R_mu in the equations of the doxygen documentation
      ri(:) = hfx_sys%xa(:,ia)

!     Identify the number of Gaussians that expand the radial part of 
!     a given orbital
      npgfi = igto%orbnl_contract(ispp%orb_index(ioa))

!     Identify the total number of Cartesian Gaussians required to expand 
!     all the atomic orbital of a given shell
      ncoi  = nco(l_i)*npgfi

!!     For debugging
!      write(6,'(a,11i5)') &
! &      'calc_prescreen_eri: Node, Nodes, io, ia, is, ioa, li, mi, npgfi, nco, ncoi = ',  &
! &      Node, Nodes, io, ia, is, ioa, l_i, m_i, npgfi, nco(l_i), ncoi 
!!     End debugging

!     Loop on all the neighbours of a given atomic orbital
!     \nu in the previous doxygen documentation
      do j = 1, numh(io)

!       Index of the neighbour
        ind = listhptr(io) + j
!       Neighbour orbital
!       This is \nu in the Equations in the doxygen documentation
        jo  = listh(ind)
!       Atom where to which neighbour orbital belongs to
        ja  = hfx_sys%iaorb(jo)
!       Species of the atom where the neighbour orbital is centered
        js  = hfx_sys%isa(ja)
!       Pointer to species type of the neighbour orbital
        jspp => species(js)
!       Pointer to the information on the Gaussian Type Orbitals expanding
!       the neighbour orbital
        jgto => hfx_gtos(js)
!       Identify the orbital index of each orbital in its atom
        joa = hfx_sys%iphorb(jo)
!       Identify the neighbour orbital angular momentum 
        l_j = lofio(js,joa)
!       Identify the neighbour orbital magnetic momentum 
        m_j = mofio(js,joa)

!       We follow only the first time that a given angular momentum
!       of a given atom appears
!       In other words the loop will continue only for different angular
!       momentum shells
        if ( m_j .ne. -l_j ) cycle

!       Identify the center of the neighbour orbital
        rj(:)  = hfx_sys%xa(:,ja)

!       Identify the relative position between the neighbour and the atom
!       whose neigbours are sought
        r_temp = rj - ri

!       Compute the actual relative position between the neighbour and
!       the atom, after performing the folding
!       hfx_sys%cell are the supercell lattice vectors in real space
        call trans_pbc(r_temp, hfx_sys%cell, hfx_sys%cell_r, r_pbc)

!       Position of the neighbour atom
!       This is R_\nu in the doxygen documentation
        rj = ri + r_pbc

!       Identify the number of Gaussians that expand the radial part of 
!       a given neighbour orbital
        npgfj = jgto%orbnl_contract(jspp%orb_index(joa))
!       Identify the total number of Cartesian Gaussians required to expand 
!       all the atomic orbital of a given shell
        ncoj = nco(l_j)*npgfj

!!       For debugging
!        write(6,'(a,13i5)') &
! &        'calc_prescreen_eri: Node, Nodes, io, j, ind, jo, maxnh, lj, mj, npgfj, ncoj, nsoli, nsolj = ', &
! &                             Node, Nodes, io, j, ind, jo, maxnh, l_j, m_j, &
! &                             npgfj, ncoj, nso(l_i), nso(l_j) 
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, npgfi  = ', &
! &                             Node, Nodes, npgfi  
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, npgfj  = ', &
! &                             Node, Nodes, npgfj  
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, l_i    = ', &
! &                             Node, Nodes, l_i
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, l_j    = ', &
! &                             Node, Nodes, l_j
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, ncoi   = ', &
! &                             Node, Nodes, ncoi
!        write(6,'(a,3i5)')                              &
! &        'calc_prescreen_eri: Node, Nodes, ncoj   = ', &
! &                             Node, Nodes, ncoj
!        write(6,'(a,2i5,3f12.5)')                       &
! &        'calc_prescreen_eri: Node, Nodes, ri     = ', &
! &                             Node, Nodes, ri(:)
!        write(6,'(a,2i5,3f12.5)')                       &
! &        'calc_prescreen_eri: Node, Nodes, rj     = ', &
! &                             Node, Nodes, rj(:)
!        do i1 = 1, npgfi
!          write(6,'(a,3i5,f12.5)')                       &
! &          'calc_prescreen_eri: Node, Nodes, ipgfi, orbnl_zeta(i1,ispp%orb_index(ioa)) = ', &
! &                               Node, Nodes, i1, igto%orbnl_zeta(i1,ispp%orb_index(ioa))
!        enddo 
!        do i1 = 1, npgfj
!          write(6,'(a,3i5,f12.5)')                       &
! &          'calc_prescreen_eri: Node, Nodes, ipgfj, orbnl_zeta(i1,jspp%orb_index(joa)) = ', &
! &                               Node, Nodes, i1, jgto%orbnl_zeta(i1,jspp%orb_index(joa))
!        enddo 
!        do i1 = 1, ncoi
!          do i2 = ioa, ioa+nso(l_i)-1
!            write(6,'(a,4i5,f12.5)')                       &
! &            'calc_prescreen_eri: Node, Nodes, incoi, sphi(1:ncoi,ioa:ioa+nso(l_i)-1) = ', &
! &                               Node, Nodes, i1, i2, igto%sphi(i1,i2)
!          enddo 
!        enddo 
!        do i1 = 1, ncoj
!          do i2 = joa, joa+nso(l_j)-1
!            write(6,'(a,4i5,f12.5)')                       &
! &            'calc_prescreen_eri: Node, Nodes, incoj, sphi(1:ncoi,joa:joa+nso(l_j)-1) = ', &
! &                               Node, Nodes, i1, i2, jgto%sphi(i1,i2)
!          enddo 
!        enddo 
!!       End debugging
 
        call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_i), &
&         1, nso(l_j), name='eri', routine='calc_prescreen_eri')
        eri(:,:,:,:) = 0.0_dp

        call calc_contract_eri2(libint_data, hfx_sys%cell, hfx_sys%cell_r, &
&         ri, rj, ri, rj, npgfi, npgfj, npgfi, npgfj,                      &
&         l_i, l_j, l_i, l_j, ncoi, ncoj, ncoi, ncoj,                      &
&         igto%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)),                   &
&         jgto%orbnl_zeta(1:npgfj, jspp%orb_index(joa)),                   &
&         igto%orbnl_zeta(1:npgfi, ispp%orb_index(ioa)),                   &
&         jgto%orbnl_zeta(1:npgfj, jspp%orb_index(joa)),                   &
&         igto%sphi(1:ncoi,ioa:ioa+nso(l_i)-1),                            &
&         jgto%sphi(1:ncoj,joa:joa+nso(l_j)-1),                            &
&         igto%sphi(1:ncoi,ioa:ioa+nso(l_i)-1),                            &
&         jgto%sphi(1:ncoj,joa:joa+nso(l_j)-1),                            &
&         hfx_opts, eri )

!       We identiy now what is the maximum value of the electron repulsion 
!       integrals eri. 
!       These eris have been computed among all the orbitals (l,m) of the shell
!       to which io and jo orbitals belong to.
!       The factor of 2.0 is a conversion from Hartrees to Rydbergs
        eri_prescreen(ind) = 2.0_dp * max(maxval(eri), -minval(eri))

        call de_alloc(eri, name='eri', routine='calc_prescreen_eri')
      end do    ! End loop on the neighbours
    end do      ! End loop on the orbitals in the unit cell

!   We identiy now what is the maximum value of eri_prescreen...
    max_eri = max(maxval(eri_prescreen),-minval(eri_prescreen))
!   ... and take the square root, since in the condition after
!   Eq.(25) of Ref. Shang:2011 appears as \f$ sqrt{M_{3}} \f$
    max_eri = sqrt(max_eri)

    if ( Node .eq. 0 ) then
      write(*,'("calc_prescreen_eri:",1X,A,E20.8)') &
&       "max_eri_prescreen", max_eri
    end if

!!   For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!    call die()
!!   End debugging


  end subroutine calc_prescreen_eri

! *****************************************************************************
!> \brief Prepare all the data required to prescreen the four center electron
!! repulsion integrals that will be computed.
!!
!! First, the neighbour lists are globalized, so all the nodes know all the 
!! neighbours of all the atoms.
!!
!! Second, the eri_prescreen matrix is allocated and computed, calling 
!! the calc_prescreen_eri subroutine.
!!
!! Third, the logical matrix um_cut is allocated and computed.
!!
!! Fourth, we compute the parameters of the quadratic functions to fit
!! the two center integrals as a function of the distance following the 
!! recipe given in Appendix C of Ref. \cite Guidon:2009
!!
!! This subroutine is called only once per molecular dynamic or conjugate
!! gradient step from siesta_forces.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] hfx_sys: data structure storing system information
!!
!! \bug When running in parallel, all the nodes perform the same operations.
!!      Potential point of optimization
!!
!! \bug um_cut is a dense matrix with the dimension of the number of orbitals
!!      of the supercell square. This can be optimized.
!!
!! POTENTIAL POINT OF OPTIMIZATION: eri_prescreen is defined as a square
!! matrix of range equal to the number of atomic orbitals in the auxiliary
!! supercell
! *****************************************************************************
  subroutine init_prescreen_eri(libint_data, hfx_opts, hfx_sys)

    use alloc         , only: de_alloc      ! Deallocation routine
    use alloc         , only: re_alloc      ! Reallocation routine
    use atmfuncs      , only: lofio         ! Returns the angular momentum
                                            !   quantum number of a given
                                            !   atomic basis orbital
    use atmfuncs      , only: mofio         ! Returns the magnetic quantum
                                            !   number of a given
                                            !   atomic basis orbital
    use atm_types     , only: nspecies      ! Number of different atomic species
    use atm_types     , only: maxn_orbnl    ! Maximum number of nl orbitals 
                                            !   (not counting different 
                                            !   "m" copies)
    use atomlist      , only: indxuo        ! Index of equivalent orbital in 
                                            !   the unit "u" cell
    use parallel      , only: IOnode        ! Input/output node
    use parallel      , only: Node          ! Local node
    use parallel      , only: Nodes         ! Total number of nodes
    use parallelsubs  , only: WhichNodeOrb  ! Given the global orbital pointer,
                                            !   this routine returns the Node 
                                            !   number where this is stored.
    use parallelsubs  , only: GlobalToLocalOrb   
                                            ! Converts an orbital index in the 
                                            !   global frame to the local frame
                                            !   if the orbital is local to this
                                            !   node. Otherwise the pointer is
                                            !   return as zero and can therefore
                                            !   be used to test whether the
                                            !   orbital is local or not.
    use listsc_module , only: listsc        ! Translates a listed neighbour   
                                            !   from unit cell to supercell
    use nao2gto_types, only: maxn_contract ! Maximum number of 
                                            !   Gaussians that enter in
                                            !   the linear combination
                                            !   to expand a given NAO
    use nao2gto_data,   only: D2Sindx
    use nao2gto_data,   only: eri_prescreen ! (phi_mu phi_nu|phi_mu phi_nu).
                                            !   It appears in Eq. (25) of
                                            !   Ref. \cite Shang:2011.
                                            !   This integral is computed as
                                            !   in Eq. (14) of
                                            !   Ref. \cite Shang:2011.
!    use nao2gto_data,   only: um_cut        ! Dense matrix that determines  
                                            !   whether two atomic orbitals in 
                                            !   the supercell are connected or
                                            !   cut
                                            !   .true.  = they are cut 
                                            !             i.e. they are not 
                                            !             neighbours
                                            !   .false. = they are not cut 
                                            !             i.e. they are 
                                            !             neighbours
    use nao2gto_data,   only: pair_dist_radii_pgf
    use nao2gto_data,   only: sfc_pgf       ! Parametrization of the screening 
                                            !    functions that are upper bound 
                                            !    to the two center integrals 
                                            !    between primitive Gaussians
                                            !    required to estimate the 
                                            !    Cauchy-Schwarz inequality.
    use nao2gto_data,   only: sfc_shell     ! Parametrization of the screening 
                                            !    functions that are upper bound 
                                            !    to the two center integrals 
                                            !    between contracted Gaussians
                                            !    required to estimate the 
                                            !    Cauchy-Schwarz inequality.
    use nao2gto_data,   only: sfc_kind      ! Similar parametrization as 
                                            !    before to screen the atoms
    use nao2gto_data,   only: list_ij       ! List of pairs of neighbour 
                                            !    orbitals where the orbital i
                                            !    is located in the unit cell
                                            !    and orbital j in the supercell
    use nao2gto_data,   only: list_kl       ! List of pairs of neighbour 
                                            !    orbitals where both the 
                                            !    orbital i and j are located 
                                            !    in the supercell
    use nao2gto_data,   only: log10_eps_schwarz
                                            ! Logarithm of the Schwarz tolerance

    use nao2gto_libint
    use nao2gto_types , only: hfx_options_type, hfx_system_type
#ifdef MPI
    use fdf
    use mpi_siesta
    use nao2gto_parallel
#endif

    implicit none

!   Arguments
    type(Libint_t)        , intent(inout) :: libint_data
!                                            Libint data structure
    type(hfx_options_type), intent(in)    :: hfx_opts
!                                            Data structure storing Hartree-Fock
!                                            exchange options
    type(hfx_system_type) , intent(in)    :: hfx_sys
!                                            Data structure storing system 
!                                            information

!   hfx_sys%nsc(3)  ! Number of unit cells required to build the 
                    !   auxiliary supercell along a given direction
!   hfx_sys%nuotot  ! Total number of atomic orbitals in the unit cell
!   hfx_sys%norb    ! Total number of atomic orbitals in the auxiliary supercell

!   Local variables
    integer  :: ncells    ! Number of unit cells required to build
                          !   the auxiliary supercell
    integer  :: nelem_ij  ! Total number of different pairs that can be 
                          !   constructed between one orbitals in the unit cell
                          !   (let us call it mu in cell 0)
                          !   and any orbital in the supercell 
                          !   (let us call it nu in cell N)
    integer  :: nelem_kl  ! Total number of different pairs that can be
                          !   constructed between one orbitals in the supercell 
                          !   (let us call it lambda in cell G)
                          !   and any orbital in the supercell 
                          !   (let us call it sigma in cell H)
    integer  :: io        ! Counter for loop on atomic orbitals
    integer  :: is        ! Atomic species to which a given orbital belongs
    integer  :: ioa       ! Orbital index of each orbital io in its atom
    integer  :: l         ! Angular quantum number of the atomic orbital
    integer  :: m         ! Magnetic quantum number of the atomic orbital
    integer  :: iu        ! Equivalent orbital in the unit cell of orbital io
    integer  :: j         ! Counter for loop on neighbors
    integer  :: ind       ! Index of the neighbor
    integer  :: ju        ! Neighbor orbital of an orbital in the unit cell
    integer  :: jo        ! Translation of the neighbour from unit cell 
                          !   to supercell
    integer  :: js        ! Chemical species of the atom to which 
                          !   the neighbor orbital belongs to
    integer  :: joa       ! Equivalent index of the neighbor orbital 
                          !   in the unit cell
    integer  :: l_j       ! Angular quantum number of the neighbor orbital
    integer  :: m_j       ! Magnetic quantum number of the neighbor orbital
    integer  :: ishell    ! Counter to loop on shells
    integer  :: jshell    ! Counter to loop on shells
    integer  :: ipgf      ! Counter to loop on primitive Gaussian functions
    integer  :: jpgf      ! Counter to loop on primitive Gaussian functions
    real(dp) :: max_eri   ! Maximum value of eri_prescreen. 
                          !   It corresponds to the variable \f$ M_{e} \f$ 
                          !   defined after Eq. (25) in Ref. \cite Shang:2011
    integer :: maxnhg 
! Variables required to globalize the neighbour list
#ifdef MPI
    integer :: MPIerror, bsizedefault
    integer :: BNode      ! Returns the node number where an orbital is stored
!    integer :: maxnhg     ! Global maximum number of neighbors over all atoms
    integer :: iio        ! Local orbital index
    integer,  dimension(:), pointer :: numhg => null()
                          ! Globalize version of numh
    integer,  dimension(:), pointer :: listhptrg => null()
                          ! Globalize version of listhptr
    integer,  dimension(:), pointer :: listhg => null()
                          ! Globalize version of listh
#endif

!-------------------------------------------------------------------------------

! Globalize the neighbor list, so all the nodes node all the list of neighbors
! of all the atoms

#ifdef MPI
    call re_alloc(numhg, 1, hfx_sys%nuotot, name='numhg', &
&     routine='setup_hfx')
    call re_alloc(listhptrg, 1, hfx_sys%nuotot, name='listhptrg', &
&     routine='setup_hfx')

!   Globalize numh
    do io = 1, hfx_sys%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        numhg(io) = hfx_sys%numh(iio)
      endif
      call MPI_Bcast(numhg(io), 1, MPI_integer, BNode, &
&       MPI_Comm_World, MPIerror)
    enddo

!   Build global listhptr
    listhptrg(1) = 0
    do io = 2, hfx_sys%nuotot
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
    enddo

!   Globalise listh
    maxnhg = listhptrg(hfx_sys%nuotot) + numhg(hfx_sys%nuotot)
    call re_alloc(listhg, 1, maxnhg, name='listhg', routine='setup_hfx')

    do io = 1, hfx_sys%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node .eq. BNode ) then
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        do jo=1,numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
&           hfx_sys%listh(hfx_sys%listhptr(iio)+1: &
&             hfx_sys%listhptr(iio)+hfx_sys%numh(iio))
        enddo
      endif

      call MPI_Bcast(listhg(listhptrg(io)+1), numhg(io), MPI_integer,&
&       BNode, MPI_Comm_World, MPIerror)
    enddo
#endif

!   Init HFX cell parameters
!   Compute the number of unit cells in the supercell
    ncells   = hfx_sys%nsc(1) * hfx_sys%nsc(2) * hfx_sys%nsc(3)
    nelem_ij = ncells * hfx_sys%nuotot * (hfx_sys%nuotot + 1) / 2
    nelem_kl = ncells*ncells*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2

!!   For debugging
!    write(6,'(a,5i5)')                                      &
! &    'init_prescreen_eri: Node, Nodes, hfx_sys%nsc = ',    &
! &    Node, Nodes, hfx_sys%nsc(:)    
!    write(6,'(a,3i5)')                                      &
! &    'init_prescreen_eri: Node, Nodes, ncells = ',         &
! &    Node, Nodes, ncells
!    write(6,'(a,3i5)')                                      &
! &    'init_prescreen_eri: Node, Nodes, hfx_sys%nuotot = ', &
! &    Node, Nodes, hfx_sys%nuotot
!    write(6,'(a,3i10)')                                     &
! &    'init_prescreen_eri: Node, Nodes, nelem_ij = ',       &
! &    Node, Nodes, nelem_ij
!    write(6,'(a,3i10)')                                     &
! &    'init_prescreen_eri: Node, Nodes, nelem_kl = ',       &
! &    Node, Nodes, nelem_kl
!    write(6,'(a,3i5)')                                       &
! &    'init_prescreen_eri: Node, Nodes, hfx_sys%norb = ',   &
! &    Node, Nodes, hfx_sys%norb
!!   End debugging
!!! --------------------------------- xmqin 2021-----------
!   Allocate D2Sindx and initialize all the entries to 0

    call re_alloc(D2Sindx, 1, hfx_sys%norb, 1, hfx_sys%norb, &
 &     name='D2Sindx', routine='init_prescreen_eri')
    D2Sindx(:,:) = 0

!   Loop over all the atomic orbitals in the associated supercell
    do io = 1, hfx_sys%norb
       iu = indxuo(io)
#ifdef MPI
!     Loop on the neighbours of the equivalent orbital in the unit cell
      do j = 1, numhg(iu)
!       Indentify the index of the neighbor in the list
        ind = listhptrg(iu) + j
!       Identify the index of the neighbor orbital 
        ju = listhg(ind)
!       Translate the listed neighbour from unit cell to supercell
        jo = listsc(io,iu,ju)
        D2Sindx(io,jo) = ind
      end do
#else
!     Do the same as before but in serial
      do j = 1, hfx_sys%numh(iu)
        ind = hfx_sys%listhptr(iu) + j
        ju = hfx_sys%listh(ind)
        jo = listsc(io,iu,ju)
        D2Sindx(io,jo) = ind
      end do
#endif
    end do
!!!--------------------------------------------------------------------------
!   Prepare data for Hartree-Fock exchange
!   to calculate ERI prescreen matrix
!    if(node.eq.0) write(6,*) "maxnh: ", hfx_sys%maxnh, maxnhg
    call re_alloc(eri_prescreen, 1, maxnhg, &
&         name='eri_prescreen', routine='init_prescreen_eri')
    eri_prescreen(:) = 0.0_dp

#ifdef MPI
    !write(*,*) "MAXNHG ", maxnhg, Node
    !write(*,*) "NUMHG ", numhg, Node
    !write(*,*) "LISTHG ", associated(listhg), Node
    !write(*,*) "LISTHPTRG ", associated(listhptrg), Node
    call calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&     maxnhg, numhg, listhptrg, listhg, max_eri)
#else
    call calc_prescreen_eri(libint_data, hfx_opts, hfx_sys, &
&     hfx_sys%maxnh, hfx_sys%numh, hfx_sys%listhptr, hfx_sys%listh, max_eri)
#endif

!   Add screenfunc method (Xinming Qin)
    if ( associated(pair_dist_radii_pgf) ) then
      deallocate(pair_dist_radii_pgf)
      pair_dist_radii_pgf => null()
    end if
    allocate(pair_dist_radii_pgf(maxn_contract, maxn_contract, maxn_orbnl, &
&     maxn_orbnl, nspecies, nspecies))

    if ( associated(sfc_kind) ) then
      deallocate(sfc_kind)
      sfc_kind => null()
    end if
    allocate(sfc_kind(nspecies, nspecies))

    if ( associated(sfc_pgf) ) then
      deallocate(sfc_pgf)
      sfc_pgf => null()
    end if
    allocate(sfc_pgf(maxn_contract, maxn_contract, maxn_orbnl, maxn_orbnl, &
&     nspecies, nspecies))

    if ( associated(sfc_shell) ) then
      deallocate(sfc_shell)
      sfc_shell => null()
    end if
    allocate(sfc_shell(maxn_orbnl,maxn_orbnl,nspecies,nspecies))

!   Initialize the screening functions
    do is = 1, nspecies
      do js = 1, nspecies
        sfc_kind(js,is)%x(:) = 0.0_dp
        do ishell = 1, maxn_orbnl
          do jshell = 1, maxn_orbnl
            sfc_shell(jshell,ishell,js,is)%x(:) = 0.0_dp
            do ipgf = 1, maxn_contract
              do jpgf = 1, maxn_contract
                 pair_dist_radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(:) = &
&                  0.0_dp
                 sfc_pgf(jpgf,ipgf,jshell,ishell,js,is)%x(:) = 0.0_dp
              end do
            end do
          end do
        end do
      end do
    end do

!   Compute the parameters that fit the logarithmic dependencies of the
!   two center integrals at large distances to a quadratic function, 
!   according to Guidon et al.
    call calc_pair_dist_radii(hfx_opts, pair_dist_radii_pgf)
    call calc_screening_functions(libint_data, hfx_opts, hfx_sys%cell, &
&     hfx_sys%cell_r, sfc_pgf, sfc_shell)

!   Build left and right orbital pair_list
!     max_element1 = ncells*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2
!     max_element2 = ncells**2*hfx_sys%nuotot*(hfx_sys%nuotot+1)/2
!   Note: element is not a simple type, hence we cannot use re_alloc
    if ( associated(list_ij%element) ) then
      deallocate(list_ij%element)
    end if
    allocate(list_ij%element(nelem_ij))
    if ( associated(list_kl%element) ) then
      deallocate(list_kl%element)
    end if
    allocate(list_kl%element(nelem_kl))

    call build_pair_list(hfx_sys%nua, hfx_sys%na, hfx_opts, hfx_sys, &
&     eri_prescreen, max_eri, list_ij)

    call build_pair_list(hfx_sys%na, hfx_sys%na, hfx_opts, hfx_sys, &
&     eri_prescreen, max_eri, list_kl)

!   Take the logarithm of the Schwarz tolerance.
!   This variable is declared in nao2gto_data.F90, 
!   defined here, 
!   and used in the calculation of the four center integrals 
!   (file nao2gto_contract.F90)
    log10_eps_schwarz = log10(hfx_opts%eps_schwarz)

! FIXME: I/O should be done in nao2gto_io.F90
! FIXME: LibFDF will soon drop MPI support
!! For debugging
    if(node.eq.0)then
      write(*,'(A,1X,I12)')  'init_prescreen_eri: nelem_ij =', nelem_ij
      write(*,'(A,1X,I12)')  'init_prescreen_eri: nelem_kl =', nelem_kl
      write(*,'(A,1X,I12)')  'init_prescreen_eri: list_ij  :',list_ij%nelement
      write(*,'(A,1X,I12)')  'init_prescreen_eri: list_kl  :',list_kl%nelement
      write(*,'(A,1X,E20.6)') 'init_prescreen_eri: max_eri =', max_eri
!      write(*,'(A,1X,E20.6," (log =",E20.6,")")') &
    endif
!! End debugging

#ifdef MPI
!    if ( IOnode ) then
!      call set_bsize_NAO2GTO(Nodes, list_kl%nelement, bsizedefault)
!      nao2gto_bsize = fdf_integer('blocksize_NAO2GTO', bsizedefault)
!!     For debugging
!      write(*,'(A,1X,I16)')     &
! &      "init_prescreen_eri: Blocksize for ERIs distribution:", nao2gto_bsize
!!     End debugging
!    end if
!    call MPI_Bcast(nao2gto_bsize, 1, MPI_integer, 0, MPI_Comm_World, MPIerror)

    call de_alloc(listhg, name='listhg', routine='setup_hfx')
    call de_alloc(listhptrg, name='listhptrg', routine='setup_hfx')
    call de_alloc(numhg, name='numhg', routine='setup_hfx')
#endif

  end subroutine init_prescreen_eri

! *****************************************************************************
! *** Private routines                                                      ***
! *****************************************************************************

! ---------------------------------------------------------------------------
!  New subroutine to find shell pair list (uv) and (mn)
!
!  Coded by xmqin, 12. 06, 2016,  different from shanghui's method.
!
!  (1) Find neighbours for all atoms in supercell using rmax ,
!      Number of neighbours  : nnia    max : maxnna = 200
!      Neighbours' index: jna(nnia)
!      PBC vectors to neighbours : xij(jna)
!      Squared distances to neighbours : r2ij(jna),
!
!      Actually, I call the original subroutine "neighbour()" of SIESTA
!      program to find atomic pair_list without additional definition
!      for pair_atomic_list '(type)'.
!      Loop over atoms to give xij and rij in batch !
!
!  (2) Find all shell pairs over orbital loop
!     (a) Hermite and trans symmetry : ishell_unit > jshell_unit
!     (b) NAO overlap criteria: rcut(io)+rcut(jo) > rij, io overlap with jo,
!     (c) Schwarz inequality included density matrix
!---------------------------------------------------------------------------
!
! Useful information:
!
! (1)  To deal with periodic boundary conditions, the unit cell is 'extended'
!      on each side.
!      There are two kinds of index for images( unit cell, atoms and obitals ) :
!      one is in the PBC supercell, the other is in the normal (or auxiliary) supercell.
!
!  PBC unit cell : Wigner-Seitz cell A(i) = -1/2a(i), 1/2a(i) , a, Ai : lattice vector.
!  PBC supercell : To extend Wigner cell on each side of direction: -nsc(i)*A/2, nsc(i)*A/2
!  Normal supercell : Direct supercell : A(i) = 1, nsc(i)*a !
!
!  The index of images should be shift between pbc supercell and normal supercell
!  because we usually use normal index.
!
!  (1) In order to take into account the 8-fold symmetry under pbc :
!
!   (u0v[R]|m[R']n[R"])        =   (v0u[-R]|m[R'-R]n[R"-R])
! = (u0v[R]|n[R"]m[R'])        =   (v0u[-R]|n[R"-R]m[R'-R])
! = (m0n[R"-R']|u[-R']v[R-R']) =   (n0m[R'-R"]|u[-R"]v[R-R"])
! = (m0n[R"-R']|v[R-R']u[-R']) =   (n0m[R'-R"]|v[R-R"]u[-R"])
!
!  The first index must be always in unit cell, so we should transfer R to 0.
!
!   pair_list by iterating in the following way.
!   we should transfer R to 0
!
!   DO iatom=1,natom
!      iatom_unit = mod(iatom-1, unit_atoms)+1
!
!      DO jatom=iatom, negbours(iatom)
!         jatom_pbc = ind1 (ind2(jatom)-ind2(iatom)+ind2(iatom_unit))
!         jatom_unit_pbc = mod(jatom-1, unit_atoms)+1
!         if (jatom_unit<iatom_unit) CYCLE
!
!         atom_ij = ncells*(iatom_unit)*(iatom_unit-1)/2
!                  +mod(jatom_pbc,unit_atoms)*iatom_unit+jatom_unit_pbc
!
!  atom_ij is the uv or mn pair list index gives u[0]v[R] or
!  m[0]n[R],  it is not equal to number of loop.
!
!  Here "_pbc" denotes images' index in pbc supercell.
!  v[0] means that we choose v[R] as reference unitcell.
!  then u, m, n must be shift to extended supercell related to v[R].
!
!!
!! if (katom+latom<=iatom+jatom)  then
!!   if ( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
! or
!----------------------------------------------------------------------------
! integer nia: Atom index of orbital (u)s into the loop .
!              To build (uv), nia = nua  for u
!                       (mn), nia = na   for m

! integer norb         : Number of orbitals in supercell
! real*8  scell(3,3)   : Supercell vectors SCELL(IXYZ,IVECT)
! integer nsc(3)       : Num. of unit cells in each supercell direction
! real*8  xa(3,na)     : Atomic positions in cartesian coordinates
! integer lasto(0:na)  : Last orbital of each atom in array iphorb
! integer iphorb(norb) : Orbital index of each orbital in its atom,
!                       where no=lasto(na)
! integer iphorb(norb)
! integer isa(na)      : Species index of each atom

!---------------------------------------------------------------------------
! *****************************************************************************
!> \brief Identify the significant shell pairs 
!!        \f$ \sqrt{\left( \mu \nu \vert \mu \nu \right)_{\rm SR}} \f$,
!!        which satisfies 
!!        \f$ \sqrt{\left( \mu \nu \vert \mu \nu \right)_{\rm SR}} \times M_{e} 
!!        \ge {\rm hfx\_opts\%eps\_pairlist} \f$.
!!        Those pairs are preselected and stored.
!!
!! It is important to note that we store shell pairs.
!! If the pair \f$ (\mu \nu) \f$ is stored, then that means the set
!! over all the possible combinations of functions with the same angular 
!! momentum of \f$ \mu \f$, and the same angular momentum of \f$ \nu \f$ will
!! be computed.
!! For example, a \f$ (ps|sd) \f$ class (or quartet) consists of 
!! 3 (for the \f$ p \f$) \f$ \times \f$ 
!! 1 (for the \f$ s \f$) \f$ \times \f$ 
!! 1 (for the \f$ s \f$) \f$ \times \f$ 
!! 6 (for the \f$ d \f$, remember that LIBINT computes using the
!! Cartesian Gaussian functions) = 18 integrals.
!!
!! According to the Schwarz inequality
!! \f{eqnarray*}{
!! \vert \left( \mu \nu \vert \lambda \sigma \right) \vert
!! \le \sqrt{\left( \mu \nu \vert \mu \nu \right) 
!!           \left( \lambda \sigma \vert \lambda \sigma \right)}
!! \le \sqrt{\left( \mu \nu \vert \mu \nu \right)} \times M_{e},
!! \f}
!! where \f$ M_{e} = \max \sqrt{\left( \mu \nu \vert \mu \nu \right)} \f$.
!!
!! Only the pair of shells that satisfy 
!! \f$ \sqrt{\left( \mu \nu \vert \mu \nu \right)} \times M_{e} \ge 
!!     {\rm hfx\_opts\%eps\_pairlist} \f$
!! will be retained. Here \f$ {\rm hfx\_opts\%eps\_pairlist} \f$ is equal to
!! tolerance1 in the Ref. \cite Shang:2011
!! Of course, in order to have a non-vanishing value of the integral,
!! orbitals \f$ \mu \f$ and \f$ \nu \f$ should overlap,
!! so \f$ r_{c}(\mu) + r_{c} (\nu) > r_{\mu \nu} \f$,
!! where \f$ r_{c} \f$ are the cutoff radii of the atomic orbitals,
!! and \f$ r_{\mu \nu} \f$ is the distance between the two centers.
!!
!! Only pairs where \f$l_{\mu} \ge l_{\nu} \f$ are stored, where
!! \f$ l_{\mu} \f$ stands for the angular momentum of the shell to which
!! atomic orbital \f$ \mu \f$ belongs. 
!! This is related with the fact that,
!! currently LIBINT has a very important restriction on the 
!! angular momentum ordering of the functions in shell quartets that 
!! it can handle. LIBINT can evaluate a shell quartet
!! \f$ (\mu \nu|\lambda \sigma) \f$ if 
!! \f$ l_{\mu} \ge l_{\nu}, l_{\lambda} \ge l_{\sigma} \f$, 
!! and \f$ l_{\lambda} + l_{\sigma} \ge l_{\mu} + l_{\nu} \f$. 
!! If one needs to compute a quartet that doesn’t conform the rule, 
!! e.g. of type \f$ (pf|sd) \f$, permutational symmetry of integrals can be 
!! utilized to compute such quartet:
!!
!! \f$ (\mu \nu | \lambda \sigma) = (\mu \nu |\sigma \lambda) = 
!! (\nu \mu| \lambda \sigma) = (\nu \mu|\sigma \lambda) = 
!! ( \lambda \sigma|\mu \nu) = (\lambda \sigma|\nu \mu) = 
!! (\sigma \lambda |\mu \nu) = (\sigma \lambda |\nu \mu) \f$. 
!!
!! In the case of \f$ (pf|sd) \f$ shell quartet, 
!! one computes quartet \f$ (ds|fp) \f$ instead, and then permutes
!! function indices back to obtain the desired
!! \f$ (pf|sd) \f$.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \param[in] nia: depending on the call, number of atoms in the unit cell
!!                 or in the supercell
!! \param[in] na:  number of atoms in the supercell
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] hfx_sys: data structure storing system information
!! \param[in] max_eri: maximum value of eri_prescreen (\f$ M_{e} \f$ in
!!                     the previous equations.
!! \param[out] list_orb: list of angular momentum shells neighbors
! *****************************************************************************
  subroutine build_pair_list(nia, na, hfx_opts, hfx_sys, &
                 eri_prescreen, max_eri, list_orb)

    use alloc         , only: re_alloc      ! Reallocation routines
    use alloc         , only: de_alloc      ! Deallocation routines
    use atmfuncs      , only: lofio         ! Returns the angular momentum 
                                            !   quantum number of a given
                                            !   atomic basis orbital
    use atmfuncs      , only: mofio         ! Returns the magnetic quantum 
                                            !   number of a given
                                            !   atomic basis orbital
    use atmfuncs      , only: rcut          ! Returns cutoff radius of a given
                                            !   atomic basis orbital
    use atomlist      , only: indxuo        ! Index of equivalent orbital in 
                                            !   the unit cell
    use atomlist      , only: lasto         ! Position of last orbital of 
                                            !   each atom
    use atomlist      , only: rmaxo         ! Maximum cutoff radius between 
                                            !   all the orbitals
    use neighbour     , only: jna=>jan      ! Atom index of neighbours
    use neighbour     , only: xij           ! Vectors from a given atom to 
                                            !   neighbours
    use neighbour     , only: r2ij          ! Square of the distance between 
                                            !   an atom and its neighbours
    use neighbour     , only: mneighb       ! Finds the neighbours of an atom
                                            !   in a cell with periodic boundary
                                            !   conditions
    use neighbour     , only: maxna=>maxnna ! Maximum number of neighbour orbits
    use sorting       , only: iorder        ! Orders elements of an integer 
                                            !   array according to some index
    use sorting       , only: order         ! Orders elements of a real array 
                                            !   according to some index
    use sorting       , only: ordix         ! Finds an order index of 
                                            !   increasing array elements
    use precision     , only: dp            ! Double precision
    use nao2gto_data  , only: subshell
    use nao2gto_data  , only: D2Sindx
    use nao2gto_index , only: indexsc
    use nao2gto_pbc   , only: trans_pbc
    use nao2gto_types , only: hfx_options_type, hfx_system_type, pair_list_type
!   For debugging
#ifdef MPI
    use mpi_siesta
#endif
!   End debugging

    implicit none

!   Arguments
    integer               , intent(in)  :: nia    ! Depending on the argument
                                                  !    in the call, it might be
                                                  !    the number of atoms in 
                                                  !    the unit cell or in the
                                                  !    supercell.
                                                  !    Known by all the nodes.
    integer               , intent(in)  :: na     ! Number of atoms in the 
                                                  !    supercell.
                                                  !    Known by all the nodes.
    type(hfx_options_type), intent(in)  :: hfx_opts
                                                  ! Data structure storing 
                                                  !    Hartree-Fock exchange  
                                                  !    options
    type(hfx_system_type) , intent(in)  :: hfx_sys 
                                                  ! Data structure storing 
                                                  !    system information
     
    real(dp)              , intent(in)  :: &
&     eri_prescreen(hfx_sys%maxnh)    !(phi_mu phi_nu|phi_mu phi_nu)
                                                  !   It appears in Eq. (25) of
                                                  !   Ref. \cite Shang:2011.
                                                  !   This integral is computed
                                                  !   as in Eq. (14) of
                                                  !   Ref. \cite Shang:2011.
    real(dp)              , intent(in)  :: max_eri
                                                  ! Maximum value of 
                                                  !   eri_prescreen
    type(pair_list_type)  , intent(out) :: list_orb
                                                  ! List of angular momentum 
                                                  !   shells neighbors

    ! Local variables
    integer  :: ncells                            ! Number of unit cells in the
                                                  !   supercell
    integer  :: n_shells
    integer  :: nnia                              ! Number of neighbours of atom
                                                  !   ia
    integer  :: is                                ! Species index of orbital io
    integer  :: js                                ! Species index of orbital jo
    integer  :: i, ja, jn, jo, jos
    integer  :: ia                                ! Counter for loop on atoms
    integer  :: io                                ! Counter for loop on orbitals
    integer  :: iuo                               ! Equivalent of orbital of io
                                                  !   in the unit cell
    integer  :: juo                               ! Equivalent of orbital of jo
                                                  !   in the unit cell
    integer  :: ioa                               ! Orbital index of each 
                                                  !   orbital io in its atom
    integer  :: joa                               ! Orbital index of each 
                                                  !   orbital jo in its atom
    integer  :: l_i                               ! Angular momentum of the 
                                                  !   orbital whose neighbours
                                                  !   are sought
    integer  :: m_i                               ! Magnetic quantum number of 
                                                  !   the orbital whose 
                                                  !   neighbours are sought
    integer  :: l_j                               ! Angular momentum of the 
                                                  !   neighbour orbital 
    integer  :: m_j                               ! Magnetic quantum number of 
                                                  !   the neighbour orbital 
    integer  :: ishell                            ! Sequential index of the 
                                                  !   angular momentum of the
                                                  !   orbital whose neighbour
                                                  !   are sought
                                                  !   in the supercell
                                                  !   1 <= ishell
                                                  !   ishell <= unit_shell*ncell
    integer  :: ishell_unit                       ! Equivalent to ishell but
                                                  !   in the unit cell
                                                  !   1 <= ishell_unit
                                                  !   ishell_unit <= unit_shell
    integer  :: jshell                            ! Sequential index of the 
                                                  !   angular momentum of the
                                                  !   neighbour orbital
                                                  !   in the supercell
                                                  !   1 <= jshell
                                                  !   jshell <= unit_shell*ncell
    integer  :: jshell_unit                       ! Equivalent to jshell but
                                                  !   in the unit cell
                                                  !   1 <= jshell_unit
                                                  !   jshell_unit <= unit_shell
    integer  :: unit_shells                       ! Total number of angular 
                                                  !   momentum shells 
                                                  !   in the unit cell,
                                                  !   i.e. if we have a 
                                                  !   unit cell with two atoms
                                                  !   and each atoms has one
                                                  !   - s shell
                                                  !   - p shell
                                                  !   then unit_shell = 4.
    real(dp) :: ri(3)                             ! Position of the atom
                                                  !   whose neighbours are 
                                                  !   sought
    real(dp) :: rj(3)                             ! Position of the neighbour 
                                                  !   atom
    real(dp) :: r_temp(3)                         ! Relative position between 
                                                  !   a given atom and 
                                                  !   its neighbours,
                                                  !   computed from the atomic
                                                  !   positions in the array XA
                                                  !   Remember that, 
                                                  !   due to the PBC,
                                                  !   a folding might have 
                                                  !   to be applied:
                                                  !   the actual neighbour 
                                                  !   might be a 
                                                  !   periodic replica of the 
                                                  !   atom included
                                                  !   in the neighbour list, 
                                                  !   located in 
                                                  !   other supercell
    real(dp) :: r_pbc(3)                          ! Actual relative position 
                                                  !   between a given atom and 
                                                  !   its neighbours, after
                                                  !   performing the folding
    real(dp) :: r2                                ! Square of the distance 
                                                  !   between a neighbour
                                                  !   and the atom
    real(dp) :: rij                               ! Distance between a neighbour
     
    integer, dimension(:), pointer :: index => null()
                                                  ! Array used to order the 
                                                  !   elements of the neighbour
                                                  !   list

!   For debugging
    integer :: MPIerror
!   End debugging

! -------------------------------------------------------------------------

!!   For debugging
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, nia = ', &
! &    Node, Nodes, nia
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, na = ', &
! &    Node, Nodes, na
!    write(6,'(a,2i5,f12.5)') &
! &    'build_pair_list: Node, Nodes, max_eri = ', &
! &    Node, Nodes, max_eri
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, maxna before = ', &
! &    Node, Nodes, maxna
!!   End debugging

!   Allocate local memory
!   hfx_sys%nuotot is the number of atomic orbitals in the unit cell
!   hfx_sys%norb   is the number of atomic orbitals in the supercell
!   hfx_sys%na     is the number of atoms in the supercell
!   hfx_sys%cell   are the supercell lattice vectors
!   hfx_sys%xa     are the atomic positions of all the atoms in the supercell
!   All these quantities are known by all the nodes

!   Compute the number of different angular momentum shells in the unit cell
    unit_shells = subshell(hfx_sys%nuotot)
!   Compute the number of unit cells in the supercell
    ncells = hfx_sys%norb / hfx_sys%nuotot
!   Initialize neighb subroutine
!    call mneighb(hfx_sys%cell, 2.0_dp*rmaxo, na, hfx_sys%xa, 0, 0, nnia )
!    call re_alloc(index, 1, maxna, name='index', routine='build_pair_list')

!!   For debugging
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, hfx_sys_nuotot = ', &
! &    Node, Nodes, hfx_sys%nuotot
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, hfx_sys_norb = ', &
! &    Node, Nodes, hfx_sys%norb
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, hfx_sys_na   = ', &
! &    Node, Nodes, hfx_sys%na
!    write(6,'(a,2i5,f12.5)') &
! &    'build_pair_list: Node, Nodes, rmaxo        = ', &
! &    Node, Nodes, rmaxo
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, unit_shells = ', &
! &    Node, Nodes, unit_shells
!    write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, ncells = ', &
! &    Node, Nodes, ncells
!        write(6,'(a,3i5)') &
! &    'build_pair_list: Node, Nodes, maxna after = ', &
! &    Node, Nodes, maxna
!!   End debugging

    list_orb%nelement = 0

!   Loop over all the atoms in the unit cell or in the supercell
!   depending on the value of nia in the call
    do ia = 1, nia 

!     Localize where this atom is centered
      ri(1:3) = hfx_sys%xa(1:3,ia)

!     Find all the neighbour atoms within a range of 
!     2.0 * maximum cutoff radii
!      call mneighb(hfx_sys%cell, 2.0_dp*rmaxo, na, hfx_sys%xa, ia, 0, nnia)

!!     For debugging
!      write(6,'(a,2i5)') &
! &      'build_pair_list: ia, nnia = ', & 
! &      ia, nnia
!!     End debugging

!     Sort the neighbours by distance
!     The neighbour list will be ordered in such a way that they will appear
!     in increasing order of the distance
!      call ordix (r2ij, 1, nnia, index)
!      call iorder(jna,  1, nnia, index)
!      call order (r2ij, 1, nnia, index)
!      call order (xij,  3, nnia, index)

!     Loop over all the atomic orbitals within the atom
      do io = lasto(ia-1) + 1, lasto(ia)

!       Identify the orbital index of each orbital in its atom
        ioa = hfx_sys%iphorb(io)
!       Identify the atomic index of the orbital 
        is  = hfx_sys%isa(hfx_sys%iaorb(io))
!       Identify the orbital angular momentum
        l_i = lofio(is,ioa)
!       Identify the magnetic angular momentum
        m_i = mofio(is,ioa)

!       We follow only the first time that a given angular momentum 
!       of a given atom appears
!       In other words the loop will continue only for different angular 
!       momentum shells
        if ( m_i .ne. -l_i ) cycle

!       Identify the sequential index of the angular momentum shell in the
!       supercell
        ishell      = subshell(io)

!       Identify the equivalent orbital in the unit cell
        iuo         = mod(io-1, hfx_sys%nuotot) + 1

!       Identify the sequential index of the angular momentum shell in the
!       supercell for the equivalent atom in the unit cell
        ishell_unit = subshell(iuo)

!!       For debugging
!        if( Node .eq. 0 ) then
!          write(6,'(a,9i7,f12.5)')                                              &
! &          'build_pair_list: ia, io, is, ioa, l_i, m_i, iuo, ishell, ishell_unit = ', & 
! &          ia, io, is, ioa, l_i, m_i, iuo, iuo, ishell_unit, rcut(is,ioa)
!        endif
!!       End debugging

!       Loop over all the neighbours
!        do jn = 1, nnia

!         Identify the atomic index of the neighbour
!          ja = jna(jn)
        do ja = 1, na
!         Identify the relative position between the neighbour and the atom
!         whose neigbours are sought
          r_temp(1:3) = hfx_sys%xa(1:3,ja) - hfx_sys%xa(1:3,ia)

!         Compute the actual relative position between the neighbour and
!         the atom, after performing the folding
!         hfx_sys%cell are the supercell lattice vectors in real space
          call trans_pbc(r_temp, hfx_sys%cell, hfx_sys%cell_r, r_pbc)

!         Square of the distance between a neighbour and the atom
          r2 = (r_pbc(1)**2) + (r_pbc(2)**2) + (r_pbc(3)**2)

!         Position of the neighbour atom
          rj = ri + r_pbc

!         Distance between a neighbour and the atom
          rij = sqrt(r2)

!!         For debugging
!          if( Node .eq. 0 ) then
!            write(6,'(a,2i7)') &
! &            '   build_pair_list: jn, ja = ', & 
! &                                 jn, ja
!            write(6,'(a,i5,3f12.5,i5,3f12.5)')   &
! &            '   build_pair_list: ia, xaia, ja, xaja = ', & 
! &            ia, hfx_sys%xa(1:3,ia), ja, hfx_sys%xa(1:3,ja)
!            write(6,'(a,2i5,7f12.5)')   &
! &            '   build_pair_list: ia, ja, r_temp, r_pbc, rij = ', & 
! &            ia, ja, r_temp, r_pbc, rij
!          endif
!!         End debugging

          do jo = lasto(ja-1)+1, lasto(ja)
!           Identify the orbital index of each orbital in its atom
            joa = hfx_sys%iphorb(jo)
!           Identify the atomic index of the orbital 
            js  = hfx_sys%isa(hfx_sys%iaorb(jo))
!           Identify the orbital angular momentum
            l_j = lofio(js,joa)
!           Identify the magnetic angular momentum
            m_j = mofio(js,joa)

!           We follow only the first time that a given angular momentum 
!           of a given atom appears
!           In other words the loop will continue only for different angular 
!           momentum shells
            if ( m_j .ne. -l_j ) cycle

!           Identify index of the neighbour orbital in the supercell
!           (jjunquer: check. might be useless)
            jos = indexsc(io, iuo, jo)

!           Identify equivalent index of the neighbour orbital in the unit cell
            juo = mod(jos-1, hfx_sys%nuotot) + 1

!           Identify the sequential index of the angular momentum shell in the
!           supercell for the neighbour orbital in the supercell
            jshell = subshell(jos)

!           Identify the sequential index of the angular momentum shell in the
!           supercell for the equivalent neighbour in the unit cell
            jshell_unit = subshell(juo)

!!           For debugging
!            if( Node .eq. 0 ) then
!              write(6,'(a,7i7,f12.5)') &
! &              '      build_pair_list: jo, jos, juo, m_j, l_j, jshell, jshell_unit, rcut = ', & 
! &                                      jo, jos, juo, m_j, l_j, jshell, jshell_unit, rcut(js,joa) 
!              write(6,'(a,3f20.12)') &
! &              '      build_pair_list: eri_prescreen(io,jo), sqrt(abs(eri_prescreen(io,jo)))*max_eri, hfx_opts%eps_pairlist = ', &
! &                     eri_prescreen(io,jo), sqrt(abs(eri_prescreen(io,jo)))*max_eri, hfx_opts%eps_pairlist 
!            endif
!!           End debugging

!           LIBINT can evaluate a shell quartet (ab|cd) if
!           \lambda(a) \ge \lambda(b).
!           Here, we compare the two angular momentum shells
            if ( jshell_unit .le. ishell_unit ) then
              if ( (rcut(is,ioa) + rcut(js,joa)) .le. rij ) cycle

!             This is the condition written two lines below Eq. (25)
!             of Ref. \cite Shang:2011, with tolerance_1 in the paper
!             begin equal to hfx_opts%eps_pairlist
!             According to the Schwarz inequality, we have
!             ( \mu \nu | \lambda \sigma ) <= 
!                      sqrt ( \mu \nu | \mu \nu ) * 
!                      sqrt ( \lambda \sigma | \lambda \sigma ) <=  
!                      sqrt ( \mu \nu | \mu \nu ) * max_preERI,
!                      where  
!                      max_preERI = max [ sqrt( \lambda\sigma | \lambda\sigma) ]
!                      is the maximum matrix element of two-center integrals.
!             We shall consider only pair of shells that pass this filter
              if ( sqrt(abs(eri_prescreen(D2Sindx(io,jo))))*max_eri .ge. &
                   hfx_opts%eps_pairlist ) then
                n_shells = ncells*(ishell_unit)*(ishell_unit-1)/2 &
&                        + ((jshell-1)/unit_shells)*ishell_unit + jshell_unit

                list_orb%nelement = list_orb%nelement + 1
                i = list_orb%nelement
                list_orb%element(i)%pair(1) = io
                list_orb%element(i)%pair(2) = jo
                list_orb%element(i)%nl_index = n_shells
                list_orb%element(i)%r1(1:3) = ri(1:3)
                list_orb%element(i)%r2(1:3) = rj(1:3)
                list_orb%element(i)%dist2 = r2
!!               For debugging
!                if (Node .eq. 0) then
!                write(6,'(a,i7)') &
! &                '         build_pair_list: list_orb%nelement = ', &
! &                list_orb%nelement
!                write(6,'(a,i7)') &
! &                '         build_pair_list: list_orb%element(i)%pair(1) = ', &
! &                list_orb%element(i)%pair(1)
!                write(6,'(a,i7)') &
! &                '         build_pair_list: list_orb%element(i)%pair(2) = ', &
! &                list_orb%element(i)%pair(2)
!                write(6,'(a,i7)') &
! &                '         build_pair_list: list_orb%element(i)%nl_index = ', &
! &                list_orb%element(i)%nl_index
!                write(6,'(a,3f12.5)') &
! &                '         build_pair_list: list_orb%element(i)%r1 = ', &
! &                list_orb%element(i)%r1
!                write(6,'(a,3f12.5)') &
! &                '         build_pair_list: list_orb%element(i)%r2 = ', &
! &                list_orb%element(i)%r2
!                write(6,'(a,f12.5)') &
! &                '         build_pair_list: list_orb%element(i)%dist2 = ', &
! &                list_orb%element(i)%dist2
!                write(6,'(a,3i5,2f12.5)') &
! &                '         build_pair_list: Node, ishell_unit, jshell_unit, rcut(is,ioa), rcut(js,joa) = ', & 
! &                Node, ishell_unit, jshell_unit, rcut(is,ioa), rcut(js,joa)
!                write(6,'(a,i5,3f12.5)') &
! &                '         build_pair_list: Node, eri_prescreen, max_eri, eps_pairlist = ', &
! &                Node, eri_prescreen(io,jo), max_eri, hfx_opts%eps_pairlist
!                write(6,'(a,3i5)') &
! &                '         build_pair_list: Node, n_shells, list_orb%nelement = ', &
! &                Node, n_shells, list_orb%nelement 
!                endif 
!!               End debugging
              end if
            end if
          end do   ! jo
        end do   ! jn
      end do   ! io
    end do   ! ia

!    call de_alloc(index, name='index', routine='build_pair_list')
!!   For debugging
!#ifdef MPI
!    call MPI_barrier(MPI_Comm_world,MPIerror)
!#endif
!    call die()
!!   End debugging

  end subroutine build_pair_list


! ***************************************************************************
!> \brief Fits the radii of the density, \f$ R_{\mu \nu}^{\rho}\f$,
!!        coming from the product of two Gaussians \f$ (\mu \nu^{\vec{a}}) \f$
!!        as a function of the distance between its centers of two Gaussians
!!        to a quadratic function.
!!        Explained in Appendix C of Ref. \cite Guidon:2009
!!
!!        Used for far-field screening
! ***************************************************************************
  subroutine calc_pair_dist_radii(hfx_opts, radii_pgf)

    use atm_types,     only : species_info
                             ! Information about the chemical species
    use atm_types,     only : species
                             ! Information about the chemical species
    use atm_types,     only : nspecies
                             ! Number of different chemical species
    use nao2gto_types, only : gto_info_type
                             ! Information about Gaussian type orbitals
    use nao2gto_types, only : hfx_options_type
                             ! Data type to store Hartree-Fock exchange options
                             !     that can be read from SIESTA input files
    use nao2gto_types, only : hfx_screen_coeff_type
                             ! Data type to store information about 
                             !    screening coefficients
    use nao2gto_data,  only : hfx_gtos
                             ! Information about Gaussian type orbitals
    use nao2gto_data,  only : pair_dist_radii_pgf
                             !  Fit the of the precomputed values
                             !     of the center of the product
                             !     densities as a function of the
                             !     inter Gaussian distance to a
                             !     two-parameter function of the
                             !     kind x(1)*rab^2 + x(2)
                             !     Described in Appendix 5C of
                             !     Ref. \cite Guidon:2009
    use nao2gto_utils, only : exp_radius_very_extended
                             ! Function that computes the radii of the product
                             !    densities of two Gaussians

    implicit none

!   Arguments
    type(hfx_options_type), intent(in) :: hfx_opts
    type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer :: radii_pgf

!   Local variables
    integer :: i             ! Counter for loop on orbitals of a given shell
                             !   and counter to loop on the grid of distances
                             !   between the center of the Gaussians
    integer :: ipgf          ! Counter for loop on Gaussians to expand the
                             !   radial part of a NAO for specie is
    integer :: jpgf          ! Counter for loop on Gaussians to expand the
                             !   radial part of a NAO for specie js
    integer :: is            ! Counter for loop on atomic species
    integer :: js            ! Counter for loop on atomic species
    integer :: io            ! Counter for loop on NAO for specie is
    integer :: jo            ! Counter for loop on NAO for specie js
    integer :: ishell        ! Number of different shells of species is
    integer :: jshell        ! Number of different shells of species js
    integer :: l_i           ! Angular momentum of the orbitals of specie is
    integer :: m_i           ! Magnetic quantum number of the orbitals
                             !    of specie is
    integer :: l_j           ! Angular momentum of the orbitals of specie js
    integer :: m_j           ! Magnetic quantum number of the orbitals
                             !    of specie js
    integer :: npgfi         ! Number of Gaussian functions required to expand
                             !    the radial part of a given shell of NAOs
                             !    for the orbital of species is
    integer :: npgfj         ! Number of Gaussian functions required to expand
                             !    the radial part of a given shell of NAOs
                             !    for the orbital of species js
    real(dp) :: radius       ! Sum of the radius of two Gaussians
                             !    (the radius of a Gaussian, is like its 
                             !    cutoff radius, and is defined as the point
                             !    where g(r)=coeff * r**l * exp(-zeta*r**2)
                             !    is smaller than a given threshold,
                             !    1.d-5 in nao2gto_transfer
                             !    subroutine
    real(dp) :: ra(3)        ! Position of the first Gaussian
                             !    It will be always assumed to be at the origin
    real(dp) :: rb(3)        ! Position of the second Gaussian
    real(dp) :: rab(3)       ! Relative position between the two Gaussians
    real(dp) :: rp(3)        ! Position P where the product of the two Gaussians
                             !    will be centered.
                             !    This is defined in Eq. (2.33) of 
                             !    Ref. \cite Fermann:2020
    real(dp) :: rap(3)       ! Relative position where the product of the
                             !    two Gaussians will be centered with respect
                             !    to Gaussian a.
                             !    This is defined in Eq. (2.33) of 
                             !    Ref. \cite Fermann:2020
    real(dp) :: rab2         ! Square of the distance between the two Gaussians
    real(dp) :: zetp         ! Sum of the exponents of the two Gaussians
    real(dp) :: ff           ! \alpha_{2} / (\alpha_{1} + \alpha_{2}),
                             ! required to define the factor K below
    real(dp) :: prefactor    ! Factor K defined in Eq. (2.37) in 
                             !    Ref. \cite Fermann:2020
    real(dp) :: table(2,0:100) 
                             ! First index:
                             !    table(1,:): distance between 
                             !    the center of the two Gaussians
                             !    table(2,:): radii of the product densities
                             !    coming from the two Gaussians
                             ! Second index:
                             !    the previous two quantities are computed over
                             !    a suitably grid of 100 points
    real(dp) :: cutoff
    real(dp) :: r1, r_max    ! Center of the product densities 
    real(dp) :: x(2)         ! Fit the of the precomputed values of the center
                             !    of the product densities as a function of 
                             !    the inter Gaussian     
                             !    distance to a two-parameter function of the 
                             !    kind x(1)*rab^2 + x(2)
    type(species_info), pointer  :: ispp => null()
                             ! Pointer to the variable containing the
                             !   information on the atomic species is
    type(species_info), pointer  :: jspp => null()
                             ! Pointer to the variable containing the
                             !   information on the atomic species js
    type(gto_info_type), pointer :: igto => null()
                             ! Pointer to the variable containing the
                             !   information on the Gaussian type orbitals
                             !   for the atomic species is
    type(gto_info_type), pointer :: jgto => null()
                             ! Pointer to the variable containing the
                             !   information on the Gaussian type orbitals
                             !   for the atomic species js
! ------------------------------------------------------------------------------

    real(dp) :: eps_new

    if( hfx_opts%eps_schwarz >1.0E-6_dp ) then
      eps_new = 1.0E-6_dp
    else
      eps_new = hfx_opts%eps_schwarz
    endif
!   Initialize some variables.
!   One of the primitive Gaussians (ra) will be always located at the origin
    ra = 0.0_dp
!   The other primitive Gaussian will change its position, but here it is
!   initialized to zero
    rb = 0.0_dp
    table = 0.0_dp

    do is = 1, nspecies
      ispp => species(is)
      igto => hfx_gtos(is)
      do js = 1, nspecies
        jspp => species(js)
        jgto => hfx_gtos(js)

        ishell = 0
        do io = 1, ispp%norbs
          l_i = ispp%orb_l(io)
          m_i = ispp%orb_m(io)
!         The loop only runs for the different shells.
!         i.e. only the first time that a given orbital in a shell appears
!         If not, it will cycle
          if ( m_i .ne. -l_i ) cycle     ! Not a normal orbital

!         Update the number of different shells
          ishell = ishell+1

!         Identify the number of Gaussian functions required to expand
!         the radial part of a given shell of NAO
          npgfi = igto%orbnl_contract(ishell)

          jshell = 0
          do jo = 1, jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j .ne. -l_j ) cycle     ! Not a normal orbital
            jshell = jshell+1
            npgfj = jgto%orbnl_contract(jshell)

!           Loop over all the Gaussian functions required to expand
!           the radial part of a given NAOs of the atom of species is
            do ipgf = 1, npgfi

!             Loop over all the Gaussian functions required to expand
!             the radial part of a given NAOs of the atom of species is
              do jpgf = 1, npgfj

!               Compute the sum of the "cutoff radius" of the two
!               primitive Gaussians
                radius = igto%pgf_radius(ipgf,ishell) + &
&                        jgto%pgf_radius(jpgf,jshell)

!               Loop over a grid of suitable distances between the center
!               of the two Gaussians.
!               Here, the initial point is when the two Gaussians are centered
!               at the origin, and the last point the sum of the two 
!               "cutoffs" of the Gaussian.
!               In between, one hundred points uniformly distributed
                do i = 0, 100
!                 Displace the second Gaussian in the x-direction,
!                 so its distance will progressively increase from 0
!                 (i.e. the same as the center of the first Gaussian)
!                 till radius 
!                 (i.e. the distances where the two Gaussians will not overlap)
                  rb(1) = 0.0_dp + 0.01_dp * radius * i

                  r_max = 0.0_dp

!                 Compute the sum of the exponents of the two Gaussians
                  zetp = igto%orbnl_zeta(ipgf,ishell) + &
 &                       jgto%orbnl_zeta(jpgf,jshell)

!                 Compute the product of the exponents of the second
!                 Gaussian with the inverse of the sum of the two exponents
!                 This is part of the exponent related with Eq. (2.37) of the
!                 notes by \cite Fermann:2020
!                 at the time of describing the Gaussian Product Theorem,
!                 where 
!                 \gamma = \alpha_1 + \alpha_2, being the alphas the exponents,
!                 and the prefactor, (K in the notes), is defined as
!                 K = \exp^{-(\alpha_{1} \alpha_{2} \overline{\vec{AB}}^{2} / 
!                             \gamma} 
!                 Since the first Gaussian is always at the origin,
!                 then \vec{A} = 0, \overline{\vec{AB}}^{2} = rab2
                  ff = jgto%orbnl_zeta(jpgf,jshell)/zetp

!                 Compute the relative positions between the two Gaussians
!                 Since the Gaussian A is always at the origin, this 
!                 equals the position of the second Gaussian
                  rab    = 0.0_dp
                  rab(1) = rb(1)

!                 Compute the square of the distance between the two Gaussians
                  rab2 = rb(1)**2

!                 Compute the prefactor K defined in Eq. (2.37) of Ref.
!                 \cite Fermann:2020
                  prefactor = exp(-igto%orbnl_zeta(ipgf,ishell)*ff*rab2)

!                 Compute the position P where the product of the two Gaussians
!                 will be centered.
!                 This is defined in Eq. (2.33) of Ref. \cite Fermann:2020
!                 Remember that the position of the first Gaussian is the origin
!                 so the term \alpha_{1} \vec{A} in the numerator cancels
                  rap(:) = ff*rab(:)

!                 Center of the Gaussian that results from the product of the
!                 two Gaussians
                  rp(:) = ra(:) + rap(:)

!                 Center of the second Gaussian (b).
!                 The first one (a) is assumed to be at the origin.
                  rb(:) = ra(:) + rab(:)

                  cutoff = 1.0_dp

!                 exp_radius_very_extended is a subroutine borrowed 
!                 from CP2K. According to Joost VandeVondele comments:
!                 but it is looking for the radius of Gaussian functions, 
!                 taking the product of two Gaussian with cartesian prefactors, !                 taking the radial part (and optionally density matrix) 
!                 also into account. 
!                 So, basically, first finding the term in a block 
!                 which contributes most to r^l * exp(-a *x^2), 
!                 and after that take into account. 
!                 So, find all the terms in the product and 
!                 compute for each of the terms their radius, 
!                 taking the largest radius found.
                  r1 = exp_radius_very_extended(l_i, l_i, l_j, l_j, &
&                      ra=ra, rb=rb, rp=rp, zetp=zetp, eps=eps_new, &
&                      prefactor=prefactor, cutoff=cutoff,epsabs=1.0d-12)
                  r_max = max(r_max,r1)

                  table(1,i) = rb(1)
                  table(2,i) = r_max

!!                 For debugging
!                  write(6,'(a,7i5,2f12.5)') &
! &                  'calc_pair_dist_radii: jpgf,ipgf,jshell,ishell,js,is,i,table = ', &
! &                  jpgf, ipgf, jshell, ishell, js, is, i, table(:,i) 
!!                 End debugging
                enddo

!               Fit the of the precomputed values of the center
!               of the product densities as a function of the inter Gaussian 
!               distance to a two-parameter function of the kind
!               x(1)*rab^2 + x(2)
                call run_powell_optimize(table,x,0.0_dp)

!               Store the two-parameters from the fitting in an array 
!               that depend on the two Gaussian functions 
!               (atom kind, shell, sets, etc).
                radii_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x

!!               For debugging
!                write(6,'(a,6i5,2f12.5)') &
! &                'calc_pair_dist_radii: jpgf,ipgf,jshell,ishell,js,is,x = ', &
! &                jpgf, ipgf, jshell, ishell, js, is, x 
!!               End debugging
              enddo ! End loop on jpgf (loop over all Gaussians required for jo)
            enddo   ! End loop on ipgf (loop over all Gaussians required for io)
          enddo     ! End loop on jo (loop over all the atomic orbitals of js)
        enddo       ! End loop on io (loop over all the atomic orbitals of is)
      enddo         ! End loop on js (loop over all species)
    enddo           ! End loop on is (loop over all species)

  end subroutine calc_pair_dist_radii

! ***************************************************************************
!> \brief Parametrization of the screening functions that are upper bound to 
!!        the two center integrals required to estimate the Cauchy-Schwarz 
!!        inequality.
!!        Explained in Appendix C of Ref. \cite Guidon:2009
!! 
!! In this Appendix, we can read how 
!! for the near-field screening one can rely on the Cauchy-Schwarz 
!! inequality
!!
!! \f{eqnarray*}{
!!   \left| \left( \mu \nu^{\vec{a}} \vert 
!!         \lambda^{\vec{b}} \sigma^{\vec{b}+\vec{c}} \right) \right|
!!         \le 
!!   \left| \left( \mu \nu^{\vec{a}} \vert \mu \nu^{\vec{a}}\right)\right|^{1/2}
!!   \times 
!!   \left| \left( \lambda \sigma^{\vec{c}} \vert \lambda \sigma^{\vec{c}} 
!!         \right) \right|^{1/2}.
!! \f}
!!  
!! This only requires two center integrals.
!! However, instead of storing or computing these two center integrals,
!! it is very efficient to instead parametrize screening functions that are
!! an upper bound to these integrals.
!! These screening functions only depend parametrically on the interatomic
!! distance \f$ \vert \vec{R}_{\mu \nu} \vert  \f$ but are different
!! for each type of Gaussian basis function (atom kind, shell, sets).
!! 
!! These fits can be easily performed once one observes that the logarithm
!! of the integral is similar to a quadratic function at larger distance
!!
!! \f{eqnarray*}{
!!   \log \left[ \left( \mu \nu \vert \mu \nu \right) (R_{\mu \nu}) \right]
!!   \approx a_{2} R_{\mu \nu}^{2} + a_{0}.
!! \f}
!!
!! This choice leads to the useful properties that the estimate 
!! decays monotonically with increasing distance,
!! and that the expression only requires the square distance between 
!! the centers. 
!! The coefficients \f$a_{0}\f$ and \f$a_{2}\f$ are calculated once and 
!! for all at the beginning of a calculation.
!! This is precisely the goal of this subroutine.
!!
!! The most important difference with respect the implementation
!! carried out in the former paper is that within the SIESTA subroutine
!! we fit the logarithm of the square root of the two center integral,
!! and not the logarithm of the two center integral as in CP2K.
!!
!! \param[in,out] libint_data: Libint data structure
!! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
!! \param[in] cell: Supercell lattice vectors in real space
!! \param[in] rcell: Reciprocal lattice of the supercell in real space
!! \param[out] coeffs_pgf: Coefficients \f$a_{0}\f$ and \f$a_{2}\f$ 
!!                        to fit the two center
!!                        integrals between primitive Gaussian functions to
!!                        a quadratic function.
!!                        The fitting coefficients are of the form
!!                        \f$ \left[ {\rm coeffs\_pgf(1)} \times R_{\mu \nu}^{2} +
!!                        {\rm coeffs\_pgf(2)} \right]  \f$
!! \param[out] coeffs_set: Coefficients \f$a_{0}\f$ and \f$a_{2}\f$ 
!!                        to fit the two center
!!                        integrals between contracted Gaussian functions to
!!                        a quadratic function.
!!                        The fitting coefficients are of the form
!!                        \f$ \left[ {\rm coeffs\_set(1)} \times R_{\mu \nu}^{2} +
!!                        {\rm coeffs\_set(2)} \right]  \f$
!!
!! \bug As it is implemented in this moment, all the nodes run the same code.
!! This can be highly improved for efficiency.
! ***************************************************************************
  subroutine calc_screening_functions( libint_data, hfx_opts, cell, rcell,  &
&                                      coeffs_pgf, coeffs_set )

    use alloc    ,     only : re_alloc
                             ! Reallocation routine
    use alloc    ,     only : de_alloc
                             ! Deallocation routine
    use atm_types,     only : species_info
                             ! Information about the chemical species
    use atm_types,     only : species
                             ! Information about the chemical species
    use atm_types,     only : nspecies
                             ! Number of different chemical species
    use nao2gto_data,  only : hfx_gtos
                             ! Information about Gaussian type orbitals
    use nao2gto_data,  only : nco
                             ! number of Cartesian Gaussians for a given 
                             !    angular momentum
    use nao2gto_data,  only : nso
                             ! number of Spherical Harmonic Gaussians for a
                             !    given angular momentum
    use nao2gto_types, only : gto_info_type
                             ! Information about Gaussian type orbitals
    use nao2gto_types, only : hfx_options_type
                             ! Data type to store Hartree-Fock exchange options
                             !     that can be read from SIESTA input files
    use nao2gto_types, only : hfx_screen_coeff_type
                             ! Data type to store information about
                             !    screening coefficients
    use nao2gto_types, only : powell_min_log
    use nao2gto_data,  only : pair_dist_radii_pgf
                             !  Fit the of the precomputed values
                             !     of the center of the product
                             !     densities as a function of the
                             !     inter Gaussian distance to a
                             !     two-parameter function of the
                             !     kind x(1)*rab^2 + x(2)
                             !     Described in Appendix 5C of
                             !     Ref. \cite Guidon:2009
    use nao2gto_contract
    use nao2gto_libint
    use nao2gto_primitive, only: calc_primitive_screen

    implicit none

!   Arguments
    type(Libint_t), intent(inout) :: libint_data
    type(hfx_options_type), intent(in) :: hfx_opts

    real(dp), intent(in)  :: cell(3,3), rcell(3,3)
    type(hfx_screen_coeff_type), dimension(:,:,:,:), pointer :: coeffs_set
    type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer :: coeffs_pgf

!   Local variables
    integer :: i             ! Counter for loop on orbitals of a given shell
                             !   and counter to loop on the grid of distances
                             !   between the center of the Gaussians
    integer :: ipgf          ! Counter for loop on Gaussians to expand the
                             !   radial part of a NAO for specie is
    integer :: jpgf          ! Counter for loop on Gaussians to expand the
                             !   radial part of a NAO for specie js
    integer :: is            ! Counter for loop on atomic species
    integer :: js            ! Counter for loop on atomic species
    integer :: io            ! Counter for loop on NAO for specie is
    integer :: jo            ! Counter for loop on NAO for specie js
    integer :: ishell        ! Number of different shells of species is
    integer :: jshell        ! Number of different shells of species js
    integer :: l_i           ! Angular momentum of the orbitals of specie is
    integer :: m_i           ! Magnetic quantum number of the orbitals 
                             !    of specie is
    integer :: l_j           ! Angular momentum of the orbitals of specie js
    integer :: m_j           ! Magnetic quantum number of the orbitals 
                             !    of specie js
    integer :: npgfi         ! Number of Gaussian functions required to expand
                             !    the radial part of a given shell of NAOs
                             !    for the orbital of species is
    integer :: npgfj         ! Number of Gaussian functions required to expand
                             !    the radial part of a given shell of NAOs
                             !    for the orbital of species js
    integer :: ncoi          ! Number of primitive Cartesian Gaussian functions 
                             !    required to expand a shell of species is
    integer :: ncoj          ! Number of primitive Cartesian Gaussian functions 
                             !    required to expand a shell of species js
    real(dp) :: zeta_i       ! Exponents of the Gaussians required to expand
                             !    the radial part of the NAOs for species is
    real(dp) :: zeta_j       ! Exponents of the Gaussians required to expand
                             !    the radial part of the NAOs for species js
    real(dp) :: radius       ! Sum of the radius of two Gaussians
                             !    (the radius of a Gaussian plays the role 
                             !    of its cutoff radius, and is defined as the 
                             !    point where g(r)=coeff* r**l * exp(-zeta*r**2)
                             !    is smaller than a given threshold,
                             !    1.d-5 in nao2gto_transfer
                             !    subroutine
    real(dp) :: ra(3)        ! Position of the first primitive Gaussian
                             !    It will be always assumed to be at the origin
    real(dp) :: rb(3)        ! Position of the second primitive Gaussian
    real(dp) :: table(2,0:100)
                             ! Table with:
                             !    the interatomic distance (table(1,:))
                             !    the log_10 ((\mu \nu \vert \mu \nu)) for
                             !    that particular distance. 
                             !    This is tabulated for 100 points 
    real(dp) :: x(2)         ! Fit the of the precomputed values of the two
                             !    center integrals as a function of
                             !    the inter Gaussian
                             !    distance to a two-parameter function of the
                             !    kind x(1)*rab^2 + x(2)

    real(dp) :: max_contraction_a, max_contraction_b, max_val, max_val_temp
    real(dp), dimension(:,:,:,:), pointer :: eri => null()
    type(species_info), pointer  :: ispp => null()
                             ! Pointer to the variable containing the 
                             !   information on the atomic species is
    type(species_info), pointer  :: jspp => null()
                             ! Pointer to the variable containing the 
                             !   information on the atomic species js
    type(gto_info_type), pointer :: igto => null()
                             ! Pointer to the variable containing the 
                             !   information on the Gaussian type orbitals 
                             !   for the atomic species is
    type(gto_info_type), pointer :: jgto => null()
                             ! Pointer to the variable containing the 
                             !   information on the Gaussian type orbitals 
                             !   for the atomic species js

! ------------------------------------------------------------------------------

!   Initialize some variables
!   One of the primitive Gaussians (ra) will be always located at the origin
    ra = 0.0_dp
!   The other primitive Gaussian will change its position, but here it is
!   initialized to zero
    rb = 0.0_dp

    table = 0.0_dp

!   Loop over the first atomic species
    do is = 1, nspecies
      ispp => species(is)
      igto => hfx_gtos(is)

!     Loop over the second atomic species
      do js = 1, nspecies
        jspp => species(js)
        jgto => hfx_gtos(js)

        ishell = 0
!       Loop over all the atomic orbitals of the first atomic specie
        do io = 1, ispp%norbs
          l_i = ispp%orb_l(io)
          m_i = ispp%orb_m(io)
!         The loop only runs for the different shells.
!         i.e. only the first time that a given orbital in a shell appears
!         If not, it will cycle
          if ( m_i /= -l_i ) cycle     ! Not a normal orbital

!         Update the number of different shells
          ishell = ishell+1

!         Identify the number of Gaussian functions required to expand
!         the radial part of a given shell of NAO
          npgfi = igto%orbnl_contract(ishell)

!         Identify the total number of primitive Cartesian Gaussian functions
!         to expand a shell 
          ncoi  = nco(l_i)*npgfi

          max_contraction_a = maxval((/(sum(abs(igto%sphi(1:ncoi,i))), &
 &            i=io,io+nso(l_i)-1)/))

!!         For debugging
!          write(6,'(a,6i5,f12.5)') &
! &          'calc_screening_functions: is, io, l_i, m_i, npgfi, ncoi, max_contraction_a = ', &  
! &          is, io, l_i, m_i, npgfi, ncoi, max_contraction_a 
!!         End debugging

          jshell = 0
!         Loop over all the atomic orbitals of the second atomic specie
          do jo = 1,jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j /= -l_j ) cycle     ! Not a normal orbital
            jshell = jshell+1
            npgfj = jgto%orbnl_contract(jshell)
            ncoj  = nco(l_j)*npgfj
            max_contraction_b =  maxval((/(sum(abs(jgto%sphi(1:ncoj,i))), &
 &                i=jo,jo+nso(l_j)-1)/))

!           Loop over all the Gaussian functions required to expand 
!           the radial part of a given NAOs of the atom of species is
            do ipgf = 1, npgfi

!             Exponent of the corresponding Gaussian
              zeta_i = igto%orbnl_zeta(ipgf,ishell)

!             Loop over all the Gaussian functions required to expand 
!             the radial part of a given NAOs of the atom of species is
              do jpgf = 1, npgfj
                zeta_j = jgto%orbnl_zeta(jpgf,jshell)
!               Compute the sum of the "cutoff radius" of the two
!               primitive Gaussians
                radius = igto%pgf_radius(ipgf,ishell) +  &
 &                       jgto%pgf_radius(jpgf,jshell)

!               Loop over a grid of suitable distances between the center
!               of the two Gaussians.
!               Here, the initial point is when the two Gaussians are centered
!               at the origin, and the last point is the sum of the two 
!               "cutoffs" of the Gaussian.
!               In between, one hundred points uniformly distributed
                do i = 0, 100
!                 Displace the second Gaussian in the x-direction,
!                 so its distance will progressively increase from 0
!                 (i.e. the same as the center of the first Gaussian)
!                 till radius
!                 (i.e. the distances where the two Gaussians will not overlap)
                  rb(1) = 0.0_dp + real(i,dp) * 0.01_dp * radius

                  max_val = 0.0_dp
                  max_val_temp = 0.0_dp

!                 Compute the four center integral between the Gaussians
!                 a and b, as a function of the distances between their centers
!                 (\mu \nu^{\vec{a}} | \mu \nu^{\vec{a}}) 
!                 in the notation of the Appendix C
!                 of Ref. \cite Guidon:2009
                  call calc_primitive_screen(libint_data, ra, rb, ra, rb, &
&                   zeta_i, zeta_j, zeta_i, zeta_j, l_i, l_j, l_i ,l_j, &
&                   max_val_temp, hfx_opts)

!                 Since the argument of the logarithm has to be a positive
!                 number, here we keep only the two center integrals if
!                 they are positive. If not, we set their values to zero.
                  max_val = max(max_val, max_val_temp)

!                 We take the square root of the two center integrals.
!                 Those are the values that finally enter the 
!                 Cauchy-Schwarz inequality
                  max_val = dsqrt(max_val)

                  max_val = max_val * max_contraction_a * max_contraction_b
                  table(1,i) = rb(1)
                  if ( max_val == 0.0_dp ) then
                    table(2,i) = powell_min_log
                  else
                    table(2,i) = log10((max_val))
                  end if
                end do

!               We fit the logarithm of the integral as a function of the 
!               interatomic distance with a quadratic function
!                
                call run_powell_optimize(table,x,powell_min_log)
                coeffs_pgf(jpgf,ipgf,jshell,ishell,js,is)%x = x

              enddo ! End loop on jpgf (loop over all Gaussians required for jo)
            enddo   ! End loop on ipgf (loop over all Gaussians required for io)
          enddo     ! End loop on jo (loop over all the atomic orbitals of js)
        enddo       ! End loop on io (loop over all the atomic orbitals of is)
      enddo         ! End loop on js (loop over all species)
    enddo           ! End loop on is (loop over all species)

    ra = 0.0_dp
    rb = 0.0_dp
    table = 0.0_dp

    do is=1,nspecies
      ispp => species(is)
      igto => hfx_gtos(is)
      do js = 1,nspecies
        jspp => species(js)
        jgto => hfx_gtos(js)

        ishell = 0
        do io = 1,ispp%norbs
           l_i = ispp%orb_l(io)
           m_i = ispp%orb_m(io)
           if (m_i .ne. -l_i) cycle     ! Not a normal orbital
           ishell = ishell+1
           npgfi = igto%orbnl_contract(ishell)
           ncoi  = nco(l_i)*npgfi

           max_contraction_a = maxval((/(sum(abs(igto%sphi(1:ncoi,i))), &
&            i=io,io+nso(l_i)-1)/))

          jshell = 0
          do jo = 1,jspp%norbs
            l_j = jspp%orb_l(jo)
            m_j = jspp%orb_m(jo)
            if ( m_j .ne. -l_j ) cycle     ! Not a normal orbital
            jshell = jshell+1
            npgfj = jgto%orbnl_contract(jshell)
            ncoj  = nco(l_j)*npgfj
            max_contraction_b = maxval((/(sum(abs(jgto%sphi(1:ncoj,i))), &
&             i=jo,jo+nso(l_j)-1)/))
            radius = igto%shell_radius(ishell)+jgto%shell_radius(jshell)

            do i=0,100
              rb(1) = 0.0_dp + real(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              max_val_temp = 0.0_dp

              call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_i), &
&               1, nso(l_j), name='eri', routine='calc_screening_functions')
              eri = 0.0_dp
              call calc_contract_eri2(libint_data, cell, rcell, ra, rb, ra, rb,&
                npgfi, npgfj, npgfi, npgfj, l_i, l_j, l_i, l_j, &
                ncoi, ncoj, ncoi, ncoj, igto%orbnl_zeta(1:npgfi, ishell), &
                jgto%orbnl_zeta(1:npgfj, jshell), &
                igto%orbnl_zeta(1:npgfi,ishell),  &
                jgto%orbnl_zeta(1:npgfj, jshell), &
                igto%sphi(1:ncoi,io:io+nso(l_i)-1), &
                jgto%sphi(1:ncoj,jo:jo+nso(l_j)-1), &
                igto%sphi(1:ncoi,io:io+nso(l_i)-1), &
                jgto%sphi(1:ncoj,jo:jo+nso(l_j)-1), &
                hfx_opts, eri)

              max_val = max(max_val, maxval(dabs(eri)))
              max_val = 2.0_dp*max_val
              max_val = dsqrt(max_val)

              table(1,i) = rb(1)
              if ( max_val == 0.0_dp ) then
                table(2,i) = powell_min_log
              else
                table(2,i) = log10((max_val))
              end if

              call de_alloc(eri, name='eri', routine='calc_screening_functions')

            end do
            call run_powell_optimize(table,x,powell_min_log)
            coeffs_set(jshell,ishell,js,is)%x = x

          enddo     ! End loop on jo (loop over all the atomic orbitals of js)
        enddo       ! End loop on io (loop over all the atomic orbitals of is)
      enddo         ! End loop on js (loop over all species)
    enddo           ! End loop on is (loop over all species)

  end subroutine calc_screening_functions

  ! ***************************************************************************
  !> \brief Driver for the Powell minimizer
  !!
  !! This routine drives the powell minimizer for the data found in a table
  !! storing function values. It constructs an approximate upper bound of
  !! the fitted function.
  !!
  !! This is done in two steps: first, we compute the symmetric weight
  !! to get a good unique initial guess, then we restart for the asymmetric
  !! one. The asymmetric function appears to have several local minima.
  !! Depending on the data to fit the loop over k can make the switch
  !! gradual, but there is seemingly not much need for it.
  !!
  !! \author Xinming Qin
  !!
  !! \par History
  !!      10.2018 Created [Xinming Qin]
  !!      11.2018 Refactored and further documented [Yann Pouillon]
  !!
  !! \param[in] table: data to fit
  !! \param[out] x: fitting coefficients of the form (x(1)*table(1)**2+x(2))
  !! \param[in] fmin: only values of table(2) larger than fmin are considered
  ! ***************************************************************************
  subroutine run_powell_optimize(table, x, fmin)

    use nao2gto_powell

    implicit none

    ! Arguments
    real(dp), intent(in)  :: table(2,0:100)
    real(dp), intent(out) :: x(2)
    real(dp), intent(in)  :: fmin

    ! Local variables
    integer :: i, k
    real(dp) :: f, large_weight, small_weight, weight
    type(opt_state_type) :: opt_state

    ! -------------------------------------------------------------------------

    ! Initial values
    x(1) = 0.0_dp
    x(2) = 0.0_dp

    do k=0,4,4

      small_weight=1.0_dp
      large_weight=small_weight*(10.0_dp**k)

      ! Init optimization run
      opt_state%state = 0
      opt_state%nvar = 2
      opt_state%iprint = 3
      opt_state%unit = 6
      opt_state%maxfun = 100000
      opt_state%rhobeg = 0.1_dp
      opt_state%rhoend = 0.000001_dp

      do

        ! Compute function
        if ( opt_state%state == 2 ) then
          opt_state%f = 0.0_dp
          do i=0,100
            f = x(1)*table(1,i)**2 +  x(2)
            if ( f > table(2,i) ) then
              weight = small_weight
            else
              weight = large_weight
            end if
            if ( table(2,i) > fmin ) then
              opt_state%f = opt_state%f + weight * (f-table(2,i))**2
            end if
          end do
        end if

        if ( opt_state%state == -1 ) exit

        call powell_optimize(opt_state%nvar, x, opt_state)
      end do

      ! Clean-up the mess
      opt_state%state = 8
      call powell_optimize(opt_state%nvar, x, opt_state)

     end do

  end subroutine run_powell_optimize

end module nao2gto_prescreen
