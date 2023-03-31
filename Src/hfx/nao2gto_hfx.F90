! *** Module: nao2gto_hfx ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Hartree-Fock exchange term of the Hamiltonian
!!
!! In the SIESTA method, the density matrix is stored using the same sparse
!! pattern as the Hamiltonian. This treatment is robust in pure DFT where the
!! effective pure-DFT potential is only determined by the electron density
!! and its gradients.
!!
!! \f$\rho(r) = \sum_{u,v}[Duv*\psi_u(r-R_u)*\psi_v(r-R_v)]\f$,
!! \f$V_{eff}[\rho(r)] /= 0\f$ only when
!! \f$\psi_u(r-R_u)\f$ overlaps with \f$\psi_v(r-R_v)\f$.
!!
!! To calculate \f$rho\f$ and \f$V_{eff}\f$, we only need to store a sparse
!! subset of \f$D_{uv}\f$ that u overlaps with v even though the actual
!! \f$D_{uv}\f$ is unknown before SCF.
!!
!! However, HFX is dependent with the density matrix of real space:
!! \f$\rho(r,r^') = \sum_{u,v} \left[ D_{uv}*\psi_u(r-R_u)*\psi_v(r^'-R_v)
!! \right]\f$
!! \f$\rho(r,r^') \approx exp[-a(r-r^'-R_u-R_v)]\f$, a property of
!! \f$\sqrt(E_g)\f$ in semiconductors.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 11.2016 Edited [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!
!! \bug Possible errors from sparse Hamiltonian and density matrix!
! *****************************************************************************
module nao2gto_hfx

  use precision

  implicit none

  private

  public :: setup_hfx

contains

  ! ***************************************************************************
  !> \brief Computes the Hartree-Fock exchange part of the Hamiltonian
  !!
  !! Parallel update HFX(u,v) to Hamiltonian,
  !! u: Local PAOs of unitcell per node,
  !! v: Global PAOs of supercell
  !!
  !! HFX(u,v) = (um|vn) * Dm (m,n)
  !! HFX(u,v) 2D structure ---> Hmat(ind) 1D sparse structure
  !!
  !! Hmat(pure DFT)<--- HFX (this subroutine)----> HFX (HSE06)
  !!
  !! \par History
  !!      - 10.2013 Created [Xinming Qin]
  !!      - 11.2016 Edited [Xinming Qin]
  !!      - 01.2018 Replaced subroutine arguments by direct use of SIESTA
  !!        variables [Yann Pouillon]
  !!
  !! \param[in,out] libint_data: initialized Libint data structure
  !! \param[in] hfx_optdata: initialized Hartree-Fock options data structure
  !! \param[in] hfx_sysdata: initialized Hartree-Fock system data structure
  !! \param[in.out] Hmat: initialized hamiltonian where to add the
  !!                      HFX contribution
  ! ***************************************************************************
  subroutine setup_hfx(libint_data, hfx_optdata, hfx_sysdata, Hmat)

    use parallel,        only: IOnode            ! Input/Output node
    use parallel,        only: Node              ! Current node
    use parallel,        only: Nodes             ! Total number of nodes
    use parallelsubs,    only: GlobalToLocalOrb  ! Converts an orbital index
                                                 !    in the global frame to the
                                                 !    local frame if the orbital
                                                 !    is local to this node.
                                                 !    Otherwise the pointer is
                                                 !    return as zero and can
                                                 !    therefore be used to test
                                                 !    whether the orbital is
                                                 !    local or not.
    use parallelsubs,    only: LocalToGlobalOrb  ! Converts an orbital index in
                                                 !    the local frame to the
                                                 !    global frame
    use parallelsubs,    only: WhichNodeOrb      ! Given the global orbital
                                                 !    pointer, this routine
                                                 !    returns the node number
                                                 !    where this is stored.
    use alloc,           only: re_alloc          ! Reallocation routine
    use alloc,           only: de_alloc          ! Deallocation routine
    use units,           only: eV                ! Conversion factor
                                                 !    from Ry to eV
#ifdef MPI
    use mpi_siesta
#endif

    use m_energies,      only: Exc
    use sparse_matrices, only: Dscf, maxnh

    use nao2gto_data,    only: eri_prescreen     ! (mu nu | mu nu), following
                                                 !    Sec. III E 1 of
                                                 !    Ref. \cite Shang:2011
    use nao2gto_data,    only: D2Sindx
    use nao2gto_data,    only: hfx_call_counter  ! Number of times that the 
                                                 !   subroutine setup_hfx 
                                                 !   is called
    use nao2gto_dm
    use nao2gto_libint,  only: Libint_t
    use nao2gto_types ,  only: hfx_options_type
    use nao2gto_types ,  only: hfx_system_type

    implicit none

    ! Arguments
    type(Libint_t), intent(inout)      :: libint_data
    type(hfx_options_type), intent(in) :: hfx_optdata
    type(hfx_system_type), intent(in)  :: hfx_sysdata
    real(dp), intent(inout) :: Hmat(hfx_sysdata%maxnh,hfx_sysdata%nspin)

    ! Local variables
    integer  :: io       ! Counter for loop on orbitals in the unit cell
    integer  :: ispin    ! Counter for loop on spins
    integer  :: jo       ! Counter for loop on neighbours
    integer  :: j        ! Counter for loop on neighbours
    integer  :: ind      ! Index of the neighbour atom
    integer  :: num_u    ! Counter for loop on atomic orbitals at the time
                         !   of computing the exchange energy
    integer  :: num_v    ! Counter for loop on atomic orbitals at the time
                         !   of computing the exchange energy
    real(dp) :: E_HFX    ! Exchange energy, added to the exchange and 
                         !   correlation energy used in SIESTA
    real(dp) :: time_start, time_end, time_dm_start, time_dm_end
                         ! Variables to control the timing
    real(dp) :: spin_factor
                         ! Factor that determines the spin degeneracy
                         !   2.0 if a non-polarized spin calculation is done
                         !   1.0 if a polarized spin calculation is done

!    integer, dimension(:,:), pointer  ::  D2Sindx => null()

!    real(dp), dimension(:,:,:), pointer :: H_EXX => null()
    real(dp), dimension(:,:), pointer :: HFX => null()
                         ! Exact exchange contribution to 
                         !   the matrix elements
                         !   This corresponds to Eq.(9) of Ref.\cite Shang:2011
                         !   First  index: number of orbitals in the unit cell
                         !   Second index: number of orbitals in the supercell
                         !   Third  index: spin
!    real(dp), dimension(:,:,:), pointer :: DM_tmp => null()
!    real(dp), dimension(:,:), pointer :: Dscf => null()
                         ! Dense density matrix
                         !   This corresponds to Eq.(10) of Ref.\cite Shang:2011
                         !   Here we apply the condition written right below
                         !   Eq. (45) of the SIESTA paper
                         !   J. M. Soler et al., 
                         !   J. Phys.: Condens. Matter 14, 2745 (2002).
                         !   \rho_{\mu \nu} \equiv \rho_{\mu^\prime \nu^\prime}
                         !   if (\nu \mu) are equivalent by translation to 
                         !   (\nu^\prime \mu^\prime)
                         !   First  index: number of orbitals in the supercell
                         !   Second index: number of orbitals in the supercell
                         !   Third  index: spin
                         !   When running in parallel this is node independent
    real(dp), dimension(:),   pointer :: Dscf_max => null()
                         !   Allocate the matrix to store the maximum value  
                         !      of the density matrix between 
                         !      angular momentum shells
                         !   First  index: number of orbitals in the supercell
                         !   Second index: number of orbitals in the supercell
#ifdef MPI
    ! Global buffers for the storage of the sparse matrix
    ! FIXME: The arrays have the save attribute in the original version,
    !        check that pointers are OK
    integer :: BNode     ! Node number where the information 
                         !   (neighbour list, density matrix),
                         !   about a given orbital is stored
    integer :: iio       ! Local index of the orbital within the node
    integer :: MPIerror  ! Type of error in the MPI subroutines
    integer :: maxnhg    ! Number of total non-zero elements (sum of numh array)
    integer :: maxndg    ! Same as maxnhg, but for the density matrix
                         !    instead of the Hamiltonian
    integer :: indg
    integer,  dimension(:),     pointer :: numhg => null()
                         ! Global version of numh
                         !    Number of nonzero elements of each row 
                         !    of sparse matrices.
                         !    All nodes have access to the full neighbour list
    integer,  dimension(:),     pointer :: listhptrg => null()
                         ! Global version of listhptr
                         !    Pointer to start of row in listh
                         !    All nodes have access to the full neighbour list
    integer,  dimension(:),     pointer :: listhg => null()
                         ! Global version of listh
                         !    List of nonzero elements of each row of 
                         !    sparse matrices
                         !    All nodes have access to the full neighbour list
!   Same as before, but for the density matrix instead of the Hamiltonian
    integer,  dimension(:),     pointer :: numdg => null()
                         ! Global version of numd
                         !    Number of nonzero elements of each row 
                         !    of sparse matrices.
                         !    All nodes have access to the full neighbour list
    integer,  dimension(:),     pointer :: listdptrg => null()
                         ! Global version of listdptr
                         !    Pointer to start of row in listd
                         !    All nodes have access to the full neighbour list
    integer,  dimension(:),     pointer :: listdg => null()
                         ! Global version of listd
                         !    List of nonzero elements of each row of 
                         !    sparse matrices
                         !    All nodes have access to the full neighbour list
    real(dp), dimension(:,:),   pointer ::  Dscfg => null()
                         ! Global version of the Density Matrix
    real(dp), dimension(:,:), pointer ::  HFXg => null()
                         ! Part of part H_EXX computed in the local node
#endif

    external :: timer

!-------------------------------------------------------------------------------

    call timer('HFX', 1)

!   Counting calls to this routine for debugging
    hfx_call_counter = hfx_call_counter + 1
    if ( Node .eq.  0 ) then
      write(6,'("setup_hfx: call #",i4.4)') hfx_call_counter
    end if

!   Globalize the neighbour list, so all the nodes do have access to it
#ifdef MPI
    call re_alloc(numhg,     1, hfx_sysdata%nuotot, name='numhg',     &
      routine='setup_hfx')
    call re_alloc(listhptrg, 1, hfx_sysdata%nuotot, name='listhptrg', &
      routine='setup_hfx')
    call re_alloc(numdg,     1, hfx_sysdata%nuotot, name='numdg',     &
      routine='setup_hfx')
    call re_alloc(listdptrg, 1, hfx_sysdata%nuotot, name='listdptrg', &
      routine='setup_hfx')

!   Globalise numh
!   Loop over all the orbitals in the unit cell.
!   All the nodes see the same number of orbitals hfx_sysdata%nuotot
    do io = 1, hfx_sysdata%nuotot
!     Identify the node where the information on this orbital is stored
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
!       Identify the local index of the orbital within the node
        call GlobalToLocalOrb(io, Node, Nodes, iio)
!       Copy the information stored in the node into the global variable
        numhg(io) = hfx_sysdata%numh(iio)
        numdg(io) = hfx_sysdata%numh(iio)
      endif
!     Broadcast the information to the rest of the nodes
      call MPI_Bcast(numhg(io), 1, MPI_integer, BNode, &
        MPI_Comm_World, MPIerror)
      call MPI_Bcast(numdg(io), 1, MPI_integer, BNode, &
        MPI_Comm_World, MPIerror)
    enddo

!   Build global listhptr
    listhptrg(1) = 0
    listdptrg(1) = 0
    do io = 2, hfx_sysdata%nuotot
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      listdptrg(io) = listdptrg(io-1) + numdg(io-1)
    enddo

!   Globalise listh
    maxnhg = listhptrg(hfx_sysdata%nuotot) + numhg(hfx_sysdata%nuotot)
    maxndg = listdptrg(hfx_sysdata%nuotot) + numdg(hfx_sysdata%nuotot)
    call re_alloc(listhg, 1, maxnhg, name='listhg', routine='setup_hfx')
    call re_alloc(listdg, 1, maxndg, name='listdg', routine='setup_hfx')

    do io = 1, hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        do jo = 1, numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
            hfx_sysdata%listh(hfx_sysdata%listhptr(iio)+1: &
              hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
          listdg(listdptrg(io)+1:listdptrg(io)+numdg(io)) = hfx_sysdata%listh( &
            hfx_sysdata%listhptr(iio)+1:hfx_sysdata%listhptr(iio)+hfx_sysdata%numh(iio))
        enddo
      endif

      call MPI_Bcast(listhg(listhptrg(io)+1), numhg(io), MPI_integer,&
        BNode, MPI_Comm_World, MPIerror)
      call MPI_Bcast(listdg(listdptrg(io)+1), numdg(io), MPI_integer,&
        BNode, MPI_Comm_World, MPIerror)
    enddo

!   Globalise Dscf
    call re_alloc(Dscfg, 1, maxndg, 1, hfx_sysdata%nspin, name='Dscfg', &
      routine='setup_hfx')
    Dscfg(:,:) = 0.0_dp

    do io = 1, hfx_sysdata%nuotot
      call WhichNodeOrb(io, Nodes, BNode)
      if ( Node == BNode ) then
        call GlobalToLocalOrb(io, Node, Nodes, iio)
        do ispin = 1, hfx_sysdata%nspin
          do jo = 1, hfx_sysdata%numh(iio)
            Dscfg(listdptrg(io)+jo,ispin) = &
              Dscf(hfx_sysdata%listhptr(iio)+jo,ispin)
          enddo
        enddo
      endif
      do ispin = 1, hfx_sysdata%nspin
        call MPI_Bcast(Dscfg(listdptrg(io)+1,ispin), numdg(io), &
          MPI_double_precision, BNode, MPI_Comm_World, MPIerror)
      enddo
    enddo
 
#endif

!!   jjunquer: For debugging
!    write(6,'(a,2i8)')' setup_hfx, Node, nuotot = ', Node, hfx_sysdata%nuotot
!    write(6,'(a,2i8)')' setup_hfx, Node, norb   = ', Node, hfx_sysdata%norb
!    write(6,'(a,2i8)')' setup_hfx, Node, nspin  = ', Node, hfx_sysdata%nspin
!    write(6,'(a,2i8)')' setup_hfx, Node, maxndg = ', Node, maxndg
!!   end jjunquer: End debugging

!   Allocate the matrix to store the exact exchange Hamiltonian matrix elements
!   This corresponds with Eq. (9) of Ref. \cite Shang:2011
!   First  dimension: number of orbitals in the unit cell
!   Second dimension: number of orbitals in the auxiliary supercell
!   Third  dimension: number of spin components
!   The dimensions are the same for all the Nodes when this is run in parallel
!    call re_alloc(H_EXX, 1, hfx_sysdata%nuotot, 1, hfx_sysdata%norb, 1, &
!      hfx_sysdata%nspin, name='H_EXX', routine='setup_hfx')
!    H_EXX(:,:,:) = 0.0_dp

!   Allocate the matrix to store the maximum value of the density matrix 
!   between angular momentum shells
!   First  dimension: number of orbitals in the auxiliary supercell
!   Second dimension: number of orbitals in the auxiliary supercell
!   The dimensions are the same for all the Nodes when this is run in parallel
!    call re_alloc(P_max, 1, hfx_sysdata%norb, 1, hfx_sysdata%norb, &
!      name='P_max', routine='setup_hfx')
!    P_max(:,:) = 0.0_dp

!    call cpu_time(time_start)
!    call cpu_time(time_dm_start)
!   Allocate the index matrix to map dense(io,jo) to sparse(ind) for reducing the memory.
 !    xmqin 2021
!    call re_alloc(D2Sindx, 1, hfx_sysdata%norb, 1, hfx_sysdata%norb, &
!     &            name='D2Sindx', routine='setup_hfx')
!    D2Sindx(:,:) = 0

#ifdef MPI
!   Transform sparse DM to full matrix
!   Here we apply the condition right below Eq. (45) of the SIESTA paper
!   J. M. Soler et al., J. Phys.: Condens. Matter 14, 2745 (2002).
!   \rho_{\mu \nu} \equiv \rho_{\mu^\prime \nu^\prime}
!   if (\nu \mu) are equivalent by translation to 
!   (\nu^\prime \mu^\prime)
!   and the DM_tmp variable is written in dense format
!   (size (number of orbitals in supercell x number of orbitals in supercell))
!    call sparse_dense(hfx_sysdata%nspin, hfx_sysdata%nuotot, &
!      hfx_sysdata%nuotot, hfx_sysdata%norb, maxndg, numdg, &
!      listdptrg, listdg, Dscfg, DM_tmp)
    call re_alloc(HFX, 1, maxnhg, 1, hfx_sysdata%nspin, name='HFX', routine='setup_hfx')
    HFX(:,:) = 0.0_dp

!   Allocate a temporal variable for the density matrix
!   First  dimension: number of orbitals in the auxiliary supercell
!   Second dimension: number of orbitals in the auxiliary supercell
!   Third  dimension: number of spin components
!   This corresponds with Eq. (10) of Ref. \cite Shang:2011

    call re_alloc(Dscf_max, 1, maxndg, name='Dscf_max', routine='setup_hfx')
    Dscf_max(:) = 0.0_dp

!    call build_D2Sindx( hfx_sysdata%nuotot,  hfx_sysdata%norb, maxndg, numdg, &
!     &                  listdptrg, listdg, D2Sindx )

!   Get the maximum value of the density matrix between shells
!    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
!      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
!      hfx_sysdata%na, hfx_sysdata%isa, maxndg, numdg, listdptrg, &
!      listdg, DM_tmp, P_max)

    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, maxndg, numdg, listdptrg, &
      listdg, D2Sindx, Dscfg, Dscf_max)

!>  Calculate ERIs and store them in RAM
!!
!! \bug Was previously done conditionally through build_hfx_potential
    call timer('ERI',1)
!    call cpu_time(time_start)
    call evaluate_eri(libint_data, hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%cell, hfx_sysdata%cell_r, &
      maxnhg, numhg, listhptrg, listhg, hfx_optdata, Dscfg, Dscf_max, HFX)
!    call cpu_time(time_end)

!    if(node.eq.0)
!    write(6,*) "ERI time :: ", time_end-time_start, " s"
    call timer('ERI',2)
    ! part H_EXX(norb,norb,nspin) per node : parallel over list_mn and list_uv
    ! u and m are local in this node, to calc H_EXX(u,m,nspin), otherwise it
    ! will be zero !
!      listdptrg, listdg, Dscfg, DM_tmp)
    call re_alloc(HFXg, 1, maxnhg, 1, hfx_sysdata%nspin, name='HFXg', routine='setup_hfx')
    HFXg(:,:) = 0.0_dp

    ! Get all H_EXX
    call MPI_AllReduce(HFX(1,1), HFXg(1,1), &
      maxnhg*hfx_sysdata%nspin, &
      MPI_double_precision, MPI_Sum, MPI_Comm_World, MPIerror)
  
    HFX(:,:) = HFXg(:,:)
    call de_alloc(HFXg, name='HFXg', routine='setup_hfx')
#else
!    call sparse_dense(hfx_sysdata%nspin, hfx_sysdata%nuotot, &
!      Dscf, DM_tmp)
!    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
!      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
!      hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%maxnh, &
!      hfx_sysdata%numh, hfx_sysdata%listhptr, hfx_sysdata%listh, &
!      DM_tmp, P_max)

    call re_alloc(HFX, 1, maxnh, 1, hfx_sysdata%nspin, name='HFX', routine='setup_hfx')
    HFX(:,:) = 0.0_dp

!   Allocate a temporal variable for the density matrix
!   First  dimension: number of orbitals in the auxiliary supercell
!   Second dimension: number of orbitals in the auxiliary supercell
!   Third  dimension: number of spin components
!   This corresponds with Eq. (10) of Ref. \cite Shang:2011
    call re_alloc(Dscf_max, 1, maxnd, name='Dscf_max', routine='setup_hfx')
    Dscf_max(:) = 0.0_dp

    call build_D2Sindx( hfx_sysdata%nuotot,  hfx_sysdata%norb, hfx_sysdata%maxnh, hfx_sysdata%numh, &
     &                  hfx_sysdata%listhptr, hfx_sysdata%listh, D2Sindx )

    call get_pmax_shell(hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%maxnh, hfx_sysdata%numh, &
      hfx_sysdata%listhptr, hfx_sysdata%listh, D2Sindx, Dscf, Dscf_max)


!> Calculate electron repulsion integrals (ERIs) 
!!
!! \bug Was previously done conditionally through build_hfx_potential
    call timer('ERI',1)
    call evaluate_eri(libint_data, hfx_sysdata%nspin, hfx_sysdata%norb, &
      hfx_sysdata%iaorb, hfx_sysdata%iphorb, hfx_sysdata%nuotot, &
      hfx_sysdata%na, hfx_sysdata%isa, hfx_sysdata%cell, hfx_sysdata%cell_r, &
      hfx_sysdata%maxnh, hfx_sysdata%numh, hfx_sysdata%listhptr, hfx_sysdata%listh,  &
      hfx_optdata, D2Sindx, Dscf_max, HFX)
    call timer('ERI',2)

#endif


!    call cpu_time(time_dm_end)

!    write(6,'(a, f12.6, a, I6)') "ERIs time = ", time_dm_end-time_dm_start, &
!                                 "(Sec), on Node =", Node

    if ( hfx_sysdata%nspin == 1 ) then
      spin_factor = 1.0_dp
    else
      spin_factor = 2.0_dp
    end if

#ifdef MPI
!   Here H_EXX and DM is global
!   Calculate EXX from dense HFX and DM.
    E_HFX = 0.0_dp
    do ispin=1,hfx_sysdata%nspin
      do io = 1, hfx_sysdata%nuotot
        do j = 1, numhg(io)
           indg = listhptrg(io) + j
           E_HFX = E_HFX - 0.25_dp*spin_factor*HFX(indg,ispin) &
      &            * Dscfg(indg,ispin)*0.25_dp
        enddo
      enddo
    enddo

#else

    E_HFX = 0.0_dp
    do ispin=1,hfx_sysdata%nspin
      do io = 1, hfx_sysdata%nuotot
        do j = 1, hfx_sysdata%numh(io)
           ind = hfx_sysdata%listhptr(io) + j
           E_HFX = E_HFX - 0.25_dp*spin_factor*HFX(ind,ispin) &
       &           * Dscf(ind,ispin)*0.25_dp
        enddo
      enddo
#endif

    Exc = Exc + E_HFX

    if ( Node .eq. 0 ) then
      write(6,'(a23,f18.6,1x,a3)') 'setup_hfx: HFX energy:',E_HFX/eV, " eV"
    endif
!    write(6,'(a,i2,a,f12.6,a)') 'setup_hfx: HFX energy (node=', Node, ") ", &
!      E_HFX/eV, " eV"


#ifdef MPI
!   Global H_EXX to local Hmat, 
!   io: global orbital, 
!   iio: local in this node 
    do ispin = 1, hfx_sysdata%nspin
      do iio = 1, hfx_sysdata%nuo
        call LocalToGlobalOrb(iio, Node, Nodes, io)
        do j = 1, hfx_sysdata%numh(iio)
           ind = hfx_sysdata%listhptr(iio) + j
           indg = listhptrg(io) + j
           Hmat(ind,ispin) = Hmat(ind,ispin) &
                           - HFX(indg,ispin)*spin_factor*0.5_dp*0.25_dp
        enddo
      enddo
    enddo

#else

    do ispin=1,hfx_sysdata%nspin
      do io=1,hfx_sysdata%nuo
        do j = 1,hfx_sysdata%numh(io)
          ind = hfx_sysdata%listhptr(io) + j
          !you need to a combine code  for ns=1, and ns=2
          Hmat(ind,ispin) = Hmat(ind,ispin) &
                          - HFX(ind,ispin)*spin_factor*0.5_dp*0.25_dp
        enddo
      enddo
    enddo
#endif

!    call cpu_time(time_end)

!    if ( Node .eq. 0 ) then
!      write(6,'(a, F12.6, a, I6)') "setup_hfx: Build HFX time =",  &
!        time_end-time_start, " (secs) on Node =", Node
!    endif

    call de_alloc(HFX, name='HFX', routine='setup_hfx')
!    call de_alloc(Dscf, name='Dscf', routine='setup_hfx')
    call de_alloc(Dscf_max, name='Dscf_max', routine='setup_hfx')
!    call de_alloc(D2Sindx, name='D2Sindx', routine='setup_hfx')
#ifdef MPI
    call de_alloc(Dscfg, name='Dscfg', routine='setup_hfx')
    call de_alloc(listdg, name='listdg', routine='setup_hfx')
    call de_alloc(listdptrg, name='listdptrg', routine='setup_hfx')
    call de_alloc(listhg, name='listhg', routine='setup_hfx')
    call de_alloc(listhptrg, name='listhptrg', routine='setup_hfx')
    call de_alloc(numdg, name='numdg', routine='setup_hfx')
    call de_alloc(numhg, name='numhg', routine='setup_hfx')
#endif

    call timer('HFX', 2)

  end subroutine setup_hfx

! *****************************************************************************
!> \brief Computes the Hartree-Fock potential corresponding to a given
!!        Hamiltonian
!!
!! We impose the same sparsity for the exact exchange matrix 
!! elements as in standard SIESTA
!!
!! The Schwarz screening is done in order to reduce the number of
!! ERIs to be computed
!! \f{eqnarray*}{
!!   \vert \left(\mu\nu^{\prime} \vert \eta^{\prime}\lambda^{\prime} 
!!   \right)_{\rm SR} \vert
!!   \le \sqrt{
!!   \left(\mu\nu^{\prime} \vert \mu \nu^{\prime} \right)_{\rm SR} 
!!  \left(\eta^{\prime}\lambda^{\prime}\vert\eta^{\prime}\lambda^{\prime}\right)
!!   _{\rm SR}
!!   }
!! \f}
!! Only when this upper bound of the ERI is smaller than 
!! \f$ \epsilon_{\rm Schwarz} \f$, the four center integral is computed.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] nspin:       Number of spin components
!! \param[in] norb:        Number of atomic orbitals in the supercell
!! \param[in] iaorb:       Atomic index of each orbital
!! \param[in] iphorb:      Orbital index of each orbital in its atom
!! \param[in] nuotot:      Number of orbitals in the unit cell
!! \param[in] na:          Number of atoms in the supercell
!! \param[in] isa:         Species index of each atom
!! \param[in] cell:        Lattice vectors of the supercell in real space
!! \param[in] rcell:       Reciprocal lattice vectors of the supercell
!! \param[in] hfx_optdata: Initialized Hartree-Fock options data structure
!! \param[in] DM_tmp:      Density matrix in dense format
!! \param[in] P_max:       Largest value of the density matrix between two
!!                         atomic orbitals of given angular momentums
!!                         in the supercell
!! \param[in,out] H_EXX:   Hartree-Fock Hamiltonian matrix elements
! *****************************************************************************
  subroutine evaluate_eri( libint_data, nspin, norb, iaorb, iphorb, nuotot, na,&
 &                         isa, cell, rcell, maxnh, numh, listhptr, listh, &
                           hfx_optdata, Dscf, Dscf_max, HFX )

    use precision,      only: i8b               ! Increased precision for
                                                !    integer kinds
    use parallel,       only: Node              ! Number of the local node
    use parallel,       only: Nodes             ! Total number of nodes
    use parallelsubs,   only: GetNodeOrbs       ! Calculates the number of 
                                                !    orbitals stored on the 
                                                !    local node.
    use parallelsubs,   only: GlobalToLocalOrb  ! Converts an orbital index 
                                                !    in the global frame to the
                                                !    local frame if the orbital
                                                !    is local to this node. 
                                                !    Otherwise the pointer is
                                                !    return as zero and can 
                                                !    therefore be used to test 
                                                !    whether the orbital is 
                                                !    local or not.
    use parallelsubs,   only: LocalToGlobalOrb  ! Converts an orbital index in 
                                                !    the local frame to the 
                                                !    global frame
    use parallelsubs,   only: WhichNodeOrb      ! Given the global orbital 
                                                !    pointer, this routine
                                                !    returns the node number
                                                !    where this is stored.
    use atm_types,      only: species           ! All the information related
                                                !    with a given chemical 
                                                !    species
    use atm_types,      only: species_info      ! Derived type for all the 
                                                !    information related with
                                                !    a given chemical species
    use atomlist,       only: indxuo            ! Index of equivalent orbital 
                                                !    in "u" cell
    use atmfuncs,       only: lofio             ! Returns total angular momentum
                                                !    quantum number of a given 
                                                !    atomic basis
    use atmfuncs,       only: mofio             ! Returns m quantum number of a
                                                !    given atomic basis
    use alloc,          only: re_alloc          ! Reallocation routine
    use alloc,          only: de_alloc          ! Deallocation routine
    use nao2gto_types
    use nao2gto_data,   only: list_ij           ! List of pair of shells that
                                                !    are preselected according 
                                                !    to the Schwarz screening
                                                !    Orbital i is within the
                                                !    home unit cell,
                                                !    Orbital j might be anywhere
                                                !    in the auxiliary supercell
    use nao2gto_data,   only: list_kl           ! List of pair of shells that
                                                !    are preselected according 
                                                !    to the Schwarz screening
                                                !    Orbitals k and l might
                                                !    be anywhere in the 
                                                !    auxiliary supercell 
    use nao2gto_data,   only: nso               ! Number of Spherical Harmonic 
                                                !    Gaussians for a given 
                                                !    angular momentum
    use nao2gto_data,   only: nco               ! Number of Cartesian Gaussians
                                                !    for a given 
                                                !    angular momentum
    use nao2gto_data,   only: hfx_gtos          ! Information about the 
                                                !    Gaussian type orbitals of
                                                !    a given species
    use nao2gto_data,   only: D2Sindx
    use nao2gto_data,   only: eri_prescreen     ! (mu nu | mu nu), following
                                                !    Sec. III E 1 of 
                                                !    Ref. \cite Shang:2011
    use nao2gto_data,   only: sfc_shell         ! 
    use nao2gto_data,   only: sfc_pgf
    use nao2gto_data,   only: pair_dist_radii_pgf
    use nao2gto_contract,  only : calc_contract_eri
    use nao2gto_libint, only: Libint_t
#ifdef MPI
    use mpi_siesta
#endif

    implicit none

!   Arguments
    type(Libint_t), intent(inout) :: libint_data
    integer,  intent(in) :: nspin         ! Number of spin components
    integer,  intent(in) :: na            ! Number of atoms in the supercell
    integer,  intent(in) :: norb          ! Number of orbitals in the supercell
                                          !   All the nodes see the same number
    integer,  intent(in) :: nuotot        ! Number of orbitals in the unit cell
                                          !   All the nodes see the same number
    integer,  intent(in) :: iaorb(norb)   ! Atomic index of each orbital  
    integer,  intent(in) :: iphorb(norb)  ! Orbital index of each orbital
                                          !   in its atom
    integer,  intent(in) :: isa(na)       ! Species index of each atom

    integer,  intent(in) :: maxnh
    integer,  intent(in) :: numh(nuotot)
    integer,  intent(in) :: listh(maxnh)
    integer,  intent(in) :: listhptr(nuotot)

      
    type(hfx_options_type), intent(in) :: hfx_optdata
                                          ! Initialized Hartree-Fock 
                                          !   options data structure
    real(dp), intent(in) :: cell(3,3)     ! Supercell lattice vectors in 
                                          !   real space (units in Bohrs)
                                          !   First  index: component
                                          !   Second index: vector index
    real(dp), intent(in) :: rcell(3,3)    ! Reciprocal lattice vectors of the
                                          !   supercell in real space 
                                          !   (units in Bohr^-1)
                                          !   First  index: component
                                          !   Second index: vector index
!    integer, intent(in)  :: D2Sindx(norb, norb)
    real(dp), intent(in) :: Dscf(maxnh,nspin)
                                          ! Density matrix in dense form
                                          ! NOTE: When running in parallel,
                                          !   this is core independent
    real(dp), intent(in) :: Dscf_max(maxnh)
                                          ! Largest value of the density matrix
                                          !   between two atomic orbitals in the
    real(dp), intent(inout) :: HFX(maxnh, nspin)
                                          ! Hartree-Fock Hamiltonian 
                                          !   matrix elements

!   Local variables
    integer  :: nsoi, nsoj, nsok, nsol, num_a, num_b, num_c, num_d,   &
                index_ij, index_kl

    integer  :: io                 ! First  orbital in the pair ij
    integer  :: jo                 ! Second orbital in the pair ij
    integer  :: ko                 ! First  orbital in the pair kl
    integer  :: lo                 ! Second orbital in the pair kl
    integer  :: is                 ! Atomic species of orbital io
    integer  :: js                 ! Atomic species of orbital jo
    integer  :: ks                 ! Atomic species of orbital ko
    integer  :: ls                 ! Atomic species of orbital lo
    integer  :: ioa                ! Orbital index of orbital io in its atom
    integer  :: joa                ! Orbital index of orbital jo in its atom
    integer  :: koa                ! Orbital index of orbital ko in its atom
    integer  :: loa                ! Orbital index of orbital lo in its atom
    integer  :: l_i                ! Orbital angular momentum of orbital io
    integer  :: l_j                ! Orbital angular momentum of orbital jo
    integer  :: l_k                ! Orbital angular momentum of orbital ko
    integer  :: l_l                ! Orbital angular momentum of orbital lo
    integer  :: npgfi              ! Number of gaussians required to expand
                                   !   the radial part of orbital io
    integer  :: npgfj              ! Number of gaussians required to expand
                                   !   the radial part of orbital jo
    integer  :: npgfk              ! Number of gaussians required to expand
                                   !   the radial part of orbital ko
    integer  :: npgfl              ! Number of gaussians required to expand
                                   !   the radial part of orbital lo
    integer  :: ncoi               ! Total number of primitive Cartesian 
                                   !   Gaussian functions required to expand
                                   !   the shell l_i or orbital io
    integer  :: ncoj               ! Total number of primitive Cartesian 
                                   !   Gaussian functions required to expand
                                   !   the shell l_j or orbital jo
    integer  :: ncok               ! Total number of primitive Cartesian 
                                   !   Gaussian functions required to expand
                                   !   the shell l_k or orbital ko
    integer  :: ncol               ! Total number of primitive Cartesian 
                                   !   Gaussian functions required to expand
                                   !   the shell l_l or orbital lo

    integer  :: ncells             ! Number of unit cells in the supercell
    integer  :: ishell             ! Index of the shell to which orbital io 
                                   !   belongs
                                   !   For instance, in Si with a DZP basis,
                                   !   there will be five shells
                                   !   shell # 1 = 3s first -zeta
                                   !   shell # 2 = 3s second-zeta
                                   !   shell # 3 = 3p first -zeta
                                   !   shell # 4 = 3p second-zeta
                                   !   shell # 5 = 3d first -zeta
                                   !   Then, given the index of a NAO,
                                   !   ishell will determine the angular 
                                   !   momentum shell to which this orbital 
                                   !   belongs
    integer  :: jshell             ! Same as ishell, but for the orbital jo
    integer  :: kshell             ! Same as ishell, but for the orbital ko
    integer  :: lshell             ! Same as ishell, but for the orbital lo

    integer(i8b) :: shell_eri_calc
    integer(i8b) :: spher_eri_calc
    integer(i8b) :: spher_eri_store

    integer(i8b) :: neris_tmp      ! Total number of the primitive eris computed
                                   !   per shell
    integer(i8b) :: tot_pgto       ! Total number of the primitive eris computed
    integer      :: i_list_ij      ! Counter for the loop on pairs of shells
                                   !   ij
    integer      :: i_list_kl      ! Counter for the loop on pairs of shells
                                   !   kl
    integer      :: list_kl_local  ! Number of pairs of orbitals kl that will 
                                   !   be handled by the local node
    integer      :: i_list_kl_local! Counter for the loop on the pair
                                   !   ok kl orbitals handled by the local node

#ifdef MPI
    integer  :: MPIerror, Request, num_loc, Status(MPI_Status_Size)
#endif
    real(dp) :: DM_max             ! P_{\rm screening} in Eq. (17) of 
                                   !   Ref. \cite Shang:2020
    real(dp) :: log10_pmax         ! Logarithm in base 10 of P_{\rm screening} 
    real(dp) :: eps_temp           ! \sqrt{(ij \vert ij) (kl \times kl) }  
                                   !   or this value multiplied by 
                                   !   P_{\rm screening}, as in Eq. (17) of
                                   !   Ref. \cite Shang:2020
    real(dp) :: nao_eri
    real(dp) :: max_val2_set       ! Calculate the upper limit of the 
                                   !   four center integral involving 
                                   !   only two orbitals of two shells
                                   !   according to the recipe given by 
                                   !   Guidon et al. in 
                                   !   J. Chem. Theory Comput. 5,3010–3021(2009)
                                   !   Eq. (49) 
                                   !   log( (mu,nu|mu,nu)(R_{mu,nu} ) 
                                   !   \approx a_2*R_{mu,nu}^2 + a_0
    real(dp) :: max_contraction_val

    real(dp) :: ri(3)              ! Position where the orbital io is centered
    real(dp) :: rj(3)              ! Position where the orbital jo is centered
    real(dp) :: rk(3)              ! Position where the orbital ko is centered
    real(dp) :: rl(3)              ! Position where the orbital lo is centered
    real(dp) :: rij2               ! Square of the distance between orbital 
                                   !    io and jo
    real(dp) :: rkl2               ! Square of the distance between orbital 
                                   !    ko and lo

!   Electron repulsion integrals between the real spherical harmonics
    real(dp), dimension(:,:,:,:), pointer :: eri => null()

    type(species_info), pointer :: ispp => null()
                                   ! Pointer to the variable with the 
                                   !    species info or orbital io
    type(species_info), pointer :: jspp => null()
                                   ! Pointer to the variable with the 
                                   !    species info or orbital jo
    type(species_info), pointer :: kspp => null()
                                   ! Pointer to the variable with the 
                                   !    species info or orbital ko
    type(species_info), pointer :: lspp => null()
                                   ! Pointer to the variable with the 
                                   !    species info or orbital lo
    type(gto_info_type), pointer :: igto => null()
                                   ! Pointer to the variable with the 
                                   !    information of the Gaussian orbitals 
                                   !    required to expand orbital io
    type(gto_info_type), pointer :: jgto => null()
                                   ! Pointer to the variable with the 
                                   !    information of the Gaussian orbitals 
                                   !    required to expand orbital jo
    type(gto_info_type), pointer :: kgto => null()
                                   ! Pointer to the variable with the 
                                   !    information of the Gaussian orbitals 
                                   !    required to expand orbital ko
    type(gto_info_type), pointer :: lgto => null()
                                   ! Pointer to the variable with the 
                                   !    information of the Gaussian orbitals 
                                   !    required to expand orbital lo

    type(hfx_screen_coeff_type), dimension(:,:), pointer :: &
      tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2

    real(dp) time1, time2, time_tot

! ------------------------------------------------------------------------------

!   Compute the number of unit cells included in the supercell
    ncells = norb / nuotot

!!   For debugging
!    write(6,'(a,3i5)')                         &
! &    'evaluate_eri: Node, Nodes, ncells   = ', Node, Nodes, ncells
!    write(6,'(a,3i5)')                         &
! &    'evaluate_eri: Node, Nodes, na       = ', Node, Nodes, na
!    write(6,'(a,3i5)')                         &
! &    'evaluate_eri: Node, Nodes, norb     = ', Node, Nodes, norb
!    write(6,'(a,3i5)')                         &
! &    'evaluate_eri: Node, Nodes, nuotot   = ', Node, Nodes, nuotot
!    write(6,*)                         &
! &    'evaluate_eri: Node, Nodes, cell     = ', Node, Nodes, cell
!    write(6,*)                         &
! &    'evaluate_eri: Node, Nodes, rcell    = ', Node, Nodes, rcell
!!   End debugging

#ifdef MPI
!   Compute the number of pairs kl that will be handled by the local node
    call GetNodeOrbs(list_kl%nelement, Node, Nodes, list_kl_local)

!!   For debugging
!    write(6,'(a,5i7)')                                                  &
! &    'evaluate_eri: Node, Nodes, list_ij%nelement, list_kl%nelement, list_kl_local = ',  &
! &                   Node, Nodes, list_ij%nelement, list_kl%nelement, list_kl_local 
!!   End debugging

#endif

!   Initialize the number of eris to be computed
    neris_tmp       = 0_i8b
    shell_eri_calc  = 0_i8b
    spher_eri_calc  = 0_i8b
    spher_eri_store = 0_i8b
    tot_pgto        = 0_i8b

    time_tot= 0.0_dp
    time1 =0.0_dp
    time2 =0.0_dp

!   Loop on all the pairs of shells ij 
!   All nodes see the same number of pairs of shells list_ij%nelement
    do i_list_ij = 1, list_ij%nelement

!     Identify the atomic orbitals in the pair ij, their centers,
!     and the square of the relative distance
      io   = list_ij%element(i_list_ij)%pair(1)
      jo   = list_ij%element(i_list_ij)%pair(2)
      ri   = list_ij%element(i_list_ij)%r1
      rj   = list_ij%element(i_list_ij)%r2
      rij2 = list_ij%element(i_list_ij)%dist2

!!     For debugging
!      write(6,'(a,3i5)')                             &
! &      'evaluate_eri: Node, Nodes, i_list_ij = ',   &
! &       Node, Nodes, i_list_ij 
!      write(6,'(a,3i5)')                             &
! &      '   evaluate_eri: Node, Nodes, pair1 = ',    &
! &          Node, Nodes, io     
!      write(6,'(a,3i5)')                             &
! &      '   evaluate_eri: Node, Nodes, pair2 = ',    &
! &          Node, Nodes, jo     
!      write(6,'(a,2i5,3f12.5)')                      &
! &      '   evaluate_eri: Node, Nodes, r1    = ',    &
! &          Node, Nodes, ri
!      write(6,'(a,2i5,3f12.5)')                      &
! &      '   evaluate_eri: Node, Nodes, r2    = ',    &
! &          Node, Nodes, rj
!      write(6,'(a,2i5,f12.5)')                       &
! &      '   evaluate_eri: Node, Nodes, dist2 = ',    &
! &          Node, Nodes, rij2
!!     End debugging

!     Extract information for the orbital io
      is = isa(iaorb(io))
      ispp => species(is)
      igto => hfx_gtos(is)
      ioa = iphorb(io)
      l_i = lofio(is, ioa)
      npgfi = igto%orbnl_contract(ispp%orb_index(ioa))
      ncoi = nco(l_i)*npgfi

!     Extract information for the orbital jo
      js = isa(iaorb(jo))
      jspp => species(js)
      jgto => hfx_gtos(js)
      joa = iphorb(jo)
      l_j = lofio(js, joa)
      npgfj = jgto%orbnl_contract(jspp%orb_index(joa))
      ncoj = nco(l_j)*npgfj

!     Determine the shells to which orbitals io and jo belongs to
      ishell = ispp%orb_index(ioa)
      jshell = jspp%orb_index(joa)

      index_ij = list_ij%element(i_list_ij)%nl_index

#ifdef MPI
      do i_list_kl_local = 1, list_kl_local
         call LocalToGlobalOrb(i_list_kl_local, Node, Nodes, i_list_kl)
#else
      do i_list_kl = 1, list_kl%nelement
#endif

!       Identify the atomic orbitals in the pair kl, their centers,
!       and the square of the relative distance
        ko = list_kl%element(i_list_kl)%pair(1)
        lo = list_kl%element(i_list_kl)%pair(2)
        rk = list_kl%element(i_list_kl)%r1
        rl = list_kl%element(i_list_kl)%r2
        rkl2 = list_kl%element(i_list_kl)%dist2

!       Extract information for the orbital ko
        ks = isa(iaorb(ko))
        kspp => species(ks)
        kgto => hfx_gtos(ks)
        koa = iphorb(ko)
        l_k = lofio(ks, koa)
        npgfk = kgto%orbnl_contract(kspp%orb_index(koa))
        ncok =nco(l_k)*npgfk

!       Extract information for the orbital lo
        ls = isa(iaorb(lo))
        lspp => species(ls)
        lgto => hfx_gtos(ls)
        loa = iphorb(lo)
        l_l = lofio(ls, loa)
        npgfl = lgto%orbnl_contract(lspp%orb_index(loa))
        ncol = nco(l_l)*npgfl

!       Determine the shells to which orbitals ko and lo belongs to
        kshell = kspp%orb_index(koa)
        lshell = lspp%orb_index(loa)

        index_kl = list_kl%element(i_list_kl)%nl_index

        if ( index_kl .le. index_ij ) then

!         Here we impose the same sparsity as in SIESTA
!         The four center integrals will be computed if and only if 
!         the four atomic orbitals are neighbors all with all

           if((D2Sindx(io,ko).eq.0).and.(D2Sindx(io,lo).eq.0).and.(D2Sindx(jo,ko).eq.0) &
               .and.(D2Sindx(jo,lo).eq.0) )  cycle

!         Calculate the maximum value of the four center integral 
!         involving only two orbitals of two shells
!         according to the recipe given by Guidon et al.
!         in J. Chem. Theory Comput. 5, 3010–3021 (2009),
!         Eq. (49) 
!         log( (mu,nu|mu,nu)(R_{mu,nu} ) \approx a_2*R_{mu,nu}^2 + a_0
          max_val2_set = sfc_shell(lshell,kshell,ls,ks)%x(1)*rkl2 + &
            sfc_shell(lshell,kshell,ls,ks)%x(2)

!         Compute the upper limit of the ERI, according to the 
!         Schwartz inequality.
!         (ij \vert ij) (kl \times kl) 
          eps_temp =  eri_prescreen(D2Sindx(io,jo))*eri_prescreen(D2Sindx(ko,lo))
!         \sqrt{ (ij \vert ij) (kl \times kl) }
!         This square root is what is used in the Schwarz inequality
          eps_temp = dsqrt(eps_temp)

!         Following Eq. (17) of Ref. \cite Shang:2020, 
!         Finally, we use the density matrix to further reduce the 
!         number of ERIs, that the maximal value of the density 
!         matrix of each shell (P_max) is calculated during every SCF cycle,
!         and then the density matrix screening is,
!         P_{\rm screening} \times 
!         \sqrt{(ij \vert ij) (kl \times kl) } .le. \epsilon_{\rm Schwarz},
!         where 
!         P_{\rm screening} = \max( \vert P^{IK}_{\rm max} \vert, 
!                                   \vert P^{IL}_{\rm max} \vert, 
!                                   \vert P^{JK}_{\rm max} \vert, 
!                                   \vert P^{JL}_{\rm max} \vert ).
!         Here four density matrix elements are needed for the maximal
!         value because of the 8-fold full permutation symmetry of the ERIs
!         is exploited in the implementation.
!         The maximal density matrix elements are chosen from the 
!         density matrix of the previous SCF cycle, 
!         which produce a stable direct SCF cycle 

          if ( hfx_optdata%DM_trunc ) then
!           Compute P_{\rm screening} in the previous equation
               DM_max = max( Dscf_max(D2Sindx(io,ko)), Dscf_max(D2Sindx(io,lo)), &
                        Dscf_max(D2Sindx(jo,ko)), Dscf_max(D2Sindx(jo,lo)) )

            if ( DM_max <= 0.0_dp ) then
              log10_pmax = log_zero
            else
              log10_pmax = log10(DM_max)
            end if

!           Compute the product of 
!           P_{\rm screening} \times \sqrt{(ij \vert ij) (kl \times kl) } 
!           as in Eq. (17) of Ref. \cite Shang:2020
            eps_temp = DM_max * eps_temp

          else
            DM_max = 1.0_dp
          endif

!         The Schwarz inequality to reduce the number of the required
!         ERIs is done here
          if ( eps_temp .gt. hfx_optdata%eps_schwarz ) then

            shell_eri_calc = shell_eri_calc + 1

            tmp_R_1 => pair_dist_radii_pgf(:,:,jshell,ishell,js,is)
            tmp_R_2 => pair_dist_radii_pgf(:,:,lshell,kshell,ls,ks)
            tmp_screen_pgf1 => sfc_pgf(:,:,jshell,ishell,js,is)
            tmp_screen_pgf2 => sfc_pgf(:,:,lshell,kshell,ls,ks)

            max_contraction_val = igto%orbnl_contraction_coeff(ishell) * &
                                  jgto%orbnl_contraction_coeff(jshell) * &
                                  kgto%orbnl_contraction_coeff(kshell) * &
                                  lgto%orbnl_contraction_coeff(lshell) * &
                                  DM_max

            call re_alloc(eri, 1, nso(l_i), 1, nso(l_j), 1, nso(l_k), &
              1, nso(l_l), name='eri', routine='evaluate_eri')
            eri(:,:,:,:) = 0.0_dp

!           The four center integral is computed here
            call calc_contract_eri( libint_data, cell, rcell, ri, rj, rk, rl, &
              npgfi, npgfj, npgfk, npgfl, l_i, l_j, l_k, l_l, &
              ncoi, ncoj, ncok, ncol, &
              igto%orbnl_zeta(1:npgfi,ispp%orb_index(ioa)), &
              jgto%orbnl_zeta(1:npgfj,jspp%orb_index(joa)), &
              kgto%orbnl_zeta(1:npgfk,kspp%orb_index(koa)), &
              lgto%orbnl_zeta(1:npgfl,lspp%orb_index(loa)), &
              igto%sphi(1:ncoi,ioa:ioa+nso(l_i)-1), &
              jgto%sphi(1:ncoj,joa:joa+nso(l_j)-1), &
              kgto%sphi(1:ncok,koa:koa+nso(l_k)-1), &
              lgto%sphi(1:ncol,loa:loa+nso(l_l)-1), &
              neris_tmp, max_contraction_val, max_val2_set, log10_pmax, &
              tmp_R_1, tmp_R_2, tmp_screen_pgf1, tmp_screen_pgf2, &
              hfx_optdata, eri )

            tot_pgto = tot_pgto + neris_tmp

            do nsoi = 1, nso(l_i)
              do nsoj = 1, nso(l_j)
                do nsok = 1, nso(l_k)
                  do nsol = 1, nso(l_l)
                    spher_eri_calc = spher_eri_calc + 1
                    if ( DM_max*dabs(eri(nsoi,nsoj,nsok,nsol)*2.0_dp) .gt. &
                         hfx_optdata%eps_stored ) then
                      spher_eri_store = spher_eri_store + 1
                      num_a = io + nsoi-1
                      num_b = jo + nsoj-1
                      num_c = ko + nsok-1
                      num_d = lo + nsol-1
!                     Here we transform the eris from Hartrees (as LIBINT)
!                     to Rydbergs (internal units in SIESTA)
                      nao_eri = eri(nsoi,nsoj,nsok,nsol) * 2.0_dp
                      
!                      write(node+100,*) num_a,num_b,num_c,num_d,nao_eri                      
!                     call cpu_time(time1)
                      call hfx_matrix(nspin, ncells, nuotot, norb, &
                              num_a, num_b, num_c, num_d, &
                              maxnh, numh, listhptr, listh,  &
                              nao_eri, Dscf, HFX )
!                     call cpu_time(time2)
!                     time_tot =time_tot + time2 -time1
                    endif
                  enddo
                enddo
              enddo
           enddo

           call de_alloc(eri, name='eri', routine='evaluate_eri')

          endif ! Schwarz inequality
        endif
      enddo !mn
    enddo   !uv

!    if(node.eq.0)
!     write(6,*) "HFX time : ", time_tot, " s"
!!   For debugging
!    write(6,'(a,i20)') &                             
! &    'evaluate_eri: neris_tmp       = ', neris_tmp
!    write(6,'(a,i20)') &
! &    'evaluate_eri: shell_eri_calc  = ', shell_eri_calc
!    write(6,'(a,i20)') &
! &    'evaluate_eri: spher_eri_calc  = ', spher_eri_calc
!    write(6,'(a,i20)') &
! &    'evaluate_eri: spher_eri_store = ', spher_eri_store
!    write(6,'(a,i20)') &
! &    'evaluate_eri: tot_pgto        = ', tot_pgto
!!   End debugging


  end subroutine evaluate_eri

! *****************************************************************************
!> \brief Computes the Hartree-Fock potential corresponding to a given
!!        Hamiltonian
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 10.2013 Created [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] nspin: ...
!! \param[in] ncells: ...
!! \param[in] nuotot: ...
!! \param[in] norb: ...
!! \param[in] io: ...
!! \param[in] jo: ...
!! \param[in] ko: ...
!! \param[in] lo: ...
!! \param[in] nao_eri: ...
!! \param[in] DM_tmp: ...
!! \param[in,out] H_EXX: ...
! *****************************************************************************
  subroutine hfx_matrix(nspin, ncells, nuotot, norb, io, jo, ko, lo, &
                  maxnh, numh, listhptr, listh, nao_eri, Dscf, HFX)

    use atomlist,      only: indxuo
    use nao2gto_index, only: indexsc
    use nao2gto_data,  only: subshell
    use nao2gto_data,   only: D2Sindx    
    implicit none

    ! Arguments
    integer,  intent(in)    :: nspin, ncells, nuotot, norb, io, jo, ko, lo
    integer,  intent(in)    :: maxnh, numh(nuotot), listh(maxnh), listhptr(nuotot)
    real(dp), intent(in)    :: nao_eri
    real(dp), intent(in)    :: Dscf(maxnh,nspin)
    real(dp), intent(inout) :: HFX(maxnh,nspin)

    ! Local variables
    integer  :: ispin, iuo, juo, kuo, luo, llo, jshell, lshell, iushell, &
                jushell, kushell, lushell, index_ij, index_kl, &
                io_trans, jo_trans, ko_trans, lo_trans
    real(dp) :: gint

    ! -------------------------------------------------------------------------

    gint = nao_eri
    iuo = io ! u is always u0
    juo = indxuo(jo)
    kuo = indxuo(ko)
    luo = indxuo(lo)
    llo = indexsc(ko, kuo, lo)

    ! num_n have to trans to play with m0 to get num_mn and
    ! campared to num_uv, so there is num_n_1
    iushell = subshell(iuo)
    jushell = subshell(juo)
    kushell = subshell(kuo)
    lushell = subshell(luo)
    jshell  = subshell(jo)
    lshell  = subshell(llo)

    index_ij = ncells*iushell*(iushell-1)/2 + &
      ((jshell-1)/subshell(nuotot))*iushell + jushell

    index_kl = ncells*kushell*(kushell-1)/2 + &
      ((lshell-1)/subshell(nuotot))*kushell + lushell

    if ( iushell  .eq. jushell  ) gint = gint * 0.5_dp
    if ( kushell  .eq. lushell  ) gint = gint * 0.5_dp
    if ( index_ij .eq. index_kl ) gint = gint * 0.5_dp

    !! HFX
    do ispin=1,nspin

      !         (u0vR|mR'nR")        =       (u0vR|nR"mR')
      ! = (v0u[-R]|m[R'-R]n[R"-R])   = (v0u[-R]|n[R"-R]m[R'-R])
      ! = (m0n[R"-R']|u[-R']v[R-R']) = (m0n[R"-R']|v[R-R']u[-R'])
      ! = (n0m[R'-R"]|u[-R"]v[R-R"]) = (n0m[R'-R"]|v[R-R"]u[-R"])

       ! 1.VEE(1[0]  2[H] | 3[G]  4[N])  (u0v[R]|m[R']n[R"])
       !   VEE(1[0]  2[H] | 4[N]  3[G])  (u0v[R])|n[R"]m[R'])

       if((D2Sindx(io, ko).ne.0) .and. (D2Sindx(jo, lo).ne.0)) then
         HFX(D2Sindx(io, ko), ispin) = HFX(D2Sindx(io, ko), ispin) &
                                +gint*Dscf(D2Sindx(jo, lo), ispin)
       endif

       if((D2Sindx(io, lo).ne.0) .and. (D2Sindx(jo, ko).ne.0)) then
         HFX(D2Sindx(io, lo), ispin) = HFX(D2Sindx(io, lo), ispin) &
                                +gint*Dscf(D2Sindx(jo, ko), ispin)
       endif

      ! 2.VEE(2[0] 1[-H]| 3[G-H] 4[N-H])  (v0u[-R]|m[R'-R]n[R"-R])
      !   VEE(2[0] 1[-H]| 4[N-H] 3[G-H])  (v0u[-R]|n[R"-R]m[R'-R])
         io_trans = indexsc( jo, juo, io )
         ko_trans = indexsc( jo, juo, ko )
         lo_trans = indexsc( jo, juo, lo )

       if((D2Sindx(juo, ko_trans).ne.0) .and. (D2Sindx(io_trans, lo_trans).ne.0)) then
         HFX(D2Sindx(juo, ko_trans), ispin) = HFX(D2Sindx(juo, ko_trans), ispin) &
                             +gint*Dscf(D2Sindx(io_trans, lo_trans), ispin)
       endif

       if((D2Sindx(juo, lo_trans).ne.0) .and. (D2Sindx(io_trans, ko_trans).ne.0)) then
         HFX(D2Sindx(juo, lo_trans), ispin) = HFX(D2Sindx(juo, lo_trans), ispin) &
                             +gint*Dscf(D2Sindx(io_trans, ko_trans), ispin)
       endif

      ! 3.VEE(3[0]  4[N-G] | 1[-G] 2[H-G]) (m0n[R"-R']|u[-R']v[R-R'])
      !   VEE(3[0]  4[N-G] |2[H-G] 1[-G] ) (m0n[R"-R']|v[R-R']u[-R'])

         io_trans = indexsc( ko, kuo, io )
         jo_trans = indexsc( ko, kuo, jo )
         lo_trans = indexsc( ko, kuo, lo )

       if((D2Sindx(kuo, io_trans).ne.0) .and. (D2Sindx(lo_trans, jo_trans).ne.0)) then
         HFX(D2Sindx(kuo, io_trans), ispin) = HFX(D2Sindx(kuo, io_trans), ispin) &
                             +gint*Dscf(D2Sindx(lo_trans, jo_trans),ispin)
       endif

       if((D2Sindx(kuo, jo_trans).ne.0) .and. (D2Sindx(lo_trans, io_trans).ne.0)) then
         HFX(D2Sindx(kuo, jo_trans), ispin) = HFX(D2Sindx(kuo, jo_trans), ispin) &
                             +gint*Dscf(D2Sindx(lo_trans, io_trans),ispin)
       endif

! 4.VEE(4[0]  3[G-N] | 1[-N] 2[H-N])  (n0m[R'-R"]|u[-R"]v[R-R"])
!   VEE(4[0]  3[G-N] | 2[H-N] 1[-N])  (n0m[R'-R"]|v[R-R"]u[-R"])

         io_trans = indexsc( lo, luo, io )
         jo_trans = indexsc( lo, luo, jo )
         ko_trans = indexsc( lo, luo, ko )

       if((D2Sindx(luo, io_trans).ne.0) .and. (D2Sindx(ko_trans, jo_trans).ne.0)) then
         HFX(D2Sindx(luo, io_trans), ispin) = HFX(D2Sindx(luo, io_trans), ispin) &
                             +gint*Dscf(D2Sindx(ko_trans, jo_trans), ispin)
       endif

       if((D2Sindx(luo, jo_trans).ne.0) .and. (D2Sindx(ko_trans, io_trans).ne.0)) then
         HFX(D2Sindx(luo, jo_trans), ispin) = HFX(D2Sindx(luo, jo_trans), ispin) &
                             +gint*Dscf(D2Sindx(ko_trans, io_trans), ispin)
       endif

    enddo

  end subroutine hfx_matrix

end module nao2gto_hfx
