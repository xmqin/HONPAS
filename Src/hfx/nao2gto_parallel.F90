module nao2gto_parallel

  implicit none

  private

  integer, public :: nao2gto_bsize

  integer, save :: NOrbLast = 0
  integer, save :: NOrbNodeLast = 0

  public :: &
&   Broadcast_NAO2GTO, &
&   GetNodeOrbs_NAO2GTO, &
&   LocalToGlobalOrb_NAO2GTO, &
&   set_bsize_NAO2GTO

contains

  ! ***************************************************************************
  !> \brief Broadcasts the relevant NAO2GTO variables to all MPI processors
  !!
  !!
  !! \par History
  !!      10.2012 Created [Xinming Qin]
  !!      12.2018 Extracted from broadcast_basis [Yann Pouillon]
  !!
  !! \note The index of cartesian GTOs, aka orbnl_index_cphi, should not
  !!       be used to calculate ERIs, thus not broadcast.
  !!
  !! \note Errors in calling DGEMM may come from wrong dimensions or missing
  !!       initializations.
  ! ***************************************************************************
  subroutine Broadcast_NAO2GTO()

#ifdef MPI
    use mpi_siesta
    use parallel, only: Node, Nodes
    use atm_types, only: maxn_orbnl, maxnorbs, nspecies
    use sparse_matrices
!    use nao2gto_common
    use nao2gto_data
    use nao2gto_types, only: gto_info_type
    use nao2gto_types, only: maxn_contract, maxnorbs_cphi, l_max, ncon_max
#endif

    implicit none

#ifdef MPI

    ! Local variables
    integer :: isp
    integer :: MPIerror
    type(gto_info_type), pointer :: gto => null()

    ! -------------------------------------------------------------------------

    if ( Nodes == 1 ) return

    ! hfx_gtos must be allocated on non-zero nodes
    if ( Node /= 0 ) then
      if ( allocated(hfx_gtos) ) deallocate(hfx_gtos)
      allocate(hfx_gtos(nspecies))
    end if

    ! HFX calculations require nco and nso from atm_types
    call MPI_Bcast(nco, l_max+1, MPI_integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast(nso, l_max+1, MPI_integer, 0, MPI_Comm_World, MPIerror)

    do isp=1,nspecies

      gto => hfx_gtos(isp)

      ! Number of contracted GTOs
      call MPI_Bcast(gto%orbnl_contract, maxn_orbnl, MPI_integer, 0, &
&       MPI_Comm_World, MPIerror)

      ! Index of spherical harmonic GTOs
      call MPI_Bcast(gto%orbnl_index_sphi, maxn_orbnl, MPI_integer, 0, &
&       MPI_Comm_World, MPIerror)

      ! Two-dimensional data broadcasts: exponents and coefficients for the
      ! spherical harmonic GTOs
      call MPI_Bcast(gto%orbnl_zeta(1,1), maxn_contract*maxn_orbnl, &
&       MPI_double_precision, 0, MPI_Comm_World, MPIerror)
      call MPI_Bcast(gto%orbnl_coefficient(1,1), maxn_contract*maxn_orbnl, &
&       MPI_double_precision, 0, MPI_Comm_World, MPIerror)

      ! Other variables
      call MPI_Bcast(gto%sphi(1,1), ncon_max*maxnorbs, MPI_double_precision, &
&       0, MPI_Comm_World, MPIerror)
      call MPI_Bcast(gto%pgf_radius(1,1), maxn_contract*maxn_orbnl, &
&       MPI_double_precision, 0, MPI_Comm_World, MPIerror)
      call MPI_Bcast(gto%shell_radius(1), maxn_orbnl, MPI_double_precision, &
&       0, MPI_Comm_World, MPIerror)
      call MPI_Bcast(gto%kind_radius, 1, MPI_double_precision, 0, &
&       MPI_Comm_World, MPIerror)
      call MPI_Bcast(gto%orbnl_contraction_coeff, maxn_orbnl, &
&       MPI_double_precision, 0, MPI_Comm_World, MPIerror)

    end do

    ! HFX options must be made available to all nodes
    call MPI_Bcast(hfx_options%DM_trunc, 1, MPI_LOGICAL, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%dump_fit_data, 1, MPI_LOGICAL, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%farfield, 1, MPI_LOGICAL, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%npts_fit, 1, MPI_INTEGER, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%potential_type, 1, MPI_INTEGER, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%omega, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%cutoff_radius, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%eps_farfield, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%eps_pairlist, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%eps_schwarz, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)
    call MPI_Bcast(hfx_options%eps_stored, 1, MPI_DOUBLE_PRECISION, 0, &
      MPI_Comm_World, MPIerror)

#endif

  end subroutine Broadcast_NAO2GTO

  !  Calculates the number of orbitals stored on the local Node.
  !
  !  Julian Gale, October 1998
  !
  !  Input :
  !
  !  integer NOrb     = The total number of orbitals in the calculation
  !  integer Node     = The local processor
  !  integer Nodes    = The total number of processors
  !
  !  Output :
  !
  !  integer NOrbNode = The number of orbitals stored on this Node - if
  !  zero
  !                     on input then calculated otherwise left unchanged
  !
  !  Modified so that value from last call is saved in case it can be
  !  re-used on the next call to save re-calculation.
  !
  subroutine GetNodeOrbs_NAO2GTO(NOrb, Node, Nodes, NOrbNode)

    use domain_decom
    use parallel, only: thisNode => Node
    use spatial, only: lspatial, nOrbPerNode

    implicit none

    ! Arguments
    integer, intent(in) :: NOrb, Node, Nodes
    integer, intent(out) :: NOrbNode

    ! Local variables
    integer Remainder, MinPerNode, RemainderBlocks, NOrbLast, &
&     NOrbNodeLast
    save NOrbLast, NOrbNodeLast

    if (NOrb .eq. NOrbLast) then

      ! Values are the same as last call - no need to recalculate
      NOrbNode = NOrbNodeLast

    else

      if (lspatial) then
        !-------------------------
        !  Spatial distribution  -
        !-------------------------
        NOrbNode = nOrbPerNode(Node+1)
      else if (use_dd) then
        if (use_dd_perm) then
          !  NOrbNode = dd_nuo
          ! We have this structure also
          NOrbNode = nOrbPerNode(Node+1)
        else
          if (thisNode /= Node) call die("Wrong use of dd GetNodeOrbs")
          NOrbNode = dd_nuo
        endif
      else
        !-----------------------------
        !  Blockcyclic distribution  -
        !-----------------------------
        ! Calculate the minimum number of orbitals per node
        MinPerNode = NOrb / (Nodes*nao2gto_bsize)

        ! Find the remainder of unassigned orbitals
        Remainder = NOrb - MinPerNode * Nodes * nao2gto_bsize

        ! Find number of complete blocks in the remainder
        RemainderBlocks = Remainder / nao2gto_bsize
        Remainder = Remainder - RemainderBlocks * nao2gto_bsize

        ! Workout the local number of orbitals
        NOrbNode = MinPerNode*nao2gto_bsize
        if (Node.lt.RemainderBlocks) NOrbNode = NOrbNode + nao2gto_bsize
        if (Node.eq.RemainderBlocks) NOrbNode = NOrbNode + Remainder
      endif

      ! Save value for next call
      NOrbNodeLast = NOrbNode

    endif   ! NOrb == NOrbLast

  end subroutine GetNodeOrbs_NAO2GTO

  !  Converts an orbital index in the local frame to the global frame
  !
  !  Julian Gale, Imperial College, December 1998
  !
  !  Input :
  !
  !  integer LOrb   = local orbital index
  !  integer Node   = local processor number
  !  integer Nodes  = global number of processors
  !
  !  From parallel.h :
  !
  !  integer BlockSize = blocking size for orbital distribution across
  !                      the nodes. Choice of value affects the
  !                      performance of the Scalapack routines
  !
  !  Output :
  !
  !  integer GOrb   = global orbital index
  !
  subroutine LocalToGlobalOrb_NAO2GTO(LOrb, Node, Nodes, GOrb)

    use domain_decom
    use parallel, only: thisNode => Node
    use spatial, only: lspatial, nL2G

    implicit none

    ! Arguments
    integer, intent(in) :: LOrb, Node, Nodes
    integer, intent(out) :: GOrb

    ! Local variables
    integer :: LEle, LBlock

    if (lspatial) then

      !-------------------------
      !  Spatial distribution  -
      !-------------------------
      GOrb = nL2G(LOrb, Node+1)

    else if (use_dd) then

      if (use_dd_perm) then
        !  GOrb = dd_invp(LOrb)
        ! We have this structures also...
        GOrb = nL2G(LOrb, Node+1)
      else
        if (thisNode /= Node) call die("Wrong use of dd LocalToGlobalOrb")
        GOrb = LOrb + llimit - 1
      endif

    else

      !-----------------------------
      !  Blockcyclic distribution  -
      !-----------------------------
      !  Find local block number
      LBlock = ((LOrb -1)/nao2gto_bsize)

      ! Substract local base line to find element number within the block
      LEle = LOrb - LBlock*nao2gto_bsize

      ! Calculate global index
      GOrb = (LBlock*Nodes + Node)*nao2gto_bsize + LEle

    endif

  end subroutine LocalToGlobalOrb_NAO2GTO

  !
  ! Finds a sensible default value for the blocksize default.
  ! When the number of orbitals is less than the blocksize
  ! typically used then lower the blocksize to ensure that
  ! some work is done on all nodes.
  !
  ! Input :
  !
  ! integer Nodes        : total number of processors
  ! integer nuotot       : total number of orbitals
  !
  ! Output :
  !
  ! integer blocksizedefault : default value of blocksize
  !
  ! Written by Julian Gale, March 2001
  !
  subroutine set_bsize_NAO2GTO(Nodes, pair_number, bsizedefault)

    implicit none

    ! Arguments
    integer, intent(in) :: Nodes, pair_number
    integer, intent(out) :: bsizedefault

    bsizedefault = ( (pair_number - 1) / Nodes) + 1

  end subroutine set_bsize_NAO2GTO

end module nao2gto_parallel
