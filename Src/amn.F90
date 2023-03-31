! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine amn( ispin )
!
!     In this subroutine we compute the overlaps between Bloch states
!     onto trial localized orbitals
!
!     These matrices will depend on three indices (see paragraph after Eq. (16) 
!     of the review by
!     N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
!     or Eq. (1.8) of the Wannier Users Guide, Version 1.2
!     $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vbar g_{n}\rangle
!     where $m$ runs between 1 and the number of bands for wannierization,
!     $n$ runs between 1 and the number of projection functions, 
!     $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
!     will be computed. 
!
!     OUTPUT: 
!     File called seedname.amn, where the overlap matrices are written in 
!     the format required by Wannier90
!   
!     BEHAVIOUR:
!     The coefficients of the wave functions are introduced
!     in amn through the matrix coeffs, computed in diagonalizeHk
!     
!     PARALLELIZATION: 
!     This subroutine is parallelized on the projectors, so each node
!     takes care on a projection.     
!
!     Implemented by J. Junquera and R. Korytar, July 2013
!

  use precision,          only: dp                  ! Real double precision type
  use parallel,           only: Nodes               ! Total number of Nodes
  use parallel,           only: Node                ! Local Node
  use parallel,           only: IONode              ! Input/output node
  use atomlist,           only: rmaxo               ! Max. cutoff atomic orbital
  use siesta_geom,        only: scell               ! Lattice vector of the
                                                    !   supercell in real space
  use siesta_geom,        only: na_s                ! Number of atoms in the
                                                    !   supercell
  use siesta_geom,        only: xa                  ! Atomic positions
  use m_siesta2wannier90, only: latvec              ! Lattice vectors in real 
                                                    !   space
  use m_siesta2wannier90, only: numproj             ! Total number of projectors
  use m_siesta2wannier90, only: numkpoints          ! Total number of k-points
                                                    !   for which the overlap of
                                                    !   the periodic part of the
                                                    !   wavefunct with a 
                                                    !   neighbour k-point will
                                                    !   be computed
  use m_siesta2wannier90, only: kpointsfrac         ! List of k points relative
                                                    !   to the reciprocal 
                                                    !   lattice vectors.
                                                    !   First  index: component
                                                    !   Second index: k-point  
                                                    !      index in the list
  use m_siesta2wannier90, only: projections         ! Trial projection functions
  use m_siesta2wannier90, only: numincbands         ! Number of bands for 
                                                    !   wannierization
                                                    !   after excluding bands  
  use m_siesta2wannier90, only: nincbands_loc       ! Number of bands for 
                                                    !   wannierization
                                                    !   after excluding bands  
                                                    !   in the local Node
  use m_siesta2wannier90, only: coeffs              ! Coefficients of the
                                                    !   wavefunctions.
                                                    !   First  index: orbital
                                                    !   Second index: band
                                                    !   Third  index: k-point
  use m_siesta2wannier90, only: Amnmat              ! Matrix of the overlaps of 
                                                    !   trial projector funtions
                                                    !   with Eigenstates of the
                                                    !   Hamiltonian
  use trialorbitalclass,  only: trialorbital        ! Derived type to define the
                                                    !    localized trial
                                                    !    orbitals
  use m_matel_registry,   only: register_in_tf_pool ! Subroutine that assigns a
                                                    !    global index to the 
                                                    !    trial projection
                                                    !    functions
  use m_new_matel,        only: new_matel           ! New MATEL implementation 
                                                    !   with the global indices 
                                                    !   of the radial functions 
                                                    !   as inputs
  use atmfuncs,           only: orb_gindex          ! Subroutine that gives
                                                    !   the global index of an
                                                    !   atomic orbital
!
! Variables for the diagonalization
!
  use atomlist,           only: no_u         ! Number of orbitals in unit cell
                                             ! NOTE: When running in parallel,
                                             !   this is core independent
  use atomlist,           only: indxuo       ! Index of equivalent orbital in 
                                             !   the unit cell
!
! Allocation/Deallocation routines
!
  use alloc,              only: re_alloc     ! Reallocation routines
  use alloc,              only: de_alloc     ! Deallocation routines

  use sys,                only: die          ! Termination routine

#ifdef MPI
  use parallelsubs,         only: GetNodeOrbs    ! Calculates the number of
                                                 !   orbitals stored on the 
                                                 !   local Node.
  use m_orderbands,       only: which_band_in_node  ! Given a node and a 
                                                    !   local index,
                                                    !   this array gives the
                                                    !   global index of the band
                                                    !   stored there
  use m_orderbands,       only: sequential_index_included_bands 
                                                    ! Sequential number of the
                                                    !   bands included for
                                                    !   wannierization
                                                    !   (the bands are listed
                                                    !   in order of incremental
                                                    !   energy)

  use mpi_siesta
#endif

  implicit none

  integer, intent(in)     :: ispin               ! Spin component

  type orbitallinkedlist
    real(dp),dimension(3)           :: center
    integer                         :: specie
    integer                         :: specieindex
    integer                         :: globalindex
    type(orbitallinkedlist),pointer :: nextitem
  end type


! Local variables
  integer  :: ik              ! Counter for loop on k-points
  integer  :: iproj           ! Counter for loop on projections
  integer  :: io              ! Counter for loop on orbital 
  integer  :: iband           ! Counter for loop on bands
  integer  :: nincbands       ! Number of bands for wannierization
  logical, save :: firsttime = .true. ! First time this subroutine is called?
  integer  :: gindex          ! Global index of the trial projector function
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexorbital ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: globalindexproj ! Global index of the neighbour atomic orbital
                              !   in the list of functions that will be
                              !   evaluated in Matel
  integer  :: indexproj       ! Index of the projector 
                              !   This index runs from 1 to the total number
                              !   of projections
  real(dp) :: trialcenter(3)  ! Position where the trial function is centered
                              !   (in Bohr)
  real(dp) :: trialrcut       ! Cutoff radius of the trial function
  real(dp) :: r12(3)          ! Relative position of the trial function with
                              !   respect to the neighbour orbital
  real(dp) :: overlap         ! Overlap between the trial function and a 
                              !   given atomic orbital
  real(dp) :: gradient(3)     ! Grandient of the overlap between 
                              !   the trial function and a given atomic orbital
  real(dp) :: phase           ! Product of the k-vector with the position
                              !   where the neighbour orbital is centered
  real(dp) :: kvector(3)      ! k-point vector in the Wannier90 grid


  complex(dp), dimension(:,:), pointer :: psiloc => null() ! Coefficients of the wave
                                             !   function (in complex format)
  complex(dp) :: exponential                 ! Exponential of exp( i * phase )
  complex(dp) :: cstar                       ! Conjugate of the coefficient
                                             !   of the wave function      
  complex(dp), parameter :: iu = cmplx(0.0_dp,1.0_dp,kind=dp) ! Imaginary unit

  type(trialorbital), save :: tf             ! Projetion function

  type(orbitallinkedlist),pointer   :: item  ! Variable that contains all the
                                             !   information of the atomic
                                             !   orbitals neighbours of a given
                                             !   trial projection function
  integer, dimension(:), save, allocatable :: projector_gindex 
                                             ! Array that gives the global index
                                             !   of a projector in the list of
                                             !   functions that will be
                                             !   evaluated by MATEL
#ifdef MPI
  integer     :: iband_global                ! Global index for a band
  integer     :: iband_sequential            ! Sequential index of the band
  integer     :: MPIerror
  complex(dp), dimension(:,:), pointer :: auxloc => null()! Temporal array for the
                                             !   the global reduction of Amnmat
#endif 

  external :: timer

! Start time counter
  call timer( 'Amn', 1 )

! Allocate memory related with the overlap matrix between the trial projection
! function and the Hamiltonian eigenstate.
! These matrices will depend on three indices (see paragraph after Eq. (16) 
! of the review by
! N. Marzari et al., Review of Modern Physics 84, 1419 (2012),
! or Eq. (1.8) of the Wannier Users Guide, Version 1.2
! $A_{m n}^{(\vec{k})} = 
!    \langle \psi_{m \vec{k}} \vbar g_{n} \rangle =
!    \sum_{\mu} \sum_{lcell} c_{\mu m}^{\ast} (\vec{k}) \times
!    exp^{i \vec{k} \cdot \left( \tau_{\mu} + \vec{R}_{lcell} \right)} \times
!    \langle \phi_{\mu} (\vec{r}-\tau_{\mu} - \vec{r}_{lcell} \vbar g_{n}\rangle
! where $m$ runs between 1 and the number of bands considered for wannierization
! $n$ runs between 1 and the number of projection functions, 
! $\vec{k}$ runs over all the k-points in the first BZ where these matrices 
! will be computed 

  nincbands = numincbands( ispin )

  call re_alloc( Amnmat,         &
 &               1, nincbands,   &
 &               1, numproj,     &
 &               1, numkpoints,  &
 &               'Amnmat',       &
 &               'Amn'           )


! Allocate memory related with a local variable where the coefficients 
! of the eigenvector at the k-point will be stored
! Only nincbands are retained for wannierization, that is why the
! second argument is made equal to nincbands
  call re_alloc( psiloc, 1, no_u, 1, nincbands, 'psiloc', 'Amn' )

! Assign a global index to the trial functions
! The same global index is assigned to a given trial function in all the nodes
  if( firsttime ) then
    allocate ( projector_gindex(numproj) )
    do iproj = 1, numproj
      tf = projections( iproj )
      call register_in_tf_pool( tf, gindex )
      projector_gindex( iproj ) = gindex
    enddo
    firsttime = .false.
  endif

kpoints:                 &
  do ik = 1, numkpoints
!   Compute the wave vector in bohr^-1 for every vector in the list
!   (done in the subroutine getkvector).
!   Remember that kpointsfrac are read from the .nnkp file in reduced units, 
!   so we have to multiply then by the reciprocal lattice vector.
    call getkvector( kpointsfrac(:,ik), kvector )

!   Initialize the local coefficient matrix for every k-point
    psiloc(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

#ifdef MPI
!   Store the local bands in this node on a complex variable
    do iband = 1, nincbands_loc
      iband_global = which_band_in_node(Node,iband)
      iband_sequential = sequential_index_included_bands(iband_global)

      do io = 1, no_u
        psiloc(io,iband_sequential) = coeffs(io,iband,ik)
      enddo

    enddo
!   Allocate workspace array for global reduction
    call re_alloc( auxloc, 1, no_u, 1, nincbands,   &
 &                 name='auxloc', routine='Amn' )
!   Global reduction of auxloc matrix
    auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    call MPI_AllReduce( psiloc(1,1), auxloc(1,1),   &
 &                      no_u*nincbands,             &
 &                      MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror )
!   After this reduction, all the nodes know the coefficients of the
!   wave function for the point ik, for all the bands and for all atomic
!   orbitals
    psiloc(:,:) = auxloc(:,:)
#else
    do iband = 1, nincbands
      do io = 1, no_u
        psiloc(io,iband) = coeffs(io,iband,ik)
      enddo
    enddo
#endif

!   Loop on the projections that will be computed in the local node
!   It would be better if each node computes all projections for
!   all its locally stored bands, instead of some projections for
!   all the bands...  In this way we will save the globalization of
!   band data.
!
#ifdef MPI
    do iproj = 1+Node, numproj, Nodes
#else
    do iproj = 1, numproj
#endif
      indexproj = iproj

!     Find the global index of the projector in the list of radial functions
!     that will be evaluated by MATEL
      globalindexproj = projector_gindex(indexproj) 

!     Find where the trial function is centered
      trialcenter = projections(indexproj)%center
  
!     Find the cutoff radius of the trial function
      trialrcut   = projections(indexproj)%rcut

!!     For debugging
!      write(6,'(a,6i5,4f12.5)')' ik, Node, nincbands, nincbands_loc, iproj, indexproj = ',  &
! &                ik, Node, nincbands, nincbands_loc, iproj,      &
! &                indexproj, trialcenter, trialrcut 
!!     End debugging

!     Find the atomic orbitals that ovelap with our radial orbital
!     centered at trialcenter and with range trialrcut
      item => get_overlapping_orbitals( scell, rmaxo, trialrcut, na_s, &
 &                                      xa, trialcenter )

OrbitalQueue:                                                        &
      do while (associated(item))
        r12 = trialcenter - item%center
        globalindexorbital = orb_gindex( item%specie, item%specieindex )
        call new_matel('S',                & ! Compute the overlap
 &                     globalindexorbital, & ! Between orbital with globalinde
 &                     globalindexproj,    & ! And projector with globalindex
 &                     r12,                & 
 &                     overlap,            & 
 &                     gradient )

        phase = -1.0_dp * dot_product( kvector, item%center )
        exponential = exp( iu * phase )

!       Loop over occupied bands
Band_loop:                                                           &
        do iband = 1, nincbands
          cstar = conjg( psiloc(item%globalindex,iband) )
          Amnmat(iband,indexproj,ik) =      & 
 &          Amnmat(iband,indexproj,ik) +    &
 &          exponential * cstar * overlap
        enddo Band_loop  ! End loop on the bands

        item => item%nextitem
      enddo OrbitalQueue

    enddo   ! Loop on projections on the local node

!   We need all the matrix Amnmat in IOnode
!   (to dump it into a file),but the results for some of the bands might
!   be computed in other nodes.
#ifdef MPI
!   Allocate workspace array for reduction
    if (IONode) then
       call re_alloc( auxloc, 1, nincbands, 1, numproj,    &
            &                 name='auxloc', routine='Amn' )
       auxloc(:,:) = cmplx(0.0_dp,0.0_dp,kind=dp)
    endif

    call MPI_Reduce( Amnmat(1,1,ik), auxloc(1,1),                       &
 &                      nincbands*numproj,                                 &
 &                      MPI_double_complex,MPI_sum,0,MPI_Comm_World,MPIerror )
    if (IONode) then
       Amnmat(:,:,ik) = auxloc(:,:)
    endif
#endif
  enddo kpoints

! Write the Amn overlap matrices in a file, in the format required
! by Wannier90
  if( IOnode ) call writeamn( ispin )

! Deallocate some of the arrays
  call de_alloc( psiloc,  'psiloc',  'Amn' )
#ifdef MPI
  if (IOnode) call de_alloc( auxloc,  'auxloc',  'Amn' )
#endif


! End time counter
  call timer( 'Amn', 2 )

  return

  contains

  function get_overlapping_orbitals( latvec, atomrcut, trialrcut, numatoms, &
 &                                   atomcoords, trialcenter )              &
 &  result(firstitem)

    use precision,          only: dp            ! Real double precision type
    use sys,                only: die                 
    use neighbour,          only: maxnna        ! Maximum number of neighbours
    use neighbour,          only: jan           ! Atom-index of neighbours
    use neighbour,          only: r2ij          ! Squared distances to neighbors
    use neighbour,          only: xij           ! Vector from a given atom
                                                !   to its neighbours
    use neighbour,          only: mneighb       ! Subroutine to compute the
                                                !   number of neighbours
    use siesta_geom,        only: isa           ! Species index of each atom
    use atomlist,           only: lasto         ! Position of last orbital 
                                                !   of each atom
    use atomlist,           only: iphorb        ! Orbital index of each  orbital
                                                !   in its atom
    use atomlist,           only: indxuo        ! Index of equivalent orbital  
                                                !   in "u" cell
    use atmfuncs,           only: rcut          ! Function that determines the
                                                !   cutoff radius of a given
                                                !   orbital of a given specie

!   This subroutine yields a list of basis orbitals that overlap 
!   with a given trial orb.
!   This subroutine bridges siesta's mneighb() with our purposes.

    implicit none
  
!   Passing variables
    real(dp),dimension(3,3),intent(in) :: latvec      ! Lattice vectors of the 
                                                      !   supercell in 
                                                      !   real space.
    real(dp)               ,intent(in) :: atomrcut    ! Maximum cutoff radius in
                                                      !   the atomic orbital 
                                                      !   basis
    integer                ,intent(in) :: numatoms    ! Number of atoms in 
                                                      !   the supercell
    real(dp),dimension(3,numatoms),intent(in) :: atomcoords  ! Atomic positions 
    real(dp),dimension(3)  ,intent(in) :: trialcenter ! Center of the 
                                                      !  trial function 
    real(dp)               ,intent(in) :: trialrcut   ! Cutoff of the 
                                                      !  trial function 

!   Passing variables
    real(dp), dimension(:,:), allocatable,save :: coords 
!   A new array containing the coordinates of the na_s atoms in the unit cell
!   plus the position of the trial function will be required to call mneighb.
!   This new array, coords, will have (3,na_s+1) dimensions

    integer          :: jneig   ! Counter for the loop on neighbours
    integer          :: nneig   ! Number of neighbors of a given trial projector
    integer          :: atom    ! Atomic index of the neighbour
    integer          :: specie  ! Atomic species of the neighbour
    integer          :: iorb    ! Counter for loop on neighbour orbitals
    integer          :: joa     ! Index for the atomic orbital within 
                                !    a given atom
    real(dp)         :: radius  ! Radius that determine the scope of the search
    real(dp)         :: rij     ! Squared distance to neighbours

    type(orbitallinkedlist), pointer     :: firstitem
    type(orbitallinkedlist), pointer     :: newitem

! Initialize mneighb
    if( .not. allocated(coords) ) then !set up static "save" variables
      allocate( coords(3,numatoms+1) )
      coords(1:3,1:numatoms) = atomcoords(1:3,1:numatoms)
    endif

    coords(:,numatoms+1) = trialcenter(:)
    radius = trialrcut + atomrcut
    call mneighb( latvec, radius, numatoms+1, coords, 0 , 0, nneig )

!
! Look for atoms that overlap with the trial orbital
!
    call mneighb( latvec, radius, numatoms+1, coords, numatoms+1, 0, nneig )
    if (nneig.gt.maxnna)                                             &
 &    call die("amn: insufficient array shapes; see mneighb(..)")

!
! Prepare output list 
!

    firstitem => null()
    do jneig = 1, nneig
      atom = jan(jneig)
      if ( atom .gt. numatoms ) cycle !it is a trial projection function
      specie = isa(atom)              !it is an atom
      rij = dsqrt( r2ij(jneig) )
      do iorb = lasto( atom-1 )+1, lasto( atom )
        joa = iphorb( iorb )
        if ( (rcut(specie,joa) + trialrcut) .gt. rij ) then
          allocate(newitem)
          newitem%center(1:3) = xij(1:3,jneig) + trialcenter(1:3)
          newitem%specie      = specie
          newitem%specieindex = joa
          newitem%globalindex = indxuo( iorb )
          newitem%nextitem    => firstitem
          firstitem => newitem
        endif
      enddo
    enddo    ! End loop on number of neighbours
  end function get_overlapping_orbitals
end subroutine amn

