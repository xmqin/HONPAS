module m_supercell

! Calculates the supercell factors by examining all atoms.
  
  use precision, only : dp
  
  implicit none
  
  public :: exact_sc_ag
  private
  
contains

  !> Calculate the exact super-cell size based on the
  !> atomic connection graph. The logic is similar to that in [[hsparse]],
  !> but instead of orbitals we deal with atoms.
  !>
  !> Two atoms are connected if any of their orbitals overlap, including
  !> the possibility of third-center KB projectors mediating an indirect overlap.
  !> For a pictorial representation of this graph,
  !> see [this page](|page|/implementation/1-auxiliary-supercell.html)
  !>
  !> Since this routine calculates explicitly the atom connections,
  !> the resulting super-cell size will be exact.
  
  subroutine exact_sc_ag(negl,ucell,na_u,isa,xa,nsc)

    use class_iSpData2D
    
    use atom_graph

#ifdef MPI
    use class_OrbitalDistribution
    use parallel, only : Nodes
    use parallelsubs, only : set_blocksizedefault
    use mpi_siesta
#endif

    !> Whether the KB projectors are neglected or not
    logical, intent(in) :: negl
    !> Unit-cell vectors
    real(dp), intent(in) :: ucell(3,3)
    !> Total number of atoms in the unit cell
    integer, intent(in) :: na_u
    !> Atomic species array
    integer, intent(in) :: isa(na_u)
    !> Atomic coordinates
    real(dp), intent(in) :: xa(3,na_u)
    !> Supercell multipliers along the three unit-cell vectors
    integer, intent(out) :: nsc(3)

    ! Calculate the atomic graph
    type(tAtomGraph) :: ag
#ifdef MPI
    type(OrbitalDistribution) :: dit
#endif
    integer, pointer :: sc(:,:)
    integer :: ia, na_l, n_nzs

#ifdef MPI
    if ( na_u > Nodes ) then ! We can only distribute something if all have something
       call set_blocksizedefault(Nodes,na_u,ia)
       call newDistribution(ia,MPI_Comm_World,dit,name='AG-dist')
       na_l = num_local_elements(dit,na_u)
       call atom_graph_generate( negl, ucell, na_u, isa, xa, ag, dit , &
            set_xijo = .false. )
    else
       call atom_graph_generate( negl, ucell, na_u, isa, xa, ag , &
            set_xijo = .false. )
    end if
#else
    na_l = na_u
    call atom_graph_generate( negl, ucell, na_u, isa, xa, ag , &
         set_xijo = .false. )
#endif

    ! Extract the largest super-cell index from the sc variable.
    n_nzs = nnzs(ag%sc_2D)
    sc => val(ag%sc_2D)

    !> Find the biggest supercell
    !> We make it symmetric on both sides.
    nsc(:) = 0
    do ia = 1 , n_nzs
       if ( abs(sc(1,ia)) > nsc(1) ) nsc(1) = abs(sc(1,ia))
       if ( abs(sc(2,ia)) > nsc(2) ) nsc(2) = abs(sc(2,ia))
       if ( abs(sc(3,ia)) > nsc(3) ) nsc(3) = abs(sc(3,ia))
    end do

    ! DEBUG
    !call atom_graph_print(ag,na_u,isa,xa)

    ! Clean-up
    call delete(ag)

#ifdef MPI
    ! Reduce the nsc
    nullify(sc) ; allocate(sc(3,1))
    call MPI_AllReduce(nsc(1),sc(1,1),3,MPI_Integer,MPI_MAX, &
         MPI_Comm_World, ia)
    nsc(:) = sc(:,1)
    deallocate(sc)
    call delete(dit)
#endif

    ! nsc now contains the one-sided range of image-cell neighbors in each of the
    ! unit-cell directions.
    ! Hence, the correct number of super-cells is:
    nsc(:) = nsc(:) * 2 + 1

  end subroutine exact_sc_ag
  
end module m_supercell
