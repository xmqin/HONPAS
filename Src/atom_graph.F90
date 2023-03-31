module atom_graph
  
  use precision,     only : dp

  use class_Sparsity
  use class_lSpData1D
  use class_iSpData2D
  use class_dSpData2D

  implicit none

  private

  ! Derived type to hold connectivity information on an atom basis

  ! Each atom has a maximum "footprint", both in orbital size and KB
  ! size. To decide whether two atoms are connected it is enough to
  ! check if their largest orbitals are connected, possibly via a KB
  ! (the largest) projector on a third atom.

  ! Note that listh contains the "unit cell" atom index. The image
  ! information is stored in the "box" array.
  ! The information in box could be used to find a tight supercell.
  ! Additionally, the array s_int stores info useful to check whether
  ! the interaction is due to direct overlap ("S" interaction) or to
  ! indirect (via KB projector) information)

  ! This atom_graph structure is also useful to perform domain-decomposition
  ! of the orbitals (not implemented yet).

  ! Example of call sequence (in siesta_init, for example):

  ! use atom_graph, only: atom_graph_generate, atom_graph_print, tAtomGraph
  ! ...
  ! type(tAtomGraph)  :: ag
  ! ...
  !  call atom_graph_generate( negl, ucell, na_u, isa, xa, ag )
  !  if (node==0)   call atom_graph_print(at_graph,na_u,isa,xa)
  !  call delete(ag)
  
  type, public :: tAtomGraph
     ! edge-based info
     type(lSpData1D) :: s_int_1D
     type(iSpData2D) :: sc_2D
     type(dSpData2D) :: xijo_2D
  end type tAtomGraph
         
  public :: atom_graph_generate, atom_graph_print

  interface delete
     module procedure delete_
  end interface delete
  public :: delete

contains

  subroutine delete_(this)
    type(tAtomGraph), intent(inout) :: this
    call delete(this%s_int_1D)
    call delete(this%sc_2D)
    call delete(this%xijo_2D)
  end subroutine delete_
  
  subroutine atom_graph_generate( negl, ucell, na_u, isa, xa, ag, dit, set_xijo)
    
    use class_OrbitalDistribution
#ifdef MPI
    use mpi_siesta, only : MPI_Comm_Self
#endif
    
    use intrinsic_missing, only : VNORM
    use atmfuncs, only : rcut, nofis, nkbfis
    use atm_types, only : nspecies, species, species_info
    use radial, only: rad_func
    use sorting
    use neighbour, only: jna=>jan, xij, r2ij, maxna => maxnna
    use neighbour, only: mneighb
    use dftu_specs, only: switch_dftu

    use sys, only : die
    use alloc, only : re_alloc, de_alloc

    integer,  intent(in) :: na_u
    integer,  intent(in) :: isa(na_u)
    real(dp), intent(in) :: ucell(3,3)  ! Unit cell
    real(dp), intent(in) :: xa(3,na_u)
    logical,  intent(in) :: negl
    type(tAtomGraph), intent(inout) :: ag

    type(OrbitalDistribution), intent(inout), target, optional :: dit
    logical, intent(in), optional :: set_xijo
    
    ! tolerance for comparing vector-coordinates
    real(dp), parameter        :: tol = 1.0d-8
    
    real(dp), allocatable :: rorbmax(:)  ! maximum ORB radius of each species
    real(dp), allocatable :: rkbmax(:)   ! maximum KB radius of each species
    real(dp), allocatable :: rdftumax(:) ! maximum DFTU radius of each species
    real(dp), allocatable :: rmax_projectors(:) ! maximum projector radius of each species
    integer :: maxnkb  = 500 ! max no. of atoms with
                             ! KB (or DFTU) projectors which 
                             ! overlap another
                             ! atom's orbitals.

    ! Local atom index and global atom index
    integer :: na_l
    integer :: lia, ia

    integer :: ikb, ind, inkb, is, isel
    integer :: ja, jnat, js, io
    integer :: ka, kna, ks, nna, nnkb, n_nzs

    type(species_info), pointer :: spp
    type(rad_func),     pointer :: pp

    real(dp) :: rci, rcj, rck, rij, rik, rjk
    real(dp) :: rmax, rmaxo, rmaxkb, rmaxdftu

    ! Local distribution used to
    type(OrbitalDistribution), pointer :: ldit
    type(Sparsity) :: sp

    integer, allocatable :: n_col(:), l_ptr(:)
    ! Pointers to the graph data structures
    integer, pointer :: l_col(:)
    logical, pointer :: s_int(:)
    integer, pointer :: sc(:,:)
    real(dp), pointer :: xijo(:,:)

    real(dp), dimension(:), pointer :: rckb => null()
    integer, dimension(:),  pointer :: index => null()
    integer, dimension(:),  pointer :: knakb => null()
    logical :: connected_s, connected_h
    logical :: third_centers ! Whether we need to consider interactions mediated
                             ! by projectors at third centers

    ! Whether the xijo array should be created
    logical :: lxijo
    
    real(dp) :: rcell(3,3)

    call delete(ag)
    lxijo = .false.
    if ( present(set_xijo) ) lxijo = set_xijo

    if ( present(dit) ) then
       ! Point to the passed distribution
       ldit => dit
    else
       ! We need to allocate a fake distribution
       nullify(ldit) ; allocate(ldit)
#ifdef MPI
       call newDistribution(na_u,MPI_Comm_Self,ldit,name='Fake dist')
#else
       call newDistribution(na_u,-1           ,ldit,name='Fake dist')
#endif
    end if

    ! Retrieve number of local elements in this distribution
    na_l = num_local_elements(ldit,na_u)
    
    call reclat(ucell,rcell,0)
    
    ! Find maximum radius of orbs and KB projectors of each specie
    allocate(rkbmax(nspecies), rorbmax(nspecies))
    allocate(rdftumax(nspecies), rmax_projectors(nspecies))
    do is = 1 , nspecies
       rorbmax(is) = 0.0_dp
       do io = 1 , nofis(is)
          rorbmax(is) = max(rorbmax(is),rcut(is,io))
       end do
       rkbmax(is) = 0.0_dp
       do ikb = 1 , nkbfis(is)
          rkbmax(is) = max(rkbmax(is),rcut(is,-ikb))
       end do
       rdftumax(is) = 0.0_dp
       if( switch_dftu ) then
          spp => species(is)
          do io = 1, spp%n_pjdftunl
             pp => spp%pjdftu(io)
             rdftumax(is) = max( rdftumax(is), pp%cutoff )
          enddo
       endif
    end do
      
    ! Find maximum range of basis orbitals and KB projectors
    rmaxo    = maxval(rorbmax (1:nspecies))
    rmaxkb   = maxval(rkbmax  (1:nspecies))
    rmaxdftu = maxval(rdftumax(1:nspecies))


    if ( negl ) then
       ! If we neglect the KB projectors in the interaction scheme
       ! we might still have to worry about the DFTU projectors
       rmax = 2._dp * (rmaxo+rmaxdftu)
       do is = 1, nspecies
          rmax_projectors(is) = rdftumax(is)
       enddo
       third_centers = switch_dftu  
    else
       rmax = 2._dp * (rmaxo+max(rmaxkb,rmaxdftu))
       do is = 1, nspecies
          rmax_projectors(is) = max(rdftumax(is),rkbmax(is))
       enddo
       third_centers = .true.
    end if

    ! Allocate local arrays that depend on parameters
    call re_alloc(knakb,1,maxnkb, name="knakb",routine="atom_graph")
    call re_alloc(rckb,1,maxnkb, name="rckb",routine="atom_graph")

    isel = 0
    ! Initialize internal data structures in neighb
    call mneighb( ucell, rmax, na_u, xa, 0, isel, nna )

    ! Initialize arrays for neighbour atoms and pointers
    allocate(n_col(na_l),l_ptr(na_l))

    ! First pass to get needed sizes
    do lia = 1 , na_l

       ! initialize number of columns
       n_col(lia) = 0

       ! get the global index
       ia = index_local_to_global(ldit,lia)

       ! Find neighbour atoms within maximum range
       call mneighb( ucell, rmax, na_u, xa, ia, isel, nna )
                                ! in case neighbor arrays have expanded
       call re_alloc(index,1,maxna,name="index",routine="atom_graph")

       ! Order neighbours in a well defined way
       call ordvec( tol, 3, nna, xij, index )
       call iorder( jna, 1, nna, index )
       call order ( r2ij, 1, nna, index )

       is  = isa(ia)
       rci = rorbmax(is)

       nnkb = 0
       if ( third_centers ) then
          do kna = 1,nna
             ka = jna(kna)
             rik = sqrt( r2ij(kna) )
             ks = isa(ka)
             ! It is only necessary to check with
             ! the *largest* projector in the third center
             rck = rmax_projectors(ks)
             if ( rci + rck > rik ) then
                call extend_projector(nnkb, kna, rck)
             end if
          end do
       end if

       ! Find atoms connected by direct overlap or
       ! through a KB projector or DFTU projector
       do jnat = 1 , nna
          connected_h = .false.
          ja = jna(jnat)
          js = isa(ja)
          rij = sqrt( r2ij(jnat) )
          rcj = rorbmax(js)
          !  Find if there is direct overlap
          if ( rci + rcj > rij ) then
             connected_h = .true.
          else
             ! Find if ja overlaps with a projector in ia's list
             do inkb = 1 , nnkb
                rck = rckb(inkb)
                kna = knakb(inkb)
                rjk = VNORM( xij(:,kna) - xij(:,jnat) )
                if ( rcj + rck > rjk ) then
                   connected_h = .true.
                   exit  ! loop over inkb
                end if
             end do
          end if
          if ( connected_h ) then
             n_col(lia) = n_col(lia) + 1
          end if
       end do

    end do

    ! Count number of non-zeroes and allocate column index
    n_nzs = sum(n_col(1:na_l))
    ! Create pointers
    l_ptr(1) = 0
    do lia = 2 , na_l
       l_ptr(lia) = l_ptr(lia-1) + n_col(lia-1)
    end do
    
    ! Create the sparsity pattern
    call newSparsity(sp,na_l,na_u,n_nzs,n_col,l_ptr, &
         name = 'Atom graph')

    ! Grab column index
    l_col => list_col(sp)

    ! Create a fake-distribution, or used the passed
    ! distribution

    ! Create the data structures
    call newlSpData1D(sp,ldit,ag%s_int_1D,name = 's_int')
    call newiSpData2D(sp,3,ldit,ag%sc_2D,name = 'sc-index', &
         sparsity_dim = 2)
    if ( lxijo ) then
       call newdSpData2D(sp,3,ldit,ag%xijo_2D,name = 'xijo', &
            sparsity_dim = 2)
    end if
    call delete(sp) ! clean up the sparsity pattern

    ! Get memory locations
    s_int => val(ag%s_int_1D)
    sc    => val(ag%sc_2D)
    if ( lxijo ) xijo => val(ag%xijo_2D)

    ! ... fill in the arrays
    do lia = 1 , na_l

       ! Reset to maintain counters
       n_col(lia) = 0
       
       ! get the global index
       ia = index_local_to_global(ldit,lia)

       ! Find neighbour atoms within maximum range
       call mneighb( ucell, rmax, na_u, xa, ia, isel, nna )
                                ! in case neighbor arrays have expanded
       call re_alloc(index,1,maxna,name="index",routine="atom_graph")

       ! Order neighbours in a well defined way
       call ordvec( tol, 3, nna, xij, index )
       call iorder( jna, 1, nna, index )
       call order ( r2ij, 1, nna, index )

       is  = isa(ia)
       rci = rorbmax(is)

       nnkb = 0
       if ( third_centers ) then
          do kna = 1,nna
             ka = jna(kna)
             rik = sqrt( r2ij(kna) )
             ks = isa(ka)
             ! It is only necessary to check with
             ! the *largest* object in the third center
             rck = rmax_projectors(ks)
             if ( rci + rck > rik ) then
                call extend_projector(nnkb, kna, rck)
             end if
          end do
       end if

       ! Find atoms connected by direct overlap or
       ! through a KB projector
       do jnat = 1 , nna
          connected_s = .false.
          connected_h = .false.
          ja = jna(jnat)
          js = isa(ja)
          rij = sqrt( r2ij(jnat) )
          rcj = rorbmax(js)
          ! Find if there is direct overlap
          if ( rci + rcj > rij ) then
             connected_s = .true.
             connected_h = .true.
          else 
             ! Find if ja overlaps with a KB/DFTU projector in ia's list
             do inkb = 1 , nnkb
                rck = rckb(inkb)
                kna = knakb(inkb)
                rjk = VNORM( xij(:,kna) - xij(:,jnat) )
                if ( rcj + rck  > rjk ) then
                   connected_h = .true.
                   exit  ! loop over inkb
                end if
             end do
          end if
          if ( connected_h ) then
             n_col(lia)    = n_col(lia) + 1
             ind           = l_ptr(lia) + n_col(lia)
             l_col(ind)    = ja
             if ( lxijo ) xijo(1:3,ind) = xij(1:3,jnat)
             sc(1:3,ind)   = reduce(xa(1:3,ia),xij(1:3,jnat), &
                  xa(1:3,ja),rcell)
             s_int(ind)    = connected_s
          end if
       end do
    end do

    ! Deallocate local arrays
    call de_alloc(index,name="index",routine="atom_graph")
    call de_alloc(knakb,name="knakb",routine="atom_graph")
    call de_alloc(rckb,name="rckb",routine="atom_graph")
    deallocate(n_col,l_ptr,rkbmax,rorbmax,rdftumax,rmax_projectors)

    ! Clean up, if the distribution has not
    ! been supplied the local distribution 
    ! needs to be deleted.
    if ( .not. present(dit) ) then
       call delete(ldit)
       deallocate(ldit)
    end if

  contains
    
    ! Possibly increase the arrays used for creating
    ! the data arrays
    ! ALSO add data to the arrays
    subroutine extend_projector(nnkb, kna, rck)
      integer, intent(inout) :: nnkb
      integer, intent(in)    :: kna
      real(dp), intent(in)   :: rck
      
      if ( nnkb == maxnkb ) then
         maxnkb = maxnkb + 10
         call re_alloc( knakb, 1, maxnkb, 'knakb', &
              'atom_graph', .true. )
         call re_alloc( rckb, 1, maxnkb, 'rckb', &
              'atom_graph', .true. )
      end if
      
      nnkb        = nnkb + 1
      knakb(nnkb) = kna
      rckb(nnkb)  = rck
      
    end subroutine extend_projector


    function reduce(xi,xij,xj,rcell) result(sc)
      real(dp), intent(in) :: xi(3), xj(3), xij(3)
      real(dp), intent(in) :: rcell(3,3)
      integer :: sc(3)
      
      real(dp) :: y(3)
      
!      print "(a,2(2x,3f10.4))", "xi, xij: ", xi, xij
!      print "(a,2(2x,3f10.4))", "xi+xij, xj: ", xi+xij, xj
      y = xi + xij - xj
!      print "(a,(2x,3f10.4))", "y: ", y
      y = matmul(y,rcell)
!      print "(a,(2x,3f10.4))", "lattice y: ", y
      sc = nint(y)
!      print "(a,3i4)", "sc: ", sc
    end function reduce
    
  end subroutine atom_graph_generate
  
  subroutine atom_graph_print(ag,na_u,isa,xa)
    
    use class_OrbitalDistribution
    
    type(tAtomGraph), intent(in) :: ag
    integer,  intent(in)    :: na_u
    integer,  intent(in)    :: isa(na_u)
    real(dp), intent(in)    :: xa(3,na_u)

    type(OrbitalDistribution), pointer :: dit
    type(Sparsity), pointer :: sp
    ! Pointers to the graph data structures
    integer, pointer :: n_col(:), l_ptr(:), l_col(:)
    logical, pointer :: s_int(:)
    integer, pointer :: sc(:,:)

    integer :: na_l, lia, ia, n_nzs, ind, ja, j

    ! The distribution
    dit => dist(ag%s_int_1D)
    sp  => spar(ag%s_int_1D)
    call attach(sp, n_col=n_col,list_ptr=l_ptr,list_col=l_col, &
         nrows = na_l , nnzs = n_nzs)

    s_int => val(ag%s_int_1D)
    sc   => val(ag%sc_2D)

    write(*,"(a,3i8)") "Atom graph: na_l, nnzs, Node: ", &
         na_l, n_nzs, dist_node(dit)

    do lia = 1 , na_l
       ia = index_local_to_global(dit,lia)

       write(*,"(/,a,i3,a,i2,i4)") "--Neighbors of atom ", &
            ia, " spec: ", isa(ia), n_col(lia)
       do j = 1 , n_col(lia)
          ind = l_ptr(lia) + j
          ja = l_col(ind)
          write(*,fmt="(4x,i3,a,i2,a,3x,3i4,1x,l1)") &
               ja, " spec: ", isa(ja), &
               " sc: ", sc(1:3,ind), s_int(ind)
       end do
    end do

  end subroutine atom_graph_print

end module atom_graph





