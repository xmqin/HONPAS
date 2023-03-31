!#define PVT_DEBUG
!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! This module contains all the different methods by 
! which the user can pivot the sparsity pattern.

module m_pivot_methods

  use precision, only : i8b

  use m_region

  implicit none

  ! The public routines
  public :: Cuthill_Mckee, rev_Cuthill_Mckee
  public :: GPS, rev_GPS
  public :: GGPS, rev_GGPS
  public :: PCG, rev_PCG
  public :: connectivity_graph
  public :: rev_connectivity_graph

#ifdef SIESTA__METIS
  public :: metis_PartGraphKway_pvt
  public :: metis_PartGraphRecursive_pvt
  public :: metis_NodeND_pvt
#endif

  public :: bandwidth, profile

  public :: sp2graphviz
  interface sp2graphviz
    module procedure sp2graphviz_sp
    module procedure sp2graphviz_lists
  end interface sp2graphviz

  private

  ! Degree determining parameters.
  ! Used for input to the idx_degree function.
  !   find lowest degree
  integer, parameter :: D_LOW = 1 
  !   find highest degree
  integer, parameter :: D_HIGH = 2
  !   find lowest sum of neighbor degree
  integer, parameter :: D_LOW_SUM = 3
  !   find highest sum of neighbor degree
  integer, parameter :: D_HIGH_SUM = 4

  ! Connect determining parameters.
  ! Used for input to the sort_connect function.
  !   find the index that has the most connecting.
  integer, parameter :: C_HIGH = 1 
  !   find the index that has the least connecting.
  integer, parameter :: C_LOW = 2

  ! Internal data-type to generate a dependency scheme of the
  ! pivot arrays
  type :: tPvtLvl
    ! The pivot table
    type(tRgn) :: pvt
    ! The level that each entry belongs to
    type(tRgn) :: lvl 
  end type tPvtLvl

  ! Linked list of PvtLvl's
  type :: tllPvtLvl
    ! Starting linked list level
    type(tPvtLvl) :: v
    ! the following level
    type(tllPvtLvl), pointer :: next => null()
  end type tllPvtLvl

  ! Linked list of levels
  ! This should be preferred over tPvtLvl as it has the
  ! same information at half the memory requirements (and probably
  ! also faster).
  type :: tLevelStructure
    ! Level number
    integer :: lvl = 0
    ! The elements in the level
    type(tRgn) :: v
    ! the next level level
    type(tLevelStructure), pointer :: next => null()
  end type tLevelStructure

#ifdef PVT_DEBUG
  character(len=*), parameter, private :: fmt_debug1_a = '(tr5,a)'
  character(len=*), parameter, private :: fmt_debug1_ai = '(tr5,a,1(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug1_aii = '(tr5,a,2(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug1_aiii = '(tr5,a,3(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug2_a = '(tr7,a)'
  character(len=*), parameter, private :: fmt_debug2_ai = '(tr7,a,1(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug2_aii = '(tr7,a,2(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug2_aiii = '(tr7,a,3(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug3_a = '(tr9,a)'
  character(len=*), parameter, private :: fmt_debug3_ai = '(tr9,a,1(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug3_aii = '(tr9,a,2(tr1,i8))'
  character(len=*), parameter, private :: fmt_debug3_aiii = '(tr9,a,3(tr1,i8))'
#endif

contains

  ! The Generalized GPS algorithm 
  ! (Progress In Electromagnetics Research, PIER 90, 121–136, 2009)
  subroutine GGPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt,priority, &
      range)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    ! The allowed range to save the node point in the new set (step II-(a))
    integer, intent(in), optional :: range

    ! The local allowed range
    integer :: lrange

    ! Temporary region used to contain the connectivity graph
    type(tRgn) :: con
    ! The level structure, contains the level of each
    type(tRgn) :: lvl

    ! The last region of the currently investigated level
    type(tRgn) :: S
    type(tRgn) :: r ! temporary region
    type(tRgn) :: renum

    ! local variables
    integer :: deg, N_lvl, depth, width, N_sub, etr_small
    integer :: i, il, j, k, etr, idx, tmp(3), target_N
    logical :: III_a, III_b

    ! The level structures sorted in linked-lists format
    integer :: N_v, N_u ! denoted 'p' and 'q' in the article, respectively
    type(tllPvtLvl), target :: vs, us
    type(tllPvtLvl), pointer :: ps
    ! Linked list of regions
    type(tRgnLL), target :: rll_list
    type(tRgnLL), pointer :: rll, rll2

    ! Sets used to contain level sets
    integer, allocatable :: lvl_set(:), h(:,:)
    logical :: suc

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'GGPS START'
#endif

    call rgn_delete(pvt)

    lrange = 0
    if ( present(range) ) lrange = range

    III_a = .false.
    III_b = .false.

    ! Find a set of pseudo-peripherals using the GPS algorithm
    call pseudo_peripheral(D_LOW,n,nnzs,n_col,l_ptr,l_col,sub,vs%v, &
        small_wd = tmp, priority = priority)
    width = tmp(1) ! the s \in L_ec where the width is the smallest
    depth = tmp(2) ! the depth of s
    etr_small = tmp(3)

#ifdef PVT_DEBUG
    write(*,*)'   primary set level [depth / width]:       ', &
        lvl_depth(vs%v%lvl),lvl_width(vs%v%lvl)
    write(*,*)'   s [min width] set level [depth / width/ etr_small]: ',depth,width,etr_small
#endif

    ! I cannot figure out (iv), it seems they should say 
    !    "has the largest depth", but I am not sure...
    !    I assume that they mean the smallest width, (width/depth above)

    ! We have now created the first 'v' point
    !   (v) -- we need to populate 'u' end points
    ps => us
    ! Get depth of current level structure (of the 'v' end)
    N_lvl = lvl_depth(vs%v%lvl)

    ! extract the nodes at the last level to search for all 
    call lvl_struct_extract(vs%v%pvt,vs%v%lvl,N_lvl,S)

    ! *** 
    !    Change, we sort according to the highest degree
    ! ***
    call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,S,r)
    S%r = rgn_pivot(sub,r%r)

    ! For the case that the first minimum width, maximum depth
    ! has not been added to the 'u' set, we add it
    if ( .not. in_rgn(S,etr_small) ) then
#ifdef PVT_DEBUG
      write(*,*)'   s [min width] set forcefully added to u end'
#endif
      call level_struct(n,nnzs,n_col,l_ptr,l_col,etr_small,sub, &
          ps%v%pvt,ps%v%lvl, priority)
      N_sub = 1
    else
      N_sub = 0
    end if
    do i = 1 , S%n
      etr = S%r(i)

#ifdef PVT_DEBUG
      write(*,'(a,i0,a,i0,a,i0)', advance='NO')' checking level structure: ',i, ' / ', S%n, &
          ' with starting index: ', etr
#endif
      call level_struct(n,nnzs,n_col,l_ptr,l_col,etr,sub,pvt,lvl,priority)
      k = lvl_depth(lvl)

      ! If the level depth is larger than the current start over 
      ! from that index (the last added element *must* be
      ! the largest level)
      if ( N_lvl < k ) then
        call die('We have apparently not found the perfect peripheral v')
      else if ( N_lvl == k .or. depth == k ) then
        ! allocate 'u' and save
        ! figure out if we should step
        if ( ps%v%pvt%n > 0 ) then
          allocate(ps%next)
          ps => ps%next
        end if
        call rgn_copy(pvt,ps%v%pvt)
        call rgn_copy(lvl,ps%v%lvl)
        N_sub = N_sub + 1
#ifdef PVT_DEBUG
        write(*,'(a)') ' ... STORED'
      else
        write(*,'(a)') ' ... DISCARDED'
#endif
      end if

    end do

#ifdef PVT_DEBUG
    write(*,*)' u set has peripherals ',N_sub
#endif

    !   (vi) --- search for other 'v' end points
    ! first populate a list with all nodes with same degree
    ! as the first 'v'

    ! To ease the search for previous found sets, we 
    ! use a list to track the already found points
    call rgn_init(con,sub%n)
    con%n = 0
    suc = rgn_push(con,vs%v%pvt%r(1)) ! push 'v' (we know it only has one)
    ps => us
    suc = rgn_push(con,ps%v%pvt%r(1))
    do while ( associated(ps%next) ) 
      ps => ps%next
      suc = rgn_push(con,ps%v%pvt%r(1))
    end do
    ! *** CHANGE 
    !     To narrow the search (we do not need a huge number of sets,
    !     it will produce enough already!)
    !     we limit the search to only look at connections within the
    !     first two levels of 'v'
    !     This will obviously increase the chance that the peripheral
    !     points exhibit a similarity.
    do i = 3 , N_lvl
      call lvl_struct_extract(vs%v%pvt,vs%v%lvl,i,S)
      do j = 1 , S%n
        if ( .not. in_rgn(con,S%r(j)) ) &
            suc = rgn_push(con,S%r(j))
      end do
    end do
    ! *** END CHANGE
    call rgn_sort(con) ! makes in_rgn faster
    deg = n_col(vs%v%pvt%r(1))
    call rgn_init(S,sub%n)
    S%n = 0 ! initialize list
    do i = 1 , sub%n
      etr = sub%r(i)
      if ( n_col(etr) /= deg ) cycle ! wrong degree
      if ( in_rgn(con,etr) ) cycle ! already found
      suc = rgn_push(S,etr)
    end do
    ! search through all with same degree and pick-out those
    ! with same depth, add them to 'v'
    N_lvl = lvl_depth(vs%v%lvl)
    ps => vs
    N_sub = 1
    do i = 1 , S%n
      etr = rgn_pivot(sub,S%r(i))

      call level_struct(n,nnzs,n_col,l_ptr,l_col,etr,sub,pvt,lvl,priority)
      k = lvl_depth(lvl)

      ! If the depth is the same as found in vs%v
      ! then we can add it
      if ( N_lvl < k ) then
        call die('We have apparently not found the perfect peripheral v')
      else if ( N_lvl == k .or. depth == k ) then
        ! add it to the list
        if ( ps%v%pvt%n > 0 ) then
          allocate(ps%next)
          ps => ps%next
        end if
        call rgn_copy(pvt,ps%v%pvt)
        call rgn_copy(lvl,ps%v%lvl)
        N_sub = N_sub + 1
      end if
    end do

#ifdef PVT_DEBUG
    write(*,*)' v has peripherals ',N_sub
#endif

    ! ### Completed algorithm I. ###

    ! At this point we have two sets, 

    ! ### Algorithm II ###

    ! To ease the operation we "attach" the 'u' at the end
    ! of 'v'
    ps => vs
    N_v = 1
    do while ( associated(ps%next) )
      N_v = N_v + 1 ; ps => ps%next
    end do
    ps%next => us

    ps => us
    N_u = 1
    do while ( associated(ps%next) )
      N_u = N_u + 1 ; ps => ps%next
    end do

#ifdef PVT_DEBUG
    write(*,*)'   attached v peripherals ',N_v
    write(*,*)'   attached u peripherals ',N_u
#endif

    ! the level set that tracks the levels each one 
    ! is assigned
    allocate(lvl_set(N_v + N_u))

    ! We initialize the pivoting and lvl arrays
    call rgn_init(pvt,sub%n)
    pvt%n = 0
    call rgn_init(lvl,sub%n)
    lvl%n = 0

    ! For each node we find the lvl_set
    !   We need a sub-copy to control which have been 
    !   added/removed
    call rgn_copy(sub,S)
    do while ( S%n > 0 )

      ! Figure out the entry
      etr = rgn_pop(S)

      ! Create the level-set
      ps => vs
      do i = 1 , N_v
        idx = rgn_pivot(ps%v%pvt,etr)
        lvl_set(i) = ps%v%lvl%r(idx)
        ps => ps%next
      end do
      do i = 1 , N_u
        idx = rgn_pivot(ps%v%pvt,etr)
        lvl_set(N_v+i) = N_lvl + 1 - ps%v%lvl%r(idx)
        ps => ps%next
      end do

#ifdef PVT_DEBUG
      write(*,'(a,1000(tr1,i0))')' II-(a) -- level set etr / ',etr,lvl_set
#endif

      ! (a) -- simple case, we add all individually
      j = minval(lvl_set(:))
      if ( all(lvl_set(:) - j <= lrange) ) then
        ! add to level lvl_set(1) and remove connectivity graph
        suc = rgn_push(pvt,etr)
        suc = rgn_push(lvl,lvl_set(1))
        call rgn_list(r,n_col(etr), &
            l_col(l_ptr(etr)+1:l_ptr(etr)+n_col(etr)) )
        ! We shall only take those in the sub-region
        if ( r%n > 0 ) &
            call rgn_intersection(r,sub,r)
        ! Remove all these from S
        call rgn_complement(r,S,S)
      end if

    end do

#ifdef PVT_DEBUG
    write(*,*)' II-(a) -- cleared same peripherals ',pvt%n
#endif

    ! ### still in algorithm II ###
    ! Now pvt contains all elements from (a)
    ! We need to handle cases (b) and (c)

    ! First create the sub-graphs
    call rgn_complement(pvt,sub,S)
    ! now S contains all elements not added to
    ! the pivoting scheme, create tree
    rll => rll_list
    N_sub = 0
    do while ( S%n > 0 ) 

      etr = rgn_pop(S)

      call rgn_init(con,1,val=etr)
      i = 0
      do while ( i < con%n ) 
        ! create connectivity graph
        i = i + 1
        etr = con%r(i)
        call rgn_list(r,n_col(etr), &
            l_col(l_ptr(etr)+1:l_ptr(etr)+n_col(etr)) )
        ! We shall only take those in the sub-region
        if ( r%n > 0 ) &
            call rgn_intersection(r,sub,r)
        ! Remove those already taken
        if ( pvt%n > 0 ) &
            call rgn_complement(pvt,r,r)
        ! Extend the connectivity graph
        call rgn_union(con,r,con)
      end do

      ! Save connectivity graph
      if ( rll%rgn%n > 0 ) then
        allocate(rll%next)
        rll => rll%next
      end if
      call rgn_copy(con,rll%rgn)
      N_sub = N_sub + 1

      ! now con' is one sub-graph
      call rgn_complement(con,S,S)

    end do

#ifdef PVT_DEBUG
    write(*,*)' created sub graphs C_i ',N_sub
#endif

    ! Sort the graphs so that rll(1), rll(2) has |V(C_1)| >= |V(C_2)|, etc.
    rll => rll_list
    do il = 1 , N_sub - 1
      do i = il + 1 , N_sub
        rll2 => rll%next
        if ( rll2%rgn%n > rll%rgn%n ) then
          ! At this point we should not have too many
          ! different regions, hence we do it "stupidly".
          call rgn_copy(rll2%rgn,r)
          call rgn_copy(rll%rgn,rll2%rgn)
          call rgn_copy(r,rll%rgn)
        end if
      end do
      rll => rll%next
    end do

#ifdef PVT_DEBUG
    write(*,*)'   sorted |V(C_1)|'
    rll => rll_list
    do il = 1 , N_sub
      write(*,'(a,i0,a,i0)') '      |C_',il,'| = ',rll%rgn%n
      rll => rll%next
    end do
#endif

    ! (c) -- figure out the level placement of the node
    !    We allocate the 'n' and 'h' variables as denoted in the
    !    article
    allocate(h(0:N_lvl,N_v+N_u))
    rll => rll_list
    do il = 1 , N_sub

      h(:,:) = 0

      !  (1) -- (2), (1) is not necessary as it will be subtracted later
      ! Loop over the elements in this sub-graph
      do i = 1 , rll%rgn%n
        etr = rll%rgn%r(i)

        ! Create the level-set for the current node
        ps => vs
        do j = 1 , N_v
          idx = rgn_pivot(ps%v%pvt,etr)
          ! the j'th element of the associated level set
          k = ps%v%lvl%r(idx)
          ! Update h, add to 'h' 1 if the element
          ! currently investigating is placed in the
          ! k'th level
          h(k,j) = h(k,j) + 1
          ps => ps%next
        end do
        do j = 1 , N_u
          idx = rgn_pivot(ps%v%pvt,etr)
          k = N_lvl + 1 - ps%v%lvl%r(idx)
          h(k,N_v+j) = h(k,N_v+j) + 1
          ps => ps%next
        end do

      end do

      !  (3) -- h_0 ^ j
      ! We now have all 'h', figure out where to add
      ! the current graph points h_0
      do j = 1 , N_v + N_u
        h(0,j) = 0
        do i = 1 , N_lvl
          h(0,j) = max(h(0,j),h(i,j))
        end do
      end do
      !  (3) -- h_0
      j = 1 ! assume j to be the first one
      do i = 1 , N_v + N_u
        if ( h(0,i) < h(0,j) ) j = i
      end do
      ! j is now the index of minimum 'h'
      ! with priority to the first one (actually, priority to 'v' set)

      !  (3) -- add all elements of rll%rgn to the
      !         levels indicated by the level set 'j'
      ! We skip to the correct level set
      ps => vs
      do i = 1 , j - 1
        ps => ps%next
      end do

      ! Add all rll%rgn to the levels as indicated by
      ! ps%lvl
      do i = 1 , rll%rgn%n
        etr = rll%rgn%r(i)

        ! figure out the level
        idx = rgn_pivot(ps%v%pvt,etr)
        ! the j'th element of the associated level set
        k = ps%v%lvl%r(idx)

        ! Add it to the level
        suc = rgn_push(pvt,etr)
        suc = rgn_push(lvl,k)

      end do

      rll => rll%next

    end do

    ! We should now have added all elements from the sub-graphs

    if ( pvt%n /= sub%n ) &
        call die('Did not collect all entries in the entry-levels')

    ! #### CLEAN UP ####

#ifdef PVT_DEBUG
    write(*,*)'   -- cleaning v'
#endif

    ! De-allocate the entries
    N_sub = 1
    do i = 1 , N_v - 1
      ps => vs
      do j = i + 1 , N_v - 1 ! we delete the "next" one
        ps => ps%next
      end do
      if ( associated(ps%next) ) then
        call rgn_delete(ps%next%v%pvt)
        call rgn_delete(ps%next%v%lvl)
        deallocate(ps%next)
        N_sub = N_sub + 1 
      end if
      nullify(ps%next)
    end do
    call rgn_delete(vs%v%pvt)
    call rgn_delete(vs%v%lvl)

#ifdef PVT_DEBUG
    write(*,*)'   -- cleaned  v', N_sub
    write(*,*)'   -- cleaning u'
#endif

    ! De-allocate the entries
    N_sub = 1
    do i = 1 , N_u - 1
      ps => us
      do j = i + 1 , N_u - 1 ! we delete the "next" one
        ps => ps%next
      end do
      if ( associated(ps%next) ) then
        call rgn_delete(ps%next%v%pvt)
        call rgn_delete(ps%next%v%lvl)
        deallocate(ps%next)
        N_sub = N_sub + 1 
      end if
      nullify(ps%next)
    end do
    call rgn_delete(us%v%pvt)
    call rgn_delete(us%v%lvl)

#ifdef PVT_DEBUG
    write(*,*)'   -- cleaned  u', N_sub
    write(*,*)'   -- cleaning rll'
#endif

    ! delete linked list for the regions
    call rgnll_delete(rll_list)

    ! ### Completed algorithm II. ###


    ! Algorithm III

    ! This algorithm will re-order the elements in pvt

    !  (i) -- It seems that this 
    !         should improve on the initial description.
    !         However, the usage is not clear to me, should it
    !         be from set 'u' or from set 'v', or any node in L_1,
    !         or any node in L_N_lvl, or any node in L_1 U l_N_lvl,
    !         or any node in G? ...
    !    ** NOTE From their other paper
    !        Wang, Guo, Shi, Progress In Electromagnetics research Symposium, 
    !        Hangzhou, China, 2008
    !       It seems like we should select the new 'v' as the node with
    !       smallest degree in L_1 U L_ec

    ! Used for calculating the degree of each node
    call rgn_init(con,n)
    con%n = 0

    !   (i)

    ! Extract the depth
    N_lvl = lvl_depth(lvl)
    ! Extract the first and last level
    call lvl_struct_extract(pvt,lvl,1,S)
    call sort_degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col,S,r)
    ! Get the node with the smallest sum
    etr = r%r(1)
    deg = degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col,etr,con)
#ifdef PVT_DEBUG
    write(*,*)'   initial level 1 degree ',etr,deg
#endif
    ! now we have a target degree for the other end-points
    call lvl_struct_extract(pvt,lvl,N_lvl,S)
    ! Now S contains all 'u' end-points
    do i = 1 , S%n
      ! Get degree of this entry
      j = degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col,S%r(i),con)
      if ( j < deg ) then
        deg = j
        etr = S%r(i)
        ! this triggers the swapping of the last indices
        III_a = .true.
      end if
    end do
    ! We should now have found the end-point with the smallest 
    ! degree (according to the LOW_SUM method)

#ifdef PVT_DEBUG
    write(*,*)'   smallest degree pseudo-peripheral ',etr,deg
#endif

    if ( etr /= r%r(1) ) then
      ! We need to swap the level structure
      do i = 1 , lvl%n
        lvl%r(i) = N_lvl + 1 - lvl%r(i)
      end do
    end if

    ! This below loop handles (ii) and (iii) simultaneously
    call rgn_init(renum,sub%n)
    renum%n = 0
    do il = 1 , N_lvl

      ! Create skipping region (all but L_il)
      call rgn_range(con,1,n) ! all
      call lvl_struct_extract(pvt,lvl,il,r) ! L_il
      call rgn_complement(r,con,con)
      ! We use a consecutive array to hold the skipping region
      call rgn_init_consecutive(n, S, con)

      ! Make sure there is space
      call rgn_init(con, sub%n)

      if ( il == 1 ) then

        ! Initialize renum
        suc = rgn_push(renum,etr)
        call rgn_consecutive_insert(S, etr)
        k     = 1
        N_sub = 1
        !  (ii) -- assign integers to the nodes in L_1
        !  (ii) -- (a-c)
        !    This is numbering the nodes in L_1
        target_N = lvl_width(lvl,il)

#ifdef PVT_DEBUG
        write(*,*)'   III (ii) start entry',etr
#endif
      else

        k = N_sub + 1 ! we assume that we have the following one
        if ( k > renum%n ) call die('GGPS: Connectivity across levels are &
            &non-existing, you must have a very special matrix.')
        ! We know that all entries in renum are ordered
        ! and their index is their associated number
        ! hence we make 'k' follow the index of the first
        ! entry for level il-1
        ! find ending index of the previous level
        N_sub = k + lvl_width(lvl,il-1) - 1

        target_N = N_sub + lvl_width(lvl,il)

      end if

#ifdef PVT_DEBUG
      write(*,*)'   III (ii-iii) current n',il, renum%n
      write(*,*)'      level width',lvl_width(lvl,il)
      write(*,*)'      k,N_sub,target ',k,N_sub,target_N
#endif

      ! Loop over all entries in the previous level
      ! As N_sub might increase, we cannot use a
      ! fixed do-loop
      i = k - 1
      do while ( renum%n < target_N ) 
        i = i + 1

        ! 1. Create connectivity graph
        if ( i <= renum%n ) then
          ! Ensure that it is clear for graph_connect (needs to be
          ! allocated)
          call rgn_init(con,sub%n)
          con%n = 0
          call graph_connect(renum%r(i),n,nnzs,n_col,l_ptr,l_col,con, &
              skip = S)
        end if

        ! 2. correction step
        ! if the current level is one and there is no connectivity
        ! we need to fix it.
        if ( N_sub < i ) then ! this *MUST* ensure i > renum%n
          ! We need to pick a new one
          ! We take the one in the current level which have not
          ! been taken yet
          call lvl_struct_extract(pvt,lvl,il,r)
          ! Remove all entries already taken
          call rgn_complement(renum,r,r)
          ! Sort the graph with respect to the LOW_SUM ( we will only
          ! add ONE element)
          call sort_degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col,r,con)
          if ( con%n > 0 ) con%n = 1
          ! now r contains a new list of entries for the current
          ! level, by allowing it to run through we will allow
          ! it to be added as noted in the algorithm
#ifdef PVT_DEBUG
          write(*,*)'   -- non-connected in lvl',il,k,i,N_sub
          write(*,*)'   -- adding non-connected entry in lvl',il,i,renum%n
#endif
          ! Ensure that the one we add will be searched for
          ! connections, note that we add one below, hence, this is
          ! will step to the one just added
          i = renum%n
        end if
        if ( con%n == 0 ) cycle
        ! Sort the graph with respect to the LOW_SUM
        call sort_degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col,con,r)
        ! add all adjacent nodes
        if ( .not. rgn_push(renum, r) ) call die('Error GGPS -- push 2')
        call rgn_consecutive_insert(S, r)

        ! if we are in the first level, we
        ! need to increase the search space
        ! to those as well (N_sub will then be renum%n)
        if ( il == 1 ) N_sub = N_sub + r%n

        ! Quick escape if the number of elements already match N_sub
        if ( renum%n == target_N ) exit

      end do

      ! Force the next loop to start from L_1
      if ( il == 1 ) N_sub = 0

    end do
    call rgn_delete(con,r,S)

    ! Check that we have collected all
    if ( renum%n /= sub%n ) then
      print *,renum%n,sub%n
      call die('GGPS: Doing Algo III (i)-(iii) we could not assert total &
          &connectivity graph.')
    end if

    !   (iv) -- create the pivot-table and create correct numbering
    if ( III_a ) then
#ifdef PVT_DEBUG
      write(*,*)'   -- swapping indices due to III(a)'
#endif
      k = renum%n
      do i = 1 , k / 2
        j = renum%r(i)
        renum%r(i) = renum%r(k+1-i)
        renum%r(k+1-i) = j
      end do
    end if

    call rgn_copy(renum,pvt)

    ! Clean-up 
    deallocate(lvl_set,h)
    call rgn_delete(con,lvl,S)
    call rgn_delete(r,renum)

    if ( pvt%n /= sub%n ) &
        call die('GGPS: Error in algorithm')

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'GGPS DONE'
#endif

  end subroutine GGPS

  ! The Generalized GPS algorithm 
  ! (Progress In Electromagnetics Research, PIER 90, 121–136, 2009)
  subroutine GGPS_new(n,nnzs,n_col,l_ptr,l_col,sub,pvt,priority, &
      range)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    ! The allowed range to save the node point in the new set (step II-(a))
    integer, intent(in), optional :: range

    ! The local allowed range
    integer :: lrange

    ! Temporary region used to contain the connectivity graph
    type(tRgn) :: con
    ! The level structure, contains the level of each
    type(tRgn) :: set, lvl

    ! The last region of the currently investigated level
    type(tRgn) :: S, skip
    type(tRgn) :: r ! temporary region
    type(tRgn) :: renum

    ! local variables
    integer :: depth, width
    integer :: k, iuv, vlvl, ulvl
    integer :: depthv, widthv, v, degreev
    integer :: depthu, widthu, u
    integer :: i, il, j, etr, idx, tmp(3)
    logical :: III_a, III_b

    type(tRgn) :: vend, uend
    integer :: n_uss, n_vss
    type(tLevelStructure) :: vs, us, xs, C
    type(tLevelStructure), allocatable :: uss(:)
    type(tLevelStructure), allocatable :: vss(:)
    type(tLevelStructure), pointer :: ps

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'GGPS START'
#endif

    call rgn_delete(pvt)

    lrange = 0
    if ( present(range) ) lrange = range

    III_a = .false.
    III_b = .false.

    ! Create skip
    call rgn_range(r, 1, n)
    call rgn_complement(sub, r, skip)
    call rgn_delete(r)

    ! Find a set of pseudo-peripherals using the GPS algorithm
    ! This will retrieve set v and width, depth, etr_small are for node v
    ! (i) -- (iii)
    call pseudo_peripheral_level_structure(D_LOW,n,nnzs,n_col,l_ptr,l_col,sub,vs, &
        small_wd = tmp, priority = priority)
    widthv = tmp(1) ! the s \in L_ec where the width is the smallest
    depthv = tmp(2) ! the depth of s
    v = tmp(3)

#ifdef PVT_DEBUG
    write(*,*)'   primary set level [depth / width]: ', &
        level_structure_depth(vs), level_structure_width(vs)
    write(*,*)'   s [min width] set level [depth / width / etr_small]: ',depthv,widthv,v
#endif

    ! (iv), find smallest width in level-structures of S (i.e. last level of L_v)
    depth = level_structure_depth(vs)
    call level_structure_level(vs, depth, S)
    widthu = huge(1)
    do i = 1 , S%n
      call rgn_list(r, 1, S%r(i:i))

      call level_structure(n,nnzs,n_col,l_ptr,l_col,us,r,skip,priority=priority)
      j = level_structure_width(us)
      if ( j < widthu ) then
        widthu = j
        depthu = level_structure_depth(us)
        u = S%r(i)
#ifdef PVT_DEBUG
        write(*,*)'   new minumim width [width / depth / etr_small] : ', &
            widthu, depthu, u
#endif
      end if

    end do


    ! this is the current depth we are targetting
    k = max(depthu, depthv)
#ifdef PVT_DEBUG
    write(*,*)'   searching for level-structures with depth [depth] : ', k
#endif

    ! We have now created the first 'v' point
    !   (v) -- we need to populate 'u' end points
    call rgn_init(uend, S%n)
    uend%n = 0
    ! also collect all uss
    allocate(uss(S%n))
    n_uss = 1
    do i = 1 , S%n
      call rgn_list(r, 1, S%r(i:i))

      call level_structure(n,nnzs,n_col,l_ptr,l_col,uss(n_uss),r,skip,priority=priority)
      j = level_structure_depth(uss(n_uss))
      if ( j == k ) then
        if ( .not. rgn_push(uend, S%r(i)) ) &
            call die('Should never happen')
#ifdef PVT_DEBUG
        write(*,*)'   populating uss [count / etr_small] : ', n_uss, S%r(i)
#endif
        n_uss = n_uss + 1
      end if

    end do
    call rgn_purge(uend)
    ! Sort u
    call rgn_sort(uend)
    n_uss = rgn_size(uend)

    ! Get all nodes with same degree as v
    !   (vi) -- we need to populate 'v' end points
    call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col, sub, pvt)
    call rgn_init(vend, sub%n)
    vend%n = 0
    ! Get degree of the root 'v' node
    degreev = n_col(vs%v%r(1))
    ! This is obviously too big, however the memory requirement of this
    ! type should be negligeble...
    allocate(vss(sub%n))
    n_vss = 1
    do i = 1, pvt%n
      ! Since the pivoting is sorted by increasing degree, we can easily
      if ( n_col(pvt%r(i)) > degreev ) exit

      call rgn_list(r, 1, pvt%r(i:i))

      ! Check degree of this node and store it if it has the same depth as v
      if ( n_col(pvt%r(i)) == degreev .and. .not. in_rgn(uend, pvt%r(i)) ) then
        call level_structure(n,nnzs,n_col,l_ptr,l_col,vss(n_vss),r,skip,priority=priority)
        j = level_structure_depth(vss(n_vss))
        if ( j == k ) then
          if ( .not. rgn_push(vend, pvt%r(i)) ) &
              call die('Should never happen')
#ifdef PVT_DEBUG
          write(*,*)'   populating vss [count / etr_small] : ', n_vss, pvt%r(i)
#endif
          n_vss = n_vss + 1
        end if
      end if

    end do
    call rgn_purge(vend)
    call rgn_sort(vend)
    n_vss = rgn_size(vend)

#ifdef PVT_DEBUG
    write(*,*)' v has peripherals ',n_vss
    write(*,*)' u has peripherals ',n_uss
#endif

    ! ### Completed algorithm I. ###

    ! At this point we have two sets, and for all the sets we have the sets

    ! ### Algorithm II ###

    ! Loop all elements and figure out their placement in the level-structure
    call init_level_structure(xs, k)
    call rgn_range(r, 1, n)
    call rgn_init(set, n)
    call rgn_init(pvt, n)
    set%n = 0
    pvt%n = 0

    ! (i) we have all Lvi = Lv
    !  -- (a) first assign all those with the same index
    step_a: do i = 1, sub%n

      ! Get currently searched element
      etr = sub%r(i)

      ! If we have already removed it, no need to check it
      if ( .not. in_rgn(r, etr) ) cycle

      ! get indices in 'v'
      vlvl = level_structure_element_lvl(vss(1), etr)
      do iuv = 2, n_vss
        if ( level_structure_element_lvl(vss(iuv), etr) /= vlvl ) then
          vlvl = 0
          cycle step_a
        end if
      end do
      ! quick skip if the level is not "uniform"
      if ( vlvl == 0 ) cycle

      ! get indices in 'u'
      j = level_structure_depth(uss(1)) + 1
      ulvl = j - level_structure_element_lvl(uss(1), etr)
      if ( ulvl /= vlvl ) cycle
      do iuv = 2, n_uss
        j = level_structure_depth(uss(iuv)) + 1

        if ( j - level_structure_element_lvl(uss(iuv), etr) /= ulvl ) then
          ulvl = 0
          cycle step_a
        end if
      end do

      ! All have the same level structure
      if ( .not. (rgn_push(set, etr) .and. rgn_push(lvl, ulvl)) ) then
        call die('GGPS: push to (II) (a)')
      end if

      ! Remove sets from r
      call rgn_consecutive_remove(r, etr)

      ! Now remove all connections
      j = l_ptr(etr)
      call rgn_list(S, n_col(etr), l_col(j+1:j+n_col(etr)))
      call rgn_consecutive_remove(r, S)

    end do step_a

    !  -- (b)
    ! split into different graphs
    if ( r%r(r%n) /= 0 ) then

      do i = 1 , sub%n

        if ( in_rgn(r, sub%r(i)) ) cycle

        ! Now figure out the full connectivity for this node

      end do

    end if

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'GGPS DONE'
#endif

  end subroutine GGPS_NEW

  subroutine rev_GGPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt,priority)
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    integer, intent(in), optional :: priority(n)

    call GGPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt, priority=priority)

    call rgn_reverse(pvt)

  end subroutine rev_GGPS

  ! The GPS algorithm 
  subroutine GPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt,priority)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! The level structure created by the algo_i_ii
    type(tLevelStructure) :: lvl
    type(tRgn) :: ipvt, r
    integer :: d, depth

#ifdef PVT_DEBUG
    write(*,fmt_debug1_ai) 'GPS START', n
#endif

    ! Ensure it is clean
    call rgn_delete(pvt)
    call rgn_init(pvt, sub%n)
    pvt%n = 0

    ! Find a set of pseudo-peripherals using the GPS algorithm
    call pseudo_peripheral_level_structure(D_LOW,n,nnzs,n_col,l_ptr,l_col, &
        sub, lvl, priority=priority)

    ! Depth of level-structure
    depth = level_structure_depth(lvl)
#ifdef PVT_DEBUG
    write(*,fmt_debug2_ai)'found left [peripheral]', lvl%v%r(1)
#endif

    ! Process the pivoting of each level
    do d = 1 , depth

      ! Order the level according to increasing degree
      call level_structure_level(lvl, d, ipvt)

#ifdef PVT_DEBUG
      write(*,fmt_debug2_aii) 'processing [level / size] ', d, ipvt%n
#endif

      call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,ipvt,r)
      if ( .not. rgn_push(pvt, r) ) call die('GPS push -- 1')

    end do

#ifdef PVT_DEBUG
    write(*,fmt_debug2_ai) 'found right [peripheral]', r%r(r%n)
#endif

    call rgn_delete(r,ipvt)
    call delete_level_structure(lvl)

    if ( pvt%n /= sub%n ) &
        call die('GPS: Error in algorithm')

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'GPS DONE'
#endif

  end subroutine GPS

  subroutine rev_GPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt,priority)
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    integer, intent(in), optional :: priority(n)

    call GPS(n,nnzs,n_col,l_ptr,l_col,sub,pvt, priority=priority)

    call rgn_reverse(pvt)

  end subroutine rev_GPS


  ! The connectivity graph
  ! This is a very simple pivoting method, entirely based on a connectivity
  ! algorithm.
  ! This may be started anywhere.
  subroutine connectivity_graph(n,nnzs,n_col,l_ptr,l_col,sub,pvt,&
      start, priority, only_sub)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The algorithm performs (only) if it knows where to start
    ! hence we can force the algorithm to start from some point
    type(tRgn), intent(in), optional :: start
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    ! Whether we should only allow to find the pivoting for the sub-graph
    ! This means that the pivoting table will not necessarily have all sub elements
    logical, intent(in), optional :: only_sub

    integer :: i, j, ptr, col
    type(tRgn) :: skip, con, r_tmp
    logical :: lonly_sub
    logical, allocatable :: in_pvt(:)

    lonly_sub = .false.
    if ( present(only_sub) ) lonly_sub = only_sub

    call rgn_init(pvt, sub%n)
    pvt%n = 0

    ! Start by initializing skip to be everything but the sub part
    call rgn_range(r_tmp, 1, n)
    call rgn_complement(sub, r_tmp, skip)
    call rgn_sort(skip)
    call rgn_delete(r_tmp)

    ! Initialize in_pvt
    allocate(in_pvt(n))
    in_pvt = .false.

    if ( present(start) ) then
      if ( .not. rgn_push(pvt, start) ) &
          call die('CG -- push 1')
      call rgn_sp_sort(pvt,n,nnzs,n_col,l_ptr,l_col, &
          start,R_SORT_MAX_FRONT)
      ! Initialize the connectivity graph (with the sorted elements)
      call rgn_copy(pvt, con)
    end if

    ! Initialize the logical array
    do i = 1, pvt%n
      in_pvt(pvt%r(i)) = .true.
    end do

    ! sort the starting configuration

    ! Continue the propagation
    do while ( pvt%n /= sub%n )

      ! Initialize room for the new columns
      if ( con%n > 0 ) then
        ! Get size
        j = 0
        do i = 1, con%n
          j = j + n_col(con%r(i))
        end do
        ! We keep the temporary region allocated.
        ! since we will use it every cycle.
        call rgn_grow(r_tmp, j)
        r_tmp%n = 0

        ! find all connections
        do i = 1, con%n
          ptr = l_ptr(con%r(i))
          do j = ptr + 1, ptr + n_col(con%r(i))
            col = l_col(j)
            if ( .not. (in_rgn(skip, col) .or. in_pvt(col)) ) then
              if ( .not. rgn_push(r_tmp, col) ) &
                  call die('CG -- push 2')
            end if
          end do
        end do

        ! Reduce to the unique values (this also sorts it)
        call rgn_uniq(r_tmp, in_place=.true.)

        ! Copy to the new connectivity region
        call copy(r_tmp, con)

      end if

      if ( con%n == 0 ) then

        ! Exit if done, or we if we shouldn't create the full graph
        if ( pvt%n == sub%n .or. lonly_sub ) exit

        ! Choose a random one
        do i = 1 , sub%n
          if ( in_rgn(skip, sub%r(i)) .or. in_pvt(sub%r(i)) ) cycle
          call rgn_grow(con, 1)
          con%n = 1
          con%r(1) = sub%r(i)
          exit
        end do

      end if

      ! Update pivoting check-indices
      do i = 1, con%n
        in_pvt(con%r(i)) = .true.
      end do

      if ( .not. rgn_push(pvt, con) ) call die('CG -- push 3')
      call rgn_sp_sort(pvt,n,nnzs,n_col,l_ptr,l_col, &
          con, R_SORT_MAX_BACK, r_logical=in_pvt)

    end do

    deallocate(in_pvt)
    call rgn_delete(skip, con, r_tmp)

    if ( pvt%n /= sub%n .and. .not. lonly_sub ) then
      call die('connectivity_graph: Error in algorithm')
    end if

  contains

    subroutine copy(from, to)
      type(tRgn), intent(in) :: from
      type(tRgn), intent(inout) :: to
      integer :: i
      call rgn_grow(to, from%n)
      do i = 1, from%n
        to%r(i) = from%r(i)
      end do
      to%n = from%n
    end subroutine copy

  end subroutine connectivity_graph

  subroutine rev_connectivity_graph(n,nnzs,n_col,l_ptr,l_col,sub,pvt,&
      start,priority, only_sub)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The algorithm performs (only) if it knows where to start
    ! hence we can force the algorithm to start from some point
    type(tRgn), intent(in), optional :: start
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    ! Whether we should only allow to find the pivoting for the sub-graph
    ! This means that the pivoting table will not necessarily have all sub elements
    logical, intent(in), optional :: only_sub

    call connectivity_graph(n,nnzs,n_col,l_ptr,l_col,sub,pvt,&
        start,priority, only_sub)

    call rgn_reverse(pvt)

  end subroutine rev_connectivity_graph


  ! The PCG algorithm 
  ! This algorithm uses the connectivity graph to sort each level.
  ! We use the peripheral search to find the end-points of the 
  ! peripheral graph. Then we add each level and sort
  ! This PCG algorithm is implemented and envisioned by Nick R. Papior.
  ! It proves more efficient than the GPS algorithm, and is a little slower.
  ! It however reduces the bandwidth by more.
  subroutine PCG(n,nnzs,n_col,l_ptr,l_col,sub,pvt,start,priority, only_sub)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The starting for the peripheral
    type(tRgn), intent(in), optional :: start
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    logical, intent(in), optional :: only_sub

    ! Another temporary level structure
    type(tLevelStructure), target :: level_s
    type(tLevelStructure) :: level_s2
    type(tLevelStructure), pointer :: ls

    ! The next starting peripheral
    type(tRgn) :: st, skip
    integer :: i, idx, st_n
    integer :: width, width2, depth, depth2
    logical :: lonly_sub

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'PCG START'
#endif

    lonly_sub = .false.
    if ( present(only_sub) ) lonly_sub = only_sub

    call rgn_range(st, 1, n)
    call rgn_complement(sub, st, skip)
    call rgn_delete(st)

    ! Ensure it is clean
    call rgn_delete(pvt)

    ! Algorithm I
    !   (i) -- pick arbitrary node with minimal degree
    if ( present(start) ) then
      call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s, start, skip, priority=priority)
    else
      idx = idx_degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col, sub, priority = priority)
      !   (ii/iii) -- search for deepest level structure in these peripherals
#ifdef PVT_DEBUG
      write(*,fmt_debug2_ai) 'initial V level structure [from]', sub%r(idx)
#endif
      call rgn_range(pvt, sub%r(idx), sub%r(idx))
      call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s, pvt, skip, priority=priority)

    end if

    search_smallest: do 

      !  (iv) -- calculate the width of this level structure
      depth = level_structure_depth(level_s)
      width = level_structure_width(level_s)

      ! extract the nodes at the last level to search for deeper ones
      call level_structure_level(level_s, depth, st)
#ifdef PVT_DEBUG
      write(*,fmt_debug2_aiii) 'checking for deeper/narrower from [depth / width / size]', depth, width, st%n
#endif

      !  (v) -- sort by the degree of the nodes (lowest to highest)
      call sort_degree(D_LOW_SUM,n,nnzs,n_col,l_ptr,l_col, st, pvt)
      call rgn_copy(pvt, st)
      call rgn_delete(pvt)

      !  (vi) -- for increasing included elements in the new
      !          level structure we try and find the smallest width structure
      !          for a subset
      st_n = st%n
      idx = 0
      do i = 1 , st_n
        
        st%n = i
        call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s2, st, skip, priority=priority)

        depth2 = level_structure_depth(level_s2)
        width2 = level_structure_width(level_s2)

#ifdef PVT_DEBUG
        write(*,fmt_debug3_aiii) 'level structure [depth / width / size]', depth2, width2, i
#endif

        ! If the level depth is larger than the current start over 
        ! from that index (the last added element *must* be
        ! the largest level)
        if ( width2 < width .or. depth2 > depth ) then
#ifdef PVT_DEBUG
          write(*,fmt_debug3_a, advance='NO') ' FOUND '
          if ( width2 < width ) write(*,'(a)', advance='NO') 'narrower '
          if ( depth2 > depth ) write(*,'(a)', advance='NO') 'deeper '
          write(*,'(a,tr1,i8)') 'level structure [size]', i
#endif
          idx = i

          ! (vii) -- switch to this new peripheral structure
          call copy_level_structure(level_s2, level_s)
          depth = depth2
          width = width2

        end if

        call delete_level_structure(level_s2)

      end do

      if ( idx > 0 ) cycle

      exit

    end do search_smallest

    call rgn_delete(st)

    ! Copy over the new stuff
    call rgn_init(pvt, sub%n)
    pvt%n = 0
    ls => level_s
    do
      if ( .not. rgn_push(pvt, ls%v) ) call die('Failure in PCG -- 1')

      if ( .not. associated(ls%next) ) exit

      ls => ls%next
    end do

    call delete_level_structure(level_s)

    if ( pvt%n /= sub%n .and. .not. lonly_sub ) &
        call die('PCG: Error in algorithm')

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'PCG DONE'
#endif

  end subroutine PCG

  subroutine rev_PCG(n,nnzs,n_col,l_ptr,l_col,sub,pvt,start,priority,only_sub)
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    type(tRgn), intent(in), optional :: start
    integer, intent(in), optional :: priority(n)
    logical, intent(in), optional :: only_sub

    call PCG(n,nnzs,n_col,l_ptr,l_col,sub,pvt,start,priority,only_sub)

    call rgn_reverse(pvt)

  end subroutine rev_PCG


#ifdef SIESTA__METIS

  subroutine metis_adjacency_graph(n,nnzs,n_col,l_ptr,l_col,sub, &
      xadj, adjncy, w, priority)
    use iso_c_binding, only: c_int
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    integer(c_int), intent(inout), allocatable :: xadj(:), adjncy(:), w(:)
    integer, intent(in), optional :: priority(n)

    integer :: io, i, ptr, ind, nc, j, nadj

    ! Allocate adjacency graphs
    allocate(xadj(0:sub%n), w(sub%n) )

    ! First count adjacencies
    xadj(0) = 0 ! + 1 : fortran-style
    do i = 1 , sub%n

      io = sub%r(i)
      ptr = l_ptr(io)
      nc = n_col(io)

      ! Count number of elements in 
      ! the sub-space
      nadj = 0
      do ind = ptr + 1 , ptr + nc
        if ( in_rgn(sub,l_col(ind)) ) then
          ! Skip "on-site" connections
          if ( l_col(ind) /= io ) then
            nadj = nadj + 1
          end if
        end if
      end do

      xadj(i) = xadj(i-1) + nadj

      if ( present(priority) ) then
        w(i) = priority(io) + 1
      else
        w(i) = 1
      end if

    end do

    ! transfer to local adjacency graph
    allocate( adjncy(xadj(sub%n)) )

    ! Create adjncy 
    nadj = 0
    do i = 1 , sub%n

      io = sub%r(i)
      ptr = l_ptr(io)
      nc = n_col(io)

      ! Count number of elements in 
      ! the sub-space
      do ind = ptr + 1 , ptr + nc
        j = rgn_pivot(sub,l_col(ind))
        if ( j > 0 ) then
          ! Skip "on-site" connections
          if ( l_col(ind) /= io ) then
            nadj = nadj + 1
            adjncy(nadj) = j - 1 ! + 1 : fortran style
          end if
        end if
      end do

      if ( nadj /= xadj(i) ) then
        print *,i,nadj, xadj(i)
        call die('metis_adjacency_graph: Error in creating &
            &adjacency graph.')
      end if

    end do

  end subroutine metis_adjacency_graph


  subroutine metis_NodeND_pvt(n,nnzs,n_col,l_ptr,l_col,sub,pvt, priority)
    use iso_c_binding, only: c_int, c_ptr, c_loc
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    integer, intent(in), optional :: priority(n)

    ! METIS variables
    integer(c_int), allocatable :: xadj(:), adjncy(:)
    integer(c_int), allocatable :: perm(:), iperm(:)
    integer(c_int), allocatable, target :: w(:)
    integer(c_int) :: iret, nvtxs, opts(100) ! In 5.0.1 it is 40, but ...
    type(c_ptr) :: wp

    interface
      integer(c_int) function METIS_SetDefaultOptions(opts) &
          bind(C, name="METIS_SetDefaultOptions")
        use iso_c_binding, only: c_int
        implicit none
        integer(c_int), dimension(*) :: opts
      end function METIS_SetDefaultOptions
      integer(c_int) function METIS_NodeND(nvtxs,xadj,adjncy,vwgt, &
          opts,perm,iperm) bind(C, name="METIS_NodeND")
        use iso_c_binding, only: c_int, c_ptr
        implicit none
        integer(c_int) :: nvtxs
        integer(c_int), dimension(*) :: xadj, adjncy, perm, iperm
        type(c_ptr), value :: vwgt
        integer(c_int), dimension(*) :: opts
      end function METIS_NodeND
    end interface

    integer :: i

    call rgn_delete(pvt)

    ! The following does C-style indexing 
    ! as the internal METIS structure is a simple offset

    !   call METIS_setdefaultoptions(opts)
    iret = METIS_setdefaultoptions(opts)
    ! METIS_OK == 1
    if ( iret /= 1 ) then
      opts(6) = 256 ! increase debug level and re-run for full dbg-lvl
      iret = METIS_setdefaultoptions(opts)
      call die('metis_NodeND: Error on initializing default options.')
    end if

    ! set options
    opts(3)  =  1 ! CTYPE == Sorted heavy-edge matching
    !    opts(5)  =  1 ! RTYPE == Greedy-based cut and volume refinement
    opts(7)  = 20 ! NITER(10) == Number of iterations
    opts(8)  =  1 ! NCUT == Number of cuts
    opts(10) =  1 ! NO2HOP == does not use 2 hop
    opts(11) =  1 ! MINCONN == Explicitly minimize the maximum connectivity
    opts(12) =  0 ! CONTIG == Forces contiguous 
    opts(13) =  1 ! COMPRESS == compress similar adjacency nodes
    opts(16) =  1 ! NSEPS(1) == tries in the separator

    ! Setup the adjacency graph
    call metis_adjacency_graph(n,nnzs,n_col,l_ptr,l_col,sub, &
        xadj, adjncy, w, priority)

    allocate(perm(sub%n), iperm(sub%n))

    ! Call metis
    wp = c_loc(w(1))
    nvtxs = sub%n
    !   call METIS_NodeND(nvtxs, xadj, adjncy, wp, opts, perm, iperm)
    iret = METIS_NodeND(nvtxs, xadj, adjncy, wp, opts, perm, iperm)
    if ( iret /= 1 ) then ! METIS_OK == 1
      print *,iret
      opts(6) = 256
      iret = METIS_NodeND(nvtxs, xadj, adjncy, wp, opts, perm, iperm)
      print *,iret
      call die('metis_NodeND: error in pivoting.')
    end if

    ! Clean-up
    deallocate(xadj,adjncy,w,iperm)

    call rgn_init(pvt,sub%n)

    ! Transfer pivoting to actual pivoting index
    do i = 1 , sub%n
      pvt%r(i) = sub%r(perm(i)+1) ! - 1 : fortran style
    end do

    ! Clean-up
    deallocate(perm)

    if ( pvt%n /= sub%n ) call die('metis_NodeND: Error in algorithm')

  end subroutine metis_NodeND_pvt

  subroutine metis_PartGraphKway_pvt(n,nnzs,n_col,l_ptr,l_col,sub,pvt, priority)
    use iso_c_binding, only: c_int, c_ptr, c_loc, c_null_ptr
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    integer, intent(in), optional :: priority(n)

    ! METIS variables
    integer(c_int), allocatable :: xadj(:), adjncy(:)
    integer(c_int), allocatable :: part(:)
    integer(c_int), allocatable, target :: w(:)
    integer(c_int) :: iret, nvtxs, ncon, nparts, old_objval, objval, opts(100) ! In 5.0.1 it is 40, but ...
    type(c_ptr) :: wp

    interface
      integer(c_int) function METIS_SetDefaultOptions(opts) &
          bind(C, name="METIS_SetDefaultOptions")
        use iso_c_binding, only: c_int
        implicit none
        integer(c_int), dimension(*) :: opts
      end function METIS_SetDefaultOptions
      integer(c_int) function METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt, &
          vsize,adjwgt,nparts,tpwgts,ubvec,opts,objval,part) bind(C, name="METIS_PartGraphKway")
        use iso_c_binding, only: c_int, c_ptr
        implicit none
        integer(c_int) :: nvtxs, ncon, nparts, objval
        integer(c_int), dimension(*) :: xadj, adjncy, opts, part
        type(c_ptr), value :: vwgt, vsize, adjwgt, tpwgts, ubvec
      end function METIS_PartGraphKway
    end interface

    ! variables for the loop
    integer :: i, ip, j
    integer :: old_bw, bw
    type(tRgn) :: next_pvt

    call rgn_delete(pvt)

    ! The following does C-style indexing 
    ! as the internal METIS structure is a simple offset

    !   call METIS_setdefaultoptions(opts)
    iret = METIS_setdefaultoptions(opts)
    ! METIS_OK == 1
    if ( iret /= 1 ) then
      opts(6) = 256 ! increase debug level and re-run for full dbg-lvl
      iret = METIS_setdefaultoptions(opts)
      call die('metis_PartGraphKway: Error on initializing default options.')
    end if

    ! set options
    opts(2)  =  0 ! OBJTYPE == edge-cut minimization
    opts(3)  =  1 ! CTYPE == Sorted heavy-edge matching
    !    opts(5)  =  1 ! RTYPE == Greedy-based cut and volume refinement
    opts(7)  = 20 ! NITER(10) == Number of iterations
    opts(8)  =  1 ! NCUTS == Number of cuts
    opts(10) =  1 ! NO2HOP == does not use 2 hop
    opts(11) =  1 ! MINCONN == Explicitly minimize the maximum connectivity
    opts(12) =  0 ! CONTIG == Forces contiguous 

    ! Allocate adjacency graphs
    call metis_adjacency_graph(n,nnzs,n_col,l_ptr,l_col,sub, &
        xadj, adjncy, w, priority)

    allocate(part(sub%n))

    ! Initialize the pivoting array
    call rgn_init(pvt,sub%n)
    call rgn_init(next_pvt,sub%n)

    iret = 1
    nvtxs = sub%n
    ncon = 1
    nparts = 1
    old_bw = huge(1)
    wp = c_loc(w(1))
    do nparts = 2 , min(nvtxs / 2 + 1, nvtxs)

      ! Call metis
      iret = METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, wp, &
          C_NULL_PTR, C_NULL_PTR, &
          nparts, & ! number of parts
          C_NULL_PTR, C_NULL_PTR, opts, objval, part)

      ! An error forces us to exit loop
      if ( iret /= 1 ) exit

      j = 0
      do ip = 0, nparts - 1
        do i = 1 , sub%n
          if ( part(i) == ip ) then
            j = j + 1
            next_pvt%r(j) = sub%r(i)
          end if
        end do
      end do

      ! Transfer pivoting to actual pivoting index
      bw = bandwidth(n,nnzs,n_col,l_ptr,l_col,next_pvt)
      if ( bw == old_bw ) then
        if ( profile(n,nnzs,n_col,l_ptr,l_col,next_pvt) &
            < profile(n,nnzs,n_col,l_ptr,l_col,pvt) ) then
          pvt%r(:) = next_pvt%r(:)
        end if
      else if ( bw < old_bw ) then
        pvt%r(:) = next_pvt%r(:)
      end if

    end do

    if ( iret /= 1 ) then ! METIS_OK == 1
      print *,iret
      call die('metis_PartGraphKway: error in pivoting.')
    end if

    ! Clean-up
    deallocate(xadj,adjncy,w,part)
    call rgn_delete(next_pvt)

    if ( pvt%n /= sub%n ) call die('metis_PartGraphKway: Error in algorithm')

  end subroutine metis_PartGraphKway_pvt

  subroutine metis_PartGraphRecursive_pvt(n,nnzs,n_col,l_ptr,l_col,sub,pvt, priority)
    use iso_c_binding, only: c_int, c_ptr, c_loc, c_null_ptr
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    integer, intent(in), optional :: priority(n)

    ! METIS variables
    integer(c_int), allocatable :: xadj(:), adjncy(:)
    integer(c_int), allocatable :: part(:)
    integer(c_int), allocatable, target :: w(:)
    integer(c_int) :: iret, nvtxs, ncon, nparts, old_objval, objval, opts(100) ! In 5.0.1 it is 40, but ...
    type(c_ptr) :: wp

    interface
      integer(c_int) function METIS_SetDefaultOptions(opts) &
          bind(C, name="METIS_SetDefaultOptions")
        use iso_c_binding, only: c_int
        implicit none
        integer(c_int), dimension(*) :: opts
      end function METIS_SetDefaultOptions
      integer(c_int) function METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt, &
          vsize,adjwgt,nparts,tpwgts,ubvec,opts,objval,part) bind(C, name="METIS_PartGraphRecursive")
        use iso_c_binding, only: c_int, c_ptr
        implicit none
        integer(c_int) :: nvtxs, ncon, nparts, objval
        integer(c_int), dimension(*) :: xadj, adjncy, opts, part
        type(c_ptr), value :: vwgt, vsize, adjwgt, tpwgts, ubvec
      end function METIS_PartGraphRecursive
    end interface

    ! variables for the loop
    integer :: i, ip, j
    integer :: old_bw, bw
    type(tRgn) :: next_pvt

    call rgn_delete(pvt)

    ! The following does C-style indexing 
    ! as the internal METIS structure is a simple offset

    !   call METIS_setdefaultoptions(opts)
    iret = METIS_setdefaultoptions(opts)
    ! METIS_OK == 1
    if ( iret /= 1 ) then
      opts(6) = 256 ! increase debug level and re-run for full dbg-lvl
      iret = METIS_setdefaultoptions(opts)
      call die('metis_PartGraphRecursive: Error on initializing default options.')
    end if

    ! set options
    opts(2)  =  0 ! OBJTYPE == edge-cut minimization
    opts(3)  =  1 ! CTYPE == Sorted heavy-edge matching
    !    opts(5)  =  1 ! RTYPE == Greedy-based cut and volume refinement
    opts(7)  = 20 ! NITER(10) == Number of iterations
    opts(8)  =  1 ! NCUTS == Number of cuts
    opts(10) =  1 ! NO2HOP == does not use 2 hop
    opts(11) =  1 ! MINCONN == Explicitly minimize the maximum connectivity
    opts(12) =  0 ! CONTIG == Forces contiguous 

    ! Allocate adjacency graphs
    call metis_adjacency_graph(n,nnzs,n_col,l_ptr,l_col,sub, &
        xadj, adjncy, w, priority)

    allocate(part(sub%n))

    ! Initialize the pivoting array
    call rgn_init(pvt,sub%n)
    call rgn_init(next_pvt,sub%n)

    iret = 1
    nvtxs = sub%n
    ncon = 1
    nparts = 1
    old_bw = huge(1)
    wp = c_loc(w(1))
    do nparts = 2 , min(nvtxs / 2 + 1, nvtxs)

      ! Call metis
      iret = METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, wp, &
          C_NULL_PTR, C_NULL_PTR, &
          nparts, & ! number of parts
          C_NULL_PTR, C_NULL_PTR, opts, objval, part)

      ! An error forces us to exit loop
      if ( iret /= 1 ) exit

      j = 0
      do ip = 0, nparts - 1
        do i = 1 , sub%n
          if ( part(i) == ip ) then
            j = j + 1
            next_pvt%r(j) = sub%r(i)
          end if
        end do
      end do

      ! Transfer pivoting to actual pivoting index
      bw = bandwidth(n,nnzs,n_col,l_ptr,l_col,next_pvt)
      if ( bw == old_bw ) then
        if ( profile(n,nnzs,n_col,l_ptr,l_col,next_pvt) &
            < profile(n,nnzs,n_col,l_ptr,l_col,pvt) ) then
          pvt%r(:) = next_pvt%r(:)
        end if
      else if ( bw < old_bw ) then
        pvt%r(:) = next_pvt%r(:)
      end if

    end do

    if ( iret /= 1 ) then ! METIS_OK == 1
      print *,iret
      call die('metis_PartGraphRecursive: error in pivoting.')
    end if

    ! Clean-up
    deallocate(xadj,adjncy,w,part)
    call rgn_delete(next_pvt)

    if ( pvt%n /= sub%n ) call die('metis_PartGraphRecursive: Error in algorithm')

  end subroutine metis_PartGraphRecursive_pvt
#endif


  subroutine breadth_first_search(n, nnzs, n_col, l_ptr, l_col, start, bfs, skip)
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)

    ! BFS search start
    type(tRgn), intent(in) :: start
    ! The BFS table
    type(tRgn), intent(inout) :: bfs
    ! Optional discarded elements
    type(tRgn), intent(in), optional :: skip

    ! The current queue of elements
    type(tRgn) :: queue, sskip

    integer :: nel
    integer :: i, ptr, el, eln

#ifdef PVT_DEBUG
    write(*,fmt_debug2_a) 'BFS START'
#endif
    ! NOT FUNCTIONAL YET

    if ( present(skip) ) then
      call rgn_init_consecutive(n, sskip, skip)
    else
      call rgn_init_consecutive(n, sskip)
    end if

    ! Actual number of elements in the level
    nel = n - rgn_size(sskip)

    ! Pre-allocate the full queue
    call rgn_init(queue, nel)
    queue%n = 0

    ! Start queue by adding all start elements
    if ( .not. rgn_push(queue, start) ) call die('Error in BFS -- 1')
    call rgn_consecutive_insert(sskip, start)

    ! Initialize the resulting breadth-first-search table
    call rgn_init(bfs, nel)
    bfs%n = 0

    ! Start creating the BFS
    do while ( rgn_size(bfs) < nel )

      if ( rgn_size(queue) == 0 ) then
        ! We have to "manually" add a fictitious point
        ! that hasn't been added yet, i.e. if the queue is
        ! empty. We have to keep filling it

        do el = 1 , n

          ! Ensure it is not skipped
          if ( in_rgn(sskip, el) ) cycle

          if ( rgn_size(queue) < rgn_size(bfs) ) then
            ! Ensure it is not in the queue
            if ( in_rgn(queue, el) ) cycle
            ! Ensure it is not already assigned
            if ( in_rgn(bfs, el) ) cycle
          else
            if ( in_rgn(bfs, el) ) cycle
            if ( in_rgn(queue, el) ) cycle
          end if

#ifdef PVT_DEBUG
          write(*,fmt_debug2_ai) 'adding to queue (non-connected) [element]', el
#endif
          if ( .not. rgn_push(queue, el) ) call die('Error in BFS -- 2')
          exit

        end do

      end if

      ! Push queued value
      el = rgn_pop(queue)
      if ( el == 0 ) call die('Error in BFS -- 3')
      if ( .not. rgn_push(bfs, el) ) call die('Error in BFS -- 4')

      ! Examine all neighbours
      ptr = l_ptr(el)
      do i = 1, n_col(el)

        eln = l_col(ptr+i)

        ! Ensure it is not skipped
        if ( in_rgn(sskip, eln) ) cycle

        if ( rgn_size(queue) < rgn_size(bfs) ) then
          ! Ensure it is not in the queue
          if ( in_rgn(queue, eln) ) cycle
          ! Ensure it is not already assigned
          if ( in_rgn(bfs, eln) ) cycle
        else
          if ( in_rgn(bfs, eln) ) cycle
          if ( in_rgn(queue, eln) ) cycle
        end if

        ! Now we can add it to the queue
        if ( .not. rgn_push(queue, eln) ) call die('Error in BFS -- 5')

      end do

    end do

    call rgn_delete(queue, sskip)

  end subroutine breadth_first_search

  subroutine level_structure(n, nnzs, n_col, l_ptr, l_col, ls, start, skip, priority, only_sub)
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)

    ! The LV table
    type(tLevelStructure), intent(inout), target :: ls
    ! Level structure search start
    type(tRgn), intent(in), optional :: start
    ! Optional discarded elements
    type(tRgn), intent(in), optional :: skip
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    ! Allow only a sub-graph
    logical, intent(in), optional :: only_sub

    ! The current queue of elements
    type(tRgn) :: queue
    logical, allocatable :: log_skip(:)

    type(tLevelStructure), pointer :: clvl
    integer :: i, nel, el, nadded
    logical :: lonly_sub

    lonly_sub = .false.
    if ( present(only_sub) ) lonly_sub = only_sub

    call delete_level_structure(ls)

    allocate(log_skip(n))
    log_skip(:) = .false.

    nel = 0
    if ( present(skip) ) then
      nel = skip%n

      ! Speed up searches in the skipped elements
      do i = 1, skip%n
        log_skip(skip%r(i)) = .true.
      end do

    end if

    ! Actual number of elements in the total level structure
    nel = n - nel

    if ( nel <= 1 ) then

      call rgn_init(ls%v, nel)
      ls%v%n = 0
      do i = 1, n
        if ( .not. log_skip(i) ) then
          if ( .not. rgn_push(ls%v, i) ) call die('level_structure: push -- 1')
        end if
      end do
      ls%lvl = 1

      deallocate(log_skip)

      return

    end if

    ! Ensure the queue has space for all elements
    call rgn_init(queue, nel)
    queue%n = 0

    ! Start queue by adding all start elements
    if ( present(start) ) then
      if ( .not. rgn_push(queue, start) ) call die('Error in LS -- 1')
    end if

    ! Ensure we have all elements in one array, to easy skip
    do i = 1, queue%n
      log_skip(queue%r(i)) = .true.
    end do

    ! Counter to count number of elements added
    nadded = 0

    ! Initialize the level structure
    ls%lvl = 1
    clvl => ls
    do

      if ( rgn_size(queue) == 0 ) then

        ! Ensure it has space
        call rgn_grow(queue, 1)
        queue%n = 0

        ! Quick escape if the graph is finished
        if ( lonly_sub ) exit

        ! We have to "manually" add a fictitious point
        ! that hasn't been added yet, i.e. if the queue is
        ! empty. We have to keep filling it...

        ! ( we could optionally do a sort_degree call
        !   to select the one with the lowest degree)

        do el = 1 , n
          ! Ensure it is not skipped
          if ( log_skip(el) ) cycle
          log_skip(el) = .true.

          if ( .not. rgn_push(queue, el) ) call die('Error in LS -- 2')
          exit

        end do

      end if

      ! Copy to the current level
      call rgn_copy(queue, clvl%v)

      ! Update number of added elements
      nadded = nadded + rgn_size(clvl%v)

      if ( nadded >= nel ) exit

      ! Prepare queue
      call rgn_grow(queue, nel - nadded)

      ! Get the next level
      call next_level(n,nnzs,n_col,l_ptr,l_col,clvl%v, queue, priority)

      ! Allocate next level and prepare
      allocate(clvl%next)
      clvl%next%lvl = clvl%lvl + 1

      clvl => clvl%next

    end do

    call rgn_delete(queue)
    deallocate(log_skip)

  contains

    subroutine next_level(n,nnzs,n_col,l_ptr,l_col, qin, qout, priority)
      integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
      integer, intent(in), optional :: priority(n)

      type(tRgn), intent(in) :: qin
      type(tRgn), intent(inout) :: qout
      integer :: i, ind, el, eln

      ! Prepare output
      qout%n = 0

      ! Start adding the queued elements
      do i = 1, rgn_size(qin)

        ! Get current row
        el = qin%r(i)
        if ( el == 0 ) call die('Error in LS -- 3')

        ! Examine all neighbours
        do ind = l_ptr(el) + 1, l_ptr(el) + n_col(el)
          eln = l_col(ind)

          ! Ensure it is not already assigned
          if ( log_skip(eln) ) cycle
          log_skip(eln) = .true.

          ! Now we can add it to the queue
          if ( .not. rgn_push(qout, eln) ) call die('Error in LS -- 5')

        end do

      end do

      call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,qout,qout,priority)

    end subroutine next_level

  end subroutine level_structure


  ! Returns the depth of the current level structure
  function level_structure_depth(level_s) result(depth)
    type(tLevelStructure), intent(in), target :: level_s
    integer :: depth
    type(tLevelStructure), pointer :: ls

    depth = 1
    ls => level_s

    do while ( associated(ls%next) )
      depth = depth + 1
      ls => ls%next
    end do

  end function level_structure_depth

  ! Returns the width of the current level structure
  function level_structure_width(level_s, ilvl) result(width)
    type(tLevelStructure), intent(in), target :: level_s
    integer, intent(in), optional :: ilvl

    type(tLevelStructure), pointer :: ls
    integer :: i, width

    ls => level_s

    if ( present(ilvl) ) then
      do i = 1, ilvl - 1
        ls => ls%next
      end do
      width = rgn_size(ls%v)
    else
      width = 0
      do
        width = max(width, rgn_size(ls%v))
        if ( .not. associated(ls%next) ) exit
        ls => ls%next
      end do
    end if

  end function level_structure_width

  ! Return the level index in the level structure in which el resides
  function level_structure_element_lvl(level_s, el) result(ilvl)
    type(tLevelStructure), intent(in), target :: level_s
    integer, intent(in) :: el
    integer :: ilvl

    type(tLevelStructure), pointer :: ls

    ! Simply go to the level and copy
    ls => level_s
    ilvl = 1
    do
      if ( in_rgn(ls%v, el) ) return
      if ( .not. associated(ls%next) ) then
        ilvl = 0
        return
      end if
      ilvl = ilvl + 1
      ls => ls%next
    end do

  end function level_structure_element_lvl


  ! Extract only a certain level from the level structure
  subroutine level_structure_level(level_s,ilvl,r)
    type(tLevelStructure), intent(in), target :: level_s
    integer, intent(in) :: ilvl
    type(tRgn), intent(inout) :: r

    integer :: i
    type(tLevelStructure), pointer :: ls

    ! Simply go to the level and copy
    ls => level_s
    do i = 1, ilvl - 1
      ls => ls%next
    end do
    call rgn_copy(ls%v, r)

  end subroutine level_structure_level


  subroutine pseudo_peripheral_level_structure(method,n,nnzs,n_col,l_ptr,l_col, &
      sub, level_s, start, &
      small_wd, priority)
    ! The method for the degree search
    integer, intent(in) :: method
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The level tree
    type(tLevelStructure), intent(inout) :: level_s
    ! An initial node(s)
    type(tRgn), intent(in), optional :: start
    ! The user can optionally obtain the smallest width and degree of
    ! the sets created by the S set in the found level set
    ! It also returns the entry that creates that structure
    integer, intent(out), optional :: small_wd(3)
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! Another temporary level structure
    type(tLevelStructure) :: level_s2

    ! Local variables
    type(tRgn) :: pvt, S, skip
    integer :: idx, i, j
    integer :: depth

#ifdef PVT_DEBUG
    write(*,fmt_debug2_a) 'PPLS START'
#endif

    call rgn_range(skip, 1, n)
    call rgn_complement(sub, skip, skip)

    ! Algorithm I
    !   (i) -- pick arbitrary node with minimal degree
    if ( present(start) ) then
      call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s, start, skip, &
          priority=priority)
    else
      idx = idx_degree(method,n,nnzs,n_col,l_ptr,l_col, sub, priority = priority)
      !   (ii/iii) -- search for deepest level structure in these peripherals
#ifdef PVT_DEBUG
      write(*,fmt_debug3_ai) 'initial V level structure [from]', sub%r(idx)
#endif
      call rgn_range(pvt, sub%r(idx), sub%r(idx))
      call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s, pvt, skip, &
          priority=priority)
    end if

    search_deepest: do 

      depth = level_structure_depth(level_s)
      
      if ( present(small_wd) ) then
        !  (iv) -- collect the depth of the smallest width structure
        small_wd(1) = huge(1)
        small_wd(2) = 0
        small_wd(3) = 0
      end if

      ! extract the nodes at the last level to search for deeper ones
      call level_structure_level(level_s, depth, S)
#ifdef PVT_DEBUG
      write(*,fmt_debug3_aii) 'checking for deeper from [depth / size]', depth, S%n
#endif

      ! Sort by the degree of the nodes (lowest to highest)
      call sort_degree(method,n,nnzs,n_col,l_ptr,l_col, S, pvt)
      call rgn_copy(pvt, S)

      do i = 1 , S%n
        call rgn_range(pvt, S%r(i), S%r(i))
        call level_structure(n, nnzs, n_col, l_ptr, l_col, level_s2, pvt, skip, &
            priority=priority)

        if ( present(small_wd) ) then
          j = level_structure_width(level_s2)
          if ( j < small_wd(1) ) then
            small_wd(1) = j
            small_wd(2) = level_structure_depth(level_s2)
            small_wd(3) = S%r(i)
          end if
        end if
        j = level_structure_depth(level_s2)
#ifdef PVT_DEBUG
        write(*,fmt_debug3_aii) 'level structure [from / depth]', S%r(i), j
#endif

        ! If the level depth is larger than the current start over 
        ! from that index (the last added element *must* be
        ! the largest level)
        if ( depth < j ) then
#ifdef PVT_DEBUG
          write(*,fmt_debug3_ai) '  deeper V set [from]', pvt%r(1)
#endif
          call copy_level_structure(level_s2, level_s)
          cycle search_deepest

        end if

      end do

      exit search_deepest

    end do search_deepest

    call delete_level_structure(level_s2)

    call rgn_delete(S, pvt, skip)

#ifdef PVT_DEBUG
    write(*,fmt_debug2_ai) 'found depth', depth
#endif

#ifdef PVT_DEBUG
    write(*,fmt_debug2_a) 'PPLS DONE'
#endif

  end subroutine pseudo_peripheral_level_structure

  recursive subroutine init_level_structure(ls, lvl)
    type(tLevelStructure), intent(inout) :: ls
    integer, intent(in) :: lvl

    call delete_level_structure(ls)
    if ( lvl <= 0 ) return
    allocate(ls%next)

    call init_level_structure(ls%next, lvl - 1)

  end subroutine init_level_structure

  recursive subroutine delete_level_structure(ls)
    type(tLevelStructure), intent(inout) :: ls

    if ( associated(ls%next) ) then
      call delete_level_structure(ls%next)
      deallocate(ls%next)
    end if

    call rgn_delete(ls%v)

  end subroutine delete_level_structure

  recursive subroutine copy_level_structure(from, to)
    type(tLevelStructure), intent(inout) :: from, to

    call delete_level_structure(to)

    ! Now copy
    to%lvl = from%lvl
    call rgn_copy(from%v, to%v)

    if ( associated(from%next) ) then
      allocate(to%next)
      call copy_level_structure(from%next, to%next)
    end if

  end subroutine copy_level_structure


  function bandwidth(n,nnzs,n_col,l_ptr,l_col,sub) result(beta)
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn) :: pvt
    integer :: beta
    integer :: i, j, ind, idx

    call rgn_init(pvt, n, val=0)
    do i = 1, sub%n
      pvt%r(sub%r(i)) = i
    end do

    beta = 0
    do i = 1 , sub%n
      idx = sub%r(i)
      do ind = l_ptr(idx) + 1 , l_ptr(idx) + n_col(idx)
        ! figure out the pivoting place
        j = pvt%r(l_col(ind))
        if ( j <= 0 ) cycle
        beta = max(beta,i-j)
      end do
    end do

    call rgn_delete(pvt)

  end function bandwidth

  function profile(n,nnzs,n_col,l_ptr,l_col,sub) result(p)
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    integer(i8b) :: p
    type(tRgn) :: pvt
    integer :: beta
    integer :: i, j, ind, idx

    call rgn_init(pvt, n, val=0)
    do i = 1, sub%n
      pvt%r(sub%r(i)) = i
    end do

    p = 0
    do i = 1 , sub%n
      idx = sub%r(i)
      beta = 0
      do ind = l_ptr(idx) + 1 , l_ptr(idx) + n_col(idx)
        ! figure out the pivoting place
        j = pvt%r(l_col(ind))
        if ( j <= 0 ) cycle
        beta = max(beta,i-j)
      end do
      p = p + beta
    end do

    call rgn_delete(pvt)

  end function profile

  ! Returns the width of the current level structure
  function lvl_width(lvl,ilvl) result(width)
    type(tRgn), intent(in) :: lvl
    integer, intent(in), optional :: ilvl
    integer :: width
    integer :: i, depth
    if ( present(ilvl) ) then
      width = count(lvl%r(1:lvl%n) == ilvl)
      return
    end if
    depth = lvl_depth(lvl)
    width = 0
    do i = 1 , depth
      width = max(width,count(lvl%r(1:lvl%n) == i))
    end do
  end function lvl_width

  ! Returns the depth of the current level structure
  function lvl_depth(lvl) result(depth)
    type(tRgn), intent(in) :: lvl
    integer :: depth
    depth = maxval(lvl%r(1:lvl%n))
  end function lvl_depth

  ! This routine creates a level structure
  subroutine level_struct(n,nnzs,n_col,l_ptr,l_col,idx,sub,pvt,lvl,priority)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! We start this pseudo-periphael from this index
    integer, intent(in) :: idx
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays,
    ! the level associated with each added index (lvl%r(1) is the level
    ! of element pvt%r(1))
    type(tRgn) :: pvt, lvl
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! Temporary region used to contain the connectivity graph
    type(tRgn) :: con, con_c, skip

    ! local variables
    integer :: i, etr, iLvl, ncon
    logical :: suc

    ! create the initial consecutive region
    call rgn_range(con, 1, n)
    call rgn_complement(sub, con, con)
    ! now create the initial consecutive region
    call rgn_init_consecutive(n, skip, con)

    ! initialize the pivoting array
    call rgn_init(pvt,sub%n)
    pvt%n = 0
    call rgn_init(lvl,sub%n)
    lvl%n = 0

    ! Counter for level
    iLvl = 0

    ! initialize loop connectivity graph
    call rgn_init(con_c,sub%n)
    con_c%n = 0

    ! initialize connectivity region
    call rgn_init(con, sub%n)
    con%n = 1
    con%r(1) = sub%r(idx)

    do while ( pvt%n < sub%n )

      ! 1. Add the current connectivity graph
      !    to the pivot table and increment level
      iLvl = iLvl + 1

      !  If the connectivity happens to be zero,
      !  we need to add an arbitrary element with lowest degree
      !  In TS this will probably be one of the worst choices, yet
      !  it is hard to select another node on another basis.
      if ( con%n == 0 ) then
        ! this limits con_c to those not chosen
        call rgn_complement(skip, sub, con_c)
        if ( con_c%n > 0 ) then
          i = idx_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,con_c, priority = priority )
          ! We will limit the addition to 1 element
          if ( .not. rgn_push(con, con_c%r(i)) ) call die('level_struct: push -- 1')
        end if
      end if

      ! Insert the current connectivity graph
      call rgn_consecutive_insert(skip, con)
      ncon = con%n
      do while ( con%n > 0 )
        ! Get index with lowest degree
        i = idx_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,con, priority = priority )
        etr = rgn_pop(con, i) ! remove entry in con
        suc = rgn_push(pvt, etr)
        suc = rgn_push(lvl, iLvl)
      end do

      ! A simple check to see if we have finished the
      ! level structure
      if ( pvt%n == sub%n ) exit

      ! 2. Create connectivity graph of all entries just added
      !    This will loop on all previously connected entries.
      !    We create a new connectivity graph and let it be
      !    added on the following loop
      ! 2a. Note that on entry con%n == 0 (due to rgn_pop)
      do i = 1 , ncon
        ! Find connections from followed entry
        ! we pick-up the entry in the order they were added
        ! to 'pvt'
        etr = pvt%r(pvt%n-ncon+i)
        call graph_connect(etr,n,nnzs,n_col,l_ptr,l_col,con_c, &
            skip = skip)

        ! Speed up big systems as we constantly add
        ! elements to the skip table, and only push new values
        call rgn_consecutive_insert(skip, con_c)
        if ( .not. rgn_push(con, con_c) ) call die('level_struct: push -- 2')

      end do

    end do

    ! Clean up
    call rgn_delete(con,con_c,skip)

  end subroutine level_struct

  ! Extract only a certain level from the level structure
  subroutine lvl_struct_extract(pvt,lvl,ilvl,r)
    type(tRgn), intent(in) :: pvt, lvl
    integer, intent(in) :: ilvl
    type(tRgn), intent(inout) :: r
    integer :: i
    logical :: suc

    ! Count number of entries for the same level
    i = count(lvl%r(1:lvl%n) == ilvl)

    call rgn_init(r,i)
    r%n = 0

    ! Copy over the entries
    do i = 1 , lvl%n
      if ( lvl%r(i) == ilvl ) then
        suc = rgn_push(r,pvt%r(i))
      end if
    end do

  end subroutine lvl_struct_extract

  subroutine pseudo_peripheral(method,n,nnzs,n_col,l_ptr,l_col,sub,lvs,&
      small_wd, priority)
    ! The method for the degree search
    integer, intent(in) :: method
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The level tree
    type(tPvtLvl), intent(inout) :: lvs
    ! The user can optionally obtain the smallest width and degree of
    ! the sets created by the S set in the found level set
    ! It also returns the entry that creates that structure
    integer, intent(out), optional :: small_wd(3)
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! Local variables
    type(tRgn) :: pvt, lvl, S
    integer :: idx, i, etr, j
    integer :: depth

    ! Algorithm I
    !   (i) -- pick arbitrary node with minimal degree
    idx = idx_degree(method,n,nnzs,n_col,l_ptr,l_col,sub, priority = priority)
    !   (ii/iii) -- search for deepest level structure in these peripherals
#ifdef PVT_DEBUG
    write(*,*)'   first v set: ',sub%r(idx), idx
#endif

    call level_struct(n,nnzs,n_col,l_ptr,l_col,idx,sub,lvs%pvt,lvs%lvl,priority)
    search_deepest: do 

      depth = lvl_depth(lvs%lvl)
      if ( present(small_wd) ) then
        !  (iv) -- collect the depth of the smallest width structure
        small_wd(1) = huge(1)
        small_wd(2) = 0
        small_wd(3) = 0
      end if

      ! extract the nodes at the last level to search for deeper ones
      call lvl_struct_extract(lvs%pvt,lvs%lvl,depth,S)

      ! Sort by the degree of the nodes (lowest to highest)
      call sort_degree(method,n,nnzs,n_col,l_ptr,l_col,S,pvt)
      ! To get the correct index look up
      S%r = rgn_pivot(sub,pvt%r)

      do i = 1 , S%n
        etr = S%r(i)
#ifdef PVT_DEBUG
        write(*,*)'     analyzing degree set: ',i,etr
#endif

        call level_struct(n,nnzs,n_col,l_ptr,l_col,etr,sub,pvt,lvl,priority)

        if ( present(small_wd) ) then
          j = lvl_width(lvl)
          if ( j < small_wd(1) ) then
            small_wd(1) = j
            j = lvl_depth(lvl)
            small_wd(2) = j
            small_wd(3) = etr
          else
            j = lvl_depth(lvl)
          end if
        else
          j = lvl_depth(lvl)
        end if

        ! If the level depth is larger than the current start over 
        ! from that index (the last added element *must* be
        ! the largest level)
        if ( depth < j ) then
#ifdef PVT_DEBUG
          write(*,*)'   changing v set: ',pvt%r(1)
#endif

          call rgn_copy(pvt,lvs%pvt)
          call rgn_copy(lvl,lvs%lvl)
          cycle search_deepest
        end if

      end do

      exit search_deepest

    end do search_deepest

    call rgn_delete(pvt,lvl,S)

  end subroutine pseudo_peripheral

  subroutine sort_degree(method,n,nnzs,n_col,l_ptr,l_col,sub,sub_sort,priority)
    ! The method of sorting 
    integer, intent(in) :: method
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The sub sorted
    type(tRgn), intent(inout) :: sub_sort
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! local variables
    type(tRgn) :: s
    integer :: i, idx, etr, nsub
    logical :: suc

    if ( sub%n <= 1 ) then
      call rgn_copy(sub,sub_sort)
      return
    end if

    ! Make copy of the sub, this is required because we can then
    ! easily POP the elements
    call rgn_copy(sub,s)

    ! This enables sub and sub_sort to be the same
    nsub = sub%n

    ! initialize the pivoting array
    call rgn_init(sub_sort,nsub)
    sub_sort%n = 0

    do i = 1 , nsub

      idx = idx_degree(method,n,nnzs,n_col,l_ptr,l_col,s, priority = priority)
      etr = rgn_pop(s, idx)
      suc = rgn_push(sub_sort,etr)
      if ( etr /= sub_sort%r(sub_sort%n) ) &
          call die('sort_degree: Error in popping')

    end do

    call rgn_delete(s)

  end subroutine sort_degree

  subroutine Cuthill_Mckee(n,nnzs,n_col,l_ptr,l_col,sub,pvt,&
      start,priority, only_sub)
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The currently indexs of the pivoted arrays
    type(tRgn), intent(inout) :: pvt
    ! The algorithm performs best if it knows where to start
    ! hence we can force the algorithm to start from some point
    type(tRgn), intent(in), optional :: start
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)
    logical, intent(in), optional :: only_sub

    ! The queue list
    type(tRgn) :: Q
    ! Temporary region used to contain the connectivity graph
    type(tRgn) :: con, skip

    ! local variables
    integer :: i, etr, idx
    logical :: suc, lonly_sub
#ifdef PVT_DEBUG
    integer :: add_zero
#endif

#ifdef PVT_DEBUG
    add_zero = 0
    write(*,fmt_debug1_a) 'Cuthill-Mckee START'
#endif

    lonly_sub = .false.
    if ( present(only_sub) ) lonly_sub = only_sub

    ! prepare lists used for the algorithm
    call rgn_init(Q,sub%n)
    Q%n = 0

    ! When we have a start immediately add the starting elements
    if ( present(start) ) then
      ! Sort the starting elements according to the lowest degree
      call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,start,skip, priority = priority)
      ! Push the starting points into the queue (sorted with priority and
      ! lowest degree)
      if ( .not. rgn_push(Q, skip) ) call die('Error in CM -- 1')
    end if

    ! initialize the pivoting array
    call rgn_init(pvt, sub%n)
    pvt%n = 0

    ! the pivoting array + the queue array, this ensures that
    ! we do not back-track already processed elements
    call rgn_range(con, 1, n)
    call rgn_complement(sub, con, con)

    ! now create the initial region used for skipping connectivities
    if ( con%n > Q%n ) then
      call rgn_init_consecutive(n, skip, con)
      call rgn_consecutive_insert(skip, Q)
    else
      call rgn_init_consecutive(n, skip, Q)
      call rgn_consecutive_insert(skip, con)
    end if

#ifdef PVT_DEBUG
    write(*,fmt_debug2_ai) 'starting queue [size]', Q%n
#endif

    ! Initialize connectivity region for graph_connect
    do while ( pvt%n < sub%n )

      ! 1. If the queue is empty we add the one with the lowest
      !    degree
      if ( Q%n == 0 ) then

        ! We will not follow the graph
        if ( lonly_sub ) exit

        ! Fake the skipping region
        ! This is just to skip the first couple of zeros
        call rgn_delete(con)
        do idx = 1, skip%n
          if ( skip%r(idx) > 0 ) then
            con%n = skip%n - idx + 1
            con%r => skip%r(idx:)
            exit
          end if
        end do
        idx = idx_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,sub, skip = con, &
            priority = priority)
        call rgn_nullify(con)
        etr = sub%r(idx)
        
#ifdef PVT_DEBUG
        write(*,fmt_debug2_ai) 'adding to queue (non-connected) [element]', etr
#endif
        ! Push the queue and the skip table
        suc = rgn_push(Q, etr)
        call rgn_consecutive_insert(skip, etr)

      end if

      ! We find the one in the queue with the highest priority
      etr = rgn_pop(Q)

      ! 2. Add it to the pivoting table
      if ( .not. rgn_push(pvt, etr) ) call die('Error in CM push -- 1')

      ! Since the sort-degree is sorting into the same array, we have
      ! to ensure con to have the correct size
      call rgn_grow(con, sub%n - pvt%n)

      ! 3. Create the connectivity graph from idx (this will remove "back" 
      !    connected entries, hence no dublicates needs to be taken into 
      !    account.)
      call graph_connect(etr,n,nnzs,n_col,l_ptr,l_col, con, skip = skip )
      call rgn_consecutive_insert(skip, con)

      ! 4. Sort the connecting elements from lowest degree to highest degree
      call sort_degree(D_LOW,n,nnzs,n_col,l_ptr,l_col,con,con, priority = priority)

#ifdef PVT_DEBUG
      if ( con%n == 0 ) then
        add_zero = add_zero + 1
        if ( add_zero > 1 ) then
          if ( mod(add_zero, 10) == 0 ) then
            write(*,'(a)', advance='NO') '|'
          else
            write(*,'(a)', advance='NO') '.'
          end if
        else
          write(*,fmt_debug2_a,advance='NO') 'empty graph .'
        end if
      else
        if ( add_zero > 0 ) write(*,*) ! new-line
        add_zero = 0
        write(*,fmt_debug2_ai) 'adding to queue [size]', con%n
      end if
#endif

      ! 5. Add all connected entries to the queue in increasing order
      if ( .not. rgn_push(Q, con) ) call die('Error in CM push -- 2')

    end do

#ifdef PVT_DEBUG
    if ( add_zero > 0 ) write(*,*) ! new-line
#endif

    call rgn_delete(Q, con, skip)

    if ( pvt%n /= sub%n .and. .not. lonly_sub ) &
        call die('CM: Error in algorithm')

#ifdef PVT_DEBUG
    write(*,fmt_debug1_a) 'Cuthill-Mckee DONE'
#endif

  end subroutine Cuthill_Mckee

  subroutine rev_Cuthill_Mckee(n,nnzs,n_col,l_ptr,l_col,sub,pvt,start,priority, only_sub)
    integer, intent(in) :: n, nnzs
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(in) :: sub
    type(tRgn), intent(inout) :: pvt
    type(tRgn), intent(in), optional :: start
    integer, intent(in), optional :: priority(n)
    logical, intent(in), optional :: only_sub

    call Cuthill_Mckee(n,nnzs,n_col,l_ptr,l_col,sub,pvt, &
        start=start,priority=priority, only_sub=only_sub)

    call rgn_reverse(pvt)

  end subroutine rev_Cuthill_Mckee


  ! Sort according to "back-connections"
  subroutine sort_connect(method,n,nnzs,n_col,l_ptr,l_col, sub, neigh, sort, priority)
    ! The method of sorting 
    integer, intent(in) :: method
    ! the dimensionality of the system
    integer, intent(in) :: n, nnzs
    ! The sparse pattern
    integer, intent(in) :: n_col(n), l_ptr(n), l_col(nnzs)
    ! The region of interest
    type(tRgn), intent(in) :: sub
    ! The neighbouring region
    type(tRgn), intent(in) :: neigh
    ! The sub sorted
    type(tRgn), intent(inout) :: sort
    ! The priority of the rows, optional
    integer, intent(in), optional :: priority(n)

    ! local variables
    type(tRgn) :: connect, sort_neigh
    integer :: i, idx, ind, nel, j

    if ( sub%n <= 1 ) then

      call rgn_copy(sub, sort)
      return

    end if

    call rgn_copy(neigh, sort_neigh)
    call rgn_sort(sort_neigh)

    ! initialize the pivoting array
    call rgn_init(connect, sub%n)
    call rgn_init(sort, sub%n)
    sort%n = 0

    ! Loop all sub elements and determine the degree that connects
    ! to the neighbouring elements

    do i = 1 , sub%n

      connect%r(i) = 0
      do ind = l_ptr(sub%r(i)) + 1 , l_ptr(sub%r(i)) + n_col(sub%r(i))
        if ( in_rgn(sort_neigh, l_col(ind)) ) &
            connect%r(i) = connect%r(i) + 1
      end do

    end do

    call rgn_delete(sort_neigh)

    ! Now do sorting according to method
    select case ( method )
    case ( C_HIGH )

      do while ( rgn_size(sort) < rgn_size(sub) )

        ! Find all elements with the highest degree
        idx = maxloc(connect%r, 1)
        nel = connect%r(idx)
        connect%r(idx) = -1
        if ( .not. rgn_push(sort, sub%r(idx)) ) call die('Programming error SC -- 1')

        do j = idx + 1, sub%n
          if ( nel == connect%r(j) ) then
            connect%r(j) = -1
            if ( .not. rgn_push(sort, sub%r(j)) ) call die('Programming error SC -- 2')
          end if
        end do

      end do

    case ( C_LOW )

      do while ( rgn_size(sort) < rgn_size(sub) )

        ! Find all elements with the highest degree
        idx = minloc(connect%r, 1)
        nel = connect%r(idx)
        connect%r(idx) = huge(1)
        if ( .not. rgn_push(sort, sub%r(idx)) ) call die('Programming error SC -- 1')

        do j = idx + 1, sub%n
          if ( nel == connect%r(j) ) then
            connect%r(j) = huge(1)
            if ( .not. rgn_push(sort, sub%r(j)) ) call die('Programming error SC -- 2')
          end if
        end do

      end do

    end select

    call rgn_delete(connect)

  end subroutine sort_connect

  ! Returns all indices that connects to the 'idx'
  ! It does not add it-self to the list
  subroutine graph_connect(idx,n,nnzs,n_col,l_ptr,l_col,con,skip)
    integer, intent(in) :: idx, n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    type(tRgn), intent(inout) :: con
    type(tRgn), intent(in), optional :: skip

    ! local variables
    integer :: ind
    logical :: suc

    ! Reset connectivity graph
    call rgn_grow(con, n_col(idx))
    con%n = 0
    if ( n_col(idx) == 0 ) return

    if ( present(skip) ) then
      do ind = l_ptr(idx) + 1 , l_ptr(idx) + n_col(idx)
        if ( l_col(ind) == idx ) cycle ! on-site
        ! skip connect if already present
        if ( in_rgn(skip,l_col(ind)) ) cycle
        ! add to the list
        suc = rgn_push(con,l_col(ind))
      end do
    else
      do ind = l_ptr(idx) + 1 , l_ptr(idx) + n_col(idx)
        if ( l_col(ind) == idx ) cycle ! on-site
        ! add to the list
        suc = rgn_push(con,l_col(ind))
      end do
    end if

  end subroutine graph_connect


  ! Returns an index dependent on the method
  !   See the D_* variables in the top to see their meaning
  function idx_degree(method,n,nnzs,n_col,l_ptr,l_col,sub,skip,priority) result(idx)
    integer, intent(in) :: method
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    ! Only look at this sub-set of the full sparsity pattern
    type(tRgn), intent(in) :: sub
    ! Skip these entries, say if they already have been added
    ! or something else
    type(tRgn), intent(in), optional :: skip
    integer, intent(in), optional :: priority(n)
    ! Returned index for the degree function
    integer :: idx, etr
    type(tRgn) :: self

    ! Local variables
    integer :: i, deg, cdeg
    logical, allocatable :: lskip(:)

    select case ( method ) 
    case ( D_LOW , D_LOW_SUM )
      deg = huge(1)
    case ( D_HIGH , D_HIGH_SUM )
      deg = 0
    end select
    select case ( method ) 
    case ( D_LOW_SUM , D_HIGH_SUM )
      ! initialize for graph_connect
      ! Allocating for maximum spots is much faster than growing everytime
      call rgn_init(self,n)
    end select

    allocate(lskip(n))
    do i = 1, sub%n
      lskip(sub%r(i)) = .false.
    end do
    if ( present(skip) ) then
      do i = 1, skip%n
        lskip(skip%r(i)) = .true.
      end do
    end if

    idx = 0
    do i = 1 , sub%n
      etr = sub%r(i)

      ! check whether we should skip it
      if ( lskip(etr) ) cycle

      ! Get the current degree (dependent on the method)
      cdeg = degree(method,n,nnzs,n_col,l_ptr,l_col,etr,self)

      select case ( method )
      case ( D_LOW , D_LOW_SUM )
        if ( cdeg < deg ) then
          idx = i
          deg = cdeg
        else if ( cdeg == deg .and. present(priority) ) then
          ! ** this should never happen if idx == 0
          if ( priority(sub%r(idx)) < priority(etr) ) then
            idx = i
          end if
        end if
      case ( D_HIGH , D_HIGH_SUM )
        if ( deg < cdeg ) then
          idx = i
          deg = cdeg
        else if ( cdeg == deg .and. present(priority) ) then
          ! ** this should never happen if idx == 0
          if ( priority(sub%r(idx)) < priority(etr) ) then
            idx = i
          end if
        end if
      end select
    end do

    deallocate(lskip)
    call rgn_delete(self)

  end function idx_degree

  function degree(method,n,nnzs,n_col,l_ptr,l_col,etr,self) result(deg)
    integer, intent(in) :: method
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    integer, intent(in) :: etr ! the entry which we want the degree from
    type(tRgn), intent(inout) :: self ! used to create connectivity graph
    integer :: deg
    integer :: j

    select case ( method )
    case ( D_LOW , D_HIGH )
      deg = n_col(etr)
    case ( D_LOW_SUM , D_HIGH_SUM )
      ! Get the connectivity graph
      call graph_connect(etr,n,nnzs,n_col,l_ptr,l_col,self)
      ! calculate degree sum
      deg = 0
      do j = 1 , self%n
        deg = deg + n_col(self%r(j))
      end do
    end select

  end function degree


  subroutine sp2graphviz_sp(file,sp,types,method,pvt)
    use class_Sparsity, only: Sparsity, attach
    character(len=*), intent(in) :: file
    type(Sparsity), intent(inout) :: sp
    ! Methods applied
    integer, intent(in), optional :: types(:)
    integer, intent(in), optional :: method
    type(tRgn), intent(in), optional :: pvt

    integer :: n, n_nzs
    integer, pointer :: ncol(:), l_ptr(:), l_col(:)

    call attach(sp,nrows_g=n,nnzs=n_nzs,n_col=ncol,list_ptr=l_ptr,list_col=l_col)

    call sp2graphviz(file,n,n_nzs,ncol,l_ptr,l_col,types,method,pvt)

  end subroutine sp2graphviz_sp

  subroutine sp2graphviz_lists(file,n,nnzs,n_col,l_ptr,l_col,types,method,pvt)
    character(len=*), intent(in) :: file
    integer, intent(in) :: n, nnzs, n_col(n), l_ptr(n), l_col(nnzs)
    ! Methods applied
    integer, intent(in), optional :: types(n)
    integer, intent(in), optional :: method
    type(tRgn), intent(in), optional :: pvt

    !                                               gray      blue      red       green     purple   turkis    orange
    character(len=7), parameter :: edges(0:6)  = (/'#D8D8D8','#5882FA','#FA5858','#58FA58','#AC58FA','#9DFCF7','#FFD681'/)
    character(len=7), parameter :: colors(0:6) = (/'#000000','#2E2EFE','#FF0000','#00FF00','#A901DB','#19A099','#FFA600'/)
    character(len=4) :: con
    integer :: i, j, k, ind, lmethod, max_types
    ! Array to hold edges in each sub-group
    integer :: nedge
    integer, allocatable :: edge(:)

    lmethod = 1
    if ( present(method) ) lmethod = method

    open(unit=555,file=trim(file),form='formatted')

    if ( present(pvt) ) then
      if ( len_trim(pvt%name) > 0 ) then
        write(555,'(a)') '// Pivoting name: '//trim(pvt%name)
      end if
    end if

    ! Inform how to process this:
    write(555,'(a)') '// This command typically produces nice graphs'

    select case ( lmethod ) 
    case ( 1 ) ! GRAPH
      write(555,'(a)') '// neato -x '//trim(file)
      write(555,'(a)') '// neato -x -Tpdf '//trim(file)//' -o graph.pdf'
      write(555,'(a)') '// neato -x -Tpng '//trim(file)//' -o graph.png'
      write(555,'(a)') 'strict graph G {'
      con = ' -- '
    case ( 2 ) ! DI-GRAPH
      write(555,'(a)') '// dot -x '//trim(file)
      write(555,'(a)') '// dot -x -Tpdf '//trim(file)//' -o digraph.pdf'
      write(555,'(a)') '// dot -x -Tpng '//trim(file)//' -o digraph.png'
      write(555,'(a)') 'strict digraph G {'
      con = ' -> '
    end select

    ! Tell graphviz to print out the edges first
    write(555,'(a)') 'outputorder = edgesfirst;'
    ! show them in increasing rank-order, from left to right
    write(555,'(a)') 'rankdir = LR;'
    ! scale the diagram to show everything (makes the image very large)
    write(555,'(a)') 'overlap = scale;'
    ! try to decrease overlaps of edges with nodes, hence paths between nodes
    ! will be clearer
    write(555,'(a)') 'splines = spline;'

    if ( present(types) ) then
      max_types = maxval(types)
    end if

    ! Allocate number of edges
    allocate(edge(maxval(n_col)))

    if ( present(types) ) then
      if ( max_types > len(edges) - 1 ) then
        write(*,*) 'Max allowed types: ',len(edges)-1
        call die('Can not represent more than 6 different types.')
      end if
      do j = 0 , max_types
        select case ( j )
        case ( 0 , 1 )
          write(555,'(a)') 'node [shape=box style=filled fillcolor="'//trim(colors(j))//'" fontcolor=white]'
        case default
          write(555,'(a)') 'node [shape=box style=filled fillcolor="'//trim(colors(j))//'"]'
        end select
        do i = 1 , n
          if ( types(i) == j ) then
            write(555,'(tr1,i0)',advance='no') get_col(i,pvt)
          end if
        end do
        write(555,'(tr1,a)') ';'
      end do
    end if
    ! Default edge
    write(555,'(a)') 'edge [color=gray] ;'

    if ( present(types) ) then

      ! print everything with respect to types
      do i = 1 , n 
        ! Always show node (sometimes it could not connect)
        if ( types(i) < 0 ) cycle
        call print_node(i,pvt)
        if ( n_col(i) <= 1 ) cycle
        do k = 0 , max_types
          if ( k /= types(i) ) cycle
          ! Create all edges with type connections
          do j = 0 , max_types

            ! get edges for this node
            nedge = 0
            do ind = l_ptr(i) + 1 , l_ptr(i) + n_col(i)
              if ( l_col(ind) <= i ) cycle
              if ( j /= types(l_col(ind)) ) cycle
              nedge = nedge + 1
              edge(nedge) = l_col(ind)
            end do

            if ( nedge == 0 ) cycle

            ! print out the edges it connects to
            if ( k == j ) then
              write(555,'(3a)') '{ edge [color="',trim(edges(k)),'"]'
            else
              write(555,'(5a)') '{ edge [color="',trim(edges(j)),';0.5:', &
                  trim(edges(k)),'"]'
            end if
            write(555,'(i0,2a)',advance='no') get_col(i,pvt),con,' {'
            do ind = 1 , nedge
              write(555,'(tr1,i0)',advance='no') get_col(edge(ind),pvt)
            end do
            write(555,'(a)') '} }'
          end do
        end do
      end do

      do k = 1 , max_types
        ! Rank them similarly
        write(555,'(tr1,a)',advance='no') '{ rank=same;'
        do i = 1 , n
          if ( types(i) == k ) then
            write(555,'(tr1,i0)',advance='no') get_col(i,pvt)
          end if
        end do
        write(555,'(tr1,a)') '}'
      end do

      ! Rank them similarly
      if ( any(types < 0) ) then
        write(555,'(tr1,a)') 'subgraph BUF {'
        write(555,'(tr1,a)') 'rank=same;'
        write(555,'(a)') 'node [shape=box style=filled fillcolor=white fontcolor=black]'
        do i = 1 , n
          if ( types(i) < 0 ) then
            write(555,'(i0)') get_col(i,pvt)
          end if
        end do
        write(555,'(tr1,a)') '}'
      end if

    else

      do i = 1 , n
        call print_node(i,pvt)
        if ( n_col(i) <= 1 ) cycle
        ! Write out the label of the node
        write(555,'(i0,2a)',advance='no') get_col(i,pvt),con,' {'
        do ind = l_ptr(i) + 1 , l_ptr(i) + n_col(i)
          if ( l_col(ind) == i ) cycle
          write(555,'(tr1,i0)',advance='no') get_col(l_col(ind),pvt)
        end do
        write(555,'(a)') '}'
      end do

    end if

    deallocate(edge)

    write(555,'(a)') '}'
    close(555)

  contains

    function get_col(col,pvt)
      integer, intent(in) :: col
      type(tRgn), intent(in), optional :: pvt
      integer :: get_col
      if ( present(pvt) ) then
        get_col = rgn_pivot(pvt,col)
      else
        get_col = col
      end if
    end function get_col

    subroutine print_node(i,pvt)
      integer, intent(in) :: i
      type(tRgn), intent(in), optional :: pvt
      ! Write out the label of the node
      write(555,'(i0,'' [label="'',i0,'' ('',i0,'')"]'')') &
          get_col(i,pvt),get_col(i,pvt),n_col(i)-1
    end subroutine print_node

  end subroutine sp2graphviz_lists

end module m_pivot_methods
