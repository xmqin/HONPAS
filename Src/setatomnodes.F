! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine setatomnodes(na,lasto,Node,Nodes)
C
C  Sets up the distribution of atoms over nodes for a 
C  spatial decomposition in parallel.
C
C  On entry :
C
C  na            = total number of atoms
C  lasto         = pointer to atom orbitals
C  Node          = local node number
C  Nodes         = total number of nodes
C
C  On exit (in modules) :
C
C  ncellpernode               = number of cells on the current node
C  ncellnodeptr(ncellpernode) = pointer to cell indices for each local cell
C  lbuffercell(ncellpernode)  = if .true. this is a buffer region cell
C  nspmin(3)                  = lower bound to real cells (not buffers)
C  nspmax(3)                  = upper bound to real cells (not buffers)
C  natomsG2L(na)              = pointer from global atom no. to local one
C  nG2L(na)                   = pointer from global orbital no. to local one
C  natomsL2G(na)              = pointer from local atom no. to global one
C  nL2G(no, nodes)            = pointer from local orbital to global one
C  natomsNode(na)             = pointer to node responsible for each atom
C  nNode(no)                  = pointer to node responsible for each orbital
C  nAtomPerNode(nodes)        = no. atoms per node
C  maxAtomPerNode             = max(nAtomsPerNode)
C  nOrbPerNode(nodes)         = no. orbitals per node
C  nCellPerNode(nodes)        = no. cells per node
C  maxCellPerNode             = max(nCellPerNode)
C 
C  Julian Gale, NRI, Curtin University, March 2004
C
C
      use alloc,     only : re_alloc, de_alloc
      use sys,       only : die
      use parallel,  only : ionode
      use precision, only : dp
      use spatial,   only : nAtomPerNode, maxAtomPerNode
     .                    , nCellPerNode, maxCellPerNode
     .                    , natomsL2G, natomsG2L, nL2G, nG2L
     .                    , nbufferx, nbuffery, nbufferz
     .                    , ncellnodeptr
     .                    , lbuffercell
     .                    , nAtomsNode, nNode, nOrbPerNode
     .                    , npgrid
     .                    , nspcell, nspcellat, nspcellatptr
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world
      use mpi_siesta, only: mpi_integer, mpi_max
#endif
      implicit none
C
C  Passed variables
C
      integer, intent(in)                      :: na
      integer, intent(in)                      :: lasto(0:na)
      integer, intent(in)                      :: Node
      integer, intent(in)                      :: Nodes
C
C  Local variables
C
      integer :: no_u, iu, LOrb, Gorb
      integer                                  :: i
      integer                                  :: ia
      integer                                  :: icx
      integer                                  :: icy
      integer                                  :: icz
      integer                                  :: ii
      integer                                  :: ind
      integer                                  :: ix
      integer                                  :: iy
      integer                                  :: iz
      integer                                  :: j
      integer                                  :: maxorb
      integer                                  :: maxxy
      integer                                  :: maxx
      integer                                  :: n
      integer                                  :: n2
      integer                                  :: n3
      integer                                  :: n5
      integer                                  :: n235(3)
      integer                                  :: nn235(3)
      integer                                  :: nadd
      integer                                  :: natomnow
      integer                                  :: natomold
      integer, pointer                         :: npgridxptr(:)
      integer, pointer                         :: npgridyptr(:)
      integer, pointer                         :: npgridzptr(:)
      integer                                  :: nspmax(3)
      integer                                  :: nspmin(3)
      integer                                  :: nspdiff
      integer                                  :: nspnobuff
      logical                                  :: ldebug
      logical                                  :: lfixx
      logical                                  :: lfixy
      logical                                  :: lfixz
      logical                                  :: lnon235
      logical                                  :: lok
      logical                                  :: lokx
      logical                                  :: loky
      logical                                  :: lokz
      logical                                  :: lpgridinput
      real(dp)                                 :: diffxy
      real(dp)                                 :: diffxz
      real(dp)                                 :: diffyz
      real(dp)                                 :: pratio
      real(dp)                                 :: ratiox
      real(dp)                                 :: ratioy
      real(dp)                                 :: ratioz
      real(dp)                                 :: ratioxy
      real(dp)                                 :: ratioxz
      real(dp)                                 :: ratioyz
      real(dp)                                 :: rnow
      real(dp)                                 :: spcelltot
      real(dp)                                 :: targetx
      real(dp)                                 :: targety
      real(dp)                                 :: targetz
      real(dp)                                 :: targetratioxy
      real(dp)                                 :: targetratioxz
      real(dp)                                 :: targetratioyz
      integer, dimension(:), pointer           :: ncellpernodelist
      integer, dimension(:), pointer           :: natompernodelist
#ifdef MPI
      integer                                  :: MPIerror
      integer, dimension(:), pointer           :: ntmp
#endif
C

      if (Nodes.eq.1) then
C**********************
C  Non-parallel case  *
C**********************
        natompernode = na
        if (natompernode.gt.maxatompernode) then
          maxatompernode = natompernode
          call re_alloc( natomsL2G, 1, maxatompernode, 'natomsL2G',
     .                   'setatomnodes' )
        endif
        do i = 1,natompernode
           natomsL2G(i) = i
        enddo
C
C  For spatial algorithm set pointers to cells
C
        ncellpernode = nspcell(1)*nspcell(2)*nspcell(3)
        if (ncellpernode.gt.maxcellpernode) then
          maxcellpernode = ncellpernode
          call re_alloc( lbuffercell, 1, maxcellpernode,
     .                   'lbuffercell', 'setatomnodes' )
          call re_alloc( ncellnodeptr, 1, maxcellpernode,
     .                   'ncellnodeptr', 'setatomnodes' )
        endif
        nspmax(1) = nspcell(1) - nbufferx
        nspmin(1) = nbufferx
        nspmax(2) = nspcell(2) - nbuffery
        nspmin(2) = nbuffery
        nspmax(3) = nspcell(3) - nbufferz
        nspmin(3) = nbufferz
        maxxy = nspcell(1)*nspcell(2)
        maxx  = nspcell(1)
        n = 0
        do iz = 1,nspcell(3)
          lokz = (iz.gt.nbufferz.and.iz.lt.(nspcell(3)-nbufferz+1))
          do iy = 1,nspcell(2)
            loky = (iy.gt.nbuffery.and.iy.lt.(nspcell(2)-nbuffery+1))
            do ix = 1,nspcell(1)
              lokx = (ix.gt.nbufferx.and.ix.lt.(nspcell(1)-nbufferx+1))
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              n = n + 1
              ncellnodeptr(n) = ind
              if (lokx.and.loky.and.lokz) then
                lbuffercell(n) = .false.
              else
                lbuffercell(n) = .true.
              endif
            enddo
          enddo
        enddo
C
C  Build lists linking atoms with parallel structure
C
        call re_alloc( natomsG2L, 1, na, 'natomsG2L', 'setatomnodes' )
        call re_alloc( natomsNode, 1, na, 'natomsNode', 'setatomnodes' )
        natomsNode(1:na) = 0
        do ii = 1,na
          natomsG2L(ii) = ii
        enddo
C
C  Build dummy lists linking orbitals with parallel structure 
C
        maxorb = lasto(na)
        call re_alloc( nL2G, 1, maxorb,1, Nodes, 'nL2G',
     &                 'setatomnodes' )
        call re_alloc( nG2L, 1, maxorb, 'nG2L', 'setatomnodes' )
        call re_alloc( nNode, 1, maxorb, 'nNode', 'setatomnodes')
        call re_alloc( nOrbPerNode, 1, Nodes, 'nOrbPerNode',
     &                 'setatomnodes' )
C
        nOrbPerNode(1) = maxorb
        do ii = 1,maxorb
          nL2G(ii,1) = ii
          nG2L(ii) = ii
          nNode(ii) = 0
        enddo
      else
C**************************
C  Spatial parallel case  *
C**************************
C
C  Set debugging flag
C
        ldebug = .false.
c        ldebug = (Node.eq.0)
C
C  Has the processor grid been previously specifed?
C
        lpgridinput = (npgrid(1)*npgrid(2)*npgrid(3).eq.Nodes)
C
C  Processor grid has already been set
C
        if (lpgridinput) goto 999
C
C  Find factors of number of processors (2,3,5)
C
        n2 = 0
        n3 = 0
        n5 = 0
        n = Nodes
        lok = .true.
        do while (lok.and.n.gt.1)
          if (mod(n,2).eq.0) then
            n2 = n2 + 1
            n = n/2
          elseif (mod(n,3).eq.0) then
            n3 = n3 + 1
            n = n/3
          elseif (mod(n,5).eq.0) then
            n5 = n5 + 1
            n = n/5
          else
            lok = .false.
          endif
        enddo
C
C  Initialise npgrid
C
        npgrid(1:3) = 1
C
C  Set flag to indicate whether there is a factor other than 2,3 or 5
C
        lnon235 = (n.ne.1)
        if (ldebug) then
          write(6,'(/,''  Spatial decomposition of processors :'',/)')
          write(6,'(''  No. of multiples of 2 = '',i4)') n2
          write(6,'(''  No. of multiples of 3 = '',i4)') n3
          write(6,'(''  No. of multiples of 5 = '',i4)') n5
          if (lnon235) then
            write(6,'(''  Other factor          = '',i4)') n
          endif
        endif
C
C  If non 2,3,5, decide which direction to fix
C
        lfixx = .false.
        lfixy = .false.
        lfixz = .false.
        if (lnon235) then
          spcelltot = (nspcell(1)-2*nbufferx)*
     .                (nspcell(2)-2*nbuffery)*
     .                (nspcell(3)-2*nbufferz)
          pratio = real(n,dp)/real(Nodes,dp)
          ratiox = pratio - (nspcell(1)-2*nbufferx)/spcelltot 
          ratioy = pratio - (nspcell(2)-2*nbuffery)/spcelltot
          ratioz = pratio - (nspcell(3)-2*nbufferz)/spcelltot
          ratiox = abs(ratiox)
          ratioy = abs(ratioy)
          ratioz = abs(ratioz)
          if (ratiox.lt.ratioy.and.ratiox.lt.ratioz) then
            lfixx = .true.
            npgrid(1) = n
          elseif (ratioy.lt.ratioz) then
            lfixy = .true.
            npgrid(2) = n
          else
            lfixz = .true.
            npgrid(3) = n
          endif
        endif
C
C  Calculate optimal target ratios for x-y and x-z
C
        targetratioxy = real(nspcell(2)-2*nbuffery,dp)/
     .    real(nspcell(1)-2*nbufferx,dp)
        targetratioxz = real(nspcell(3)-2*nbufferz,dp)/
     .    real(nspcell(1)-2*nbufferx,dp)
        targetratioyz = real(nspcell(3)-2*nbufferz,dp)/
     .    real(nspcell(2)-2*nbuffery,dp)
C
C  Group factors in proportion to numbers of cells
C
        n235(1) = n5
        n235(2) = n3
        n235(3) = n2
        nn235(1) = 5
        nn235(2) = 3
        nn235(3) = 2
        do ind = 1,3
          ii = n235(ind)
          n = nn235(ind)
          do i = 1,ii
            if (lfixx) then
              ratioyz = real(npgrid(3),dp)/real(npgrid(2),dp)
              diffyz = targetratioyz - ratioyz
              if (diffyz.lt.0.0_dp) then
                npgrid(2) = n*npgrid(2)
              else
                npgrid(3) = n*npgrid(3)
              endif
            elseif (lfixy) then
              ratioxz = real(npgrid(3),dp)/real(npgrid(1),dp)
              diffxz = targetratioxz - ratioxz
              if (diffxz.lt.0.0_dp) then
                npgrid(1) = n*npgrid(1)
              else
                npgrid(3) = n*npgrid(3)
              endif
            elseif (lfixz) then
              ratioxy = real(npgrid(2),dp)/real(npgrid(1),dp)
              diffxy = targetratioxy - ratioxy
              if (diffxy.lt.0.0_dp) then
                npgrid(1) = n*npgrid(1)
              else
                npgrid(2) = n*npgrid(2)
              endif
            else
              ratioxy = real(npgrid(2),dp)/real(npgrid(1),dp)
              ratioxz = real(npgrid(3),dp)/real(npgrid(1),dp)
              diffxy = targetratioxy - ratioxy
              diffxz = targetratioxz - ratioxz
              if (diffxy.lt.0.0_dp.and.diffxz.lt.0.0_dp) then
                npgrid(1) = n*npgrid(1)
              elseif (diffxy.gt.diffxz) then
                npgrid(2) = n*npgrid(2)
              else
                npgrid(3) = n*npgrid(3)
              endif
            endif
          enddo
        enddo

 999    continue

        if (ldebug) then
          write(6,'(/,''  Processor grid = '',3(i4,1x))') 
     .      (npgrid(n),n=1,3)
        endif
C
C  Calculate target number of atoms per node in each direction
C
        targetx = real(na,dp)/real(npgrid(1),dp)
        targety = real(na,dp)/real(npgrid(2),dp)
        targetz = real(na,dp)/real(npgrid(3),dp)
C
C  Allocate local memory
C
        nullify( npgridxptr )
        call re_alloc( npgridxptr, 0, npgrid(1), 'npgridxptr',
     &                 'setatomnodes' )
        nullify( npgridyptr )
        call re_alloc( npgridyptr, 0, npgrid(2), 'npgridyptr',
     &                 'setatomnodes' )
        nullify( npgridzptr )
        call re_alloc( npgridzptr, 0, npgrid(3), 'npgridzptr',
     &                 'setatomnodes' )
C
C
C Now count cells, and divide them in x, y, and z 
C according to number of nodes available in each direction

        npgridxptr(0) = nbufferx
        npgridyptr(0) = nbuffery
        npgridzptr(0) = nbufferz
C
        maxxy = nspcell(1)*nspcell(2)
        maxx  = nspcell(1)

C  Assign cells to nodes based on cumulative numbers of atoms : X
C
        natomnow = 0
        natomold = 0
        npgridxptr(1:npgrid(1)) = 0
        ii = 1
        do ix = nbufferx+1,nspcell(1) - nbufferx
          if (ii == npgrid(1)) exit
          ! If we've reached the last cell, then all other atoms must
          ! go onto next node. Go over each plane of x
          do iy = nbuffery+1,nspcell(2) - nbuffery
            do iz = nbufferz+1,nspcell(3) - nbufferz
              ! Within each plane of x, calculate the index of each cell
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              ! and accumulate the number of atoms in each cell
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          ! natomnow contains all atoms up till now
          rnow = natomnow - natomold
          ! rnow contains all atoms on this plane of x, plus
          ! any previous planes not yet assigned to a node.
          if (rnow/targetx > 0.90_dp) then
            npgridxptr(ii) = ix
            ii = ii + 1
            natomold = natomnow
          endif
        enddo
       ! Put all remaining atoms onto next node; in case we still
       ! have nodes left over, leave them empty.
        npgridxptr(ii:npgrid(1)) = nspcell(1) - nbufferx
       ! npgridxptr contains, indexed by processor grid in the x
       ! direction, the first x cell that processor is responsible for.

C Now repeat the above algorithm for y & z directions.
C
C  Assign cells to nodes based on cumulative numbers of atoms : Y
C
        natomnow = 0
        natomold = 0
        npgridyptr(1:npgrid(2)) = 0
        ii = 1
        do iy = nbuffery+1,nspcell(2) - nbuffery
          if (ii == npgrid(2)) exit
          do iz = nbufferz+1,nspcell(3) - nbufferz
            do ix = nbufferx+1,nspcell(1) - nbufferx
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          rnow = natomnow - natomold
          if (rnow/targety > 0.90_dp) then
            npgridyptr(ii) = iy
            ii = ii + 1
            natomold = natomnow
          endif
        enddo
        npgridyptr(ii:npgrid(2)) = nspcell(2) - nbuffery
C
C  Assign cells to nodes based on cumulative numbers of atoms : Z
C
        natomnow = 0
        natomold = 0
        npgridzptr(1:npgrid(3)) = 0
        ii = 1
        do iz = nbufferz+1,nspcell(3) - nbufferz
          if (ii == npgrid(3)) exit
          do ix = nbufferx+1,nspcell(1) - nbufferx
            do iy = nbuffery+1,nspcell(2) - nbuffery
              ind = (iz-1)*maxxy + (iy-1)*maxx + ix
              natomnow = natomnow + nspcellat(ind)
            enddo
          enddo
          rnow = natomnow - natomold
          if (rnow/targetz > 0.90_dp) then
            npgridzptr(ii) = iz
            ii = ii + 1
            natomold = natomnow
          endif
        enddo
        npgridzptr(ii:npgrid(3)) = nspcell(3) - nbufferz
C
C  Build lists linking cells to nodes
C
        n = 0
        ncellpernode = 0
        do iz = 1,npgrid(3)
          do iy = 1,npgrid(2)
            do ix = 1,npgrid(1)
              n = n + 1
              if (Node.eq.n-1) then
C Calculate maximum & minimum number of cells that this node
C will get in each direction:
                nspmax(1) = npgridxptr(ix) 
                nspmin(1) = npgridxptr(ix-1)
                nspmax(2) = npgridyptr(iy)
                nspmin(2) = npgridyptr(iy-1)
                nspmax(3) = npgridzptr(iz)
                nspmin(3) = npgridzptr(iz-1)
C
C Adjust these maxima & minima for buffering
                nspnobuff = nspcell(1) - 2*nbufferx
                nspdiff = min(nspnobuff-nspmax(1)+nspmin(1),nbufferx)
                if (nspdiff.gt.0) nspmin(1) = nspmin(1) - nspdiff
                nspdiff = min(nspnobuff-nspmax(1)+nspmin(1),nbufferx)
                if (nspdiff.gt.0) nspmax(1) = nspmax(1) + nspdiff
                nspnobuff = nspcell(2) - 2*nbuffery
                nspdiff = min(nspnobuff-nspmax(2)+nspmin(2),nbuffery)
                if (nspdiff.gt.0) nspmin(2) = nspmin(2) - nspdiff
                nspdiff = min(nspnobuff-nspmax(2)+nspmin(2),nbuffery)
                if (nspdiff.gt.0) nspmax(2) = nspmax(2) + nspdiff
                nspnobuff = nspcell(3) - 2*nbufferz
                nspdiff = min(nspnobuff-nspmax(3)+nspmin(3),nbufferz)
                if (nspdiff.gt.0) nspmin(3) = nspmin(3) - nspdiff
                nspdiff = min(nspnobuff-nspmax(3)+nspmin(3),nbufferz)
                if (nspdiff.gt.0) nspmax(3) = nspmax(3) + nspdiff
C
C Ensure these maxima & minima are not outside the bounds of the cells.
                nspmax(1) = min(nspmax(1),nspcell(1)-1)
                nspmin(1) = max(nspmin(1),1)
                nspmax(2) = min(nspmax(2),nspcell(2)-1)
                nspmin(2) = max(nspmin(2),1)
                nspmax(3) = min(nspmax(3),nspcell(3)-1)
                nspmin(3) = max(nspmin(3),1)
C
C nadd is the number of cells we've allocated for this node
                nadd = (npgridzptr(iz) - npgridzptr(iz-1) + 2*nbufferz)*
     .                 (npgridyptr(iy) - npgridyptr(iy-1) + 2*nbuffery)*
     .                 (npgridxptr(ix) - npgridxptr(ix-1) + 2*nbufferx)
                if (ncellpernode+nadd.gt.maxcellpernode) then
                  maxcellpernode = ncellpernode + nadd
                  call re_alloc( lbuffercell, 1, maxcellpernode,
     .                           'lbuffercell', 'setatomnodes' )
                  call re_alloc( ncellnodeptr,1,maxcellpernode,
     .                           'ncellnodeptr', 'setatomnodes' )
                endif
C
C The lok[xyz] variable below are only necessary when using buffers.
                do icz = npgridzptr(iz-1)+1-nbufferz,
     .              npgridzptr(iz)+nbufferz
                  lokz = (icz.gt.npgridzptr(iz-1).and.
     .              icz.le.npgridzptr(iz))
                  do icy = npgridyptr(iy-1)+1-nbuffery,
     .                npgridyptr(iy)+nbuffery
                    loky = (icy.gt.npgridyptr(iy-1).and.
     .                icy.le.npgridyptr(iy))
                    do icx = npgridxptr(ix-1)+1-nbufferx,
     .                  npgridxptr(ix)+nbufferx
                      lokx = (icx.gt.npgridxptr(ix-1).and.
     .                  icx.le.npgridxptr(ix))
                      ind = (icz-1)*maxxy + (icy-1)*maxx + icx
                      ncellpernode = ncellpernode + 1
                      ncellnodeptr(ncellpernode) = ind
                      if (lokx.and.loky.and.lokz) then
                        lbuffercell(ncellpernode) = .false.
                      else
                        lbuffercell(ncellpernode) = .true.
                      endif
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
        
        if (ionode) then
          nullify( ncellpernodelist )
          call re_alloc( ncellpernodelist, 1, nodes , 
     &                   'ncellpernodelist', 'setatomnodes' )
        else
          nullify( ncellpernodelist )
          call re_alloc( ncellpernodelist, 0, 0, 
     &                   'ncellpernodelist', 'setatomnodes' )
        endif
#ifdef MPI
        call mpi_gather(ncellpernode,1,mpi_integer,
     .       ncellpernodelist,1,mpi_integer,
     .       0,mpi_comm_world,mpierror)
#else
        ncellpernodelist=ncellpernode
#endif
        if (ionode) then
           write(6,*)
           do i=0,nodes-1
              write(6,'(''  Cells per Processor = '',2(i8,1x))') 
     .             i,ncellpernodelist(i+1)
           enddo
           if (any(ncellpernodelist==0)) then
              write(6,*)
              write(6,*) "Warning: Bad load balancing: ",
     .                   "no cells allocated to node",
     .                   minloc(ncellpernodelist)-1
              write(6,*) "Try specifying a number of nodes that is ",
     .                   "an exact factor of the number of cells:"
     .                   ,product(nspcell)
              call die()
           endif
        endif
        call de_alloc( ncellpernodelist, 'ncellpernodelist',
     &                 'setatomnodes' )
C
C  Find size for atom link array and then allocate
C
        natompernode = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            natompernode = natompernode + nspcellat(ind)
          endif
        enddo
        if (natompernode.gt.maxatompernode) then
          maxatompernode = natompernode
          call re_alloc( natomsL2G, 1, maxatompernode, 'natomsL2G',
     &                   'setatomnodes' )
        endif
C
C  Build lists linking atoms to nodes
C
        natompernode = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              natompernode = natompernode + 1
              natomsL2G(natompernode) = nspcellatptr(n,ind)
            enddo
          endif
        enddo

        if (ionode) then
          nullify( natompernodelist )
          call re_alloc( natompernodelist, 1, nodes , 
     &                   'natompernodelist', 'setatomnodes' )
        else
          nullify( natompernodelist )
          call re_alloc( natompernodelist, 0, 0 , 
     &                   'natompernodelist', 'setatomnodes' )
        endif
#ifdef MPI
        call mpi_gather(natompernode,1,mpi_integer,
     .       natompernodelist,1,mpi_integer,
     .       0,mpi_comm_world,mpierror)
#else
        natompernodelist=natompernode
#endif
        if (ionode) then
          write(6,*)
          do i = 0,nodes-1
            write(6,'(''  Atoms per Processor = '',2(i8,1x))') 
     .             i,natompernodelist(i+1)
          enddo
          if (any(natompernodelist==0)) then
             call die("Bad load balancing. Nodes with no atoms")
          endif
        endif
        call de_alloc( natompernodelist, 'natompernodelist',
     &                 'setatomnodes' )
C
C  Build lists linking atoms with parallel structure
C
        call re_alloc( natomsG2L, 1, na, 'natomsG2L', 'setatomnodes' )
        call re_alloc( natomsNode, 1, na, 'natomsNode', 'setatomnodes' )
        natomsNode(1:na) = 0
        natomsG2L(1:na) = 0
        j = 0
        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              ia = nspcellatptr(n,ind)
              j = j + 1
              natomsG2L(ia) = j
              natomsNode(ia) = Node
            enddo
          endif
        enddo
C
C  Globalise data
C
        maxorb = lasto(na)
#ifdef MPI
        nullify( ntmp )
        call re_alloc( ntmp, 1, max(maxorb,Nodes), 'ntmp',
     &                 'setatomnodes' )
        ntmp(1:na) = natomsNode(1:na)
        call MPI_AllReduce(ntmp,natomsNode,na,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        ntmp(1:na) = natomsG2L(1:na)
        call MPI_AllReduce(ntmp,natomsG2L,na,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
#endif
C
C  Build lists linking orbitals with parallel structure 
C
        ! Instead of maxorb, nL2G's first dimension could
        ! be max(no_l)
        call re_alloc( nL2G, 1, maxorb,1,Nodes, 'nL2G', 'setatomnodes' )
        call re_alloc( nG2L, 1, maxorb, 'nG2L', 'setatomnodes' )
        call re_alloc( nNode, 1, maxorb, 'nNode', 'setatomnodes' )
        call re_alloc( nOrbPerNode, 1, Nodes, 'nOrbPerNode',
     &                 'setatomnodes' )
C
        nOrbPerNode(1:Nodes) = 0
        nNode(1:maxorb) = 0
        nL2G(1:maxorb,1:Nodes) = 0
        nG2L(1:maxorb) = 0

        do ii = 1,ncellpernode
          if (.not.lbuffercell(ii)) then
            ind = ncellnodeptr(ii)
            do n = 1,nspcellat(ind)
              ia = nspcellatptr(n,ind)
              do j = lasto(ia-1)+1,lasto(ia)
                nOrbPerNode(Node+1) = nOrbPerNode(Node+1) + 1
                nL2G(nOrbPerNode(Node+1),Node+1) = j
                nG2L(j) = nOrbPerNode(Node+1)
                nNode(j) = Node
              enddo         
            enddo         
          endif
        enddo
C
C  Globalise data
C
#ifdef MPI
!
!       This is a crude use of allreduce with mpi_max.
!       It is based on the fact that "non-assigned" slots are zero,
!       so upon the reduction the global array will contain the
!       data from all processors. It would be more clear to use
!       mpi_all_gather...

        ntmp(1:maxorb) = nNode(1:maxorb)
        call MPI_AllReduce(ntmp,nNode,maxorb,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        ! 
        do ii = 1,Nodes
          ntmp(1:maxorb) = nL2G(1:maxorb,ii)
          call MPI_AllReduce(ntmp,nL2G(1,ii),maxorb,MPI_integer,
     .      MPI_max,MPI_Comm_World,MPIerror)
        enddo

        ! This globalization seems useless, since nG2L should
        ! only be used by the local node.
        ntmp(1:maxorb) = nG2L(1:maxorb)
        call MPI_AllReduce(ntmp,nG2L,maxorb,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)

        ntmp(1:Nodes) = nOrbPerNode(1:Nodes)
        call MPI_AllReduce(ntmp,nOrbPerNode,Nodes,MPI_integer,MPI_max,
     .    MPI_Comm_World,MPIerror)
        call de_alloc( ntmp,  'ntmp', 'setatomnodes' )
#endif
         no_u = maxorb
         if (node==0) then
            call io_assign(iu)
            open(unit=iu,file="SPATIAL_INDEXES",form="formatted")
            write(iu,*) "nl2g"
            do Lorb = 1, no_u	
               write(iu, "(i6,8i9)") Lorb, (nL2G(Lorb,j),j=1,Nodes)
            enddo
            write(iu,*) "nNode"
            do Gorb = 1, no_u	
               write(iu, "(i6,i5)") Gorb, nNode(Gorb)
            enddo
            write(iu,*) "nG2L (node 0)"
            do Gorb = 1, no_u	
               write(iu, "(i6,i5)") Gorb, nG2L(Gorb)
            enddo
            write(iu,*) "nOrbPerNode"
            do j = 1, Nodes
               write(iu, "(i6,i8)") j-1, nOrbPerNode(j)
            enddo
            call io_close(iu)
         endif


C
C  Free local memory
C
        call de_alloc( npgridxptr, 'npgridxptr', 'setatomnodes' )
        call de_alloc( npgridyptr, 'npgridyptr', 'setatomnodes' )
        call de_alloc( npgridzptr, 'npgridzptr', 'setatomnodes' )
      endif
      return
      end
