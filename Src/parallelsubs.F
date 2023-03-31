! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module parallelsubs
C
C  Parallelisation related subroutines
C
      use parallel,     only : BlockSize, ProcessorY
      use parallel,     only : thisNode => Node
      use spatial,      only : lspatial, nNode, nL2G, nG2L,
     &                         nOrbPerNode
      use domain_decom, only : use_dd, use_dd_perm, dd_nuo, dd_perm,
     &                         ulimit, llimit, dd_invp, dd_nnode
      use sys,          only : die
      implicit none

      logical, public, save :: pexsi_dist 
      integer, public, save :: pexsi_bs

      contains

      subroutine set_processorY(n)
      integer, intent(in) :: n
      ProcessorY = n
      end subroutine set_processorY

      subroutine GetNodeOrbs( NOrb, Node, Nodes, NOrbNode)
      use parallel, only: thisNode => Node
C
C  Calculates the number of orbitals stored on the local Node.
C
C  Julian Gale, October 1998
C
C  Input :
C
C  integer NOrb     = The total number of orbitals in the calculation
C  integer Node     = The local processor
C  integer Nodes    = The total number of processors
C
C  Output :
C
C  integer NOrbNode = The number of orbitals stored on this Node - if zero
C                     on input then calculated otherwise left unchanged
C
C  Modified so that value from last call is saved in case it can be
C  re-used on the next call to save re-calculation.
C

C Passed arguments
      integer NOrb, NOrbNode, Node, Nodes

C Local variables
      integer Remainder, MinPerNode, RemainderBlocks, NOrbLast,
     .  NOrbNodeLast

C Saved variables
      save NOrbLast, NOrbNodeLast

      data
     .  NOrbLast / 0 /
     .  NOrbNodeLast / 0 /

!      if (thisNode /= Node)
!     $       call message("Non-local use of GetNodeOrbs")

      if (pexsi_dist) then
         remainder = Norb - pexsi_bs*Nodes
         if (Node == Nodes-1) then
            NorbNode = pexsi_bs + remainder
         else
            NorbNode = pexsi_bs
         endif
         
         RETURN
      endif

      if (NOrb .eq. NOrbLast) then
C Values are the same as last call - no need to recalculate
        NOrbNode = NOrbNodeLast

      else

        if (lspatial) then
C-------------------------
C  Spatial distribution  -
C-------------------------
          NOrbNode = nOrbPerNode(Node+1)

        else if (use_dd) then
           if (use_dd_perm) then
              !  NOrbNode = dd_nuo
              ! We have this structure also
              NOrbNode = nOrbPerNode(Node+1)
           else
              if (thisNode /= Node)
     $             call die("Wrong use of dd GetNodeOrbs")
              NOrbNode = dd_nuo
           endif
        else
C-----------------------------
C  Blockcyclic distribution  -
C-----------------------------
C Calculate the minimum number of orbitals per node
          MinPerNode = NOrb / (Nodes*BlockSize)

C Find the remainder of unassigned orbitals
          Remainder = NOrb - MinPerNode * Nodes * BlockSize

C Find number of complete blocks in the remainder
          RemainderBlocks = Remainder / BlockSize
          Remainder = Remainder - RemainderBlocks * BlockSize

C Workout the local number of orbitals
          NOrbNode = MinPerNode*BlockSize
          if (Node.lt.RemainderBlocks) NOrbNode = NOrbNode + BlockSize
          if (Node.eq.RemainderBlocks) NOrbNode = NOrbNode + Remainder
        endif

C Save value for next call
        NOrbNodeLast = NOrbNode

      endif

      end subroutine GetNodeOrbs

      subroutine GlobalToLocalOrb( GOrb, Node, Nodes, LOrb)
      use parallel, only: thisNode => Node
C
C  Converts an orbital index in the global frame to the local frame
C  if the orbital is local to this node. Otherwise the pointer is
C  return as zero and can therefore be used to test whether the
C  orbital is local or not.
C
C  Julian Gale, Imperial College, November 1998
C
C  Input :
C
C  integer GOrb   = global orbital index
C  integer Node   = local processor number
C  integer Nodes  = global number of processors
C
C  From parallel.h :
C
C  integer BlockSize = blocking size for orbital distribution across
C                      the nodes. Choice of value affects the 
C                      performance of the Scalapack routines
C
C  Output :
C
C  integer LOrb   = local orbital index
C

      integer GOrb, Node, Nodes, LOrb

C  Local variables
      integer OrbCheck, GBlock, LEle, LBlock

!      if (thisNode /= Node)
!     $       call message("Non-local use of GlobalToLocalOrb")

      if (pexsi_dist) then
         call WhichNodeOrb(Gorb,Nodes,OrbCheck)
         if (OrbCheck == Node) then
            lorb = gorb - pexsi_bs*Node
         else
            lorb = 0
         endif
         
      else if (lspatial) then 
C-------------------------
C  Spatial distribution  -
C-------------------------
        if (Node.eq.nNode(GOrb)) then
          LOrb = nG2L(GOrb)
        else
          LOrb = 0
        endif
      else if (use_dd) then
        if (use_dd_perm) then
          ! LOrb = dd_perm(GOrb)
           LOrb = nG2L(GOrb)
        else
          if (GOrb.ge.llimit .and. GOrb.lt.ulimit) then
            LOrb = GOrb - llimit + 1
          else
            LOrb = 0
          endif
        endif
      else
C-----------------------------
C  Blockcyclic distribution  -
C-----------------------------
C  Find global block number
        GBlock = ((GOrb -1)/BlockSize)

C  Substract global base line to find element number within the block
        LEle = GOrb - GBlock*BlockSize

C  Find the block number on the local node
        LBlock = ((GBlock - Node)/Nodes)

C  Generate the local orbital pointer based on the local block number
        LOrb = LEle + LBlock*Blocksize

C  Check that this is consistent - if it is not then this
C  local orbital is not on this node and so we return 0
C  to indicate this.
        OrbCheck = (LBlock*Nodes + Node)*BlockSize + LEle
        if (OrbCheck.ne.GOrb) LOrb = 0
      endif

      end  subroutine GlobalToLocalOrb


      subroutine LocalToGlobalOrb( LOrb, Node, Nodes, GOrb)
      use parallel, only: thisNode => Node
C
C  Converts an orbital index in the local frame to the global frame
C
C  Julian Gale, Imperial College, December 1998
C
C  Input :
C
C  integer LOrb   = local orbital index
C  integer Node   = local processor number
C  integer Nodes  = global number of processors
C
C  From parallel.h :
C
C  integer BlockSize = blocking size for orbital distribution across
C                      the nodes. Choice of value affects the 
C                      performance of the Scalapack routines
C
C  Output :
C
C  integer GOrb   = global orbital index
C

      integer LOrb, Node, Nodes, GOrb

C  Local variables
      integer LEle, LBlock

!      if (thisNode /= Node)
!     $       call message("Non-local use of LocalToGlobalOrb")

      if (pexsi_dist) then
         gorb = pexsi_bs*Node + lorb

      else if (lspatial) then 
C-------------------------
C  Spatial distribution  - 
C-------------------------
        GOrb = nL2G(LOrb,Node+1)

      else if (use_dd) then
        if (use_dd_perm) then
          !  GOrb = dd_invp(LOrb)
          ! We have this structures also...
           GOrb = nL2G(LOrb,Node+1)
        else
           if (thisNode /= Node)
     $          call die("Wrong use of dd LocalToGlobalOrb")

          GOrb = LOrb + llimit - 1
        endif
      else   
C-----------------------------
C  Blockcyclic distribution  -
C-----------------------------
C  Find local block number
        LBlock = ((LOrb -1)/BlockSize)

C  Substract local base line to find element number within the block
        LEle = LOrb - LBlock*BlockSize

C  Calculate global index
        GOrb = (LBlock*Nodes + Node)*BlockSize + LEle
      endif
      return
      end  subroutine LocalToGlobalOrb

      subroutine WhichNodeOrb( GOrb, Nodes, Node)
C
C  Given the global orbital pointer, this routine
C  returns the Node number where this is stored.
C
C  Julian Gale, Imperial College, January 1999
C
C  Input :
C
C  integer GOrb   = global orbital index
C  integer Nodes  = total number of Nodes
C
C  From parallel.h :
C
C  integer BlockSize = blocking size for orbital distribution across
C                      the nodes. Choice of value affects the 
C                      performance of the Scalapack routines
C
C  Output :
C
C  integer Node   = Node where this orbital is stored locally
C
      integer GOrb, Node, Nodes

C  Local variables
      integer GBlock

      if (pexsi_dist) then
         Node = (Gorb-1)/pexsi_bs
         ! Check for overflow: if gorb is in the "remainder" set, 
         ! the above will be incorrect (e.g. 4 (or even 5) in a {0,1,2,3} set)
         ! (Consider the case norbs=13, pexsi_bs=2, nodes=5)
         if (Node > Nodes-1) Node = Nodes-1

      else if (lspatial) then 
C-------------------------
C  Spatial distribution  - 
C-------------------------
        Node = nNode(GOrb)
      else if (use_dd) then
         if (.not. use_dd_perm) then
            call die("WhichNodeOrb not ready in DD")
         endif
        Node = dd_nnode(GOrb)
        !! Node = nNode(GOrb)
      else   
C-----------------------------
C  Blockcyclic distribution  -
C-----------------------------
C  Find global block number
        GBlock = ((GOrb -1)/BlockSize)

C  Find the Node number that has this block
        Node = mod(GBlock,Nodes)
      endif

      return
      end  subroutine WhichNodeOrb


      subroutine WhichMeshNode( GMesh, Nxyz, Nodes, Node)
C
C  Given the global mesh pointer, this routine
C  returns the Node number where this is stored.
C
C  Julian Gale, Imperial College, March 1999
C
C  Input :
C
C  integer GMesh  = global mesh index
C  integer Nxyz(3)= number of grid points along each axis
C  integer Nodes  = total number of Nodes
C
C  Output :
C
C  integer Node   = Node where this orbital is stored locally
C

C  Passed variables
      integer GMesh, Nxyz(3), Node, Nodes

C  Local variables
      integer 
     .  Gy, Gz, Mesh(3)

C  Trap zero mesh case
      if (Nxyz(1)*Nxyz(2)*Nxyz(3).eq.0) then
        Node = 0
        return
      endif

C  Find grid indices along each axis
      Gz = ((GMesh-1)/(Nxyz(1)*Nxyz(2))) + 1
      Gy = GMesh - (Gz-1)*Nxyz(1)*Nxyz(2)
      Gy = ((Gy-1)/Nxyz(1))+1
      Mesh(1) = 1
      Mesh(2) = Gy
      Mesh(3) = Gz

C  Call subroutine to find node number from global XYZ index
      call WhichMeshXYZNode( Mesh, Nxyz, Nodes, Node)

      end subroutine WhichMeshNode

      subroutine HowManyMeshPerNode(Nxyz, Node, Nodes, NMeshPN, NxyzL)
C
C  Given the total dimensions of the grid of points it works
C  out how many of the points are on the current node.
C
C  Julian Gale, Imperial College, March 1999
C
C  Input :
C
C  integer Nxyz(3) = dimensions of grid of points
C  integer Node    = local processor number
C  integer Nodes   = global number of processors
C
C  From parallel.h :
C
C  integer ProcessorY = second dimension of processor grid - 
C                      must be a factor of the number of processors
C                      being used. 
C
C  Output :
C
C  integer NMeshPN  = total number of local mesh points
C  integer NxyzL(3) = local dimensions of mesh points on this node
C

      implicit none
C  Passed variables
      integer Nxyz(3), Node, Nodes, NMeshPN, NxyzL(3)

C  Local variables
      integer 
     .  BlockSizeY, BlockSizeZ, NRemY, NRemZ, Py, Pz,
     .  ProcessorZ

C  Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      ProcessorZ = Nodes/ProcessorY

C  Find processor grid location
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1

C  Set blocking sizes
      BlockSizeY = (Nxyz(2)/ProcessorY)
      BlockSizeZ = (Nxyz(3)/ProcessorZ)

      NRemY = Nxyz(2) - BlockSizeY*ProcessorY
      NRemZ = Nxyz(3) - BlockSizeZ*ProcessorZ
      if (Py-1.lt.NRemY) BlockSizeY = BlockSizeY + 1
      if (Pz-1.lt.NRemZ) BlockSizeZ = BlockSizeZ + 1

C  Trap zero grid case
      if (BlockSizeY.eq.0.or.BlockSizeZ.eq.0) then
        NMeshPN = 0
        NxyzL(1) = Nxyz(1)
        NxyzL(2) = Nxyz(2)
        NxyzL(3) = Nxyz(3)
        return
      endif

C  Assign blocksizes as local grid dimensions
      NxyzL(1) = Nxyz(1)
      NxyzL(2) = BlockSizeY
      NxyzL(3) = BlockSizeZ

C  Calculate total number of grid points by multiplying dimensions
      NMeshPN = NxyzL(1)*NxyzL(2)*NxyzL(3)

      return
      end  subroutine HowManyMeshPerNode


      subroutine GlobalToLocalMesh( GMesh, Nxyz, Node, Nodes, LMesh)
C
C  Converts an orbital index in the global frame to the local frame
C  if the orbital is local to this node. Otherwise the pointer is
C  return as zero and can therefore be used to test whether the
C  orbital is local or not.
C
C  Julian Gale, Imperial College, January 1999
C
C  Input :
C
C  integer GMesh   = global mesh point index
C  integer Nxyz(3) = dimensions of grid of points
C  integer Node    = local processor number
C  integer Nodes   = global number of processors
C
C  Output :
C
C  integer LMesh   = local mesh point index
C

C  Passed variables
      integer GMesh, Nxyz(3), Node, Nodes, LMesh

C  Local variables
      integer 
     .  Mesh(3), MeshL(3), NxyzL(3), Gx, Gy, Gz, NMeshPN

C  Trap zero mesh size case
      if (Nxyz(1)*Nxyz(2)*Nxyz(3).eq.0) then
        LMesh = GMesh
        return
      endif

C  Find grid indices along each axis
      Gz = ((GMesh-1)/(Nxyz(1)*Nxyz(2))) + 1
      Gx = GMesh - (Gz-1)*Nxyz(1)*Nxyz(2)
      Gy = ((Gx-1)/Nxyz(1))+1
      Gx = Gx - (Gy-1)*Nxyz(1)
      Mesh(1) = Gx
      Mesh(2) = Gy
      Mesh(3) = Gz

C  Call subroutine to find local mesh pointer based on global XYZ grid
      call GlobalToLocalXYZMesh( Mesh, Nxyz, Node, Nodes, MeshL)

C  If MeshL = 0 then this point is not on this node -> return
      if (MeshL(1)*MeshL(2)*MeshL(3).eq.0) then
        LMesh = 0
        return
      endif

C  Find dimensions of local mesh so that we can work out pointer
      call HowManyMeshPerNode(Nxyz, Node, Nodes, NMeshPN, NxyzL)

C  Generate the local mesh pointer based on the local element and blocks
      LMesh = (MeshL(3)-1)*(Nxyz(1)*NxyzL(2)) + 
     .        (MeshL(2)-1)*Nxyz(1) + Gx

      end  subroutine GlobalToLocalMesh


      subroutine WhichMeshXYZNode( Mesh, Nxyz, Nodes, Node)
C
C  Given the global mesh pointer, this routine
C  returns the Node number where this is stored.
C
C  Julian Gale, Imperial College, March 1999
C
C  Input :
C
C  integer Mesh(3) = global mesh point by grid reference
C  integer Nxyz(3) = number of grid points along each axis
C  integer Nodes   = total number of Nodes
C
C  From parallel.h :
C
C  integer ProcessorY = second dimension of processor grid - 
C                      must be a factor of the number of processors
C                      being used. 
C
C  Output :
C
C  integer Node   = Node where this orbital is stored locally
C

C  Passed variables
      integer Mesh(3), Nxyz(3), Node, Nodes

C  Local variables
      integer 
     .  ProcessorZ, PGy, PGz, BlockSizeY, BlockSizeZ,
     .  GSplitY, GSplitZ, NRemY, NRemZ

C  Trap zero mesh size case
      if (Nxyz(1)*Nxyz(2)*Nxyz(3).eq.0) then
        Node = 0
        return
      endif

C  Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

      ProcessorZ = Nodes/ProcessorY

C  Find blocksizes along axes
      BlockSizeY = (Nxyz(2)/ProcessorY)
      BlockSizeZ = (Nxyz(3)/ProcessorZ)
      NRemY = Nxyz(2) - BlockSizeY*ProcessorY
      NRemZ = Nxyz(3) - BlockSizeZ*ProcessorZ

C  Find grid point where the change from BlockSize+1 to BlockSize happens
      GSplitY = NRemY*(BlockSizeY + 1)
      GSplitZ = NRemZ*(BlockSizeZ + 1)

C  Find processor grid coordinates for this point
      if (Mesh(2).le.GSplitY) then
        PGy = (Mesh(2) - 1)/(BlockSizeY + 1) + 1
      else
        PGy = (Mesh(2) - 1 - GSplitY)/(BlockSizeY) + NRemY + 1
      endif
      if (Mesh(3).le.GSplitZ) then
        PGz = (Mesh(3) - 1)/(BlockSizeZ + 1) + 1
      else
        PGz = (Mesh(3) - 1 - GSplitZ)/(BlockSizeZ) + NRemZ + 1
      endif
      
C  Calculate processor number for this grid point
      Node = (PGy - 1)*ProcessorZ + PGz - 1

      return
      end  subroutine WhichMeshXYZNode


      subroutine GlobalToLocalXYZMesh( Mesh, Nxyz, Node, Nodes, LMesh)
C
C  Converts an orbital index in the global frame to the local frame
C  if the orbital is local to this node. Otherwise the pointer is
C  return as zero and can therefore be used to test whether the
C  orbital is local or not.
C
C  Julian Gale, Imperial College, January 1999
C
C  Input :
C
C  integer Mesh(3) = global mesh point grid reference
C  integer Nxyz(3) = dimensions of grid of points
C  integer Node    = local processor number
C  integer Nodes   = global number of processors
C
C  From parallel.h :
C
C  integer ProcessorY = second dimension of processor grid - 
C                      must be a factor of the number of processors
C                      being used. 
C
C  Output :
C
C  integer LMesh(3) = local mesh point grid reference
C

C  Passed variables
      integer Mesh(3), Nxyz(3), Node, Nodes, LMesh(3)

C  Local variables
      integer 
     .  ProcessorZ, Py, Pz, PGy, PGz, BlockSizeY, BlockSizeZ,
     .  GSplitY, GSplitZ, NRemY, NRemZ

C  Trap zero mesh size case
      if (Nxyz(1)*Nxyz(2)*Nxyz(3).eq.0) then
        LMesh(1) = Mesh(1)
        LMesh(2) = Mesh(2)
        LMesh(3) = Mesh(3)
        return
      endif

C  Check that ProcessorY is a factor of the number of processors
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      ProcessorZ = Nodes/ProcessorY

C  Find processor grid location
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1

C  Find blocksizes along axes
      BlockSizeY = (Nxyz(2)/ProcessorY)
      BlockSizeZ = (Nxyz(3)/ProcessorZ)
      NRemY = Nxyz(2) - BlockSizeY*ProcessorY
      NRemZ = Nxyz(3) - BlockSizeZ*ProcessorZ

C  Find grid point where the change from BlockSize+1 to BlockSize happens
      GSplitY = NRemY*(BlockSizeY + 1)
      GSplitZ = NRemZ*(BlockSizeZ + 1)

C  Find processor grid coordinates for this point
      if (Mesh(2).le.GSplitY) then
        PGy = (Mesh(2) - 1)/(BlockSizeY + 1) + 1
      else
        PGy = (Mesh(2) - 1 - GSplitY)/(BlockSizeY) + NRemY + 1
      endif
      if (Mesh(3).le.GSplitZ) then
        PGz = (Mesh(3) - 1)/(BlockSizeZ + 1) + 1
      else
        PGz = (Mesh(3) - 1 - GSplitZ)/(BlockSizeZ) + NRemZ + 1
      endif
      
C  If this node isn't the local one then set pointer to zero and return
      if (Py.ne.PGy.or.Pz.ne.PGz) then
        LMesh(1) = 0
        LMesh(2) = 0
        LMesh(3) = 0
        return
      endif

C  Get local mesh pointer by subtracting baseline
      LMesh(1) = Mesh(1)
      if (Mesh(2).le.GSplitY) then
        LMesh(2) = Mesh(2) - (PGy - 1)*(BlockSizeY + 1)
      else
        LMesh(2) = Mesh(2) - (PGy - 1)*BlockSizeY - NRemY
      endif
      if (Mesh(3).le.GSplitZ) then
        LMesh(3) = Mesh(3) - (PGz - 1)*(BlockSizeZ + 1)
      else
        LMesh(3) = Mesh(3) - (PGz - 1)*BlockSizeZ - NRemZ
      endif

      end subroutine GlobalToLocalXYZMesh

      subroutine set_processorYdefault(Nodes,procYdefault)
C
C Finds a sensible default value for the processor grid in the Y
C direction to try to get an equal split in X and Y. Currently only
C looks for factors of 2, 3 and 5.
C
C Input :
C
C integer Nodes        : total number of processors
C
C Output :
C
C integer procYdefault : default value of Y grid parameter
C
C Written by Julian Gale, November 1999
C
C Passed arguments
      integer
     .  Nodes, procYdefault
C Local variables
      integer
     .  Nx, Ny, Nrem
      logical
     .  factor

C Initialise values
      Nx = 1
      Ny = 1
      Nrem = Nodes
      factor = .true.

C Loop looking for factors
      do while (factor.and.Nrem.gt.1)
        factor = .false.
        if (mod(Nrem,2).eq.0) then
          Nrem = Nrem/2
          factor = .true.
          if (Nx.gt.Ny) then
            Ny = 2*Ny
          else
            Nx = 2*Nx
          endif
        endif
        if (mod(Nrem,3).eq.0) then
          Nrem = Nrem/3
          factor = .true.
          if (Nx.gt.Ny) then
            Ny = 3*Ny
          else
            Nx = 3*Nx
          endif
        endif
        if (mod(Nrem,5).eq.0) then
          Nrem = Nrem/5
          factor = .true.
          if (Nx.gt.Ny) then
            Ny = 5*Ny
          else
            Nx = 5*Nx
          endif
        endif
      enddo

C Choose default value as lowest of Nx and Ny
      if (Nx.gt.Ny) then
        procYdefault = Ny
      else
        procYdefault = Nx
      endif

      end subroutine set_processorYdefault

#ifdef MPI
      subroutine set_blocksizedefault(Nodes,nuotot,bs)
C
C Finds a sensible default value for the blocksize default.
C When the number of orbitals is less than the blocksize
C typically used then lower the blocksize to ensure that 
C some work is done on all nodes.
C
C Input :
C
C integer Nodes        : total number of processors
C integer nuotot       : total number of orbitals
C
C Output :
C
C integer bs : default value of blocksize
C
C Written by Julian Gale, March 2001
C Modified by Alberto Garcia, October 2013
C
C Passed arguments
      integer, intent(in)  :: Nodes, nuotot
      integer, intent(out) :: bs

      integer  :: n1, nn

      ! Scalapack routine for block-cyclic distributions
      integer, external :: numroc

C Compare number of orbitals against sensible number
! Note that this number could be optimized...

      if (nuotot.gt.24*Nodes) then
        bs = 24
      else if (nuotot.lt.Nodes) then
        bs = 1
      else

         bs = nuotot/Nodes + 1
         do
            ! Check the number of orbitals handled
            ! by the first and last processors
            ! Avoid idle processors, and too-high
            ! imbalances

            n1 = numroc(nuotot,bs,0,0,nodes)
            nn = numroc(nuotot,bs,nodes-1,0,nodes)

            if (nn == 0) then
               bs = bs - 1
            else if (dble(n1)/nn > 2) then
               bs = bs - 1
            else
               exit
            endif
         enddo
               
      endif

      end subroutine set_blocksizedefault
#endif
      end module parallelsubs
