! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine cellxc( irel, ider, cell, NMesh, NSpan, maxp, mtype, 
     .                   XMesh, nspin, Dens, EX, EC, DX, DC, VXC, 
     .                   DVXCDN, stressl )

C *******************************************************************
C Finds total exchange-correlation energy and potential in a
C   periodic cell.
C This version implements the Local (spin) Density Approximation and
C   the Generalized-Gradient-Aproximation with the 'explicit mesh 
C   functional' approach of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
C Wrtten by J.M.Soler using algorithms developed by 
C   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996 - Aug.1997
C Parallel version written by J.Gale. June 1999.
C Argument DVXCDN added by J.Junquera. November 2000.
C Adapted for multiple functionals in the same run by J.Gale 2005
C ************************* INPUT ***********************************
C integer irel         : Relativistic exchange? (0=>no, 1=>yes)
C integer ider         : Return dVxc/drho in DVXCDN?
C                        0=>no, 1=>yes (available only for LDA)
C real*8  cell(3,3)    : Unit cell vectors cell(ixyz,ivector)
C integer NMesh(3)     : Number of mesh divisions of each vector
C integer NSpan(3)     : Physical dimensions of arrays XMesh, Dens and
C                        VXC (or memory span between array elements)
C                        See usage section for more information
C integer maxp         : Physical dimension of XMesh, Dens, and VXC
C integer mtype        : Mesh type:
C                        0 => Uniform mesh
C                        1 => Adaptive mesh, given in cartesian coord
C                        2 => Adaptive mesh, given in cell-vector coord
C real    XMesh(3,maxp): Mesh point coordinates (not used if mtype=0)
C                        When mtype=2, cartesian coordinates are
C                        Xcart(ix,im) = Sum_iv(cell(ix,iv)*XMesh(iv,ip))
C                        Notice single precision in this version
C integer nspin        : nspin=1 => unpolarized; nspin=2 => polarized;
C                        nspin=4 => non-collinear polarization
C real    Dens(maxp,nspin) : Total (nspin=1) or spin (nspin=2) electron
C                        density at mesh points, ordered as
C                        IP = I1+NSpan(1)*((I2-1)+NSpan(2)*(I3-1)),
C                        with I1=1,...,NMesh(1), etc 
C                        For non-collinear polarization, the density
C                        matrix is given by: Dens(1)=D11, Dens(2)=D22,
C                        Dens(3)=Real(D12), Dens(4)=Im(D12)
C                        Notice single precision in this version
C ************************* OUTPUT **********************************
C real*8  EX              : Total exchange energy
C real*8  EC              : Total correlation energy
C real*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C real*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C real    VXC(maxp,nspin) : (Spin) exch-corr potential
C                           Notice single precision in this version
C real    DVXCDN(maxp,nspin,nspin) : Derivatives of exchange-correlation
C                           potential respect the charge density
C                           Not used unless ider=1. Available only for LDA
C real*8  stressl(3,3)    : xc contribution to the stress tensor,
C                           assuming constant density (not charge),
C                           i.e. r->r' => rho'(r') = rho(r)
C                           For plane-wave and grid (finite diff) basis
C                           sets, density rescaling gives an extra term
C                           (not included) (DX+DC-EX-EC)/cell_volume for
C                           the diagonal elements of stress. For other
C                           basis sets, the extra term is, in general:
C                           IntegralOf(v_xc * d_rho/d_strain) / cell_vol
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUnit below
C Stress in EUnit/Bohr**3
C ************************ USAGE ************************************
C Typical calls for different array dimensions:
C     parameter ( maxp = 1000000 )
C     DIMENSION NMesh(3), Dens(maxp,2), VXC(maxp,2)
C     Find cell vectors
C     Find density at N1*N2*N3 mesh points (less than maxp) and place 
C       them consecutively in array Dens
C     NMesh(1) = N1
C     NMesh(2) = N2
C     NMesh(3) = N3
C     CALL cellXC( 0, 0, cell, NMesh, NMesh, maxp, 0, XMesh,
C    .             2, Dens, EX, EC, DX, DC, VXC, DVXCDN, stress )
C Or alternatively:
C     parameter ( M1=100, M2=100, M3=100 )
C     DIMENSION NMesh(3), NSpan(3), Dens(M1,M2,M3,2), VXC(M1,M2,M3,2)
C     DATA NSpan / M1, M2, M3 /
C     Find cell vectors
C     Find Dens at N1*N2*N3 mesh points
C     NMesh(1) = N1
C     NMesh(2) = N2
C     NMesh(3) = N3
C     CALL cellXC( 0, 0, cell, NMesh, NSpan, M1*M2*M3, 0, XMesh,
C    .             2, Dens, EX, EC, DX, DC, VXC, DVXCDN, stress )
C ********* BEHAVIOUR ***********************************************
C - Notice that XMesh is not used if mtype=0, and DVXCDN is not
C   used if ider=0. This means that you do not even need to dimension 
C   them in the calling program. See usage section for calling examples.
C - If FNCTNL='LDA', Dens and VXC may be the same physical array.
C - Stops and prints a warning if arguments functl or XCauth are not
C   one of the allowed possibilities.
C - Since the exchange and correlation part is usually a small fraction
C   of a typical electronic structure calculation, this routine has
C   been coded with emphasis on simplicity and functionality, not in
C   eficiency. The parallel version was written by J.Gale.
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC, RECLAT, TIMER, VOLCEL
C *******************************************************************

      use precision,    only : dp, grid_p
      use xcmod,        only : nXCfunc, XCfunc, XCauth
      use xcmod,        only : XCweightX, XCweightC
      use sys,          only :die
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
C  Modules
      use parallel,     only : Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use mpi_siesta
      use mesh,         only : nsm, nmeshg
#endif

      implicit none

C Argument types and dimensions
      integer,       intent(in)   :: ider
      integer,       intent(in)   :: irel
      integer,       intent(in)   :: maxp
      integer,       intent(in)   :: mtype
      integer,       intent(in)   :: NMesh(3)
      integer,       intent(in)   :: NSpan(3)
      integer,       intent(in)   :: nspin
      real(dp),      intent(in)   :: cell(3,3)
      real(dp),      intent(out)  :: DC
      real(dp),      intent(out)  :: DX
      real(dp),      intent(out)  :: EC
      real(dp),      intent(out)  :: EX
      real(dp),      intent(inout)  :: stressl(3,3)

C If you change next line, please change also the argument
C explanations in the header comments
!!! Pending (AG)
*     real(dp)
      real(grid_p),          intent(in)   :: Dens(maxp,nspin)
      real(grid_p),          intent(out)  :: DVXCDN(maxp,nspin,nspin)
      real(grid_p),          intent(out)  :: VXC(maxp,nspin)
      real(grid_p),          intent(in)   :: XMesh(3,*)

C Fix the order of the numerical derivatives
C NN is the number of points used in each coordinate and direction,
C i.e. a total of 6*NN neighbour points is used to find the gradients
      integer,       parameter    :: NN = 3

C Fix energy unit:  EUnit=1.0 => Hartrees,
C                   EUnit=0.5 => Rydbergs,
C                   EUnit=0.03674903 => eV
      real(dp),      parameter    :: EUnit = 0.5_dp

C Fix switch to skip points with zero density
      logical,       parameter    :: skip0 = .true.
      real(dp),      parameter    :: denmin = 1.0e-15_dp

C Parameter mspin must be equal or larger than nspin
      integer,       parameter    :: mspin = 4

#ifdef MPI
C MPI related variables
      integer 
     .  MPIerror, ProcessorZ, BlockSizeY, BlockSizeZ, MeshNsm(3),
     .  BlockSizeYMax, BlockSizeZMax, NRemY, NRemZ,
     .  Py, Pz
      integer 
     .  NSwap, NDummy, IDen1, IDen2, IOut, INN, MeshDum(3), 
     .  Requests(4), Statuss(MPI_Status_Size,4), Tag, Pleft2, Pright2,
     .  Pleft3, Pright3, PzMod, PyMod, NMeshPN, JPNN(3,-NN:NN),
     .  nmax, nrleft, nrright, nsleft, nsright
      real(dp)
     .  Buffer(4),Buffer2(4)
      logical 
     .  NoLocalPoints

      real(grid_p), dimension(:,:,:), pointer :: bdens, bvxc
      real(grid_p), dimension(:,:),   pointer :: bbuffer, dbuffer
      integer,      dimension(:,:),   pointer :: PMesh

#endif

C Local variables and arrays
      logical           GGA, GGAfunctl
      integer           I1, I2, I3, IC, IN, IP, IS, IX,
     .                  J1, J2, J3, JN, JP(3,-NN:NN), JS, JX,
     .                  KS, NPG, nf
      real(dp)          D(mspin), DECDD(mspin), DECDGD(3,mspin),
     .                  dentot, DEXDD(mspin), DEXDGD(3,mspin),
     .                  DGDM(-NN:NN), DGIDFJ(3,3,-NN:NN),
     .                  DMDX(3,3), DVol, DXDM(3,3),
     .                  DVCDN(mspin*mspin), DVXDN(mspin*mspin),
     .                  EPSC, EPSX, F1, F2, GD(3,mspin),
     .                  VOLCEL, volume, stress(3,3)
      external          GGAXC, LDAXC, RECLAT, VOLCEL

C Start time counter (intended only for debugging and development)
      call timer( 'cellXC', 1 )

C Set GGA switch
      GGA = .false.
      do nf = 1,nXCfunc
        if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
          GGA = .true.
        else
          if ( XCfunc(nf).ne.'LDA' .and. XCfunc(nf).ne.'lda' .and.
     .         XCfunc(nf).ne.'LSD' .and. XCfunc(nf).ne.'lsd' ) then
            write(6,*) 'cellXC: Unknown functional ', XCfunc(nf)
            call die()
          endif
        endif
      enddo

C Check ider
      if (ider.ne.0 .and. GGA)
     $     call die('cellXC: ider=1 available only for LDA')

C Check value of mspin
      if (mspin.lt.nspin) then
        write(6,*) 'cellXC: parameter mspin must be at least ', nspin
        call die()
      endif

#ifdef MPI
C If GGA and the number of processors is greater than 1 then we need
C to exchange border densities so that the density derivatives can
C be calculated.
      if (GGA.and.Nodes.gt.1) then

C Work out processor grid dimensions
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        ProcessorZ = Nodes/ProcessorY

C Work out grid coordinates of local node in processor grid
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1

C Work out maximum blocksizes
        MeshNSM(1) = NMeshG(1)/NSM
        MeshNSM(2) = NMeshG(2)/NSM
        MeshNSM(3) = NMeshG(3)/NSM
        BlockSizeY = (MeshNsm(2)/ProcessorY)
        BlockSizeZ = (MeshNsm(3)/ProcessorZ)
        NRemY = MeshNsm(2) - BlockSizeY*ProcessorY
        NRemZ = MeshNsm(3) - BlockSizeZ*ProcessorZ
        if (NRemY.gt.0) then
          BlockSizeYMax = BlockSizeY + 1
        else
          BlockSizeYMax = BlockSizeY
        endif
        if (NRemZ.gt.0) then
          BlockSizeZMax = BlockSizeZ + 1
        else
          BlockSizeZMax = BlockSizeZ
        endif
        BlockSizeY = BlockSizeY*NSM
        BlockSizeZ = BlockSizeZ*NSM
        BlockSizeYMax = BlockSizeYMax*NSM
        BlockSizeZMax = BlockSizeZMax*NSM

C Check that there are enough points on each Node to cover the finite
C difference interval - if not then exit!
        if (NN.gt.BlockSizeY.or.NN.gt.BlockSizeZ) then
          if (Node.eq.0) then
            write(6,'(''  ERROR - number of fine mesh points per '',
     .        ''Node must be greater than finite difference order '')')
          endif
          call die()
        endif

C Check that this Node was some points 
        call HowManyMeshPerNode(MeshNSM, Node, Nodes, NDummy, MeshDum)
        NoLocalPoints = (NDummy.eq.0)
        if (NoLocalPoints) then
          EX = 0.0_dp
          EC = 0.0_dp
          DX = 0.0_dp
          DC = 0.0_dp
          if (GGA) then
            do IS = 1,nspin
              do I3 = 0,NMesh(3)-1
                do I2 = 0,NMesh(2)-1
                  do I1 = 0,NMesh(1)-1
                    IP = 1+I1+NSpan(1)*I2+NSpan(1)*NSpan(2)*I3
                    VXC(IP,IS) = 0.0_grid_p
                  enddo
                enddo
              enddo
            enddo
          endif
          goto 991
        endif

C Allocate memory for border densities and potential

        nullify( bdens )
        call re_alloc( bdens,
     &                 1, NSpan(1)*max(BlockSizeYMax,BlockSizeZMax),
     &                 1, 4*NN, 1, nspin, name='bdens',
     &                 routine='cellxc' )
        nullify( bvxc )
        call re_alloc( bvxc,
     &                 1,NSpan(1)*max(BlockSizeYMax,BlockSizeZMax),
     &                 1, 4*NN, 1, nspin, name='bvxc',
     &                 routine='cellxc' )
        nullify( PMesh )
        call re_alloc( PMesh, 1, 3, 1, 4, name='PMesh',
     &                 routine='cellxc' )

C Exchange density information in the Z-direction
        if (NMesh(3).ne.NMeshG(3)) then

C Find Nodes to exchange with
          PzMod = Pz - 1
          if (PzMod.lt.1) PzMod = PzMod + ProcessorZ
          Pleft3 = (Py - 1)*ProcessorZ + PzMod - 1
          PzMod = Pz + 1
          if (PzMod.gt.ProcessorZ) PzMod = PzMod - ProcessorZ
          Pright3 = (Py - 1)*ProcessorZ + PzMod - 1

C Find size of mesh at nodes to be exchanged with
          call HowManyMeshPerNode(MeshNsm,Pleft3,Nodes,NMeshPN,
     .      PMesh(1:,3))
          call HowManyMeshPerNode(MeshNsm,Pright3,Nodes,NMeshPN,
     .      PMesh(1:,4))
          do ix = 1,3
            PMesh(ix,3) = PMesh(ix,3)*NSM
            PMesh(ix,4) = PMesh(ix,4)*NSM
          enddo

C Find sizes of transfer arrays
          nrleft = PMesh(1,3)*PMesh(2,3)*nn*nspin
          nrright = PMesh(1,4)*PMesh(2,4)*nn*nspin
          nsleft = NSpan(1)*NSpan(2)*nn*nspin
          nsright = NSpan(1)*NSpan(2)*nn*nspin

C Allocate arrays for data storage
          nmax = max(nrleft,nrright)
          nullify( bbuffer )
          call re_alloc( bbuffer, 1, nmax, 1, 2, name='bbuffer',
     &                   routine='cellxc' )
          nmax = max(nsleft,nsright)
          nullify( dbuffer )
          call re_alloc( dbuffer, 1, nmax, 1, 2, name='dbuffer',
     &                   routine='cellxc' )

C Build buffer transfer data array
          nsleft = 0
          nsright = 0
          do in = 1,nn
            do is = 1,nspin

C Collect data to be sent
              NSwap = NSpan(1)*NSpan(2)
              IDen1 = NSwap*(NMesh(3)-in)
              dbuffer(nsright+1:nsright+NSwap,2) = 
     .          dens(IDen1+1:IDen1+NSwap,is)
              nsright = nsright + NSwap

              IDen1 = NSwap*(in-1)
              dbuffer(nsleft+1:nsleft+NSwap,1) = 
     .          dens(IDen1+1:IDen1+NSwap,is)
              nsleft = nsleft + NSwap

            enddo
          enddo

C Transfers
          call MPI_IRecv(bbuffer(1,1),nrleft,MPI_grid_real,
     .      Pleft3,1,MPI_Comm_World,Requests(1),MPIerror)
          call MPI_IRecv(bbuffer(1,2),nrright,MPI_grid_real,
     .      Pright3,2,MPI_Comm_World,Requests(2),MPIerror)
          call MPI_ISend(dbuffer(1,1),nsleft,MPI_grid_real,
     .      Pleft3,2,MPI_Comm_World,Requests(3),MPIerror)
          call MPI_ISend(dbuffer(1,2),nsright,MPI_grid_real,
     .      Pright3,1,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
          call MPI_WaitAll(4,Requests,Statuss,MPIerror)

C Barrier - just in case
          call MPI_Barrier(MPI_Comm_World,MPIerror)

C Place transferred data in correct arrays
          nrleft = 0
          nrright = 0
          do in = 1,nn
            do is = 1,nspin

C Sort received data
              NSwap = PMesh(1,3)*PMesh(2,3)
              bdens(1:NSwap,2*nn+in,is) = 
     .          bbuffer(nrleft+1:nrleft+NSwap,1)
              nrleft = nrleft + NSwap

              NSwap = PMesh(1,4)*PMesh(2,4)
              bdens(1:NSwap,3*nn+in,is) = 
     .          bbuffer(nrright+1:nrright+NSwap,2)
              nrright = nrright + NSwap

            enddo
          enddo

C Free arrays for data storage
          call de_alloc( dbuffer, name='dbuffer' )
          call de_alloc( bbuffer, name='bbuffer' )

        endif

C Exchange density information in the Y-direction
        if (NMesh(2).ne.NMeshG(2)) then

C Find Nodes to exchange with
          PyMod = Py - 1
          if (PyMod.lt.1) PyMod = PyMod + ProcessorY
          Pleft2 = (PyMod - 1)*ProcessorZ + Pz - 1
          PyMod = Py + 1
          if (PyMod.gt.ProcessorY) PyMod = PyMod - ProcessorY
          Pright2 = (PyMod - 1)*ProcessorZ + Pz - 1

C Find size of mesh at nodes to be exchanged with
          call HowManyMeshPerNode(MeshNsm,Pleft2,Nodes,NMeshPN,
     .      PMesh(1:,1))
          call HowManyMeshPerNode(MeshNsm,Pright2,Nodes,NMeshPN,
     .      PMesh(1:,2))
          do ix = 1,3
            PMesh(ix,1) = PMesh(ix,1)*NSM
            PMesh(ix,2) = PMesh(ix,2)*NSM
          enddo

C Find size of transfer arrays
          nrleft = PMesh(1,1)*PMesh(3,1)*nn*nspin
          nrright = PMesh(1,2)*PMesh(3,2)*nn*nspin
          nsleft = NSpan(1)*NMesh(3)*nn*nspin
          nsright = NSpan(1)*NMesh(3)*nn*nspin

C Allocate transfer arrays
          nmax = max(nrleft,nrright)
          nullify( bbuffer )
          call re_alloc( bbuffer, 1, nmax, 1, 2, name='bbuffer',
     &                   routine='cellxc' )
          nmax = max(nsleft,nsright)
          nullify( dbuffer )
          call re_alloc( dbuffer, 1, nmax, 1, 2, name='dbuffer',
     &                   routine='cellxc' )

C Fill data array for sending
          nsleft = 0
          nsright = 0
          do in = 1,nn
            do is = 1,nspin

C Post density sends
              NSwap = NSpan(1)
              IDen1 = (NSpan(2)-in)*NSpan(1)
              do ix = 1,NMesh(3)
                dbuffer(nsright+1:nsright+NSwap,2) =
     .            dens(IDen1+1:IDen1+NSwap,is)
                nsright = nsright + NSwap
                IDen1 = IDen1 + NSpan(1)*NSpan(2)
              enddo
              IDen1 = (in-1)*NSpan(1)
              do ix = 1,NMesh(3)
                dbuffer(nsleft+1:nsleft+NSwap,1) =
     .            dens(IDen1+1:IDen1+NSwap,is)
                nsleft = nsleft + NSwap
                IDen1 = IDen1 + NSpan(1)*NSpan(2)
              enddo

            enddo
          enddo

C Transfers
          call MPI_IRecv(bbuffer(1,1),nrleft,MPI_grid_real,
     .      Pleft2,1,MPI_Comm_World,Requests(1),MPIerror)
          call MPI_IRecv(bbuffer(1,2),nrright,MPI_grid_real,
     .      Pright2,2,MPI_Comm_World,Requests(2),MPIerror)
          call MPI_ISend(dbuffer(1,1),nsleft,MPI_grid_real,
     .      Pleft2,2,MPI_Comm_World,Requests(3),MPIerror)
          call MPI_ISend(dbuffer(1,2),nsright,MPI_grid_real,
     .      Pright2,1,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
          call MPI_WaitAll(4,Requests,Statuss,MPIerror)

C Barrier - just in case
          call MPI_Barrier(MPI_Comm_World,MPIerror)

C Place data in required array location
          nrleft = 0
          nrright = 0
          do in = 1,nn
            do is = 1,nspin

C Post density receives
              do ix = 1,PMesh(3,1)
                NSwap = PMesh(1,1)
                IDen1 = (ix-1)*PMesh(1,1)
                bdens(IDen1+1:IDen1+NSwap,in,is) =
     .            bbuffer(nrleft+1:nrleft+NSwap,1)
                nrleft = nrleft + NSwap
              enddo
              do ix = 1,PMesh(3,2)
                NSwap = PMesh(1,2)
                IDen1 = (ix-1)*PMesh(1,2)
                bdens(IDen1+1:IDen1+NSwap,nn+in,is) =
     .            bbuffer(nrright+1:nrright+NSwap,2)
                nrright = nrright + NSwap
              enddo

            enddo
          enddo

C Free arrays for data storage
          call de_alloc( dbuffer, name='dbuffer' )
          call de_alloc( bbuffer, name='bbuffer' )

        endif

      else

        NoLocalPoints = .false.

      endif

C Barrier to ensure that all densities are in place
  991 call MPI_Barrier(MPI_Comm_World,MPIerror)
      if (NoLocalPoints) goto 992
#endif

C Find weights of numerical derivation from Lagrange interp. formula
      do IN = -NN,NN
        F1 = 1.0_dp
        F2 = 1.0_dp
        do JN = -NN,NN
          if (JN.ne.IN .and. JN.ne.0) F1 = F1 * (0  - JN)
          if (JN.ne.IN)               F2 = F2 * (IN - JN)
        enddo
        DGDM(IN) = F1 / F2
      enddo
      DGDM(0) = 0.0_dp

C Find total number of mesh points
#ifdef MPI
      NPG = NMeshG(1) * NMeshG(2) * NMeshG(3)
#else
      NPG = NMesh(1) * NMesh(2) * NMesh(3)
#endif

C Find Jacobian matrix dx/dmesh for uniform mesh
      if (mtype.eq.0) then

C Find mesh cell volume 
        DVol = VOLCEL( cell ) / DBLE(NPG)

        if (GGA) then

C Find mesh unit vectors and reciprocal mesh vectors
          do IC = 1,3
            do IX = 1,3
#ifdef MPI
              DXDM(IX,IC) = cell(IX,IC) / NMeshG(IC)
#else
              DXDM(IX,IC) = cell(IX,IC) / NMesh(IC)
#endif
            enddo
          enddo
          call reclat( DXDM, DMDX, 0 )

C Find the weights for the derivative d(gradF(i))/d(F(j)) of
C the gradient at point i with respect to the value at point j
          do IN = -NN,NN
            do IC = 1,3
              do IX = 1,3
                DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
              enddo
            enddo
          enddo

        endif

      endif

C Initialize output
      EX = 0.0_dp
      EC = 0.0_dp
      DX = 0.0_dp
      DC = 0.0_dp
      VXC(1:maxp,1:nspin) = 0.0_grid_p
      if (ider.eq.1) then
        DVXCDN(1:maxp,1:nspin,1:nspin) = 0.0_grid_p
      endif

      if (GGA) then
#ifdef MPI
C Initialise buffer regions of Vxc
        if (Nodes.gt.1) then
          do is = 1,nspin
            if (NMesh(2).ne.NMeshG(2)) then
              NSwap = NSpan(1)*BlockSizeZMax
              do in = 1,2*nn
                do ix = 1,NSwap
                  bvxc(ix,in,is) = 0.0_grid_p
                enddo
              enddo
            endif
            if (NMesh(3).ne.NMeshG(3)) then
              NSwap = NSpan(1)*BlockSizeYMax
              do in = 1,2*nn
                do ix = 1,NSwap
                  bvxc(ix,2*nn+in,is) = 0.0_grid_p
                enddo
              enddo
            endif
          enddo
        endif
#endif
      endif
      stress(1:3,1:3) = 0.0_dp

C Loop on mesh points
      do I3 = 0,NMesh(3)-1
      do I2 = 0,NMesh(2)-1
      do I1 = 0,NMesh(1)-1

C Find mesh index of this point
        IP = 1 + I1 + NSpan(1) * I2 + NSpan(1) * NSpan(2) * I3

C Skip point if density=0
        if (skip0) then
          dentot = 0.0_dp
          do IS = 1,MIN(nspin,2)
            dentot = dentot + MAX(0.0_grid_p,Dens(IP,IS))
          enddo
          if (dentot .lt. denmin) then
            do IS = 1,nspin
              VXC(IP,IS) = 0.0_grid_p
            enddo
            goto 210
          endif
        endif

C Find mesh indexes of neighbour points
C Note : a negative index indicates a point from the buffer region
        if (GGA .or. mtype.ne.0) then

C X-direction
          do IN = -NN,NN
            J1 = MOD( I1+IN+100*NMesh(1), NMesh(1) )
            JP(1,IN) = 1+J1+NSpan(1)*I2+NSpan(1)*NSpan(2)*I3
          enddo

C Y-direction
#ifdef MPI
          if (NMesh(2).eq.NMeshG(2)) then
#endif
            do IN = -NN,NN
              J2 = MOD( I2+IN+100*NMesh(2), NMesh(2) )
              JP(2,IN) = 1+I1+NSpan(1)*J2+NSpan(1)*NSpan(2)*I3
            enddo
#ifdef MPI
          else
            do IN = -NN,NN
              J2 = I2+IN
              if (J2.lt.0) then
C Out of block to the left - negative index
                IOut = -J2
                JP(2,IN) = -(1+I1+PMesh(1,1)*I3)
                JPNN(2,IN) = J2
              elseif (J2.gt.(NMesh(2)-1)) then
C Out of block to the right - negative index
                IOut = J2 - NMesh(2) + 1
                JP(2,IN) = -(1+I1+PMesh(1,2)*I3)
                JPNN(2,IN) = IOut
              else
C In block - positive index
                JP(2,IN) = 1+I1+NSpan(1)*J2+NSpan(1)*NSpan(2)*I3
              endif
            enddo
          endif
#endif

C Z-direction
#ifdef MPI
          if (NMesh(3).eq.NMeshG(3)) then
#endif
            do IN = -NN,NN
              J3 = MOD( I3+IN+100*NMesh(3), NMesh(3) )
              JP(3,IN) = 1+I1+NSpan(1)*I2+NSpan(1)*NSpan(2)*J3
            enddo
#ifdef MPI
          else
            do IN = -NN,NN
              J3 = I3+IN
              if (J3.lt.0) then
C Out of block to the left - negative index
                IOut = -J3
                JP(3,IN) = -(1+I1+PMesh(1,3)*I2)
                JPNN(3,IN) = J3
              elseif (J3.gt.(NMesh(3)-1)) then
C Out of block to the right - negative index
                IOut = J3 - NMesh(3) + 1
                JP(3,IN) = -(1+I1+PMesh(1,4)*I2)
                JPNN(3,IN) = IOut
              else
C In block - positive index
                JP(3,IN) = 1+I1+NSpan(1)*I2+NSpan(1)*NSpan(2)*J3
              endif
            enddo
          endif
#endif

        endif

C Find Jacobian matrix dx/dmesh for adaptative mesh
        if (mtype .ne. 0) then

C Find dx/dmesh
          do IC = 1,3
            do IX = 1,3
              DXDM(IX,IC) = 0.0_dp
              do IN = -NN,NN
                if (mtype .eq. 1) then
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                          XMesh(IX,JP(IC,IN)) * DGDM(IN)
                else
                  DXDM(IX,IC) = DXDM(IX,IC) +
     .                   ( cell(IX,1) * XMesh(1,JP(IC,IN)) +
     .                     cell(IX,2) * XMesh(2,JP(IC,IN)) +
     .                     cell(IX,3) * XMesh(3,JP(IC,IN)) ) * DGDM(IN)
                endif
              enddo
            enddo
          enddo

C Find inverse of matrix dx/dmesh
          call reclat( DXDM, DMDX, 0 )

C Find differential of volume = determinant of Jacobian matrix
          DVol = VOLCEL( DXDM )

C Find the weights for the derivative d(gradF(i))/d(F(j)), of
C the gradient at point i with respect to the value at point j
          if (GGA) then
            do IN = -NN,NN
              do IC = 1,3
                do IX = 1,3
                  DGIDFJ(IX,IC,IN) = DMDX(IX,IC) * DGDM(IN)
                enddo
              enddo
            enddo
          endif

        endif

C  Find density and gradient of density at this point
        do IS = 1,nspin
          D(IS) = Dens(IP,IS)
        enddo
C Test to ensure that densities are always > 0 added to 
C avoid floating point errors in ggaxc. JDG & JMS
        do IS = 1,min(nspin,2)
          D(IS) = max(0.0_dp,D(IS))
*         D(IS) = max(denmin,D(IS))
        enddo
        if (GGA) then
#ifdef MPI
          if (Nodes.eq.1) then
#endif
            do IS = 1,nspin
              do IX = 1,3
                GD(IX,IS) = 0.0_dp
                do IN = -NN,NN
                  GD(IX,IS) = GD(IX,IS) +
     .                      DGIDFJ(IX,1,IN) * Dens(JP(1,IN),IS) +
     .                      DGIDFJ(IX,2,IN) * Dens(JP(2,IN),IS) +
     .                      DGIDFJ(IX,3,IN) * Dens(JP(3,IN),IS)
                enddo
              enddo
            enddo
#ifdef MPI
          else
            do IS = 1,nspin
              do IX = 1,3
                GD(IX,IS) = 0.0_dp
                do IN = -NN,NN
                  GD(IX,IS) = GD(IX,IS) +
     .              DGIDFJ(IX,1,IN) * Dens(JP(1,IN),IS)
                enddo
                do IN = -NN,NN
                  if (JP(2,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     .                DGIDFJ(IX,2,IN) * Dens(JP(2,IN),IS)
                  else
                    INN = JPNN(2,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,2,IN) * bdens(-JP(2,IN),-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,2,IN) * bdens(-JP(2,IN),NN+INN,IS)
                    endif
                  endif
                enddo
                do IN = -NN,NN
                  if (JP(3,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     .                DGIDFJ(IX,3,IN) * Dens(JP(3,IN),IS)
                  else
                    INN = JPNN(3,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,3,IN) * bdens(-JP(3,IN),2*NN-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     .                  DGIDFJ(IX,3,IN) * bdens(-JP(3,IN),3*NN+INN,IS)
                    endif
                  endif
                enddo
              enddo
            enddo
          endif
#endif
        endif

C Loop over all functionals
        do nf = 1,nXCfunc

C Is this a GGA?
          if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
            GGAfunctl = .true.
          else
            GGAfunctl = .false.
          endif

C Find exchange and correlation energy densities and their 
C derivatives with respect to density and density gradient
          if (GGAfunctl) then
            call ggaxc( XCauth(nf), irel, nspin, D, GD,
     .                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
          else
            call ldaxc( XCauth(nf), irel, nspin, D, EPSX, EPSC, DEXDD, 
     .                  DECDD, DVXDN, DVCDN )
          endif

C Scale return values by weight for this functional
          EPSX = XCweightX(nf)*EPSX
          EPSC = XCweightC(nf)*EPSC
          do IS = 1,nspin
            DEXDD(IS) = XCweightX(nf)*DEXDD(IS)
            DECDD(IS) = XCweightC(nf)*DECDD(IS)
          enddo
          if (GGAfunctl) then
            do IS = 1,nspin
              DEXDGD(1:3,IS) = XCweightX(nf)*DEXDGD(1:3,IS)
              DECDGD(1:3,IS) = XCweightC(nf)*DECDGD(1:3,IS)
            enddo
          else
            DVXDN(1:nspin*nspin) = XCweightX(nf)*DVXDN(1:nspin*nspin)
            DVCDN(1:nspin*nspin) = XCweightC(nf)*DVCDN(1:nspin*nspin)
          endif

C Add contributions to exchange-correlation energy and its
C derivatives with respect to density at all points
          do IS = 1,min(nspin,2)
            EX = EX + DVol * D(IS) * EPSX
            EC = EC + DVol * D(IS) * EPSC
            DX = DX + DVol * D(IS) * EPSX
            DC = DC + DVol * D(IS) * EPSC
          enddo
          do IS = 1,nspin
            DX = DX - DVol * D(IS) * DEXDD(IS)
            DC = DC - DVol * D(IS) * DECDD(IS)
            if (GGAfunctl) then
              VXC(IP,IS) = VXC(IP,IS) + DEXDD(IS) + DECDD(IS)
#ifdef MPI
              if (Nodes.eq.1) then
#endif
                do IN = -NN,NN
                  do IC = 1,3
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(IC,IN),IS) *
     .                        DEXDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                      DC = DC - DVol * Dens(JP(IC,IN),IS) *
     .                        DECDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                      VXC(JP(IC,IN),IS) = VXC(JP(IC,IN),IS) + 
     .                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,IC,IN)
                    enddo
                  enddo
                enddo
#ifdef MPI
              else
                do IN = -NN,NN

C X-direction
                  do IX = 1,3
                    DX = DX - DVol * Dens(JP(1,IN),IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,1,IN)
                    DC = DC - DVol * Dens(JP(1,IN),IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,1,IN)
                    VXC(JP(1,IN),IS) = VXC(JP(1,IN),IS) + 
     .                (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,1,IN)
                  enddo

C Y-direction
                  if (JP(2,IN).gt.0) then
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(2,IN),IS) *
     .                      DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      DC = DC - DVol * Dens(JP(2,IN),IS) *
     .                      DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      VXC(JP(2,IN),IS) = VXC(JP(2,IN),IS) + 
     .                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,2,IN)
                    enddo
                  else
                    INN = JPNN(2,IN)
                    if (INN.lt.0) then
                      do IX = 1,3
                        DX = DX - DVol * bdens(-JP(2,IN),-INN,IS) *
     .                        DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        DC = DC - DVol * bdens(-JP(2,IN),-INN,IS) *
     .                        DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        bvxc(-JP(2,IN),-INN,IS) = (DEXDGD(IX,IS) + 
     .                        DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     .                        bvxc(-JP(2,IN),-INN,IS)
                      enddo
                    else
                      do IX = 1,3
                        DX = DX - DVol * bdens(-JP(2,IN),NN+INN,IS) *
     .                        DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        DC = DC - DVol * bdens(-JP(2,IN),NN+INN,IS) *
     .                        DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        bvxc(-JP(2,IN),NN+INN,IS) = (DEXDGD(IX,IS) +
     .                      DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     .                  bvxc(-JP(2,IN),NN+INN,IS)
                      enddo
                    endif
                  endif

C Z-direction
                  if (JP(3,IN).gt.0) then
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(3,IN),IS) *
     .                  DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      DC = DC - DVol * Dens(JP(3,IN),IS) *
     .                  DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      VXC(JP(3,IN),IS) = VXC(JP(3,IN),IS) + 
     .                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,3,IN)
                    enddo
                  else
                    INN = JPNN(3,IN)
                    if (INN.lt.0) then
                      do IX = 1,3
                        DX = DX - DVol * bdens(-JP(3,IN),2*NN-INN,IS) *
     .                    DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        DC = DC - DVol * bdens(-JP(3,IN),2*NN-INN,IS) *
     .                    DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        bvxc(-JP(3,IN),2*NN-INN,IS) = (DEXDGD(IX,IS) + 
     .                    DECDGD(IX,IS)) * DGIDFJ(IX,3,IN) + 
     .                    bvxc(-JP(3,IN),2*NN-INN,IS)
                      enddo
                    else
                      do IX = 1,3
                        DX = DX - DVol * bdens(-JP(3,IN),3*NN+INN,IS) *
     .                    DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        DC = DC - DVol * bdens(-JP(3,IN),3*NN+INN,IS) *
     .                    DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        bvxc(-JP(3,IN),3*NN+INN,IS) = (DEXDGD(IX,IS) + 
     .                    DECDGD(IX,IS))*DGIDFJ(IX,3,IN) +
     .                    bvxc(-JP(3,IN),3*NN+INN,IS)
                      enddo
                    endif
                  endif

                enddo
              endif

#endif
            else
              VXC(IP,IS) = VXC(IP,IS) + DEXDD(IS) + DECDD(IS)
              if (ider .eq. 1) then
                do JS = 1, nspin
                  KS = JS + (IS-1)*nspin
                  DVXCDN(IP,JS,IS) = DVXCDN(IP,JS,IS) + 
     .              DVXDN(KS) + DVCDN(KS)
                enddo
              endif
            endif
          enddo

C Add contribution to stress due to change in gradient of density
C originated by the deformation of the mesh with strain
          if (GGAfunctl) then
            do JX = 1,3
              do IX = 1,3
                do IS = 1,nspin
                  stress(IX,JX) = stress(IX,JX) - DVol * GD(IX,IS) *
     .                             ( DEXDGD(JX,IS) + DECDGD(JX,IS) )
                enddo
              enddo
            enddo
          endif

C End of loop over functionals
        enddo

  210   continue

      enddo
      enddo
      enddo


#ifdef MPI
C Return buffer region contributions to VXC to their correct nodes
  992 if (GGA.and.Nodes.gt.1) then
        if (NoLocalPoints) goto 993

C Exchange density information in the Z-direction
        if (NMesh(3).ne.NMeshG(3)) then

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post VXC receive
              NSwap = NSpan(1)*NSpan(2)
              Tag = 8*in + 2*is
              call MPI_IRecv(bdens(1,2*nn+in,is),NSwap,MPI_grid_real,
     .             Pleft3,Tag+1,MPI_Comm_World,Requests(1),MPIerror)
              call MPI_IRecv(bdens(1,3*nn+in,is),NSwap,MPI_grid_real,
     .             Pright3,Tag+2,MPI_Comm_World,Requests(2),MPIerror)

C Post VXC send
              NSwap = PMesh(1,4)*PMesh(2,4)
              call MPI_ISend(bvxc(1,3*nn+in,is),NSwap,MPI_grid_real,
     .          Pright3,Tag+1,MPI_Comm_World,Requests(3),MPIerror)
              NSwap = PMesh(1,3)*PMesh(2,3)
              call MPI_ISend(bvxc(1,2*nn+in,is),NSwap,MPI_grid_real,
     .          Pleft3,Tag+2,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
              call MPI_WaitAll(4,Requests,Statuss,MPIerror)
              call MPI_Barrier(MPI_Comm_World,MPIerror)

            enddo

          enddo

        endif

C Exchange density information in the Y-direction
        if (NMesh(2).ne.NMeshG(2)) then

C Loop over points of buffer region
          do in = 1,nn

C Loop over spins
            do is = 1,nspin

C Post VXC receives
              NSwap = NSpan(1)*NSpan(3)
              Tag = 8*in + 2*is
              call MPI_IRecv(bdens(1,in,is),NSwap,MPI_grid_real,
     .          Pleft2,Tag+1,MPI_Comm_World,Requests(1),MPIerror)
              call MPI_IRecv(bdens(1,nn+in,is),NSwap,MPI_grid_real,
     .          Pright2,Tag+2,MPI_Comm_World,Requests(2),MPIerror)

C Post VXC sends
              NSwap = PMesh(1,2)*PMesh(3,2)
              call MPI_ISend(bvxc(1,nn+in,is),NSwap,MPI_grid_real,
     .          Pright2,Tag+1,MPI_Comm_World,Requests(3),MPIerror)
              NSwap = PMesh(1,1)*PMesh(3,1)
              call MPI_ISend(bvxc(1,in,is),NSwap,MPI_grid_real,
     .          Pleft2,Tag+2,MPI_Comm_World,Requests(4),MPIerror)

C Wait for transfers to complete
              call MPI_WaitAll(4,Requests,Statuss,MPIerror)
              call MPI_Barrier(MPI_Comm_World,MPIerror)

            enddo

          enddo

        endif

C Wait for end of transfers
  993   call MPI_Barrier(MPI_Comm_World,MPIerror)
        if (NoLocalPoints) goto 999

C Return buffer VXC values to main array
        if (NMesh(2).ne.NMeshG(2)) then
          do is = 1,nspin
            do in = 1,nn
              IDen1 = (in-1)*NSpan(1) 
              IDen2 = (NMesh(2)-in)*NSpan(1) 
              do ix = 1,NMesh(3)
                do ic = 1,NMesh(1)
                  VXC(IDen1+ic,is) = VXC(IDen1+ic,is) +
     .              bdens((ix-1)*NSpan(1)+ic,in,is)
                  VXC(IDen2+ic,is) = VXC(IDen2+ic,is) +
     .              bdens((ix-1)*NSpan(1)+ic,nn+in,is)
                enddo
                IDen1 = IDen1 + NSpan(1)*NSpan(2)
                IDen2 = IDen2 + NSpan(1)*NSpan(2)
              enddo
            enddo
          enddo
        endif

        if (NMesh(3).ne.NMeshG(3)) then
          do is = 1,nspin
            do in = 1,nn
              IDen1 = NSpan(1)*NSpan(2)*(in-1)
              IDen2 = NSpan(1)*NSpan(2)*(NMesh(3)-in)
              do ix = 1,NSpan(1)*NMesh(2)
                VXC(IDen1+ix,is) = VXC(IDen1+ix,is) +
     .            bdens(ix,2*nn+in,is)
                VXC(IDen2+ix,is) = VXC(IDen2+ix,is) +
     .            bdens(ix,3*nn+in,is)
              enddo
            enddo
          enddo
        endif

      endif
#endif

C Add contribution to stress from the change of volume with strain and
C divide by volume to get correct stress definition (dE/dStrain)/Vol
      volume = VOLCEL( cell )
      do JX = 1,3
        stress(JX,JX) = stress(JX,JX) + EX + EC
        do IX = 1,3
          stress(IX,JX) = stress(IX,JX) / volume
        enddo
      enddo
      
C Divide by energy unit
      EX = EX / EUnit
      EC = EC / EUnit
      DX = DX / EUnit
      DC = DC / EUnit
      do IS = 1,nspin
        do I3 = 0,NMesh(3)-1
        do I2 = 0,NMesh(2)-1
        do I1 = 0,NMesh(1)-1
          IP = 1 + I1 + NSpan(1) * I2 + NSpan(1) * NSpan(2) * I3
          VXC(IP,IS) = VXC(IP,IS) / EUnit
        enddo
        enddo
        enddo
      enddo
      do JX = 1,3
        do IX = 1,3
          stressl(IX,JX) = stressl(IX,JX) + (stress(IX,JX) / EUnit)
        enddo
      enddo

      if (ider.eq.1 .and. .not.GGA) then
        do IS = 1,nspin
        do JS = 1,nspin
          do I3 = 0,NMesh(3)-1
          do I2 = 0,NMesh(2)-1
          do I1 = 0,NMesh(1)-1
            IP = 1 + I1 + NSpan(1) * I2 + NSpan(1) * NSpan(2) * I3
            DVXCDN(IP,JS,IS) = DVXCDN(IP,JS,IS) / EUnit
          enddo
          enddo
          enddo
        enddo
        enddo
      endif

C Re-entry point for nodes with no mesh points stored locally
  999 continue

#ifdef MPI
C Collect together quantities
      Buffer(1) = EX
      Buffer(2) = EC
      Buffer(3) = DX
      Buffer(4) = DC
      call MPI_AllReduce(Buffer,Buffer2,4,MPI_double_precision,MPI_sum,
     .  MPI_Comm_World,MPIerror)
      EX = Buffer2(1)
      EC = Buffer2(2)
      DX = Buffer2(3)
      DC = Buffer2(4)

C Free memory
      if (GGA.and.Nodes.gt.1.and..not.NoLocalPoints) then
        call de_alloc( bdens, name='bdens' )
        call de_alloc( bvxc, name='bvxc' )
        call de_alloc( PMesh, name='PMesh' )
      endif
#endif

C Stop time counter
      call timer( 'cellXC', 2 )

      end
