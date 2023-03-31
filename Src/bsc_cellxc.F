! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine bsc_cellxc( irel, ider, cell,
     $                   NMesh, NSpan, maxp, mtype, 
     &                   XMesh, nspin, Dens, EX, EC, DX, DC, VXC, 
     &                   DVXCDN, stressl )

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
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      use mesh,         only : nsm, meshLim
      use siestaxc,     only : ggaxc, ldaxc
      use siestaxc,     only : getXC
#ifdef MPI
C  Modules
      use mesh,         only : nmeshg
      use parallel,     only : Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use mpi_siesta
      use moreMeshSubs, only : distExtMeshData, gathExtMeshData
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
      integer  :: IOut, INN, JPNN(3,-NN:NN), MPIerror
#endif

C Local variables and arrays
      logical           GGA, GGAfunctl
      integer           I1, I2, I3, IC, IN, IP, IS, IX,
     &                  J1, J2, J3, JN, JP(3,-NN:NN), JS, JX,
     &                  KS, NPG, nf
      real(dp)          D(mspin), DECDD(mspin), DECDGD(3,mspin),
     &                  dentot, DEXDD(mspin), DEXDGD(3,mspin),
     &                  DGDM(-NN:NN), DGIDFJ(3,3,-NN:NN),
     &                  DMDX(3,3), DVol, DXDM(3,3),
     &                  DVCDN(mspin*mspin), DVXDN(mspin*mspin),
     &                  EPSC, EPSX, F1, F2, GD(3,mspin),
     &                  VOLCEL, volume, stress(3,3)
      external          RECLAT, VOLCEL
      
      integer           BS(3), iDistr
      real(grid_p), pointer :: bdensX(:,:,:), bdensY(:,:,:),
     &                         bdensZ(:,:,:), bvxcX(:,:,:),
     &                         bvxcY(:,:,:), bvxcZ(:,:,:)

      integer                 :: nXCfunc
      character(len=20)       :: XCauth(10), XCfunc(10)
      real(dp)                :: XCweightX(10), XCweightC(10)

#ifdef DEBUG
      call write_debug( '    PRE cellxc' )
#endif

C Start time counter (intended only for debugging and development)
#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 3 )
#endif
      call timer( 'cellXC', 1 )

      call getXC( nXCfunc, XCfunc, XCauth, XCweightX, XCweightC)

      GGA = .false.
      do nf = 1,nXCfunc
         if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga' ) then
            GGA = .true.
         endif
         if ( XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
            call die("BSC's cellxc cannot handle VDW functionals")
         endif
      enddo
      
      if (ider.ne.0 .and. GGA) then
         call die('BSC cellXC: ider=1 available only for LDA')
      endif

C Check value of mspin
      if (mspin.lt.nspin) then
        write(6,*) 'cellXC: parameter mspin must be at least ', nspin
        call die()
      endif

      BS(1) = (meshLim(2,1)-meshLim(1,1)+1)*NSM
      BS(2) = (meshLim(2,2)-meshLim(1,2)+1)*NSM
      BS(3) = (meshLim(2,3)-meshLim(1,3)+1)*NSM
#ifdef MPI
C If GGA and the number of processors is greater than 1 then we need
C to exchange border densities so that the density derivatives can
C be calculated.
      if (GGA.and.Nodes.gt.1) then
        if (NN.gt.BS(1).or.NN.gt.BS(2).or.NN.gt.BS(3)) then
          if (Node.eq.0) then
            write(6,'(''  ERROR - number of fine mesh points per '',
     &      ''Node must be greater than finite difference order '')')
          endif
          call die()
        endif

        iDistr = 3
        if (BS(1).ne.NMeshG(1)) then
          nullify(bdensX)
          call re_alloc( bdensX, 1, BS(2)*BS(3), 1, 2*NN,
     &                   1, nspin, 'bdensX', 'cellxc' )
          nullify(bvxcX)
          call re_alloc( bvxcX, 1, BS(2)*BS(3), 1, 2*NN,
     &                   1, nspin, 'bvxcX', 'cellxc' )
          call distExtMeshData( iDistr, 1, BS(2)*BS(3), NSM, NN, NSPIN,
     &                          maxp, NMeshG, dens, bdensX )
        endif

        if (BS(2).ne.NMeshG(2)) then
          nullify(bdensY)
          call re_alloc( bdensY, 1, BS(1)*BS(3), 1, 2*NN,
     &                   1, nspin, 'bdensY', 'cellxc' )
          nullify(bvxcY)
          call re_alloc( bvxcY, 1, BS(1)*BS(3), 1, 2*NN,
     &                   1, nspin, 'bvxcY', 'cellxc' )
          call distExtMeshData( iDistr, 2, BS(1)*BS(3), NSM, NN, NSPIN,
     &                          maxp, NMeshG, dens, bdensY )
        endif

        if (BS(3).ne.NMeshG(3)) then
          nullify(bdensZ)
          call re_alloc( bdensZ, 1, BS(1)*BS(2), 1, 2*NN,
     &                   1, nspin, 'bdensZ', 'cellxc' )
          nullify(bvxcZ)
          call re_alloc( bvxcZ, 1, BS(1)*BS(2), 1, 2*NN,
     &                   1, nspin, 'bvxcZ', 'cellxc' )
          call distExtMeshData( iDistr, 3, BS(1)*BS(2), NSM, NN, NSPIN,
     &                          maxp, NMeshG, dens, bdensZ )
        endif
      endif
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
          call RECLAT( DXDM, DMDX, 0 )

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

#ifdef MPI
C Initialise buffer regions of Vxc
      if (GGA.and.Nodes.gt.1) then
        if (BS(1).ne.NMeshG(1)) bvxcX = 0.0_grid_p
        if (BS(2).ne.NMeshG(2)) bvxcY = 0.0_grid_p
        if (BS(3).ne.NMeshG(3)) bvxcZ = 0.0_grid_p
      endif
#endif
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
#ifdef MPI
          if (NMesh(1).eq.NMeshG(1)) then
#endif
            do IN = -NN,NN
              J1 = MOD( I1+IN+100*BS(1), BS(1) )
              JP(1,IN) = 1+J1+BS(1)*I2+BS(1)*BS(2)*I3
            enddo
#ifdef MPI
          else
            do IN = -NN,NN
              J1 = I1+IN
              if (J1.lt.0) then
C Out of block to the left - negative index
                IOut = -J1
                JP(1,IN) = -(1+I2+BS(2)*I3)
                JPNN(1,IN) = J1
              elseif (J1.gt.(BS(1)-1)) then
C Out of block to the right - negative index
                IOut = J1 - BS(1) + 1
                JP(1,IN) = -(1+I2+BS(2)*I3)
                JPNN(1,IN) = IOut
              else
C In block - positive index
                JP(1,IN) = 1+J1+BS(1)*I2+BS(1)*BS(2)*I3
              endif
            enddo
          endif
#endif


C Y-direction
#ifdef MPI
          if (NMesh(2).eq.NMeshG(2)) then
#endif
            do IN = -NN,NN
              J2 = MOD( I2+IN+100*BS(2), BS(2) )
              JP(2,IN) = 1+I1+BS(1)*J2+BS(1)*BS(2)*I3
            enddo
#ifdef MPI
          else
            do IN = -NN,NN
              J2 = I2+IN
              if (J2.lt.0) then
C Out of block to the left - negative index
                IOut = -J2
                JP(2,IN) = -(1+I1+BS(1)*I3)
                JPNN(2,IN) = J2
              elseif (J2.gt.(BS(2)-1)) then
C Out of block to the right - negative index
                IOut = J2 - BS(2) + 1
                JP(2,IN) = -(1+I1+BS(1)*I3)
                JPNN(2,IN) = IOut
              else
C In block - positive index
                JP(2,IN) = 1+I1+BS(1)*J2+BS(1)*BS(2)*I3
              endif
            enddo
          endif
#endif

C Z-direction
#ifdef MPI
          if (NMesh(3).eq.NMeshG(3)) then
#endif
            do IN = -NN,NN
              J3 = MOD( I3+IN+100*BS(3), BS(3) )
              JP(3,IN) = 1+I1+BS(1)*I2+BS(1)*BS(2)*J3
            enddo
#ifdef MPI
          else
            do IN = -NN,NN
              J3 = I3+IN
              if (J3.lt.0) then
C Out of block to the left - negative index
                IOut = -J3
                JP(3,IN) = -(1+I1+BS(1)*I2)
                JPNN(3,IN) = J3
              elseif (J3.gt.(BS(3)-1)) then
C Out of block to the right - negative index
                IOut = J3 - BS(3) + 1
                JP(3,IN) = -(1+I1+BS(1)*I2)
                JPNN(3,IN) = IOut
              else
C In block - positive index
                JP(3,IN) = 1+I1+BS(1)*I2+BS(1)*BS(2)*J3
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
     &                          XMesh(IX,JP(IC,IN)) * DGDM(IN)
                else
                  DXDM(IX,IC) = DXDM(IX,IC) +
     &                   ( cell(IX,1) * XMesh(1,JP(IC,IN)) +
     &                     cell(IX,2) * XMesh(2,JP(IC,IN)) +
     &                     cell(IX,3) * XMesh(3,JP(IC,IN)) ) * DGDM(IN)
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
     &                      DGIDFJ(IX,1,IN) * Dens(JP(1,IN),IS) +
     &                      DGIDFJ(IX,2,IN) * Dens(JP(2,IN),IS) +
     &                      DGIDFJ(IX,3,IN) * Dens(JP(3,IN),IS)
                enddo
              enddo
            enddo
#ifdef MPI
          else
            do IS = 1,nspin
              do IX = 1,3
                GD(IX,IS) = 0.0_dp
                do IN = -NN,NN

                  if (JP(1,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     &                DGIDFJ(IX,1,IN) * Dens(JP(1,IN),IS)
                  else
                    INN = JPNN(1,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,1,IN) * bdensX(-JP(1,IN),-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,1,IN) * bdensX(-JP(1,IN),NN+INN,IS)
                    endif
                  endif

                enddo
                do IN = -NN,NN
                  if (JP(2,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     &                DGIDFJ(IX,2,IN) * Dens(JP(2,IN),IS)
                  else
                    INN = JPNN(2,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,2,IN) * bdensY(-JP(2,IN),-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,2,IN) * bdensY(-JP(2,IN),NN+INN,IS)
                    endif
                  endif
                enddo
                do IN = -NN,NN
                  if (JP(3,IN).gt.0) then
                    GD(IX,IS) = GD(IX,IS) +
     &                DGIDFJ(IX,3,IN) * Dens(JP(3,IN),IS)
                  else
                    INN = JPNN(3,IN)
                    if (INN.lt.0) then
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,3,IN) * bdensZ(-JP(3,IN),-INN,IS)
                    else
                      GD(IX,IS) = GD(IX,IS) +
     &                  DGIDFJ(IX,3,IN) * bdensZ(-JP(3,IN),NN+INN,IS)
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
     &                  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
          else
            call ldaxc( XCauth(nf), irel, nspin, D, EPSX, EPSC, DEXDD, 
     &                  DECDD, DVXDN, DVCDN )
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
     &                        DEXDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                      DC = DC - DVol * Dens(JP(IC,IN),IS) *
     &                        DECDGD(IX,IS) * DGIDFJ(IX,IC,IN)
                      VXC(JP(IC,IN),IS) = VXC(JP(IC,IN),IS) + 
     &                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,IC,IN)
                    enddo
                  enddo
                enddo
#ifdef MPI
              else
                do IN = -NN,NN

C X-direction
                  if (JP(1,IN).gt.0) then
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(1,IN),IS) *
     &                      DEXDGD(IX,IS) * DGIDFJ(IX,1,IN)
                      DC = DC - DVol * Dens(JP(1,IN),IS) *
     &                      DECDGD(IX,IS) * DGIDFJ(IX,1,IN)
                      VXC(JP(1,IN),IS) = VXC(JP(1,IN),IS) + 
     &                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,1,IN)
                    enddo
                  else
                    INN = JPNN(1,IN)
                    if (INN.lt.0) then
                      do IX = 1,3
                        DX = DX - DVol * bdensX(-JP(1,IN),-INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,1,IN)
                        DC = DC - DVol * bdensX(-JP(1,IN),-INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,1,IN)
                        bvxcX(-JP(1,IN),-INN,IS) = (DEXDGD(IX,IS) + 
     &                       DECDGD(IX,IS))*DGIDFJ(IX,1,IN) +
     &                       bvxcX(-JP(1,IN),-INN,IS)
                      enddo
                    else
                      do IX = 1,3
                        DX = DX - DVol * bdensX(-JP(1,IN),NN+INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,1,IN)
                        DC = DC - DVol * bdensX(-JP(1,IN),NN+INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,1,IN)
                        bvxcX(-JP(1,IN),NN+INN,IS) = (DEXDGD(IX,IS) +
     &                       DECDGD(IX,IS))*DGIDFJ(IX,1,IN) +
     &                       bvxcX(-JP(1,IN),NN+INN,IS)
                      enddo
                    endif
                  endif
C Y-direction
                  if (JP(2,IN).gt.0) then
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(2,IN),IS) *
     &                      DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      DC = DC - DVol * Dens(JP(2,IN),IS) *
     &                      DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                      VXC(JP(2,IN),IS) = VXC(JP(2,IN),IS) + 
     &                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,2,IN)
                    enddo
                  else
                    INN = JPNN(2,IN)
                    if (INN.lt.0) then
                      do IX = 1,3
                        DX = DX - DVol * bdensY(-JP(2,IN),-INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        DC = DC - DVol * bdensY(-JP(2,IN),-INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        bvxcY(-JP(2,IN),-INN,IS) = (DEXDGD(IX,IS) + 
     &                       DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     &                       bvxcY(-JP(2,IN),-INN,IS)
                      enddo
                    else
                      do IX = 1,3
                        DX = DX - DVol * bdensY(-JP(2,IN),NN+INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        DC = DC - DVol * bdensY(-JP(2,IN),NN+INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,2,IN)
                        bvxcY(-JP(2,IN),NN+INN,IS) = (DEXDGD(IX,IS) +
     &                       DECDGD(IX,IS))*DGIDFJ(IX,2,IN) +
     &                       bvxcY(-JP(2,IN),NN+INN,IS)
                      enddo
                    endif
                  endif

C Z-direction
                  if (JP(3,IN).gt.0) then
                    do IX = 1,3
                      DX = DX - DVol * Dens(JP(3,IN),IS) *
     &                  DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      DC = DC - DVol * Dens(JP(3,IN),IS) *
     &                  DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                      VXC(JP(3,IN),IS) = VXC(JP(3,IN),IS) + 
     &                  (DEXDGD(IX,IS)+DECDGD(IX,IS))*DGIDFJ(IX,3,IN)
                    enddo
                  else
                    INN = JPNN(3,IN)
                    if (INN.lt.0) then
                      do IX = 1,3
                        DX = DX - DVol * bdensZ(-JP(3,IN),-INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        DC = DC - DVol * bdensZ(-JP(3,IN),-INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        bvxcZ(-JP(3,IN),-INN,IS) = (DEXDGD(IX,IS) + 
     &                       DECDGD(IX,IS)) * DGIDFJ(IX,3,IN) + 
     &                       bvxcZ(-JP(3,IN),-INN,IS)
                      enddo
                    else
                      do IX = 1,3
                        DX = DX - DVol * bdensZ(-JP(3,IN),NN+INN,IS) *
     &                       DEXDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        DC = DC - DVol * bdensZ(-JP(3,IN),NN+INN,IS) *
     &                       DECDGD(IX,IS) * DGIDFJ(IX,3,IN)
                        bvxcZ(-JP(3,IN),NN+INN,IS) = (DEXDGD(IX,IS) + 
     &                       DECDGD(IX,IS))*DGIDFJ(IX,3,IN) +
     &                       bvxcZ(-JP(3,IN),NN+INN,IS)
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
     &              DVXDN(KS) + DVCDN(KS)
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
     &                             ( DEXDGD(JX,IS) + DECDGD(JX,IS) )
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
      if (GGA.and.Nodes.gt.1) then
        if (BS(1).ne.NMeshG(1)) then
          call gathExtMeshData( iDistr, 1, BS(2)*BS(3), NSM, NN, NSPIN,
     &                          maxp, NMeshG, bvxcX, VXC )
        endif
        if (BS(2).ne.NMeshG(2)) then
          call gathExtMeshData( iDistr, 2, BS(1)*BS(3), NSM, NN, NSPIN,
     &                          maxp, NMeshG, bvxcY, VXC )
        endif
        if (BS(3).ne.NMeshG(3)) then
          call gathExtMeshData( iDistr, 3, BS(1)*BS(2), NSM, NN, NSPIN,
     &                          maxp, NMeshG, bvxcZ, VXC )
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
          IP = 1 + I1 + BS(1) * I2 + BS(1) * BS(2) * I3
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
            IP = 1 + I1 + BS(1) * I2 + BS(1) * BS(2) * I3
            DVXCDN(IP,JS,IS) = DVXCDN(IP,JS,IS) / EUnit
          enddo
          enddo
          enddo
        enddo
        enddo
      endif

#ifdef MPI
C Free memory
      if (GGA.and.Nodes.gt.1) then
        if (BS(1).ne.NMeshG(1)) then
          call de_alloc( bdensX, 'bdensX', 'cellxc' )
          call de_alloc( bvxcX, 'bvxcX', 'cellxc' )
        endif
        if (BS(2).ne.NMeshG(2)) then
          call de_alloc( bdensY, 'bdensY', 'cellxc' )
          call de_alloc( bvxcY, 'bvxcY', 'cellxc' )
        endif
        if (BS(3).ne.NMeshG(3)) then
          call de_alloc( bvxcZ, 'bvxcZ', 'cellxc' )
          call de_alloc( bdensZ, 'bdensZ', 'cellxc' )
        endif
      endif
#endif

C Stop time counter
#ifdef _TRACE_
      call MPI_Barrier( MPI_Comm_World, MPIerror )
      call MPItrace_event( 1000, 0 )
#endif
      call timer( 'cellXC', 2 )

      end subroutine bsc_cellxc
