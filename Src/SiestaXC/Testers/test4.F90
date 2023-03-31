PROGRAM siestaXCtest4

  ! Compares the energy and potential calculated by atomXC and cellXC.
  ! J.M.Soler. Sept.2009

  ! Used module procedures
  USE siestaXC, only: atomXC
  USE siestaXC, only: cellXC
  USE siestaXC, only: setXC

  ! Used module parameters
  USE siestaXC, only: dp
  USE siestaXC, only: gp => grid_p

! Used MPI types
#ifdef MPI
  USE mpi_siesta, only: MPI_Double_Precision
  USE mpi_siesta, only: MPI_Max
  USE mpi_siesta, only: MPI_Sum
  USE mpi_siesta, only: MPI_Comm_World
#endif

  implicit none

  ! Tester parameters
  integer, parameter:: irel  =  0 ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nSpin =  2 ! Number of spin components
  integer, parameter:: nfTot = 19 ! Number of functionals
  integer, parameter:: nr = 501   ! Number of radial points
  integer, parameter:: nx = 60    ! Number of grid points per lattice vector
  integer, parameter:: n1cut = 8  ! Cutoff parameter
  integer, parameter:: n2cut = 2  ! Cutoff parameter:
                                  !    fCut(r)=(1-(r/rMax)**n1cut)**n2cut
  real(dp),parameter:: dWidth = 2._dp ! Width of density distribution, in Bohr
  real(dp),parameter:: Qtot = 10._dp  ! Integral of density distribution
  real(dp),parameter:: spinPol= 2._dp ! Integral of densUp - densDown
  real(dp),parameter:: rMax = 12._dp  ! Cutoff radius, in Bohr
  real(dp),parameter:: rBuff = 3._dp  ! Radial buffer of zero density, in Bohr
  real(dp),parameter:: densMin  = 1.e-9_dp  ! Min. density to proceed

  ! List of functionals to be tested
!  integer, parameter:: nf = nfTot-4   ! Number of tested functionals
!  integer:: indexf(nf) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,  18/) 
!                         ! Indexes from list below (only one VDW allowed)

  ! Same to test a single functional
  integer, parameter:: nf = 1        ! Number of tested functionals
  integer:: indexf(nf) = (/18/)      ! Indexes from list below

  ! All functionals available
  !                  1,           2,          3,           4,   
  !                  5,           6,          7,           8, 
  !                  9,          10,         11,          12,
  !                 13,          14,         15,          16,
  !                 17,          18,         19
  character(len=3):: &
    func(nfTot)=(/'LDA',       'LDA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'GGA',       'GGA', &
                  'GGA',       'GGA',       'VDW',       'VDW', &
                  'VDW',       'VDW',       'VDW'       /)
  character(len=10):: &
    auth(nfTot)=(/'PZ        ','PW92      ','PW91      ','PBE       ', &
                  'RPBE      ','revPBE    ','LYP       ','WC        ', &
                  'PBEJsJrLO ','PBEJsJrHEG','PBEGcGxLO ','PBEGcGxHEG', &
                  'PBESOL    ','AM05      ','DRSLL     ','LMKLL     ', &
                  'C09       ','BH        ','VV        '/) 

  ! Tester variables and arrays
  integer :: cellMesh(3) = (/nx,nx,nx/)
  integer :: i1, i1max, i2, i2max, i3, i3max, ir, irmax, &
             lb1, lb2, lb3, mSpin, myNode, nNodes, one, two, ub1, ub2, ub3
  real(dp):: atomDens(nr,nSpin), atomEc, atomEx, atomDc, atomDx, &
             atomVxc(nr,nSpin), avgDiffVxc, &
             cell(3,3), cellEc, cellEx, cellDc, cellDx, &
             d0, d0s(nSpin), diffVxc(nSpin), dr, dx, Ecut, kCut, &
             latConst, maxDiffVxc, pi, r, recCell(3,3), rMesh(nr), &
             stress(3,3), sumDiffVxc, tmp, Vxc(nSpin), &
             wc(nfTot), wr, wx(nfTot), x(3), x0(3)
  real(gp),allocatable:: cellDens(:,:,:,:), cellVxc(:,:,:,:)

#ifdef MPI
  ! MPI-related variables
  integer :: MPIerror
  integer :: nLarger, nxNode
#endif

#ifdef MPI
  ! Initialize MPI and get myNode and nNodes
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
#else
  myNode = 0
  nNodes = 1
#endif

  ! Initialize hybrid XC functional with all tested functionals
  mSpin = min( nSpin, 2 )
  wx = 1._dp / nf
  wc = 1._dp / nf
  call setXC( nf, func(indexf), auth(indexf), wx(indexf), wc(indexf) )

  ! Find radial mesh points and gaussian density
  pi = acos(-1._dp)
  d0 = Qtot / (2*pi*dWidth**2)**1.5_dp    ! Total density at origin
  if (nSpin==1) then
    d0s(1) = d0
  else
    one = 1   ! A silly thing to satisfy the compiler when nSpin=1
    two = 2
    d0s(one) = d0 * (Qtot + spinPol) / Qtot / 2 ! Spin up density at origin
    d0s(two) = d0 * (Qtot - spinPol) / Qtot / 2 ! Spin down density at origin
  end if
  dr = rmax / (nr-1)                      ! Interval between radial points
  do ir = 1,nr
    rMesh(ir) = dr * (ir-1)               ! Radial point values
    atomDens(ir,:) = DensOfR( d0s(:), rMesh(ir) )
  end do

  ! Find exchange and correlation energy and potential from radial density
  call atomXC( irel, nr, nr, rMesh, nSpin, atomDens, &
               atomEx, atomEc, atomDx, atomDc, atomVxc )

  ! Define fcc unit cell, such that a sphere of radius rMax+rBuff fits in it
  latConst = (rMax+rBuff) * 2*sqrt(2._dp)
  cell(:,1) = (/ 0.0_dp, 0.5_dp, 0.5_dp /)
  cell(:,2) = (/ 0.5_dp, 0.0_dp, 0.5_dp /)
  cell(:,3) = (/ 0.5_dp, 0.5_dp, 0.0_dp /)
  cell(:,:) = cell(:,:) * latConst

  ! Define reciprocal unit cell
  recCell(:,1) = (/-1.0_dp, 1.0_dp, 1.0_dp /)
  recCell(:,2) = (/ 1.0_dp,-1.0_dp, 1.0_dp /)
  recCell(:,3) = (/ 1.0_dp, 1.0_dp,-1.0_dp /)
  recCell(:,:) = recCell(:,:) * 2*pi/latConst
  kCut = cellMesh(1) * sqrt(sum(recCell(:,1)**2)) / 2  ! Max. wave vector
  Ecut = kCut**2                                       ! Mesh cutoff, in Ry
  dx = pi / kCut                                 ! Dist. between mesh planes

  ! Find the box of mesh points own by my processor
#ifdef MPI
  ! Do simplest thing: divide only along first axis
  nxNode = Nx / nNodes          ! Points per node along first vector
  nLarger = nx - nxNode*nNodes  ! Number of nodes with one more point
  if (myNode<nLarger) then      ! My node has nx+1 points
    lb1 = (nxNode+1)*(myNode-1)
    ub1 = (nxNode+1)*(myNode-1) - 1
  else                          ! My node has nx points
    lb1 = (nxNode+1)*nLarger + nxNode*(myNode-nLarger)
    ub1 = (nxNode+1)*nLarger + nxNode*(myNode-nLarger+1) - 1
  end if
#else
  ! All points belong to the only processor
  lb1 = 0
  ub1 = nx-1
#endif
  lb2 = 0
  lb3 = 0
  ub2 = nx-1
  ub3 = nx-1

  ! Allocate arrays for density and potential
  allocate( cellDens(lb1:ub1,lb2:ub2,lb3:ub3,nSpin), &
             cellVxc(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) )

  ! Find density at mesh points
  x0(:) = sum(cell,2) / 2     ! Center of cell
  do i3 = lb3,ub3
  do i2 = lb2,ub2
  do i1 = lb1,ub1
    x(:) = i1*cell(:,1)/cellMesh(1) &   ! Mesh point position
         + i2*cell(:,2)/cellMesh(2) &
         + i3*cell(:,3)/cellMesh(3)
    r = sqrt( sum((x-x0)**2) )          ! Distance to center of cell
    cellDens(i1,i2,i3,:) = DensOfR( d0s(:), r )
  end do ! i1
  end do ! i2
  end do ! i3

  ! Find exchange and correlation energy and potential from density in cell
  call cellXC( irel, cell, cellMesh, lb1, ub1, lb2, ub2, lb3, ub3, nSpin, &
               cellDens, cellEx, cellEc, cellDx, cellDc, stress, cellVxc )

  ! Print parameters
  if (myNode==0) then
    print'(/,a,10(a3,4x))', 'funcs= ', func(indexf)
    print  '(a,10(a6,1x))', 'auths= ', auth(indexf)
    print'(a,3f12.6)', 'dr, dx, Ecut = ', dr, dx, Ecut
    print'(a,2f12.6)', 'rMax, rBuff = ', rMax, rBuff
  end if

  ! Compare energies
  if (myNode==0) then
    print'(a,3f15.9)', 'atomEx, cellEx, diff =', atomEx, cellEx, atomEx-cellEx
    print'(a,3f15.9)', 'atomEc, cellEc, diff =', atomEc, cellEc, atomEc-cellEc
    print'(a,3f15.9)', 'atomDx, cellDx, diff =', atomDx, cellDx, atomDx-cellDx
    print'(a,3f15.9)', 'atomDc, cellDc, diff =', atomDc, cellDc, atomDc-cellDc
  end if

  ! Write potentials
  if (nNodes==1) then
    open( unit=44, file='atomVxc.out' )
    do ir = 1,nr
      if (sum(atomDens(ir,1:mSpin)) < densMin) cycle
      write(44,'(3f15.9)') rMesh(ir), atomVxc(ir,1), atomVxc(ir,mSpin)
    end do
    close( unit=44 )
    open( unit=44, file='cellVxc.out' )
    do i1 = cellMesh(1)/2,0,-1
      i2 = cellMesh(2) / 2
      i3 = cellMesh(3) / 2
      if (sum(cellDens(i1,i2,i3,1:mSpin)) < densMin) cycle
      x(:) = i1*cell(:,1)/cellMesh(1) &   ! Mesh point position
           + i2*cell(:,2)/cellMesh(2) &
           + i3*cell(:,3)/cellMesh(3)
      write(44,'(3f15.9)') &
        sqrt(sum((x-x0)**2)), cellVxc(i1,i2,i3,1), cellVxc(i1,i2,i3,mSpin)
    end do
    close( unit=44 )
  end if

  ! Compare potentials
  sumDiffVxc = 0
  maxDiffVxc = 0
  do i3 = lb3,ub3
  do i2 = lb2,ub2
  do i1 = lb1,ub1
    if (sum(cellDens(i1,i2,i3,1:mSpin)) < densMin) cycle
    x(:) = cell(:,1)*i1/cellMesh(1) &
         + cell(:,2)*i2/cellMesh(2) &
         + cell(:,3)*i3/cellMesh(3)
    r = sqrt(sum((x-x0)**2))
    if (r>=rMax) cycle
    ! Simplest thing: a linear interpolation of atomVxc (requires large nr)
    ir = r/dr + 1
    wr = (r - ir*dr) / dr
    Vxc(:) = atomVxc(ir,:)*(1-wr) + atomVxc(ir+1,:)*wr
    diffVxc(:) = abs(cellVxc(i1,i2,i3,:)-Vxc(:))
    if (maxval(diffVxc(:)) > maxDiffVxc) then
      i1max = i1
      i2max = i2
      i3max = i3
      irmax = ir
    end if
    sumDiffVxc = sumDiffVxc + sum(diffVxc(:)**2)
    maxDiffVxc = max( maxDiffVxc, maxval(diffVxc(:)) )
  end do ! i1
  end do ! i2
  end do ! i3
#ifdef MPI
  ! Find sumDiffVxc and maxDiffVxc accross all processor nodes
  tmp = sumDiffVxc
  call MPI_AllReduce( tmp, sumDiffVxc, 1, MPI_double_precision, &
                      MPI_Sum, MPI_Comm_World, MPIerror )
  tmp = maxDiffVxc
  call MPI_AllReduce( tmp, maxDiffVxc, 1, MPI_double_precision, &
                      MPI_Max, MPI_Comm_World, MPIerror )
#endif
  avgDiffVxc = sumDiffVxc / nSpin / nx**3
  if (myNode==0) then
    print'(a,2f15.9)', 'avgDiffVxc, maxDiffVxc = ', avgDiffVxc, maxDiffVxc
!    print'(a,4i6)', 'i1max,i2max,i3max, irmax = ', i1max, i2max, i3max, irmax
  end if

! Finalize MPI
#ifdef MPI
  call MPI_Finalize( MPIerror )
#endif

CONTAINS

FUNCTION DensOfR( d0, r )

  ! Returns a radial density distribution

  implicit none
  real(dp),intent(in):: d0(nSpin)  ! Density at center of charge distribution
  real(dp),intent(in):: r          ! Distance to center of charge distribution
  real(dp)           :: DensOfR(nSpin)  ! Electron density

  ! Use a simple gaussian distribution
  DensOfR = d0 * exp(-r**2/2/dWidth**2)

  ! Impose a smooth radial cutoff
  DensOfR = DensOfR * ( 1 - (r/rMax)**n1cut )**n2cut

END FUNCTION DensOfR

END PROGRAM siestaXCtest4

