PROGRAM siestaXCtest3

  ! Tester of cellXC. Compares the potential and stress with the numerical
  ! (finite-difference) derivatives of the energy. J.M.Soler. Sept.2009

  ! Used module procedures
  USE siestaXC, only: cellXC
  USE siestaXC, only: myMeshBox
  USE siestaXC, only: nfft
  USE siestaXC, only: setMeshDistr
  USE siestaXC, only: setXC

  ! Used module parameters
  USE siestaXC, only: dp
  USE siestaXC, only: gp => grid_p

! Used MPI procedures and types
#ifdef MPI
  USE mpi_siesta, only: MPI_AllReduce
  USE mpi_siesta, only: MPI_Comm_World
  USE mpi_siesta, only: MPI_Double_Precision
  USE mpi_siesta, only: MPI_Sum
#endif

  implicit none

  ! Tester parameters
  integer, parameter:: irel  = 0         ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nSpin = 2         ! Number of spin components
  integer, parameter:: nfTot = 19        ! Number of functionals available
  integer, parameter:: nRan =  6         ! Number of random points for test
  real(dp),parameter:: latConst = 10._dp ! Lattice constant, in Bohr
  real(dp),parameter:: Ecut = 30._dp     ! Planewave cutoff of the mesh
  real(dp),parameter:: Qtot = 10._dp     ! Integral of density distribution
  real(dp),parameter:: D1byD0 = 0.99_dp  ! Fractional change of density
  real(dp),parameter:: spinPol= 2._dp    ! Integral of densUp - densDown
  real(dp),parameter:: latCons= 10._dp   ! Lattice constant in Bohr
  real(dp),parameter:: asym = 0.5_dp     ! Orthorhombic asymetry 
  real(dp),parameter:: deltaDens = 1.e-6_dp ! Used for numerical derivatives
  real(dp),parameter:: deltaStrain = 1.e-4_dp ! Used for numerical derivatives

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

  ! A few random numbers
  real(dp):: ran(nRan) = (/0.749218032_dp, 0.928517579_dp, 0.043866380_dp, &
                           0.669289084_dp, 0.392851044_dp, 0.502184791_dp /)

  ! Tester variables and arrays
  logical :: pointIsMine
  integer :: i, i1, i2, i3, ic, iDelta, ip, iRan, iSpin, lb1, lb2, lb3, &
             myBox(2,3), myDistr, myNode, &
             nMesh(3), nNodes, np, ub1, ub2, ub3
  real(dp):: cell(3,3), cell0(3,3), d0(nSpin), Dc, Dc0, &
             dEdDens, dEdStrain(3,3), dVol, Dx, Dx0, dxMax, &
             Ec, Ec0, Ex, Ex0, kCut, pi, &
             strain(3,3), stress(3,3), stress0(3,3), &
             tmp, volume, VxcSpin, wc(nfTot), wx(nfTot), x(3)
  real(gp),allocatable:: dens(:,:,:,:), dens0(:,:,:,:), &
                         Vxc(:,:,:,:), Vxc0(:,:,:,:)

#ifdef MPI
  ! MPI-related variables
  integer:: MPIerror
#endif

  ! Initialize hybrid XC functional with all tested functionals
  wx = 1._dp / nf
  wc = 1._dp / nf
  call setXC( nf, func(indexf), auth(indexf), wx(indexf), wc(indexf) )

  ! Define a near-orthorhombic unit cell
  ! This is just to assume orthorhombicity in calculating nMesh from Ecut
  cell(:,:) = 0
  cell(1,1) = 1 - asym
  cell(2,2) = 1
  cell(3,3) = 1 + asym
  cell(:,1) = cell(:,1) + (/ 0.000_dp,  0.120_dp, -0.175_dp /)
  cell(:,2) = cell(:,2) + (/ 0.221_dp,  0.000_dp,  0.030_dp /)
  cell(:,3) = cell(:,3) + (/ 0.116_dp, -0.231_dp,  0.000_dp /)
  cell = cell * latConst
  volume = ( cell(2,1)*cell(3,2) - cell(3,1)*cell(2,2) ) * cell(1,3)  &
         + ( cell(3,1)*cell(1,2) - cell(1,1)*cell(3,2) ) * cell(2,3)  &
         + ( cell(1,1)*cell(2,2) - cell(2,1)*cell(1,2) ) * cell(3,3)

#ifdef MPI
  ! Initialize MPI and get myNode and nNodes
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
#else
  myNode = 0
  nNodes = 1
#endif

  ! Find the number of mesh points in each direction
  pi = acos(-1._dp)
  kCut = sqrt(Ecut)
  dxMax = pi / kCut
  do ic = 1,3
    nMesh(ic) = ceiling( cell(ic,ic) / dxMax )
    call nfft( nMesh(ic) )
  end do

  ! Print parameters
  if (myNode==0) then
    print'(/,a,10(a3,4x))', 'funcs= ', func(indexf)
    print  '(a,10(a6,1x))', 'auths= ', auth(indexf)
    print'(a,3i6,f12.6)', 'nMesh, Ecut =', nMesh, Ecut
  end if

  ! Find the box of mesh points own by my processor
  call setMeshDistr( myDistr, nMesh )
  call myMeshBox( nMesh, myDistr, myBox )
  lb1 = myBox(1,1)
  ub1 = myBox(2,1)
  lb2 = myBox(1,2)
  ub2 = myBox(2,2)
  lb3 = myBox(1,3)
  ub3 = myBox(2,3)

  ! Allocate arrays for density and potential
  allocate( dens(lb1:ub1,lb2:ub2,lb3:ub3,nSpin), &
           dens0(lb1:ub1,lb2:ub2,lb3:ub3,nSpin), &
             Vxc(lb1:ub1,lb2:ub2,lb3:ub3,nSpin), &
            Vxc0(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) )

  ! Find density at mesh points
  np = product(nMesh)
  dVol = volume / np
  if (nSpin==1) then
    d0(1) = Qtot / volume
  else
    d0(1) = (Qtot+spinPol) / 2 / volume
    d0(2) = (Qtot-spinPol) / 2 / volume
  end if
  do i3 = lb3,ub3
  do i2 = lb2,ub2
  do i1 = lb1,ub1
    x(:) = cell(:,1)*i1/nMesh(1) &   ! Mesh point position
         + cell(:,2)*i2/nMesh(2) &
         + cell(:,3)*i3/nMesh(3)
    dens(i1,i2,i3,:) = d0(:) * (1 + D1byD0 * cos( 2*pi*x(1)/cell(1,1) ) &
                                           * cos( 2*pi*x(2)/cell(2,2) ) &
                                           * cos( 2*pi*x(3)/cell(3,3) ) )
  end do ! i1
  end do ! i2
  end do ! i3

  ! Find exchange and correlation energy and potential from density in cell
  call cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, nSpin, &
               dens, Ex, Ec, Dx, Dc, stress, Vxc )

  ! Print total energies
  if (myNode==0) then
    print'(a,2f15.9)', 'Ex, Ec =', Ex, Ec
    print'(a,2f15.9)', 'Dx, Dc =', Dx, Dc
  end if

  ! Store unperturbed magnitudes
  cell0 = cell
  dens0 = dens
  Ex0 = Ex
  Ec0 = Ec
  Dx0 = Dx
  Dc0 = Dc
  stress0 = stress
  Vxc0 = Vxc

  ! Compare Vxc with numerical derivative of d(Ex+Ec)/dDens
  do iRan = 1,nRan

    ! Choose one spin component
    if (ran(iRan)< 0.5_dp) iSpin = 1
    if (ran(iRan)>=0.5_dp) iSpin = 2
    ran(iRan) = 2*mod(ran(iRan),0.5_dp)  ! ran(iRan) again in range (0:1)

    ! Choose one mesh point. Notice that index ranges are (0:nMesh-1)
    ip = np * ran(iRan)
    i3 = ip / nMesh(1) / nMesh(2)  ! Undo ip=i1+nMesh(1)*i2+nMesh(1)*nMesh(2)*i3
    ip = ip - nMesh(1) * nMesh(2) * i3
    i2 = ip / nMesh(1)
    i1 = ip - nMesh(1) * i2
    ! Find if point i1,i2,i3 is in my box
    pointIsMine = (lb1<=i1 .and. i1<=ub1 .and. &
                   lb2<=i2 .and. i2<=ub2 .and. &
                   lb3<=i3 .and. i3<=ub3)

    ! Find energies with dens(i1,i2,i3,iSpin)+/-delta
    dEdDens = 0
    do iDelta = -1,1,2
      dens = dens0
      if (pointIsMine) &
        dens(i1,i2,i3,iSpin) = dens0(i1,i2,i3,iSpin) + iDelta*deltaDens
      call cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, nSpin, &
                   dens, Ex, Ec, Dx, Dc, stress, Vxc )
      dEdDens = dEdDens + iDelta*(Ex+Ec)/(2*deltaDens)
    end do ! iDelta

#ifdef MPI
    ! Send Vxc(i1,i2,i3,iSpin) to other nodes
    if (pointIsMine) then
      VxcSpin = Vxc0(i1,i2,i3,iSpin)
    else
      VxcSpin = 0
    end if
    tmp = VxcSpin
    call MPI_AllReduce( tmp, VxcSpin, 1, MPI_double_precision, &
                        MPI_Sum, MPI_Comm_World, MPIerror )
#else
    VxcSpin = Vxc0(i1,i2,i3,iSpin)
#endif
    if (myNode==0) then
      if (iRan==1) print'(4a6,3a15)', &
        'i1','i2','i3','iSpin','Vxc','dExc/dDens','diff'
      print'(4i6,3f15.9)', &
        i1, i2, i3, iSpin, VxcSpin, dEdDens/dVol, VxcSpin-dEdDens/dVol
    end if

  end do ! iRan

  ! Compare stress with numerical derivative of dE/dStrain
  dens = dens0
  do i2 = 1,3
    do i1 = 1,3
      dEdStrain(i1,i2) = 0
      do iDelta = -1,1,2
        strain(:,:) = 0
        forall(i=1:3) strain(i,i) = 1
        strain(i1,i2) = strain(i1,i2) + iDelta*deltaStrain
        cell = matmul( strain, cell0 )
        call cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, nSpin, &
                     dens, Ex, Ec, Dx, Dc, stress, Vxc )
        dEdStrain(i1,i2) = dEdStrain(i1,i2) + iDelta*(Ex+Ec)/(2*deltaStrain)
      end do ! iDelta
    end do ! i1
  end do ! i2
  if (myNode==0) then
    print'(a,/,(3f15.9))', 'xc stress * volume =', stress0*volume
    print'(a,/,(3f15.9))', 'dExc/dStrain =', dEdStrain
    print'(a,/,(3f15.9))', 'diff =', (stress0*volume-dEdStrain)
  end if

! Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

END PROGRAM siestaXCtest3

