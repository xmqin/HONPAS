PROGRAM siestaXCtest2

  ! Compares potential and the numerical derivative of the energy 
  ! calculated by atomXC. J.M.Soler. Sept.2009

  ! Used module procedures
  USE siestaXC, only: atomXC
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
  integer, parameter:: irel  =  0 ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nSpin =  2 ! Number of spin components
  integer, parameter:: nfTot = 19 ! Number of functionals
  integer, parameter:: nr    = 101 ! Number of radial points
  integer, parameter:: n1cut =  8 ! Cutoff parameter
  integer, parameter:: n2cut =  2 ! Cutoff parameter:
                                  !    fCut(r)=(1-(r/rMax)**n1cut)**n2cut
  real(dp),parameter:: dWidth = 2._dp ! Width of density distribution, in Bohr
  real(dp),parameter:: Qtot = 10._dp  ! Integral of density distribution
  real(dp),parameter:: spinPol= 2._dp ! Integral of densUp - densDown
  real(dp),parameter:: rMax = 20._dp  ! Cutoff radius, in Bohr
  real(dp),parameter:: deltaDens = 1.e-8_dp  ! Finite diff. change
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
  integer :: iDelta, ir, irmax, ismax, iSpin, one, two
  real(dp):: avgDiffVxc, dDensdr, dens(nr,nSpin), dens0(nr,nSpin), &
             d0tot, d0(nSpin), dEdDens, dDens, diffVxc, &
             Dc, Dc0, dr, dVol, Dx, Dx0, Ec, Ec0, Ex, Ex0, &
             kf, kg, maxDiffVxc, pi, r, rMesh(nr), &
             Vxc(nr,nSpin), Vxc0(nr,nSpin), wc(nfTot), wr, wx(nfTot)

#ifdef MPI
  ! Initialize MPI, even though this test is intended to be run serially
  integer:: MPIerror, myNode, nNodes
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
#endif

  ! Initialize hybrid XC functional with all tested functionals
  wx = 1._dp / nf
  wc = 1._dp / nf
  call setXC( nf, func(indexf), auth(indexf), wx(indexf), wc(indexf) )

  ! Find radial mesh points and gaussian density
  pi = acos(-1._dp)
  d0tot = Qtot / (2*pi*dWidth**2)**1.5_dp    ! Total density at origin
  if (nSpin==1) then
    d0(1) = d0tot
  else
    one = 1   ! A silly thing to satisfy the compiler when nSpin=1
    two = 2
    d0(one) = d0tot * (Qtot + spinPol) / Qtot / 2 ! Spin up density at origin
    d0(two) = d0tot * (Qtot - spinPol) / Qtot / 2 ! Spin down density at origin
  end if
  dr = rmax / (nr-1)                      ! Interval between radial points
  do ir = 1,nr
    rMesh(ir) = dr * (ir-1)               ! Radial point values
    dens0(ir,:) = DensOfR( d0(:), rMesh(ir) )
  end do

  ! Find exchange and correlation energy and potential from radial density
  call atomXC( irel, nr, nr, rMesh, nSpin, dens0, Ex0, Ec0, Dx0, Dc0, Vxc0 )

  ! Print parameters
  print'(/,a,10(a3,4x))', 'funcs= ', func(indexf)
  print  '(a,10(a6,1x))', 'auths= ', auth(indexf)
  print'(a,2f12.6)', 'dr, rMax = ', dr, rMax

  ! Calculate finite-difference derivatives
  open( unit=44, file='test2.Vxc' )
  avgDiffVxc = 0
  maxDiffVxc = 0
  do iSpin = 1,nSpin
    do ir = 2,nr
      if (dens0(ir,iSpin)<densMin) cycle ! Do nothing if dens=0
      dVol = 4*pi*rMesh(ir)**2 * dr
      dDens = min( deltaDens, dens0(ir,iSpin)/100 )
      dEdDens = 0
      do iDelta = -1,1,2
        dens = dens0
        dens(ir,iSpin) = dens0(ir,iSpin) + iDelta*dDens
        call atomXC( irel, nr, nr, rMesh, nSpin, dens, Ex, Ec, Dx, Dc, Vxc )
        dEdDens = dEdDens + iDelta * (Ex+Ec) / (2*dDens) / dVol
      end do ! iDelta
      diffVxc = Vxc0(ir,iSpin) - dEdDens
      avgDiffVxc = avgDiffVxc + diffVxc**2
      if (abs(diffVxc) > maxDiffVxc) then
        maxDiffVxc = abs(diffVxc)
        irMax = ir
        isMax = iSpin
      end if
      kf = (3*pi**2*sum(dens0(ir,:)))**(1.0_dp/3)
      dDensdr = sum(dens0(ir+1,:)-dens0(ir-1,:))/(rMesh(ir+1)-rMesh(ir-1))
      kg = abs(dDensdr)/sum(dens0(ir,:))
      if (ir==2) then
        print'(a5,a9,4a15)','iSpin','r','dens','Vxc','dExc/dDens','diff'
        print'(i5,f9.3,4f15.9)', &
          ispin, rMesh(1), dens0(1,iSpin), Vxc0(1,iSpin)
!        print'(a5,a9,6a15)','iSpin','r','dens','kf','kg','Vxc', &
!                            'dExc/dDens','diff'
!        print'(i5,f9.3,6f15.9)', &
!          ispin, rMesh(1), dens0(1,iSpin), kf, kg, Vxc0(1,iSpin)
      end if
      print'(i5,f9.3,4f15.9)', &
        ispin, rMesh(ir), dens0(ir,iSpin), Vxc0(ir,iSpin), dEdDens, diffVxc
!      print'(i5,f9.3,6f15.9)', &
!        ispin, rMesh(ir), dens0(ir,iSpin), kf, kg, &
!        Vxc0(ir,iSpin), dEdDens, diffVxc
      write(44,'(f9.3,4f15.9)') &
        rMesh(ir), dens0(ir,iSpin), Vxc0(ir,iSpin), dEdDens, diffVxc
    end do ! ir
  end do ! iSpin
  avgDiffVxc = sqrt( avgDiffVxc / nSpin / nr )
  print'(a,2f15.9)', 'avgDiffVxc, maxDiffVxc = ', avgDiffVxc, maxDiffVxc
!  print'(a,2i6)', 'irMax, iSpinMax = ', irmax, ismax
  close( unit=44 )

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

END PROGRAM siestaXCtest2

