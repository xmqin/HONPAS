! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program protoNEB

! 'Prototype' nudget elastic band (NEB) using siesta as a subroutine.
! But not really a NEB: it only finds the energy of the NH3 molecule for 
! several geometries, to illustrate the use of siesta as a subroutine. 
! J.M.Soler & J.Nara, May.2014

  use fsiesta, only: siesta_forces ! Calculate atomic forces, energy, and stress
  use fsiesta, only: siesta_launch ! Start a siesta process
  use fsiesta, only: siesta_quit   ! Finish siesta process
  use fsiesta, only: siesta_units  ! Set physical units
#ifdef MPI
  use mpi_interfaces, only: MPI_AllReduce
  use mpi_interfaces, only: MPI_Bcast
  use mpi_interfaces, only: MPI_Comm_World
  use mpi_interfaces, only: MPI_Comm_Rank
  use mpi_interfaces, only: MPI_Comm_Size
  use mpi_interfaces, only: MPI_Comm_Split
  use mpi_interfaces, only: MPI_Double_Precision
  use mpi_interfaces, only: MPI_Finalize
  use mpi_interfaces, only: MPI_Init
  use mpi_interfaces, only: MPI_Min
#endif

  implicit none

! Internal parameters and variables
  integer,parameter:: dp = kind(1.d0)
  integer,parameter:: na = 4             ! number of atoms
  integer,parameter:: np = 3             ! number of points (nudges)
  integer:: ia, ih, ip, ix, maxLen, MPIerror, myNode, myPoint, &
            mySiestaComm, mySiestaNode, nNodes, nSiestaNodes, numLen
  real(dp):: cell(3,3), cellSize, energy(np), fa(3,na,np), &
             phi, pi, rNH, theta, thetaMax, thetaMin, xa(3,na,np)
  real(dp):: recvBuff(3*na*np), sendBuff(3*na*np)
  character(len=20):: numName, nodeString
  character(len=12):: myLabel

! Open output file
  close(unit=6)
  open(unit=6, file='protoNEB.out')

! Set internal parameters
  pi = acos(-1._dp)
  thetaMin = pi/3           ! min. angle between NH bonds and symmetry axis
  thetaMax = pi/2           ! max. angle between NH bonds and symmetry axis
  rNH = 1.0_dp              ! N-H bond length (Ang)
  cellSize = 10._dp         ! cell size (Ang)

! Set NH3 geometries
  xa = 0                    ! nitrogen atom (atom 4) at origin
  do ip = 1,np              ! loop on geometries (points)
    theta = thetaMin + (ip-1)*(thetaMax-thetaMin)/(np-1)
    do ih = 1,3             ! loop on Hydrogen atoms
      phi = 2*pi*(ih-1)/3._dp
      xa(1,ih,ip) = rNH*sin(theta)*cos(phi)
      xa(2,ih,ip) = rNH*sin(theta)*sin(phi)
      xa(3,ih,ip) = rNH*cos(theta)
    enddo
  enddo
  cell = 0                  ! cell unit vectors
  do ix = 1,3
    cell(ix,ix) = cellSize
  enddo

! Initialize MPI and set siesta communicators
#ifdef MPI
  call MPI_Init( MPIerror )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, nNodes, MPIerror )
  if (mod(nNodes,np)/=0) &
    stop 'protoNEB ERROR: #Cores must be a multiple of #Geometries'
  nSiestaNodes = nNodes/np
  myPoint = myNode/nSiestaNodes+1
  call MPI_Comm_Split( MPI_Comm_World, myPoint, 0, mySiestaComm, MPIerror )
  call MPI_Comm_Rank( mySiestaComm, mySiestaNode, MPIerror )
#else
  stop 'protoNEB ERROR: this version must run with MPI'
#endif

! Set siesta system labels
  do ip = 1,np
    ! Find name of this node's number, and its name length
    write(numName,*) myPoint
    numName = adjustl(numName)
    numLen = len_trim(numName)

    ! Set node number string, e.g. 00, 01, ... 63
    maxLen = 2
    nodeString = '00'
    nodeString(maxLen-numLen+1:maxLen) = trim(numName)
    nodeString(maxLen+1:) = ' '
    myLabel = 'NH3-point'//nodeString(1:maxLen)
  enddo
    
! Start siesta processes and set the physical units
  if (myNode==0) write(6,*) 'protoNEB: calling siesta_launch'
  call siesta_launch( myLabel, mpi_comm=mySiestaComm )
  ! alternatively, to let siesta create its communicators:
  ! call siesta_launch( myLabel, nSiestaNodes ) 
  call siesta_units( 'Ang', 'eV' )

! Find energy and forces with my geometry
  call siesta_forces( myLabel, na, xa(:,:,myPoint), cell, &
                      energy(myPoint), fa(:,:,myPoint) )

! Gather energies and forces at all geometries
#ifdef MPI
  sendBuff = 1.0e100_dp
  sendBuff(myPoint) = energy(myPoint)
  call MPI_AllReduce( sendBuff, recvBuff, np, MPI_double_precision, &
                      MPI_Min, MPI_Comm_World, MPIerror )
  energy = recvBuff(1:np)
  sendBuff = 1.0e100_dp
  sendBuff(3*na*(myPoint-1)+1:3*na*myPoint) = reshape(fa(:,:,myPoint),(/3*na/))
  call MPI_AllReduce( sendBuff, recvBuff, 3*4*np, MPI_double_precision, &
                      MPI_Min, MPI_Comm_World, MPIerror )
  fa = reshape(recvBuff,(/3,na,np/))
#endif

! Print energies and forces
  if (myNode==0) then
    do ip = 1,np
      write(6,'(/a,i4)')   'protoNEB: Point', ip
      write(6,'(a,f12.6)') 'energy (eV) =', energy(ip)
      write(6,'(a,/,(i4,3f12.6))') 'forces (eV/Ang) =', (ia,fa(:,ia,ip),ia=1,na)
    enddo
  endif

! Stop siesta processes
  call siesta_quit( myLabel )

! Finalize MPI and close output file
#ifdef MPI
  call MPI_Finalize( MPIerror )
#endif
  close(unit=6)

end program protoNEB

