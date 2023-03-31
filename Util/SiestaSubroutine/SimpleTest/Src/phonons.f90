! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program phonons

! Compute the force-constant matrix in parallel
! This version uses MPI and siesta as a subroutine. It must be compiled
! together with siesta.
! A. Garcia (2011)

  use mpi
  use fsiesta
  use sys, only: die

  implicit none
  integer,parameter:: dp = kind(1.d0)

  integer,parameter :: na = 3
  integer :: error, ia, i
  integer :: myNode, MPI_Comm, Nodes
  real(dp):: e, fa(3,na), xa(3,na), faplus(3,na), faminus(3,na)
  real(dp):: xaflat(3*na), fc(3*na,3*na), fcbuff(3*na)
  real(dp), parameter :: delta_u = 0.05

!
! Original coordinates
!
  data xa / 0.0, 0.0, 0.0, &
            0.7, 0.7, 0.0, &
           -0.7, 0.7, 0.0 /

! Initialize MPI and get my node's index
  call MPI_Init( error )
  call MPI_Comm_Size( MPI_Comm_World, Nodes, error )
  if (Nodes /= 9) call die("You need 9 nodes!")

  call MPI_Comm_Rank( MPI_Comm_World, myNode, error )

! Split the communicator to have as many new contexts as nodes 
! (This is only for demonstration purposes -- in real life one
!  might want to have each team comprise a few processors)

! In this example one should have 9 processors, or else
! things would not work. In real life one should decide beforehand
! how to distribute the load.

  call MPI_Comm_Split( MPI_Comm_World, myNode, myNode, MPI_Comm, error )

! Set physical units
  call siesta_units( 'Ang', 'eV' )

! Launch a siesta process in each node
! Note that each process is using the same fdf file. This is
! not optimal in the current version, since every processor
! will attempt to create files with the same names!
! There should be a more flexible framework for this

  call siesta_launch( 'h2o', mpi_comm= MPI_Comm )
  if (myNode==0) print*, 'siesta launched'

!
!  Compute a column of the force constant matrix in each processor
!  (Each processor takes care of a coordinate)
  
! Find forces for + displacement
  xaflat = reshape(xa,(/3*na/))
  xaflat(1+myNode) = xaflat(1+myNode) + delta_u
  xa = reshape(xaflat,(/3,3/))
  call siesta_forces( 'h2o', na, xa, energy=e, fa=faplus )
  print *, "Node ", myNode, " Eplus: ", e

! Find forces for - displacement
  xaflat(1+myNode) = xaflat(1+myNode) - 2* delta_u
  xa = reshape(xaflat,(/3,3/))
  call siesta_forces( 'h2o', na, xa, energy=e, fa=faminus )
  print *, "Node ", myNode, " Eminus: ", e
  
  fa = 0.5*(faplus-faminus)/delta_u
  fcbuff(1:3*na) = reshape(fa,(/3*na/))

! Done computing with Siesta
  call siesta_quit( 'h2o' )

! Gather all the pieces in the root node (of the original communicator!!)

  call mpi_gather(fcbuff(1:3*na),3*na,MPI_double_precision, &
                  fc,3*na,MPI_double_precision, 0, MPI_Comm_World, error)

  if (myNode==0) then
    do i = 1, 3*na     
       print "(9g12.4)", fc(i,:)
    enddo
 endif

! Finalize MPI
  call MPI_Finalize( error )

end program phonons

