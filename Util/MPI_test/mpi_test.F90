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
program mpi_test

use parallel, only: Node, Nodes

#ifdef MPI
use mpi_siesta, only: mpi_init, mpi_comm_rank, mpi_comm_size
use mpi_siesta, only: mpi_comm_world
#endif

use precision, only: dp
use m_mpi_utils, only: broadcast, globalize_sum

#ifdef MPI
integer :: mpierror
#endif

real(dp) :: a(1:3,0:2), b(1:3,0:2)

#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      node = 0
      nodes = 1
#endif


if (Node == 0) then
   print *, "Working with ", nodes, " nodes."
   call random_number(a)
   print *, a
   print *, "----------------"
endif

call broadcast(a(:,:))

if (Node == 1)  print *, a
      
call globalize_sum(a,b)

if (Node == 0) then
   print *, "-After reduction---------------"
   print *, b
   print *, "-Should equal-------------"
   print *, nodes*a
endif

end program mpi_test
