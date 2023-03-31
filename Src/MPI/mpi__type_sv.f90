      MODULE MPI__type_SV
!
!
!       Michael Hennecke
!       A Fortran 90 interface to MPI version 1.1
!       RZ Uni Karlsruhe, Internal Report 63/96
!
!         Auxiliary module template MPI__type_SV
!
!         03-Oct-1996, hennecke@rz.uni-karlsruhe.de (v0.9c beta)
!
!       Permission is granted to copy and distribute this file
!       or modified versions of this file for no fee, provided the 
!       copyright notice and this permission notice are preserved 
!       on all copies.
!
!       (C) 1996  Michael Hennecke, RZ Universitaet Karlsruhe
!
!
        USE TIMER_MPI_M, only: timer_mpi

        IMPLICIT NONE
        PRIVATE

!       ... A.9  Fortran Bindings for Point-to-Point Communication  ...

!!!     PUBLIC :: MPI_PACK
!!!     INTERFACE MPI_PACK
!!!       MODULE PROCEDURE MPI_PACK_T
!!!     END INTERFACE ! MPI_PACK

!!!     PUBLIC :: MPI_UNPACK
!!!     INTERFACE MPI_UNPACK
!!!       MODULE PROCEDURE MPI_UNPACK_T
!!!     END INTERFACE ! MPI_UNPACK

!       ... A.10  Fortran Bindings for Collective Communication  ...

        PUBLIC :: MPI_GATHER
        INTERFACE MPI_GATHER
          MODULE PROCEDURE MPI_GATHER_T
        END INTERFACE ! MPI_GATHER

        PUBLIC :: MPI_GATHERV
        INTERFACE MPI_GATHERV
          MODULE PROCEDURE MPI_GATHERV_T
        END INTERFACE ! MPI_GATHERV

        PUBLIC :: MPI_SCATTER
        INTERFACE MPI_SCATTER
          MODULE PROCEDURE MPI_SCATTER_T
        END INTERFACE ! MPI_SCATTER

        PUBLIC :: MPI_SCATTERV
        INTERFACE MPI_SCATTERV
          MODULE PROCEDURE MPI_SCATTERV_T
        END INTERFACE ! MPI_SCATTERV

        PUBLIC :: MPI_ALLGATHER
        INTERFACE MPI_ALLGATHER
          MODULE PROCEDURE MPI_ALLGATHER_T
        END INTERFACE ! MPI_ALLGATHER

        PUBLIC :: MPI_ALLGATHERV
        INTERFACE MPI_ALLGATHERV
          MODULE PROCEDURE MPI_ALLGATHERV_T
        END INTERFACE ! MPI_ALLGATHERV

        PUBLIC :: MPI_ALLTOALL
        INTERFACE MPI_ALLTOALL
          MODULE PROCEDURE MPI_ALLTOALL_T
        END INTERFACE ! MPI_ALLTOALL

        PUBLIC :: MPI_ALLTOALLV
        INTERFACE MPI_ALLTOALLV
          MODULE PROCEDURE MPI_ALLTOALLV_T
        END INTERFACE ! MPI_ALLTOALLV

        PUBLIC :: MPI_REDUCE
        INTERFACE MPI_REDUCE
          MODULE PROCEDURE MPI_REDUCE_T
        END INTERFACE ! MPI_REDUCE

        PUBLIC :: MPI_ALLREDUCE
        INTERFACE MPI_ALLREDUCE
          MODULE PROCEDURE MPI_ALLREDUCE_T
        END INTERFACE ! MPI_ALLREDUCE

        PUBLIC :: MPI_REDUCE_SCATTER
        INTERFACE MPI_REDUCE_SCATTER
          MODULE PROCEDURE MPI_REDUCE_SCATTER_T
        END INTERFACE ! MPI_REDUCE_SCATTER

        PUBLIC :: MPI_SCAN
        INTERFACE MPI_SCAN
          MODULE PROCEDURE MPI_SCAN_T
        END INTERFACE ! MPI_SCAN

      CONTAINS

!       ... A.9  Fortran Bindings for Point-to-Point Communication  ...

        SUBROUTINE MPI_PACK_T(                                          &
     &      INBUF, INCOUNT, DATATYPE, OUTBUF, OUTCOUNT, POSITION,       &
     &      COMM, IERROR)
          type, INTENT(IN)  :: INBUF
          INTEGER, INTENT(IN)  :: INCOUNT
          INTEGER, INTENT(IN)  :: DATATYPE
          type, INTENT(OUT) :: OUTBUF(*)
          INTEGER, INTENT(IN)  :: OUTCOUNT
          INTEGER, INTENT(INOUT) :: POSITION
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_PACK
          CALL     MPI_PACK(                                            &
     &      INBUF, INCOUNT, DATATYPE, OUTBUF, OUTCOUNT, POSITION,       &
     &      COMM, IERROR)
        END SUBROUTINE MPI_PACK_T
        
        SUBROUTINE MPI_UNPACK_T(                                        &
     &      INBUF, INSIZE, POSITION, OUTBUF, OUTCOUNT, DATATYPE,        &
     &      COMM, IERROR)
          type, INTENT(IN)  :: INBUF
          INTEGER, INTENT(IN)  :: INSIZE
          INTEGER, INTENT(INOUT) :: POSITION
          type, INTENT(OUT) :: OUTBUF(*)
          INTEGER, INTENT(IN)  :: OUTCOUNT
          INTEGER, INTENT(IN)  :: DATATYPE
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_UNPACK
          CALL     MPI_UNPACK(                                          &
     &      INBUF, INSIZE, POSITION, OUTBUF, OUTCOUNT, DATATYPE,        &
     &      COMM, IERROR)
        END SUBROUTINE MPI_UNPACK_T
        
!      ... A.10  Fortran Bindings for Collective Communication  ...
        
        SUBROUTINE MPI_GATHER_T(                                        &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, ROOT, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNT
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: ROOT
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_GATHER
          call timer_mpi('MPI_GATHER',1)
          CALL     MPI_GATHER(                                          &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, ROOT, COMM, IERROR) 
          call timer_mpi('MPI_GATHER',2)
        END SUBROUTINE MPI_GATHER_T
        
        SUBROUTINE MPI_GATHERV_T(                                       &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,          &
     &      DISPLS, RECVTYPE, ROOT, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNTS(*)
          INTEGER, INTENT(IN)  :: DISPLS(*)
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: ROOT
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_GATHERV
          CALL     MPI_GATHERV(                                         &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,          &
     &      DISPLS, RECVTYPE, ROOT, COMM, IERROR) 
        END SUBROUTINE MPI_GATHERV_T
                  
        SUBROUTINE MPI_SCATTER_T(                                       &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, ROOT, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNT
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: ROOT
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_SCATTER
          CALL     MPI_SCATTER(                                         &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, ROOT, COMM, IERROR) 
        END SUBROUTINE MPI_SCATTER_T
        
        SUBROUTINE MPI_SCATTERV_T(                                      &
     &      SENDBUF, SENDCOUNTS, DISPLS, SENDTYPE, RECVBUF,             &
     &      RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNTS(*)
          INTEGER, INTENT(IN)  :: DISPLS(*)
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNT
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: ROOT
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_SCATTERV
          CALL     MPI_SCATTERV(                                        &
     &      SENDBUF, SENDCOUNTS, DISPLS, SENDTYPE, RECVBUF,             &
     &      RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR) 
        END SUBROUTINE MPI_SCATTERV_T
        
        SUBROUTINE MPI_ALLGATHER_T(                                     &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNT
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_ALLGATHER
          call timer_mpi('MPI_ALLGATHER',1)
          CALL     MPI_ALLGATHER(                                       &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, COMM, IERROR) 
          call timer_mpi('MPI_ALLGATHER',2)
        END SUBROUTINE MPI_ALLGATHER_T
        
        SUBROUTINE MPI_ALLGATHERV_T(                                    &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,          &
     &      DISPLS, RECVTYPE, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNTS(*)
          INTEGER, INTENT(IN)  :: DISPLS(*)
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_ALLGATHERV
          CALL     MPI_ALLGATHERV(                                      &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,          &
     &      DISPLS, RECVTYPE, COMM, IERROR) 
        END SUBROUTINE MPI_ALLGATHERV_T
        
        SUBROUTINE MPI_ALLTOALL_T(                                      &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNT
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNT
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_ALLTOALL
          CALL     MPI_ALLTOALL(                                        &
     &      SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,           &
     &      RECVTYPE, COMM, IERROR) 
        END SUBROUTINE MPI_ALLTOALL_T
        
        SUBROUTINE MPI_ALLTOALLV_T(                                     &
     &      SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE, RECVBUF,            &
     &      RECVCOUNTS, RDISPLS, RECVTYPE, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          INTEGER, INTENT(IN)  :: SENDCOUNTS(*)
          INTEGER, INTENT(IN)  :: SDISPLS(*)
          INTEGER, INTENT(IN)  :: SENDTYPE
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNTS(*)
          INTEGER, INTENT(IN)  :: RDISPLS(*)
          INTEGER, INTENT(IN)  :: RECVTYPE
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_ALLTOALLV
          CALL     MPI_ALLTOALLV(                                       &
     &      SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE, RECVBUF,            &
     &      RECVCOUNTS, RDISPLS, RECVTYPE, COMM, IERROR) 
        END SUBROUTINE MPI_ALLTOALLV_T
        
        SUBROUTINE MPI_REDUCE_T(                                        &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: COUNT
          INTEGER, INTENT(IN)  :: DATATYPE
          INTEGER, INTENT(IN)  :: OP
          INTEGER, INTENT(IN)  :: ROOT
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_REDUCE
          call timer_mpi('MPI_REDUCE',1)
          CALL     MPI_REDUCE(                                          &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR) 
          call timer_mpi('MPI_REDUCE',2)
        END SUBROUTINE MPI_REDUCE_T
        
        SUBROUTINE MPI_ALLREDUCE_T(                                     &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: COUNT
          INTEGER, INTENT(IN)  :: DATATYPE
          INTEGER, INTENT(IN)  :: OP
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_ALLREDUCE
          call timer_mpi('MPI_ALLREDUCE',1)
          CALL     MPI_ALLREDUCE(                                       &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR) 
          call timer_mpi('MPI_ALLREDUCE',2)
        END SUBROUTINE MPI_ALLREDUCE_T
        
        SUBROUTINE MPI_REDUCE_SCATTER_T(                                &
     &      SENDBUF, RECVBUF, RECVCOUNTS, DATATYPE, OP, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: RECVCOUNTS(*)
          INTEGER, INTENT(IN)  :: DATATYPE
          INTEGER, INTENT(IN)  :: OP
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_REDUCE_SCATTER
          CALL     MPI_REDUCE_SCATTER(                                  &
     &      SENDBUF, RECVBUF, RECVCOUNTS, DATATYPE, OP, COMM, IERROR) 
        END SUBROUTINE MPI_REDUCE_SCATTER_T
        
        SUBROUTINE MPI_SCAN_T(                                          &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR) 
          type, INTENT(IN)  :: SENDBUF
          type, INTENT(OUT) :: RECVBUF(*)
          INTEGER, INTENT(IN)  :: COUNT
          INTEGER, INTENT(IN)  :: DATATYPE
          INTEGER, INTENT(IN)  :: OP
          INTEGER, INTENT(IN)  :: COMM
          INTEGER, INTENT(OUT) :: IERROR 
          EXTERNAL MPI_SCAN
          CALL     MPI_SCAN(                                            &
     &      SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR) 
        END SUBROUTINE MPI_SCAN_T
       
      END MODULE MPI__type_SV
