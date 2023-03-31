! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE MPI_INTERFACES
!
!
!       Michael Hennecke
!       A Fortran 90 interface to MPI version 1.1
!       RZ Uni Karlsruhe, Internal Report 63/96
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
!       Modified by Alberto Garcia, wdpgaara@lg.ehu.es                       
!
        USE MPI__INCLUDE

!       ... generic overloads for <choice> argument routines ...

        USE MPI__logical_V   
        USE MPI__logical_S
        USE MPI__character_V  
        USE MPI__character_S

        include "V_S.uses"

!       ... these are for two different-rank <choice> arguments ...

        USE MPI__logical_VS   
        USE MPI__logical_SV
        USE MPI__character_VS 
        USE MPI__character_SV

        include "VS.uses"

!       ... this is for *ALL* combinations of type/rank of SENDRECV ...
!       USE MPI__sendrecv

        IMPLICIT NONE
!
!       ... A.9  Fortran Bindings for Point-to-Point Communication  ...

        INTERFACE
          SUBROUTINE MPI_GET_COUNT(STATUS, DATATYPE, COUNT, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: COUNT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GET_COUNT
          
          SUBROUTINE MPI_WAIT(REQUEST, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(INOUT) :: REQUEST
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_WAIT
          
          SUBROUTINE MPI_TEST(REQUEST, FLAG, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(INOUT) :: REQUEST
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TEST
          
          SUBROUTINE MPI_REQUEST_FREE(REQUEST, IERROR)
            INTEGER, INTENT(INOUT) :: REQUEST
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_REQUEST_FREE
          
          SUBROUTINE MPI_WAITANY(                                       &
     &        COUNT, ARRAY_OF_REQUESTS, INDEX, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: INDEX
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_WAITANY
          
          SUBROUTINE MPI_TESTANY(                                       &
     &        COUNT, ARRAY_OF_REQUESTS, INDEX, FLAG, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: INDEX
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TESTANY

          SUBROUTINE MPI_WAITALL(                                       &
     &        COUNT, ARRAY_OF_REQUESTS, ARRAY_OF_STATUSES, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_WAITALL
          
          SUBROUTINE MPI_TESTALL(                                       &
     &        COUNT, ARRAY_OF_REQUESTS, FLAG, ARRAY_OF_STATUSES, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TESTALL
          
          SUBROUTINE MPI_WAITSOME(                                      &
     &        INCOUNT, ARRAY_OF_REQUESTS, OUTCOUNT, ARRAY_OF_INDICES,   &
     &        ARRAY_OF_STATUSES, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: INCOUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: OUTCOUNT
            INTEGER, INTENT(OUT) :: ARRAY_OF_INDICES(*)
            INTEGER, INTENT(OUT) :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_WAITSOME
          
          SUBROUTINE MPI_TESTSOME(                                      &
     &        INCOUNT, ARRAY_OF_REQUESTS, OUTCOUNT, ARRAY_OF_INDICES,   &
     &        ARRAY_OF_STATUSES, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: INCOUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: OUTCOUNT
            INTEGER, INTENT(OUT) :: ARRAY_OF_INDICES(*)
            INTEGER, INTENT(OUT) :: ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TESTSOME
          
          SUBROUTINE MPI_IPROBE(SOURCE, TAG, COMM, FLAG, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: SOURCE
            INTEGER, INTENT(IN)  :: TAG
            INTEGER, INTENT(IN)  :: COMM
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_IPROBE
          
          SUBROUTINE MPI_PROBE(SOURCE, TAG, COMM, STATUS, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: SOURCE
            INTEGER, INTENT(IN)  :: TAG
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_PROBE
          
          SUBROUTINE MPI_CANCEL(REQUEST, IERROR)
            INTEGER, INTENT(IN)  :: REQUEST
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CANCEL
          
          SUBROUTINE MPI_TEST_CANCELLED(STATUS, FLAG, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: STATUS(MPI_STATUS_SIZE)
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TEST_CANCELLED
          
          SUBROUTINE MPI_START(REQUEST, IERROR)
            INTEGER, INTENT(INOUT) :: REQUEST
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_START
          
          SUBROUTINE MPI_STARTALL(COUNT, ARRAY_OF_REQUESTS, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(INOUT) :: ARRAY_OF_REQUESTS(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_STARTALL
          
          SUBROUTINE MPI_TYPE_CONTIGUOUS(                               &
     &        COUNT, OLDTYPE, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: OLDTYPE
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_CONTIGUOUS
          
          SUBROUTINE MPI_TYPE_VECTOR(                                   &
     &        COUNT, BLOCKLENGTH, STRIDE, OLDTYPE, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: BLOCKLENGTH
            INTEGER, INTENT(IN)  :: STRIDE
            INTEGER, INTENT(IN)  :: OLDTYPE
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_VECTOR
          
          SUBROUTINE MPI_TYPE_HVECTOR(                                  &
     &        COUNT, BLOCKLENGTH, STRIDE, OLDTYPE, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: BLOCKLENGTH
            INTEGER, INTENT(IN)  :: STRIDE
            INTEGER, INTENT(IN)  :: OLDTYPE
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_HVECTOR
          
          SUBROUTINE MPI_TYPE_INDEXED(                                  &
     &        COUNT, ARRAY_OF_BLOCKLENGTHS, ARRAY_OF_DISPLACEMENTS,     &
     &        OLDTYPE, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: ARRAY_OF_BLOCKLENGTHS(*)
            INTEGER, INTENT(IN)  :: ARRAY_OF_DISPLACEMENTS(*)
            INTEGER, INTENT(IN)  :: OLDTYPE
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_INDEXED
          
          SUBROUTINE MPI_TYPE_HINDEXED(                                 &
     &        COUNT, ARRAY_OF_BLOCKLENGTHS, ARRAY_OF_DISPLACEMENTS,     &
     &        OLDTYPE, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: ARRAY_OF_BLOCKLENGTHS(*)
            INTEGER, INTENT(IN)  :: ARRAY_OF_DISPLACEMENTS(*)
            INTEGER, INTENT(IN)  :: OLDTYPE
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_HINDEXED
          
          SUBROUTINE MPI_TYPE_STRUCT(                                   &
     &        COUNT, ARRAY_OF_BLOCKLENGTHS, ARRAY_OF_DISPLACEMENTS,     &
     &        ARRAY_OF_TYPES, NEWTYPE, IERROR)
            INTEGER, INTENT(IN)  :: COUNT
            INTEGER, INTENT(IN)  :: ARRAY_OF_BLOCKLENGTHS(*)
            INTEGER, INTENT(IN)  :: ARRAY_OF_DISPLACEMENTS(*)
            INTEGER, INTENT(IN)  :: ARRAY_OF_TYPES(*)
            INTEGER, INTENT(OUT) :: NEWTYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_STRUCT
          
          SUBROUTINE MPI_TYPE_EXTENT(DATATYPE, EXTENT, IERROR)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: EXTENT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_EXTENT
          
          SUBROUTINE MPI_TYPE_SIZE(DATATYPE, SIZE, IERROR)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: SIZE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_SIZE
          
          SUBROUTINE MPI_TYPE_COUNT(DATATYPE, COUNT, IERROR)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: COUNT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_COUNT
          
          SUBROUTINE MPI_TYPE_LB( DATATYPE, DISPLACEMENT, IERROR)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: DISPLACEMENT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_LB
          
          SUBROUTINE MPI_TYPE_UB( DATATYPE, DISPLACEMENT, IERROR)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: DISPLACEMENT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_UB
          
          SUBROUTINE MPI_TYPE_COMMIT(DATATYPE, IERROR)
            INTEGER, INTENT(INOUT) :: DATATYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_COMMIT
          
          SUBROUTINE MPI_TYPE_FREE(DATATYPE, IERROR)
            INTEGER, INTENT(INOUT) :: DATATYPE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TYPE_FREE
          
          SUBROUTINE MPI_GET_ELEMENTS(STATUS, DATATYPE, COUNT, IERROR)
            USE MPI__INCLUDE, ONLY: MPI_STATUS_SIZE
            INTEGER, INTENT(IN)  :: STATUS(MPI_STATUS_SIZE)
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(OUT) :: COUNT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GET_ELEMENTS
          
          SUBROUTINE MPI_PACK_SIZE(                                     &
     &        INCOUNT, DATATYPE, COMM, SIZE, IERROR)
            INTEGER, INTENT(IN)  :: INCOUNT
            INTEGER, INTENT(IN)  :: DATATYPE
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: SIZE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_PACK_SIZE
        END INTERFACE
          
!       ... A.10  Fortran Bindings for Collective Communication  ...
          
        INTERFACE
          SUBROUTINE MPI_BARRIER(COMM, IERROR) 
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_BARRIER
          
          SUBROUTINE MPI_OP_CREATE( FUNCTION, COMMUTE, OP, IERROR) 
            EXTERNAL FUNCTION 
            LOGICAL, INTENT(IN)  :: COMMUTE 
            INTEGER, INTENT(OUT) :: OP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_OP_CREATE
          
          SUBROUTINE MPI_OP_FREE( OP, IERROR) 
            INTEGER, INTENT(IN)  :: OP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_OP_FREE
        END INTERFACE

!       ... A.11  Fortran Bindings for Groups, Contexts, etc.  ...

        INTERFACE
          SUBROUTINE MPI_GROUP_SIZE(GROUP, SIZE, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(OUT) :: SIZE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_SIZE
          
          SUBROUTINE MPI_GROUP_RANK(GROUP, RANK, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(OUT) :: RANK
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_RANK
          
          SUBROUTINE MPI_GROUP_TRANSLATE_RANKS(                         &
     &        GROUP1, N, RANKS1, GROUP2, RANKS2, IERROR)
            INTEGER, INTENT(IN)  :: GROUP1
            INTEGER, INTENT(IN)  :: N
            INTEGER, INTENT(IN)  :: RANKS1(*)
            INTEGER, INTENT(IN)  :: GROUP2
            INTEGER, INTENT(OUT) :: RANKS2(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_TRANSLATE_RANKS
          
          SUBROUTINE MPI_GROUP_COMPARE(GROUP1, GROUP2, RESULT, IERROR)
            INTEGER, INTENT(IN)  :: GROUP1
            INTEGER, INTENT(IN)  :: GROUP2
            INTEGER, INTENT(OUT) :: RESULT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_COMPARE
          
          SUBROUTINE MPI_COMM_GROUP(COMM, GROUP, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: GROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_GROUP
          
          SUBROUTINE MPI_GROUP_UNION(GROUP1, GROUP2, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP1
            INTEGER, INTENT(IN)  :: GROUP2
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_UNION
          
          SUBROUTINE MPI_GROUP_INTERSECTION(                            &
     &        GROUP1, GROUP2, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP1
            INTEGER, INTENT(IN)  :: GROUP2
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_INTERSECTION
          
          SUBROUTINE MPI_GROUP_DIFFERENCE(                              &
     &        GROUP1, GROUP2, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP1
            INTEGER, INTENT(IN)  :: GROUP2
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_DIFFERENCE
          
          SUBROUTINE MPI_GROUP_INCL(GROUP, N, RANKS, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(IN)  :: N
            INTEGER, INTENT(IN)  :: RANKS(*)
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_INCL
          
          SUBROUTINE MPI_GROUP_EXCL(GROUP, N, RANKS, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(IN)  :: N
            INTEGER, INTENT(IN)  :: RANKS(*)
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_EXCL
          
          SUBROUTINE MPI_GROUP_RANGE_INCL(                              &
     &        GROUP, N, RANGES, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(IN)  :: N
            INTEGER, INTENT(IN)  :: RANGES(3,*)
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_RANGE_INCL
          
          SUBROUTINE MPI_GROUP_RANGE_EXCL(                              &
     &        GROUP, N, RANGES, NEWGROUP, IERROR)
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(IN)  :: N
            INTEGER, INTENT(IN)  :: RANGES(3,*)
            INTEGER, INTENT(OUT) :: NEWGROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_RANGE_EXCL
          
          SUBROUTINE MPI_GROUP_FREE(GROUP, IERROR)
            INTEGER, INTENT(INOUT) :: GROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GROUP_FREE
          
          SUBROUTINE MPI_COMM_SIZE(COMM, SIZE, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: SIZE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_SIZE
          
          SUBROUTINE MPI_COMM_RANK(COMM, RANK, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: RANK
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_RANK
          
          SUBROUTINE MPI_COMM_COMPARE(COMM1, COMM2, RESULT, IERROR)
            INTEGER, INTENT(IN)  :: COMM1
            INTEGER, INTENT(IN)  :: COMM2
            INTEGER, INTENT(OUT) :: RESULT
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_COMPARE
          
          SUBROUTINE MPI_COMM_DUP(COMM, NEWCOMM, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: NEWCOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_DUP
          
          SUBROUTINE MPI_COMM_CREATE(COMM, GROUP, NEWCOMM, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: GROUP
            INTEGER, INTENT(OUT) :: NEWCOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_CREATE
          
          SUBROUTINE MPI_COMM_SPLIT(COMM, COLOR, KEY, NEWCOMM, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: COLOR
            INTEGER, INTENT(IN)  :: KEY
            INTEGER, INTENT(OUT) :: NEWCOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_SPLIT
          
          SUBROUTINE MPI_COMM_FREE(COMM, IERROR)
            INTEGER, INTENT(INOUT) :: COMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_FREE
          
          SUBROUTINE MPI_COMM_TEST_INTER(COMM, FLAG, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: IERROR
            LOGICAL :: FLAG 
          END SUBROUTINE MPI_COMM_TEST_INTER
          
          SUBROUTINE MPI_COMM_REMOTE_SIZE(COMM, SIZE, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: SIZE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_REMOTE_SIZE
          
          SUBROUTINE MPI_COMM_REMOTE_GROUP(COMM, GROUP, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: GROUP
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_COMM_REMOTE_GROUP
          
          SUBROUTINE MPI_INTERCOMM_CREATE(                              &
     &        LOCAL_COMM, LOCAL_LEADER, PEER_COMM, REMOTE_LEADER,       &
     &        TAG, NEWINTERCOMM, IERROR)
            INTEGER, INTENT(IN)  :: LOCAL_COMM
            INTEGER, INTENT(IN)  :: LOCAL_LEADER
            INTEGER, INTENT(IN)  :: PEER_COMM
            INTEGER, INTENT(IN)  :: REMOTE_LEADER
            INTEGER, INTENT(IN)  :: TAG
            INTEGER, INTENT(OUT) :: NEWINTERCOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_INTERCOMM_CREATE
          
          SUBROUTINE MPI_INTERCOMM_MERGE(                               &
     &        INTERCOMM, HIGH, NEWINTRACOMM, IERROR)
            INTEGER, INTENT(IN)  :: INTERCOMM
            LOGICAL, INTENT(IN)  :: HIGH 
            INTEGER, INTENT(OUT) :: NEWINTRACOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_INTERCOMM_MERGE
          
          SUBROUTINE MPI_KEYVAL_CREATE(                                 &
     &        COPY_FN, DELETE_FN, KEYVAL, EXTRA_STATE, IERROR)
            EXTERNAL COPY_FN, DELETE_FN 
            INTEGER, INTENT(OUT) :: KEYVAL
            INTEGER, INTENT(IN)  :: EXTRA_STATE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_KEYVAL_CREATE
          
          SUBROUTINE MPI_KEYVAL_FREE(KEYVAL, IERROR)
            INTEGER, INTENT(INOUT) :: KEYVAL
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_KEYVAL_FREE
          
          SUBROUTINE MPI_ATTR_PUT(COMM, KEYVAL, ATTRIBUTE_VAL, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: KEYVAL
            INTEGER, INTENT(IN)  :: ATTRIBUTE_VAL
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ATTR_PUT
          
          SUBROUTINE MPI_ATTR_GET(                                      &
     &        COMM, KEYVAL, ATTRIBUTE_VAL, FLAG, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: KEYVAL
            INTEGER, INTENT(OUT) :: ATTRIBUTE_VAL
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ATTR_GET
          
          SUBROUTINE MPI_ATTR_DELETE(COMM, KEYVAL, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: KEYVAL
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ATTR_DELETE
        END INTERFACE
          
!       ... A.12  Fortran Bindings for Process Topologies  ...

        INTERFACE
          SUBROUTINE MPI_CART_CREATE(                                   &
     &        COMM_OLD, NDIMS, DIMS, PERIODS, REORDER, COMM_CART,       &
     &        IERROR)
            INTEGER, INTENT(IN)  :: COMM_OLD
            INTEGER, INTENT(IN)  :: NDIMS
            INTEGER, INTENT(IN)  :: DIMS(*)
            LOGICAL, INTENT(IN)  :: PERIODS(*)
            LOGICAL, INTENT(IN)  :: REORDER 
            INTEGER, INTENT(OUT) :: COMM_CART
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_CREATE
          
          SUBROUTINE MPI_DIMS_CREATE(NNODES, NDIMS, DIMS, IERROR)
            INTEGER, INTENT(IN)  :: NNODES
            INTEGER, INTENT(IN)  :: NDIMS
            INTEGER, INTENT(INOUT) :: DIMS(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_DIMS_CREATE
          
          SUBROUTINE MPI_GRAPH_CREATE(                                  &
     &        COMM_OLD, NNODES, INDEX, EDGES, REORDER, COMM_GRAPH,      &
     &        IERROR)
            INTEGER, INTENT(IN)  :: COMM_OLD
            INTEGER, INTENT(IN)  :: NNODES
            INTEGER, INTENT(IN)  :: INDEX(*)
            INTEGER, INTENT(IN)  :: EDGES(*)
            LOGICAL, INTENT(IN)  :: REORDER 
            INTEGER, INTENT(OUT) :: COMM_GRAPH
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPH_CREATE

          SUBROUTINE MPI_TOPO_TEST(COMM, STATUS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: STATUS
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_TOPO_TEST
          
          SUBROUTINE MPI_GRAPHDIMS_GET(COMM, NNODES, NEDGES, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: NNODES
            INTEGER, INTENT(OUT) :: NEDGES
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPHDIMS_GET
          
          SUBROUTINE MPI_GRAPH_GET(                                     &
     &        COMM, MAXINDEX, MAXEDGES, INDEX, EDGES, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: MAXINDEX
            INTEGER, INTENT(IN)  :: MAXEDGES
            INTEGER, INTENT(OUT) :: INDEX(*)
            INTEGER, INTENT(OUT) :: EDGES(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPH_GET
          
          SUBROUTINE MPI_CARTDIM_GET(COMM, NDIMS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: NDIMS
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CARTDIM_GET
          
          SUBROUTINE MPI_CART_GET(                                      &
     &        COMM, MAXDIMS, DIMS, PERIODS, COORDS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: MAXDIMS
            INTEGER, INTENT(OUT) :: DIMS(*)
            LOGICAL, INTENT(OUT) :: PERIODS(*) 
            INTEGER, INTENT(OUT) :: COORDS(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_GET
          
          SUBROUTINE MPI_CART_RANK(COMM, COORDS, RANK, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: COORDS(*)
            INTEGER, INTENT(OUT) :: RANK
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_RANK
          
          SUBROUTINE MPI_CART_COORDS(                                   &
     &        COMM, RANK, MAXDIMS, COORDS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: RANK
            INTEGER, INTENT(IN)  :: MAXDIMS
            INTEGER, INTENT(OUT) :: COORDS(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_COORDS
          
          SUBROUTINE MPI_GRAPH_NEIGHBORS_COUNT(                         &
     &        COMM, RANK, NNEIGHBORS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: RANK
            INTEGER, INTENT(OUT) :: NNEIGHBORS
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPH_NEIGHBORS_COUNT
          
          SUBROUTINE MPI_GRAPH_NEIGHBORS(                               &
     &        COMM, RANK, MAXNEIGHBORS, NEIGHBORS, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: RANK
            INTEGER, INTENT(IN)  :: MAXNEIGHBORS
            INTEGER, INTENT(OUT) :: NEIGHBORS(*)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPH_NEIGHBORS
          
          SUBROUTINE MPI_CART_SHIFT(                                    &
     &        COMM, DIRECTION, DISP, RANK_SOURCE, RANK_DEST, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: DIRECTION
            INTEGER, INTENT(IN)  :: DISP
            INTEGER, INTENT(OUT) :: RANK_SOURCE
            INTEGER, INTENT(OUT) :: RANK_DEST
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_SHIFT
          
          SUBROUTINE MPI_CART_SUB(COMM, REMAIN_DIMS, NEWCOMM, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            LOGICAL, INTENT(IN)  :: REMAIN_DIMS(*) 
            INTEGER, INTENT(OUT) :: NEWCOMM
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_SUB
          
          SUBROUTINE MPI_CART_MAP(                                      &
     &        COMM, NDIMS, DIMS, PERIODS, NEWRANK, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: NDIMS
            INTEGER, INTENT(IN)  :: DIMS(*)
            LOGICAL, INTENT(IN)  :: PERIODS(*) 
            INTEGER, INTENT(OUT) :: NEWRANK
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_CART_MAP
          
          SUBROUTINE MPI_GRAPH_MAP(                                     &
     &        COMM, NNODES, INDEX, EDGES, NEWRANK, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: NNODES
            INTEGER, INTENT(IN)  :: INDEX(*)
            INTEGER, INTENT(IN)  :: EDGES(*)
            INTEGER, INTENT(OUT) :: NEWRANK
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GRAPH_MAP
        END INTERFACE
          
!       ... A.13  Fortran Bindings for Environmental Inquiry  ...

        INTERFACE
          SUBROUTINE MPI_GET_PROCESSOR_NAME(NAME, RESULTLEN, IERROR)
            CHARACTER*(*), INTENT(OUT) :: NAME
            INTEGER, INTENT(OUT) :: RESULTLEN
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_GET_PROCESSOR_NAME
          
          SUBROUTINE MPI_ERRHANDLER_CREATE(FUNCTION, ERRHANDLER, IERROR)
            EXTERNAL FUNCTION 
            INTEGER, INTENT(OUT) :: ERRHANDLER
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERRHANDLER_CREATE
          
          SUBROUTINE MPI_ERRHANDLER_SET(COMM, ERRHANDLER, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: ERRHANDLER
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERRHANDLER_SET
          
          SUBROUTINE MPI_ERRHANDLER_GET(COMM, ERRHANDLER, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(OUT) :: ERRHANDLER
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERRHANDLER_GET
          
          SUBROUTINE MPI_ERRHANDLER_FREE(ERRHANDLER, IERROR)
            INTEGER, INTENT(IN)  :: ERRHANDLER  !  bug?
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERRHANDLER_FREE
          
          SUBROUTINE MPI_ERROR_STRING(                                  &
     &        ERRORCODE, STRING, RESULTLEN, IERROR)
            INTEGER, INTENT(IN)  :: ERRORCODE
            CHARACTER*(*), INTENT(OUT) :: STRING 
            INTEGER, INTENT(OUT) :: RESULTLEN
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERROR_STRING
          
          SUBROUTINE MPI_ERROR_CLASS(ERRORCODE, ERRORCLASS, IERROR)
            INTEGER, INTENT(IN)  :: ERRORCODE
            INTEGER, INTENT(OUT) :: ERRORCLASS
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ERROR_CLASS
          
!!! AG: Removed WTIME and WTICK specifications.
!!!     They are both real*8 in all implementations
!!!          
          SUBROUTINE MPI_INIT(IERROR)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_INIT
          
          SUBROUTINE MPI_FINALIZE(IERROR)
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_FINALIZE
          
          SUBROUTINE MPI_INITIALIZED(FLAG, IERROR)
            LOGICAL, INTENT(OUT) :: FLAG 
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_INITIALIZED
          
          SUBROUTINE MPI_ABORT(COMM, ERRORCODE, IERROR)
            INTEGER, INTENT(IN)  :: COMM
            INTEGER, INTENT(IN)  :: ERRORCODE
            INTEGER, INTENT(OUT) :: IERROR 
          END SUBROUTINE MPI_ABORT
        END INTERFACE
          
!       ... A.14  Fortran Bindings for Profiling  ...

        INTERFACE
          SUBROUTINE MPI_PCONTROL(LEVEL)
            INTEGER, INTENT(IN) :: LEVEL
!           ... maybe more arguments ...
          END SUBROUTINE MPI_PCONTROL
        END INTERFACE

        public

      END MODULE MPI_INTERFACES
