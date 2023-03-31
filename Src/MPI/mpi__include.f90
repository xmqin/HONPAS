      MODULE MPI__INCLUDE
!
!
!       Michael Hennecke
!       A Fortran 90 interface to MPI version 1.1
!       RZ Uni Karlsruhe, Internal Report 63/96
!
!         Auxiliary module MPI__INCLUDE
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
        IMPLICIT NONE

        INCLUDE 'mpif.h'

        public :: mpi_real
        public :: mpi_complex
        public :: mpi_double_complex
        public :: mpi_double_precision
        public :: mpi_2double_precision
        public :: mpi_integer, mpi_character, mpi_logical
        public :: mpi_integer8
        public :: mpi_packed
        public :: mpi_maxloc, mpi_sum, mpi_max, mpi_lor
        public :: mpi_comm_world
        public :: mpi_comm_self

        ! For threaded MPI
        public :: mpi_thread_single
        public :: mpi_thread_funneled


        ! Added by Toby White, <tow21@cam.ac.uk>; 24/03/2005
        ! All of these are in mpich-1.1 and should be visible.
        ! (See mpi-report-1.1/node169
        
        ! MPI return codes
        public :: MPI_SUCCESS
        public :: MPI_ERR_BUFFER
        public :: MPI_ERR_COUNT
        public :: MPI_ERR_TYPE
        public :: MPI_ERR_TAG
        public :: MPI_ERR_COMM 
        public :: MPI_ERR_RANK 
        public :: MPI_ERR_REQUEST
        public :: MPI_ERR_ROOT 
        public :: MPI_ERR_GROUP 
        public :: MPI_ERR_OP 
        public :: MPI_ERR_TOPOLOGY  
        public :: MPI_ERR_DIMS 
        public :: MPI_ERR_ARG 
        public :: MPI_ERR_UNKNOWN
        public :: MPI_ERR_TRUNCATE  
        public :: MPI_ERR_OTHER 
        public :: MPI_ERR_INTERN 
        public :: MPI_ERR_PENDING 
        public :: MPI_ERR_IN_STATUS  
        public :: MPI_ERR_LASTCODE
        
        ! assorted constants
        public :: MPI_BOTTOM
        public :: MPI_PROC_NULL
        public :: MPI_ANY_SOURCE
        public :: MPI_ANY_TAG
        public :: MPI_UNDEFINED
        public :: MPI_BSEND_OVERHEAD
        public :: MPI_KEYVAL_INVALID
        public :: MPI_IN_PLACE

        ! MPI_status
        public :: MPI_STATUS_SIZE
        public :: MPI_STATUSES_IGNORE
        public :: MPI_SOURCE
        public :: MPI_TAG
        public :: MPI_ERROR

        ! error handling specifiers
        public :: MPI_ERRORS_ARE_FATAL
        public :: MPI_ERRORS_RETURN
        
        ! Maximum sizes
        public :: MPI_MAX_PROCESSOR_NAME
        public :: MPI_MAX_ERROR_STRING

        ! Reserved communicators
        !public :: MPI_COMM_WORLD
        !public :: MPI_COMM_SELF

        ! datatypes for reduction functions
        public :: MPI_2REAL
        !public :: MPI_2DOUBLE_PRECISION
        public :: MPI_2INTEGER

        ! special datatypes for constructing derived datatypes
        public :: MPI_UB
        public :: MPI_LB

        ! results of communicator and group comparisons
        public :: MPI_IDENT
        public :: MPI_CONGRUENT
        public :: MPI_SIMILAR
        public :: MPI_UNEQUAL

        ! environmental inquiry keys
        public :: MPI_TAG_UB
        public :: MPI_IO
        public :: MPI_HOST
        public :: MPI_WTIME_IS_GLOBAL

        ! collective operations
        !public :: MPI_MAX
        public :: MPI_MIN
        !public :: MPI_SUM
        public :: MPI_PROD
        !public :: MPI_MAXLOC
        public :: MPI_MINLOC
        public :: MPI_BAND
        public :: MPI_BOR
        public :: MPI_BXOR
        public :: MPI_LAND
        !public :: MPI_LOR
        public :: MPI_LXOR

        ! Null handles
        public :: MPI_GROUP_NULL
        public :: MPI_COMM_NULL
        public :: MPI_DATATYPE_NULL
        public :: MPI_REQUEST_NULL
        public :: MPI_OP_NULL
        public :: MPI_ERRHANDLER_NULL

        ! Empty group
        public :: MPI_GROUP_EMPTY

        ! Topology types
        public :: MPI_GRAPH
        public :: MPI_CART

        private

      END MODULE MPI__INCLUDE
