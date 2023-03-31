! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
!     Auxiliary file for IBM P4 machines 
!     when using the PESSL library in parallel runs.
!
      SUBROUTINE CPUTIM (TIME)

      DOUBLE PRECISION TIME

      TIME = MCLOCK()*0.01D0
      END
!
!     Add these two routines here if you are using pessl.
!     (Thanks to Vladimir Timochevski and Fedwa El-Mellouhi)
!
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     $                     LLD, INFO )
*
*  -- ScaLAPACK tools routine (version 1.5) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  DESCINIT initializes the descriptor vector with the 8 input arguments
*  M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  DESC    (output) INTEGER array of dimension DLEN_.
*          The array descriptor of a distributed matrix to be set.
*
*  M       (global input) INTEGER
*          The number of rows in the distributed matrix. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns in the distributed matrix. N >= 0.
*
*  MB      (global input) INTEGER
*          The blocking factor used to distribute the rows of the
*          matrix. MB >= 1.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  IRSRC   (global input) INTEGER
*          The process row over which the first row of the matrix is
*          distributed. 0 <= IRSRC < NPROW.
*
*  ICSRC   (global input) INTEGER
*          The process column over which the first column of the
*          matrix is distributed. 0 <= ICSRC < NPCOL.
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation on the matrix. The context itself is global.
*
*  LLD     (local input)  INTEGER
*          The leading dimension of the local array storing the local
*          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Note
*  ====
*
*  If the routine can recover from an erroneous input argument, it will
*  return an acceptable descriptor vector.  For example, if LLD = 0 on
*  input, DESC(LLD_) will contain the smallest leading dimension
*  required to store the specified M-by-N distributed matrix, INFO
*  will be set  -9 in that case.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, local_PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.MAX( 1, NUMROC( M, MB, MYROW, IRSRC,
     $                                NPROW ) ) ) THEN
         INFO = -9
      END IF
*
      IF( INFO.NE.0 )
     $   CALL local_PXERBLA( ICTXT, 'DESCINIT', -INFO )
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
      DESC( M_ )  = MAX( 0, M )
      DESC( N_ )  = MAX( 0, N )
      DESC( MB_ ) = MAX( 1, MB )
      DESC( NB_ ) = MAX( 1, NB )
      DESC( RSRC_ ) = MAX( 0, MIN( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = MAX( 0, MIN( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = MAX( LLD, MAX( 1, NUMROC( DESC( M_ ), DESC( MB_ ),
     $                              MYROW, DESC( RSRC_ ), NPROW ) ) )
*
      RETURN
*
*     End DESCINIT
*
      END
!
!     Local copy of pxerbla, called local_pxerbla
!     IBM's PESSL apparently does its own error handling, and does
!     not use the pxerbla model of the reference implementation.
!
      SUBROUTINE local_PXERBLA( ICTXT, SRNAME, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 20, 1995
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFO
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SRNAME
*     ..
*
*  Purpose
*  =======
*
*  PXERBLA is an error handler for the ScaLAPACK routines. It is called
*  by a ScaLAPACK routine if an input parameter has an invalid value.
*  A message is printed.  Installers may consider modifying this routine
*  in order to call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  SRNAME  (global input) CHARACTER*(*)
*          The name of the routine which called PXERBLA.
*
*  INFO    (global input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
*
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number', I4, ' had an illegal value' )
*
*     End of PXERBLA
*
      END
