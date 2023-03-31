      SUBROUTINE PDSYEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
     $                    DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 14, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYEVD computes  all the eigenvalues and eigenvectors
*  of a real symmetric matrix A by calling the recommended sequence
*  of ScaLAPACK routines.
*
*  In its present form, PDSYEVD assumes a homogeneous system and makes
*  no checks for consistency of the eigenvalues or eigenvectors across
*  the different processes.  Because of this, it is possible that a
*  heterogeneous system may return incorrect results without any error
*  messages.
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;     (NOT IMPLEMENTED YET)
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
*          global dimension (N, N), local dimension ( LLD_A,
*          LOCc(JA+N-1) )
*          On entry, the symmetric matrix A.  If UPLO = 'U', only the
*          upper triangular part of A is used to define the elements of
*          the symmetric matrix.  If UPLO = 'L', only the lower
*          triangular part of A is used to define the elements of the
*          symmetric matrix.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          If INFO=0, the eigenvalues in ascending order.
*
*  Z       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          Z contains the orthonormal eigenvectors
*          of the symmetric matrix A.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace required.
*
*  LWORK   (local input) INTEGER
*          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
*          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
*          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK array.  The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = 7*N + 8*NPCOL + 2
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The algorithm failed to compute the INFO/(N+1) th
*                eigenvalue while working on the submatrix lying in
*                global rows and columns mod(INFO,N+1).
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices sub( A ), sub( Z ) must verify
*  some alignment properties, namely the following expression
*  should be true:
*  ( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
*    IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
*    with IROFFA = MOD( IA-1, MB_A )
*     and ICOFFA = MOD( JA-1, NB_A ).
*
*  Further Details
*  ======= =======
*
*  Contributed by Francoise Tisseur, University of Manchester.
*
*  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
*              Conquer Algorithm for the Symmetric Eigenvalue Problem
*              on Distributed Memory Architectures",
*              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
*              (see also LAPACK Working Note 132)
*                http://www.netlib.org/lapack/lawns/lawn132.ps
*
*  =====================================================================
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IACOL, IAROW, ICOFFA, ICOFFZ, ICTXT, IINFO,
     $                   INDD, INDE, INDE2, INDTAU, INDWORK, INDWORK2,
     $                   IROFFA, IROFFZ, ISCALE, LIWMIN, LLWORK,
     $                   LLWORK2, LWMIN, MYCOL, MYROW, NB, NP, NPCOL,
     $                   NPROW, NQ, OFFSET, TRILWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
*     ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANSY
      EXTERNAL           LSAME, INDXG2P, NUMROC, PDLAMCH, PDLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DSCAL, PCHK1MAT,
     $                   PDLARED1D, PDLASCL, PDLASET, PDORMTR, PDSTEDC,
     $                   PDSYTRD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Quick return
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Test the input arguments.
*
      CALL BLACS_GRIDINFO( DESCZ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            NB = DESCA( NB_ )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFZ = MOD( IZ-1, DESCZ( MB_ ) )
            ICOFFZ = MOD( JZ-1, DESCZ( NB_ ) )
            IAROW = INDXG2P( IA, NB, MYROW, DESCA( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
            LQUERY = ( LWORK.EQ.-1 )
            TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
            LWMIN = MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DBLE( LWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( IROFFA.NE.IROFFZ .OR. ICOFFA.NE.ICOFFZ ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
               INFO = -( 1200+CSRC_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set up pointers into the WORK array
*
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDE2 = INDD + N
      INDWORK = INDE2 + N
      LLWORK = LWORK - INDWORK + 1
      INDWORK2 = INDD
      LLWORK2 = LWORK - INDWORK2 + 1
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
      ANRM = PDLANSY( 'M', UPLO, N, A, IA, JA, DESCA, WORK( INDWORK ) )
*
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PDLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
*
*     Reduce symmetric matrix to tridiagonal form.
*
*
      CALL PDSYTRD( UPLO, N, A, IA, JA, DESCA, WORK( INDD ),
     $              WORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     Copy the values of D, E to all processes.
*
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDD ), W,
     $                WORK( INDWORK ), LLWORK )
*
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDE2 ), WORK( INDE ),
     $                WORK( INDWORK ), LLWORK )
*
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
*
      IF( UPPER ) THEN
         OFFSET = 1
      ELSE
         OFFSET = 0
      END IF
      CALL PDSTEDC( 'I', N, W, WORK( INDE+OFFSET ), Z, IZ, JZ, DESCZ,
     $              WORK( INDWORK2 ), LLWORK2, IWORK, LIWORK, INFO )
*
      CALL PDORMTR( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK2 ),
     $              LLWORK2, IINFO )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of PDSYEVD
*
      END
      SUBROUTINE PDSYNGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                     DESCB, SCALE, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, IBTYPE, INFO, JA, JB, LWORK, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*
*  =======
*
*  PDSYNGST reduces a complex Hermitian-definite generalized
*  eigenproblem to standard form.
*
*  PDSYNGST performs the same function as PDHEGST, but is based on
*  rank 2K updates, which are faster and more scalable than
*  triangular solves (the basis of PDSYNGST).
*
*  PDSYNGST calls PDHEGST when UPLO='U', hence PDHENGST provides
*  improved performance only when UPLO='L', IBTYPE=1.
*
*  PDSYNGST also calls PDHEGST when insufficient workspace is
*  provided,  hence PDSYNGST provides improved
*  performance only when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*  In the following sub( A ) denotes A( IA:IA+N-1, JA:JA+N-1 ) and
*  sub( B ) denotes B( IB:IB+N-1, JB:JB+N-1 ).
*
*  If IBTYPE = 1, the problem is sub( A )*x = lambda*sub( B )*x,
*  and sub( A ) is overwritten by inv(U**H)*sub( A )*inv(U) or
*  inv(L)*sub( A )*inv(L**H)
*
*  If IBTYPE = 2 or 3, the problem is sub( A )*sub( B )*x = lambda*x or
*  sub( B )*sub( A )*x = lambda*x, and sub( A ) is overwritten by
*  U*sub( A )*U**H or L**H*sub( A )*L.
*
*  sub( B ) must have been previously factorized as U**H*U or L*L**H by
*  PDPOTRF.
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
*  IBTYPE   (global input) INTEGER
*          = 1: compute inv(U**H)*sub( A )*inv(U) or
*               inv(L)*sub( A )*inv(L**H);
*          = 2 or 3: compute U*sub( A )*U**H or L**H*sub( A )*L.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored and sub( B ) is
*                  factored as U**H*U;
*          = 'L':  Lower triangle of sub( A ) is stored and sub( B ) is
*                  factored as L*L**H.
*
*  N       (global input) INTEGER
*          The order of the matrices sub( A ) and sub( B ).  N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
*          the leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix, and its strictly
*          lower triangular part is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains
*          the lower triangular part of the matrix, and its strictly
*          upper triangular part is not referenced.
*
*          On exit, if INFO = 0, the transformed matrix, stored in the
*          same format as sub( A ).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)). On entry,
*          this array contains the local pieces of the triangular factor
*          from the Cholesky factorization of sub( B ), as returned by
*          PDPOTRF.
*
*  IB      (global input) INTEGER
*          B's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JB      (global input) INTEGER
*          B's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  SCALE   (global output) DOUBLE PRECISION
*          Amount by which the eigenvalues should be scaled to
*          compensate for the scaling performed in this routine.
*          At present, SCALE is always returned as 1.0, it is
*          returned here to allow for future enhancement.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                  dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( NB * ( NP0 +1 ), 3 * NB )
*
*          When IBTYPE = 1 and UPLO = 'L', PDSYNGST provides improved
*          performance when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*          where NB = MB_A = NB_A,
*          NP0 = NUMROC( N, NB, 0, 0, NPROW ),
*          NQ0 = NUMROC( N, NB, 0, 0, NPROW ),
*
*          NUMROC ia a ScaLAPACK tool functions
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the
*          optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONEHALF, ONE, MONE
      PARAMETER          ( ONEHALF = 0.5D0, ONE = 1.0D0, MONE = -1.0D0 )
      INTEGER            DLEN_, CTXT_, MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( DLEN_ = 9, CTXT_ = 2, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, INDAA, INDG, INDR, INDRT, IROFFA,
     $                   IROFFB, J, K, KB, LWMIN, LWOPT, MYCOL, MYROW,
     $                   NB, NP0, NPCOL, NPK, NPROW, NQ0, POSTK
*     ..
*     .. Local Arrays ..
      INTEGER            DESCAA( DLEN_ ), DESCG( DLEN_ ),
     $                   DESCR( DLEN_ ), DESCRT( DLEN_ ), IDUM1( 2 ),
     $                   IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           LSAME, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   PDGEMM, PDLACPY, PDSYGST, PDSYMM, PDSYR2K,
     $                   PDTRSM, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      SCALE = 1.0D0
*
      NB = DESCA( MB_ )
*
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $              NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $              NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $              NPCOL )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            NP0 = NUMROC( N, NB, 0, 0, NPROW )
            NQ0 = NUMROC( N, NB, 0, 0, NPCOL )
            LWMIN = MAX( NB*( NP0+1 ), 3*NB )
            IF( IBTYPE.EQ.1 .AND. .NOT.UPPER ) THEN
               LWOPT = 2*NP0*NB + NQ0*NB + NB*NB
            ELSE
               LWOPT = LWMIN
            END IF
            WORK( 1 ) = DBLE( LWOPT )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IBTYPE.LT.1 .OR. IBTYPE.GT.3 ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( N.LT.0 ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -9
            ELSE IF( ICOFFB.NE.0 .OR. IBCOL.NE.IACOL ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( MB_ ) ) THEN
               INFO = -( 1100+MB_ )
            ELSE IF( DESCB( NB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 1100+NB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1100+CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
            END IF
         END IF
         IDUM1( 1 ) = IBTYPE
         IDUM2( 1 ) = 1
         IF( UPPER ) THEN
            IDUM1( 2 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'L' )
         END IF
         IDUM2( 2 ) = 2
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IB,
     $                  JB, DESCB, 11, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYNGST', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
      IF( IBTYPE.NE.1 .OR. UPPER .OR. LWORK.LT.LWOPT ) THEN
         CALL PDSYGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                 DESCB, SCALE, INFO )
         RETURN
      END IF
*
      CALL DESCSET( DESCG, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCR, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCRT, NB, N, NB, NB, IAROW, IACOL, ICTXT, NB )
      CALL DESCSET( DESCAA, NB, NB, NB, NB, IAROW, IACOL, ICTXT, NB )
*
      INDG = 1
      INDR = INDG + DESCG( LLD_ )*NB
      INDAA = INDR + DESCR( LLD_ )*NB
      INDRT = INDAA + DESCAA( LLD_ )*NB
*
      DO 30 K = 1, N, NB
*
         KB = MIN( N-K+1, NB )
         POSTK = K + KB
         NPK = N - POSTK + 1
*
*
         CALL PDLACPY( 'A', N-POSTK+1, KB, B, POSTK+IB-1, K+JB-1, DESCB,
     $                 WORK( INDG ), POSTK, 1, DESCG )
         CALL PDLACPY( 'A', N-POSTK+1, KB, A, POSTK+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), POSTK, 1, DESCR )
         CALL PDLACPY( 'A', KB, K-1, A, K+IA-1, JA, DESCA,
     $                 WORK( INDRT ), 1, 1, DESCRT )
*
         CALL PDLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), K, 1, DESCR )
         CALL PDTRSM( 'Right', 'L', 'N', 'N', NPK, KB, MONE, B, K+IB-1,
     $                K+JB-1, DESCB, WORK( INDG ), POSTK, 1, DESCG )
*
         CALL PDSYMM( 'Right', 'L', NPK, KB, ONEHALF, A, K+IA-1, K+JA-1,
     $                DESCA, WORK( INDG ), POSTK, 1, DESCG, ONE,
     $                WORK( INDR ), POSTK, 1, DESCR )
*
         CALL PDSYR2K( 'Lower', 'No T', NPK, KB, ONE, WORK( INDG ),
     $                 POSTK, 1, DESCG, WORK( INDR ), POSTK, 1, DESCR,
     $                 ONE, A, POSTK+IA-1, POSTK+JA-1, DESCA )
*
         CALL PDGEMM( 'No T', 'No Conj', NPK, K-1, KB, ONE,
     $                WORK( INDG ), POSTK, 1, DESCG, WORK( INDRT ), 1,
     $                1, DESCRT, ONE, A, POSTK+IA-1, JA, DESCA )
*
         CALL PDSYMM( 'Right', 'L', NPK, KB, ONE, WORK( INDR ), K, 1,
     $                DESCR, WORK( INDG ), POSTK, 1, DESCG, ONE, A,
     $                POSTK+IA-1, K+JA-1, DESCA )
*
         CALL PDTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, K-1,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, K+IA-1, JA,
     $                DESCA )
*
         CALL PDLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDAA ), 1, 1, DESCAA )
*
         IF( MYROW.EQ.DESCAA( RSRC_ ) .AND. MYCOL.EQ.DESCAA( CSRC_ ) )
     $        THEN
            DO 20 I = 1, KB
               DO 10 J = 1, I
                  WORK( INDAA+J-1+( I-1 )*DESCAA( LLD_ ) )
     $               = WORK( INDAA+I-1+( J-1 )*DESCAA( LLD_ ) )
   10          CONTINUE
   20       CONTINUE
         END IF
*
         CALL PDTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1,
     $                1, DESCAA )
*
         CALL PDTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', KB, KB, ONE,
     $                B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1, 1,
     $                DESCAA )
*
         CALL PDLACPY( 'L', KB, KB, WORK( INDAA ), 1, 1, DESCAA, A,
     $                 K+IA-1, K+JA-1, DESCA )
*
         CALL PDTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', NPK, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, POSTK+IA-1,
     $                K+JA-1, DESCA )
*
         DESCR( CSRC_ ) = MOD( DESCR( CSRC_ )+1, NPCOL )
         DESCG( CSRC_ ) = MOD( DESCG( CSRC_ )+1, NPCOL )
         DESCRT( RSRC_ ) = MOD( DESCRT( RSRC_ )+1, NPROW )
         DESCAA( RSRC_ ) = MOD( DESCAA( RSRC_ )+1, NPROW )
         DESCAA( CSRC_ ) = MOD( DESCAA( CSRC_ )+1, NPCOL )
   30 CONTINUE
*
      WORK( 1 ) = DBLE( LWOPT )
*
      RETURN
      END
      INTEGER          FUNCTION PJLAENV( ICTXT, ISPEC, NAME, OPTS, N1,
     $                 N2, N3, N4 )
*
*  -- ScaLAPACK test routine (version 1.6.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ICTXT, ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*
*  =======
*
*  PJLAENV is called from the ScaLAPACK symmetric and Hermitian
*  tailored eigen-routines to choose
*  problem-dependent parameters for the local environment.  See ISPEC
*  for a description of the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (global input) INTEGER
*          Specifies the parameter to be returned as the value of
*          PJLAENV.
*          = 1: the data layout blocksize;
*          = 2: the panel blocking factor;
*          = 3: the algorithmic blocking factor;
*          = 4: execution path control;
*          = 5: maximum size for direct call to the LAPACK routine
*
*  NAME    (global input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (global input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (global input) INTEGER
*  N2      (global input) INTEGER
*  N3      (global input) INTEGER
*  N4      (global input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*          At present, only N1 is used, and it (N1) is used only for
*          'TTRD'
*
* (PJLAENV) (global or local output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if PJLAENV = -k, the k-th argument had an illegal
*          value.
*
*          Most parameters set via a call to PJLAENV must be identical
*          on all processors and hence PJLAENV will return the same
*          value to all procesors (i.e. global output).  However some,
*          in particular, the panel blocking factor can be different
*          on each processor and hence PJLAENV can return different
*          values on different processors (i.e. local output).
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling PJLAENV from
*  the ScaLAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by PJLAENV is checked for validity
*      in the calling subroutine.  For example, PJLAENV is used to
*      retrieve the optimal blocksize for STRTRI as follows:
*
*      NB = PJLAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  PJLAENV is patterned after ILAENV and keeps the same interface in
*  anticipation of future needs, even though PJLAENV is only sparsely
*  used at present in ScaLAPACK.  Most ScaLAPACK codes use the input
*  data layout blocking factor as the algorithmic blocking factor -
*  hence there is no need or opportunity to set the algorithmic or
*  data decomposition blocking factor.
*
*  pXYYtevx.f and pXYYtgvx.f and pXYYttrd.f are the only codes which
*  call PJLAENV in this release.  pXYYtevx.f and pXYYtgvx.f redistribute
*  the data to the best data layout for each transformation.  pXYYttrd.f
*  uses a data layout blocking factor of 1 and a
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, GLOBAL, SNAME
      CHARACTER          C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*8        SUBNAM
      INTEGER            I, IC, IDUMM, IZ, LLD_, MB_, MSZ, M_, NB, NB_,
     $                   N_
      REAL               BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   RSRC_
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR
*     ..
*
*
*     .. External Subroutines ..
      EXTERNAL           IGAMX2D
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*
*
      GO TO ( 10, 10, 10, 10, 10 )ISPEC
*
*     Invalid value for ISPEC
*
      PJLAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      PJLAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.100 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 2: 2 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 3: 4 )
      C3 = SUBNAM( 5: 7 )
      C4 = C3( 2: 3 )
*
*     This is to keep ftnchek happy
*
      IF( ( N2+N3+N4 )*0.NE.0 ) THEN
         C4 = OPTS
         C3 = C4
      END IF
*
      GO TO ( 50, 60, 70, 80, 90 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  data layout block size
*     (global - all processes must use the same value)
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'LLT' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 1
            ELSE
               NB = 1
            END IF
         ELSE IF( C3.EQ.'GST' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BCK' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRS' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      END IF
*
*
      PJLAENV = NB
      GLOBAL = .TRUE.
      GO TO 100
*
   60 CONTINUE
*
*     ISPEC = 2:  panel blocking factor (Used only in PxyyTTRD)
*     (local - different processes may use different values)
*
      NB = 16
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      PJLAENV = NB
      GLOBAL = .FALSE.
      GO TO 100
*
*
   70 CONTINUE
*
*     ISPEC = 3:  algorithmic blocking factor (Used only in PxyyTTRD)
*     (global - all processes must use the same value)
*
      NB = 1
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 16
            ELSE
               NB = 16
            END IF
         END IF
      END IF
      PJLAENV = NB
      GLOBAL = .TRUE.
      GO TO 100
*
   80 CONTINUE
*
*     ISPEC = 4:  Execution path options (Used only in PxyyTTRD)
*     (global - all processes must use the same value)
*
      PJLAENV = -4
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
*           V and H interleaved (default is not interleaved)
            IF( N1.EQ.1 ) THEN
               PJLAENV = 1
            END IF
*
*           Two ZGEMMs (default is one ZGEMM)
            IF( N1.EQ.2 ) THEN
               PJLAENV = 0
            END IF
*           Balanced Update (default is minimum communication update)
            IF( N1.EQ.3 ) THEN
               PJLAENV = 0
            END IF
         END IF
      END IF
      GLOBAL = .TRUE.
      GO TO 100
*
   90 CONTINUE
*
*     ISPEC = 5:  Minimum size to justify call to parallel code
*     (global - all processes must use the same value)
*
      MSZ = 0
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               MSZ = 100
            ELSE
               MSZ = 100
            END IF
         END IF
      END IF
      PJLAENV = MSZ
      GLOBAL = .TRUE.
      GO TO 100
*
  100 CONTINUE
*
      IF( GLOBAL ) THEN
         CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, PJLAENV, 1, IDUMM,
     $                 IDUMM, -1, -1, IDUMM )
      END IF
*
*
*
      RETURN
*
*     End of PJLAENV
*
      END
      SUBROUTINE PDSTEDC( COMPZ, N, D, E, Q, IQ, JQ, DESCQ, WORK, LWORK,
     $                    IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 13, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, IQ, JQ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*  PDSTEDC computes all eigenvalues and eigenvectors of a
*  symmetric tridiagonal matrix in parallel, using the divide and
*  conquer algorithm.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See DLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.    (NOT IMPLEMENTED YET)
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original dense symmetric
*                  matrix also.  On entry, Z contains the orthogonal
*                  matrix used to reduce the original matrix to
*                  tridiagonal form.            (NOT IMPLEMENTED YET)
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  D       (global input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in descending order.
*
*  E       (global input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Q       (local output) DOUBLE PRECISION array,
*          local dimension ( LLD_Q, LOCc(JQ+N-1))
*          Q  contains the orthonormal eigenvectors of the symmetric
*          tridiagonal matrix.
*          On output, Q is distributed across the P processes in block
*          cyclic format.
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace needed.
*
*  LWORK   (local input/output) INTEGER,
*          the dimension of the array WORK.
*          LWORK = 6*N + 2*NP*NQ
*          NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
*          NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
*
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK array.  The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = 2 + 7*N + 8*NPCOL
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The algorithm failed to compute the INFO/(N+1) th
*                eigenvalue while working on the submatrix lying in
*                global rows and columns mod(INFO,N+1).
*
*  Further Details
*  ======= =======
*
*  Contributed by Francoise Tisseur, University of Manchester.
*
*  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
*              Conquer Algorithm for the Symmetric Eigenvalue Problem
*              on Distributed Memory Architectures",
*              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
*              (see also LAPACK Working Note 132)
*                http://www.netlib.org/lapack/lawns/lawn132.ps
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICOFFQ, IIQ, IPQ, IQCOL, IQROW, IROFFQ, JJQ,
     $                   LDQ, LIWMIN, LWMIN, MYCOL, MYROW, NB, NP,
     $                   NPCOL, NPROW, NQ
      DOUBLE PRECISION   ORGNRM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   DLANST
      EXTERNAL           INDXG2P, LSAME, NUMROC, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DLASCL, DSTEDC,
     $                   INFOG2L, PDLAED0, PDLASRT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Test the input parameters.
*
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      LDQ = DESCQ( LLD_ )
      NB = DESCQ( NB_ )
      NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
      NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IQ, JQ, DESCQ, 8, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            IROFFQ = MOD( IQ-1, DESCQ( MB_ ) )
            ICOFFQ = MOD( JQ-1, DESCQ( NB_ ) )
            IQROW = INDXG2P( IQ, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            IQCOL = INDXG2P( JQ, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = 6*N + 2*NP*NQ
            LIWMIN = 2 + 7*N + 8*NPCOL
*
            WORK( 1 ) = DBLE( LWMIN )
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
            IF( .NOT.LSAME( COMPZ, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( IROFFQ.NE.ICOFFQ .OR. ICOFFQ.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -12
            END IF
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCQ( CTXT_ ), 'PDSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return
*
      IF( N.EQ.0 )
     $   GO TO 10
      CALL INFOG2L( IQ, JQ, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     $              IQROW, IQCOL )
      IF( N.EQ.1 ) THEN
         IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL )
     $      Q( 1 ) = ONE
         GO TO 10
      END IF
*
*     If N is smaller than the minimum divide size NB, then
*     solve the problem with the serial divide and conquer
*     code locally.
*
      IF( N.LE.NB ) THEN
         IF( ( MYROW.EQ.IQROW ) .AND. ( MYCOL.EQ.IQCOL ) ) THEN
            IPQ = IIQ + ( JJQ-1 )*LDQ
            CALL DSTEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( N+1 ) + N
               GO TO 10
            END IF
         END IF
         GO TO 10
      END IF
*
*     If P=NPROW*NPCOL=1, solve the problem with DSTEDC.
*
      IF( NPCOL*NPROW.EQ.1 ) THEN
         IPQ = IIQ + ( JJQ-1 )*LDQ
         CALL DSTEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK, IWORK,
     $                LIWORK, INFO )
         GO TO 10
      END IF
*
*     Scale matrix to allowable range, if necessary.
*
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.NE.ZERO ) THEN
         CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
         CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, N-1, 1, E, N-1, INFO )
      END IF
*
      CALL PDLAED0( N, D, E, Q, IQ, JQ, DESCQ, WORK, IWORK, INFO )
*
*     Sort eigenvalues and corresponding eigenvectors
*
      CALL PDLASRT( 'I', N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, IWORK,
     $              LIWORK, INFO )
*
*           Scale back.
*
      IF( ORGNRM.NE.ZERO )
     $   CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
*
   10 CONTINUE
*
      IF( LWORK.GT.0 )
     $   WORK( 1 ) = DBLE( LWMIN )
      IF( LIWORK.GT.0 )
     $   IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of PDSTEDC
*
      END
      SUBROUTINE PDLAED0( N, D, E, Q, IQ, JQ, DESCQ, WORK, IWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INFO, IQ, JQ, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAED0 computes all eigenvalues and corresponding eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  D       (global input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in descending order.
*
*  E       (global input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Q       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Q, LOCc(JQ+N-1))
*          Q  contains the orthonormal eigenvectors of the symmetric
*          tridiagonal matrix.
*          On output, Q is distributed across the P processes in block
*          cyclic format.
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*
*  WORK    (local workspace ) DOUBLE PRECISION array, dimension (LWORK)
*          LWORK = 6*N + 2*NP*NQ, with
*          NP = NUMROC( N, MB_Q, MYROW, IQROW, NPROW )
*          NQ = NUMROC( N, NB_Q, MYCOL, IQCOL, NPCOL )
*          IQROW = INDXG2P( IQ, NB_Q, MYROW, RSRC_Q, NPROW )
*          IQCOL = INDXG2P( JQ, MB_Q, MYCOL, CSRC_Q, NPCOL )
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          LIWORK = 2 + 7*N + 8*NPCOL
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The algorithm failed to compute the INFO/(N+1) th
*                eigenvalue while working on the submatrix lying in
*                global rows and columns mod(INFO,N+1).
*
*  =====================================================================
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ID, IDCOL, IDROW, IID, IINFO, IIQ, IM1, IM2,
     $                   IPQ, IQCOL, IQROW, J, JJD, JJQ, LDQ, MATSIZ,
     $                   MYCOL, MYROW, N1, NB, NBL, NBL1, NPCOL, NPROW,
     $                   SUBPBS, TSUBPBS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, DGERV2D,
     $                   DGESD2D, DSTEQR, INFOG2L, PDLAED1, PXERBLA
*     ..
*     .. External Functions ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Test the input parameters.
*
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( DESCQ( NB_ ).GT.N .OR. N.LT.2 )
     $   INFO = -1
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCQ( CTXT_ ), 'PDLAED0', -INFO )
         RETURN
      END IF
*
      NB = DESCQ( NB_ )
      LDQ = DESCQ( LLD_ )
      CALL INFOG2L( IQ, JQ, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     $              IQROW, IQCOL )
*
*     Determine the size and placement of the submatrices, and save in
*     the leading elements of IWORK.
*
      TSUBPBS = ( N-1 ) / NB + 1
      IWORK( 1 ) = TSUBPBS
      SUBPBS = 1
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.1 ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
*     Divide the matrix into TSUBPBS submatrices of size at most NB
*     using rank-1 modifications (cuts).
*
      DO 40 I = NB + 1, N, NB
         IM1 = I - 1
         D( IM1 ) = D( IM1 ) - ABS( E( IM1 ) )
         D( I ) = D( I ) - ABS( E( IM1 ) )
   40 CONTINUE
*
*     Solve each submatrix eigenproblem at the bottom of the divide and
*     conquer tree. D is the same on each process.
*
      DO 50 ID = 1, N, NB
         CALL INFOG2L( IQ-1+ID, JQ-1+ID, DESCQ, NPROW, NPCOL, MYROW,
     $                 MYCOL, IID, JJD, IDROW, IDCOL )
         MATSIZ = MIN( NB, N-ID+1 )
         IF( MYROW.EQ.IDROW .AND. MYCOL.EQ.IDCOL ) THEN
            IPQ = IID + ( JJD-1 )*LDQ
            CALL DSTEQR( 'I', MATSIZ, D( ID ), E( ID ), Q( IPQ ), LDQ,
     $                   WORK, INFO )
            IF( INFO.NE.0 ) THEN
               CALL PXERBLA( DESCQ( CTXT_ ), 'DSTEQR', -INFO )
               RETURN
            END IF
            IF( MYROW.NE.IQROW .OR. MYCOL.NE.IQCOL ) THEN
               CALL DGESD2D( DESCQ( CTXT_ ), MATSIZ, 1, D( ID ), MATSIZ,
     $                       IQROW, IQCOL )
            END IF
         ELSE IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL ) THEN
            CALL DGERV2D( DESCQ( CTXT_ ), MATSIZ, 1, D( ID ), MATSIZ,
     $                    IDROW, IDCOL )
         END IF
   50 CONTINUE
*
      IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL ) THEN
         CALL DGEBS2D( DESCQ( CTXT_ ), 'A', ' ', N, 1, D, N )
      ELSE
         CALL DGEBR2D( DESCQ( CTXT_ ), 'A', ' ', N, 1, D, N, IQROW,
     $                 IQCOL )
      END IF
*
*     Successively merge eigensystems of adjacent submatrices
*     into eigensystem for the corresponding larger matrix.
*
*     while ( SUBPBS > 1 )
*
   60 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         IM2 = SUBPBS - 2
         DO 80 I = 0, IM2, 2
            IF( I.EQ.0 ) THEN
               NBL = IWORK( 2 )
               NBL1 = IWORK( 1 )
               IF( NBL1.EQ.0 )
     $            GO TO 70
               ID = 1
               MATSIZ = MIN( N, NBL*NB )
               N1 = NBL1*NB
            ELSE
               NBL = IWORK( I+2 ) - IWORK( I )
               NBL1 = NBL / 2
               IF( NBL1.EQ.0 )
     $            GO TO 70
               ID = IWORK( I )*NB + 1
               MATSIZ = MIN( NB*NBL, N-ID+1 )
               N1 = NBL1*NB
            END IF
*
*     Merge lower order eigensystems (of size N1 and MATSIZ - N1)
*     into an eigensystem of size MATSIZ.
*
            CALL PDLAED1( MATSIZ, N1, D( ID ), ID, Q, IQ, JQ, DESCQ,
     $                    E( ID+N1-1 ), WORK, IWORK( SUBPBS+1 ), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = IINFO*( N+1 ) + ID
            END IF
*
   70       CONTINUE
            IWORK( I / 2+1 ) = IWORK( I+2 )
   80    CONTINUE
         SUBPBS = SUBPBS / 2
*
         GO TO 60
      END IF
*
*     end while
*
   90 CONTINUE
      RETURN
*
*     End of PDLAED0
*
      END
      SUBROUTINE PDLASRT( ID, N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, 
     $                    IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 22, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, IQ, JQ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLASRT Sort the numbers in D in increasing order and the
*  corresponding vectors in Q.
*
*  Arguments
*  =========
*
*  ID      (global input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order. (NOT IMPLEMENTED YET)
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( Q ). N >= 0.
*
*  D       (global input/output) DOUBLE PRECISION array, dimmension (N)
*          On exit, the number in D are sorted in increasing order.
*
*  Q       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_Q, LOCc(JQ+N-1) ). This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be copied from.
*
*  IQ      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( Q ).
*
*  JQ      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( Q ).
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*          dimension (LWORK)
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK = MAX( N, NP * ( NB + NQ ))
*          where
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW ),
*          NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                  dimension (LIWORK)
*
*  LIWORK (local or global input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = N + 2*NB + 2*NPCOL
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CL, COL, DUMMY, I, ICTXT, IID, IIQ, INDCOL,
     $                   INDX, INDXC, INDXG, IPQ, IPQ2, IPW, IPWORK, J,
     $                   JJQ, K, L, LDQ, LEND, LIWMIN, LWMIN, MYCOL,
     $                   MYROW, NB, ND, NP, NPCOL, NPROW, NQ, PSQ, QCOL,
     $                   QTOT, SBUF
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      EXTERNAL           INDXG2L, INDXG2P, LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA, DCOPY,
     $                   DGERV2D, DGESD2D, DLACPY, DLAPST
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IQ, JQ, DESCQ, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            LDQ = DESCQ( LLD_ )
            NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = MAX( N, NP*( NB+NQ ) )
            LIWMIN = N + 2*( NB+NPCOL )
            IF( .NOT.LSAME( ID, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN ) THEN
               INFO = -9
            ELSE IF( LIWORK.LT.LIWMIN ) THEN
               INFO = -11
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDLASRT', -INFO )
         RETURN
      END IF
*
*     Set Pointers
*
      INDXC = 1
      INDX = INDXC + N
      INDXG = INDX
      INDCOL = INDXG + NB
      QTOT = INDCOL + NB
      PSQ = QTOT + NPCOL
*
      IID = 1
      IPQ2 = 1
      IPW = IPQ2 + NP*NQ
*
      DUMMY = 0
      IIQ = INDXG2L( IQ, NB, DUMMY, DUMMY, NPROW )
*
*     Sort the eigenvalues in D
*
      CALL DLAPST( 'I', N, D, IWORK( INDX ), INFO )
*
      DO 10 L = 0, N - 1
         WORK( IID+L ) = D( IWORK( INDX+L ) )
         IWORK( INDXC-1+IWORK( INDX+L ) ) = IID + L
   10 CONTINUE
      CALL DCOPY( N, WORK, 1, D, 1 )
*
      ND = 0
   20 CONTINUE
      IF( ND.LT.N ) THEN
         LEND = MIN( NB, N-ND )
         J = JQ + ND
         QCOL = INDXG2P( J, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )
         K = 0
         DO 30 L = 0, LEND - 1
            I = JQ - 1 + IWORK( INDXC+ND+L )
            CL = INDXG2P( I, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )
            IWORK( INDCOL+L ) = CL
            IF( MYCOL.EQ.CL ) THEN
               IWORK( INDXG+K ) = IWORK( INDXC+ND+L )
               K = K + 1
            END IF
   30    CONTINUE
*
         IF( MYCOL.EQ.QCOL ) THEN
            DO 40 CL = 0, NPCOL - 1
               IWORK( QTOT+CL ) = 0
   40       CONTINUE
            DO 50 L = 0, LEND - 1
               IWORK( QTOT+IWORK( INDCOL+L ) ) = IWORK( QTOT+
     $            IWORK( INDCOL+L ) ) + 1
   50       CONTINUE
            IWORK( PSQ ) = 1
            DO 60 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
   60       CONTINUE
            DO 70 L = 0, LEND - 1
               CL = IWORK( INDCOL+L )
               I = JQ + ND + L
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IIQ + ( JJQ-1 )*LDQ
               IPWORK = IPW + ( IWORK( PSQ+CL )-1 )*NP
               CALL DCOPY( NP, Q( IPQ ), 1, WORK( IPWORK ), 1 )
               IWORK( PSQ+CL ) = IWORK( PSQ+CL ) + 1
   70       CONTINUE
            IWORK( PSQ ) = 1
            DO 80 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
   80       CONTINUE
            DO 90 L = 0, K - 1
               I = IWORK( INDXG+L )
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IPQ2 + ( JJQ-1 )*NP
               IPWORK = IPW + ( IWORK( PSQ+MYCOL )-1 )*NP
               CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
               IWORK( PSQ+MYCOL ) = IWORK( PSQ+MYCOL ) + 1
   90       CONTINUE
            DO 100 CL = 1, NPCOL - 1
               COL = MOD( MYCOL+CL, NPCOL )
               SBUF = IWORK( QTOT+COL )
               IF( SBUF.NE.0 ) THEN
                  IPWORK = IPW + ( IWORK( PSQ+COL )-1 )*NP
                  CALL DGESD2D( DESCQ( CTXT_ ), NP, SBUF,
     $                          WORK( IPWORK ), NP, MYROW, COL )
               END IF
  100       CONTINUE
*
         ELSE
*
            IF( K.NE.0 ) THEN
               CALL DGERV2D( DESCQ( CTXT_ ), NP, K, WORK( IPW ), NP,
     $                       MYROW, QCOL )
               DO 110 L = 0, K - 1
                  I = JQ - 1 + IWORK( INDXG+L )
                  JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
                  IPQ = 1 + ( JJQ-1 )*NP
                  IPWORK = IPW + L*NP
                  CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
  110          CONTINUE
            END IF
         END IF
         ND = ND + NB
         GO TO 20
      END IF
      CALL DLACPY( 'Full', NP, NQ, WORK, NP, Q( IIQ ), LDQ )
*
*     End of PDLASRT
*
      END
      SUBROUTINE PDLAED1( N, N1, D, ID, Q, IQ, JQ, DESCQ, RHO, WORK,
     $                    IWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ID, INFO, IQ, JQ, N, N1
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAED1 computes the updated eigensystem of a diagonal
*  matrix after modification by a rank-one symmetric matrix,
*  in parallel.
*
*    T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)
*
*     where Z = Q'u, u is a vector of length N with ones in the
*     N1 and N1 + 1 th elements and zeros elsewhere.
*
*     The eigenvectors of the original matrix are stored in Q, and the
*     eigenvalues are in D.  The algorithm consists of three stages:
*
*        The first stage consists of deflating the size of the problem
*        when there are multiple eigenvalues or if there is a zero in
*        the Z vector.  For each such occurence the dimension of the
*        secular equation problem is reduced by one.  This stage is
*        performed by the routine PDLAED2.
*
*        The second stage consists of calculating the updated
*        eigenvalues. This is done by finding the roots of the secular
*        equation via the routine SLAED4 (as called by PDLAED3).
*        This routine also calculates the eigenvectors of the current
*        problem.
*
*        The final stage consists of computing the updated eigenvectors
*        directly using the updated eigenvalues.  The eigenvectors for
*        the current problem are multiplied with the eigenvectors from
*        the overall problem.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*
*  N1      (input) INTEGER
*          The location of the last eigenvalue in the leading
*          sub-matrix.
*          min(1,N) <= N1 <= N.
*
*  D       (global input/output) DOUBLE PRECISION array, dimension (N)
*          On entry,the eigenvalues of the rank-1-perturbed matrix.
*          On exit, the eigenvalues of the repaired matrix.
*
*  ID      (global input) INTEGER
*          Q's global row/col index, which points to the beginning
*          of the submatrix which is to be operated on.
*
*  Q       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Q, LOCc(JQ+N-1))
*          Q  contains the orthonormal eigenvectors of the symmetric
*          tridiagonal matrix.
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  RHO    (input) DOUBLE PRECISION
*         The subdiagonal entry used to create the rank-1 modification.
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension 6*N + 2*NP*NQ
*
*  IWORK   (local workspace/output) INTEGER array,
*          dimension 7*N + 8*NPCOL + 2
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The algorithm failed to compute the ith eigenvalue.
*
*  =====================================================================
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, COLTYP, IBUF, ICTOT, ICTXT, IDLMDA, IIQ,
     $                   INDCOL, INDROW, INDX, INDXC, INDXP, INDXR, INQ,
     $                   IPQ, IPQ2, IPSM, IPU, IPWORK, IQ1, IQ2, IQCOL,
     $                   IQQ, IQROW, IW, IZ, J, JC, JJ2C, JJC, JJQ, JNQ,
     $                   K, LDQ, LDQ2, LDU, MYCOL, MYROW, NB, NN, NN1,
     $                   NN2, NP, NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCQ2( DLEN_ ), DESCU( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DESCINIT, INFOG1L,
     $                   INFOG2L, PDGEMM, PDLAED2, PDLAED3, PDLAEDZ,
     $                   PDLASET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*
*     Test the input parameters.
*
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ID.GT.DESCQ( N_ ) ) THEN
         INFO = -4
      ELSE IF( N1.GE.N ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCQ( CTXT_ ), 'PDLAED1', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     The following values are  integer pointers which indicate
*     the portion of the workspace used by a particular array
*     in PDLAED2 and PDLAED3.
*
      ICTXT = DESCQ( CTXT_ )
      NB = DESCQ( NB_ )
      LDQ = DESCQ( LLD_ )
*
      CALL INFOG2L( IQ-1+ID, JQ-1+ID, DESCQ, NPROW, NPCOL, MYROW, MYCOL,
     $              IIQ, JJQ, IQROW, IQCOL )
*
      NP = NUMROC( N, DESCQ( MB_ ), MYROW, IQROW, NPROW )
      NQ = NUMROC( N, DESCQ( NB_ ), MYCOL, IQCOL, NPCOL )
*
      LDQ2 = MAX( NP, 1 )
      LDU = LDQ2
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IPQ2 = IW + N
      IPU = IPQ2 + LDQ2*NQ
      IBUF = IPU + LDU*NQ
*     (IBUF est de taille 3*N au maximum)
*
      ICTOT = 1
      IPSM = ICTOT + NPCOL*4
      INDX = IPSM + NPCOL*4
      INDXC = INDX + N
      INDXP = INDXC + N
      INDCOL = INDXP + N
      COLTYP = INDCOL + N
      INDROW = COLTYP + N
      INDXR = INDROW + N
*
      CALL DESCINIT( DESCQ2, N, N, NB, NB, IQROW, IQCOL, ICTXT, LDQ2,
     $               INFO )
      CALL DESCINIT( DESCU, N, N, NB, NB, IQROW, IQCOL, ICTXT, LDU,
     $               INFO )
*
*     Form the z-vector which consists of the last row of Q_1 and the
*     first row of Q_2.
*
      IPWORK = IDLMDA
      CALL PDLAEDZ( N, N1, ID, Q, IQ, JQ, LDQ, DESCQ, WORK( IZ ),
     $              WORK( IPWORK ) )
*
*     Deflate eigenvalues.
*
      IPQ = IIQ + ( JJQ-1 )*LDQ
      CALL PDLAED2( ICTXT, K, N, N1, NB, D, IQROW, IQCOL, Q( IPQ ), LDQ,
     $              RHO, WORK( IZ ), WORK( IW ), WORK( IDLMDA ),
     $              WORK( IPQ2 ), LDQ2, WORK( IBUF ), IWORK( ICTOT ),
     $              IWORK( IPSM ), NPCOL, IWORK( INDX ), IWORK( INDXC ),
     $              IWORK( INDXP ), IWORK( INDCOL ), IWORK( COLTYP ),
     $              NN, NN1, NN2, IQ1, IQ2 )
*
*
*     Solve Secular Equation.
*
      IF( K.NE.0 ) THEN
         CALL PDLASET( 'A', N, N, ZERO, ONE, WORK( IPU ), 1, 1, DESCU )
         CALL PDLAED3( ICTXT, K, N, NB, D, IQROW, IQCOL, RHO,
     $                 WORK( IDLMDA ), WORK( IW ), WORK( IZ ),
     $                 WORK( IPU ), LDQ2, WORK( IBUF ), IWORK( INDX ),
     $                 IWORK( INDCOL ), IWORK( INDROW ), IWORK( INDXR ),
     $                 IWORK( INDXC ), IWORK( ICTOT ), NPCOL, INFO )
*
*     Compute the updated eigenvectors.
*
         IQQ = MIN( IQ1, IQ2 )
         IF( NN1.GT.0 ) THEN
            INQ = IQ - 1 + ID
            JNQ = JQ - 1 + ID + IQQ - 1
            CALL PDGEMM( 'N', 'N', N1, NN, NN1, ONE, WORK( IPQ2 ), 1,
     $                   IQ1, DESCQ2, WORK( IPU ), IQ1, IQQ, DESCU,
     $                   ZERO, Q, INQ, JNQ, DESCQ )
         END IF
         IF( NN2.GT.0 ) THEN
            INQ = IQ - 1 + ID + N1
            JNQ = JQ - 1 + ID + IQQ - 1
            CALL PDGEMM( 'N', 'N', N-N1, NN, NN2, ONE, WORK( IPQ2 ),
     $                   N1+1, IQ2, DESCQ2, WORK( IPU ), IQ2, IQQ,
     $                   DESCU, ZERO, Q, INQ, JNQ, DESCQ )
         END IF
*
         DO 10 J = K + 1, N
            JC = IWORK( INDX+J-1 )
            CALL INFOG1L( JQ-1+JC, NB, NPCOL, MYCOL, IQCOL, JJC, COL )
            CALL INFOG1L( JC, NB, NPCOL, MYCOL, IQCOL, JJ2C, COL )
            IF( MYCOL.EQ.COL ) THEN
               IQ2 = IPQ2 + ( JJ2C-1 )*LDQ2
               INQ = IPQ + ( JJC-1 )*LDQ
               CALL DCOPY( NP, WORK( IQ2 ), 1, Q( INQ ), 1 )
            END IF
   10    CONTINUE
      END IF
*
   20 CONTINUE
      RETURN
*
*     End of PDLAED1
*
      END
      SUBROUTINE DLAPST( ID, N, D, INDX, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            INDX( * )
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*  DLAPST is a modified version of the LAPACK routine DLASRT.
*
*  Define a permutation INDX that sorts the numbers in D
*  in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input)  DOUBLE PRECISION array, dimension (N)
*          The array to be sorted.
*
*  INDX    (ouput) INTEGER array, dimension (N).
*          The permutation which sorts the array D.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, ITMP, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAPST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, N
         INDX( I ) = I
   10 CONTINUE
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   20 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 40 I = START + 1, ENDD
               DO 30 J = I, START + 1, -1
                  IF( D( INDX( J ) ).GT.D( INDX( J-1 ) ) ) THEN
                     ITMP = INDX( J )
                     INDX( J ) = INDX( J-1 )
                     INDX( J-1 ) = ITMP
                  ELSE
                     GO TO 40
                  END IF
   30          CONTINUE
   40       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 60 I = START + 1, ENDD
               DO 50 J = I, START + 1, -1
                  IF( D( INDX( J ) ).LT.D( INDX( J-1 ) ) ) THEN
                     ITMP = INDX( J )
                     INDX( J ) = INDX( J-1 )
                     INDX( J-1 ) = ITMP
                  ELSE
                     GO TO 60
                  END IF
   50          CONTINUE
   60       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( INDX( START ) )
         D2 = D( INDX( ENDD ) )
         I = ( START+ENDD ) / 2
         D3 = D( INDX( I ) )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   70       CONTINUE
   80       CONTINUE
            J = J - 1
            IF( D( INDX( J ) ).LT.DMNMX )
     $         GO TO 80
   90       CONTINUE
            I = I + 1
            IF( D( INDX( I ) ).GT.DMNMX )
     $         GO TO 90
            IF( I.LT.J ) THEN
               ITMP = INDX( I )
               INDX( I ) = INDX( J )
               INDX( J ) = ITMP
               GO TO 70
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
  100       CONTINUE
  110       CONTINUE
            J = J - 1
            IF( D( INDX( J ) ).GT.DMNMX )
     $         GO TO 110
  120       CONTINUE
            I = I + 1
            IF( D( INDX( I ) ).LT.DMNMX )
     $         GO TO 120
            IF( I.LT.J ) THEN
               ITMP = INDX( I )
               INDX( I ) = INDX( J )
               INDX( J ) = ITMP
               GO TO 100
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 20
      RETURN
*
*     End of DLAPST
*
      END
      SUBROUTINE PDLAED2( ICTXT, K, N, N1, NB, D, DROW, DCOL, Q, LDQ,
     $                    RHO, Z, W, DLAMDA, Q2, LDQ2, QBUF, CTOT, PSM,
     $                    NPCOL, INDX, INDXC, INDXP, INDCOL, COLTYP, NN,
     $                    NN1, NN2, IB1, IB2 )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            DCOL, DROW, IB1, IB2, ICTXT, K, LDQ, LDQ2, N,
     $                   N1, NB, NN, NN1, NN2, NPCOL
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            COLTYP( * ), CTOT( 0: NPCOL-1, 4 ),
     $                   INDCOL( N ), INDX( * ), INDXC( * ), INDXP( * ),
     $                   PSM( 0: NPCOL-1, 4 )
      DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ),
     $                   Q2( LDQ2, * ), QBUF( * ), W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAED2 sorts the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny entry in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  ICTXT  (global input) INTEGER
*         The BLACS context handle, indicating the global context of
*         the operation on the matrix. The context itself is global.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  N1     (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) < N1 < N.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  DROW   (global input) INTEGER
*          The process row over which the first row of the matrix D is
*          distributed. 0 <= DROW < NPROW.
*
*  DCOL   (global input) INTEGER
*          The process column over which the first column of the
*          matrix D is distributed. 0 <= DCOL < NPCOL.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (N1,N1)
*         and (N1+1, N1+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,NQ).
*
*  RHO    (global input/output) DOUBLE PRECISION
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         PDLAED3.
*
*  Z      (global input) DOUBLE PRECISION array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  DLAMDA (global output) DOUBLE PRECISION array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED3 to form the secular equation.
*
*  W      (global output) DOUBLE PRECISION array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to SLAED3.
*
*  Q2     (output) DOUBLE PRECISION array, dimension (LDQ2, NQ)
*         A copy of the first K eigenvectors which will be used by
*
*  LDQ2    (input) INTEGER
*         The leading dimension of the array Q2.
*
*  QBUF   (workspace) DOUBLE PRECISION array, dimension 3*N
*
*  CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  PSM    (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  NPCOL   (global input) INTEGER
*          The total number of columns over which the distributed
*           submatrix is distributed.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDXC  (output) INTEGER array, dimension (N)
*         The permutation used to arrange the columns of the deflated
*         Q matrix into three groups:  the first group contains non-zero
*         elements only at and above N1, the second contains
*         non-zero elements only below N1, and the third is dense.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  INDCOL (workspace) INTEGER array, dimension (N)
*
*  COLTYP (workspace/output) INTEGER array, dimension (N)
*         During execution, a label which will indicate which of the
*         following types a column in the Q2 matrix is:
*         1 : non-zero in the upper half only;
*         2 : dense;
*         3 : non-zero in the lower half only;
*         4 : deflated.
*
*  NN     (global output) INTEGER, the order of matrix U, (PDLAED1).
*  NN1    (global output) INTEGER, the order of matrix Q1, (PDLAED1).
*  NN2    (global output) INTEGER, the order of matrix Q2, (PDLAED1).
*  IB1    (global output) INTEGER, pointeur on Q1, (PDLAED1).
*  IB2    (global output) INTEGER, pointeur on Q2, (PDLAED1).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0,
     $                   TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, CT, I, IAM, IE1, IE2, IMAX, INFO, J, JJQ2,
     $                   JJS, JMAX, JS, K2, MYCOL, MYROW, N1P1, N2, NJ,
     $                   NJCOL, NJJ, NP, NPROCS, NPROW, PJ, PJCOL, PJJ
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX, INDXG2L, INDXL2G, NUMROC
      DOUBLE PRECISION   DLAPY2, PDLAMCH
      EXTERNAL           IDAMAX, INDXG2L, INDXL2G, NUMROC, PDLAMCH,
     $                   DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_PINFO, DCOPY, DGERV2D,
     $                   DGESD2D, DLAPST, DROT, DSCAL, INFOG1L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, SQRT
*     ..
*     .. External Functions ..
*     ..
*     .. Local Arrays ..
      INTEGER            PTT( 4 )
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NUMROC( N, NB, MYROW, DROW, NPROW )
*
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1.  Since z is the concatenation of
*     two normalized vectors, norm2(z) = sqrt(2).
*
      T = ONE / SQRT( TWO )
      CALL DSCAL( N, T, Z, 1 )
*
*     RHO = ABS( norm(z)**2 * RHO )
*
      RHO = ABS( TWO*RHO )
*
*     Calculate the allowable deflation tolerance
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = PDLAMCH( ICTXT, 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         GO TO 220
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
*
      CALL DLAPST( 'I', N, D, INDX, INFO )
*
      DO 10 I = 1, N1
         COLTYP( I ) = 1
   10 CONTINUE
      DO 20 I = N1P1, N
         COLTYP( I ) = 3
   20 CONTINUE
      COL = DCOL
      DO 40 I = 1, N, NB
         DO 30 J = 0, NB - 1
            IF( I+J.LE.N )
     $         INDCOL( I+J ) = COL
   30    CONTINUE
         COL = MOD( COL+1, NPCOL )
   40 CONTINUE
*
      K = 0
      K2 = N + 1
      DO 50 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N )
     $         GO TO 80
         ELSE
            PJ = NJ
            GO TO 60
         END IF
   50 CONTINUE
   60 CONTINUE
      J = J + 1
      NJ = INDX( J )
      IF( J.GT.N )
     $   GO TO 80
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( PJ )
         C = Z( NJ )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = DLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) )
     $         COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4
            CALL INFOG1L( NJ, NB, NPCOL, MYCOL, DCOL, NJJ, NJCOL )
            CALL INFOG1L( PJ, NB, NPCOL, MYCOL, DCOL, PJJ, PJCOL )
            IF( INDCOL( PJ ).EQ.INDCOL( NJ ) .AND. MYCOL.EQ.NJCOL ) THEN
               CALL DROT( NP, Q( 1, PJJ ), 1, Q( 1, NJJ ), 1, C, S )
            ELSE IF( MYCOL.EQ.PJCOL ) THEN
               CALL DGESD2D( ICTXT, NP, 1, Q( 1, PJJ ), NP, MYROW,
     $                       NJCOL )
               CALL DGERV2D( ICTXT, NP, 1, QBUF, NP, MYROW, NJCOL )
               CALL DROT( NP, Q( 1, PJJ ), 1, QBUF, 1, C, S )
            ELSE IF( MYCOL.EQ.NJCOL ) THEN
               CALL DGESD2D( ICTXT, NP, 1, Q( 1, NJJ ), NP, MYROW,
     $                       PJCOL )
               CALL DGERV2D( ICTXT, NP, 1, QBUF, NP, MYROW, PJCOL )
               CALL DROT( NP, QBUF, 1, Q( 1, NJJ ), 1, C, S )
            END IF
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T
            K2 = K2 - 1
            I = 1
   70       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 70
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 60
   80 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
*
*     Count up the total number of the various types of columns, then
*     form a permutation which positions the four column types into
*     four uniform groups (although one or more of these groups may be
*     empty).
*
      DO 100 J = 1, 4
         DO 90 I = 0, NPCOL - 1
            CTOT( I, J ) = 0
   90    CONTINUE
         PTT( J ) = 0
  100 CONTINUE
      DO 110 J = 1, N
         CT = COLTYP( J )
         COL = INDCOL( J )
         CTOT( COL, CT ) = CTOT( COL, CT ) + 1
  110 CONTINUE
*
*     PSM(*) = Position in SubMatrix (of types 1 through 4)
*
      DO 120 COL = 0, NPCOL - 1
         PSM( COL, 1 ) = 1
         PSM( COL, 2 ) = 1 + CTOT( COL, 1 )
         PSM( COL, 3 ) = PSM( COL, 2 ) + CTOT( COL, 2 )
         PSM( COL, 4 ) = PSM( COL, 3 ) + CTOT( COL, 3 )
  120 CONTINUE
      PTT( 1 ) = 1
      DO 140 I = 2, 4
         CT = 0
         DO 130 J = 0, NPCOL - 1
            CT = CT + CTOT( J, I-1 )
  130    CONTINUE
         PTT( I ) = PTT( I-1 ) + CT
  140 CONTINUE
*
*     Fill out the INDXC array so that the permutation which it induces
*     will place all type-1 columns first, all type-2 columns next,
*     then all type-3's, and finally all type-4's.
*
      DO 150 J = 1, N
         JS = INDXP( J )
         COL = INDCOL( JS )
         CT = COLTYP( JS )
         I = INDXL2G( PSM( COL, CT ), NB, COL, DCOL, NPCOL )
         INDX( J ) = I
         INDXC( PTT( CT ) ) = I
         PSM( COL, CT ) = PSM( COL, CT ) + 1
         PTT( CT ) = PTT( CT ) + 1
  150 CONTINUE
*
*
      DO 160 J = 1, N
         JS = INDXP( J )
         JJS = INDXG2L( JS, NB, J, J, NPCOL )
         COL = INDCOL( JS )
         IF( COL.EQ.MYCOL ) THEN
            I = INDX( J )
            JJQ2 = INDXG2L( I, NB, J, J, NPCOL )
            CALL DCOPY( NP, Q( 1, JJS ), 1, Q2( 1, JJQ2 ), 1 )
         END IF
  160 CONTINUE
*
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      CALL DCOPY( N, D, 1, Z, 1 )
      DO 170 J = K + 1, N
         JS = INDXP( J )
         I = INDX( J )
         D( I ) = Z( JS )
  170 CONTINUE
*
      PTT( 1 ) = 1
      DO 190 I = 2, 4
         CT = 0
         DO 180 J = 0, NPCOL - 1
            CT = CT + CTOT( J, I-1 )
  180    CONTINUE
         PTT( I ) = PTT( I-1 ) + CT
  190 CONTINUE
*
*
      IB1 = INDXC( 1 )
      IE1 = IB1
      IB2 = INDXC( PTT( 2 ) )
      IE2 = IB2
      DO 200 I = 2, PTT( 3 ) - 1
         IB1 = MIN( IB1, INDXC( I ) )
         IE1 = MAX( IE1, INDXC( I ) )
  200 CONTINUE
      DO 210 I = PTT( 2 ), PTT( 4 ) - 1
         IB2 = MIN( IB2, INDXC( I ) )
         IE2 = MAX( IE2, INDXC( I ) )
  210 CONTINUE
      NN1 = IE1 - IB1 + 1
      NN2 = IE2 - IB2 + 1
      NN = MAX( IE1, IE2 ) - MIN( IB1, IB2 ) + 1
  220 CONTINUE
      RETURN
*
*     End of PDLAED2
*
      END
      SUBROUTINE PDLAED3( ICTXT, K, N, NB, D, DROW, DCOL, RHO, DLAMDA,
     $                    W, Z, U, LDU, BUF, INDX, INDCOL, INDROW,
     $                    INDXR, INDXC, CTOT, NPCOL, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            DCOL, DROW, ICTXT, INFO, K, LDU, N, NB, NPCOL
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( 0: NPCOL-1, 4 ), INDCOL( * ),
     $                   INDROW( * ), INDX( * ), INDXC( * ), INDXR( * )
      DOUBLE PRECISION   BUF( * ), D( * ), DLAMDA( * ), U( LDU, * ),
     $                   W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAED3 finds the roots of the secular equation, as defined by the
*  values in D, W, and RHO, between 1 and K.  It makes the
*  appropriate calls to SLAED4
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  ICTXT  (global input) INTEGER
*         The BLACS context handle, indicating the global context of
*         the operation on the matrix. The context itself is global.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  D      (input/output) DOUBLE PRECISION array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  DROW   (global input) INTEGER
*          The process row over which the first row of the matrix D is
*          distributed. 0 <= DROW < NPROW.
*
*  DCOL   (global input) INTEGER
*          The process column over which the first column of the
*          matrix D is distributed. 0 <= DCOL < NPCOL.
*
*  Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (N1,N1)
*         and (N1+1, N1+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,NQ).
*
*  RHO    (global input/output) DOUBLE PRECISION
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         PDLAED3.
*
*  DLAMDA (global output) DOUBLE PRECISION array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED3 to form the secular equation.
*
*  W      (global output) DOUBLE PRECISION array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to SLAED3.
*
*  Z      (global input) DOUBLE PRECISION array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  U     (global output) DOUBLE PRECISION array
*         global dimension (N, N), local dimension (LDU, NQ).
*         Q  contains the orthonormal eigenvectors of the symmetric
*         tridiagonal matrix.
*
*  LDU    (input) INTEGER
*         The leading dimension of the array U.
*
*  QBUF   (workspace) DOUBLE PRECISION array, dimension 3*N
*
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDCOL (workspace) INTEGER array, dimension (N)
*
*
*  INDROW (workspace) INTEGER array, dimension (N)
*
*
*  INDXR (workspace) INTEGER array, dimension (N)
*
*
*  INDXC (workspace) INTEGER array, dimension (N)
*
*  CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  NPCOL   (global input) INTEGER
*          The total number of columns over which the distributed
*           submatrix is distributed.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute the ith eigenvalue.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, GI, I, IINFO, IIU, IPD, IU, J, JJU, JU,
     $                   KK, KL, KLC, KLR, MYCOL, MYKL, MYKLR, MYROW,
     $                   NPROW, PDC, PDR, ROW
      DOUBLE PRECISION   AUX, TEMP
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           INDXG2L, DLAMC3, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     $                   DGERV2D, DGESD2D, DLAED4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      ROW = DROW
      COL = DCOL
      DO 20 I = 1, N, NB
         DO 10 J = 0, NB - 1
            IF( I+J.LE.N ) THEN
               INDROW( I+J ) = ROW
               INDCOL( I+J ) = COL
            END IF
   10    CONTINUE
         ROW = MOD( ROW+1, NPROW )
         COL = MOD( COL+1, NPCOL )
   20 CONTINUE
*
      MYKL = CTOT( MYCOL, 1 ) + CTOT( MYCOL, 2 ) + CTOT( MYCOL, 3 )
      KLR = MYKL / NPROW
      IF( MYROW.EQ.DROW ) THEN
         MYKLR = KLR + MOD( MYKL, NPROW )
      ELSE
         MYKLR = KLR
      END IF
      PDC = 1
      COL = DCOL
   30 CONTINUE
      IF( MYCOL.NE.COL ) THEN
         PDC = PDC + CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
         COL = MOD( COL+1, NPCOL )
         GO TO 30
      END IF
      PDR = PDC
      KL = KLR + MOD( MYKL, NPROW )
      ROW = DROW
   40 CONTINUE
      IF( MYROW.NE.ROW ) THEN
         PDR = PDR + KL
         KL = KLR
         ROW = MOD( ROW+1, NPROW )
         GO TO 40
      END IF
*
      DO 50 I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
         Z( I ) = ONE
   50 CONTINUE
      IF( MYKLR.GT.0 ) THEN
         KK = PDR
         DO 80 I = 1, MYKLR
            CALL DLAED4( K, KK, DLAMDA, W, BUF, RHO, BUF( K+I ), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = KK
            END IF
*
*     ..Compute part of z
*
            DO 60 J = 1, KK - 1
               Z( J ) = Z( J )*( BUF( J ) /
     $                  ( DLAMDA( J )-DLAMDA( KK ) ) )
   60       CONTINUE
            Z( KK ) = Z( KK )*BUF( KK )
            DO 70 J = KK + 1, K
               Z( J ) = Z( J )*( BUF( J ) /
     $                  ( DLAMDA( J )-DLAMDA( KK ) ) )
   70       CONTINUE
            KK = KK + 1
   80    CONTINUE
*
         IF( MYROW.NE.DROW ) THEN
            CALL DCOPY( K, Z, 1, BUF, 1 )
            CALL DGESD2D( ICTXT, K+MYKLR, 1, BUF, K+MYKLR, DROW, MYCOL )
         ELSE
            IPD = 2*K + 1
            CALL DCOPY( MYKLR, BUF( K+1 ), 1, BUF( IPD ), 1 )
            IF( KLR.GT.0 ) THEN
               IPD = MYKLR + IPD
               ROW = MOD( DROW+1, NPROW )
               DO 100 I = 1, NPROW - 1
                  CALL DGERV2D( ICTXT, K+KLR, 1, BUF, K+KLR, ROW,
     $                          MYCOL )
                  CALL DCOPY( KLR, BUF( K+1 ), 1, BUF( IPD ), 1 )
                  DO 90 J = 1, K
                     Z( J ) = Z( J )*BUF( J )
   90             CONTINUE
                  IPD = IPD + KLR
                  ROW = MOD( ROW+1, NPROW )
  100          CONTINUE
            END IF
         END IF
      END IF
*
      IF( MYROW.EQ.DROW ) THEN
         IF( MYCOL.NE.DCOL .AND. MYKL.NE.0 ) THEN
            CALL DCOPY( K, Z, 1, BUF, 1 )
            CALL DCOPY( MYKL, BUF( 2*K+1 ), 1, BUF( K+1 ), 1 )
            CALL DGESD2D( ICTXT, K+MYKL, 1, BUF, K+MYKL, MYROW, DCOL )
         ELSE IF( MYCOL.EQ.DCOL ) THEN
            IPD = 2*K + 1
            COL = DCOL
            KL = MYKL
            DO 120 I = 1, NPCOL - 1
               IPD = IPD + KL
               COL = MOD( COL+1, NPCOL )
               KL = CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
               IF( KL.NE.0 ) THEN
                  CALL DGERV2D( ICTXT, K+KL, 1, BUF, K+KL, MYROW, COL )
                  CALL DCOPY( KL, BUF( K+1 ), 1, BUF( IPD ), 1 )
                  DO 110 J = 1, K
                     Z( J ) = Z( J )*BUF( J )
  110             CONTINUE
               END IF
  120       CONTINUE
            DO 130 I = 1, K
               Z( I ) = SIGN( SQRT( -Z( I ) ), W( I ) )
  130       CONTINUE
*
         END IF
      END IF
*
*     Diffusion
*
      IF( MYROW.EQ.DROW .AND. MYCOL.EQ.DCOL ) THEN
         CALL DCOPY( K, Z, 1, BUF, 1 )
         CALL DCOPY( K, BUF( 2*K+1 ), 1, BUF( K+1 ), 1 )
         CALL DGEBS2D( ICTXT, 'All', ' ', 2*K, 1, BUF, 2*K )
      ELSE
         CALL DGEBR2D( ICTXT, 'All', ' ', 2*K, 1, BUF, 2*K, DROW, DCOL )
         CALL DCOPY( K, BUF, 1, Z, 1 )
      END IF
*
*     Copy of D at the good place
*
      KLC = 0
      KLR = 0
      DO 140 I = 1, K
         GI = INDX( I )
         D( GI ) = BUF( K+I )
         COL = INDCOL( GI )
         ROW = INDROW( GI )
         IF( COL.EQ.MYCOL ) THEN
            KLC = KLC + 1
            INDXC( KLC ) = I
         END IF
         IF( ROW.EQ.MYROW ) THEN
            KLR = KLR + 1
            INDXR( KLR ) = I
         END IF
  140 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      IF( MYKL.NE.0 ) THEN
         DO 180 J = 1, MYKL
            KK = INDXC( J )
            JU = INDX( KK )
            JJU = INDXG2L( JU, NB, J, J, NPCOL )
            CALL DLAED4( K, KK, DLAMDA, W, BUF, RHO, AUX, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = KK
            END IF
            IF( K.EQ.1 .OR. K.EQ.2 ) THEN
               DO 150 I = 1, KLR
                  KK = INDXR( I )
                  IU = INDX( KK )
                  IIU = INDXG2L( IU, NB, J, J, NPROW )
                  U( IIU, JJU ) = BUF( KK )
  150          CONTINUE
               GO TO 180
            END IF
*
            DO 160 I = 1, K
               BUF( I ) = Z( I ) / BUF( I )
  160       CONTINUE
            TEMP = DNRM2( K, BUF, 1 )
            DO 170 I = 1, KLR
               KK = INDXR( I )
               IU = INDX( KK )
               IIU = INDXG2L( IU, NB, J, J, NPROW )
               U( IIU, JJU ) = BUF( KK ) / TEMP
  170       CONTINUE
*
  180    CONTINUE
      END IF
*
  190 CONTINUE
*
      RETURN
*
*     End of PDLAED3
*
      END
      SUBROUTINE PDLAEDZ( N, N1, ID, Q, IQ, JQ, LDQ, DESCQ, Z, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ID, IQ, JQ, LDQ, N, N1
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * )
      DOUBLE PRECISION   Q( LDQ, * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAEDZ Form the z-vector which consists of the last row of Q_1
*  and the first row of Q_2.
*  =====================================================================
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            COL, I, IBUF, ICTXT, IIQ, IIZ1, IIZ2, IQCOL,
     $                   IQROW, IZ, IZ1, IZ1COL, IZ1ROW, IZ2, IZ2COL,
     $                   IZ2ROW, J, JJQ, JJZ1, JJZ2, MYCOL, MYROW, N2,
     $                   NB, NBLOC, NPCOL, NPROW, NQ1, NQ2, ZSIZ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     $                   DGERV2D, DGESD2D, INFOG2L
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      ICTXT = DESCQ( CTXT_ )
      NB = DESCQ( NB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( ID, ID, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     $              IQROW, IQCOL )
      N2 = N - N1
*
*     Form z1 which consist of the last row of Q1
*
      CALL INFOG2L( IQ-1+( ID+N1-1 ), JQ-1+ID, DESCQ, NPROW, NPCOL,
     $              MYROW, MYCOL, IIZ1, JJZ1, IZ1ROW, IZ1COL )
      NQ1 = NUMROC( N1, NB, MYCOL, IZ1COL, NPCOL )
      IF( ( MYROW.EQ.IZ1ROW ) .AND. ( NQ1.NE.0 ) ) THEN
         CALL DCOPY( NQ1, Q( IIZ1, JJZ1 ), LDQ, WORK, 1 )
         IF( MYROW.NE.IQROW .OR. MYCOL.NE.IQCOL )
     $      CALL DGESD2D( ICTXT, NQ1, 1, WORK, NQ1, IQROW, IQCOL )
      END IF
*
*     Proc (IQROW, IQCOL) receive the parts of z1
*
      IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL ) THEN
         COL = IZ1COL
         DO 20 I = 0, NPCOL - 1
            NQ1 = NUMROC( N1, NB, COL, IZ1COL, NPCOL )
            IF( NQ1.GT.0 ) THEN
               IF( IZ1ROW.NE.IQROW .OR. COL.NE.IQCOL ) THEN
                  IBUF = N1 + 1
                  CALL DGERV2D( ICTXT, NQ1, 1, WORK( IBUF ), NQ1,
     $                          IZ1ROW, COL )
               ELSE
                  IBUF = 1
               END IF
               IZ1 = 0
               IZ = I*NB + 1
               NBLOC = ( NQ1-1 ) / NB + 1
               DO 10 J = 1, NBLOC
                  ZSIZ = MIN( NB, NQ1-IZ1 )
                  CALL DCOPY( ZSIZ, WORK( IBUF+IZ1 ), 1, Z( IZ ), 1 )
                  IZ1 = IZ1 + NB
                  IZ = IZ + NB*NPCOL
   10          CONTINUE
            END IF
            COL = MOD( COL+1, NPCOL )
   20    CONTINUE
      END IF
*
*     Form z2 which consist of the first row of Q2
*
      CALL INFOG2L( IQ-1+( ID+N1 ), JQ-1+( ID+N1 ), DESCQ, NPROW, NPCOL,
     $              MYROW, MYCOL, IIZ2, JJZ2, IZ2ROW, IZ2COL )
      NQ2 = NUMROC( N2, NB, MYCOL, IZ2COL, NPCOL )
      IF( ( MYROW.EQ.IZ2ROW ) .AND. ( NQ2.NE.0 ) ) THEN
         CALL DCOPY( NQ2, Q( IIZ2, JJZ2 ), LDQ, WORK, 1 )
         IF( MYROW.NE.IQROW .OR. MYCOL.NE.IQCOL )
     $      CALL DGESD2D( ICTXT, NQ2, 1, WORK, NQ2, IQROW, IQCOL )
      END IF
*
*     Proc (IQROW, IQCOL) receive the parts of z2
*
      IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL ) THEN
         COL = IZ2COL
         DO 40 I = 0, NPCOL - 1
            NQ2 = NUMROC( N2, NB, COL, IZ2COL, NPCOL )
            IF( NQ2.GT.0 ) THEN
               IF( IQROW.NE.IZ2ROW .OR. IQCOL.NE.COL ) THEN
                  IBUF = 1 + N2
                  CALL DGERV2D( ICTXT, NQ2, 1, WORK( IBUF ), NQ2,
     $                          IZ2ROW, COL )
               ELSE
                  IBUF = 1
               END IF
               IZ2 = 0
               IZ = NB*I + N1 + 1
               NBLOC = ( NQ2-1 ) / NB + 1
               DO 30 J = 1, NBLOC
                  ZSIZ = MIN( NB, NQ2-IZ2 )
                  CALL DCOPY( ZSIZ, WORK( IBUF+IZ2 ), 1, Z( IZ ), 1 )
                  IZ2 = IZ2 + NB
                  IZ = IZ + NB*NPCOL
   30          CONTINUE
            END IF
            COL = MOD( COL+1, NPCOL )
   40    CONTINUE
      END IF
*
*     proc(IQROW,IQCOL) broadcast Z=(Z1,Z2)
*
      IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL ) THEN
         CALL DGEBS2D( ICTXT, 'All', ' ', N, 1, Z, N )
      ELSE
         CALL DGEBR2D( ICTXT, 'All', ' ', N, 1, Z, N, IQROW, IQCOL )
      END IF
*
      RETURN
*
*     End of PDLAEDZ
*
      END
      SUBROUTINE PZHEEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
     $                    DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK,
     $                    LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 25, 2002
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( * ), WORK( * ), Z( * )
*     
*
*  Purpose
*  =======
*
*  PZHEEVD computes all the eigenvalues and eigenvectors of a Hermitian
*  matrix A by using a divide and conquer algorithm.
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;    (NOT IMPLEMENTED YET)
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (local input/workspace) block cyclic COMPLEX*16 array,
*          global dimension (N, N), local dimension ( LLD_A,
*          LOCc(JA+N-1) )
*
*          On entry, the symmetric matrix A.  If UPLO = 'U', only the
*          upper triangular part of A is used to define the elements of
*          the symmetric matrix.  If UPLO = 'L', only the lower
*          triangular part of A is used to define the elements of the
*          symmetric matrix.
*
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          If DESCA( CTXT_ ) is incorrect, PZHEEV cannot guarantee
*          correct error reporting.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          If INFO=0, the eigenvalues in ascending order.
*
*  Z       (local output) COMPLEX*16 array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          Z contains the orthonormal eigenvectors of the matrix A.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) COMPLEX*16 array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace needed for the
*          computation.
*
*  LWORK   (local input) INTEGER
*          If eigenvectors are requested:
*            LWORK = N + ( NP0 + MQ0 + NB ) * NB,
*          with  NP0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPROW )
*                MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine calculates the size for all
*          work arrays. Each of these values is returned in the first
*          entry of the corresponding work array, and no error message
*          is issued by PXERBLA.
*
*  RWORK   (local workspace/output) DOUBLE PRECISION array,
*          dimension (LRWORK)
*          On output RWORK(1) returns the real workspace needed to
*          guarantee completion.  If the input parameters are incorrect,
*          RWORK(1) may also be incorrect.
*
*  LRWORK  (local input) INTEGER
*          Size of RWORK array.
*          LRWORK >= 1 + 9*N + 3*NP*NQ,
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
*          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On output IWORK(1) returns the integer workspace needed.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = 7*N + 8*NPCOL + 2
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
*                converge in PDLAED3.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices sub( A ), sub( Z ) must verify
*  some alignment properties, namely the following expression
*  should be true:
*  ( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
*    IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
*    with IROFFA = MOD( IA-1, MB_A )
*     and ICOFFA = MOD( JA-1, NB_A ).
*
*  Further Details
*  ======= =======
*
*  Contributed by Francoise Tisseur, University of Manchester.
*
*  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
*              Conquer Algorithm for the Symmetric Eigenvalue Problem
*              on Distributed Memory Architectures",
*              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
*              (see also LAPACK Working Note 132)
*                http://www.netlib.org/lapack/lawns/lawn132.ps
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION               ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY
      INTEGER            CSRC_A, I, IACOL, IAROW, ICOFFA, IINFO, IIZ,
     $                   INDD, INDE, INDE2, INDRWORK, INDTAU, INDWORK,
     $                   INDZ, IPR, IPZ, IROFFA, IROFFZ, ISCALE, IZCOL,
     $                   IZROW, J, JJZ, LDR, LDZ, LIWMIN, LLRWORK,
     $                   LLWORK, LRWMIN, LWMIN, MB_A, MYCOL, MYROW, NB,
     $                   NB_A, NN, NP0, NPCOL, NPROW, NQ, NQ0, OFFSET,
     $                   RSRC_A
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCRZ( 9 ), IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      DOUBLE PRECISION   PZLANHE, PDLAMCH
      EXTERNAL           LSAME, INDXG2L, INDXG2P, NUMROC, PZLANHE,
     $                   PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCINIT, INFOG2L,
     $                   PZELGET, PZHETRD, PCHK2MAT, PZLASCL, PZLASET,
     $                   PZUNMTR, PDLARED1D, PDLASET, PDSTEDC, PXERBLA,
     $                   DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, ICHAR, MAX, MIN, MOD, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      INFO = 0
*
*     Quick return
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Test the input arguments.
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE 
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         CALL CHK1MAT( N, 2, N, 2, IZ, JZ, DESCZ, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            LOWER = LSAME( UPLO, 'L' )
            NB_A = DESCA( NB_ )
            MB_A = DESCA( MB_ )
            NB = NB_A
            RSRC_A = DESCA( RSRC_ )
            CSRC_A = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A )
            ICOFFA = MOD( JA-1, NB_A )
            IAROW = INDXG2P( IA, NB_A, MYROW, RSRC_A, NPROW )
            IACOL = INDXG2P( JA, MB_A, MYCOL, CSRC_A, NPCOL )
            NP0 = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ0 = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
            IROFFZ = MOD( IZ-1, MB_A )
            CALL INFOG2L( IZ, JZ, DESCZ, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIZ, JJZ, IZROW, IZCOL )
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
*
*           Compute the total amount of space needed
*
            NN = MAX( N, NB, 2 )
            NQ = NUMROC( NN, NB, 0, 0, NPCOL )
            LWMIN = N + ( NP0+NQ+NB )*NB
            LRWMIN = 1 + 9*N + 3*NP0*NQ0
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DCMPLX( LWMIN )
            RWORK( 1 ) = DBLE( LRWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -14
            ELSE IF( LRWORK.LT.LRWMIN .AND. LRWORK.NE.-1 ) THEN
               INFO = -16
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -4
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFA.NE.IROFFZ ) THEN
               INFO = -10
            ELSE IF( IAROW.NE.IZROW ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
               INFO = -( 1200+N_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            END IF
         END IF
         IF( LOWER ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'U' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IZ,
     $                  JZ, DESCZ, 11, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PZHEEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Set up pointers into the WORK array
*
      INDTAU = 1
      INDWORK = INDTAU + N
      LLWORK = LWORK - INDWORK + 1
*
*     Set up pointers into the RWORK array
*
      INDE = 1
      INDD = INDE + N
      INDE2 = INDD + N
      INDRWORK = INDE2 + N
      LLRWORK = LRWORK - INDRWORK + 1
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
*
      ANRM = PZLANHE( 'M', UPLO, N, A, IA, JA, DESCA,
     $       RWORK( INDRWORK ) )
*
*
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PZLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
*
*     Reduce Hermitian matrix to tridiagonal form.
*
      CALL PZHETRD( UPLO, N, A, IA, JA, DESCA, RWORK( INDD ),
     $              RWORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     Copy the values of D, E to all processes
*
*     Here PxLARED1D is used to redistribute the tridiagonal matrix.
*     PxLARED1D, however, doesn't yet workMx Mawith arbritary matrix
*     distributions so we have PxELGET as a backup.
*
      OFFSET = 0
      IF( IA.EQ.1 .AND. JA.EQ.1 .AND. RSRC_A.EQ.0 .AND. CSRC_A.EQ.0 )
     $     THEN
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDD ), W,
     $                   RWORK( INDRWORK ), LLRWORK )
*
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDE2 ),
     $                   RWORK( INDE ), RWORK( INDRWORK ), LLRWORK )
         IF( .NOT.LOWER )
     $      OFFSET = 1
      ELSE
         DO 10 I = 1, N
            CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA-1,
     $                    DESCA )
            W( I ) = DBLE( WORK( INDWORK ) )
   10    CONTINUE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA,
     $                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   20       CONTINUE
         ELSE
            DO 30 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA, I+JA-1,
     $                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   30       CONTINUE
         END IF
      END IF
*
*     Call PDSTEDC to compute eigenvalues and eigenvectors.
*
      INDZ = INDE + N
      INDRWORK = INDZ + NP0*NQ0
      LLRWORK = LRWORK - INDRWORK + 1
      LDR = MAX( 1, NP0 )
      CALL DESCINIT( DESCRZ, DESCZ( M_ ), DESCZ( N_ ), DESCZ( MB_ ),
     $               DESCZ( NB_ ), DESCZ( RSRC_ ), DESCZ( CSRC_ ),
     $               DESCZ( CTXT_ ), LDR, INFO )
      CALL PZLASET( 'Full', N, N, CZERO, CONE, Z, IZ, JZ, DESCZ )
      CALL PDLASET( 'Full', N, N, ZERO, ONE, RWORK( INDZ ), 1, 1,
     $              DESCRZ )
      CALL PDSTEDC( 'I', N, W, RWORK( INDE+OFFSET ), RWORK( INDZ ), IZ,
     $              JZ, DESCRZ, RWORK( INDRWORK ), LLRWORK, IWORK,
     $              LIWORK, IINFO )
*
      LDZ = DESCZ( LLD_ )
      LDR = DESCRZ( LLD_ )
      IIZ = INDXG2L( IZ, NB, MYROW, MYROW, NPROW )
      JJZ = INDXG2L( JZ, NB, MYCOL, MYCOL, NPCOL )
      IPZ = IIZ + ( JJZ-1 )*LDZ
      IPR = INDZ - 1 + IIZ + ( JJZ-1 )*LDR
      DO 50 J = 0, NQ0 - 1
         DO 40 I = 0, NP0 - 1
            Z( IPZ+I+J*LDZ ) = RWORK( IPR+I+J*LDR )
   40    CONTINUE
   50 CONTINUE
*
*     Z = Q * Z
*
      CALL PZUNMTR( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 ) = DCMPLX( LWMIN )
      RWORK( 1 ) = DBLE( LRWMIN )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of PZHEEVD
*
      END
      SUBROUTINE PZHENGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                     DESCB, SCALE, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, IBTYPE, INFO, JA, JB, LWORK, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*16         A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*
*  =======
*
*  PZHENGST reduces a complex Hermitian-definite generalized
*  eigenproblem to standard form.
*
*  PZHENGST performs the same function as PZHEGST, but is based on
*  rank 2K updates, which are faster and more scalable than
*  triangular solves (the basis of PZHENGST).
*
*  PZHENGST calls PZHEGST when UPLO='U', hence PZHENGST provides
*  improved performance only when UPLO='L', IBTYPE=1.
*
*  PZHENGST also calls PZHEGST when insufficient workspace is
*  provided,  hence PZHENGST provides improved
*  performance only when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*  In the following sub( A ) denotes A( IA:IA+N-1, JA:JA+N-1 ) and
*  sub( B ) denotes B( IB:IB+N-1, JB:JB+N-1 ).
*
*  If IBTYPE = 1, the problem is sub( A )*x = lambda*sub( B )*x,
*  and sub( A ) is overwritten by inv(U**H)*sub( A )*inv(U) or
*  inv(L)*sub( A )*inv(L**H)
*
*  If IBTYPE = 2 or 3, the problem is sub( A )*sub( B )*x = lambda*x or
*  sub( B )*sub( A )*x = lambda*x, and sub( A ) is overwritten by
*  U*sub( A )*U**H or L**H*sub( A )*L.
*
*  sub( B ) must have been previously factorized as U**H*U or L*L**H by
*  PZPOTRF.
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
*  IBTYPE   (global input) INTEGER
*          = 1: compute inv(U**H)*sub( A )*inv(U) or
*               inv(L)*sub( A )*inv(L**H);
*          = 2 or 3: compute U*sub( A )*U**H or L**H*sub( A )*L.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored and sub( B ) is
*                  factored as U**H*U;
*          = 'L':  Lower triangle of sub( A ) is stored and sub( B ) is
*                  factored as L*L**H.
*
*  N       (global input) INTEGER
*          The order of the matrices sub( A ) and sub( B ).  N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
*          the leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix, and its strictly
*          lower triangular part is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains
*          the lower triangular part of the matrix, and its strictly
*          upper triangular part is not referenced.
*
*          On exit, if INFO = 0, the transformed matrix, stored in the
*          same format as sub( A ).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)). On entry,
*          this array contains the local pieces of the triangular factor
*          from the Cholesky factorization of sub( B ), as returned by
*          PZPOTRF.
*
*  IB      (global input) INTEGER
*          B's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JB      (global input) INTEGER
*          B's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  SCALE   (global output) DOUBLE PRECISION
*          Amount by which the eigenvalues should be scaled to
*          compensate for the scaling performed in this routine.
*          At present, SCALE is always returned as 1.0, it is
*          returned here to allow for future enhancement.
*
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( NB * ( NP0 +1 ), 3 * NB )
*
*          When IBTYPE = 1 and UPLO = 'L', PZHENGST provides improved
*          performance when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*          where NB = MB_A = NB_A,
*          NP0 = NUMROC( N, NB, 0, 0, NPROW ),
*          NQ0 = NUMROC( N, NB, 0, 0, NPROW ),
*
*          NUMROC ia a ScaLAPACK tool functions
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the
*          optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*
*
*     .. Parameters ..
      COMPLEX*16         ONEHALF, ONE, MONE
      DOUBLE PRECISION   RONE
      PARAMETER          ( ONEHALF = ( 0.5D0, 0.0D0 ),
     $                   ONE = ( 1.0D0, 0.0D0 ),
     $                   MONE = ( -1.0D0, 0.0D0 ), RONE = 1.0D0 )
      INTEGER            DLEN_, CTXT_, MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( DLEN_ = 9, CTXT_ = 2, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, INDAA, INDG, INDR, INDRT, IROFFA,
     $                   IROFFB, J, K, KB, LWMIN, LWOPT, MYCOL, MYROW,
     $                   NB, NP0, NPCOL, NPK, NPROW, NQ0, POSTK
*     ..
*     .. Local Arrays ..
      INTEGER            DESCAA( DLEN_ ), DESCG( DLEN_ ),
     $                   DESCR( DLEN_ ), DESCRT( DLEN_ ), IDUM1( 2 ),
     $                   IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           LSAME, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   PXERBLA, PZGEMM, PZHEGST, PZHEMM, PZHER2K,
     $                   PZLACPY, PZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      SCALE = 1.0D0
*
      NB = DESCA( MB_ )
*
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $              NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $              NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $              NPCOL )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            NP0 = NUMROC( N, NB, 0, 0, NPROW )
            NQ0 = NUMROC( N, NB, 0, 0, NPCOL )
            LWMIN = MAX( NB*( NP0+1 ), 3*NB )
            IF( IBTYPE.EQ.1 .AND. .NOT.UPPER ) THEN
               LWOPT = 2*NP0*NB + NQ0*NB + NB*NB
            ELSE
               LWOPT = LWMIN
            END IF
            WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IBTYPE.LT.1 .OR. IBTYPE.GT.3 ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( N.LT.0 ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -9
            ELSE IF( ICOFFB.NE.0 .OR. IBCOL.NE.IACOL ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( MB_ ) ) THEN
               INFO = -( 1100+MB_ )
            ELSE IF( DESCB( NB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 1100+NB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1100+CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
            END IF
         END IF
         IDUM1( 1 ) = IBTYPE
         IDUM2( 1 ) = 1
         IF( UPPER ) THEN
            IDUM1( 2 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'L' )
         END IF
         IDUM2( 2 ) = 2
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IB,
     $                  JB, DESCB, 11, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZHENGST', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
      IF( IBTYPE.NE.1 .OR. UPPER .OR. LWORK.LT.LWOPT ) THEN
         CALL PZHEGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                 DESCB, SCALE, INFO )
         RETURN
      END IF
*
      CALL DESCSET( DESCG, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCR, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCRT, NB, N, NB, NB, IAROW, IACOL, ICTXT, NB )
      CALL DESCSET( DESCAA, NB, NB, NB, NB, IAROW, IACOL, ICTXT, NB )
*
      INDG = 1
      INDR = INDG + DESCG( LLD_ )*NB
      INDAA = INDR + DESCR( LLD_ )*NB
      INDRT = INDAA + DESCAA( LLD_ )*NB
*
      DO 30 K = 1, N, NB
*
         KB = MIN( N-K+1, NB )
         POSTK = K + KB
         NPK = N - POSTK + 1
*
*
         CALL PZLACPY( 'A', N-POSTK+1, KB, B, POSTK+IB-1, K+JB-1, DESCB,
     $                 WORK( INDG ), POSTK, 1, DESCG )
         CALL PZLACPY( 'A', N-POSTK+1, KB, A, POSTK+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), POSTK, 1, DESCR )
         CALL PZLACPY( 'A', KB, K-1, A, K+IA-1, JA, DESCA,
     $                 WORK( INDRT ), 1, 1, DESCRT )
*
         CALL PZLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), K, 1, DESCR )
         CALL PZTRSM( 'Right', 'L', 'N', 'N', NPK, KB, MONE, B, K+IB-1,
     $                K+JB-1, DESCB, WORK( INDG ), POSTK, 1, DESCG )
*
         CALL PZHEMM( 'Right', 'L', NPK, KB, ONEHALF, A, K+IA-1, K+JA-1,
     $                DESCA, WORK( INDG ), POSTK, 1, DESCG, ONE,
     $                WORK( INDR ), POSTK, 1, DESCR )
*
         CALL PZHER2K( 'Lower', 'No T', NPK, KB, ONE, WORK( INDG ),
     $                 POSTK, 1, DESCG, WORK( INDR ), POSTK, 1, DESCR,
     $                 RONE, A, POSTK+IA-1, POSTK+JA-1, DESCA )
*
         CALL PZGEMM( 'No T', 'No Conj', NPK, K-1, KB, ONE,
     $                WORK( INDG ), POSTK, 1, DESCG, WORK( INDRT ), 1,
     $                1, DESCRT, ONE, A, POSTK+IA-1, JA, DESCA )
*
         CALL PZHEMM( 'Right', 'L', NPK, KB, ONE, WORK( INDR ), K, 1,
     $                DESCR, WORK( INDG ), POSTK, 1, DESCG, ONE, A,
     $                POSTK+IA-1, K+JA-1, DESCA )
*
         CALL PZTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, K-1,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, K+IA-1, JA,
     $                DESCA )
*
         CALL PZLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDAA ), 1, 1, DESCAA )
*
         IF( MYROW.EQ.DESCAA( RSRC_ ) .AND. MYCOL.EQ.DESCAA( CSRC_ ) )
     $        THEN
            DO 20 I = 1, KB
               DO 10 J = 1, I
                  WORK( INDAA+J-1+( I-1 )*DESCAA( LLD_ ) )
     $               = DCONJG( WORK( INDAA+I-1+( J-1 )*
     $               DESCAA( LLD_ ) ) )
   10          CONTINUE
   20       CONTINUE
         END IF
*
         CALL PZTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1,
     $                1, DESCAA )
*
         CALL PZTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', KB, KB, ONE,
     $                B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1, 1,
     $                DESCAA )
*
         CALL PZLACPY( 'L', KB, KB, WORK( INDAA ), 1, 1, DESCAA, A,
     $                 K+IA-1, K+JA-1, DESCA )
*
         CALL PZTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', NPK, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, POSTK+IA-1,
     $                K+JA-1, DESCA )
*
         DESCR( CSRC_ ) = MOD( DESCR( CSRC_ )+1, NPCOL )
         DESCG( CSRC_ ) = MOD( DESCG( CSRC_ )+1, NPCOL )
         DESCRT( RSRC_ ) = MOD( DESCRT( RSRC_ )+1, NPROW )
         DESCAA( RSRC_ ) = MOD( DESCAA( RSRC_ )+1, NPROW )
         DESCAA( CSRC_ ) = MOD( DESCAA( CSRC_ )+1, NPCOL )
   30 CONTINUE
*
      WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
*
      RETURN
      END
*
#ifndef HAS_DLAMC3
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
***********************************************************
#endif
