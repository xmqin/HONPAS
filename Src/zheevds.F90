  subroutine zheevds( jobz, uplo, n, a, lda, z, ldz, w, work, lwork, &
                      rwork, lrwork, iwork, liwork, info )

  !  -- lapack driver routine (version 3.0) --
  !     univ. of tennessee, univ. of california berkeley, nag ltd.,
  !     courant institute, argonne national lab, and rice university
  !     june 30, 1999
  !     special siesta version to allow specific storage of eigenvectors
  !     jdg august 2004

  ! scalar arguments
  character       :: jobz
  character       :: uplo

  integer         :: info
  integer         :: lda
  integer         :: ldz
  integer         :: liwork
  integer         :: lrwork
  integer         :: lwork
  integer         :: n

  ! array arguments ..
  integer         :: iwork( * )
  double precision:: rwork( * ), w( * )
  complex*16      :: a( lda, * ), z( ldz, * ), work( * )


!  Purpose
!  =======
!
!  ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a
!  complex Hermitian matrix A.  If eigenvectors are desired, it uses a
!  divide and conquer algorithm.
!
!  The divide and conquer algorithm makes very mild assumptions about
!  floating point arithmetic. It will work on machines with a guard
!  digit in add/subtract, or on those binary machines without guard
!  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!  Cray-2. It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.
!          If N <= 1,                LWORK must be at least 1.
!          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
!          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  RWORK   (workspace/output) DOUBLE PRECISION array,
!                                         dimension (LRWORK)
!          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!
!  LRWORK  (input) INTEGER
!          The dimension of the array RWORK.
!          If N <= 1,                LRWORK must be at least 1.
!          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
!          If JOBZ  = 'V' and N > 1, LRWORK must be at least
!                         1 + 5*N + 2*N**2.
!
!          If LRWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the RWORK array,
!          returns this value as the first entry of the RWORK array, and
!          no error message related to LRWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          If N <= 1,                LIWORK must be at least 1.
!          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Jeff Rutter, Computer Science Division, University of California
!     at Berkeley, USA
!
!  =====================================================================
!
  !------------------------------------------------------------- Input Variables
  double precision   :: zero
  double precision   :: one
  parameter          ( zero = 0.0d0, one = 1.0d0 )
  complex*16         cone
  parameter          ( cone = ( 1.0d0, 0.0d0 ) )
  !------------------------------------------------------------- Local Variables
  logical            :: lower
  logical            :: lquery
  logical            :: wantz

  integer            :: iinfo
  integer            :: imax
  integer            :: inde
  integer            :: indrwk
  integer            :: indtau
  integer            :: indwk2
  integer            :: indwrk
  integer            :: iscale
  integer            :: liopt
  integer            :: liwmin
  integer            :: llrwk
  integer            :: llwork
  integer            :: llwrk2
  integer            :: lopt
  integer            :: lropt
  integer            :: lrwmin
  integer            :: lwmin

  double precision   :: anrm
  double precision   :: bignum
  double precision   :: eps
  double precision   :: rmax
  double precision   :: rmin
  double precision   :: safmin
  double precision   :: sigma
  double precision   :: smlnum

  ! external functions ..
  logical            :: lsame
  double precision   :: dlamch
  double precision   :: zlanhe

  external           :: lsame
  external           :: dlamch
  external           :: zlanhe

  ! external subroutines ..
  external           :: dscal
  external           :: dsterf
  external           :: xerbla
  external           :: zhetrd
  external           :: zlacpy
  external           :: zlascl
  external           :: zstedc
  external           :: zunmtr

  ! intrinsic functions ..
  intrinsic          dble, int, max, sqrt
  !--------------------------------------------------------------------- BEGIN
  ! executable statements test the input parameters.
  wantz = lsame( jobz, 'V' )
  lower = lsame( uplo, 'L' )
  lquery = ( lwork.eq.-1 .or. lrwork.eq.-1 .or. liwork.eq.-1 )
  info = 0
  IF( n .le. 1 ) THEN
    lwmin  = 1
    lrwmin = 1
    liwmin = 1
    lopt   = lwmin
    lropt  = lrwmin
    liopt  = liwmin
  ELSE
    if( wantz ) then
      lwmin  = 2*n 
      lrwmin = 1 + 5*n + 2*n**2
      liwmin = 3 + 5*n
    else
      lwmin  = n + 1
      lrwmin = n
      liwmin = 1
    end if
    lopt  = lwmin
    lropt = lrwmin
    liopt = liwmin
  ENDIF

  if( .not.( wantz .or. lsame( jobz, 'N' ) ) ) then
    info = -1
  else if( .not.( lower .or. lsame( uplo, 'U' ) ) ) then
    info = -2
  else if( n.lt.0 ) then
    info = -3
  else if( lda.lt.max( 1, n ) ) then
    info = -5
  else if( lwork.lt.lwmin .and. .not.lquery ) then
    info = -8
  else if( lrwork.lt.lrwmin .and. .not.lquery ) then
    info = -10
  else if( liwork.lt.liwmin .and. .not.lquery ) then
    info = -12
  end if
  if( info.eq.0 ) then
    work( 1 ) = lopt
    rwork( 1 ) = lropt
    iwork( 1 ) = liopt
  end if
  if( info.ne.0 ) then
    call xerbla( 'ZHEEVDS', -info )
    return
  else if( lquery ) then
    return
  end if
  ! quick return if possible
  if( n.eq.0 ) then
    return
  endif
  if( n.eq.1 ) then
    w( 1 ) = a( 1, 1 )
    if( wantz ) then
      a( 1, 1 ) = cone
    endif
    return
  end if

  ! get machine constants.
  safmin = dlamch( 'Safe minimum' )
  eps = dlamch( 'Precision' )
  smlnum = safmin / eps
  bignum = one / smlnum
  rmin = sqrt( smlnum )
  rmax = sqrt( bignum )

  ! scale matrix to allowable range, if necessary.
  anrm = zlanhe( 'M', uplo, n, a, lda, rwork )
  iscale = 0
  if( anrm.gt.zero .and. anrm.lt.rmin ) then
    iscale = 1
    sigma = rmin / anrm
  else if( anrm.gt.rmax ) then
    iscale = 1
    sigma = rmax / anrm
  end if
  if( iscale.eq.1 ) then
    call zlascl( uplo, 0, 0, one, sigma, n, n, a, lda, info )
  endif

  ! call zhetrd to reduce hermitian matrix to tridiagonal form.
  inde   = 1
  indtau = 1
  indwrk = indtau + n
  indrwk = inde + n
  indwk2 = indwrk
  llwork = lwork - indwrk + 1
  llwrk2 = lwork - indwk2 + 1
  llrwk  = lrwork - indrwk + 1
  call zhetrd( uplo, n, a, lda, w, rwork( inde ), work( indtau ), &
               z, llwork, iinfo )
  lopt = max( dble( lopt ), dble( n )+dble( z(1,1) ) )

  ! for eigenvalues only, call dsterf.  for eigenvectors, first call
  ! zstedc to generate the eigenvector matrix, work(indwrk), of the
  ! tridiagonal matrix, then call zunmtr to multiply it to the
  ! householder transformations represented as householder vectors in a.
  if( .not.wantz ) then
    call dsterf( n, w, rwork( inde ), info )
  else
    call zstedc( 'I', n, w, rwork( inde ), z, ldz,               &
                 work( indwk2 ), llwrk2, rwork( indrwk ), llrwk, &
                 iwork, liwork, info )

    call zunmtr( 'L', uplo, 'N', n, n, a, lda, work( indtau ), &
                 z, ldz, work( indwk2 ), llwrk2, iinfo )
    lopt = max( lopt, n+n**2+int( work( indwk2 ) ) )
  end if
  ! if matrix was scaled, then rescale eigenvalues appropriately.
  if( iscale.eq.1 ) then
    if( info.eq.0 ) then
      imax = n
    else
      imax = info - 1
    end if
    call dscal( imax, one / sigma, w, 1 )
  end if
  work( 1 )  = lopt
  rwork( 1 ) = lropt
  iwork( 1 ) = liopt
  return

  !----------------------------------------------------------------------- END
  end subroutine zheevds
