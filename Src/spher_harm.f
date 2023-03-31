! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module spher_harm
      use precision
      use sys
      use alloc, only: re_alloc, de_alloc

      implicit none

      real(dp), pointer, private :: Y(:) => null()
      real(dp), pointer, private :: DYDR(:,:) => null()
      INTEGER,           private :: MAX_LM=-1

      public :: rlylm, ylmexp, ylmylm, lofilm
      public :: reset_spher_harm
      public :: gauleg !! for phirphi...
      interface ylmexp
        module procedure ylmexp_1, ylmexp_2
      end interface

      private

      CONTAINS

      subroutine rlylm( LMAX, R, RLY, GRLY )
      integer, intent(in)   :: lmax
      real(dp), intent(in)  :: r(3)
      real(dp), intent(out) :: rly(0:)
!      real(dp), intent(out) :: grly(3,0:)   !!! Not accepted...
      real(dp), intent(out) :: grly(1:,0:)

C FINDS REAL SPHERICAL HARMONICS MULTIPLIED BY R**L: RLY=R**L*YLM,
C AND THEIR GRADIENTS GRLY, AT POINT R:
C    YLM = C * PLM( COS(THETA) ) * SIN(M*PHI)   FOR   M <  0
C    YLM = C * PLM( COS(THETA) ) * COS(M*PHI)   FOR   M >= 0
C WITH (THETA,PHI) THE POLAR ANGLES OF R, C A POSITIVE NORMALIZATION
C CONSTANT AND PLM ASSOCIATED LEGENDRE POLYNOMIALS.
C THE ORDER OF THE Y'S IS THAT IMPLICIT IN THE NESTED LOOPS
C    DO L=0,LMAX
C      DO M=-L,L
C WITH A UNIFIED INDEX ILM=1,...,LMAX**2 INCREASING BY ONE UNIT IN THE
C  INNER LOOP.
C WRITTEN BY J.M.SOLER. AUG/96
C *********** INPUT ***************************************************
C INTEGER LMAX : Maximum angular momentum quantum number required.
C REAL*8  R(3) : Position at which Y and GY are required.
C *********** OUTPUT **************************************************
C REAL*8 RLY(LMAX*LMAX)    : Real spherical harmonics times r**l,
C                             at point R, as explained above.
C REAL*8 GRLY(3,LMAX*LMAX) : Gradient of the RLY's at point R.
C *********** UNITS ***************************************************
C Units of R are arbitrary. Units of RLY and GRLY are related to those
C  of R in the obvious way.
C *********************************************************************


C Dimension parameters for internal variables
      INTEGER MAXL, MAXLP1
      PARAMETER ( MAXL = 10, MAXLP1 = MAXL+1 )

C Other internal parameters
      REAL(DP) TINY, ZERO, HALF, ONE, TWO, THREE, SIX
      PARAMETER (TINY=1.e-14_dp, ZERO=0._dp, HALF=0.5_dp, ONE=1._dp,
     .           TWO=2._dp, THREE=3._dp, SIX=6._dp)

      real(dp), parameter :: pi = 3.14159265358979323846264338327950_dp
      real(dp), parameter :: pi4 = 4 * pi

C Internal variables
      INTEGER
     .  I, ILM, ILM0, L, LMXMX, M
      REAL(DP)
     .  C(0:MAXLP1*MAXLP1), COSM, COSMM1, COSPHI,
     .  ZP(0:MAXLP1,0:MAXLP1), FAC, GY(3),
     .  P(0:MAXLP1,0:MAXLP1),
     .  RL(-1:MAXL), RSIZE, RX, RY, RZ, RXY,
     .  SINM, SINMM1, SINPHI, YY
      SAVE LMXMX, C
      DATA LMXMX /-1/

C Evaluate normalization constants once and for all
      IF (LMAX.GT.LMXMX) THEN
         IF (LMAX.GT.MAXL) call die('YLM: MAXL too small')
         DO 20 L=0,LMAX
            ILM0=L*L+L
            DO 15 M=0,L
               FAC=(2*L+1)/pi4
               DO 10 I=L-M+1,L+M
                  FAC=FAC/I
   10          CONTINUE
               C(ILM0+M)=SQRT(FAC)
C              Next line because YY's are real combinations of M and -M
               IF (M.NE.0) C(ILM0+M)=C(ILM0+M)*SQRT(TWO)
               C(ILM0-M)=C(ILM0+M)
   15       CONTINUE
   20    CONTINUE
         LMXMX=LMAX
      ENDIF

C Initalize to zero
      DO 25 ILM = 0,(LMAX+1)*(LMAX+1)-1
        RLY(ILM) = ZERO
        GRLY(1,ILM) = ZERO
        GRLY(2,ILM) = ZERO
        GRLY(3,ILM) = ZERO
   25 CONTINUE

C Explicit formulas up to L=2
      IF (LMAX.LE.2) THEN

         RLY(0) = C(0)

C        Label 999 is the exit point
         IF (LMAX.EQ.0) GOTO 999

         RLY(1)    = -(C(1)*R(2))
         GRLY(2,1) = -C(1)

         RLY(2)    =  C(2)*R(3)
         GRLY(3,2) =  C(2)

         RLY(3)    = -(C(3)*R(1))
         GRLY(1,3) = -C(3)

         IF (LMAX.EQ.1) GOTO 999

         RLY(4)    =  C(4)*SIX*R(1)*R(2)
         GRLY(1,4) =  C(4)*SIX*R(2)
         GRLY(2,4) =  C(4)*SIX*R(1)

         RLY(5)    = (-C(5))*THREE*R(2)*R(3)
         GRLY(2,5) = (-C(5))*THREE*R(3)
         GRLY(3,5) = (-C(5))*THREE*R(2)

         RLY(6)    =  C(6)*HALF*(TWO*R(3)*R(3)-R(1)*R(1)-R(2)*R(2))
         GRLY(1,6) = (-C(6))*R(1)
         GRLY(2,6) = (-C(6))*R(2)
         GRLY(3,6) =  C(6)*TWO*R(3)

         RLY(7)    = (-C(7))*THREE*R(1)*R(3)
         GRLY(1,7) = (-C(7))*THREE*R(3)
         GRLY(3,7) = (-C(7))*THREE*R(1)

         RLY(8)    =  C(8)*THREE*(R(1)*R(1)-R(2)*R(2))
         GRLY(1,8) =  C(8)*SIX*R(1)
         GRLY(2,8) = (-C(8))*SIX*R(2)

         GOTO 999
      ENDIF

C Special case for R=0
      RSIZE = SQRT( R(1)*R(1)+R(2)*R(2)+R(3)*R(3) )
      IF ( RSIZE .LT. TINY ) THEN
        RLY(0) = C(0)
        GRLY(2,1) = -C(1)
        GRLY(3,2) =  C(2)
        GRLY(1,3) = -C(3)
        GOTO 999
      ENDIF

C Avoid z axis
      RX = R(1) / RSIZE
      RY = R(2) / RSIZE
      RZ = R(3) / RSIZE
      RXY = SQRT(RX*RX+RY*RY)
      IF (RXY.LT.TINY) THEN
        RX = TINY
        RXY = SQRT(RX*RX+RY*RY)
      ENDIF

C General algorithm based on routine PLGNDR of 'Numerical Recipes'
C See also J.M.Soler notes of 23/04/96.
C     Find associated Legendre polynomials and their derivative
      DO 60 M=LMAX,0,-1
         P(M,M+1)=ZERO
         P(M,M)=ONE
         FAC=ONE
         DO 30 I=1,M
            P(M,M)=-(P(M,M)*FAC*RXY)
            FAC=FAC+TWO
   30    CONTINUE
         P(M+1,M)=RZ*(2*M+1)*P(M,M)
         DO 40 L=M+2,LMAX
            P(L,M)=(RZ*(2*L-1)*P(L-1,M)-(L+M-1)*P(L-2,M))/(L-M)
   40    CONTINUE
         DO 50 L=M,LMAX
            ZP(L,M)=-((M*P(L,M)*RZ/RXY+P(L,M+1))/RXY)
   50   CONTINUE
   60 CONTINUE
C     Find spherical harmonics and their gradient
      RL(-1) = ZERO
      RL(0)  = ONE
      DO 70 L = 1,LMAX
        RL(L) = RL(L-1)*RSIZE
   70 CONTINUE
      COSPHI=RX/RXY
      SINPHI=RY/RXY
      COSM=ONE
      SINM=ZERO
      DO 90 M=0,LMAX
        DO 80 L=M,LMAX
          ILM=L*L+L-M
          GY(1)=(-ZP(L,M))*RX *RZ *SINM - P(L,M)*M*COSM*SINPHI/RXY
          GY(2)=(-ZP(L,M))*RY *RZ *SINM + P(L,M)*M*COSM*COSPHI/RXY
          GY(3)= ZP(L,M)*RXY*RXY*SINM
          YY=C(ILM)*P(L,M)*SINM
          RLY(ILM)=RL(L)*YY
          GRLY(1,ILM)= RX*L*RL(L-1)*YY + RL(L)*GY(1)*C(ILM)/RSIZE
          GRLY(2,ILM)= RY*L*RL(L-1)*YY + RL(L)*GY(2)*C(ILM)/RSIZE
          GRLY(3,ILM)= RZ*L*RL(L-1)*YY + RL(L)*GY(3)*C(ILM)/RSIZE
          
          ILM=L*L+L+M
          GY(1)=(-ZP(L,M))*RX *RZ *COSM + P(L,M)*M*SINM*SINPHI/RXY
          GY(2)=(-ZP(L,M))*RY *RZ *COSM - P(L,M)*M*SINM*COSPHI/RXY
          GY(3)= ZP(L,M)*RXY*RXY*COSM
          YY=C(ILM)*P(L,M)*COSM
          RLY(ILM)=RL(L)*YY
          GRLY(1,ILM)= RX*L*RL(L-1)*YY + RL(L)*GY(1)*C(ILM)/RSIZE
          GRLY(2,ILM)= RY*L*RL(L-1)*YY + RL(L)*GY(2)*C(ILM)/RSIZE
          GRLY(3,ILM)= RZ*L*RL(L-1)*YY + RL(L)*GY(3)*C(ILM)/RSIZE
   75     CONTINUE
   80   CONTINUE
        COSMM1=COSM
        SINMM1=SINM
        COSM=COSMM1*COSPHI-SINMM1*SINPHI
        SINM=COSMM1*SINPHI+SINMM1*COSPHI
   90 CONTINUE

  999 CONTINUE
      END SUBROUTINE RLYLM

      SUBROUTINE YLMYLM( ILM1, ILM2, R, YY, DYYDR )
      implicit none
      INTEGER, intent(in)       ::       ILM1, ILM2
      REAL(DP), intent(in)      ::       R(3)
      REAL(DP), intent(out)     ::       YY, DYYDR(3)

C Returns the product of two real spherical harmonics (SH),
C  times r**l each, and its gradient.
C *********** INPUT ***************************************************
C INTEGER ILM1, ILM2 : Combined SH indexes. ILM = L*L+L+M+1
C REAL*8  R(3)       : Vector towards the (theta,phi) direction.
C *********** OUTPUT **************************************************
C REAL*8  YY       : Product of the two SH.
C REAL*8  DYYDR(3) : Derivative (gradient) of YY with respect to R
C *********************************************************************
C Written by J.M.Soler. Feb' 96.
C *********************************************************************
      INTEGER           I, L, LM

      L  = MAX( LOFILM(ILM1), LOFILM(ILM2) )
      LM = (L+1)*(L+1)

      if (LM.gt.MAX_LM) then

        MAX_LM = LM
        call re_alloc( Y, 0, MAX_LM, 'Y', 'spher_harm' )
        call re_alloc( DYDR, 1, 3, 0, MAX_LM, 'DYDR', 'spher_harm' )
      endif

      CALL RLYLM( L, R, Y, DYDR )

      YY = Y(ILM1-1) * Y(ILM2-1)
      do I = 1,3
        DYYDR(I) = DYDR(I,ILM1-1)*Y(ILM2-1) + Y(ILM1-1)*DYDR(I,ILM2-1)
      enddo

      END SUBROUTINE YLMYLM

      INTEGER FUNCTION LOFILM( ILM )
      integer, intent(in) :: ilm

C Finds the total angular momentum quantum number L which corresponds
C to the combined index ILM == (L,M) implicit in the nested loop
C    ILM = 0
C    DO L=0,LMAX
C      DO M=-L,L
C        ILM = ILM + 1
C with ILM = 1,2,...,L**2
C Written by J.M.Soler. April 1996.

      INTEGER JLM, MAXL
      PARAMETER ( MAXL = 100 )

      if ( ILM .LE. 0 )  call die('LOFILM: ILM not allowed')
      JLM = 0
      DO 10 LOFILM = 0,MAXL
        JLM = JLM + 2*LOFILM + 1
        IF ( JLM .GE. ILM ) RETURN
   10 CONTINUE
      call die('LOFILM: ILM too large')
      END function lofilm

!
!     There are now two implementations of ylmexp, with
!     one and two indexes for "func"
!
      SUBROUTINE YLMEXP_2( LMAX, RLYLM_F, FUNC, IS, IO, IR1, NR, RMAX,
     .                   NY, ILM, FLM )
C Makes a radial times spherical-harmonic expansion of a function.

      integer, intent(in)          :: lmax
      interface
         subroutine rlylm_f(lmax,rvec,yy,grady)
         use precision, only: dp
         integer, intent(in)   :: lmax
         real(dp), intent(in)  :: rvec(3)
         real(dp), intent(out) :: yy(0:)
         real(dp), intent(out) :: grady(1:,0:)
         end subroutine rlylm_f

         subroutine func(i1,i2,rvec,yy,grady)
         use precision, only: dp
         integer, intent(in)   :: i1, i2
         real(dp), intent(in)  :: rvec(3)
         real(dp), intent(out) :: yy, grady(3)
         end subroutine func
       end interface

       integer, intent(in)        :: is, io
       integer, intent(in)        :: ir1, nr
       real(dp), intent(in)       :: rmax

       integer, intent(out)       :: ny
       integer, intent(out)       :: ilm(:)
       real(dp), dimension(ir1:,:), intent(out)  :: flm

C Written by J.M.Soler. September 1995.
C ************************* INPUT ***********************************
C INTEGER  LMAX                     : Max. ang. momentum quantum number
C EXTERNAL RLYLM_F(Lmax,Rvec,yy,GradY) : Returns spherical harmonics,
C                                     times r**l
C EXTERNAL FUNC(IS,IO,Rvec,F,GradF) : Function to be expanded.
C INTEGER  IS, IO                   : Indexes to call FUNC.
C INTEGER  IR1                      : First radial index.
C INTEGER  NR                       : Last radial index.
C REAL*8   RMAX                     : Maximum radius: R(IR)=RMAX*IR/NR
C ************************* OUTPUT **********************************
C INTEGER  NY            : Number of spherical harmonics required.
C INTEGER  ILM(NY)       : Spherical harmonics indexes: ILM=L*L+L+M+1.
C REAL*8   FLM(IR1:NR,NY): Radial expansion for each spherical harmonic:
C   F(Rvec) = Sum_iy( FLM(Rmod,iy) * YY(Rdir,ILM(iy)) ),
C   with   Rmod = RMAX*IR/NR
C ************************* UNITS ***********************************
C RMAX must be in the same units of the argument Rvec in FUNC.
C FLM is in the same units of the argument F in FUNC.
C ************************* BEHAVIOUR *******************************
C If function FUNC contains angular components with L higher than LMAX,
C   they will corrupt those computed with lower L.
C *******************************************************************


C Tolerance for FLM -------------------------------------------------
      REAL(DP) FTOL
      PARAMETER ( FTOL = 1.e-12_dp )
C -------------------------------------------------------------------

C Dimension parameters for internal variables -----------------------
      INTEGER MAXL, MAXLM, MAXSP
      PARAMETER ( MAXL  = 8 )
      PARAMETER ( MAXLM = (MAXL+1) * (MAXL+1) )
      PARAMETER ( MAXSP = (MAXL+1) * (2*MAXL+1) )
C -------------------------------------------------------------------

C Declare internal variables ----------------------------------------
      INTEGER
     .  IM, IR, ISP, IZ, JLM, JR, JY, NLM, NSP
      REAL(DP)
     .  DFDX(3), DYDX(3,MAXLM), 
     .  F(MAXSP), FY, PHI, PI, R, RX(3), THETA, 
     .  W(MAXL+1), WSP,
     .  X(3,MAXSP), YY(MAXLM), YSP(MAXSP,MAXLM),
     .  Z(MAXL+1)
      EXTERNAL CHKDIM
*     EXTERNAL TIMER
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER( 'YLMEXP', 1 )
C -------------------------------------------------------------------

C Check maximum angular momentum ------------------------------------
      CALL CHKDIM( 'YLMEXP', 'MAXL', MAXL, LMAX, 1 )
C -------------------------------------------------------------------

C Find special points and weights for gaussian quadrature -----------
      CALL GAULEG( -1._dp, 1._dp, Z, W, LMAX+1 )
C -------------------------------------------------------------------

C Find weighted spherical harmonics at special points ---------------
      PI = 4._dp * ATAN(1._dp)
      NLM = (LMAX+1)**2
      NSP = 0
      DO 30 IZ = 1,LMAX+1
        THETA = ACOS( Z(IZ) )
        DO 20 IM = 0,2*LMAX
          NSP = NSP + 1
          PHI = IM * 2._dp * PI / (2*LMAX+1)
          X(1,NSP) = SIN(THETA) * COS(PHI)
          X(2,NSP) = SIN(THETA) * SIN(PHI)
          X(3,NSP) = COS(THETA)
          CALL RLYLM_F( LMAX, X(1,NSP), YY, DYDX )
          WSP = W(IZ) * (2._dp*PI) / (2*LMAX+1)
          DO 10 JLM = 1,NLM
            YSP(NSP,JLM) = WSP * YY(JLM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C -------------------------------------------------------------------
      ILM = 0

C Expand FUNC in spherical harmonics at each radius -----------------
      NY = 0

C     Loop on radial points
      DO 80 IR = IR1,NR
        R = IR * RMAX / NR

C       Find function at points on a sphere of radius R
        DO 40 ISP = 1,NSP
          RX(1) = R * X(1,ISP)
          RX(2) = R * X(2,ISP)
          RX(3) = R * X(3,ISP)
          CALL FUNC( IS, IO, RX, F(ISP), DFDX )
   40   CONTINUE

C       Expand F(R) in spherical harmonics
        DO 70 JLM = 1,NLM
!          FY = DDOT(NSP,F,1,YSP(1,JLM),1)
          FY = dot_product(F(1:nsp),YSP(1:nsp,JLM))
          IF ( ABS(FY) .GT. FTOL ) THEN
C           Find JY corresponding to JLM
            DO 50 JY = 1,NY
              IF ( ILM(JY) .EQ. JLM ) THEN
                FLM(IR,JY) = FY
                GOTO 70
              ENDIF
   50       CONTINUE

C           New spherical harmonic required.
            NY = NY + 1
            ILM(NY) = JLM
            DO 60 JR = IR1,NR
              FLM(JR,NY) = 0._dp
   60       CONTINUE
            FLM(IR,NY) = FY
          ENDIF
   70   CONTINUE

   80 CONTINUE
C -------------------------------------------------------------------

C Special case for zero function ------------------------------------
      IF (NY .EQ. 0) THEN
        NY = 1
        ILM(1) = 1
        DO 90 IR = IR1,NR
          FLM(IR,1) = 0._dp
   90   CONTINUE
      ENDIF
C -------------------------------------------------------------------

C Stop time counter -------------------------------------------------
*     CALL TIMER( 'YLMEXP', 2 )
C -------------------------------------------------------------------

      end subroutine ylmexp_2
!
!     Version with a single global index to call func
!
      SUBROUTINE YLMEXP_1( LMAX, RLYLM_F, FUNC, IG, IR1, NR, RMAX,
     .                   NY, ILM, FLM )
C Makes a radial times spherical-harmonic expansion of a function.

      integer, intent(in)          :: lmax
      interface
         subroutine rlylm_f(lmax,rvec,yy,grady)
         use precision, only: dp
         integer, intent(in)   :: lmax
         real(dp), intent(in)  :: rvec(3)
         real(dp), intent(out) :: yy(0:)
         real(dp), intent(out) :: grady(1:,0:)
         end subroutine rlylm_f

         subroutine func(ig,rvec,yy,grady)
         use precision, only: dp
         integer, intent(in)   :: ig
         real(dp), intent(in)  :: rvec(3)
         real(dp), intent(out) :: yy, grady(3)
         end subroutine func
       end interface

       integer, intent(in)        :: ig
       integer, intent(in)        :: ir1, nr
       real(dp), intent(in)       :: rmax

       integer, intent(out)       :: ny
       integer, intent(out)       :: ilm(:)
       real(dp), dimension(ir1:,:), intent(out)  :: flm

C Written by J.M.Soler. September 1995.
C ************************* INPUT ***********************************
C INTEGER  LMAX                     : Max. ang. momentum quantum number
C EXTERNAL RLYLM_F(Lmax,Rvec,yy,GradY) : Returns spherical harmonics,
C                                     times r**l
C EXTERNAL FUNC(IS,IO,Rvec,F,GradF) : Function to be expanded.
C INTEGER  IG                       : Global Index to call FUNC.
C INTEGER  IR1                      : First radial index.
C INTEGER  NR                       : Last radial index.
C REAL*8   RMAX                     : Maximum radius: R(IR)=RMAX*IR/NR
C ************************* OUTPUT **********************************
C INTEGER  NY            : Number of spherical harmonics required.
C INTEGER  ILM(NY)       : Spherical harmonics indexes: ILM=L*L+L+M+1.
C REAL*8   FLM(IR1:NR,NY): Radial expansion for each spherical harmonic:
C   F(Rvec) = Sum_iy( FLM(Rmod,iy) * YY(Rdir,ILM(iy)) ),
C   with   Rmod = RMAX*IR/NR
C ************************* UNITS ***********************************
C RMAX must be in the same units of the argument Rvec in FUNC.
C FLM is in the same units of the argument F in FUNC.
C ************************* BEHAVIOUR *******************************
C If function FUNC contains angular components with L higher than LMAX,
C   they will corrupt those computed with lower L.
C *******************************************************************


C Tolerance for FLM -------------------------------------------------
      REAL(DP) FTOL
      PARAMETER ( FTOL = 1.e-12_dp )
C -------------------------------------------------------------------

C Dimension parameters for internal variables -----------------------
      INTEGER MAXL, MAXLM, MAXSP
      PARAMETER ( MAXL  = 8 )
      PARAMETER ( MAXLM = (MAXL+1) * (MAXL+1) )
      PARAMETER ( MAXSP = (MAXL+1) * (2*MAXL+1) )
C -------------------------------------------------------------------

C Declare internal variables ----------------------------------------
      INTEGER
     .  IM, IR, ISP, IZ, JLM, JR, JY, NLM, NSP
      REAL(DP)
     .  DFDX(3), DYDX(3,MAXLM), 
     .  F(MAXSP), FY, PHI, PI, R, RX(3), THETA, 
     .  W(MAXL+1), WSP,
     .  X(3,MAXSP), YY(MAXLM), YSP(MAXSP,MAXLM),
     .  Z(MAXL+1)
      EXTERNAL CHKDIM
*     EXTERNAL TIMER
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER( 'YLMEXP', 1 )
C -------------------------------------------------------------------

C Check maximum angular momentum ------------------------------------
      CALL CHKDIM( 'YLMEXP', 'MAXL', MAXL, LMAX, 1 )
C -------------------------------------------------------------------

C Find special points and weights for gaussian quadrature -----------
      CALL GAULEG( -1._dp, 1._dp, Z, W, LMAX+1 )
C -------------------------------------------------------------------

C Find weighted spherical harmonics at special points ---------------
      PI = 4._dp * ATAN(1._dp)
      NLM = (LMAX+1)**2
      NSP = 0
      DO 30 IZ = 1,LMAX+1
        THETA = ACOS( Z(IZ) )
        DO 20 IM = 0,2*LMAX
          NSP = NSP + 1
          PHI = IM * 2._dp * PI / (2*LMAX+1)
          X(1,NSP) = SIN(THETA) * COS(PHI)
          X(2,NSP) = SIN(THETA) * SIN(PHI)
          X(3,NSP) = COS(THETA)
          CALL RLYLM_F( LMAX, X(1,NSP), YY, DYDX )
          WSP = W(IZ) * (2._dp*PI) / (2*LMAX+1)
          DO 10 JLM = 1,NLM
            YSP(NSP,JLM) = WSP * YY(JLM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C -------------------------------------------------------------------
      ILM = 0

C Expand FUNC in spherical harmonics at each radius -----------------
      NY = 0

C     Loop on radial points
      DO 80 IR = IR1,NR
        R = IR * RMAX / NR

C       Find function at points on a sphere of radius R
        DO 40 ISP = 1,NSP
          RX(1) = R * X(1,ISP)
          RX(2) = R * X(2,ISP)
          RX(3) = R * X(3,ISP)
          CALL FUNC( IG, RX, F(ISP), DFDX )
   40   CONTINUE

C       Expand F(R) in spherical harmonics
        DO 70 JLM = 1,NLM
!          FY = DDOT(NSP,F,1,YSP(1,JLM),1)
          FY = dot_product(F(1:nsp),YSP(1:nsp,JLM))
          IF ( ABS(FY) .GT. FTOL ) THEN
C           Find JY corresponding to JLM
            DO 50 JY = 1,NY
              IF ( ILM(JY) .EQ. JLM ) THEN
                FLM(IR,JY) = FY
                GOTO 70
              ENDIF
   50       CONTINUE

C           New spherical harmonic required.
            NY = NY + 1
            ILM(NY) = JLM
            DO 60 JR = IR1,NR
              FLM(JR,NY) = 0._dp
   60       CONTINUE
            FLM(IR,NY) = FY
          ENDIF
   70   CONTINUE

   80 CONTINUE
C -------------------------------------------------------------------

C Special case for zero function ------------------------------------
      IF (NY .EQ. 0) THEN
        NY = 1
        ILM(1) = 1
        DO 90 IR = IR1,NR
          FLM(IR,1) = 0._dp
   90   CONTINUE
      ENDIF
C -------------------------------------------------------------------

C Stop time counter -------------------------------------------------
*     CALL TIMER( 'YLMEXP', 2 )
C -------------------------------------------------------------------

      end subroutine ylmexp_1

      subroutine gauleg(X1,X2,X,W,N)
!
!     Abscissas and weights for Gauss-Legendre quadrature formula
!
      integer, intent(in)   :: N
      real(dp), intent(in)  :: X1,X2
      real(dp), intent(out) :: X(N),W(N)

      real(dp), PARAMETER :: EPS=3.e-14_dp

      integer m, i, j
      real(dp) xm, xl, z, p1, p2, p3, pp, z1

      M = (N+1)/2 
      XM = 0.5_DP*(X2+X1)
      XL = 0.5_DP*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654_DP*(I-.25_DP)/(N+.5_DP))
1       CONTINUE
          P1=1._DP
          P2=0._DP
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2._DP*J-1._DP)*Z*P2-(J-1._DP)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1._DP)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        W(I)=2._DP*XL/((1._DP-Z*Z)*PP*PP)
        X(N+1-I)=XM+XL*Z
        W(N+1-I)=W(I)
12    CONTINUE

      end subroutine gauleg

      SUBROUTINE RESET_SPHER_HARM( )
      implicit none
      if (MAX_LM.gt.0) then
        MAX_LM = -1
        call de_alloc( Y,    'Y', 'spher_harm' )
        call de_alloc( DYDR, 'DYDR', 'spher_harm' )
      endif
      END SUBROUTINE RESET_SPHER_HARM

      end module spher_harm
