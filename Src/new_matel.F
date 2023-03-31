! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!
      module m_new_matel
      public :: new_matel
      CONTAINS
      SUBROUTINE new_MATEL( OPERAT, IG1, IG2, R12, S12, DSDR )
C *******************************************************************
C Finds two-center matrix elements between 'atomic orbitals' 
C with finite radial and angular momentum cutoffs.
C Written by J.M.Soler. April 1995.
C Matrix elements of the position operator by DSP. June, 1999
C Electrostatic interaction added by JMS. July 2002.
C Introduction of unified 'global' indexes by AG. July 2011.
C ************************* INPUT ***********************************
C CHARACTER OPERAT : Operator to be used. The valid options are:
C   'S' => Unity (overlap). Uppercase required for all values.
C   'T' => -Laplacian
C   'U' => 1/|r'-r| (with evaluate returning charge distributions)
C   'X' => x, returning <phi1(r-R12)|x|phi2(r)> (origin on second atom)
C   'Y' => y, returning <phi1(r-R12)|y|phi2(r)>
C   'Z' => z, returning <phi1(r-R12)|z|phi2(r)>
C INTEGER IG1    : Global index of 1st function (must be positive)
C INTEGER IG2    : Gloabal index of 2nd function
C                    Indexes IG1, IG2 are used only to call 
C                    routines LCUT, RCUT and EVALUATE (see below), and 
C                    may have other meanings within those routines
C REAL*8  R12(3) : Vector from first to second atom
C ************************* OUTPUT **********************************
C REAL*8 S12      : Matrix element between orbitals.
C REAL*8 DSDR(3)  : Derivative (gradient) of S12 with respect to R12.
C ************************* ROUTINES CALLED *************************
C The following functions must exist:
C
C INTEGER FUNCTION LCUT(IG)
C   Returns the maximum angular momentum of the functions
C Input:
C   INTEGER IG : Global function index
C
C REAL*8 FUNCTION RCUT(IG)
C   Returns cutoff radius of the function
C Input:
C   INTEGER IG : Global function index
C
C SUBROUTINE EVALUATE(IG,R,PHI,GRPHI)
C   Returns the value of the functions to be integrated
C Input:
C   INTEGER IG : Global function index
C   REAL*8  R(3) : Position with respect to atom
C Output:
C   REAL*8  PHI      : Value of orbital at point R
C   REAL*8  GRPHI(3) : Gradient of PHI at point R  
C ************************* UNITS ***********************************
C Length units are arbitrary, but must be consistent in MATEL, RCUT
C   and EVALUATE. The laplacian unit is (length unit)**(-2).
C ************************* BEHAVIOUR *******************************
C 1) If |R12| > RCUT(IS1,IO1) + RCUT(IS2,IO2), returns U(Rmax)*Rmax/R
C    for OPERAT='U', and exactly zero in all other cases.
C 2) If (IS1.LE.0 .OR. IS2.LE.0) all internal tables are erased for
C    reinitialization and nothing is calculated.
C *******************************************************************
C    6  10        20        30        40        50        60        7072

C Modules -----------------------------------------------------------
      use precision, only : dp
      use parallel, only: node
      USE ALLOC
      USE m_matel_registry, ONLY: EVALUATE, LCUT, RCUT
      USE m_matel_registry, ONLY: EVALUATE_X, EVALUATE_Y, EVALUATE_Z
      use interpolation, only: spline, splint
      use m_errorf, only: derf
      use spher_harm, only: rlylm, ylmexp, ylmylm, lofilm
      use spher_harm, only: reset_spher_harm
      use sys,     only: die
      use m_radfft
      use siesta_options, only: matel_NRTAB
C -------------------------------------------------------------------

C Argument types and dimensions -------------------------------------
      IMPLICIT          NONE
      CHARACTER         OPERAT
      INTEGER           IG1, IG2
      real(dp)          DSDR(3), R12(3), S12
C -------------------------------------------------------------------

C Internal precision parameters  ------------------------------------
C  NRTAB is the number of radial points for matrix-element tables.
C  NQ is the number of radial points in reciprocal space.
C  EXPAND is a factor to expand some array sizes
C  Q2CUT is the required planewave cutoff for orbital expansion
C    (in Ry if lengths are in Bohr).
C  GWBYDR is the width of a gaussian used to neutralize the charge
C    distributions for operat='U', in units of the radial interval
C  FFTOL is the tolerance for considering equal the radial part of
C    two orbitals.
C  TINY is a small number to avoid a zero denominator
      INTEGER, save :: NRTAB = 1024
      INTEGER,          PARAMETER :: NQ     =  1024
      real(dp),         PARAMETER :: EXPAND =  1.20_dp
      real(dp),         PARAMETER :: Q2CUT  =  2.50e3_dp
      real(dp),         PARAMETER :: GWBYDR =  1.5_dp
      real(dp),         PARAMETER :: FFTOL  =  1.e-8_dp
      real(dp),         PARAMETER :: TINY   =  1.e-12_dp
      CHARACTER(LEN=*), PARAMETER :: MYNAME =  'MATEL'
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------
      INTEGER ::
     &  I, IF1, IF2, IFF, IFFY, IFLM1, IFLM2, 
     &  IG, IOPER, IQ, IR, IX,
     &  JF1, JF2, JFF, JFFR, JFFY, JFLM1, JFLM2, JLM, 
     &  JG1, JG2, JR,
     &  L, L1, L2, L3, LMAX,
     &  NILM, NILM1, NILM2, NJLM1, NJLM2
      INTEGER, SAVE ::
     &  MF=0, MFF=0, MFFR=0, MFFY=0, 
     &  NF=0, NFF=0, NFFR=0, NFFY=0, NR=NQ

      INTEGER, POINTER, SAVE ::
     &  IFFR(:), ILM(:), ILMFF(:), INDF(:,:), INDFF(:,:,:),
     &  INDFFR(:), INDFFY(:), NLM(:,:)

      real(dp) ::
     &  C, CH(2), CPROP, DFFR0, DFFRMX, DSRDR, ERRF,
     &  FFL(0:NQ), FFQ(0:NQ), FR(0:NQ,2), FQ(0:NQ,2), GAUSS, 
     &  Q, R, SR, VR(0:NQ,2), VQ(0:NQ,2), X12(3)

      real(dp), SAVE ::
     &  DQ, DR, DRTAB, GWIDTH, PI, QMAX, RMAX

      real(dp), POINTER, SAVE ::
     &  CFFR(:), DYDR(:,:), F(:,:), FFR(:,:,:), FFY(:,:), Y(:)

      LOGICAL ::
     &  FAR, FOUND, PROPOR

      LOGICAL, SAVE ::
     &  NULLIFIED=.FALSE.

      TYPE(allocDefaults) ::
     &  OLDEFS

      EXTERNAL  PROPOR, TIMER
C -------------------------------------------------------------------
C Start time counter 
*     CALL TIMER( MYNAME, 1 )

C Set allocation defaults 
      CALL ALLOC_DEFAULT( OLD=OLDEFS, COPY=.TRUE., SHRINK=.FALSE. )

C Nullify pointers 
      IF (.NOT.NULLIFIED) THEN
        NULLIFY( IFFR, ILM, ILMFF, INDF, INDFF, INDFFR, INDFFY, NLM,
     &           CFFR, DYDR, F, FFR, FFY, Y )
        CALL RE_ALLOC( INDF, 1, 1, 1, 1, 'INDF', MYNAME )
        CALL RE_ALLOC( INDFFY, 0,MFF, 'INDFFY', MYNAME )
        INDFFY(0) = 0
        NULLIFIED = .TRUE.

!       Set only at specific times. If matel_nrtab is changed the new
!       value will not take effect unless there is a re-initialization,
!       to avoid corruption of the tables.
        NRTAB = matel_NRTAB

      ENDIF

C Check if tables must be re-initialized 
      IF ( IG1.LE.0 .OR. IG2.LE.0 ) THEN
        CALL RESET_SPHER_HARM( )
        CALL RESET_RADFFT( )
        CALL DE_ALLOC( IFFR, 'IFFR', MYNAME )
        CALL DE_ALLOC( ILM, 'ILM', MYNAME )
        CALL DE_ALLOC( ILMFF, 'ILMFF', MYNAME )
        CALL DE_ALLOC( INDF, 'INDF', MYNAME )
        CALL DE_ALLOC( INDFF, 'INDFF', MYNAME )
        CALL DE_ALLOC( INDFFR, 'INDFFR', MYNAME )
        CALL DE_ALLOC( INDFFY, 'INDFFY', MYNAME )
        CALL DE_ALLOC( NLM, 'NLM', MYNAME )
        CALL DE_ALLOC( CFFR, 'CFFR', MYNAME )
        CALL DE_ALLOC( DYDR, 'DYDR', MYNAME )
        CALL DE_ALLOC( F, 'F', MYNAME )
        CALL DE_ALLOC( FFR, 'FFR', MYNAME )
        CALL DE_ALLOC( FFY, 'FFY', MYNAME )
        CALL DE_ALLOC( Y, 'Y', MYNAME )
        MF   = 0
        MFF  = 0
        MFFR = 0
        MFFY = 0
        NF   = 0
        NFF  = 0
        NFFR = 0
        NFFY = 0
        NULLIFIED = .FALSE.
        GOTO 900
      ENDIF

C Check argument OPERAT 
      IF ( OPERAT .EQ. 'S' ) THEN
        IOPER = 1
      ELSEIF ( OPERAT .EQ. 'T' ) THEN
        IOPER = 2
      ELSEIF ( OPERAT .EQ. 'U' ) THEN
        IOPER = 3
      ELSEIF ( OPERAT .EQ. 'X' ) THEN
        IOPER = 4
      ELSEIF ( OPERAT .EQ. 'Y' ) THEN
        IOPER = 5
      ELSEIF ( OPERAT .EQ. 'Z' ) THEN
        IOPER = 6
      ELSE
        call die('MATEL: Invalid value of argument OPERAT')
      ENDIF

C Check size of orbital index table
      IF ( MAX(IG1,IG2).GT.UBOUND(INDF,1)   .OR.
     &     MAX(IOPER,3).GT.UBOUND(INDF,2) ) THEN
        CALL RE_ALLOC( INDF, 1, MAX(UBOUND(INDF,1),IG1,IG2),
     &                       1,MAX(UBOUND(INDF,2),IOPER,3),
     &                'INDF', MYNAME)
        CALL RE_ALLOC( NLM,  1,MAX(UBOUND(INDF,1),IG1,IG2),
     &                       1,MAX(UBOUND(INDF,2),IOPER,3),
     &                'NLM', MYNAME )
      ENDIF

C Find radial expansion of each orbital -----------------------------
      DO I = 1,2
        IF (I .EQ. 1) THEN
          IG = IG1
          FOUND = (INDF(IG,1) .NE. 0)
        ELSE
          IG = IG2
          FOUND = (INDF(IG,IOPER) .NE. 0)
        ENDIF
        IF (.NOT.FOUND) THEN
*         CALL TIMER( 'MATEL1', 1 )
          PI = 4._dp * ATAN(1._dp)
C         Factor 2 below is because we will expand the product of
C         two orbitals
          QMAX = 2._dp * SQRT( Q2CUT )
          DQ = QMAX / NQ
          DR = PI / QMAX
          NR = NQ
          RMAX = NR * DR
          DRTAB = RMAX / NRTAB
          IF ( RCUT(IG) .GT. RMAX ) THEN
            call die('MATEL: NQ too small for required cutoff.')
          ENDIF

C         Reallocate arrays
          L = LCUT( IG )
          NILM = (L+2)**2
          IF (NF+NILM .GT. MF) MF = EXPAND * (NF+NILM)
          CALL RE_ALLOC( F, 0,NQ, 1,MF, 'F', MYNAME )
          CALL RE_ALLOC( ILM,     1,MF, 'ILM', MYNAME )
          CALL RE_ALLOC( INDFF,   1,MF, 1,MF, 1, MAX(IOPER,3),
     &                   'INDFF', MYNAME )

C         Expand orbital in spherical harmonics
          IF ((I.EQ.1) .OR. (IOPER.LE.3)) THEN
            CALL YLMEXP( L, RLYLM, EVALUATE, IG, 0, NQ, RMAX,
     &                   NILM, ILM(NF+1:), F(:,NF+1:) )
            INDF(IG,1) = NF+1
            INDF(IG,2) = NF+1
            INDF(IG,3) = NF+1
            NLM(IG,1) = NILM
            NLM(IG,2) = NILM
            NLM(IG,3) = NILM
          ELSE
            IF(IOPER.EQ.4) THEN
              CALL YLMEXP( L+1, RLYLM, EVALUATE_X, IG, 0, NQ, RMAX,
     &                     NILM, ILM(NF+1:), F(:,NF+1:) )
            ELSEIF(IOPER.EQ.5) THEN
              CALL YLMEXP( L+1, RLYLM, EVALUATE_Y, IG, 0, NQ, RMAX,
     &                     NILM, ILM(NF+1:), F(:,NF+1:) )
            ELSE
              CALL YLMEXP( L+1, RLYLM, EVALUATE_Z, IG, 0, NQ, RMAX,
     &                     NILM, ILM(NF+1:), F(:,NF+1:) )
            ENDIF
            INDF(IG,IOPER) = NF+1
            NLM(IG,IOPER) = NILM
          ENDIF

C         Store orbital in k-space
          DO JLM = 1,NILM
            NF = NF + 1
            L = LOFILM( ILM(NF) )
            CALL RADFFT( L, NQ, RMAX, F(0:NQ,NF), F(0:NQ,NF) )
*           F(NQ,NF) = 0._dp
          ENDDO

*         CALL TIMER( 'MATEL1', 2 )
        ENDIF
      ENDDO
C -------------------------------------------------------------------

C Find radial expansion of overlap ----------------------------------
      IF1 = INDF(IG1,1)
      IF2 = INDF(IG2,IOPER)

      IF ( INDFF(IF1,IF2,IOPER) .EQ. 0 ) THEN
*       CALL TIMER( 'MATEL2', 1 )

        NILM1 = NLM(IG1,1)
        NILM2 = NLM(IG2,IOPER)
        DO IFLM1 = IF1,IF1+NILM1-1
        DO IFLM2 = IF2,IF2+NILM2-1

C         Check interaction range
          IF (RCUT(IG1)+RCUT(IG2) .GT. RMAX) THEN
            call die('MATEL: NQ too small for required cutoff.')
          ENDIF

C         Special case for operat='U'
          IF (OPERAT.EQ.'U') THEN

C           Check that charge distributions are spherical
            IF (ILM(IFLM1).NE.1 .OR. ILM(IFLM2).NE.1) THEN
              CALL DIE('matel: ERROR: operat=U for L=0 only')
            ENDIF

C           Find L=0, Q=0 components of charge distributions
C           The actual charges are these times Y00=1/sqrt(4*pi)
            CH(1) = (2._dp*PI)**1.5_dp * F(0,IFLM1)
            CH(2) = (2._dp*PI)**1.5_dp * F(0,IFLM2)

C           Add gaussian neutralizing charge and find potential
            GWIDTH = GWBYDR * DR
            DO IQ = 1,NQ
              Q = IQ * DQ
              GAUSS = EXP( -0.5_dp*(Q*GWIDTH)**2 )
              FQ(IQ,1) = F(IQ,IFLM1) - F(0,IFLM1) * GAUSS
              FQ(IQ,2) = F(IQ,IFLM2) - F(0,IFLM2) * GAUSS
              VQ(IQ,1) = FQ(IQ,1) * 4._dp*PI/(Q*Q)
              VQ(IQ,2) = FQ(IQ,2) * 4._dp*PI/(Q*Q)
            ENDDO
            FQ(0,1) = 0._dp
            FQ(0,2) = 0._dp
            VQ(0,1) = 0._dp
            VQ(0,2) = 0._dp

C           Return to real space
            L = 0
            CALL RADFFT( L, NQ, NQ*PI/RMAX, FQ(0:NQ,1), FR(0:NQ,1) )
            CALL RADFFT( L, NQ, NQ*PI/RMAX, FQ(0:NQ,2), FR(0:NQ,2) )
            CALL RADFFT( L, NQ, NQ*PI/RMAX, VQ(0:NQ,1), VR(0:NQ,1) )
            CALL RADFFT( L, NQ, NQ*PI/RMAX, VQ(0:NQ,2), VR(0:NQ,2) )

C           Subtract neutralizing charge
            DO IR = 0,NR
              R = IR * DR
              GAUSS = EXP(-0.5_dp*(R/GWIDTH)**2) / 
     &                (2._dp*PI*GWIDTH**2)**1.5_dp
              IF (IR.EQ.0) THEN
                ERRF = SQRT(2._dp/PI) / GWIDTH
              ELSE
                ERRF = DERF(SQRT(0.5_dp)*R/GWIDTH) / R
              ENDIF
              FR(IR,1) = FR(IR,1) + CH(1) * GAUSS
              FR(IR,2) = FR(IR,2) + CH(2) * GAUSS
              VR(IR,1) = VR(IR,1) + CH(1) * ERRF
              VR(IR,2) = VR(IR,2) + CH(2) * ERRF
            ENDDO

C           Select NRTAB out of NQ points
            DO IR = 0,NRTAB
              JR = IR * NR / NRTAB
              FR(IR,1) = FR(JR,1)
              FR(IR,2) = FR(JR,2)
              VR(IR,1) = VR(JR,1)
              VR(IR,2) = VR(JR,2)
            ENDDO

          ELSE
            DO IQ = 0,NQ
              FQ(IQ,1) = F(IQ,IFLM1)
              FQ(IQ,2) = F(IQ,IFLM2)
            ENDDO
          ENDIF

C         Find orbitals convolution by multiplication in k-space
          C = ( 2.0_dp * PI )**1.5_dp
          DO IQ = 0,NQ
            Q = IQ * DQ
            FFQ(IQ) = C * FQ(IQ,1) * FQ(IQ,2)
            IF ( OPERAT .EQ. 'T' ) THEN
              FFQ(IQ) = FFQ(IQ) * Q*Q
            ELSEIF (OPERAT.EQ.'U' .AND. IQ.GT.0) THEN
              FFQ(IQ) = FFQ(IQ) * 4._dp*PI/(Q*Q)
            ENDIF
          ENDDO

C         Loop on possible values of l quantum number of product
          L1 = LOFILM( ILM(IFLM1) )
          L2 = LOFILM( ILM(IFLM2) )
          CALL RE_ALLOC( IFFR, 0,L1+L2, 'IFFR', MYNAME )
          CALL RE_ALLOC( CFFR, 0,L1+L2, 'CFFR', MYNAME )
          DO L3 = ABS(L1-L2), L1+L2, 2

C           Return to real space
            CALL RADFFT( L3, NQ, NQ*PI/RMAX, FFQ, FFL )
*           FFL(NQ) = 0._dp
            IF (MOD(ABS(L1-L2-L3)/2,2) .NE. 0) THEN
              DO IR = 0,NR
                FFL(IR) = - FFL(IR)
              ENDDO
            ENDIF

C           Divide by R**L
            IF (L3 .NE. 0) THEN
              DO IR = 1,NR
                R = IR * DR
                FFL(IR) = FFL(IR) / R**L3
              ENDDO
C             Parabolic extrapolation to R=0
              FFL(0) = ( 4.0_dp * FFL(1) - FFL(2) ) / 3.0_dp
            ENDIF

C           Select NRTAB out of NR points
            IF (MOD(NR,NRTAB) .NE. 0)
     &        CALL DIE('matel ERROR: NQ must be multiple of NRTAB')
            DO IR = 0,NRTAB
              JR = IR * NR / NRTAB
              FFL(IR) = FFL(JR)
            ENDDO

C           Find if radial function is already in table
            FOUND = .FALSE.
            SEARCH: DO JG1 = 1,SIZE(INDF,1)
                    DO JG2 = 1,SIZE(INDF,1)
              JF1 = INDF(JG1,1)
              JF2 = INDF(JG2,IOPER)
              IF (JF1.NE.0 .AND. JF2.NE.0) THEN
                NJLM1 = NLM(JG1,1)
                NJLM2 = NLM(JG2,IOPER)
                DO JFLM1 = JF1,JF1+NJLM1-1
                DO JFLM2 = JF2,JF2+NJLM2-1
                  JFF = INDFF(JFLM1,JFLM2,IOPER)
                  IF (JFF .NE. 0) THEN
                    DO JFFY = INDFFY(JFF-1)+1, INDFFY(JFF)
                      JFFR = INDFFR(JFFY)
                      IF ( PROPOR(NRTAB,FFL(1),FFR(1,1,JFFR),
     &                            FFTOL,CPROP)           ) THEN
                        FOUND = .TRUE.
                        IFFR(L3) = JFFR
                        CFFR(L3) = CPROP
                        EXIT SEARCH
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
                ENDDO
              ENDIF
            ENDDO
            ENDDO SEARCH

C           Store new radial function
            IF (.NOT.FOUND) THEN
              NFFR = NFFR + 1
              IF (NFFR .GT. MFFR) THEN
                MFFR = EXPAND * NFFR
                CALL RE_ALLOC( FFR, 0,NRTAB, 1, 2, 1, MFFR, 'FFR',
     &                         MYNAME )
              ENDIF
              IFFR(L3) = NFFR
              CFFR(L3) = 1._dp
              DO IR = 0,NRTAB
                FFR(IR,1,NFFR) = FFL(IR)
              ENDDO

C             Add neutralizing-charge corrections
              IF (OPERAT .EQ. 'U') THEN
                DO IR = 0,NRTAB
                  IF (IR.EQ.0) THEN
                    ERRF = 1._dp / SQRT(PI) / GWIDTH
                  ELSE
                    R = IR * DRTAB
                    ERRF = DERF(0.5_dp*R/GWIDTH) / R
                  ENDIF
                  FFR(IR,1,NFFR) = 
     &                FFR(IR,1,NFFR) 
     &              - CH(1) * CH(2) * ERRF
     &              + CH(1) * VR(IR,2) + CH(2) * VR(IR,1)
     &              - 2._dp*PI*GWIDTH**2 * CH(1) * FR(IR,2)
     &              - 2._dp*PI*GWIDTH**2 * CH(2) * FR(IR,1)
                ENDDO
              ENDIF

C             Setup spline interpolation
!!            Force derivative, rather than second derivative, to zero
!!              DFFR0 = HUGE(1.0_dp)
              DFFR0 = 0.0_dp
              DFFRMX = 0.0_dp
              CALL SPLINE( RMAX/NRTAB, FFR(0:NRTAB,1,NFFR), NRTAB+1, 
     &                     DFFR0, DFFRMX, FFR(0:NRTAB,2,NFFR) )
            ENDIF
          ENDDO

C         Reallocate some arrays
          NILM = (L1+L2+1)**2
          IF (NFFY+NILM .GT. MFFY) THEN
            MFFY = EXPAND * (NFFY+NILM)
            CALL RE_ALLOC( FFY, 1, 1, 1,MFFY, 'FFY', MYNAME )
            CALL RE_ALLOC( ILMFF, 1,MFFY, 'ILMFF', MYNAME )
            CALL RE_ALLOC( INDFFR, 1,MFFY, 'INDFFR', MYNAME )
          ENDIF
          CALL RE_ALLOC( Y, 1, NILM, 'Y', MYNAME, .FALSE. )
          CALL RE_ALLOC( DYDR, 1, 3, 1, NILM, 'DYDR', MYNAME, .FALSE. )

C         Expand the product of two spherical harmonics (SH) also in SH
          CALL YLMEXP( L1+L2, RLYLM, YLMYLM, ILM(IFLM1), ILM(IFLM2),
     &                 1, 1, 1.0_dp, NILM, ILMFF(NFFY+1:),
     &                 FFY(:,NFFY+1:))

C         Loop on possible lm values of orbital product
          DO I = 1,NILM
            NFFY = NFFY + 1
            JLM = ILMFF(NFFY)
            L3 = LOFILM( JLM )
            INDFFR(NFFY) = IFFR(L3)
            FFY(1,NFFY) = FFY(1,NFFY) * CFFR(L3)
          ENDDO

          NFF = NFF + 1
          IF (NFF .GT. MFF) THEN
            MFF = EXPAND * NFF
            CALL RE_ALLOC( INDFFY, 0,MFF, 'INDFFY', MYNAME )
          ENDIF
          INDFFY(NFF) = NFFY
          INDFF(IFLM1,IFLM2,IOPER) = NFF
        ENDDO
        ENDDO

*       CALL TIMER( 'MATEL2', 2 )
      ENDIF
C End of initialization section -------------------------------------

C Find value of matrix element and its gradient ---------------------
*     CALL TIMER( 'MATEL3', 1 )

C     Initialize output
      S12 = 0.0_dp
      DSDR(1) = 0.0_dp
      DSDR(2) = 0.0_dp
      DSDR(3) = 0.0_dp

C     Avoid R12=0
      X12(1) = R12(1)
      X12(2) = R12(2)
      X12(3) = R12(3)
      R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      IF (R .LT. TINY) THEN
        X12(3) = TINY
        R = SQRT( X12(1)*X12(1) + X12(2)*X12(2) + X12(3)*X12(3) )
      ENDIF

C     Find if orbitals are far (out of range)
      IF (R .GT. RCUT(IG1)+RCUT(IG2) ) THEN
        FAR = .TRUE.
        IF (OPERAT.EQ.'U') THEN
          IF1 = INDF(IG1,1)
          IF2 = INDF(IG2,IOPER)
          IFF = INDFF(IF1,IF2,IOPER)
          IFFY = INDFFY(IFF)
          JFFR = INDFFR(IFFY)
          S12 = FFR(NRTAB,1,JFFR) * RMAX / R
          DO IX = 1,3
            DSDR(IX) = - S12 * R12(IX) / R**2
          ENDDO
        ENDIF
      ELSE
        FAR = .FALSE.
      ENDIF

C     Find spherical harmonics times R**L
      IF (.NOT.FAR) THEN
        IF1 = INDF(IG1,1)
        IF2 = INDF(IG2,IOPER)
        NILM1 = NLM(IG1,1)
        NILM2 = NLM(IG2,IOPER)
        DO IFLM1 = IF1,IF1+NILM1-1
        DO IFLM2 = IF2,IF2+NILM2-1
          IFF = INDFF(IFLM1,IFLM2,IOPER)
          LMAX = 0
          DO IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
            JLM = ILMFF(IFFY)
            LMAX = MAX( LMAX, LOFILM(JLM) )
          ENDDO
          CALL RLYLM( LMAX, X12, Y, DYDR )

C         Interpolate radial functions and obtain SH expansion
          DO IFFY = INDFFY(IFF-1)+1, INDFFY(IFF)
            JFFR = INDFFR(IFFY)
            CALL SPLINT( RMAX/NRTAB, FFR(0:NRTAB,1,JFFR),
     &                   FFR(0:NRTAB,2,JFFR), NRTAB+1, R, SR, DSRDR )
            JLM = ILMFF(IFFY)
            S12 = S12 + SR * FFY(1,IFFY) * Y(JLM)
            DO IX = 1,3
              DSDR(IX) = DSDR(IX) +
     &                   DSRDR * FFY(1,IFFY) * Y(JLM) * X12(IX) / R +
     &                   SR * FFY(1,IFFY) * DYDR(IX,JLM)
            ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDIF

*     CALL TIMER( 'MATEL3', 2 )
C -------------------------------------------------------------------

C Restore allocation defaults 
  900 CONTINUE
      CALL ALLOC_DEFAULT( RESTORE=OLDEFS )

C Stop time counter
*     CALL TIMER( MYNAME, 2 )

!------------------------ Internal procedures

      END SUBROUTINE new_MATEL
      end module m_new_matel
