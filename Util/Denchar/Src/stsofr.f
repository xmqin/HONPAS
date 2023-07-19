
      SUBROUTINE STSOFR( NA, NO, NUO, MAXNA, NSPIN, 
     .                  ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .                  RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .                  ARMUNI, IUNITCD, RMAXO,
     .                  EMIN, NE, DELTAE, ETA, XSTS, BFUNC )
C **********************************************************************
C Computes the local density of states as a function of energy
C at a given point in real space, for the simulation of STS spectra
C Coded by P. Ordejon. July 2004
C **********************************************************************

      use precision
      USE FDF
      USE ATMFUNCS
      USE CHEMICAL
      USE LISTSC_MODULE, ONLY: LISTSC

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NUO, IUNITCD, NK,
     .  NSPIN, MAXNA, NWF(NK), NUMWF, 
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  INDW(NK,NUMWF), NE, BFUNC
      real(dp), INTENT(IN) ::
     .  XA(3,NA), ARMUNI, RMAXO,
     .  EMIN, DELTAE, ETA, XSTS(3)

      real(dp), INTENT(IN) ::
     . CELL(3,3), 
     . RPSI(NUO,NK,NUMWF,NSPIN), IPSI(NUO,NK,NUMWF,NSPIN),
     . E(NK,NUMWF,NSPIN), K(NK,3)

C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NUO              : Total number of orbitals in Unit Cell
C INTEGER MAXNA            : Maximum number of neighbours of any atom
C INTEGER NSPIN            : Number of different spin polarizations
C                            Nspin = 1 => unpolarized, Nspin = 2 => polarized
C INTEGER ISA(NA)          : Species index of each atom
C INTEGER IPHORB(NO)       : Orital index of each orbital in its atom
C INTEGER INDXUO(NO)       : Equivalent orbital in unit cell
C INTEGER LASTO(0:NA)      : Last orbital of each atom in array iphorb
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C                            (in bohr)
C REAL*8  CELL(3,3)        : Supercell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C REAL*8 RPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (real part)
C REAL*8 IPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (imag part)
C REAL*8 E(NK,NUMWF,NSPIN) : Eigen energies
C INTEGER INDW(NUMWF,NK)   : Index of the wavefunctions
C INTEGER NWF(NK)          : Number of wavefncts to print for each k-point
C INTEGER NUMWF            : Max num of wavefncts to print a given k-point
C INTEGER NK               : Number of k-points
C REAL*8 K(NK,3)           : k-points
C REAL*8  ARMUNI           : Conversion factor for the charge density
C INTEGER IUNITCD          : Unit of the charge density
C REAL*8  RMAXO            : Maximum range of basis orbitals
C REAL*8 EMIN              : Minimum energy in STS curve
C INTEGER NE               : Number of points in STS curve
C REAL*8 DELTAE            : Energy step between points in STS curve
C REAL*8 ETA               : Lorenzian or gaussian width
C REAL*8 XSTS(3)           : Position where STS is calculated
C INTEGER BFUNC            : Broadening function
C                            1 = Lorentzian   2 = Gaussian
C **********************************************************************

      INTEGER, DIMENSION(:), ALLOCATABLE :: JNA

      REAL, DIMENSION(:,:), allocatable :: RWF, IMWF, MWF, PWF

      real(dp), DIMENSION(:), ALLOCATABLE :: R2IJ, LDOS

      real(dp), DIMENSION(:,:), ALLOCATABLE ::  XIJ

      INTEGER
     .  IA, ISEL, NNA, UNIT, IE
 
      INTEGER
     .  I, J, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND, IX, IY, IZ, NX, NY, NZ, IWF, 
     .  INDWF, IZA(NA), IK

      real(dp)
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, GRPHIMU(3),
     .  OCELL(3,3), PHASE, SI, CO, PI, ENER, PSI2

      real(dp)
     .  RWAVE, RWAVEUP, RWAVEDN,
     .  IWAVE, IWAVEUP, IWAVEDN
 
      LOGICAL FIRST

      CHARACTER
     .  SNAME*40, FNAME*60, 
     .  PASTE*60, CHAR1*10, CHAR2*10, 
     .  EXT*20, EXT2*25

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE, NEIGHB

C **********************************************************************
C INTEGER IA               : Atom whose neighbours are needed.
C                            A routine initialization must be done by
C                            a first call with IA = 0
C                            If IA0=0, point X0 is used as origin instead
C INTEGER ISEL             : Single-counting switch (0=No, 1=Yes). If ISEL=1,
C                            only neighbours with JA.LE.IA are included in JNA
C INTEGER NNA              : Number of non-zero orbitals at a point in 
C                            real space
C INTEGER JNA(MAXNA)       : Atom index of neighbours. The neighbours
C                            atoms might be in the supercell
C REAL*8  XIJ(3,MAXNA)     : Vectors from point in real space to orbitals
C REAL*8  R2IJ(MAXNA)      : Squared distance to atomic orbitals
C REAL*8  XPO(3)           : Coordinates of the point of the plane respect
C                            we are going to calculate the neighbours orbitals
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************


C Allocate some variables ---------------------------------------------
      IF (NSPIN .NE. 1 .AND. NSPIN .NE. 2) THEN
        WRITE(6,*)'BAD NUMBER NSPIN IN STS.F'
        WRITE(6,*)'NSPIN = ',NSPIN
        WRITE(6,*)'IT MUST BE 1 => NON POLARIZED, OR 2 = > POLARIZED'
        STOP
      ENDIF

      PI = 4.0D0 * ATAN(1.0D0)

      ALLOCATE(LDOS(NE))
      CALL MEMORY('A','D',SIZE(LDOS),'stsofr')
C Initiallize ldos vector ---------------------------------------------

      DO IE = 1,NE
        LDOS(IE) = 0.0d0
      ENDDO

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        CALL MEMORY('D','I',SIZE(JNA),'stsofr')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'stsofr')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'stsofr')
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'stsofr')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'stsofr')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'stsofr')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2


C Assign the point where STS will be computed

      DO IX = 1,3
        XPO(IX) = XSTS(IX)
      ENDDO

C Open file to store sts spectrum   -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(SNAME,'.STS')
      CALL IO_ASSIGN(UNIT)
      OPEN(UNIT = UNIT, FILE = FNAME, STATUS = 'UNKNOWN',
     .     FORM = 'FORMATTED')
      REWIND(UNIT)

      WRITE(6,'(A)')
      WRITE(6,'(A)')
     .  '   Generating STS spectra at point'
      WRITE(6,'(3F10.5)')
     .  XPO(1),XPO(2),XPO(3)
      WRITE(6,'(A,I6)')
     .  '   Please, wait....'
      WRITE(6,'(A)')

C Localize non-zero orbitals each point in real space ---------------
     
      IA   = 0
      ISEL = 0
      NNA  = MAXNA

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )


C Loop over wavefunctions 

      DO IK  = 1, NK
      DO IWF = 1, NWF(IK)

        INDWF = INDW(IK,IWF)


C Initialize the value of the wave function  -----------------------
        RWAVE  = 0.D0
        RWAVEUP   = 0.D0
        RWAVEDN   = 0.D0
        IWAVE  = 0.D0
        IWAVEUP   = 0.D0
        IWAVEDN   = 0.D0

C Loop over Non-zero orbitals ------------------------------------------ 
        DO 110 IAT1 = 1, NNA
          IF( R2IJ(IAT1) .GT. RMAX2 ) EXIT

          IAVEC1   = JNA(IAT1)
          IS1      = ISA(IAVEC1)
          XVEC1(1) = -XIJ(1,IAT1)
          XVEC1(2) = -XIJ(2,IAT1)
          XVEC1(3) = -XIJ(3,IAT1)

C Calculate phase of Bloch wavefunctions

          PHASE = K(IK,1)*(XPO(1)+XIJ(1,IAT1))+
     .            K(IK,2)*(XPO(2)+XIJ(2,IAT1))+
     .            K(IK,3)*(XPO(3)+XIJ(3,IAT1))

          SI=SIN(PHASE)
          CO=COS(PHASE)

          DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
            IPHI1 = IPHORB(IO)
            IUO   = INDXUO(IO)
            CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

            IF ( NSPIN .EQ. 1 ) THEN
              RWAVE  = RWAVE  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,1)*CO - IPSI(IUO,IK,IWF,1)*SI)
              IWAVE  = IWAVE  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,1)*SI + IPSI(IUO,IK,IWF,1)*CO)
            ELSEIF (NSPIN .EQ. 2) THEN 
              RWAVEUP  = RWAVEUP  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,1)*CO - IPSI(IUO,IK,IWF,1)*SI)
              IWAVEUP  = IWAVEUP  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,1)*SI + IPSI(IUO,IK,IWF,1)*CO)
              RWAVEDN  = RWAVEDN  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,2)*CO - IPSI(IUO,IK,IWF,2)*SI)
              IWAVEDN  = IWAVEDN  + PHIMU * 
     .        (RPSI(IUO,IK,IWF,2)*SI + IPSI(IUO,IK,IWF,2)*CO)
            ENDIF



 120      ENDDO
 110    ENDDO


C Loop over energies to calculate the local density of states 

            ENER = EMIN

            DO IE = 1,NE

              IF (NSPIN .EQ. 1) THEN
                PSI2 = RWAVE**2+IWAVE**2
                IF (BFUNC .EQ. 1) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (ETA/PI) / ((ENER - E(IK,IWF,1))**2 + ETA**2)
                ELSE IF (BFUNC .EQ. 2) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (1.0d0/(ETA*DSQRT(PI))) * 
     .             DEXP(-((ENER - E(IK,IWF,1))/ETA)**2)
                ELSE 
                  STOP 'Wrong convolution function in stsofr.f'
                ENDIF
              ELSE IF (NSPIN .EQ. 2) THEN
                PSI2 = RWAVEUP**2+IWAVEUP**2
                IF (BFUNC .EQ. 1) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (ETA/PI) / ((ENER - E(IK,IWF,1))**2 + ETA**2)
                ELSE IF (BFUNC .EQ. 2) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (1.0d0/(ETA*DSQRT(PI))) * 
     .             DEXP(-((ENER - E(IK,IWF,1))/ETA)**2)
                ELSE 
                  STOP 'Wrong convolution function in stsofr.f'
                ENDIF
                PSI2 = RWAVEDN**2+IWAVEDN**2
                IF (BFUNC .EQ. 1) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (ETA/PI) / ((ENER - E(IK,IWF,2))**2 + ETA**2)
                ELSE IF (BFUNC .EQ. 2) THEN
                  LDOS(IE) = LDOS(IE) + (ARMUNI/NK) * 2.0 * PSI2 *
     .             (1.0d0/(ETA*DSQRT(PI))) * 
     .             DEXP(-((ENER - E(IK,IWF,2))/ETA)**2)
                ELSE 
                  STOP 'Wrong convolution function in stsofr.f'
                ENDIF
              ELSE 
                STOP 'Wrong NSPIN'
              ENDIF

              ENER = ENER + DELTAE
            ENDDO
           
      ENDDO
      ENDDO


C Print out LDOS -------------------------------------------------
      WRITE(UNIT,'(A)')  '# LDOS vs. Energy curve'
      WRITE(UNIT,'(A)')  '# '
      WRITE(UNIT,'(A,3e13.5,A)')
     .             '# LDOS computed on point ',(XSTS(I),I=1,3),' bohr'
      IF (IUNITCD .EQ. 1) THEN
        WRITE(UNIT,'(A)') '# in units of elecron/(bohr)**3/eV'
      ELSE IF (IUNITCD .EQ. 2) THEN
        WRITE(UNIT,'(A)') '# in units of elecron/(ang)**3/eV'
      ELSE IF (IUNITCD .EQ. 3) THEN
        WRITE(UNIT,'(A)') '# in units of elecron/unitcell/eV'
      ELSE 
        STOP 'Wrong IUNITCD in stsofr.f'
      ENDIF
      WRITE(UNIT,'(A)') '# '
      WRITE(UNIT,'(A)') '# Energy in eV)'
      WRITE(UNIT,'(A)') '# '
      WRITE(UNIT,'(A)') '#       Energy           LDOS'
      WRITE(UNIT,'(A)') '# '

      ENER = EMIN
      DO IE = 1,NE
        WRITE(UNIT,'(2f15.5)') ENER,LDOS(IE)
        ENER = ENER + DELTAE
      ENDDO

      WRITE(6,'(A)')
      WRITE(6,'(A)')
     .  '   Your STS output file is:'
      WRITE(6,'(A,A)') '   ',FNAME

      CALL IO_CLOSE(UNIT)


      CALL MEMORY('D','D',SIZE(LDOS),'stsofr')
      DEALLOCATE(LDOS)


      RETURN    
      END
