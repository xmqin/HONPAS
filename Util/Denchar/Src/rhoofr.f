
      SUBROUTINE RHOOFR( NA, NO, NUO, MAXND, MAXNA, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO, 
     .                   XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DATM, 
     .                   IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .                   DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD, ISCALE, RMAXO )
C **********************************************************************
C Compute the density of charge at the points of a plane or a 3D grid
C in real space
C Coded by J. Junquera November'98
C Modified by P. Ordejon to include 3D capabilities, June 2003
C **********************************************************************

      use precision
      USE FDF
      USE ATMFUNCS
      USE CHEMICAL
      USE LISTSC_MODULE, ONLY: LISTSC
      use planed, only: plane

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NUO, IOPTION, NPX, NPY, NPZ, ISCALE, IUNITCD,
     .  IDIMEN, MAXND, NSPIN, MAXNA,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  NUMD(NUO), LISTDPTR(NUO), LISTD(MAXND)

      real(dp), INTENT(IN) ::
     .  XA(3,NA), XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  RMAXO

      real(dp), INTENT(IN) ::
     . CELL(3,3), DSCF(MAXND,NSPIN),
     . DATM(NO)

C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NUO              : Total number of orbitals in Unit Cell
C INTEGER MAXND            : Maximum number
C                            of basis orbitals interacting, either directly
C                            or through a KB projector, with any orbital
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
C INTEGER NUMD(NUO)        : Control vector of Density Matrix
C                            (number of non-zero elements of eaach row)
C INTEGER LISTDPTR(NUO)    : Pointer to where each row of listh starts - 1
C                            The reason for pointing to the element before
C                            the first one is so that when looping over the
C                            elements of a row there is no need to shift by
C                            minus one.
C INTEGER LISTD(MAXND)     : Control vector of Density Matrix
C                            (list of non-zero elements of each row)
C REAL*8  DSCF(MAXND,NSPIN): Density Matrix
C REAL*8  DATM(NO)          : Neutral atom charge density of each orbital
C INTEGER IDIMEN           : Specify if the run is to plot quantities
C                            in a plane or in a 3D grid (2 or 3, respect)
C INTEGER IOPTION          : Option to generate the plane
C                            1 = Normal vector
C                            2 = Two vectors contained in the plane
C                            3 = Three points in the plane
C INTEGER NPX,NPY          : Number of points generated along x and y 
C                            direction in a system of reference in which
C                            the third components od the points of the plane is
C                            zero (Plane Reference Frame; PRF)
C REAL*8  XMIN, XMAX       : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX       : Limits of the plane in the PRF for y-direction
C REAL*8  NORMAL(3)        : Components of the normal vector used to define 
C                            the plane
C REAL*8  DIRVER1(3)       : Components of the first vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  DIRVER2(3)       : Components of the second vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only used if ioption = 3)
C INTEGER IUNITCD          : Unit of the charge density
C INTEGER ISCALE           : Unit if the points of the plane
C REAL*8  ARMUNI           : Conversion factor for the charge density
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER
     .  NPLAMAX, NAPLA, NAINCELL

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  INDICES, JNA

      REAL, DIMENSION(:,:), allocatable :: RHO
      REAL, DIMENSION(:), allocatable   :: DRHO

      real(dp), DIMENSION(:), ALLOCATABLE ::
     .   R2IJ, DISCF, DIATM, DIUP, DIDOWN

      real(dp), DIMENSION(:,:), ALLOCATABLE ::
     .   PLAPO, POINRE, XAPLA, XIJ, XAINCELL

      INTEGER
     .  NPO, IA, ISEL, NNA, UNIT1, UNIT2, UNIT3, UNIT4
 
      INTEGER
     .  I, J, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND, IX, IY, IZ, NX, NY, NZ, IZA(NA)

      real(dp)
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, PHINU, GRPHIMU(3), GRPHINU(3), 
     .  OCELL(3,3)

      real(dp)
     .  DENCHAR, DENCHAR0, DENUP, DENDOWN

      LOGICAL FIRST

      CHARACTER
     .  SNAME*30, FNAMESCF*38, FNAMEDEL*38, FNAMEUP*38, FNAMEDOWN*38,
     .  PASTE*38

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE
     .  NEIGHB, WROUT

C **********************************************************************
C INTEGER NPLAMAX          : Maximum number of points in the plane
C REAL*8  PLAPO(NPLAMAX,3) : Coordinates of the points of the plane in PRF
C REAL*8  POINRE(NPLAMAX,3): Coordinates of the points of the plane in Lattice
C                            Reference Frame
C INTEGER NAPLA            : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)      : Indices of the atoms whose coordinates will 
C                            be rotated from the lattice reference frame 
C                            to the in-plane reference frame
C REAL*8  XAPLA(3,NA)      : Atomic coordinates in plane reference frame
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
C REAL*8  DENCHAR          : Self-Consistent Density of Charge at a given
C                            point in real space
C REAL*8  DENCHAR0         : Harris Density of Charge at a given point
C                            in real space
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************

C Allocate some variables ---------------------------------------------
      NPLAMAX = NPX * NPY * NPZ

      ALLOCATE(PLAPO(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'rhoofr')

      ALLOCATE(POINRE(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'rhoofr')

      ALLOCATE(INDICES(NA))
      CALL MEMORY('A','I',NA,'rhoofr')

      ALLOCATE(XAPLA(3,NA))
      CALL MEMORY('A','D',3*NA,'rhoofr')

      ALLOCATE(XAINCELL(3,NA))
      CALL MEMORY('A','D',3*NA,'rhoofr')

      ALLOCATE(DISCF(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIATM(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIUP(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

      ALLOCATE(DIDOWN(NO))
      CALL MEMORY('A','D',NO,'rhoofr')

C Build the plane ------------------------------------------------------
          CALL PLANE( NA, NPLAMAX, IDIMEN, IOPTION, 
     .                XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, 
     .                NPX, NPY, NPZ, COORPO, NORMAL, 
     .                DIRVER1, DIRVER2, 
     .                XA, NAPLA, INDICES, ISCALE,
     .                POINRE, PLAPO, XAPLA )   
C End building of the plane --------------------------------------------

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        CALL MEMORY('D','I',SIZE(JNA),'rhoofr')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'rhoofr')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'rhoofr')
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'rhoofr')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'rhoofr')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'rhoofr')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

C Open files to store charge density -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')

      IF (NSPIN .EQ. 1) THEN
        IF (IDIMEN .EQ. 2) THEN
          FNAMESCF = PASTE(SNAME,'.CON.SCF')
          FNAMEDEL = PASTE(SNAME,'.CON.DEL')
          CALL IO_ASSIGN(UNIT1)
          OPEN(UNIT = UNIT1, FILE = FNAMESCF, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT1)
          CALL IO_ASSIGN(UNIT2)
          OPEN(UNIT = UNIT2, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT2)
        ELSEIF (IDIMEN .EQ. 3) THEN
          FNAMESCF = PASTE(SNAME,'.RHO.cube')
          FNAMEDEL = PASTE(SNAME,'.DRHO.cube')
          CALL IO_ASSIGN(UNIT1)
          OPEN(UNIT = UNIT1, FILE = FNAMESCF, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT1)
          CALL IO_ASSIGN(UNIT2)
          OPEN(UNIT = UNIT2, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT2)
        ENDIF
      ELSEIF (NSPIN .EQ. 2) THEN
        IF (IDIMEN .EQ. 2) THEN
          FNAMESCF = PASTE(SNAME,'.CON.MAG' )
          FNAMEDEL = PASTE(SNAME,'.CON.DEL' )
          FNAMEUP  = PASTE(SNAME,'.CON.UP'  )
          FNAMEDOWN= PASTE(SNAME,'.CON.DOWN')
          CALL IO_ASSIGN(UNIT1)
          OPEN(UNIT = UNIT1, FILE = FNAMESCF, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT1)
          CALL IO_ASSIGN(UNIT2)
          OPEN(UNIT = UNIT2, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT2)
          CALL IO_ASSIGN(UNIT3)
          OPEN(UNIT = UNIT3, FILE = FNAMEUP, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT3)
          CALL IO_ASSIGN(UNIT4)
          OPEN(UNIT = UNIT4, FILE = FNAMEDOWN, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT4)
        ELSE IF (IDIMEN .EQ. 3) THEN
          FNAMEUP = PASTE(SNAME,'.RHO.UP.cube' )
          FNAMEDOWN = PASTE(SNAME,'.RHO.DOWN.cube' )
          FNAMEDEL = PASTE(SNAME,'.DRHO.cube' )
          CALL IO_ASSIGN(UNIT1)
          OPEN(UNIT = UNIT1, FILE = FNAMEUP, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT1)
          CALL IO_ASSIGN(UNIT2)
          OPEN(UNIT = UNIT2, FILE = FNAMEDOWN, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT2)
          CALL IO_ASSIGN(UNIT3)
          OPEN(UNIT = UNIT3, FILE = FNAMEDEL, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNIT3)
        ENDIF
      ELSE
        WRITE(6,*)'BAD NUMBER NSPIN IN RHOOFR.F'
        WRITE(6,*)'NSPIN = ',NSPIN
        WRITE(6,*)'IT MUST BE 1 => NON POLARIZED, OR 2 = > POLARIZED'
        STOP
      ENDIF

      IF (IDIMEN .EQ. 2) THEN
C Select all atoms in list to be printed out
        NAINCELL=NAPLA
        DO IA=1,NAPLA
          DO IX=1,3
            XAINCELL(IX,IA)=XAPLA(IX,IA)
          ENDDO
        ENDDO
      ELSE IF (IDIMEN .EQ.3) THEN
        DO IX = 1,3
          DO IY = 1,3
            OCELL(IX,IY)=0.D0
          ENDDO
        ENDDO
C   Determine cell size
        OCELL(1,1) = DABS(XMAX-XMIN)
        OCELL(2,2) = DABS(YMAX-YMIN)
        OCELL(3,3) = DABS(ZMAX-ZMIN)
C   Determine atoms which are within the plotting box
C   so that only they will be printed
        NAINCELL=0
        DO IA=1,NA
          IF ( (XAPLA(1,IA).LT.XMIN*1.1) .OR. (XAPLA(1,IA).GT.XMAX*1.1)
     .    .OR. (XAPLA(2,IA).LT.YMIN*1.1) .OR. (XAPLA(2,IA).GT.YMAX*1.1)
     .    .OR. (XAPLA(3,IA).LT.ZMIN*1.1) .OR. (XAPLA(3,IA).GT.ZMAX*1.1))
     .    GOTO 90
          NAINCELL=NAINCELL+1
          IZA(NAINCELL) = ATOMIC_NUMBER(ISA(IA))
          DO IX=1,3
            XAINCELL(IX,NAINCELL)=XAPLA(IX,IA)
          ENDDO
90        CONTINUE
        ENDDO
        IF (NSPIN .EQ. 1) THEN
          WRITE(UNIT1,*) FNAMESCF
          WRITE(UNIT1,*) FNAMESCF
          WRITE(UNIT1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
          WRITE(UNIT1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
          WRITE(UNIT1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
          WRITE(UNIT1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
          DO IA = 1,NAINCELL
          WRITE(UNIT1,'(i5,4f12.6)') IZA(IA),0.0,
     .                               (XAINCELL(IX,IA),IX=1,3)
          ENDDO
          WRITE(UNIT2,*) FNAMEDEL
          WRITE(UNIT2,*) FNAMEDEL
          WRITE(UNIT2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
          WRITE(UNIT2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
          WRITE(UNIT2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
          WRITE(UNIT2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
          DO IA = 1,NAINCELL
          WRITE(UNIT2,'(i5,4f12.6)') IZA(IA),0.0,
     .                               (XAINCELL(IX,IA),IX=1,3)
          ENDDO
        ELSE IF (NSPIN .EQ. 2) THEN
          WRITE(UNIT1,*) FNAMEUP
          WRITE(UNIT1,*) FNAMEUP
          WRITE(UNIT1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
          WRITE(UNIT1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
          WRITE(UNIT1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
          WRITE(UNIT1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
          DO IA = 1,NAINCELL
            WRITE(UNIT1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                 (XAINCELL(IX,IA),IX=1,3)
          ENDDO
          WRITE(UNIT2,*) FNAMEDOWN
          WRITE(UNIT2,*) FNAMEDOWN
          WRITE(UNIT2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
          WRITE(UNIT2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
          WRITE(UNIT2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
          WRITE(UNIT2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
          DO IA = 1,NAINCELL
            WRITE(UNIT2,'(i5,4f12.6)') IZA(IA),0.0,
     .                                 (XAINCELL(IX,IA),IX=1,3)
          ENDDO
          WRITE(UNIT3,*) FNAMEDEL
          WRITE(UNIT3,*) FNAMEDEL
          WRITE(UNIT3,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
          WRITE(UNIT3,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
          WRITE(UNIT3,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
          WRITE(UNIT3,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
          DO IA = 1,NAINCELL
            WRITE(UNIT3,'(i5,4f12.6)') IZA(IA),0.0,
     .                                 (XAINCELL(IX,IA),IX=1,3)
          ENDDO
        ENDIF
      ENDIF

      CALL WROUT(IDIMEN, .TRUE., .FALSE., IOPTION, NORMAL, COORPO, 
     .             DIRVER1, DIRVER2,
     .             NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .             IUNITCD,  NA, NAINCELL, INDICES, XAINCELL )
   

      WRITE(6,'(A)')
      WRITE(6,'(A)')
     .    '   Generating Charge Density values on the grid now....'
      WRITE(6,'(A)')
     .    '   Please, wait....'
      WRITE(6,'(A)')


C Allocate space for density in 3D-grid

       ALLOCATE(RHO(NPX*NPY*NPZ,NSPIN))
       CALL MEMORY('A','S',NPX*NPY*NPZ*NSPIN,'RHOOFR')
       ALLOCATE(DRHO(NPX*NPY*NPZ))
       CALL MEMORY('A','S',NPX*NPY*NPZ,'RHOOFR')

      
C Loop over all points in real space -----------------------------------
C      DO 100 NPO = 1, NPX*NPY*NPZ
      NPO = 0
      DO 102 NZ = 1,NPZ
      DO 101 NY = 1,NPY
      DO 100 NX = 1,NPX
        NPO = NPO + 1

C Initialize the density of charge at each point -----------------------
        DENCHAR  = 0.D0
        DENCHAR0 = 0.D0
        DENUP    = 0.D0
        DENDOWN  = 0.D0

C Localize non-zero orbitals at each point in real space ---------------
        DO IX = 1,3
          XPO(IX) = POINRE(NPO,IX)
        ENDDO
     
        IA   = 0
        ISEL = 0
        NNA  = MAXNA

        CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .               NNA, JNA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
         DO 110 IAT1 = 1, NNA
           IF( R2IJ(IAT1) .GT. RMAX2 ) EXIT

           IAVEC1   = JNA(IAT1)
           IS1      = ISA(IAVEC1)
           XVEC1(1) = -XIJ(1,IAT1)
           XVEC1(2) = -XIJ(2,IAT1)
           XVEC1(3) = -XIJ(3,IAT1)


           DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
             IPHI1 = IPHORB(IO)
             IUO   = INDXUO(IO)
             CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

C Copy full row IUO of Density Matrix to DI(j) --------------------------
             IF (IO . EQ. IUO) THEN
               DO 130 IN = 1, NUMD(IUO)
                 IND = IN + LISTDPTR(IUO)
                 J = LISTD(IND)
                 IF (NSPIN .EQ. 1) THEN
                   DISCF(J) = DSCF(IND,1)
                 ELSEIF (NSPIN .EQ. 2) THEN
                   DIUP(J)   = DSCF(IND,1)
                   DIDOWN(J) = DSCF(IND,2)
                 ENDIF
 130           ENDDO
             ELSE
               DO 135 IN = 1, NUMD(IUO)
                 IND = IN + LISTDPTR(IUO)
                 J = LISTSC( IO, IUO, LISTD(IND) )
                 IF (NSPIN .EQ. 1) THEN
                   DISCF(J) = DSCF(IND,1)
                 ELSEIF (NSPIN .EQ. 2) THEN
                   DIUP(J)   = DSCF(IND,1)
                   DIDOWN(J) = DSCF(IND,2)
                 ENDIF
 135           ENDDO
             ENDIF

             DENCHAR0 = DENCHAR0 + PHIMU*PHIMU*DATM(IUO)
C            WRITE(6,*) IO,DATM(IO)

C Loop again over neighbours -------------------------------------------
             DO 140 IAT2 = 1, NNA
               IAVEC2   = JNA(IAT2)
               IS2      = ISA(IAVEC2)
               XVEC2(1) = -XIJ(1,IAT2)
               XVEC2(2) = -XIJ(2,IAT2)
               XVEC2(3) = -XIJ(3,IAT2)
               
               DO 150 JO = LASTO(IAVEC2-1) + 1, LASTO(IAVEC2)
                 IPHI2 = IPHORB(JO)
                 CALL PHIATM( IS2, IPHI2, XVEC2, PHINU, GRPHINU )
                 
                 IF ( NSPIN .EQ. 1 ) THEN
                   DENCHAR  = DENCHAR  + PHINU*PHIMU*DISCF(JO)
                 ELSEIF (NSPIN .EQ. 2) THEN 
                   DENUP    = DENUP    + PHINU*PHIMU*DIUP(JO)
                   DENDOWN  = DENDOWN  + PHINU*PHIMU*DIDOWN(JO)
                 ENDIF
                 
 150           ENDDO

 140         ENDDO
 120       ENDDO
 110     ENDDO
 
         IF (IDIMEN .EQ. 2) THEN
           IF ( NSPIN .EQ. 1 ) THEN
             WRITE(UNIT1,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),DENCHAR*ARMUNI
             WRITE(UNIT2,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),(DENCHAR-DENCHAR0)*ARMUNI 
           ELSEIF ( NSPIN .EQ. 2 ) THEN
             WRITE(UNIT1,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),(DENUP-DENDOWN)*ARMUNI
             WRITE(UNIT2,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),
     .            (DENUP+DENDOWN-DENCHAR0)*ARMUNI
             WRITE(UNIT3,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),DENUP*ARMUNI
             WRITE(UNIT4,'(3F12.5)')
     .            PLAPO(NPO,1),PLAPO(NPO,2),DENDOWN*ARMUNI
           ENDIF
         ELSE IF (IDIMEN .EQ. 3) THEN
           IF (NSPIN .EQ. 1) THEN
             RHO(NPO,1) = DENCHAR*ARMUNI
             DRHO(NPO) = (DENCHAR-DENCHAR0)*ARMUNI
           ELSE IF (NSPIN .EQ.2) THEN
             RHO(NPO,1) = DENUP*ARMUNI
             RHO(NPO,2) = DENDOWN*ARMUNI
             DRHO(NPO) = (DENUP+DENDOWN-DENCHAR0)*ARMUNI
           ENDIF
         ENDIF



         IF (IDIMEN .EQ. 2) THEN
           IF ( MOD(NPO,NPX) .EQ. 0 ) THEN
CC-AG        Use * for single line separation,
CC-AG        since (/) gives two...
CC-AG             WRITE(UNIT1,'(/)')
CC-AG             WRITE(UNIT2,'(/)')
                  WRITE(UNIT1,*)
                  WRITE(UNIT2,*)
             IF ( NSPIN .EQ. 2 ) THEN
               WRITE(UNIT3,*)
               WRITE(UNIT4,*)
             ENDIF
           ENDIF
         ENDIF


C End x loop
 100  ENDDO  
C End y and z loops
 101  ENDDO  
 102  ENDDO  


      IF (IDIMEN .EQ. 3) THEN
        IF (NSPIN .EQ. 1) THEN
          DO NX=1,NPX
            DO NY=1,NPY
              WRITE(UNIT1,'(6e13.5)')
     .          (RHO(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
              WRITE(UNIT2,'(6e13.5)')
     .          (DRHO(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY),NZ=1,NPZ)
            ENDDO
          ENDDO
        ELSE IF (NSPIN .EQ. 2) THEN
          DO NX=1,NPX
            DO NY=1,NPY
              WRITE(UNIT1,'(6e13.5)')
     .          (RHO(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
              WRITE(UNIT2,'(6e13.5)')
     .          (RHO(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
              WRITE(UNIT3,'(6e13.5)')
     .          (DRHO(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY),NZ=1,NPZ)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      WRITE(6,'(A)')
      WRITE(6,'(A)')
     .    '   Your output files are:'
      IF (NSPIN .EQ. 1) THEN
        WRITE(6,'(A,A)') '   ',FNAMESCF
        WRITE(6,'(A,A)') '   ',FNAMEDEL
      ELSE IF (NSPIN .EQ. 2) THEN
        WRITE(6,'(A,A)') '   ',FNAMESCF
        WRITE(6,'(A,A)') '   ',FNAMEDEL
        WRITE(6,'(A,A)') '   ',FNAMEUP
        WRITE(6,'(A,A)') '   ',FNAMEDOWN
      ENDIF


      CALL IO_CLOSE(UNIT1)
      CALL IO_CLOSE(UNIT2)
      if (NSPIN .EQ. 2) CALL IO_CLOSE(UNIT3)
      IF (IDIMEN .EQ. 2 .AND. NSPIN .EQ. 2) CALL IO_CLOSE(UNIT4)
     
          
      RETURN    
      END
