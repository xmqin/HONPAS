
      SUBROUTINE WAVOFR( NA, NO, NUO, MAXNA, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .                   RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .                   IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .                   DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD, ISCALE, RMAXO )
C **********************************************************************
C Compute the wave functions at the points of a plane or a 3D grid
C in real space
C Coded by P. Ordejon, from Junquera's rhoofr. July 2003
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
     .  IDIMEN, NSPIN, MAXNA, NK, NWF(NK), NUMWF, 
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  INDW(NK,NUMWF)
      real(dp), INTENT(IN) ::
     .  XA(3,NA), XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  RMAXO

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
C INTEGER IDIMEN           : Specify if the run is to plot quantities
C                            in a plane or in a 3D grid (2 or 3, respect)
C INTEGER IOPTION          : Option to generate the plane
C                            1 = Normal vector
C                            2 = Two vectors contained in the plane
C                            3 = Three points in the plane
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C REAL*8  XMIN, XMAX       : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX       : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only used if ioption = 3)
C REAL*8  NORMAL(3)        : Components of the normal vector used to define 
C                            the plane
C REAL*8  DIRVER1(3)       : Components of the first vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  DIRVER2(3)       : Components of the second vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  ARMUNI           : Conversion factor for the charge density
C INTEGER IUNITCD          : Unit of the charge density
C INTEGER ISCALE           : Unit if the points of the plane
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER
     .  NPLAMAX, NAPLA, NAINCELL

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  INDICES, JNA

      REAL, DIMENSION(:,:), allocatable :: RWF, IMWF, MWF, PWF

      real(dp), DIMENSION(:), ALLOCATABLE ::
     .   R2IJ

      real(dp), DIMENSION(:,:), ALLOCATABLE ::
     .   PLAPO, POINRE, XAPLA, XIJ, XAINCELL

      INTEGER
     .  NPO, IA, ISEL, NNA, UNITRE1, UNITRE2, UNITIM1, UNITIM2,
     .  UNITPH1, UNITPH2, UNITMO1, UNITMO2
 
      INTEGER
     .  I, J, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND, IX, IY, IZ, NX, NY, NZ, IWF, 
     .  INDWF, IZA(NA), IK

      real(dp)
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, GRPHIMU(3),
     .  OCELL(3,3), PHASE, SI, CO, PI

      real(dp)
     .  RWAVE, RWAVEUP, RWAVEDN,
     .  IWAVE, IWAVEUP, IWAVEDN,
     .  MWAVE, MWAVEUP, MWAVEDN,
     .  PWAVE, PWAVEUP, PWAVEDN
 
      LOGICAL FIRST

      CHARACTER
     .  SNAME*40, FNAMEWFRE*60, FNAMEWFIM*60, 
     .  FNAMEWFURE*60, FNAMEWFUIM*60, FNAMEWFDRE*60, FNAMEWFDIM*60, 
     .  FNAMEWFMO*60, FNAMEWFPH*60,
     .  FNAMEWFUMO*60, FNAMEWFUPH*60, FNAMEWFDMO*60, FNAMEWFDPH*60,
     .  PASTE*60, CHAR1*10, CHAR2*10, ITOCHAR*10, 
     .  EXT*20, EXT2*25

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE,
     .  NEIGHB, WROUT, ITOCHAR

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
C REAL RWF(NPO,NSPIN)      : Wave fnctn at each point of the grid (real part)
C REAL IMWF(NPO,NSPIN)     : Wave fnctn at each point of the grid (imag part)
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************

C Allocate some variables ---------------------------------------------


      PI = 4.0D0 * ATAN(1.0D0)

      NPLAMAX = NPX * NPY * NPZ

      ALLOCATE(PLAPO(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'wavofr')

      ALLOCATE(POINRE(NPLAMAX,3))
      CALL MEMORY('A','D',3*NPLAMAX,'wavofr')

      ALLOCATE(INDICES(NA))
      CALL MEMORY('A','I',NA,'wavofr')

      ALLOCATE(XAPLA(3,NA))
      CALL MEMORY('A','D',3*NA,'wavofr')

      ALLOCATE(XAINCELL(3,NA))
      CALL MEMORY('A','D',3*NA,'wavofr')

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
        CALL MEMORY('D','I',SIZE(JNA),'wavofr')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'wavofr')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'wavofr')
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'wavofr')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'wavofr')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'wavofr')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

C Loop over wavefunctions to write out

      DO IK  = 1, NK
      DO IWF = 1,NWF(IK)


        INDWF = INDW(IK,IWF)
        CHAR1 = ITOCHAR(INDWF)
        CHAR2 = ITOCHAR(IK)
C Open files to store wave functions -----------------------------------
        SNAME = FDF_STRING('SystemLabel','siesta')

        IF (NSPIN .EQ. 1) THEN
          IF (IDIMEN .EQ. 2) THEN
            FNAMEWFRE = PASTE(SNAME,'.CON.K')
            FNAMEWFRE = PASTE(FNAMEWFRE,CHAR2)
            FNAMEWFRE = PASTE(FNAMEWFRE,'.WF')
            FNAMEWFRE = PASTE(FNAMEWFRE,CHAR1)
            FNAMEWFPH = PASTE(FNAMEWFRE,'.PHASE')
            FNAMEWFMO = PASTE(FNAMEWFRE,'.MOD')
            FNAMEWFIM = PASTE(FNAMEWFRE,'.IMAG')
            FNAMEWFRE = PASTE(FNAMEWFRE,'.REAL')
          ELSEIF (IDIMEN .EQ. 3) THEN
            FNAMEWFRE = PASTE(SNAME,'.K')
            FNAMEWFRE = PASTE(FNAMEWFRE,CHAR2)
            FNAMEWFRE = PASTE(FNAMEWFRE,'.WF')
            FNAMEWFRE = PASTE(FNAMEWFRE,CHAR1)
            FNAMEWFPH = PASTE(FNAMEWFRE,'.PHASE')
            FNAMEWFPH = PASTE(FNAMEWFPH,'.cube')
            FNAMEWFMO = PASTE(FNAMEWFRE,'.MOD')
            FNAMEWFMO = PASTE(FNAMEWFMO,'.cube')
            FNAMEWFIM = PASTE(FNAMEWFRE,'.IMAG')
            FNAMEWFIM = PASTE(FNAMEWFIM,'.cube')
            FNAMEWFRE = PASTE(FNAMEWFRE,'.REAL')
            FNAMEWFRE = PASTE(FNAMEWFRE,'.cube')
          ENDIF
          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FNAMEWFRE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FNAMEWFIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FNAMEWFMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FNAMEWFPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)
        ELSEIF (NSPIN .EQ. 2) THEN
          IF (IDIMEN .EQ. 2) THEN
            FNAMEWFURE = PASTE(SNAME,'.CON.K')
            FNAMEWFURE = PASTE(FNAMEWFURE,CHAR2)
            FNAMEWFURE = PASTE(FNAMEWFURE,'.WF')
            FNAMEWFURE = PASTE(FNAMEWFURE,CHAR1)
            FNAMEWFDRE = PASTE(FNAMEWFURE,'.DOWN')
            FNAMEWFDPH = PASTE(FNAMEWFDRE,'.PHASE')
            FNAMEWFDMO = PASTE(FNAMEWFDRE,'.MOD')
            FNAMEWFDIM = PASTE(FNAMEWFDRE,'.IMAG')
            FNAMEWFDRE = PASTE(FNAMEWFDRE,'.REAL')
            FNAMEWFURE = PASTE(FNAMEWFURE,'.UP')
            FNAMEWFUPH = PASTE(FNAMEWFURE,'.PHASE')
            FNAMEWFUMO = PASTE(FNAMEWFURE,'.MOD')
            FNAMEWFUIM = PASTE(FNAMEWFURE,'.IMAG')
            FNAMEWFURE = PASTE(FNAMEWFURE,'.REAL')
          ELSE IF (IDIMEN .EQ. 3) THEN
            FNAMEWFURE = PASTE(SNAME,'.K')
            FNAMEWFURE = PASTE(FNAMEWFURE,CHAR2)
            FNAMEWFURE = PASTE(FNAMEWFURE,'.WF')
            FNAMEWFURE = PASTE(FNAMEWFURE,CHAR1)
            FNAMEWFDPH = PASTE(FNAMEWFURE,'.DOWN.PHASE')
            FNAMEWFDPH = PASTE(FNAMEWFDPH,'.cube')
            FNAMEWFDMO = PASTE(FNAMEWFURE,'.DOWN.MOD')
            FNAMEWFDMO = PASTE(FNAMEWFDMO,'.cube')
            FNAMEWFDIM = PASTE(FNAMEWFURE,'.DOWN.IMAG')
            FNAMEWFDIM = PASTE(FNAMEWFDIM,'.cube')
            FNAMEWFDRE = PASTE(FNAMEWFURE,'.DOWN.REAL')
            FNAMEWFDRE = PASTE(FNAMEWFDRE,'.cube')
            FNAMEWFUPH = PASTE(FNAMEWFURE,'.UP.PHASE')
            FNAMEWFUPH = PASTE(FNAMEWFUPH,'.cube')
            FNAMEWFUMO = PASTE(FNAMEWFURE,'.UP.MOD')
            FNAMEWFUMO = PASTE(FNAMEWFUMO,'.cube')
            FNAMEWFUIM = PASTE(FNAMEWFURE,'.UP.IMAG')
            FNAMEWFUIM = PASTE(FNAMEWFUIM,'.cube')
            FNAMEWFURE = PASTE(FNAMEWFURE,'.UP.REAL')
            FNAMEWFURE = PASTE(FNAMEWFURE,'.cube')
          ENDIF

          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FNAMEWFURE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITRE2)
          OPEN(UNIT = UNITRE2, FILE = FNAMEWFDRE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE2)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FNAMEWFUIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITIM2)
          OPEN(UNIT = UNITIM2, FILE = FNAMEWFDIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM2)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FNAMEWFUMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITMO2)
          OPEN(UNIT = UNITMO2, FILE = FNAMEWFDMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO2)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FNAMEWFUPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)
          CALL IO_ASSIGN(UNITPH2)
          OPEN(UNIT = UNITPH2, FILE = FNAMEWFDPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH2)
        ELSE
          WRITE(6,*)'BAD NUMBER NSPIN IN WAVOFR.F'
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
          NAINCELL=0
          DO IA=1,NA
            IF ((XAPLA(1,IA).LT.XMIN*1.1).OR.(XAPLA(1,IA).GT.XMAX*1.1)
     .     .OR. (XAPLA(2,IA).LT.YMIN*1.1).OR.(XAPLA(2,IA).GT.YMAX*1.1)
     .     .OR. (XAPLA(3,IA).LT.ZMIN*1.1).OR.(XAPLA(3,IA).GT.ZMAX*1.1))
     .      GOTO 90
            NAINCELL=NAINCELL+1
            IZA(NAINCELL) = ATOMIC_NUMBER(ISA(IA))
            DO IX=1,3
              XAINCELL(IX,NAINCELL)=XAPLA(IX,IA)
            ENDDO
90          CONTINUE
          ENDDO

          IF (NSPIN .EQ. 1) THEN
            WRITE(UNITRE1,*) FNAMEWFRE
            WRITE(UNITRE1,*) FNAMEWFRE
            WRITE(UNITRE1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITRE1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITRE1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITRE1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITRE1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITIM1,*) FNAMEWFIM
            WRITE(UNITIM1,*) FNAMEWFIM
            WRITE(UNITIM1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITIM1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITIM1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITIM1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITIM1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITMO1,*) FNAMEWFMO
            WRITE(UNITMO1,*) FNAMEWFMO
            WRITE(UNITMO1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITMO1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITMO1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITMO1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITMO1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITPH1,*) FNAMEWFPH
            WRITE(UNITPH1,*) FNAMEWFPH
            WRITE(UNITPH1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITPH1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITPH1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITPH1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITPH1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
          ELSE IF (NSPIN .EQ. 2) THEN
            WRITE(UNITRE1,*) FNAMEWFURE
            WRITE(UNITRE1,*) FNAMEWFURE
            WRITE(UNITRE1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITRE1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITRE1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITRE1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITRE1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITRE2,*) FNAMEWFDRE
            WRITE(UNITRE2,*) FNAMEWFDRE
            WRITE(UNITRE2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITRE2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITRE2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITRE2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITRE2,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITIM1,*) FNAMEWFUIM
            WRITE(UNITIM1,*) FNAMEWFUIM
            WRITE(UNITIM1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITIM1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITIM1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITIM1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITIM1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITIM2,*) FNAMEWFDIM
            WRITE(UNITIM2,*) FNAMEWFDIM
            WRITE(UNITIM2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITIM2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITIM2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITIM2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITIM2,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITMO1,*) FNAMEWFUMO
            WRITE(UNITMO1,*) FNAMEWFUMO
            WRITE(UNITMO1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITMO1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITMO1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITMO1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITMO1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITMO2,*) FNAMEWFDMO
            WRITE(UNITMO2,*) FNAMEWFDMO
            WRITE(UNITMO2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITMO2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITMO2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITMO2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITMO2,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITPH1,*) FNAMEWFUPH
            WRITE(UNITPH1,*) FNAMEWFUPH
            WRITE(UNITPH1,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITPH1,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITPH1,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITPH1,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITPH1,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
            WRITE(UNITPH2,*) FNAMEWFDPH
            WRITE(UNITPH2,*) FNAMEWFDPH
            WRITE(UNITPH2,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(UNITPH2,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(UNITPH2,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(UNITPH2,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(UNITPH2,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
          ENDIF
        ENDIF


        CALL WROUT(IDIMEN, .FALSE., .TRUE., IOPTION, NORMAL, COORPO,
     .             DIRVER1, DIRVER2,
     .             NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .             IUNITCD,  NA, NAINCELL, INDICES, XAINCELL )

        WRITE(6,'(A)')
        WRITE(6,'(A)')
     .    '   Generating values of wavefunction on the grid'
        WRITE(6,'(A,I6)')
     .    '   for k-point ',IK
        WRITE(6,'(A,I6)')
     .    '   for eigenstate number ',INDW(IK,IWF)
        WRITE(6,'(A,I6)')
     .    '   Please, wait....'
        WRITE(6,'(A)')


C Allocate space for wave functions in 3D-grid

        IF (.NOT.ALLOCATED(RWF)) THEN
          ALLOCATE(RWF(NPX*NPY*NPZ,NSPIN))
          CALL MEMORY('A','D',NPX*NPY*NPZ*NSPIN,'WAVOFR')
          ALLOCATE(IMWF(NPX*NPY*NPZ,NSPIN))
          CALL MEMORY('A','D',NPX*NPY*NPZ*NSPIN,'WAVOFR')
          ALLOCATE(MWF(NPX*NPY*NPZ,NSPIN))
          CALL MEMORY('A','D',NPX*NPY*NPZ*NSPIN,'WAVOFR')
          ALLOCATE(PWF(NPX*NPY*NPZ,NSPIN))
          CALL MEMORY('A','D',NPX*NPY*NPZ*NSPIN,'WAVOFR')
        ENDIF

C Loop over wavefunctions to write out
      
C Loop over all points in real space -----------------------------------
C      DO 100 NPO = 1, NPX*NPY*NPZ
        NPO = 0
        DO 102 NZ = 1,NPZ
        DO 101 NY = 1,NPY
        DO 100 NX = 1,NPX
          NPO = NPO + 1


C Initialize the wave function at each point -----------------------
          RWAVE  = 0.D0
          RWAVEUP   = 0.D0
          RWAVEDN   = 0.D0
          IWAVE  = 0.D0
          IWAVEUP   = 0.D0
          IWAVEDN   = 0.D0
          MWAVE  = 0.D0
          MWAVEUP   = 0.D0
          MWAVEDN   = 0.D0
          PWAVE  = 0.D0
          PWAVEUP   = 0.D0
          PWAVEDN   = 0.D0

C Localize non-zero orbitals at each point in real space ---------------
          DO IX = 1,3
            XPO(IX) = POINRE(NPO,IX)
          ENDDO
     
          IA   = 0
          ISEL = 0
          NNA  = MAXNA

          CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .                 NNA, JNA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
           DO 110 IAT1 = 1, NNA
             IF( R2IJ(IAT1) .GT. RMAX2 ) EXIT

             IAVEC1   = JNA(IAT1)
             IS1      = ISA(IAVEC1)
             XVEC1(1) = -XIJ(1,IAT1)
             XVEC1(2) = -XIJ(2,IAT1)
             XVEC1(3) = -XIJ(3,IAT1)

C We take atom 1 as center of coordinates

C             PHASE = K(IK,1)*XA(1,IAVEC1)+
C     .               K(IK,2)*XA(2,IAVEC1)+
C     .               K(IK,3)*XA(3,IAVEC1)

             PHASE = K(IK,1)*(XPO(1)+XIJ(1,IAT1))+
     .               K(IK,2)*(XPO(2)+XIJ(2,IAT1))+
     .               K(IK,3)*(XPO(3)+XIJ(3,IAT1))

             SI=SIN(PHASE)
             CO=COS(PHASE)

             DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
               IPHI1 = IPHORB(IO)
               IUO   = INDXUO(IO)
               CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

               IF ( NSPIN .EQ. 1 ) THEN
                 RWAVE  = RWAVE  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,1)*CO - IPSI(IUO,IK,IWF,1)*SI)
                 IWAVE  = IWAVE  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,1)*SI + IPSI(IUO,IK,IWF,1)*CO)
               ELSEIF (NSPIN .EQ. 2) THEN 
                 RWAVEUP  = RWAVEUP  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,1)*CO - IPSI(IUO,IK,IWF,1)*SI)
                 IWAVEUP  = IWAVEUP  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,1)*SI + IPSI(IUO,IK,IWF,1)*CO)
                 RWAVEDN  = RWAVEDN  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,2)*CO - IPSI(IUO,IK,IWF,2)*SI)
                 IWAVEDN  = IWAVEDN  + PHIMU * 
     .            (RPSI(IUO,IK,IWF,2)*SI + IPSI(IUO,IK,IWF,2)*CO)
               ENDIF

C calculate module and phase of wavefunction

               IF ( NSPIN .EQ. 1 ) THEN
                 MWAVE = DSQRT(RWAVE**2 + IWAVE**2)
                 IF (DABS(IWAVE) .LT. 1.D-6) IWAVE = DABS(IWAVE)
                 IF (DABS(RWAVE) .LT. 1.D-12) THEN
                   IF (IWAVE .GT. 0.0D0) PWAVE = PI/2.0D0
                   IF (IWAVE .LT. 0.0D0) PWAVE = -PI/2.0D0
                 ELSE
                   PWAVE = DATAN(IWAVE/RWAVE)
                 ENDIF
                 IF (RWAVE .LT. 0.0D0 .AND. IWAVE .GE. 0.0d0)
     .              PWAVE = PWAVE + PI
                 IF (RWAVE .LT. 0.0D0 .AND. IWAVE .LT. 0.0D0)
     .              PWAVE = PWAVE - PI
               ELSEIF (NSPIN .EQ. 2) THEN 
                 MWAVEUP = DSQRT(RWAVEUP**2 + IWAVEUP**2)

                 IF (DABS(IWAVEUP) .LT. 1.D-6) IWAVEUP = DABS(IWAVEUP)
                 IF (DABS(RWAVEUP) .LT. 1.D-12) THEN
                   IF (IWAVEUP .GT. 0.0D0) PWAVEUP = PI/2.0D0
                   IF (IWAVEUP .LT. 0.0D0) PWAVEUP = -PI/2.0D0
                 ELSE
                   PWAVEUP = DATAN(IWAVEUP/RWAVEUP)
                 ENDIF
                 IF (RWAVEUP .LT. 0.0D0 .AND. IWAVEUP .GE. 0.0d0)
     .              PWAVEUP = PWAVEUP + PI
                 IF (RWAVEUP .LT. 0.0D0 .AND. IWAVEUP .LT. 0.0D0)
     .              PWAVEUP = PWAVEUP - PI

                 MWAVEDN = DSQRT(RWAVEDN**2 + IWAVEDN**2)

                 IF (DABS(IWAVEDN) .LT. 1.D-6) IWAVEDN = DABS(IWAVEDN)
                 IF (DABS(RWAVEDN) .LT. 1.D-12) THEN
                   IF (IWAVEDN .GT. 0.0D0) PWAVEDN = PI/2.0D0
                   IF (IWAVEDN .LT. 0.0D0) PWAVEDN = -PI/2.0D0
                 ELSE
                   PWAVEDN = DATAN(IWAVEDN/RWAVEDN)
                 ENDIF
                 IF (RWAVEDN .LT. 0.0D0 .AND. IWAVEDN .GE. 0.0d0)
     .              PWAVEDN = PWAVEDN + PI
                 IF (RWAVEDN .LT. 0.0D0 .AND. IWAVEDN .LT. 0.0D0)
     .              PWAVEDN = PWAVEDN - PI

               ENDIF

 120         ENDDO
 110       ENDDO
 
           IF (IDIMEN .EQ. 2) THEN
             IF ( NSPIN .EQ. 1 ) THEN
               WRITE(UNITRE1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),RWAVE*SQRT(ARMUNI)
               WRITE(UNITIM1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),IWAVE*SQRT(ARMUNI)
               WRITE(UNITMO1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),MWAVE*SQRT(ARMUNI)
               WRITE(UNITPH1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),PWAVE
             ELSEIF ( NSPIN .EQ. 2 ) THEN
               WRITE(UNITRE1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),RWAVEUP*SQRT(ARMUNI)
               WRITE(UNITRE2,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),RWAVEDN*SQRT(ARMUNI)
               WRITE(UNITIM1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),IWAVEUP*SQRT(ARMUNI)
               WRITE(UNITIM2,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),IWAVEDN*SQRT(ARMUNI)
               WRITE(UNITMO1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),MWAVEUP*SQRT(ARMUNI)
               WRITE(UNITMO2,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),MWAVEDN*SQRT(ARMUNI)
               WRITE(UNITPH1,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),PWAVEUP
               WRITE(UNITPH2,'(3F12.5)')
     .              PLAPO(NPO,1),PLAPO(NPO,2),PWAVEDN
             ENDIF
           ELSE IF (IDIMEN .EQ. 3) THEN
             IF (NSPIN .EQ. 1) THEN
               RWF(NPO,1) = RWAVE*SQRT(ARMUNI)
               IMWF(NPO,1) = IWAVE*SQRT(ARMUNI)
               MWF(NPO,1) = MWAVE*SQRT(ARMUNI)
               PWF(NPO,1) = PWAVE
             ELSE IF (NSPIN .EQ.2) THEN
               RWF(NPO,1) = RWAVEUP*SQRT(ARMUNI)
               RWF(NPO,2) = RWAVEDN*SQRT(ARMUNI)
               IMWF(NPO,1) = IWAVEUP*SQRT(ARMUNI)
               IMWF(NPO,2) = IWAVEDN*SQRT(ARMUNI)
               MWF(NPO,1) = MWAVEUP*SQRT(ARMUNI)
               MWF(NPO,2) = MWAVEDN*SQRT(ARMUNI)
               PWF(NPO,1) = PWAVEUP
               PWF(NPO,2) = PWAVEDN
             ENDIF
           ENDIF



           IF (IDIMEN .EQ. 2) THEN
             IF ( MOD(NPO,NPX) .EQ. 0 ) THEN
               WRITE(UNITRE1,*)
               WRITE(UNITIM1,*)
               IF ( NSPIN .EQ. 2 ) THEN
                 WRITE(UNITRE2,*)
                 WRITE(UNITIM2,*)
               ENDIF
             ENDIF
           ENDIF


C End x loop
 100    ENDDO  
C End y and z loops
 101    ENDDO  
 102    ENDDO  

        IF (IDIMEN .EQ. 3) THEN
          IF (NSPIN .EQ. 1) THEN
            DO NX=1,NPX
              DO NY=1,NPY
                WRITE(UNITRE1,'(6e13.5)')
     .            (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITIM1,'(6e13.5)')
     .            (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITMO1,'(6e13.5)')
     .            (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITPH1,'(6e13.5)')
     .            (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
              ENDDO
            ENDDO
          ELSE IF (NSPIN .EQ. 2) THEN
            DO NX=1,NPX
              DO NY=1,NPY
                WRITE(UNITRE1,'(6e13.5)')
     .            (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITRE2,'(6e13.5)')
     .            (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                WRITE(UNITIM1,'(6e13.5)')
     .            (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITIM2,'(6e13.5)')
     .            (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                WRITE(UNITMO1,'(6e13.5)')
     .            (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITMO2,'(6e13.5)')
     .            (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                WRITE(UNITPH1,'(6e13.5)')
     .            (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                WRITE(UNITPH2,'(6e13.5)')
     .            (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
              ENDDO
            ENDDO
          ENDIF

          WRITE(6,'(A)')
          WRITE(6,'(A)')
     .      '   Your output files are:'
          IF (NSPIN .EQ. 1) THEN
            WRITE(6,'(A,A)') '   ',FNAMEWFRE
            WRITE(6,'(A,A)') '   ',FNAMEWFIM
            WRITE(6,'(A,A)') '   ',FNAMEWFMO
            WRITE(6,'(A,A)') '   ',FNAMEWFPH
          ELSE IF (NSPIN .EQ. 2) THEN
            WRITE(6,'(A,A)') '   ',FNAMEWFURE
            WRITE(6,'(A,A)') '   ',FNAMEWFUIM
            WRITE(6,'(A,A)') '   ',FNAMEWFDRE
            WRITE(6,'(A,A)') '   ',FNAMEWFDIM
            WRITE(6,'(A,A)') '   ',FNAMEWFUMO
            WRITE(6,'(A,A)') '   ',FNAMEWFUPH
            WRITE(6,'(A,A)') '   ',FNAMEWFDMO
            WRITE(6,'(A,A)') '   ',FNAMEWFDPH
          ENDIF
        ENDIF



        CALL IO_CLOSE(UNITRE1)
        CALL IO_CLOSE(UNITIM1)
        IF (NSPIN .EQ. 2) CALL IO_CLOSE(UNITRE2)
        IF (NSPIN .EQ. 2) CALL IO_CLOSE(UNITIM2)
        CALL IO_CLOSE(UNITMO1)
        CALL IO_CLOSE(UNITPH1)
        IF (NSPIN .EQ. 2) CALL IO_CLOSE(UNITMO2)
        IF (NSPIN .EQ. 2) CALL IO_CLOSE(UNITPH2)
     
          
      ENDDO
      ENDDO

      CALL MEMORY('D','D',SIZE(RWF),'wavofr')
      DEALLOCATE(RWF)
      CALL MEMORY('D','D',SIZE(IMWF),'wavofr')
      DEALLOCATE(IMWF)
      CALL MEMORY('D','D',SIZE(MWF),'wavofr')
      DEALLOCATE(MWF)
      CALL MEMORY('D','D',SIZE(PWF),'wavofr')
      DEALLOCATE(PWF)

      RETURN    
      END
