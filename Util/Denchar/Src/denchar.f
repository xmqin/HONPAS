
      PROGRAM DENCHAR

C **********************************************************************
C Reads density matrix from SIESTA and calculates the charge density
C at the points of a plane in real space, or at a 3D grid of points
C Coded by J. Junquera 11/98
C Modified by J. Junquera 07/01 
C 3D and wavefunction capabilities coded by P. Ordejon, June 2003
C Modified to handle complex wavefunctions and multiple k-points
C by P. Ordejon, July 2004
C Modified to include the local density of states at a point
C by P. Ordejon, July 2004
C Modified to use the more efficient WFSX format for wave functions
C by A. Garcia, May 2012
C
C Version: 2.1
C **********************************************************************
C
C  Modules
C
      USE PRECISION
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
      use parallel, only: nodes, node, ionode

      IMPLICIT NONE

      INTEGER
     .   NO_U, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NSC(3), NK, NUMWF, NE
      integer :: ns_dummy

      INTEGER
     .  IDIMEN, IOPTION, NPX, NPY, NPZ, IUNITCD, ISCALE, BFUNC

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR, NWF

      INTEGER, DIMENSION(:,:), ALLOCATABLE ::
     .  INDW


      real(dp)
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO

      real(dp)
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3),
     .  DIRVER2(3), ARMUNI, EMIN, DELTAE, ETA, XSTS(3)

      real(dp), DIMENSION(:,:,:,:), ALLOCATABLE ::
     .   RPSI,IPSI

      real(dp), DIMENSION(:,:,:), ALLOCATABLE ::
     .   E

      real(dp), DIMENSION(:,:), ALLOCATABLE ::
     .   XA, DSCF, K

      real(dp), DIMENSION(:), ALLOCATABLE ::
     .   DATM

      CHARACTER
     .  FILEIN*20, FILEOUT*20

      LOGICAL 
     .  FOUND, CHARGE, WAVES, STS

      EXTERNAL
     .  IODM, READPLA, REDATA_DENCHAR, REINIT, RHOOFR, VOLCEL
      external :: readwavesx

      DATA NORMAL /0.D0,0.D0,1.D0/
      DATA COORPO /1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/
      DATA DIRVER1 /1.D0,0.D0,0.D0/
      DATA DIRVER2 /0.D0,1.D0,0.D0/

C ****** READ FROM SIESTA **********************************************
C INTEGER NO_U                : Total number of orbitals in the unit cell
C INTEGER NO_S                : Total number of orbitals in the supercell
C INTEGER NA_S                : Total number of atoms in the supercell
C INTEGER NSPIN               : Number of different spin polarizations
C                               Nspin = 1 => Unpolarized, Nspin = 2 => Polarized
C INTEGER MAXND               : Maximum number
C                               of basis orbitals interacting, either directly
C                               or through a KB projector, with any orbital
C INTEGER MAXNA               : Maximum number of neighbours of any atom
C INTEGER NSC(3)              : Num. of unit cells in each supercell direction
C INTEGER NUMWF               : Max num of wavefncts to print for a given k-po.
C INTEGER NWF(NK)             : Num of wavefncts to print for each k-point
C INTEGER ISA(MAXA)           : Species index of each atom in the supercell
C INTEGER LASTO(0:MAXA)       : Last orbital of each atom in array iphorb
C INTEGER IPHORB(MAXO)        : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)        : Equivalent orbital in unit cell
C INTEGER NUMD(NO_U)          : Number of nonzero elements of each row of the
C                               Hamiltonian matrix between atomic orbitals
C INTEGER LISTD(MAXND)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTDPTR(NO_U)      : Pointer to where each row of listh starts - 1
C                               The reason for pointing to the element before
C                               the first one is so that when looping over the
C                               elements of a row there is no need to shift by
C                               minus one.
C INTEGER NK                  : Number of k-points in wave functions file
C INTEGER INDW(NK,NUMWF)      : Indices of wavefuncs to print for each k-point
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)
C REAL*8  DATM(NO_S)          : Occupations of basis orbitals in free atom
C REAL*8  DSCF(MAXND,NSPIN)   : Density Matrix (DM)
C REAL*8  RPSI(NO_U,NK,NUMWF,NSPIN): Wave function coefficients (real part)
C REAL*8  IPSI(NO_U,NK,NUMWF,NSPIN): Wave function coefficients (imag part)
C REAL*8  E(NK,NUMWF,NSPIN)        : Wave function energies
C REAL*8  K(NK,3)             : K-points
C ****** INFORMATION OF THE POINT, PLANE OR 3D GRID ***********************
C INTEGER IDIMEN              : Specifies 2D or 3D mode
C LOGICAL CHARGE              : Should charge density be computed?
C LOGICAL WAVES               : Should wave functions be computed?
C LOGICAL STS                 : Should STS simulation be done?
C INTEGER IOPTION             : Option to generate the plane or 3D grid
C                               1 = Normal vector to xy plane (ie, z direction)
C                               2 = Two vectors belonging to the xy plane
C                               3 = Three points of the xy plane
C                               4 = Three atomic indices define the xy plane
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) directions in a system of 
C                               reference in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER NE                  : Number of energies for STS curve
C REAL*8 EMIN                 : Minimum energy of STS curve
C REAL*8 DELTAE               : Step in STS curve
C REAL*8 ETA                  : Lorenzian width used in STS curve
C INTEGER IUNITCD             : Units of the charge density
C INTEGER ISCALE              : Units of the points of the plane or 3D grid
C REAL*8  XMIN, XMAX          : Limits of the plane in the PRF for x-axis
C REAL*8  YMIN, YMAX          : Limits of the plane in the PRF for y-axis
C REAL*8  ZMIN, ZMAX          : Limits of the or 3D grid in the PRF for z-axis
C REAL*8  COORPO(3,3)         : Coordinates of the three points used 
C                               to define the xy plane
C REAL*8  NORMAL(3)           : Components of the normal vector 
C                               used to define the xy plane (z-direction)
C REAL*8  DIRVER1(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8  DIRVER2(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8 XSTS(3)              : Coordinates of points where STS spectrum
C                               is calculated
C INTEGER BFUNC               : Broadening function for STS spectra:
C                               1 = Lorentzian, 2 = Gaussian
C REAL*8  ARMUNI              : Conversion factors for the charge density
C ****** INTERNAL VARIABLES ********************************************
C LOGICAL FOUND               : Has DM been found in disk?
C                               (Only when task = 'read')
C **********************************************************************
C
      nodes = 1
      node = 0
      ionode = .true.

C Set up fdf -----------------------------------------------------------
      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINIT( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN, IDIMEN,
     .            CHARGE, WAVES, STS )

C Allocate some variables ----------------------------------------------
      ALLOCATE(XA(3,NA_S))
      CALL MEMORY('A','D',3*NA_S,'denchar')

      ALLOCATE(LASTO(0:NA_S))
      CALL MEMORY('A','D',NA_S+1,'denchar')

      ALLOCATE(ISA(NA_S))
      CALL MEMORY('A','D',NA_S,'denchar')

      ALLOCATE(IPHORB(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

      ALLOCATE(INDXUO(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

      ALLOCATE(DATM(NO_S))
      CALL MEMORY('A','D',NO_S,'denchar')

C Read some variables from SIESTA --------------------------------------
      CALL REDATA_DENCHAR( NO_S, NA_S, NO_U, MAXND, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO, DATM )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns_dummy)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))


C Allocate variables

C If this is a charge calculation, allocate space for DM
      IF (CHARGE) THEN
        ALLOCATE(LISTDPTR(NO_U))
        CALL MEMORY('A','I',NO_U,'denchar')
        LISTDPTR(:) = 0

        ALLOCATE(NUMD(NO_U))
        CALL MEMORY('A','I',NO_U,'denchar')
        NUMD(:) = 0

C Allocate some other variables ----------------------------------------
        IF (.NOT.ALLOCATED(LISTD)) THEN
          ALLOCATE(LISTD(MAXND))
          CALL MEMORY('A','I',MAXND,'denchar')
        ENDIF
        
        IF (ALLOCATED(DSCF)) THEN
          CALL MEMORY('D','D',SIZE(DSCF),'denchar')
          DEALLOCATE(DSCF)
        ENDIF
        ALLOCATE(DSCF(MAXND,NSPIN))
        CALL MEMORY('A','D',MAXND*NSPIN,'denchar')

C Read Density Matrix from files ---------------------------------------
        CALL IODM('READ', MAXND, NO_U, NSPIN,
     .            NUMD, LISTDPTR, LISTD, DSCF, FOUND )
        IF (.NOT. FOUND) THEN
          WRITE(6,*)' DENSITY MATRIX NOT FOUND              '
          WRITE(6,*)' CHECK YOU HAVE COPIED IT FROM THE       '
          WRITE(6,*)' DIRECTORY WHERE YOU HAVE RUN SIESTA   '
          STOP
        ENDIF 
      ENDIF

      IF (WAVES .OR. STS) THEN

C Call readwaves just to find out number of wavefunctions to print
C or to calculate STS spectra.
C Allocate temporary space for eigenvalues and eigenfunctions
        NUMWF = 1
        NK = 1

        ALLOCATE(NWF(NK))
        CALL MEMORY('A','I',NK,'denchar')
        ALLOCATE(INDW(NK,NUMWF))
        CALL MEMORY('A','I',NUMWF*NK,'denchar')
        ALLOCATE(K(NK,3))
        CALL MEMORY('A','D',NK*3,'denchar')
        ALLOCATE(E(NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NUMWF*NSPIN*NK,'denchar')
        ALLOCATE(RPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'denchar')
        ALLOCATE(IPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'denchar')

!!!!!!!!        CALL READWAVES(NSPIN,NO_U,0,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)
        CALL READWAVESX(NSPIN,NO_U,0,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)

C deallocate temporary space 
        CALL MEMORY('D','I',SIZE(NWF),'denchar')
        DEALLOCATE(NWF)
        CALL MEMORY('D','I',SIZE(INDW),'denchar')
        DEALLOCATE(INDW)
        CALL MEMORY('D','D',SIZE(K),'denchar')
        DEALLOCATE(K)
        CALL MEMORY('D','D',SIZE(E),'denchar')
        DEALLOCATE(E)
        CALL MEMORY('D','D',SIZE(RPSI),'denchar')
        DEALLOCATE(RPSI)
        CALL MEMORY('D','D',SIZE(IPSI),'denchar')
        DEALLOCATE(IPSI)

C allocate space for eigenenergies and (complex) eigenfunctions 
        ALLOCATE(NWF(NK))
        CALL MEMORY('A','I',NK,'denchar')
        ALLOCATE(INDW(NK,NUMWF))
        CALL MEMORY('A','I',NUMWF*NK,'denchar')
        ALLOCATE(K(NK,3))
        CALL MEMORY('A','D',NK*3,'denchar')
        ALLOCATE(E(NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NUMWF*NSPIN*NK,'denchar')
        ALLOCATE(RPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'denchar')
        ALLOCATE(IPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'denchar')

C call readwaves again to actually read wavefunctions
!!!!        CALL READWAVES(NSPIN,NO_U,1,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)
        CALL READWAVESX(NSPIN,NO_U,1,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)
      ENDIF

C Find out parameters of STS simulation
      IF (STS) THEN
        CALL READSTS(VOLUME, EMIN, NE, DELTAE, ETA, XSTS, BFUNC,
     .               IUNITCD, ARMUNI)
      ENDIF

      IF (CHARGE .OR. WAVES) THEN
C Read option to generate the plane or 3D-grid -------------------------
        CALL READPLA( NA_S, XA, VOLUME, IDIMEN,
     .                IOPTION, IUNITCD, ISCALE, NPX, NPY, NPZ,
     .                XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .                COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                ARMUNI )
      ENDIF

C Form Density Matrix for Neutral and Isolated Atoms -------------------
      IF (CHARGE) THEN
      
        CALL RHOOFR( NA_S, NO_S, NO_U, MAXND, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO,
     .               XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DATM,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO )
      ENDIF

      IF (WAVES) THEN
        CALL WAVOFR( NA_S, NO_S, NO_U, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .               RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO )
      ENDIF

      IF (STS) THEN
        CALL STSOFR( NA_S, NO_S, NO_U, MAXNA, NSPIN,
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .               RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .               ARMUNI, IUNITCD, RMAXO,
     .               EMIN, NE, DELTAE, ETA, XSTS, BFUNC )
      ENDIF

      END PROGRAM DENCHAR
