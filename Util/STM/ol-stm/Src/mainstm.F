
      PROGRAM MAINSTM

C **********************************************************************
C Reads eigenfncts and eigenvalues from SIESTA and computes STM images
C
C Coded by P. Ordejon, November 2004
C
C Modified by N. Lorente, August 2005
C
C **********************************************************************
C
C  Modules
C
      USE PRECISION
      USE PARALLEL
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF

      IMPLICIT NONE

      INTEGER ns,
     .   NO_U, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NSC(3), NK, NUMWF

      INTEGER
     .  NPX, NPY, NPZ, IUNITCD, IX, IV, NSCX, NSCY

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR, NWF

      INTEGER, DIMENSION(:,:), ALLOCATABLE ::
     .   INDW

      DOUBLE PRECISION
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO, UCELL(3,3), V0, EMAX, EMIN

      DOUBLE PRECISION
     .  ZMIN, ZMAX, ZREF, ZREF2,
     .  ARMUNI

      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE ::
     .   RPSI,IPSI

      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::
     .   E

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .   XA, DSCF, K

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .   DATM, DSCFNA

      CHARACTER
     .  FILEIN*20, FILEOUT*20

      LOGICAL 
     .  FOUND

      EXTERNAL
     .  IODM, REDATA, REINITSTM, VOLCEL


C ****** READ FROM SIESTA **********************************************
! integer ns                  : Total number of species
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
C REAL*8  UCELL(3,3)          : Unit cell vectors UCELL(IXYZ,IVECT)
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
C ****** INFORMATION OF THE PLANE OR 3D GRID ***************************
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) directions in a system of 
C                               reference in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER IUNITCD             : Units of the charge density
C REAL*8  ZREF                : Position of reference plane for wf. extrapol.
C REAL*8  ZREF2                : Position of reference plane for wf. extrapol.
C REAL*8  ZMIN, ZMAX          : Limits of the or 3D grid in the PRF for z-axis
C INTEGER NSCX, NSCY          : Number of cells to plot in cube file (in x
C                               and y directions)
C REAL*8  ARMUNI              : Conversion factors for the charge density
C ****** INTERNAL VARIABLES ********************************************
C REAL*8  DSCFNA(MAXND)       : Density Matrix for Neutral Atoms
C LOGICAL FOUND               : Has DM been found in disk?
C                               (Only when task = 'read')
C **********************************************************************


C Set up fdf -----------------------------------------------------------
      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

      WRITE(6,*)
      WRITE(6,*) 'STM Simulation program '
      WRITE(6,*) 'P. Ordejon and N. Lorente, Nov. 04'
      WRITE(6,*)
      WRITE(6,*) 'Reading information from SIESTA'
      WRITE(6,*)


C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINITSTM( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN)

C Allocate some variables ----------------------------------------------
      ALLOCATE(XA(3,NA_S))
      CALL MEMORY('A','D',3*NA_S,'mainstm')

      ALLOCATE(LASTO(0:NA_S))
      CALL MEMORY('A','D',NA_S+1,'mainstm')

      ALLOCATE(ISA(NA_S))
      CALL MEMORY('A','D',NA_S,'mainstm')

      ALLOCATE(IPHORB(NO_S))
      CALL MEMORY('A','D',NO_S,'mainstm')

      ALLOCATE(INDXUO(NO_S))
      CALL MEMORY('A','D',NO_S,'mainstm')

      ALLOCATE(DATM(NO_S))
      CALL MEMORY('A','D',NO_S,'mainstm')

C Read some variables from SIESTA --------------------------------------
      CALL REDATA( NO_S, NA_S, NO_U, MAXND, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO, DATM )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))

C Calculate unit cell vectors
      DO IX=1,3
        DO IV=1,3
          UCELL(IX,IV) = CELL(IX,IV)/NSC(IV)
        ENDDO
      ENDDO


C Allocate variables



C call readwaves just to find out number of wavefunctions to print
C allocate temporary space for eigenvalues and eigenfunctions
        NUMWF = 1
        NK = 1

        ALLOCATE(NWF(NK))
        CALL MEMORY('A','I',NK,'mainstm')
        ALLOCATE(INDW(NK,NUMWF))
        CALL MEMORY('A','I',NUMWF*NK,'mainstm')
        ALLOCATE(K(NK,3))
        CALL MEMORY('A','D',NK*3,'mainstm')
        ALLOCATE(E(NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NUMWF*NSPIN*NK,'mainstm')
        ALLOCATE(RPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'mainstm')
        ALLOCATE(IPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'mainstm')

        CALL READWAVES(NSPIN,NO_U,0,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)

C deallocate temporary space 
        CALL MEMORY('D','I',SIZE(NWF),'mainstm')
        DEALLOCATE(NWF)
        CALL MEMORY('D','I',SIZE(INDW),'mainstm')
        DEALLOCATE(INDW)
        CALL MEMORY('D','D',SIZE(K),'mainstm')
        DEALLOCATE(K)
        CALL MEMORY('D','D',SIZE(E),'mainstm')
        DEALLOCATE(E)
        CALL MEMORY('D','D',SIZE(RPSI),'mainstm')
        DEALLOCATE(RPSI)
        CALL MEMORY('D','D',SIZE(IPSI),'mainstm')
        DEALLOCATE(IPSI)

C allocate space for eigenenergies and (complex) eigenfunctions 
        ALLOCATE(NWF(NK))
        CALL MEMORY('A','I',NK,'mainstm')
        ALLOCATE(INDW(NK,NUMWF))
        CALL MEMORY('A','I',NUMWF*NK,'mainstm')
        ALLOCATE(K(NK,3))
        CALL MEMORY('A','D',NK*3,'mainstm')
        ALLOCATE(E(NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NUMWF*NSPIN*NK,'mainstm')
        ALLOCATE(RPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'mainstm')
        ALLOCATE(IPSI(NO_U,NK,NUMWF,NSPIN))
        CALL MEMORY('A','D',NO_U*NUMWF*NSPIN*NK,'mainstm')

C call readwaves again to actually read wavefunctions
        CALL READWAVES(NSPIN,NO_U,1,NWF,NUMWF,NK,RPSI,IPSI,E,K,INDW)

C Read option for STM 
      CALL READSTM( VOLUME, 
     .     IUNITCD, NPX, NPY, NPZ, ZREF, ZREF2, ZMIN, ZMAX, EMAX, EMIN,
     .      NSCX, NSCY, ARMUNI )

C Call routine to calculate the value of the potential at vacuum

      CALL VACPOT(V0)

C Call routine to compute STM image

        CALL STM   ( NA_S, NO_S, NO_U, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .               RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .               ZREF,ZREF2, ZMIN, ZMAX, NPX, NPY, NPZ, NSCX, NSCY,
     .               V0, EMAX, EMIN,
     .               ARMUNI, IUNITCD, RMAXO )


      END PROGRAM MAINSTM
