
      SUBROUTINE STM   ( NA, NO, NUO, MAXNA, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .                   RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .                   ZREF, ZMIN, ZMAX, NPX, NPY, NPZ, NSCX, NSCY,
     .                   V0, EMAX, EMIN,
     .                   ARMUNI, IUNITCD, RMAXO )

C **********************************************************************
C Simulate STM images in the Tersoff-Hamann approximation, by
C extrapolating the wavefunctions into vacuum
C
C Coded by P. Ordejon and N. Lorente,  November 2004
C **********************************************************************

      USE ATMFUNCS
      USE FDF
      USE CHEMICAL


      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NUO, NPX, NPY, NPZ, IUNITCD,
     .  NSPIN, MAXNA, NK, NWF(NK), NUMWF, 
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  INDW(NK,NUMWF), NSCX, NSCY

      DOUBLE PRECISION, INTENT(IN) ::
     .  XA(3,NA), ZMIN, ZMAX, ZREF,
     .  ARMUNI, RMAXO, V0, EMAX, EMIN

      DOUBLE PRECISION, INTENT(IN) ::
     . CELL(3,3), 
     . RPSI(NUO,NK,NUMWF,NSPIN), IPSI(NUO,NK,NUMWF,NSPIN),
     . E(NK,NUMWF,NSPIN), K(NK,3)

      DOUBLE PRECISION ::
     . UCELL(3,3), VOLCEL

      EXTERNAL ::
     . VOLCEL
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
C REAL*8  UCELL(3,3)       : Unit cell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C REAL*8 RPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (real part)
C REAL*8 IPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (imag part)
C REAL*8 E(NK,NUMWF,NSPIN) : Eigen energies in eV
C INTEGER INDW(NUMWF,NK)   : Index of the wavefunctions
C INTEGER NWF(NK)          : Number of wavefncts to print for each k-point
C INTEGER NUMWF            : Max num of wavefncts to print a given k-point
C INTEGER NK               : Number of k-points
C REAL*8 K(NK,3)           : k-points
c REAL*8 ZREF              : Position of reference plane for wf. estrapol.
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction for the STM scan
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C INTEGER NSCX, NSCY       : Number of cells in x and y direction to plot
C                            in cube file
C REAL*8  V0               : Value of the potential at the vacuum region in eV
C REAL*8  EMAX             : Maximum value for the energy window for STM in eV
C REAL*8  EMIN             : Minimum value for the energy window for STM in eV
C REAL*8  ARMUNI           : Conversion factor for the charge density
C INTEGER IUNITCD          : Unit of the charge density
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  JNA

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     .   R2IJ

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     .   XIJ

      INTEGER
     .  IA, ISEL, NNA, I, J, IN, IAT1, IO, IUO, IAVEC1, 
     .  IS1, IPHI1, NX, NY, NZ, IWF, IK, ISPIN, UNITRE1,
     .  IX, IY, IZ, NSX, NSY, NAU

      DOUBLE PRECISION
     .  DOT, RMAX, XPO(3), RMAX2, XVEC1(3),
     .  PHIMU, GRPHIMU(3),
     .  PHASE, SI, CO, ENER, PMIKR, SIMIKR, COMIKR, USAVE, VC, VU

      DOUBLE PRECISION, ALLOCATABLE :: RHO(:,:,:)

      DOUBLE COMPLEX
     .  CWAVE, EXPPHI, EXMIKR

      DOUBLE COMPLEX, ALLOCATABLE :: CW(:,:), CWE(:,:,:)
 
      LOGICAL FIRST

      CHARACTER
     .   SNAME*40, FNAME*60, PASTE*60

      EXTERNAL
     .  NEIGHB, IO_ASSIGN, IO_CLOSE, PASTE

C      CHARACTER
C     .  SNAME*40, FNAMEWFRE*60, FNAMEWFIM*60, 
C     .  FNAMEWFURE*60, FNAMEWFUIM*60, FNAMEWFDRE*60, FNAMEWFDIM*60, 
C     .  FNAMEWFMO*60, FNAMEWFPH*60,
C     .  FNAMEWFUMO*60, FNAMEWFUPH*60, FNAMEWFDMO*60, FNAMEWFDPH*60,
C     .  PASTE*60, CHAR1*10, CHAR2*10, ITOCHAR*10, 
C     .  EXT*20, EXT2*25
C
C      EXTERNAL
C     .  IO_ASSIGN, IO_CLOSE, PASTE, PLANE,
C     .  NEIGHB, WROUT, ITOCHAR

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


C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        CALL MEMORY('D','I',SIZE(JNA),'stm')
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        CALL MEMORY('D','D',SIZE(R2IJ),'stm')
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        CALL MEMORY('D','D',SIZE(XIJ),'stm')
        DEALLOCATE(XIJ)
      ENDIF

      ALLOCATE(JNA(MAXNA))
      CALL MEMORY('A','I',MAXNA,'stm')
      ALLOCATE(R2IJ(MAXNA))
      CALL MEMORY('A','D',MAXNA,'stm')
      ALLOCATE(XIJ(3,MAXNA))
      CALL MEMORY('A','D',3*MAXNA,'stm')

      ALLOCATE(CW(0:NPX-1,0:NPY-1))
      CALL MEMORY('A','Z',NPX*NPY,'stm')
      ALLOCATE(CWE(0:NPX-1,0:NPY-1,0:NPZ-1))
      CALL MEMORY('A','Z',NPX*NPY*NPZ,'stm')
      ALLOCATE(RHO(0:NPX-1,0:NPY-1,0:NPZ-1))
      CALL MEMORY('A','D',NPX*NPY*NPZ,'stm')

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

      IF (NSPIN .GT. 2)  STOP 'stm: WRONG NSPIN'

C Check that cell is orthorombic

      IF (UCELL(3,1) /= 0.0D0 .OR. UCELL(3,2) /= 0.0D0 .OR.
     .    UCELL(1,3) /= 0.0D0 .OR. UCELL(2,3) /= 0.0D0) THEN
        WRITE(6,*) 'error: the code only accepts monoclinic cells'
        WRITE(6,*) '       with Z as the vertical axis'
        STOP
      ENDIF


C Initialize density

      RHO = 0.0D0

C Loop over k-points and wavefunctions to include in the STM image

      DO IK  = 1, NK
      WRITE(6,*) 'stm:  Processing kpoin ',IK
      WRITE(6,*) '     --------------------------------'
      DO IWF = 1,NWF(IK)

        WRITE(6,*) 'stm:     Processing w.f. ',iwf

C Check that we have a bound state (E below vacuum level)
        DO ISPIN = 1,NSPIN

          ENER = E(IK,IWF,ISPIN)
          IF (ENER .LT. EMIN .OR. ENER .GT. EMAX) GOTO 99

          IF (E(IK,IWF,ISPIN) .GT. V0) THEN
            WRITE(6,*) 'ERROR: ENERGY EIGENVALUE ',IWF,
     .      ' FOR K-POINT ', IK, 'FOR SPIN ',ISPIN
            WRITE(6,*) '       IS ABOVE VACUUM LEVEL'
           STOP
          ENDIF

C Calculate the value of the w.f. at each point of the reference plane

C Loop over all points in real space -----------------------------------

          DO 101 NY = 1,NPY
          DO 100 NX = 1,NPX


C Initialize the wave function at each point -----------------------
            CWAVE   = (0.0D0, 0.0D0)

C Determine position of current point in the reference plane
            XPO(1) = (NX-1)*UCELL(1,1)/NPX + (NY-1)*UCELL(1,2)/NPY 
            XPO(2) = (NX-1)*UCELL(2,1)/NPX + (NY-1)*UCELL(2,2)/NPY 
            XPO(3) = ZREF

C Phase to cancel the phase of the wave function: -i.k.r
            PMIKR = -(K(IK,1)*XPO(1) + K(IK,2)*XPO(2) + K(IK,3)*XPO(3))
            SIMIKR=DSIN(PMIKR)
            COMIKR=DCOS(PMIKR)
            EXMIKR=DCMPLX(COMIKR,SIMIKR)

C Localize non-zero orbitals at each point in real space ---------------
     
            IA   = 0
            ISEL = 0
            NNA  = MAXNA

            CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .                   NNA, JNA, XIJ, R2IJ, FIRST )

C Loop over Non-zero orbitals ------------------------------------------ 
            DO 110 IAT1 = 1, NNA
              IF( R2IJ(IAT1) .GT. RMAX2 ) GOTO 110

              IAVEC1   = JNA(IAT1)
              IS1      = ISA(IAVEC1)
              XVEC1(1) = -XIJ(1,IAT1)
              XVEC1(2) = -XIJ(2,IAT1)
              XVEC1(3) = -XIJ(3,IAT1)

C XPO + XIJ(IAT1) is just the absolute position of atom IAT1

              PHASE = K(IK,1)*(XPO(1)+XIJ(1,IAT1))+
     .                K(IK,2)*(XPO(2)+XIJ(2,IAT1))+
     .                K(IK,3)*(XPO(3)+XIJ(3,IAT1))

              SI=DSIN(PHASE)
              CO=DCOS(PHASE)
              EXPPHI=DCMPLX(CO,SI)

              DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
                IPHI1 = IPHORB(IO)
                IUO   = INDXUO(IO)
                CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

                CWAVE  = CWAVE  + PHIMU * 
     .          DCMPLX(RPSI(IUO,IK,IWF,ISPIN),IPSI(IUO,IK,IWF,ISPIN)) *
     .          EXPPHI * EXMIKR

 120          ENDDO
 110        ENDDO

            CW(NX-1,NY-1)  = CWAVE * SQRT(ARMUNI)

C End x loop
100       ENDDO  
C End y loop
101       ENDDO  

C Call routine to extrapolate wave function and compute STM image

          ENER = E(IK,IWF,ISPIN)
          CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                     CW,ENER,K(IK,1),CWE)

          RHO = RHO + DREAL(CWE*DCONJG(CWE))

C End loop on spin
99      ENDDO
C End kpoint and wavefunctions loops
      ENDDO
      ENDDO

C Normalize for number of k-points

      RHO = RHO/NK
C ...


C Write charge density in cube format
C Check if lattice vectors in xy plane are orthogonal

      DOT = 0.0
      DO IX=1,3
        DOT = DOT+UCELL(IX,1)*UCELL(IX,2)
      ENDDO


      call io_assign(unitre1)
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(SNAME,'.STM.cube')
      IF (DABS(DOT) .GT. 1.0D-2) THEN
        WRITE(6,*)
        WRITE(6,*) 'stm: WARNING: The cell is not orthorombic, so the'
        WRITE(6,*) '     results can not be plotted in cube format'
        WRITE(6,*)
        GOTO 200
      ELSE
        WRITE(6,*)
        WRITE(6,*) 'stm: writing cube format file ',FNAME
        WRITE(6,*)
        WRITE(6,*) '     ',NSCX,' x ',NSCY,' cells in cube plot'
      ENDIF

C Calculate number of atoms in unit cell
      VC = VOLCEL(CELL)
      VU = VOLCEL(UCELL)
      NAU = NA / IDNINT(VC/VU)

      open(unitre1,file=FNAME,form='formatted',status='unknown')
      WRITE(UNITRE1,*) 'STM'
      WRITE(UNITRE1,*) 'STM'
      WRITE(UNITRE1,'(i5,4f12.6)') NAU, 0.0, 0.0, ZMIN
      WRITE(UNITRE1,'(i5,4f12.6)') NPX*NSCX,(UCELL(1,J)/(NPX-1),J=1,3)
      WRITE(UNITRE1,'(i5,4f12.6)') NPY*NSCY,(UCELL(2,J)/(NPY-1),J=1,3)
      WRITE(UNITRE1,'(i5,4f12.6)') NPZ,0.0,0.0,((ZMAX-ZMIN)/(NPZ-1))

      DO IA = 1, NAU
        WRITE(UNITRE1,'(i5,4f12.6)') ATOMIC_NUMBER(ISA(IA)),0.0, 
     .                               (XA(IX,IA),IX=1,3)
      ENDDO

      DO NSX = 1,NSCX
      DO NX=0,NPX-1
        DO NSY = 1,NSCY
        DO NY=0,NPY-1
          WRITE(UNITRE1,'(6e13.5)')
     .            (RHO(NX,NY,NZ),NZ=0,NPZ-1)
        ENDDO
        ENDDO
      ENDDO
      ENDDO


200   CONTINUE

      call io_close(unitre1)

C Write charge density in Siesta format
C

      call io_assign(unitre1)
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = PASTE(SNAME,'.STM.siesta')
      WRITE(6,*)
      WRITE(6,*) 'stm: writing SIESTA format file', FNAME
      WRITE(6,*)
      open(unitre1,file=FNAME,form='unformatted',
     .         status='unknown')
C write the correct range of the z axis: temporarily override ucell
      USAVE = UCELL(3,3)
      UCELL(3,3) = abs(ZMAX-ZMIN)
      WRITE(unitre1) UCELL
C restore ucell
      UCELL(3,3)=USAVE
      WRITE(unitre1) NPX, NPY, NPZ, 1

      DO IZ=0,NPZ-1
        DO IY=0,NPY-1
C         WRITE(unitre1) (RHO(IND+IX),IX=1,NPX)
          WRITE(unitre1) (REAL(RHO(IX,IY,NPZ-1-IZ)),IX=0,NPX-1)
        ENDDO
      ENDDO

      call io_close(unitre1)



C CLOSE ALLOCATABLE ARRAYS


      DEALLOCATE(RHO)
      CALL MEMORY('D','D',NPX*NPY*NPZ,'stm')
      DEALLOCATE(CWE)
      CALL MEMORY('D','Z',NPX*NPY*NPZ,'stm')
      DEALLOCATE(CW)
      CALL MEMORY('D','Z',NPX*NPY,'stm')

      RETURN    
      END
