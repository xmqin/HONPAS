! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE STM   ( NA, NO, NUO, MAXNA, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .                   RPSI, IPSI, E, INDW, NWF, NUMWF, NK, K,
     .                   ZREF, ZMIN, ZMAX, NPX, NPY, NPZ, V0,
     .                   ARMUNI, IUNITCD, RMAXO )

C **********************************************************************
C Simulate STM images in the Tersoff-Hamann approximation, by
C extrapolating the wavefunctions into vacuum
C Coded by P. Ordejon and N. Lorente,  November 2004
C **********************************************************************

      USE ATMFUNCS


      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NUO, NPX, NPY, NPZ, IUNITCD,
     .  NSPIN, MAXNA, NK, NWF(NK), NUMWF, 
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA),
     .  INDW(NK,NUMWF)

      DOUBLE PRECISION, INTENT(IN) ::
     .  XA(3,NA), ZMIN, ZMAX, ZREF,
     .  ARMUNI, RMAXO, V0

      DOUBLE PRECISION, INTENT(IN) ::
     . CELL(3,3), UCELL(3,3),
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
C REAL*8  UCELL(3,3)       : Unit cell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C REAL*8 RPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (real part)
C REAL*8 IPSI(NUO,NK,NUMWF,NSPIN): Wave function coefficients (imag part)
C REAL*8 E(NK,NUMWF,NSPIN) : Eigen energies
C INTEGER INDW(NUMWF,NK)   : Index of the wavefunctions
C INTEGER NWF(NK)          : Number of wavefncts to print for each k-point
C INTEGER NUMWF            : Max num of wavefncts to print a given k-point
C INTEGER NK               : Number of k-points
C REAL*8 K(NK,3)           : k-points
c REAL*8 ZREF              : Position of reference plane for wf. estrapol.
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction for the STM scan
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C REAL*8  V0               : Value of the potential at the vacuum region
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
     .  IS1, IPHI1, NX, NY, NZ, IWF, IK, ISPIN

      DOUBLE PRECISION
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  PHIMU, GRPHIMU(3),
     .  PHASE, SI, CO, ENER, PMIKR, SIMIKR, COMIKR

      DOUBLE COMPLEX
     .  CWAVE, CWAVEUP, CWAVEDN, EXPPHI, EXMIKR

      DOUBLE COMPLEX, ALLOCATABLE :: CW(:,:), CWUP(:,:), CWDN(:,:),
     .                      CWE(:,:,:), CWEUP(:,:,:), CWEDN(:,:,:)
 
      LOGICAL FIRST

      EXTERNAL
     .  NEIGHB

C      CHARACTER
C     .  SNAME*40, FNAMEWFRE*60, FNAMEWFIM*60, 
C     .  FNAMEWFURE*60, FNAMEWFUIM*60, FNAMEWFDRE*60, FNAMEWFDIM*60, 
C     .  FNAMEWFMO*60, FNAMEWFPH*60,
C     .  FNAMEWFUMO*60, FNAMEWFUPH*60, FNAMEWFDMO*60, FNAMEWFDPH*60,
C     .  CHAR1*10, CHAR2*10, ITOCHAR*10, 
C     .  EXT*20, EXT2*25
C
C      EXTERNAL
C     .  IO_ASSIGN, IO_CLOSE, PLANE,
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

      IF (NSPIN .EQ. 1) THEN
        ALLOCATE(CW(0:NPX-1,0:NPY-1))
        CALL MEMORY('A','Z',NPX*NPY,'stm')
        ALLOCATE(CWE(0:NPX-1,0:NPY-1,0:NPZ-1))
        CALL MEMORY('A','Z',NPX*NPY*NPZ,'stm')
      ELSE IF (NSPIN .EQ. 2) THEN
        ALLOCATE(CWUP(0:NPX-1,0:NPY-1))
        CALL MEMORY('A','Z',NPX*NPY,'stm')
        ALLOCATE(CWEUP(0:NPX-1,0:NPY-1,0:NPZ-1))
        CALL MEMORY('A','Z',NPX*NPY*NPZ,'stm')
        ALLOCATE(CWDN(0:NPX-1,0:NPY-1))
        CALL MEMORY('A','Z',NPX*NPY,'stm')
        ALLOCATE(CWEDN(0:NPX-1,0:NPY-1,0:NPZ-1))
        CALL MEMORY('A','Z',NPX*NPY*NPZ,'stm')
      ELSE
        STOP 'stm: WRONG NSPIN'
      ENDIF

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

C Loop over k-points and wavefunctions to include in the STM image

      DO IK  = 1, NK
      DO IWF = 1,NWF(IK)

C Check that we have a bound state (E below vacuum level)
C       DO ISPIN = 1,NSPIN
C         IF (E(IK,IWF,ISPIN) .GT. V0) THEN
C          WRITE(6,*) 'ERROR: ENERGY EIGENVALUE ',IWF,
C    .     ' FOR K-POINT ', IK, 'FOR SPIN ',ISPIN
C          WRITE(6,*) '       IS ABOVE VACUUM LEVEL'
C          STOP
C         ENDIF
C       ENDDO


C Calculate the value of the w.f. at each point of the reference plane

C Loop over all points in real space -----------------------------------

        DO 101 NY = 1,NPY
        DO 100 NX = 1,NPX


C Initialize the wave function at each point -----------------------
          CWAVE   = (0.0D0, 0.0D0)
          CWAVEUP = (0.0D0, 0.0D0)
          CWAVEDN = (0.0D0, 0.0D0)

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
     .                 NNA, JNA, XIJ, R2IJ, FIRST )

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
     .              K(IK,2)*(XPO(2)+XIJ(2,IAT1))+
     .              K(IK,3)*(XPO(3)+XIJ(3,IAT1))

            SI=DSIN(PHASE)
            CO=DCOS(PHASE)
            EXPPHI=DCMPLX(CO,SI)

            DO 120 IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
              IPHI1 = IPHORB(IO)
              IUO   = INDXUO(IO)
              CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )

              IF ( NSPIN .EQ. 1 ) THEN
                CWAVE  = CWAVE  + PHIMU * 
     .          DCMPLX(RPSI(IUO,IK,IWF,1),IPSI(IUO,IK,IWF,1)) *
     .          EXPPHI * EXMIKR
              ELSEIF (NSPIN .EQ. 2) THEN 
                CWAVEUP  = CWAVEUP  + PHIMU * 
     .          DCMPLX(RPSI(IUO,IK,IWF,1),IPSI(IUO,IK,IWF,1))*
     .          EXPPHI * EXMIKR
                CWAVEDN  = CWAVEDN  + PHIMU * 
     .          DCMPLX(RPSI(IUO,IK,IWF,2),IPSI(IUO,IK,IWF,2))*
     .          EXPPHI * EXMIKR
              ENDIF

 120        ENDDO
 110      ENDDO

          IF ( NSPIN .EQ. 1 ) THEN
            CW(NX-1,NY-1)  = CWAVE * SQRT(ARMUNI)
          ELSEIF (NSPIN .EQ. 2) THEN 
            CWUP(NX-1,NY-1)  = CWAVEUP * SQRT(ARMUNI)
            CWDN(NX-1,NY-1)  = CWAVEDN * SQRT(ARMUNI)
          ENDIF

C End x loop
 100    ENDDO  
C End y loop
 101    ENDDO  

C Call routine to extrapolate wave function and compute STM image

        IF ( NSPIN .EQ. 1 ) THEN
          ENER = E(IK,IWF,1)
          CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                     CW,ENER,K(IK,1),CWE)
        ELSEIF (NSPIN .EQ. 2) THEN 
          ENER = E(IK,IWF,1)
          CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                     CWUP,ENER,K(IK,1),CWEUP)
          ENER = E(IK,IWF,2)
          CALL EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                     CWDN,ENER,K(IK,1),CWEDN)
        ENDIF
          
C End kpoint and wavefunctions loops
      ENDDO
      ENDDO

      IF (NSPIN .EQ. 1) THEN
        DEALLOCATE(CWE)
        CALL MEMORY('D','Z',NPX*NPY*NPZ,'stm')
        DEALLOCATE(CW)
        CALL MEMORY('D','Z',NPX*NPY,'stm')
      ELSE IF (NSPIN .EQ. 2) THEN
        DEALLOCATE(CWEDN)
        DEALLOCATE(CWDN)
        DEALLOCATE(CWEUP)
        DEALLOCATE(CWUP)
        CALL MEMORY('D','Z',NPX*NPY*NPZ,'stm')
        CALL MEMORY('D','Z',NPX*NPY,'stm')
        CALL MEMORY('D','Z',NPX*NPY*NPZ,'stm')
        CALL MEMORY('D','Z',NPX*NPY,'stm')
      ENDIF


      RETURN    
      END
