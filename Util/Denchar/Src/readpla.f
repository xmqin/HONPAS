
      SUBROUTINE READPLA(MAXA, XA, VOLUME, IDIMEN,
     .                  IOPTION, IUNITCD, ISCALE, NPX, NPY, NPZ,
     .                  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .                  COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                  ARMUNI ) 

C **********************************************************************
C Read the data file to prepare the plane or 3D grid 
C in which we are going to  calculate the charge density
C or wavefunctions
C
C Coded by J. Junquera, November 98
C Modified by P. Ordejon to include 3D and wavefunction capabilities
C **********************************************************************

      use precision
      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  MAXA, IDIMEN
      
      real(dp), INTENT(IN) ::
     .  XA(3,MAXA), VOLUME

      INTEGER, INTENT(OUT) ::
     .  IOPTION, NPX, NPY, NPZ, ISCALE, IUNITCD
     
      real(dp), INTENT(OUT) ::
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), 
     .  ARMUNI

C ****** INPUT *********************************************************
C INTEGER MAXA           : Maximum number of atoms
C REAL*8  XA(3,MAXA)     : Atomic coordinates
C REAL*8 VOLUME          : Volumen of unit cell (in bohr**3)
C INTEGER IDIMEN         : Specify if the run is to plot quantities
C                          in a plane or in a 3D grid (2 or 3, respect)
C ****** OUTPUT ********************************************************
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane
C                          3 = Three points of the plane
C                          4 = Three atomic indices
C INTEGER IUNITCD        : Units for the electron density
C                          IUNITCD = 1 => Ele/(bohr)**3
C                          IUNITCD = 2 => Ele/(Ang)**3
C                          IUNITCD = 3 => Ele/(unitcell)
C INTEGER ISCALE         : Units for the atomic positions
C                          (ISCALE = 1 => Bohrs, ISCALE = 2 => Ang)
C INTEGER NPX, NPY, NPZ  : Number of points generated along x and y
C                          (and z, for 3D-grids) directions ina a system 
C                          of reference in which the third component of 
C                          the points of the plane is zero (Plane 
C                          Reference Frame; PRF)
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX     : Limits of the z-direction in the PRF (for 3D-grids)
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define 
C                          the plane. COORPO(POINT,IX)
C REAL*8  NORMAL(3)      : Components of the normal vector used to define
C                          the plane 
C REAL*8  DIRVER1(3)     : Components of the first vector 
C                          contained in the xy plane 
C REAL*8  DIRVER2(3)     : Components of the first vector 
C                          contained in the xy plane 
C REAL*8  ARMUNI         : Conversion factors for the charge density
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER 
     .  OGP*22, OGP_DEFECT*22,
     .  CPF*22, CPF_DEFECT*22,
     .  UCD*22, UCD_DEFECT*22

      INTEGER
     .  IUNIT, IX, JX, NPX_DEFECT, NPY_DEFECT, NPZ_DEFECT,
     .  IND1, IND2, IND3

      real(dp)
     .  ORIGIN(3), XDIR(3)

      LOGICAL 
     .  LEQI, COLIN

      EXTERNAL 
     .  LEQI, COLINEAR

      DATA ORIGIN /0.D0,0.D0,0.D0/
      DATA XDIR   /1.D0,0.D0,0.D0/
      DATA IND1   /1/
      DATA IND2   /2/
      DATA IND3   /3/
      DATA COLIN  /.FALSE./


      CPF_DEFECT = 'Bohr'

      CPF = FDF_STRING('Denchar.CoorUnits',CPF_DEFECT)

      IF (LEQI(CPF,'bohr')) then
        ISCALE = 1
      ELSEIF (LEQI(CPF,'ang')) then
        ISCALE = 2
      ELSE
        WRITE(6,*)'readpla: ERROR Denchar.CoorUnits must be Ang or Bohr'
        STOP
      ENDIF

      UCD_DEFECT = 'Ele/bohr**3'
      UCD = FDF_STRING('Denchar.DensityUnits',UCD_DEFECT)
      IF (LEQI(UCD,'ele/bohr**3')) then
        IUNITCD = 1
      ELSEIF (LEQI(UCD,'ele/ang**3')) then
        IUNITCD = 2
      ELSEIF (LEQI(UCD,'ele/unitcell')) then
        IUNITCD = 3
      ELSE
       WRITE(6,'(A)')' readpla: ERROR   Wrong Option in Units of      '
       WRITE(6,'(A)')' readpla:  Charge Density                       '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:' 
       WRITE(6,'(A)')' readpla:                                       '
       WRITE(6,'(A)')' readpla:      - Ele/bohr**3                    '
       WRITE(6,'(A)')' readpla:      - Ele/ang**3                     '
       WRITE(6,'(A)')' readpla:      - Ele/unitcell                   '
       STOP
      ENDIF

      NPX_DEFECT = 50
      NPY_DEFECT = 50
      NPZ_DEFECT = 50
      NPX = FDF_INTEGER('Denchar.NumberPointsX',NPX_DEFECT)
      NPY = FDF_INTEGER('Denchar.NumberPointsY',NPY_DEFECT)
      IF (IDIMEN .EQ. 2) THEN
        NPZ = 1
      ELSE IF (IDIMEN .EQ. 3) THEN
        NPZ = FDF_INTEGER('Denchar.NumberPointsZ',NPZ_DEFECT)
      ENDIF

      XMIN = FDF_PHYSICAL('Denchar.MinX',-3.D0,'Bohr')
      XMAX = FDF_PHYSICAL('Denchar.MaxX', 3.D0,'Bohr')
      YMIN = FDF_PHYSICAL('Denchar.MinY',-3.D0,'Bohr')
      YMAX = FDF_PHYSICAL('Denchar.MaxY', 3.D0,'Bohr')
      IF (IDIMEN .EQ. 2) THEN
        ZMIN = 0.D0
        ZMAX = 0.D0
      ELSE IF (IDIMEN .EQ. 3) THEN
        ZMIN = FDF_PHYSICAL('Denchar.MinZ',-3.D0,'Bohr')
        ZMAX = FDF_PHYSICAL('Denchar.MaxZ', 3.D0,'Bohr')
      ENDIF

      OGP_DEFECT = 'NormalVector'
      OGP = FDF_STRING('Denchar.PlaneGeneration',OGP_DEFECT)
      IF (LEQI(OGP,'normalvector')) then
        IOPTION = 1
      ELSEIF (LEQI(OGP,'twolines')) then
        IOPTION = 2
      ELSEIF (LEQI(OGP,'threepoints')) then
        IOPTION = 3
      ELSEIF (LEQI(OGP,'threeatomicindices')) then
        IOPTION = 4
      ELSE
       WRITE(6,'(A)')' readpla: ERROR Wrong Option to Generate Plane  '
       WRITE(6,'(A)')' readpla:  You must choose one of the following:'
       WRITE(6,'(A)')' readpla:                                       '
       WRITE(6,'(A)')' readpla:      - NormalVector                   '
       WRITE(6,'(A)')' readpla:      - TwoLines                       '
       WRITE(6,'(A)')' readpla:      - ThreePoints                    '
       WRITE(6,'(A)')' readpla:      - ThreeAtomicIndices             '
       STOP
      ENDIF

      IF ( FDF_BLOCK('Denchar.CompNormalVector',IUNIT) ) THEN
        READ(IUNIT,*)(NORMAL(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Comp2Vectors',IUNIT) ) THEN
        READ(IUNIT,*)(DIRVER1(IX),IX=1,3)
        READ(IUNIT,*)(DIRVER2(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Coor3Points',IUNIT) ) THEN
        READ(IUNIT,*)(COORPO(1,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(2,IX),IX=1,3)
        READ(IUNIT,*)(COORPO(3,IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.Indices3Atoms',IUNIT) ) THEN
        READ(IUNIT,*)IND1, IND2, IND3
      ENDIF



      IF ( IOPTION .EQ. 4 ) THEN
        DO IX = 1,3
          COORPO(1,IX) = XA(IX,IND1)
        ENDDO
        DO IX = 1,3
          COORPO(2,IX) = XA(IX,IND2)
        ENDDO
        DO IX = 1,3
          COORPO(3,IX) = XA(IX,IND3)
        ENDDO
      ENDIF

C Check if the three points are colinear -------------------------------
      IF ((IOPTION .EQ. 3) .OR. (IOPTION .EQ. 4))THEN
         CALL COLINEAR( COORPO, COLIN )
         IF(COLIN) THEN
           WRITE(6,*)'The coordinates of the three points are colinear'
           WRITE(6,*)'and do not define a plane' 
           WRITE(6,*)'Please, check these coordinates in the input file'
           STOP
         ENDIF
      ENDIF
 

      IF ( FDF_BLOCK('Denchar.PlaneOrigin',IUNIT) ) THEN
        READ(IUNIT,*)(ORIGIN(IX),IX=1,3)
      ENDIF

      IF ( FDF_BLOCK('Denchar.X_Axis',IUNIT) ) THEN
        READ(IUNIT,*)(XDIR(IX),IX=1,3)
      ENDIF

      IF (IOPTION .LT. 3) THEN
        DO IX = 1,3      
          COORPO(1,IX) = ORIGIN(IX)
        ENDDO
        IF(IOPTION .EQ. 1) THEN
          DO IX = 1,3
            COORPO(2,IX) = XDIR(IX)
          ENDDO
        ENDIF
      ENDIF

C Scale points coordinates
C   Iscale = 1 => Do nothing
C   Iscale = 2 => Multiply by 1./0.529177 (Ang --> Bohr)

      IF( (ISCALE .EQ. 2) .AND. (IOPTION .NE. 4) ) THEN
        DO IX = 1,3
          DO JX = 1,3
            COORPO(JX,IX) = 1.D0 / 0.529177D0 * COORPO(JX,IX)
          ENDDO
          ORIGIN(IX)  = 1.D0 / 0.529177D0 * ORIGIN(IX)
          DIRVER1(IX) = 1.D0 / 0.529177D0 * DIRVER1(IX)
          DIRVER2(IX) = 1.D0 / 0.529177D0 * DIRVER2(IX)
        ENDDO
      ENDIF 

C Units of Charge Density
C   Iunitcd = 1 => Do nothing
C   Iunitcd = 2 => Multiply by (1.d0 / 0.529177d0) **3 (bohr**3 --> Ang**3)
C   Iunitcd = 3 => Multiply by volume unit cell (in bohrs**3) 

      IF (IUNITCD .EQ. 1) THEN
        ARMUNI = 1.D0
      ELSEIF( IUNITCD .EQ. 2 ) THEN
        ARMUNI = (1.D0 / 0.529177D0)**3 
      ELSEIF( IUNITCD .EQ. 3 ) THEN
        ARMUNI = VOLUME
      ENDIF

      END

