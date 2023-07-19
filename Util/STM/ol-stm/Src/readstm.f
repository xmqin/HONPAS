
      SUBROUTINE READSTM(VOLUME, 
     .      IUNITCD, NPX, NPY, NPZ, ZREF, ZREF2, ZMIN, ZMAX, EMAX, EMIN,
     .      NSCX, NSCY, ARMUNI ) 

C **********************************************************************
C Read the data file with the input for the STM calculation
C
C Coded by P. Ordejon, November 2004
C Based on readpla. (by Junquera and Ordejon)
C
C Modified by N. Lorente, August 2005
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) ::
     .  VOLUME

      INTEGER, INTENT(OUT) ::
     .  NPX, NPY, NPZ, IUNITCD, NSCX, NSCY
     
      DOUBLE PRECISION, INTENT(OUT) ::
     .  ZREF, zref2, ZMIN, ZMAX, ARMUNI, EMIN, EMAX


C ****** INPUT *********************************************************
C REAL*8 VOLUME          : Volumen of unit cell (in bohr**3)
C ****** OUTPUT ********************************************************
C INTEGER IUNITCD        : Units for the electron density
C                          IUNITCD = 1 => Ele/(bohr)**3
C                          IUNITCD = 2 => Ele/(Ang)**3
C                          IUNITCD = 3 => Ele/(unitcell)
C INTEGER NPX, NPY, NPZ  : Number of points generated along x, y and z
C REAL*8 ZREF            : Position of reference plane for wf. extrapol.(Bohr)
C REAL*8 ZREF2            : Position of second reference plane.(Bohr)
C REAL*8  ZMIN, ZMAX     : Limits of the z-direction (in Bohr)
C INTEGER NSCX, NSCY     : Number of cells in x and y direction to plot
C                          in cube file
C REAL*8  EMAX, EMIN     : Energy limits for STM scan (en eV!!!)
C REAL*8  ARMUNI         : Conversion factors for the charge density
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER 
     .  UCD*22, UCD_DEFECT*22

      INTEGER
     .  IUNIT, NPX_DEFECT, NPY_DEFECT, NPZ_DEFECT

      LOGICAL 
     .  LEQI

      EXTERNAL 
     .  LEQI


C READ UNITS OF CHARGE DENSITY TO USE

      UCD_DEFECT = 'Ele/bohr**3'
      UCD = FDF_STRING('STM.DensityUnits',UCD_DEFECT)
      IF (LEQI(UCD,'ele/bohr**3')) then
        IUNITCD = 1
      ELSEIF (LEQI(UCD,'ele/ang**3')) then
        IUNITCD = 2
      ELSEIF (LEQI(UCD,'ele/unitcell')) then
        IUNITCD = 3
      ELSE
       WRITE(6,'(A)')' readstm: ERROR   Wrong Option in Units of      '
       WRITE(6,'(A)')' readstm:  Charge Density                       '
       WRITE(6,'(A)')' readstm:  You must choose one of the following:' 
       WRITE(6,'(A)')' readstm:                                       '
       WRITE(6,'(A)')' readstm:      - Ele/bohr**3                    '
       WRITE(6,'(A)')' readstm:      - Ele/ang**3                     '
       WRITE(6,'(A)')' readstm:      - Ele/unitcell                   '
       STOP
      ENDIF

      NPX_DEFECT = 50
      NPY_DEFECT = 50
      NPZ_DEFECT = 50
      NPX = FDF_INTEGER('STM.NumberPointsX',NPX_DEFECT)
      NPY = FDF_INTEGER('STM.NumberPointsY',NPY_DEFECT)
      NPZ = FDF_INTEGER('STM.NumberPointsZ',NPZ_DEFECT)

      IF ((2*(NPX/2) .NE. NPX) .OR. (NPX .LT. 2)) THEN
       WRITE(6,'(A)')' readstm: ERROR   NPX must be positive and even'
       STOP
      ENDIF

      IF ((2*(NPY/2) .NE. NPY) .OR. (NPY .LT. 2)) THEN
       WRITE(6,'(A)')' readstm: ERROR   NPY must be positive and even'
       STOP
      ENDIF

      EMIN = FDF_PHYSICAL('STM.Emin',-1.0d10,'eV')
      EMAX = FDF_PHYSICAL('STM.Emax',1.0d10,'eV')

      ZREF = FDF_PHYSICAL('STM.RefZ',1.0d40,'Bohr')
      ZREF2 = FDF_PHYSICAL('STM.RefZ2',ZREF,'Bohr')
      ZMIN = FDF_PHYSICAL('STM.MinZ',1.0d40,'Bohr')
      ZMAX = FDF_PHYSICAL('STM.MaxZ',1.0d40,'Bohr')

      IF (ZREF .GT. 0.99999d40) THEN
        WRITE(6,*) 'ERROR: You must specify STM.REFZ in input'
        STOP
      ENDIF
      IF (ZMIN .GT. 0.99999d40) THEN
        WRITE(6,*) 'ERROR: You must specify STM.MINZ in input'
        STOP
      ENDIF
      IF (ZMAX .GT. 0.99999d40) THEN
        WRITE(6,*) 'ERROR: You must specify STM.MAXZ in input'
        STOP
      ENDIF

C Number of cells to plot in cube file in x and y directions
      NSCX = FDF_INTEGER('STM.NumberCellsX',1)
      NSCY = FDF_INTEGER('STM.NumberCellsY',1)

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

