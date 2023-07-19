
      SUBROUTINE WROUT(IDIMEN, CHARGE, WAVES, IOPTION, NORMAL, COORPO, 
     .                 DIRVER1, DIRVER2, 
     .                 NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, 
     .                 ZMIN, ZMAX, IUNITCD,
     .                 MAXATOM, NAPLA, INDICES, XAPLA )

C **********************************************************************
C Dump input data to ouput
C Modified to make general writeout is done only to standard output, 
C not to the individual data files.
C
C Written by J. Junquera Feb '99
C Modified by P. Ordejon, June 2003
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      LOGICAL
     .  CHARGE, WAVES

      INTEGER
     .  IOPTION,  NPX, NPY, NPZ, IUNITCD, MAXATOM, NAPLA, IDIMEN,
     .  INDICES(MAXATOM)

      DOUBLE PRECISION
     .  NORMAL(3), COORPO(3,3), DIRVER1(3), DIRVER2(3), XAPLA(3,MAXATOM)

      DOUBLE PRECISION
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX

C **************  INPUT  ***********************************************
C INTEGER IDIMEN         : 2D or 3D run
C LOGICAL CHARGE         : Are we writting charge output?
C LOGICAL WAVES          : Are we writting wavefunctions output?
C INTEGER IDIMEN         : 2D or 3D run
C INTEGER IOPTION        : Option to generate the plane
C                          1 = Normal vector
C                          2 = Two vectors belonging to the plane

C                          3 = Three points of the plane
C                          4 = Three atomic indices
C REAL*8  NORMAL(3)      : Components of the normal vector
C REAL*8  COORPO(3,3)    : Coordinates of the three points used to define
C                          the plane 
C REAL*8  DIRVER(3)      : Components of two vector contained in the plane
C INTEGER NPX,NPY,NPZ    : Number of points along x and y and z 
C REAL*8  XMIN, XMAX     : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX     : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX     : Limits of the grid z-direction
C INTEGER IUNITCD        : Unit of the charge density
C INTEGER MAXATOM        : Total number of atoms in supercell
C INTEGER NAPLA          : Number of atoms whose coordiantes has been rotated   
C INTEGER INDICES(MAXATOM): Indices of tha atoms whose coordinates has 
C                           been roated
C REAL*8  XAPLA(3,MAXATOM): Atomic coordiantes in the in-plane reference frame
C **********************************************************************

C ***************  INTERNAL VARIABLES **********************************
      CHARACTER
     .  SNAME*30

      INTEGER
     .  IX, IP, IA, UNIT1

      LOGICAL, SAVE :: FRSTME

      DATA FRSTME /.TRUE./



C Open files to store charge density -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')

      UNIT1 = 6


C Write general information only if called for the first time

      IF (FRSTME) THEN


        WRITE(UNIT1,'(A)')
     .    '                          ************************       '
        WRITE(UNIT1,'(A)')
     .    '                          *  WELCOME TO DENCHAR  *       '
        WRITE(UNIT1,'(A)')
     .    '                          ************************       '

        WRITE(UNIT1,'(A,A)')
     .    '  You are running DENCHAR for system: ',SNAME
        WRITE(UNIT1,'(A)')
     .    '  '

        WRITE(UNIT1,'(A)')
        IF (IDIMEN .EQ. 2) THEN
          WRITE(UNIT1,'(A)')
     .   '  You have chosen the 2D mode. Values of the functions'
          WRITE(UNIT1,'(A)')
     .    '  will be given in a 2D grid'
        ELSE IF (IDIMEN .EQ. 3) THEN
          WRITE(UNIT1,'(A)')
     .   '  You have chosen the 3D mode. Values of the functions'
          WRITE(UNIT1,'(A)')
     .    '  will be given in a 3D grid, in Gaussian Cube format'
        ENDIF

        WRITE(UNIT1,'(A)')
        WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the x-direction : ',
     .    '  ', NPX
        WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the y-direction : ',
     .    '  ', NPY
        IF (IDIMEN .EQ. 3)
     .    WRITE(UNIT1,'(A,/,A,I5)')
     .    '  Number of points in the z-direction : ',
     .    '  ', NPZ
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the x-component of the window : ',
     .    '  ', XMIN,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the x-component of the window : ',
     .    '  ', XMAX,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the y-component of the window : ',
     .    '  ', YMIN,' bohrs'
        WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the y-component of the window : ',
     .    '  ', YMAX,' bohrs'
        IF (IDIMEN .EQ. 3) THEN
          WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Minimum value of the z-component of the window : ',
     .    '  ', ZMIN,' bohrs'
          WRITE(UNIT1,'(A,/,A,F12.5,A)')
     .    '  Maximum value of the z-component of the window : ',
     .    '  ', ZMAX,' bohrs'
        ENDIF

        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A,/A)')
     .    '  The options you have chosen to generate the plane',
     .    '  are the following: '

        IF( IOPTION .EQ. 1 ) THEN

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : NormalVector'
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the normal vector : ',
     .    '  ',(NORMAL(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Origin of the plane : ',
     .    '  ',(COORPO(1,IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Another point to define the X direction : ',
     .    '  ',(COORPO(2,IX),IX=1,3)

        ELSEIF( IOPTION .EQ. 2 ) THEN 
        
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : TwoLines'
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the first vector inside the plane :',
     .    '  ',(DIRVER1(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Components of the second vector inside the plane:',
     .    '  ',(DIRVER2(IX),IX=1,3)
          WRITE(UNIT1,'(A,/,A,3F12.5)')
     .    '  Origin of the plane : ',
     .    '  ',(COORPO(1,IX),IX=1,3)

        ELSEIF( IOPTION .EQ. 3 ) THEN 

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : ThreePoints'
          WRITE(UNIT1,'(A)')
     .    '  Coordinates of three points in the plane : '
          DO IP = 1,3
            WRITE(UNIT1,'(A,3F12.5)')
     .    '  '  ,(COORPO(IP,IX),IX=1,3)
          ENDDO
      
        ELSEIF( IOPTION .EQ. 4 ) THEN 

          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Option to generate the plane : ThreeAtomicIndices'
          WRITE(UNIT1,'(A)')
     .    '  Position of the three atoms : '
          DO IP = 1,3
            WRITE(UNIT1,'(A,3F12.5)')
     .    '  '  ,(COORPO(IP,IX),IX=1,3)
          ENDDO

        ENDIF

        IF ( IUNITCD .EQ. 1) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/(bohr**3)'
        ELSEIF ( IUNITCD .EQ. 2) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/(angstrom**3)'
        ELSEIF( IUNITCD .EQ. 3) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A,/,A)')
     .    '  Unit of the charge density in output files : ',
     .    '  Electrons/unit cell'
        ENDIF


        IF( NAPLA .NE. 0) THEN
          WRITE(UNIT1,'(A)')
     .    '  '
          WRITE(UNIT1,'(A)')
     .    '  Atomic coordinates in the in-plane reference frame'
          WRITE(UNIT1,'(A,19(1H ),A)')
     .    '  Atomic Index','Atomic coordinates'
          DO IA = 1, NAPLA
            WRITE(UNIT1,'(A,I14,5X,3F15.4)')
     .      '',INDICES(IA), (XAPLA(IX,INDICES(IA)),IX=1,3)
          ENDDO
        ENDIF

        FRSTME = .FALSE.

      ENDIF

      IF (CHARGE) THEN
        
        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A)')
     .    '  You are now computing charge density on the grid'

      ENDIF

      IF (WAVES) THEN
        
        WRITE(UNIT1,'(A)')
     .    '  '
        WRITE(UNIT1,'(A)')
     .    '  You are now computing Wave Functions on the grid'

      ENDIF


      END
       
      
