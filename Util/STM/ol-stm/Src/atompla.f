      SUBROUTINE ATOMPLA( NA, ORIGIN, XA, MROT, IDIMEN, ISCALE,
     .                   NATINPLA, INDICES, XAPLANE )
C **********************************************************************
C Find the coordinates of some selected atoms on the 
C plane-reference-frame. Atoms in plane will have the third coordinate 
C equal to zero.
C For 3D-grid, all the atoms are selected
C Coded by J. Junquera May '99
C Modified by P. Ordejon to include 3D capability; June 2003
C **********************************************************************

      USE FDF
      USE PARSE
      USE SYS

      IMPLICIT NONE

      INTEGER
     .  NA, NATINPLA, INDICES(NA), ISCALE

      DOUBLE PRECISION
     .  ORIGIN(3), XA(3,NA), MROT(3,3), XAPLANE(3,NA)

C ******* INPUT ********************************************************
C INTEGER NA             : Number of atoms
C REAL*8  ORIGIN(3)      : Origin of the plane reference frame
C REAL*8  XA(3,NA)       : Atomic coordinates in lattice reference frame
C REAL*8  MROT(3,3)      : Rotation matrix from the lattice reference frame 
C                          to the in-plane reference frame
C INTEGER IDIMEN         : Determines if run is plane or 3D-grid
C INTEGER ISCALE         : Units of the atomic coordinates
C                          ISCALE = 1 => bohrs, ISCALE = 2 => Ang
C ******* OUTPUT *******************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)    : Atomic indices of the atoms whose coordinates
C                          will be rotated
C REAL*8  XAPLANE(3,NA)  : Atomic coordinates in plane reference frame
C **********************************************************************

      CHARACTER
     .  LINE*150

      INTEGER 
     .  NATINPL_DEFECT, IDIMEN, IUNIT, IAT, IX
 
      DOUBLE PRECISION
     .  VAUX1(3), VAUX2(3)

      LOGICAL 
     .  ATINPLA

      TYPE(PARSED_LINE), POINTER :: P
      TYPE(BLOCK), POINTER       :: BP

      EXTERNAL 
     .  MATVECT

C **********************************************************************
C INTEGER NATINPLA       : Number of atoms whose coordinates will be 
C                          rotated from the lattice reference frame to 
C                          the in-plane reference frame
C REAL*8  VAUX1(3)       : Auxiliar vector
C REAL*8  VAUX2(3)       : Auxiliar vector
C LOGICAL ATINPLA        : Is block Denchar.AtomsInPlane present in fdf?
C **********************************************************************

      IF (IDIMEN .EQ. 2) THEN

C Read fdf data block 'Denchar.AtomsInPlane' ---------------------------
        NATINPL_DEFECT = 0
        NATINPLA = 0
  
        NULLIFY(BP)
        IF ( .NOT. FDF_BLOCK('Denchar.AtomsInPlane',BP) )  GOTO 2000

        LOOP: DO
          IF (.NOT. FDF_BLINE(BP,LINE)) EXIT LOOP
          P=>DIGEST(LINE)
          IF (.NOT. MATCH(P,"I") ) 
     .         CALL DIE("Wrong format in Denchar.AtomsInPlane")
          NATINPLA = NATINPLA + 1
          INDICES(NATINPLA) = INTEGERS(P,1) 
          CALL DESTROY(P)
        ENDDO LOOP
        CALL DESTROY(BP) 
 2000   CONTINUE

      ELSE IF (IDIMEN .EQ. 3) THEN
        NATINPLA = NA
        DO IAT = 1, NA
          INDICES(IAT) = IAT
        ENDDO
      ELSE
        CALL DIE("Wrong IDIMEN in ATOMPLA")
      ENDIF
     

C Rotate the coordinates -----------------------------------------------
        DO IAT = 1, NATINPLA

          DO IX = 1,3
            VAUX1(IX) = XA(IX,INDICES(IAT)) - ORIGIN(IX)
          ENDDO

          CALL MATVECT(MROT, VAUX1, VAUX2)

          DO IX = 1,3
            XAPLANE(IX,INDICES(IAT)) = VAUX2(IX)
          ENDDO
  
C        IF (ISCALE .EQ. 2) THEN
C          DO IX = 1,3
C           XAPLANE(IX,INDICES(IAT))=XAPLANE(IX,INDICES(IAT))*0.529177D0
C          ENDDO
C        ENDIF

        ENDDO

        RETURN
        END
