      SUBROUTINE READWAVES(NSPIN,NORB,IFLAG,NWF,NUMWF,NK,
     .                     RPSI,IPSI,E,K,IND)

C Reads the wavefunctions and energies from a file written by Siesta
C
C P. Ordejon, July 2003
C
C Modified for complex wavefunctions and multiple k-points,
C P. Ordejon, July 2004
C **************** INPUT ********************************************
C INTEGER NSPIN     : Number of spin components
C INTEGER NORB      : Number of basis orbitals
C INTEGER IFLAG     : 0=only read and return num of wavefcts and k-po.
C                     1=actually read wavefunctions
C **************** INPUT OR OUTPUT **********************************
C INTEGER NK        : Number of k-points to read 
C                     input(output) if IFLAG=0(1)
C INTEGER NUMWF     : Max num of wavefncts to read for a given k-po.
C INTEGER NWF(NK)   : Number of wavefunctions to read for each k-po.
C                     input(output) if IFLAG=0(1)
C **************** OUTPUT *******************************************
C REAL*8 RPSI(NORB,NK,NUMWF,NSPIN): Wavefunctions (real part)
C REAL*8 IPSI(NORB,NK,NUMWF,NSPIN): Wavefunctions (imag part)
C REAL*8 E(NK,NUMWF,NSPIN)        : Eigenvalues
C REAL*8 K(NK,3)                  : K-points
C INTEGER IND(NK,NUMWF)           : List of indexes of wavefunctions
C *******************************************************************

C Modules

      use precision
      use fdf


      IMPLICIT NONE

      INTEGER IFLAG, NSPIN, NORB, NUMWF, NK
      INTEGER IND(NK,NUMWF), NWF(NK)
      real(dp) RPSI(NORB,NK,NUMWF,NSPIN), 
     .                 IPSI(NORB,NK,NUMWF,NSPIN),
     .                 E(NK,NUMWF,NSPIN), K(NK,3)

C INTERNAL VARIABLES .............
      INTEGER UNIT, NSP, NUO, ISPIN, IISPIN, IWF, IIWF, IORB
      INTEGER IDUMB, NUMBERWF, NUMK, IK, IIK
      real(dp) REPSI,IMPSI

      CHARACTER PASTE*33
      CHARACTER, SAVE :: SNAME*30, FNAME*33
      CHARACTER CHDUMB*20

      SAVE UNIT
      EXTERNAL PASTE
C ..................


      IF (IFLAG .EQ. 0) THEN
        SNAME = FDF_STRING('SystemLabel','siesta')
        FNAME = PASTE(SNAME,'.WFS')

        CALL IO_ASSIGN(UNIT)
        OPEN (UNIT, FILE=FNAME, FORM='unformatted', STATUS='unknown')


        READ(UNIT) NK
C        IF (NK .NE. 1) THEN
C          WRITE(6,*) 'Wavefunctions file contains more then 1 k-point'
C          WRITE(6,*) 'DENCHAR can only handle the Gamma point!!'
C          STOP
C        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        DO IK = 1, NK
          READ(UNIT) 
          READ(UNIT) 
          READ(UNIT) NUMBERWF
          IF (NUMBERWF .GT. NUMWF) NUMWF = NUMBERWF
          DO IWF = 1, NUMBERWF
            READ(UNIT)
            READ(UNIT)
            DO IORB = 1, NORB
              READ(UNIT)
            ENDDO
          ENDDO
        ENDDO

        REWIND(UNIT)
        
        RETURN

      ELSE IF (IFLAG .EQ. 1) THEN

        READ(UNIT) NUMK
        IF (NK .NE. NUMK) THEN
          WRITE(6,*) 'Error in number of k-points'
          STOP
        ENDIF
        READ(UNIT) NSP
        IF (NSP .NE. NSPIN) THEN
          WRITE(6,*) 'NSPIN is not consistent between data files!'
          STOP
        ENDIF
        READ(UNIT) NUO
        IF (NUO .NE. NORB) THEN
          WRITE(6,*) 'Nr. of orbs is not consistent between data files!'
          STOP
        ENDIF

        DO IK = 1, NK
          DO ISPIN = 1,NSPIN
            READ(UNIT) IIK,K(IK,1),K(IK,2),K(IK,3)
            IF (IK .NE. IIK) THEN
              WRITE(6,*) 'Inconsistent order of kpoints in WFS file!'
              STOP
            ENDIF
            READ(UNIT) IISPIN
            IF (IISPIN .NE. ISPIN) THEN
              WRITE(6,*) 'Inconsistent order of spins in WFS file!'
              STOP
            ENDIF
            READ(UNIT) NWF(IK)
   

            DO IWF = 1,NWF(IK)
              READ(UNIT) IND(IK,IWF)
              READ(UNIT) E(IK,IWF,ISPIN)
              DO IORB = 1, NORB
                READ(UNIT) IDUMB,CHDUMB,IDUMB,IDUMB,IDUMB,CHDUMB,
     .                     REPSI,IMPSI
                RPSI(IORB,IK,IWF,ISPIN)=REPSI
                IPSI(IORB,IK,IWF,ISPIN)=IMPSI
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        CLOSE (UNIT)
        CALL IO_CLOSE(UNIT)
      ELSE
        WRITE(6,*) 'IFLAG must be either 0 or 1 in READWAVE!!'
        STOP
      ENDIF
           

      RETURN

      END
        

