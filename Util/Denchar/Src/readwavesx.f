! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE READWAVESX(NSPIN,NORB,IFLAG,NWF,NUMWF,NK,
     .                     RPSI,IPSI,E,K,IND)

C Reads the wavefunctions and energies from a WFSX file written by Siesta
C
C P. Ordejon, July 2003
C
C Modified for complex wavefunctions and multiple k-points,
C P. Ordejon, July 2004
C Conversion to WFSX format,
C A. Garcia, May 2012
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

      use precision, only: dp, sp
      use fdf


      IMPLICIT NONE

      INTEGER IFLAG, NSPIN, NORB, NUMWF, NK
      INTEGER IND(NK,NUMWF), NWF(NK)
      real(dp) RPSI(NORB,NK,NUMWF,NSPIN), 
     .                 IPSI(NORB,NK,NUMWF,NSPIN),
     .                 E(NK,NUMWF,NSPIN), K(NK,3)

C INTERNAL VARIABLES .............
      INTEGER UNIT, NSP, NUO, ISPIN, IISPIN, IWF, IIWF, j
      INTEGER IDUMB, NUMBERWF, NUMK, IK, IIK
      logical :: gamma_flag
      real(sp), allocatable :: psi(:,:)

      CHARACTER :: SNAME*30, FNAME*33

      SAVE UNIT
C ..................


      IF (IFLAG .EQ. 0) THEN

        NUMWF = 0   ! Initialize

        SNAME = FDF_STRING('SystemLabel','siesta')
        FNAME = trim(sname) // '.WFSX'

        CALL IO_ASSIGN(UNIT)
        OPEN (UNIT, FILE=FNAME, FORM='unformatted', STATUS='old')


        READ(UNIT) NK

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
        READ(UNIT)  ! Orbital labels

        DO IK = 1, NK
          READ(UNIT)   ! ik, k, weight
          READ(UNIT)   ! ispin
          READ(UNIT) NUMBERWF
          IF (NUMBERWF .GT. NUMWF) NUMWF = NUMBERWF
          DO IWF = 1, NUMBERWF
            READ(UNIT)  ! indwf
            READ(UNIT)  ! eigenvalue
            READ(UNIT)  ! actual wf data in a single record
          ENDDO
        ENDDO

        REWIND(UNIT)
        
        RETURN

      ELSE IF (IFLAG .EQ. 1) THEN

        READ(UNIT) NUMK, GAMMA_FLAG
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

        if (gamma_flag) then
         allocate(psi(1,nuo))
        else
         allocate(psi(2,nuo))
        endif

        read(UNIT)  ! Orbital labels

        DO IK = 1, NK
          DO ISPIN = 1,NSPIN
            READ(UNIT) IIK,K(IK,1),K(IK,2),K(IK,3)
            IF (IK .NE. IIK) THEN
              WRITE(6,*) 'Inconsistent order of kpoints in WFSX file!'
              STOP
            ENDIF
            READ(UNIT) IISPIN
            IF (IISPIN .NE. ISPIN) THEN
              WRITE(6,*) 'Inconsistent order of spins in WFSX file!'
              STOP
            ENDIF
            READ(UNIT) NWF(IK)
   

            DO IWF = 1,NWF(IK)
              READ(UNIT) IND(IK,IWF)
              READ(UNIT) E(IK,IWF,ISPIN)
              read(unit) (psi(1:,j), j=1,nuo)
              if (gamma_flag) then
                 RPSI(1:nuo,IK,IWF,ISPIN)=dble(psi(1,1:nuo))
                 IPSI(1:nuo,IK,IWF,ISPIN)=0.d0
              else
                 RPSI(1:nuo,IK,IWF,ISPIN)=dble(psi(1,1:nuo))
                 IPSI(1:nuo,IK,IWF,ISPIN)=dble(psi(2,1:nuo))
              endif
            ENDDO
          ENDDO
        ENDDO

        CLOSE (UNIT)
        CALL IO_CLOSE(UNIT)
      ELSE
        WRITE(6,*) 'IFLAG must be either 0 or 1 in READWAVESX!!'
        STOP
      ENDIF
           

      RETURN

      END
        

