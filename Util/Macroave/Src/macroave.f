! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      PROGRAM MACROAVE
C **********************************************************************
C The MACROAVE program implements the macroscopic average technique,
C introduced by A. Baldereschi and coworkers
C (A. Baldereschi, S. Baroni, and R. Resta, Phys. Rev. Lett. 61, 734 (1988)).
C This is an extremely powerful method that relates
C microscopic quantities, typical outputs of first-principles codes,
C with macroscopic magnitudes, needed to perform electrostatic analysis.
C Within this methodology, we will be able of washing out all the
C wiggles of the rapidly-varying functions of position (resembling
C the underlying atomic structure) of the microscopic quantities,
C blowing up only the macroscopic features. 
C It can be used to compute band offsets, work functions, effective
C charges, and high frequency dielectric constants, among others
C interesting physical properties. 
C Ref: L. Colombo, R. Resta and S. Baroni  Phys Rev B  44, 5572 (1991)
C Coded by P. Ordejon and J. Junquera, April 1999
C Modified by J. Junquera, November 2001
C **********************************************************************

C Modules --------------------------------------------------------------
      USE DEFS_BASIS
      USE DEFS_COMMON
      use interpolation, only: spline, splint
      use m_fft_gpfa,only : fft_gpfa_ez

      IMPLICIT NONE

C ********* PARAMETERS *************************************************
C INTEGER   NP    : Parameter needed to define maximum number of points 
C                   for FFT grid.
C                   Number of points = (2**NP)
C INTEGER   N     : Maximum number of complex data point for FFT. 
C                   MUST be a power of 2
C REAL*8   HARTREE: Conversion factor from Hartrees to Rydbergs
C                   1 hartree = 2 Ry
C REAL*8   RYDBERG: Conversion factor from Rydbergs to eV
C                   1 Ry = 13.6058 eV
C **********************************************************************

      INTEGER NP, N
      PARAMETER(NP = 12)
      PARAMETER(N = 2**NP) 

      DOUBLE PRECISION HARTREE, RYDBERG
      PARAMETER(HARTREE = 2.D0)
      PARAMETER(RYDBERG = 13.6058D0)

C ********* VARIABLES **************************************************
      CHARACTER
     .   CODE*10, INTERP*10

      CHARACTER
     .   SNAME*15, INPDATA*15, FNAMERHO*24, 
     .   FNAMEPLAVE*26, FNAMEDELV*25

      LOGICAL
     .  SIESTA, ABINIT, POTENTIAL, CHARGE, TOTALCHARGE,
     .  FOUND

      logical :: linear =.false., splin=.false., poly=.false.

      INTEGER
     .   NATOMS, NSPIN, NSM, MESH(3), NCONV, NPOINTS, NPT, 
     .   NZ

      INTEGER 
     .  IS, IP, IX, II, I, IJ, IG, J, IA,
     .  UNIT1, UNIT2, UNIT3, UNIT4

      REAL, ALLOCATABLE :: 
     .  RHOS(:,:) 

      DOUBLE PRECISION, ALLOCATABLE :: 
     .  RHO(:,:) 

      DOUBLE PRECISION
     .   CELL(3,3), DCELL(3,3)

      DOUBLE PRECISION
     .   L, SUR, DS, LENGTH, CONVFAC, QREN, QTOT,
     .   LAV1, LAV2, DELTA, VOL, VOLCEL, SURPLA

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::
     .   Z, RHOZ, D2RHOZ, DRHODZ 
      
      DOUBLE PRECISION
     .   DATA(2*N), TH(2*N), V(2*N), X, GSQ, YP1, YPN,
     .   VREEC(N), VIMEC(N), RE(N), IM(N), PHI

      COMPLEX*8
     .  A, B, C
  
C ABINIT variables
      INTEGER
     .  NGFFT13(3), NSPPOL, FFORM, RDWR, UNITFI

      DOUBLE PRECISION
     .  RPRIMD(3,3)

      TYPE(HDR_TYPE) :: HDR
C end ABINIT variables

      EXTERNAL
     .   IO_ASSIGN, IO_CLOSE, THETAFT, VOLCEL 

C *********************************************************************
C CHARACTER CODE       : First principles-code used to generate the
C                        electrostatic potential at the points of a grid
C                        in real space. It is read from file 'macroave.in'
C CHARACTER SNAME      : System Label 
C                        (If code = ABINIT, then SNAME = FNAMERHO)
C CHARACTER INPDATA    : Calculate the band offset from the charge
C                        density or from the electrostatic potential?
C CHARACTER FNAMERHO   : Name of the file where the electrostatic
C                        potential at the mesh is stored
C CHARACTER FNAMEPLAVE : Name of the file where the planar average of the 
C                        electronic potential will be stored
C CHARACTER FNAMEDELV  : Name of the file where the profile
C                        of the electrostatic potential will be stored
C LOGICAL   SIESTA     : Have you used SIESTA to get the electrostatic potential
C                        or the charge density?  
C LOGICAL   ABINIT     : Have you used ABINIT to get the electrostatic potential
C                        or the charge density?  
C LOGICAL   LINEAR     : Linear interpolation to get the charge 
C                        density/potential in the FFT grid.
C LOGICAL   SPLIN      : Cubic spline interpolation to get the charge 
C                        density/potential in the FFT grid.
C LOGICAL   POTENTIAL  : We are going to compute the band offset from
C                        the electrostatic potential 
C LOGICAL   CHARGE     : We are going to compute the band offset from
C                        the charge density
C LOGICAL   FOUND      : Were data found? (only when task in iorho ='read')
C INTEGER   NATOMS     : Number of atoms in unit cell
C INTEGER   NSPIN      : Number of spin polarizations (1 or 2)
C INTEGER   MESH(3)    : Number of mesh divisions of each lattice vectors,
C                        INCLUDING subgrid
C INTEGER   NSM        : Number of sub-mesh points per mesh point
C                        (not used in this version)
C INTEGER   NPOINTS    : Number of mesh subdivisions in the normal direction 
C                        to the interface
C INTEGER   NCONV      : Number of convolutions required to calculate the
C                        macroscopic average
C INTEGER   NPT        : Total number of mesh points (included subpoints)
C REAL*8    CELL(3,3)  : Unit cell lattice vectors (a.u.) CELL(IXYZ,IVECT)
C REAL*8    DS         : Differential area per point of the mesh
C REAL*8    SUR        : Area of a plane parallel to the interface
C REAL*8    LENGTH     : Distance between two planes parallel to the interface
C REAL*8    L          : Length of the cell in the direction nomal 
C                        to the interface (a.u.)
C REAL*8    LAV1       : Linear period of the electrostatic potential
C                        in the bulklike region for the first material   
C REAL*8    LAV2       : Linear period of the electrostatic potential 
C                        in the bulklike region for the second material   
C REAL*8    CONVFAC    : Conversion factor for the output units
C REAL*8    QTOT       : Total electronic charge in the unit cell
C REAL*8    QREN       : Total electronic charge calculated from
C                        the input data file
C REAL*4    Z(MESH(3)) : Z-coordinate of the planes where the elec. density 
C                        is averaged
C REAL*4    RHO        : Electron density
C                        Notice single precision in this version
C REAL*8    RHOZ       : Planar average of the electron density
C REAL*8    DATA(2*N)  : Fourier coefficients of the planar average density
C REAL*8    TH(2*N)    : Fourier coefficients of the step functions
C REAL*8    V(2*N)     : Fourier coefficients of the potential
C REAL*4    VREEC(N)   : Real part of the electronic potential in real space
C REAL*4    VIMEC(N)   : Imag. part of the electronic potential in real space
C REAL*4    VRENC(N)   : Real part of the nuclear potential in real space
C REAL*4    VIMNC(N)   : Imag. part of the nuclear potential in real space
C *********************************************************************

C Reading input data from a file ---------------------------------------
      CALL IO_ASSIGN(UNIT1)
        OPEN(UNIT=UNIT1, FILE='macroave.in', STATUS='OLD')
          READ(UNIT1,'(A)')CODE
          READ(UNIT1,'(A)')INPDATA
          READ(UNIT1,'(A)')SNAME
          READ(UNIT1,*)NCONV
          READ(UNIT1,*)LAV1
          READ(UNIT1,*)LAV2
          READ(UNIT1,*)QTOT
          READ(UNIT1,*)INTERP
      CALL IO_CLOSE(UNIT1)

C Which code has been used to get the electrostatic potential? ---------
      IF ( CODE .EQ. 'siesta' .OR. CODE .EQ. 'SIESTA' .OR.
     .     CODE .EQ. 'Siesta' ) THEN
             SIESTA = .TRUE.
             ABINIT = .FALSE.
      ELSE IF ( CODE .EQ. 'abinit' .OR. CODE .EQ. 'ab-init' .OR.
     .          CODE .EQ. 'ABINIT' .OR. CODE .EQ. 'AB-INIT' .OR.
     .          CODE .EQ. 'Abinit' .OR. CODE .EQ. 'Ab-init') THEN
             SIESTA = .FALSE.
             ABINIT = .TRUE.
      ELSE 
        WRITE(6,*) 'macroave: Unknown code ', CODE 
        STOP
      ENDIF
             
C Are we going to compute the band offset from the charge density or
C from the electrostatic potential? ------------------------------------ 
      IF ( INPDATA .EQ. 'potential' .OR. INPDATA .EQ. 'POTENTIAL' .OR.
     .     INPDATA .EQ. 'Potential' ) THEN
             POTENTIAL   = .TRUE.
             CHARGE      = .FALSE.
             TOTALCHARGE = .FALSE.
      ELSE IF ( INPDATA .EQ. 'charge' .OR. INPDATA .EQ. 'CHARGE' .OR.
     .          INPDATA .EQ. 'Charge' ) THEN
             POTENTIAL   = .FALSE.
             CHARGE      = .TRUE. 
             TOTALCHARGE = .FALSE.
      ELSE IF ( INPDATA .EQ. 'totalcharge' .OR. 
     .          INPDATA .EQ. 'Totalcharge' .OR.
     .          INPDATA .EQ. 'TotalCharge' .OR.
     .          INPDATA .EQ. 'TOTALCHARGE' ) THEN
             POTENTIAL   = .FALSE.
             CHARGE      = .FALSE.
             TOTALCHARGE = .TRUE.
      ELSE 
        WRITE(6,*) 'macroave: Unknown input data  ', INPDATA
        STOP
      ENDIF

C What kind of interpolation will we use to get the charge density/
C potential in a FFT grid? ---------------------------------------------
          
      IF ( INTERP .EQ. 'linear' .OR. INTERP .EQ. 'Linear' .OR.
     .     INTERP .EQ. 'LINEAR' ) THEN
         LINEAR = .true.
      ELSE IF ( INTERP .EQ. 'spline' .OR. INTERP .EQ. 'Spline' .OR.
     .          INTERP .EQ. 'SPLINE' ) THEN
         SPLIN = .true.
      ! New option 'poly' with a cleaner implmentation
      ! It can replace 'linear' below, or provide a new order
      ! (interface to be decided)
      ELSE IF ( INTERP .EQ. 'poly') THEN
         POLY = .true.
      ENDIF

C Reading charge density from a file -----------------------------------
      IF ( SIESTA ) THEN
         IF (POTENTIAL) THEN
           FNAMERHO = TRIM(SNAME)//'.VH'
         ELSEIF (CHARGE) THEN 
           FNAMERHO = TRIM(SNAME)//'.RHO'
         ELSEIF (TOTALCHARGE) THEN 
           FNAMERHO = TRIM(SNAME)//'.TOCH'
         ENDIF
      ELSE IF ( ABINIT ) THEN
         FNAMERHO = SNAME
      ENDIF

      IF (SIESTA) THEN
        NSM   = 1 
        NPT   = 0
        NSPIN = 0
        CALL IORHO( 'READ', FNAMERHO, DCELL, MESH, NSM, NPT, NSPIN,
     .               RHOS, FOUND )
        IF (FOUND) THEN
          ALLOCATE( RHOS(NPT,NSPIN) )
          ALLOCATE( RHO(NPT,NSPIN) )
          CALL IORHO( 'READ', FNAMERHO, DCELL, MESH, NSM, NPT, NSPIN,
     .                 RHOS, FOUND )
          DO I = 1, 3
            DO J = 1, 3
              CELL(J,I) = DCELL(J,I)
            ENDDO
          ENDDO
C Transform the density or the potential read from SIESTA 
C from a single precision variable to a double precision variable
          DO IS = 1, NSPIN
            DO IP = 1, NPT
              RHO(IP,IS) = RHOS(IP,IS) * 1.0D0
            ENDDO  
          ENDDO 

        ELSE
          WRITE(6,*)'macroave: ERROR: file not found: ', FNAMERHO
          STOP
        ENDIF

      ELSE IF (ABINIT) THEN
        CALL IO_ASSIGN(UNIT2)
        OPEN(UNIT=UNIT2, FILE=FNAMERHO, FORM='UNFORMATTED',STATUS='OLD')
        RDWR = 1
        CALL HDR_IO(FFORM,HDR,RDWR,UNIT2)
C       For debugging 
c        rdwr=4 ; unitfi=6
c        call hdr_io(fform,hdr,rdwr,unitfi)

        DO I = 1, 3
          MESH(I) = HDR%NGFFT(I)
          DO J = 1, 3
            CELL(J,I) = HDR%RPRIMD(J,I)
          ENDDO
        ENDDO
        NSPIN = HDR%NSPPOL
        NPT = MESH(1) * MESH(2) * MESH(3)
        ALLOCATE( RHO(NPT,NSPIN) )

        DO IS = 1, NSPIN
          READ(UNIT2) (RHO(IP,IS),IP=1,NPT)
C Units for the potential in Ab-init are in Hartrees, 
C so we transform them into Ry. No transformation is
C needed for the charge density
C (it is directly read in electrons/bohr**3).
          IF (POTENTIAL) THEN 
            DO IP = 1, NPT 
              RHO(IP,IS) = RHO(IP,IS) * HARTREE
            ENDDO
          ENDIF
        ENDDO
        CALL IO_CLOSE(UNIT2)
      ENDIF

C Initialize some variables (we suppose orthorombic cells) -------------

      L  = CELL(3,3)
      SUR = SURPLA( CELL )
      VOL = VOLCEL( CELL )
      DS = SUR/( MESH(1) * MESH(2) )
      LENGTH = L/DBLE(N)
      NPOINTS = MESH(3)

      ALLOCATE(Z(NPOINTS+1))
      ALLOCATE(RHOZ(NPOINTS+1))
      ALLOCATE(D2RHOZ(NPOINTS+1))
      ALLOCATE(DRHODZ(N))

      RHOZ(1:NPOINTS+1)   = 0.D0
      D2RHOZ(1:NPOINTS+1) = 0.D0
      DRHODZ(1:N)         = 0.D0

      IF (POTENTIAL) THEN 
         CONVFAC = RYDBERG
      ELSE IF (CHARGE) THEN
         CONVFAC = 1.0D0
      ELSE IF (TOTALCHARGE) THEN
         CONVFAC = 1.0D0
      ENDIF 
  

C Loop over all points and calculate the planar average ----------------
C Warning: The planar average is only done for the first component of
C          RHO. Spin polarization is not implemented yet ---------------
      DO IP = 1, NPT
        NZ = (IP-1) / (MESH(1)*MESH(2)) + 1
        RHOZ(NZ) =  RHOZ(NZ) + RHO(IP,1)*DS
      ENDDO

      DO IP = 1, NPOINTS
        RHOZ(IP) = RHOZ(IP) / SUR
      ENDDO

      DO IP  = 1, NPOINTS
        Z(IP) = (IP-1)*CELL(3,3)/DBLE(NPOINTS)
      ENDDO 

C Calculate electrostatic potential or electronic charge density -------
C in fft grid, interpolating the planar average calculated before ------
      IF (SPLIN) THEN  
        Z(NPOINTS+1)    = L
        RHOZ(NPOINTS+1) = RHOZ(1)
        DELTA = L/DBLE(NPOINTS)
        YP1 = ( RHOZ(2) - RHOZ(NPOINTS) )   / (2.0D0*DELTA)
        YPN = YP1
        CALL SPLINE(DELTA, RHOZ, NPOINTS+1, YP1, YPN, D2RHOZ)
        I = 0
        DO II = 1, 2*N-1, 2
          I = I + 1
          X = (I-1)*L/DBLE(N)
          CALL SPLINT( DELTA, RHOZ, D2RHOZ, NPOINTS+1, X, DATA(II), 
     .                 DRHODZ(I) )
          DATA(II+1) = 0.D0
        ENDDO
      ELSE IF (LINEAR) THEN 
        I = 0
        DO II = 1,2*N-1,2
          I = I + 1
          X = (I-1)*L/DBLE(N)
          DO IJ = 1, NPOINTS
             IF (X .EQ. Z(IJ)) THEN
                DATA(II) = RHOZ(IJ)
                DATA(II+1) = 0.D0
                GOTO 20
             ENDIF
             IF (Z(IJ) .GT. X) THEN
                DATA(II) = RHOZ(IJ-1) +
     .                     (X-Z(IJ-1))*(RHOZ(IJ)-RHOZ(IJ-1))/
     .                     (Z(IJ)-Z(IJ-1))
                DATA(II+1) = 0.D0
                GOTO 20
             ENDIF
          ENDDO
          ! X > Z(NPOINTS)
          DATA(II)=RHOZ(NPOINTS) +
     .             (X-Z(NPOINTS))*(RHOZ(1)-RHOZ(NPOINTS))/
     .             (Z(NPOINTS)-Z(NPOINTS-1))
          DATA(II+1) = 0.D0
 20       CONTINUE
        ENDDO
      ELSE IF (POLY) THEN
        Z(NPOINTS+1)    = L
        RHOZ(NPOINTS+1) = RHOZ(1)
        I = 0
        do II = 1,2*N-1,2
           I = I + 1
           X = (I-1)*L/DBLE(N)
           ! quadratic polynomial for now
           ! interface to be decided
           call dpnint1(2,z,rhoz,npoints+1,x,data(ii),.true.)
           data(ii+1) = 0.d0
        enddo
      ENDIF

C  Renormalize the charge density ---------------------------------------
      IF (CHARGE .OR. TOTALCHARGE) THEN
        QREN = 0.D0
        DO IP = 1, 2*N-1, 2
          QREN = QREN + DATA(IP)*LENGTH*SUR
        ENDDO
        DO IP = 1, 2*N-1, 2
          IF (CHARGE) THEN
            DATA(IP) = DATA(IP) * QTOT/QREN
          ELSEIF(TOTALCHARGE) THEN
            DATA(IP) = DATA(IP) - QREN/VOL
          ENDIF
        ENDDO
        QREN = 0.D0
        DO IP = 1, 2*N-1, 2
          QREN = QREN + DATA(IP)*LENGTH*SUR
        ENDDO
C For debugging
c        WRITE(6,*)' QREN = ', QREN
      ENDIF 
C ...
 
C Print planar average of the electrostatic potential or ---------------
C the electronic charge density ----------------------------------------
      FNAMEPLAVE = TRIM(SNAME)//'.PAV'
      CALL IO_ASSIGN(UNIT3)
        OPEN(UNIT=UNIT3, FILE=FNAMEPLAVE,STATUS='UNKNOWN') 
          I = 0
          DO II = 1, 2*N-1, 2
            I = I+1
            X=(I-1)*L/DBLE(N)
c            WRITE(UNIT3,'(3F20.12)')X,
c     .           DATA(II)*CONVFAC,DATA(II+1)*CONVFAC
            WRITE(UNIT3,'(2F20.12)')X,
     .           DATA(II)*CONVFAC
          ENDDO
      CALL IO_CLOSE(UNIT3)
C ...


C Calculate fourier transform of the electrostatic potential or 
C the electronic density -----------------------------------------------
      CALL fft_gpfa_ez(DATA,N,1)       
C ...

C Calculate macroscopic average of the electrostatic potential or the
C electronic charge density taking the convolution with two step functions.
C In Fourier space, it is a product of the Fourier transform components -
C The decompositions in the sum over II is due to the spetial way in which
C the data are stored in the fft subroutine ( see Fig. 12.2.2, in 
C 'Numerical Recipes, The Art of Scientific Computing' 
C by W.H. Press, S.A. Teukolsky, W.T. Veterling and B.P. Flannery, 
C Cambridge U.P. 1987 and 1992.


      CALL THETAFT(N,L,LAV1,TH)

      DO II = 1, N+1, 2
        A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
        B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
        C = A*B
        DATA(II) = REAL(C)*L/DBLE(N)
        DATA(II+1) = AIMAG(C)*L/DBLE(N)
      ENDDO

      DO II = N+3, 2*N-1, 2
        A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
        B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
        C = A*B
        DATA(II) = REAL(C)*L/DBLE(N)
        DATA(II+1) = AIMAG(C)*L/DBLE(N)
      ENDDO


      IF (NCONV .EQ. 2) THEN  
        CALL THETAFT(N,L,LAV2,TH)

        DO II = 1, N+1, 2
          A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
          B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
          C = A*B
          DATA(II) = REAL(C)*L/DBLE(N)
          DATA(II+1) = AIMAG(C)*L/DBLE(N)
c          IF ( POISON ) THEN
c            IG = (II-1) / 2
c            GSQ= (2.D0*PI*IG/L)**2
c            IF(GSQ .GT. 0.D0) THEN
c              V(II) = DATA(II) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
c              V(II+1) = DATA(II+1) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
c            ELSE
c              V(II) = 0.D0
c              V(II+1) = 0.D0
c            ENDIF
c          ENDIF
        ENDDO

        DO II = N+3, 2*N-1, 2
          A = DATA(II)*(1.D0,0.D0) + DATA(II+1)*(0.D0,1.D0)
          B = TH(II)*(1.D0,0.D0) + TH(II+1)*(0.D0,1.D0)
          C = A*B
          DATA(II) = REAL(C)*L/DBLE(N)
          DATA(II+1) = AIMAG(C)*L/DBLE(N)
c          IF ( POISON ) THEN
c            IG = (-2*N+II-1) / 2
c            GSQ= (2.D0*PI*IG/L)**2
c            IF(GSQ .GT. 0.D0) THEN
c              V(II) = DATA(II) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
c              V(II+1) = DATA(II+1) * (4.D0*PI/GSQ) * HARTREE * RYDBERG
c            ELSE
c              V(II) = 0.D0
c              V(II+1) = 0.D0
c            ENDIF
c          ENDIF
        ENDDO

      ENDIF 
C ...

C Transform average electronic density and potential to real space -----
C The decompositions in the sum over J is due to the spetial way in which
C the data are stored in the fft subroutine ( see Fig. 12.2.2, in 
C 'Numerical Recipes, The Art of Scientific Computing' 
C by W.H. Press, S.A. Teukolsky, W.T. Veterling and B.P. Flannery, 
C Cambridge U.P. 1987 and 1992.

      DO II = 1, N
        RE(II) = 0.D0
        IM(II) = 0.D0
c        IF ( POISON ) THEN
c          VREEC(II) = 0.D0
c          VIMEC(II) = 0.D0
c        ENDIF
        DO J = 1, N+1, 2 
          PHI = -2.D0 * PI * (II-1) * ( (J-1)/2 ) / DBLE(N)
          RE(II)=RE(II)+(1.D0/DBLE(N))*(DATA(J)*COS(PHI)
     .                         -DATA(J+1)*SIN(PHI))
          IM(II)=IM(II)+(1.D0/DBLE(N))*(DATA(J)*SIN(PHI)
     .                         +DATA(J+1)*COS(PHI))
c          IF ( POISON ) THEN
c            VREEC(II)=VREEC(II)+(1.D0/DBLE(N))*(V(J)*COS(PHI)
c     .                           -V(J+1)*SIN(PHI))
c            VIMEC(II)=VIMEC(II)+(1.D0/DBLE(N))*(V(J)*SIN(PHI)
c     .                           +V(J+1)*COS(PHI))
c          ENDIF
        ENDDO

        DO J = N+3, 2*N-1, 2
          PHI = -2.0D0 * PI * (II-1) * ((-2*N+J-1)/2) / DBLE(N)
          RE(II)=RE(II)+(1.D0/DBLE(N))*(DATA(J)*COS(PHI)
     .                           -DATA(J+1)*SIN(PHI))
          IM(II)=IM(II)+(1.D0/DBLE(N))*(DATA(J)*SIN(PHI)
     .                         +DATA(J+1)*COS(PHI))
c          IF ( POISON ) THEN
c            VREEC(II)=VREEC(II)+(1.D0/DBLE(N))*(V(J)*COS(PHI)
c     .                           -V(J+1)*SIN(PHI))
c            VIMEC(II)=VIMEC(II)+(1.D0/DBLE(N))*(V(J)*SIN(PHI)
c     .                           +V(J+1)*COS(PHI))
c          ENDIF
        ENDDO
      ENDDO
C ...

C Print averaged electronic charge density and potential ---------------
      FNAMEDELV = TRIM(SNAME)//'.MAV'
      CALL IO_ASSIGN(UNIT4)
        OPEN(UNIT=UNIT4, FILE=FNAMEDELV, STATUS='UNKNOWN') 
          DO I = 1, N
            X=(I-1)*L/DBLE(N)
c            WRITE(UNIT4,'(3F20.5)')X,
c     .           RE(I)*CONVFAC ,IM(I)*CONVFAC
            WRITE(UNIT4,'(2F20.12)')X,
     .           RE(I)*CONVFAC 
          ENDDO
      CALL IO_CLOSE(UNIT4)
C ...

C Print electrostatic potential ----------------------------------------
c      IF (POISON) THEN
c        FNAMEVEC = TRIM(SNAME)//'.VEC'
c        CALL IO_ASSIGN(UNIT5)
c        OPEN(UNIT=UNIT5, FILE=FNAMEVEC, STATUS='UNKNOWN')
c        DO I = 1, N
c          X=(I-1)*L/DBLE(N)
c          WRITE(UNIT5,'(3F20.12)')X, VREEC(I), VIMEC(I)
c        ENDDO
c        CALL IO_CLOSE(UNIT5)
c      ENDIF
C ...

      DEALLOCATE(Z)
      DEALLOCATE(RHOZ)
      DEALLOCATE(DRHODZ)
      DEALLOCATE(D2RHOZ)

      END

