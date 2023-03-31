! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE WAVOFR( NA, NO, no_u, MAXNA, NSPIN, nspin_blocks,
     .                   non_coll,
     .                   ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .                   wf_unit, NK, gamma_wfsx,
     .                   IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .                   ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .                   DIRVER1, DIRVER2, 
     .                   ARMUNI, IUNITCD, ISCALE, RMAXO,
     $                   kpoint_selected, wf_selected )
C **********************************************************************
C Compute the wave functions at the points of a plane or a 3D grid
C in real space
C     Coded by P. Ordejon, from Junquera's rhoofr. July 2003
!     Extended to NC spin by Alberto Garcia, 2019      
C **********************************************************************

      use precision
      USE FDF
      USE ATMFUNCS, only: phiatm
      USE CHEMICAL, only: atomic_number
      use planed, only: plane

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  NA, NO, NO_U, IOPTION, NPX, NPY, NPZ, ISCALE, IUNITCD,
     .  IDIMEN, NSPIN, nspin_blocks, MAXNA, NK,
     .  ISA(NA), IPHORB(NO), INDXUO(NO), LASTO(0:NA)
      real(dp), INTENT(IN) ::
     .  XA(3,NA), XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3), DIRVER2(3), ARMUNI,
     .  RMAXO
      logical, intent(in) :: non_coll, gamma_wfsx
      integer, intent(in) :: wf_unit
      integer, intent(in) :: kpoint_selected, wf_selected
      real(dp), INTENT(IN) :: CELL(3,3)


C ****** INPUT *********************************************************
C INTEGER NA               : Total number of atoms in Supercell
C INTEGER NO               : Total number of orbitals in Supercell
C INTEGER NO_U             : Total number of orbitals in Unit Cell
C INTEGER MAXNA            : Maximum number of neighbours of any atom
C INTEGER NSPIN            : Number of different spin polarizations: 1, 2, 4, 8
! integer nspin_blocks     : blocks in WFSX file       
! logical non_coll         : NC/SOC wfs data?
C INTEGER ISA(NA)          : Species index of each atom
C INTEGER IPHORB(NO)       : Orital index of each orbital in its atom
C INTEGER INDXUO(NO)       : Equivalent orbital in unit cell
C INTEGER LASTO(0:NA)      : Last orbital of each atom in array iphorb
C REAL*8  XA(3,NA)         : Atomic positions in cartesian coordinates
C                            (in bohr)
C REAL*8  CELL(3,3)        : Supercell vectors CELL(IXYZ,IVECT)
C                            (in bohr)
C INTEGER NK               : Number of k-points
C INTEGER IDIMEN           : Specify if the run is to plot quantities
C                            in a plane or in a 3D grid (2 or 3, respect)
C INTEGER IOPTION          : Option to generate the plane
C                            1 = Normal vector
C                            2 = Two vectors contained in the plane
C                            3 = Three points in the plane
C INTEGER NPX,NPY,NPZ      : Number of points along x and y and z
C REAL*8  XMIN, XMAX       : Limits of the plane in the PRF for x-direction
C REAL*8  YMIN, YMAX       : Limits of the plane in the PRF for y-direction
C REAL*8  ZMIN, ZMAX       : Limits of the z-direction
C REAL*8  COORPO(3,3)      : Coordinates of the three points used to define
C                            the plane (Only used if ioption = 3)
C REAL*8  NORMAL(3)        : Components of the normal vector used to define 
C                            the plane
C REAL*8  DIRVER1(3)       : Components of the first vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C REAL*8  DIRVER2(3)       : Components of the second vector contained 
C                            in the plane
C                            (Only used if ioption = 2)
C INTEGER IUNITCD          : Unit of the charge density
C INTEGER ISCALE           : Unit if the points of the plane
C REAL*8  RMAXO            : Maximum range of basis orbitals
C **********************************************************************

      logical :: filtering_k_points, filtering_wfs
      
      INTEGER  NPLAMAX, NAPLA, NAINCELL
      INTEGER, DIMENSION(:), ALLOCATABLE ::  INDICES, JNA

      REAL(SP), ALLOCATABLE :: wf_single(:,:)
      COMPLEX(DP), ALLOCATABLE :: wf(:,:)
      complex(dp), allocatable :: CWAVE(:)

      real(dp) :: k(3)
      REAL, DIMENSION(:,:), allocatable :: RWF, IMWF, MWF, PWF

      real(dp), DIMENSION(:), ALLOCATABLE ::  R2IJ
      real(dp), DIMENSION(:,:), ALLOCATABLE ::
     .   PLAPO, POINRE, XAPLA, XIJ, XAINCELL

      INTEGER
     .  NPO, IA, ISEL, NNA, UNITRE1, UNITRE2, UNITIM1, UNITIM2,
     .  UNITPH1, UNITPH2, UNITMO1, UNITMO2
 
      INTEGER
     .  I, J, IN, IAT1, IAT2, JO, IO, IUO, IAVEC, IAVEC1, IAVEC2, 
     .  IS1, IS2, IPHI1, IPHI2, IND, IX, IY, IZ, NX, NY, NZ, IWF, 
     .  INDWF, IZA(NA), IK, is

      integer :: spinor_comps, ispin, iwf_orig
      integer :: idummy, number_of_wfns, spinor_dim
      real(dp) :: ener
      complex(dp) :: expphi
      
      real(dp)
     .  RMAX, XPO(3), RMAX2, XVEC1(3),
     .  XVEC2(3), PHIMU, GRPHIMU(3),
     .  OCELL(3,3), PHASE, SI, CO, PI

      real(dp)
     .  RWAVE, RWAVEUP, RWAVEDN,
     .  IWAVE, IWAVEUP, IWAVEDN,
     .  MWAVE, MWAVEUP, MWAVEDN,
     .  PWAVE, PWAVEUP, PWAVEDN
 
      LOGICAL FIRST, empty_domain

      CHARACTER
     .  SNAME*40, FN_RE*60, FN_IM*60, 
     .  FN_URE*60, FN_UIM*60, FN_DRE*60, FN_DIM*60, 
     .  FN_MO*60, FN_PH*60,
     .  FN_UMO*60, FN_UPH*60, FN_DMO*60, FN_DPH*60,
     .  CHAR1*10, CHAR2*10, EXT*20, EXT2*25

      EXTERNAL IO_ASSIGN, IO_CLOSE, NEIGHB, WROUT

C **********************************************************************
C INTEGER NPLAMAX          : Maximum number of points in the plane
C REAL*8  PLAPO(NPLAMAX,3) : Coordinates of the points of the plane in PRF
C REAL*8  POINRE(NPLAMAX,3): Coordinates of the points of the plane in Lattice
C                            Reference Frame
C INTEGER NAPLA            : Number of atoms whose coordinates will be rotated
C INTEGER INDICES(NA)      : Indices of the atoms whose coordinates will 
C                            be rotated from the lattice reference frame 
C                            to the in-plane reference frame
C REAL*8  XAPLA(3,NA)      : Atomic coordinates in plane reference frame
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
C REAL RWF(NPO,NSPIN)      : Wave fnctn at each point of the grid (real part)
C REAL IMWF(NPO,NSPIN)     : Wave fnctn at each point of the grid (imag part)
C INTEGER IZA(NA)          : Atomic number of each atom
C **********************************************************************

      filtering_k_points = (kpoint_selected /= 0 )
      filtering_wfs = (wf_selected /= 0 )
      
C     Allocate some variables ---------------------------------------------

      PI = 4.0D0 * ATAN(1.0D0)

      select case ( NSPIN )
      case ( 1, 2, 4, 8 )
! ok
      case default
        WRITE(6,*)'BAD NUMBER NSPIN IN WAVOFR.F'
        WRITE(6,*)'NSPIN = ',NSPIN
        WRITE(6,*)'IT MUST BE 1, 2, 4 or 8'
        STOP
      end select


      NPLAMAX = NPX * NPY * NPZ

      ALLOCATE(PLAPO(NPLAMAX,3))
      ALLOCATE(POINRE(NPLAMAX,3))
      ALLOCATE(INDICES(NA))
      ALLOCATE(XAPLA(3,NA))

      ALLOCATE(XAINCELL(3,NA))

C Build the plane ------------------------------------------------------
          CALL PLANE( NA, NPLAMAX, IDIMEN, IOPTION, 
     .                XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, 
     .                NPX, NPY, NPZ, COORPO, NORMAL, 
     .                DIRVER1, DIRVER2, 
     .                XA, NAPLA, INDICES, ISCALE,
     .                POINRE, PLAPO, XAPLA )   

C Initialize neighbour subroutine --------------------------------------
      IA = 0
      ISEL = 0
      RMAX = RMAXO
      NNA  = MAXNA
      IF (ALLOCATED(JNA)) THEN
        DEALLOCATE(JNA)
      ENDIF
      IF (ALLOCATED(R2IJ)) THEN
        DEALLOCATE(R2IJ)
      ENDIF
      IF (ALLOCATED(XIJ)) THEN
        DEALLOCATE(XIJ)
      ENDIF
      ALLOCATE(JNA(MAXNA))
      ALLOCATE(R2IJ(MAXNA))
      ALLOCATE(XIJ(3,MAXNA))

      FIRST = .TRUE.
      DO I = 1,3
        XPO(I) = 0.D0
      ENDDO

      CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .             NNA, JNA, XIJ, R2IJ, FIRST )
 
      FIRST = .FALSE.
      RMAX2 =  RMAXO**2

! Stream over wavefunctions in file

      ! The first dimension of wf_single is the number of real numbers per orbital
      ! to be read from the WFSX file:
      ! 1 for real wfs, 2 for complex, and four for the two spinor components
      ! wf is a complex array which holds either a wfn or a two-component spinor.

      if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(1:no_u,2))
        spinor_comps = 2
      else
        spinor_comps = 1
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           allocate(wf(1:no_u,1))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(1:no_u,1))
        endif
      endif
      allocate(CWAVE(spinor_comps))

      k_loop: DO IK  = 1, NK
         ! ** Put filters in these loops
         spin_loop: do ispin = 1, nspin_blocks

         read(wf_unit) idummy, k(1:3)
            if (idummy /= ik) then
               write(6,*) "ik index mismatch in WFS file"
               WRITE(6,*) "ik in file, ik: ", idummy, ik
            endif
         read(wf_unit) idummy
            if (idummy /= ispin) then
               write(6,*) "ispin index mismatch in WFS file"
               WRITE(6,*) "ispin in file, ispin: ", idummy, ispin
            endif
         read(wf_unit) number_of_wfns

         if (filtering_k_points .and. (ik /= kpoint_selected)) then
            write(0,*) 'k-point ', ik, ' skipped'
            do IWF = 1, number_of_wfns
               read(wf_unit) 
               read(wf_unit) 
               read(wf_unit) 
            enddo
            CYCLE k_loop
         endif
         
         WRITE(6,"(a,i0,a,i0)") 'Processing kpoint ',IK,
     $              ' nwf: ', number_of_wfns
         WRITE(6,*) '     --------------------------------'

         wf_loop: DO IWF = 1, number_of_wfns

            read(wf_unit) iwf_orig
            ! Note that we mean the *original* index
            if (filtering_wfs .and. (iwf_orig /= wf_selected)) then
               read(wf_unit) 
               read(wf_unit)
               CYCLE wf_loop
            endif
               
            if (iwf_orig /= iwf) then
               ! The file holds a subset of wfs, with the original indexes...
               WRITE(6,*) 'Original wf index: ', iwf_orig
            endif
            read(wf_unit) ener

            read(wf_unit) (wf_single(:,io), io=1,no_u)
            ! Use a double precision complex form in what follows
            if ( non_coll) then
               wf(:,1) = cmplx(wf_single(1,:), wf_single(2,:), kind=dp)
               wf(:,2) = cmplx(wf_single(3,:), wf_single(4,:), kind=dp)
            else
               if (gamma_wfsx) then
                  wf(:,1) = cmplx(wf_single(1,:), 0.0_sp, kind=dp)
               else
                  wf(:,1) = cmplx(wf_single(1,:),wf_single(2,:),kind=dp)
               endif
            endif

            write(CHAR1,"(i0)") Iwf_Orig
            write(CHAR2,"(i0)") IK
            
!     Open files to store wave functions -----------------
            SNAME = FDF_STRING('SystemLabel','siesta')

         IF (NSPIN .EQ. 1) THEN

          IF (IDIMEN .EQ. 2) THEN
            FN_RE = TRIM(SNAME)//'.CON.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FN_PH = TRIM(FN_RE)//'.PHASE'
            FN_MO = TRIM(FN_RE)//'.MOD'
            FN_IM = TRIM(FN_RE)//'.IMAG'
            FN_RE = TRIM(FN_RE)//'.REAL'
          ELSEIF (IDIMEN .EQ. 3) THEN
            FN_RE = TRIM(SNAME)//'.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FN_PH = TRIM(FN_RE)//'.PHASE.cube'
            FN_MO = TRIM(FN_RE)//'.MOD.cube'
            FN_IM = TRIM(FN_RE)//'.IMAG.cube'
            FN_RE = TRIM(FN_RE)//'.REAL.cube'
          ENDIF

          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FN_RE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FN_IM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FN_MO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FN_PH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)

       ELSE ! nspin >= 2
          ! We will reuse 'up' and 'down' for the spinor components
           IF (IDIMEN .EQ. 2) THEN
            FN_URE = TRIM(SNAME)//'.CON.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FN_DRE = TRIM(FN_URE)//'.DOWN'
            FN_DPH = TRIM(FN_DRE)//'.PHASE'
            FN_DMO = TRIM(FN_DRE)//'.MOD'
            FN_DIM = TRIM(FN_DRE)//'.IMAG'
            FN_DRE = TRIM(FN_DRE)//'.REAL'

            FN_URE = TRIM(FN_URE)//'.UP'
            FN_UPH = TRIM(FN_URE)//'.PHASE'
            FN_UMO = TRIM(FN_URE)//'.MOD'
            FN_UIM = TRIM(FN_URE)//'.IMAG'
            FN_URE = TRIM(FN_URE)//'.REAL'

           ELSE IF (IDIMEN .EQ. 3) THEN
            FN_URE = TRIM(SNAME)//'.K' //
     $            trim(char2) // '.WF.' // trim(CHAR1)
            FN_DPH = TRIM(FN_URE)//'.DOWN.PHASE.cube'
            FN_DMO = TRIM(FN_URE)//'.DOWN.MOD.cube'
            FN_DIM = TRIM(FN_URE)//'.DOWN.IMAG.cube'
            FN_DRE = TRIM(FN_URE)//'.DOWN.REAL.cube'
            FN_UPH = TRIM(FN_URE)//'.UP.PHASE.cube'
            FN_UMO = TRIM(FN_URE)//'.UP.MOD.cube'
            FN_UIM = TRIM(FN_URE)//'.UP.IMAG.cube'
            FN_URE = TRIM(FN_URE)//'.UP.REAL.cube'
           ENDIF

          CALL IO_ASSIGN(UNITRE1)
          OPEN(UNIT = UNITRE1, FILE = FN_URE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE1)
          CALL IO_ASSIGN(UNITRE2)
          OPEN(UNIT = UNITRE2, FILE = FN_DRE, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITRE2)
          CALL IO_ASSIGN(UNITIM1)
          OPEN(UNIT = UNITIM1, FILE = FN_UIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM1)
          CALL IO_ASSIGN(UNITIM2)
          OPEN(UNIT = UNITIM2, FILE = FN_DIM, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITIM2)
          CALL IO_ASSIGN(UNITMO1)
          OPEN(UNIT = UNITMO1, FILE = FN_UMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO1)
          CALL IO_ASSIGN(UNITMO2)
          OPEN(UNIT = UNITMO2, FILE = FN_DMO, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITMO2)
          CALL IO_ASSIGN(UNITPH1)
          OPEN(UNIT = UNITPH1, FILE = FN_UPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH1)
          CALL IO_ASSIGN(UNITPH2)
          OPEN(UNIT = UNITPH2, FILE = FN_DPH, STATUS = 'UNKNOWN',
     .         FORM = 'FORMATTED')
          REWIND(UNITPH2)
        ENDIF

        IF (IDIMEN .EQ. 2) THEN
C Select all atoms in list to be printed out
          NAINCELL=NAPLA
          DO IA=1,NAPLA
            DO IX=1,3
              XAINCELL(IX,IA)=XAPLA(IX,IA)
            ENDDO
          ENDDO

        ELSE IF (IDIMEN .EQ.3) THEN
          DO IX = 1,3
            DO IY = 1,3
              OCELL(IX,IY)=0.D0
            ENDDO
          ENDDO
C     Determine cell size
          ! Note that we are *always* using an orthorhombic box!
          OCELL(1,1) = DABS(XMAX-XMIN)
          OCELL(2,2) = DABS(YMAX-YMIN)
          OCELL(3,3) = DABS(ZMAX-ZMIN)
C   Determine atoms which are within the plotting box
          NAINCELL=0
          DO IA=1,NA
            IF ((XAPLA(1,IA).LT.XMIN*1.1).OR.(XAPLA(1,IA).GT.XMAX*1.1)
     .     .OR. (XAPLA(2,IA).LT.YMIN*1.1).OR.(XAPLA(2,IA).GT.YMAX*1.1)
     .     .OR. (XAPLA(3,IA).LT.ZMIN*1.1).OR.(XAPLA(3,IA).GT.ZMAX*1.1))
     .      GOTO 90
            NAINCELL=NAINCELL+1
            IZA(NAINCELL) = ATOMIC_NUMBER(ISA(IA))
            DO IX=1,3
              XAINCELL(IX,NAINCELL)=XAPLA(IX,IA)
            ENDDO
90          CONTINUE
          ENDDO

          IF (NSPIN .EQ. 1) THEN
            call write_cube_header(unitre1,fn_re)
            call write_cube_header(unitim1,fn_im)
            call write_cube_header(unitmo1,fn_mo)
            call write_cube_header(unitph1,fn_ph)
      
          ELSE IF (NSPIN .EQ. 2) THEN
             if (ispin == 1 ) then
                call write_cube_header(unitre1,fn_ure)
                call write_cube_header(unitim1,fn_uim)
                call write_cube_header(unitmo1,fn_umo)
                call write_cube_header(unitph1,fn_uph)
             else
                call write_cube_header(unitre2,fn_dre)
                call write_cube_header(unitim2,fn_dim)
                call write_cube_header(unitmo2,fn_dmo)
                call write_cube_header(unitph2,fn_dph)
             endif
          ELSE ! 4 or 8
            call write_cube_header(unitre1,fn_ure)
            call write_cube_header(unitim1,fn_uim)
            call write_cube_header(unitmo1,fn_umo)
            call write_cube_header(unitph1,fn_uph)
            
            call write_cube_header(unitre2,fn_dre)
            call write_cube_header(unitim2,fn_dim)
            call write_cube_header(unitmo2,fn_dmo)
            call write_cube_header(unitph2,fn_dph)

          ENDIF
        ENDIF


        CALL WROUT(IDIMEN, .FALSE., .TRUE., IOPTION, NORMAL, COORPO,
     .             DIRVER1, DIRVER2,
     .             NPX, NPY, NPZ, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .             IUNITCD,  NA, NAINCELL, INDICES, XAINCELL )

        WRITE(6,'(A)')
        WRITE(6,'(A,i0,a,i0)')
     .    '   Generating grid values for wf... k-point:',IK,
     $       ' wf #:', iwf_orig


C Allocate space for wave functions in 3D-grid

        spinor_dim = max(nspin,2)
        IF (.NOT.ALLOCATED(RWF)) THEN
          ALLOCATE(RWF(NPX*NPY*NPZ,spinor_dim))
          ALLOCATE(IMWF(NPX*NPY*NPZ,spinor_dim))
          ALLOCATE(MWF(NPX*NPY*NPZ,spinor_dim))
          ALLOCATE(PWF(NPX*NPY*NPZ,spinor_dim))
        ENDIF

      
C Loop over all points in real space -----------------------------------

        empty_domain = .true.
        NPO = 0
        DO 102 NZ = 1,NPZ
        DO 101 NY = 1,NPY
        DO 100 NX = 1,NPX
          NPO = NPO + 1

          ! Initialize the wave function at each point -----------------------
          CWAVE(:) = 0.0_dp

C Localize non-zero orbitals at each point in real space ---------------
          DO IX = 1,3
            XPO(IX) = POINRE(NPO,IX)
          ENDDO
     
          IA   = 0
          ISEL = 0
          NNA  = MAXNA

          CALL NEIGHB( CELL, RMAX, NA, XA, XPO, IA, ISEL, 
     .                 NNA, JNA, XIJ, R2IJ, FIRST )

          if (nna /= 0) empty_domain = .false.

          
C Loop over Non-zero orbitals ------------------------------------------ 
           DO IAT1 = 1, NNA
             IF( R2IJ(IAT1) .GT. RMAX2 ) EXIT

             IAVEC1   = JNA(IAT1)
             IS1      = ISA(IAVEC1)
             XVEC1(1) = -XIJ(1,IAT1)
             XVEC1(2) = -XIJ(2,IAT1)
             XVEC1(3) = -XIJ(3,IAT1)

             PHASE = K(1)*(XPO(1)+XIJ(1,IAT1))+
     .               K(2)*(XPO(2)+XIJ(2,IAT1))+
     .               K(3)*(XPO(3)+XIJ(3,IAT1))

             SI=SIN(PHASE)
             CO=COS(PHASE)
             EXPPHI=CMPLX(CO,SI,dp)
             
             DO IO = LASTO(IAVEC1-1) + 1, LASTO(IAVEC1)
               IPHI1 = IPHORB(IO)
               IUO   = INDXUO(IO)
               CALL PHIATM( IS1, IPHI1, XVEC1, PHIMU, GRPHIMU )
               !     Note implicit loop over spinor components
               CWAVE(:) = CWAVE(:) + PHIMU * WF(iuo,:) * EXPPHI 
            ENDDO       ! orbitals
         ENDDO ! atoms

           IF ( NSPIN .EQ. 1 ) THEN

              rwave = real(cwave(1), dp)
              iwave = aimag(cwave(1))

              call mod_and_phase(rwave,iwave,mwave,pwave)

      
           ELSEIF (NSPIN .EQ. 2) THEN

              if (ispin == 1) then
                 rwaveup = real(cwave(1), dp)
                 iwaveup = aimag(cwave(1))
                 call mod_and_phase(rwaveup,iwaveup,mwaveup,pwaveup)
              else
                 rwavedn = real(cwave(1), dp)
                 iwavedn = aimag(cwave(1))
                 call mod_and_phase(rwavedn,iwavedn,mwavedn,pwavedn)
              endif

           ELSE ! 4 or 8

              rwaveup = real(cwave(1), dp)
              iwaveup = aimag(cwave(1))
              rwavedn = real(cwave(2), dp)
              iwavedn = aimag(cwave(2))
              call mod_and_phase(rwaveup,iwaveup,mwaveup,pwaveup)
              call mod_and_phase(rwavedn,iwavedn,mwavedn,pwavedn)

           ENDIF

           IF (IDIMEN .EQ. 2) THEN
              IF ( NSPIN .EQ. 1 ) THEN
                 call write_plapo(unitre1,rwave)
                 call write_plapo(unitim1,iwave)
                 call write_plapo(unitmo1,mwave)
                 call write_plapo(unitph1,pwave)
             ELSEIF ( NSPIN .EQ. 2 ) THEN
                if (ispin == 1) then              
                 call write_plapo(unitre1,rwaveup)
                 call write_plapo(unitim1,iwaveup)
                 call write_plapo(unitmo1,mwaveup)
                 call write_plapo(unitph1,pwaveup)
                else
                 call write_plapo(unitre2,rwavedn)
                 call write_plapo(unitim2,iwavedn)
                 call write_plapo(unitmo2,mwavedn)
                 call write_plapo(unitph2,pwavedn)
                endif
             ELSE ! 4 or 8

               call write_plapo(unitre1,rwaveup)
               call write_plapo(unitim1,iwaveup)
               call write_plapo(unitmo1,mwaveup)
               call write_plapo(unitph1,pwaveup)
               
               call write_plapo(unitre2,rwavedn)
               call write_plapo(unitim2,iwavedn)
               call write_plapo(unitmo2,mwavedn)
               call write_plapo(unitph2,pwavedn)
               
             ENDIF
              
          ELSE IF (IDIMEN .EQ. 3) THEN
             IF (NSPIN .EQ. 1) THEN
                RWF(NPO,1) = RWAVE
                IMWF(NPO,1) = IWAVE
                MWF(NPO,1) = MWAVE
                PWF(NPO,1) = PWAVE
             ELSE IF (NSPIN .EQ.2) THEN
               if (ispin == 1) then
                  RWF(NPO,1) = RWAVEUP
                  IMWF(NPO,1) = IWAVEUP
                  MWF(NPO,1) = MWAVEUP
                  PWF(NPO,1) = PWAVEUP
               else
                  RWF(NPO,2) = RWAVEDN
                  IMWF(NPO,2) = IWAVEDN
                  MWF(NPO,2) = MWAVEDN
                  PWF(NPO,2) = PWAVEDN
               endif
             ELSE ! 4 or 8

               RWF(NPO,1) = RWAVEUP
               IMWF(NPO,1) = IWAVEUP
               MWF(NPO,1) = MWAVEUP
               PWF(NPO,1) = PWAVEUP

               RWF(NPO,2) = RWAVEDN
               IMWF(NPO,2) = IWAVEDN
               MWF(NPO,2) = MWAVEDN
               PWF(NPO,2) = PWAVEDN

            ENDIF
         ENDIF

C End x loop
 100    ENDDO  
C End y and z loops
 101    ENDDO  
 102    ENDDO  

        if (empty_domain) then
           call die("This domain is not covered by any atoms!")
        endif
        
        IF (IDIMEN .EQ. 3) THEN
! Write cube file data; Z dim changes fastest

           DO NX=1,NPX
              DO NY=1,NPY

                if ( NSPIN == 1 .or.
     &              (NSPIN == 2 .and. ispin == 1 ) .or.
     &              NSPIN >= 4 ) then
                  WRITE(UNITRE1,'(6e13.5)')
     $                (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                  WRITE(UNITIM1,'(6e13.5)')
     $                (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                  WRITE(UNITMO1,'(6e13.5)')
     $                (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                  WRITE(UNITPH1,'(6e13.5)')
     $                (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,1),NZ=1,NPZ)
                end if
                 
                if ( (NSPIN == 2 .and. ispin == 2) .or.
     &              NSPIN >= 4 ) then
                  WRITE(UNITRE2,'(6e13.5)')
     $                (RWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                  WRITE(UNITIM2,'(6e13.5)')
     $                (IMWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                  WRITE(UNITMO2,'(6e13.5)')
     $                (MWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                  WRITE(UNITPH2,'(6e13.5)')
     $                (PWF(NX+(NY-1)*NPX+(NZ-1)*NPX*NPY,2),NZ=1,NPZ)
                end if
              ENDDO
           ENDDO

          WRITE(6,'(A)')
          WRITE(6,'(A)')
     .      '   Your output files are:'
          IF (NSPIN .EQ. 1) THEN
            WRITE(6,'(A,A)') '   ',FN_RE
            WRITE(6,'(A,A)') '   ',FN_IM
            WRITE(6,'(A,A)') '   ',FN_MO
            WRITE(6,'(A,A)') '   ',FN_PH
          ELSE ! 2, 4 or 8
            WRITE(6,'(A,A)') '   ',FN_URE
            WRITE(6,'(A,A)') '   ',FN_UIM
            WRITE(6,'(A,A)') '   ',FN_DRE
            WRITE(6,'(A,A)') '   ',FN_DIM
            WRITE(6,'(A,A)') '   ',FN_UMO
            WRITE(6,'(A,A)') '   ',FN_UPH
            WRITE(6,'(A,A)') '   ',FN_DMO
            WRITE(6,'(A,A)') '   ',FN_DPH
            if (nspin >= 4) write(6,"(a)")
     $           "... up and down for spinor components"
          ENDIF
        ENDIF ! 3D


        CALL IO_CLOSE(UNITRE1)
        CALL IO_CLOSE(UNITIM1)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITRE2)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITIM2)
        CALL IO_CLOSE(UNITMO1)
        CALL IO_CLOSE(UNITPH1)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITMO2)
        IF (NSPIN >= 2) CALL IO_CLOSE(UNITPH2)
     
          
      ENDDO  wf_loop ! wfn
      ENDDO  spin_loop ! spin_block
      ENDDO  k_loop ! k-point

      DEALLOCATE(RWF)
      DEALLOCATE(IMWF)
      DEALLOCATE(MWF)
      DEALLOCATE(PWF)

      CONTAINS
      
        subroutine write_plapo(lun,w)
        integer, intent(in) :: lun
        real(dp), intent(in) :: w
        ! rest by host association
        
        WRITE(lun,'(3F12.5)') PLAPO(NPO,1),PLAPO(NPO,2),w
        IF ( MOD(NPO,NPX) .EQ. 0 ) WRITE(lun,*)
        end subroutine write_plapo
      
       subroutine write_cube_header(lun,fname)
       integer, intent(in) :: lun
       character(len=*), intent(in) :: fname
       ! rest by host association
      
            WRITE(LUN,*) FNAME
            WRITE(LUN,*) FNAME
            WRITE(LUN,'(i5,4f12.6)') NAINCELL, XMIN, YMIN, ZMIN
            WRITE(LUN,'(i5,4f12.6)') NPX,(OCELL(1,J)/(NPX-1),J=1,3)
            WRITE(LUN,'(i5,4f12.6)') NPY,(OCELL(2,J)/(NPY-1),J=1,3)
            WRITE(LUN,'(i5,4f12.6)') NPZ,(OCELL(3,J)/(NPZ-1),J=1,3)
            DO IA = 1,NAINCELL
              WRITE(LUN,'(i5,4f12.6)') IZA(IA),0.0,
     .                                   (XAINCELL(IX,IA),IX=1,3)
            ENDDO
       end subroutine write_cube_header

      END

      subroutine mod_and_phase(rw,iw,mw,pw)
      integer, parameter :: dp = selected_real_kind(10,100)
      real(dp), intent(in) :: rw
      real(dp), intent(inout) :: iw
      real(dp), intent(out) :: mw, pw

      real(dp), parameter :: pi = 3.1141592653589_dp
      
      MW = SQRT(RW**2 + IW**2)
      IF (ABS(IW) .LT. 1.D-6) IW = ABS(IW)
      IF (ABS(RW) .LT. 1.D-12) THEN
         IF (IW .GT. 0.0D0) PW = PI/2.0D0
         IF (IW .LT. 0.0D0) PW = -PI/2.0D0
      ELSE
         PW = ATAN(IW/RW)
      ENDIF
      IF (RW .LT. 0.0D0 .AND. IW .GE. 0.0d0)
     .     PW = PW + PI
      IF (RW .LT. 0.0D0 .AND. IW .LT. 0.0D0)
     .     PW = PW - PI
      end
