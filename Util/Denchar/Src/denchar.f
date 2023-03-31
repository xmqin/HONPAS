! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      PROGRAM DENCHAR

C **********************************************************************
C Reads density matrix from SIESTA and calculates the charge density
C at the points of a plane in real space, or at a 3D grid of points
C Coded by J. Junquera 11/98
C Modified by J. Junquera 07/01 
C 3D and wavefunction capabilities coded by P. Ordejon, June 2003
C Modified to handle complex wavefunctions and multiple k-points
C by P. Ordejon, July 2004
C Modified to use the more efficient WFSX format for wave functions
C by A. Garcia, May 2012
C Added more files in 3D mode (A. Garcia, Nov 2016)
! Added support for NC/SOC wfs and streamlined by A. Garcia (March-June 2019)
C **********************************************************************
C
C  Modules
C
      USE PRECISION
      USE BASIS_IO
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
      use parallel, only: nodes, node
      use m_getopts

      IMPLICIT NONE

      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels

      INTEGER
     .   NO_U, NO_S, NA_S, NSPIN, MAXND, MAXNA,
     .   NSC(3), NK
      integer :: ns_dummy

      INTEGER
     .  IDIMEN, IOPTION, NPX, NPY, NPZ, IUNITCD, ISCALE

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO, 
     .  NUMD, LISTD, LISTDPTR

      real(dp)
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO

      real(dp)
     .  XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .  COORPO(3,3), NORMAL(3), DIRVER1(3),
     .  DIRVER2(3), ARMUNI, EMIN

      real(dp), DIMENSION(:,:), ALLOCATABLE     :: XA, DSCF, K
      real(dp), DIMENSION(:), ALLOCATABLE       :: DATM

      CHARACTER
     .  FILEIN*20, FILEOUT*20, sname*30

      LOGICAL 
     .  FOUND, CHARGE, WAVES, ionode, debug
      logical :: gamma_wfsx, non_coll
      integer :: nspin_wfsx, nspin_blocks, wf_unit
      integer :: no_u_wfsx
      integer :: kpoint_selected, wf_selected
      
      EXTERNAL
     .  IODM, READPLA, REDATA_DENCHAR, REINIT, RHOOFR, VOLCEL

      DATA NORMAL /0.D0,0.D0,1.D0/
      DATA COORPO /1.D0,0.D0,0.D0,0.D0,1.D0,0.D0,0.D0,0.D0,1.D0/
      DATA DIRVER1 /1.D0,0.D0,0.D0/
      DATA DIRVER2 /0.D0,1.D0,0.D0/

C ****** READ FROM SIESTA **********************************************
C INTEGER NO_U                : Total number of orbitals in the unit cell
C INTEGER NO_S                : Total number of orbitals in the supercell
C INTEGER NA_S                : Total number of atoms in the supercell
C INTEGER NSPIN               : Number of different spin polarizations
C                               Nspin = 1 => Unpolarized, Nspin = 2 => Polarized
C INTEGER MAXND               : Maximum number
C                               of basis orbitals interacting, either directly
C                               or through a KB projector, with any orbital
C INTEGER MAXNA               : Maximum number of neighbours of any atom
C INTEGER NSC(3)              : Num. of unit cells in each supercell direction
C INTEGER NUMWF               : Max num of wavefncts to print for a given k-po.
C INTEGER NWF(NK)             : Num of wavefncts to print for each k-point
C INTEGER ISA(MAXA)           : Species index of each atom in the supercell
C INTEGER LASTO(0:MAXA)       : Last orbital of each atom in array iphorb
C INTEGER IPHORB(MAXO)        : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)        : Equivalent orbital in unit cell
C INTEGER NUMD(NO_U)          : Number of nonzero elements of each row of the
C                               Hamiltonian matrix between atomic orbitals
C INTEGER LISTD(MAXND)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C INTEGER LISTDPTR(NO_U)      : Pointer to where each row of listh starts - 1
C                               The reason for pointing to the element before
C                               the first one is so that when looping over the
C                               elements of a row there is no need to shift by
C                               minus one.
C INTEGER NK                  : Number of k-points in wave functions file
C INTEGER INDW(NK,NUMWF)      : Indices of wavefuncs to print for each k-point
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)
C REAL*8  DATM(NO_S)          : Occupations of basis orbitals in free atom
C REAL*8  DSCF(MAXND,NSPIN)   : Density Matrix (DM)
C ****** INFORMATION OF THE POINT, PLANE OR 3D GRID ***********************
C INTEGER IDIMEN              : Specifies 2D or 3D mode
C LOGICAL CHARGE              : Should charge density be computed?
C LOGICAL WAVES               : Should wave functions be computed?
C INTEGER IOPTION             : Option to generate the plane or 3D grid
C                               1 = Normal vector to xy plane (ie, z direction)
C                               2 = Two vectors belonging to the xy plane
C                               3 = Three points of the xy plane
C                               4 = Three atomic indices define the xy plane
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) directions in a system of 
C                               reference in which the third component of the 
C                               points of the plane is zero 
C                               (Plane Reference Frame; PRF)
C INTEGER IUNITCD             : Units of the charge density
C INTEGER ISCALE              : Units of the points of the plane or 3D grid
C REAL*8  XMIN, XMAX          : Limits of the plane in the PRF for x-axis
C REAL*8  YMIN, YMAX          : Limits of the plane in the PRF for y-axis
C REAL*8  ZMIN, ZMAX          : Limits of the or 3D grid in the PRF for z-axis
C REAL*8  COORPO(3,3)         : Coordinates of the three points used 
C                               to define the xy plane
C REAL*8  NORMAL(3)           : Components of the normal vector 
C                               used to define the xy plane (z-direction)
C REAL*8  DIRVER1(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8  DIRVER2(3)          : Components of the first vector contained 
C                               in the plane
C REAL*8  ARMUNI              : Conversion factors for the charge density
C ****** INTERNAL VARIABLES ********************************************
C LOGICAL FOUND               : Has DM been found in disk?
C                               (Only when task = 'read')
C **********************************************************************
C
!     Process options
!
      kpoint_selected = 0
      wf_selected = 0
      n_opts = 0
      do
         call getopts('dhk:w:',
     $        opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
         case ('d')
            debug = .true.
         case ('k')
            read(opt_arg,*) kpoint_selected
         case ('w')
            read(opt_arg,*) wf_selected
         case ('h')
            write(0,*) " See manual "
            STOP
         case ('?',':')
            write(0,*) "Invalid option: ", opt_arg(1:1)
            write(0,*) "Use -h option for manual"
            write(0,*) ""
            !call manual()
            STOP
         end select
      enddo

      nargs = command_argument_count()
      nlabels = nargs - n_opts + 1
      if (nlabels == 1)  then
        call get_command_argument(n_opts,value=FILEIN,status=iostat)
        if ( iostat /= 0 ) then
           stop "Cannot get name of fdf file in command line"
        end if
      else
        FILEIN = 'stdin'
      endif

      nodes = 1
      node = 0
      ionode = .true.

C Set up fdf -----------------------------------------------------------
      !      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINIT( NO_S, NA_S, NO_U, MAXND, MAXNA, NSPIN, IDIMEN,
     .     CHARGE, WAVES )

C Allocate some variables ----------------------------------------------
      ALLOCATE(XA(3,NA_S))
      ALLOCATE(LASTO(0:NA_S))
      ALLOCATE(ISA(NA_S))
      ALLOCATE(IPHORB(NO_S))
      ALLOCATE(INDXUO(NO_S))
      ALLOCATE(DATM(NO_S))

C Read some variables from SIESTA --------------------------------------
      CALL REDATA_DENCHAR( NO_S, NA_S, NO_U, MAXND, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO, DATM )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns_dummy)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))

C If this is a charge calculation, allocate space for DM
      IF (CHARGE) THEN
        ALLOCATE(LISTDPTR(NO_U))
        LISTDPTR(:) = 0

        ALLOCATE(NUMD(NO_U))
        NUMD(:) = 0

C Allocate some other variables ----------------------------------------
        IF (.NOT.ALLOCATED(LISTD)) THEN
          ALLOCATE(LISTD(MAXND))
        ENDIF
        
        IF (ALLOCATED(DSCF)) THEN
          DEALLOCATE(DSCF)
        ENDIF
        ALLOCATE(DSCF(MAXND,NSPIN))

C Read Density Matrix from files ---------------------------------------
        CALL IODM('READ', MAXND, NO_U, NSPIN,
     .            NUMD, LISTDPTR, LISTD, DSCF, FOUND )
        IF (.NOT. FOUND) THEN
          WRITE(6,*)' DENSITY MATRIX NOT FOUND              '
          WRITE(6,*)' CHECK YOU HAVE COPIED IT FROM THE       '
          WRITE(6,*)' DIRECTORY WHERE YOU HAVE RUN SIESTA   '
          STOP
        ENDIF 
      ENDIF

      IF (WAVES) THEN

      ! Read header of WFSX file

      SNAME = FDF_STRING('SystemLabel','siesta')
      CALL IO_ASSIGN(wf_unit)
      OPEN (wf_unit, FILE=trim(SNAME)//'.WFSX', FORM='unformatted',
     $       STATUS='unknown',position='rewind')

      read(wf_unit) nk, gamma_wfsx
      read(wf_unit) nspin_wfsx

      ! Non-collinear or SOC files have a single "spin" block as opposed to collinear-spin
      ! files, which contain two spin blocks per k section.
      non_coll = (nspin_wfsx >= 4)
      nspin_blocks = nspin_wfsx
      if (non_coll) nspin_blocks = 1

      read(wf_unit) no_u_wfsx
      if (no_u /= no_u_wfsx) then
         call die("Mismatch in no_u in WFSX and DIM/PLD files")
      endif
      read(wf_unit)   ! orbital labels

      ENDIF

      IF (CHARGE .OR. WAVES) THEN
C Read option to generate the plane or 3D-grid -------------------------
        CALL READPLA( NA_S, XA, VOLUME, IDIMEN,
     .                IOPTION, IUNITCD, ISCALE, NPX, NPY, NPZ,
     .                XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     .                COORPO, NORMAL, DIRVER1, DIRVER2, 
     .                ARMUNI )
      ENDIF

C Form Density Matrix for Neutral and Isolated Atoms -------------------
      IF (CHARGE) THEN
      
        CALL RHOOFR( NA_S, NO_S, NO_U, MAXND, MAXNA, NSPIN, 
     .               ISA, IPHORB, INDXUO, LASTO,
     .               XA, CELL, NUMD, LISTD, LISTDPTR, DSCF, DATM,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO )
      ENDIF

      IF (WAVES) THEN
        CALL WAVOFR( NA_S, NO_S, NO_U, MAXNA, NSPIN_wfsx, nspin_blocks,
     $               non_coll,
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL,
     .               wf_unit, NK, gamma_wfsx,
     .               IDIMEN, IOPTION, XMIN, XMAX, YMIN, YMAX, 
     .               ZMIN, ZMAX, NPX, NPY, NPZ, COORPO, NORMAL, 
     .               DIRVER1, DIRVER2, 
     .               ARMUNI, IUNITCD, ISCALE, RMAXO,
     .               kpoint_selected, wf_selected)
      ENDIF

      END PROGRAM DENCHAR
