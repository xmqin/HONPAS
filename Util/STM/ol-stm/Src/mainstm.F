! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      PROGRAM MAINSTM

! Reads eigenfunctions and eigenvalues from SIESTA and computes the
!     LDOS needed to simulate STM images or STS curves.
!     Input:
!     - WFSX file
!     - .DIM and .PLD files      
!     - .ion files
!     - .VH file
!
!     Output for STM:
!     - A 'gridfunc' file, which for STM simulations will have information
!     about all the planes desired. (It marks the contents as "non-periodic"
! along Z.
!     Output for STS:
!     - A 'gridfunc' file with a single plane perpendicular to Z, and using
!     the 'spin' dimension to store values for different energies.
!     - An auxiliary file with the values of the energy used in sampling.
      
! Coded by P. Ordejon, November 2004
! Modified by N. Lorente, August 2005
! Re-structured, adapted to NC/SOC, and extended to STS: Alberto Garcia, March 2019


      USE PRECISION, only:dp
      USE BASIS_IO, only: read_basis_ascii
      USE LISTSC_MODULE, ONLY: LISTSC_INIT
      USE FDF
      use broadening, only: width, gaussian_broadener

      IMPLICIT NONE

      INTEGER ns,
     .   NO_U, NO_S, NA_S, NSPIN, MAXNA,
     .   NSC(3), NK, no_u_in_wfs_file, idummy

      INTEGER
     .  NPX, NPY, NPZ, IUNITCD, IX, IV

      INTEGER, DIMENSION(:), ALLOCATABLE ::
     .  ISA, LASTO, IPHORB, INDXUO

      real(dp), DIMENSION(:), ALLOCATABLE  ::  WK
      INTEGER, DIMENSION(:,:), ALLOCATABLE ::  INDW

      REAL(DP)
     .   CELL(3,3), VOLUME, VOLCEL, RMAXO, UCELL(3,3), V0, EMAX, EMIN

      REAL(DP)  ZMIN, ZMAX, ZREF, ARMUNI

      REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE ::  RPSI,IPSI
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE ::   E
      REAL(DP), DIMENSION(:,:), ALLOCATABLE ::   XA, K

      CHARACTER  FILEIN*20, FILEOUT*20, SNAME*40

      INTEGER :: wf_unit
      logical :: gamma, non_coll
      integer :: nspin_blocks

      integer :: ne
      logical :: want_sts
      EXTERNAL :: REDATA, REINITSTM, VOLCEL

C ****** READ FROM SIESTA **********************************************
! integer ns                  : Total number of species
C INTEGER NO_U                : Total number of orbitals in the unit cell
C INTEGER NO_S                : Total number of orbitals in the supercell
C INTEGER NA_S                : Total number of atoms in the supercell
C INTEGER NSPIN               : Number of different spin polarizations (1, 2, 4 [NC/SOC])
C INTEGER MAXNA               : Maximum number of neighbours of any atom   !** Ambiguous: range needed
C INTEGER NSC(3)              : Num. of unit cells in each supercell direction
C INTEGER ISA(MAXA)           : Species index of each atom in the supercell
C INTEGER LASTO(0:MAXA)       : Last orbital of each atom in array iphorb
C INTEGER IPHORB(MAXO)        : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)        : Equivalent orbital in unit cell
C INTEGER NK                  : Number of k-points in wave functions file
C REAL*8  CELL(3,3)           : Supercell vectors CELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  UCELL(3,3)          : Unit cell vectors UCELL(IXYZ,IVECT)
C                               (units in bohrs)
C REAL*8  VOLUME              : Volumen of unit cell (in bohr**3)
C REAL*8  RMAXO               : Maximum range of basis orbitals
C REAL*8  XA(3,NA_S)          : Atomic coordinates in cartesian coordinates
C                               (units in bohrs)

C ****** INFORMATION FOR 3D GRID ***************************
C INTEGER NPX, NPY, NPZ       : Number of points generated along x and y
C                               (and z for 3D grids) 
C INTEGER IUNITCD             : Units of the charge density
C REAL*8  ZREF                : Position of reference plane for wf. extrapol.
C REAL*8  ZMIN, ZMAX          : Limits of the 3D grid for z-axis
C REAL*8  ARMUNI              : Conversion factors for the charge density


      FILEIN  = 'stdin'
      FILEOUT = 'out.fdf'
      CALL FDF_INIT(FILEIN,FILEOUT)

      WRITE(6,*)
      WRITE(6,*) 'STM/STS Simulation program' //
     $           ' with optional wfs projection'
      WRITE(6,*) 'P. Ordejon and N. Lorente, Nov. 04'
      WRITE(6,*) 'A. Garcia, March 2019'
      WRITE(6,*)
      WRITE(6,*) 'Reading information from SIESTA'
      WRITE(6,*)


C Read some variables from SIESTA to define the limits of some arrays --
      CALL REINITSTM( NO_S, NA_S, NO_U, IDUMMY, MAXNA, NSPIN)

      ALLOCATE(XA(3,NA_S))
      ALLOCATE(LASTO(0:NA_S))
      ALLOCATE(ISA(NA_S))
      ALLOCATE(IPHORB(NO_S))
      ALLOCATE(INDXUO(NO_S))

C Read some variables from SIESTA --------------------------------------
      CALL REDATA( NO_S, NA_S, NO_U, IDUMMY, NSPIN,
     .             ISA, IPHORB, INDXUO, LASTO,
     .             CELL, NSC, XA, RMAXO )

C Read the information about the basis set -----------------------------
      CALL READ_BASIS_ASCII(ns)

C Initialize listsc ----------------------------------------------------
      CALL LISTSC_INIT( NSC, NO_U )

C Calculate the volume of the unit cell --------------------------------
      VOLUME = VOLCEL( CELL ) / (NSC(1) * NSC(2) * NSC(3))

C Calculate unit cell vectors
      DO IX=1,3
        DO IV=1,3
          UCELL(IX,IV) = CELL(IX,IV)/NSC(IV)
        ENDDO
      ENDDO

      ! Read heading info from WFSX file

      SNAME = FDF_STRING('SystemLabel','siesta')

      CALL IO_ASSIGN(wf_unit)
      OPEN (wf_unit, FILE=trim(SNAME)//'.WFSX', FORM='unformatted',
     $       STATUS='unknown',position='rewind')

      read(wf_unit) nk, gamma
      read(wf_unit) nspin

      ! Non-collinear or SOC files have a single "spin" block as opposed to collinear-spin
      ! files, which contain two spin blocks per k section.
      non_coll = (nspin >= 4)
      nspin_blocks = nspin
      if (non_coll) nspin_blocks = 1

      read(wf_unit) no_u_in_wfs_file
      if (no_u /= no_u_in_wfs_file) then
         call die("Mismatch in no_u in WFSX and DIM/PLD files")
      endif
      read(wf_unit)   ! orbital labels

      ! The file is now positioned to start reading wfs in routines 'stm' or 'sts'

      CALL READSTM( VOLUME, 
     .     IUNITCD, NPX, NPY, NPZ, ZREF, ZMIN, ZMAX, EMAX, EMIN,
     .     ARMUNI )

! Call routine to calculate the value of the potential at vacuum

      CALL VACPOT(V0)

! STS calculations can be done on several planes perpendicular to z
! at the same time. But note that the resulting data file might be large.
! It can be argued that a single point, or line, might be enough. But      
! we can easily reuse the code for 'stm' in the 'sts' routine. In 
! postprocessing the extra information can be filtered.
! Besides, for the "projection" technique we need a full XY-unit cell,
! since we need to fourier-transform the wavefunctions.      

      ne = fdf_get('STS.NumberOfPoints',1)
      want_sts = (ne > 1)
      
      if (want_sts) then

           ! Later: use function pointers for broadener
           width = fdf_get('STS.broadening',0.2_dp)  ! in eV
          
           CALL STS ( NA_S, NO_S, NO_U, MAXNA, nspin,
     $               nspin_blocks, non_coll,
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .               wf_unit, NK, gamma,
     .               ZREF, ZMIN, ZMAX, NPX, NPY, NPZ,
     .               V0, EMAX, EMIN, NE, gaussian_broadener,
     .               ARMUNI, IUNITCD, RMAXO )
        else

           CALL STM ( NA_S, NO_S, NO_U, MAXNA, nspin,
     $               nspin_blocks, non_coll,
     .               ISA, IPHORB, INDXUO, LASTO, XA, CELL, UCELL,
     .               wf_unit, NK, gamma,
     .               ZREF, ZMIN, ZMAX, NPX, NPY, NPZ,
     .               V0, EMAX, EMIN,
     .               ARMUNI, IUNITCD, RMAXO )

        endif
        
        call io_close(wf_unit)

      END PROGRAM MAINSTM
