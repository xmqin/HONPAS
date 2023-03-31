! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE RANGER( MODE, NX, CELL, RANGE, NA, XA,
     .                   NAMOVE, IAMOVE,
     .                   IA0, ISC, X0,
     .                   NNA, JAN, XIJ, R2IJ )

C ********************************************************************
C Finds the neighbours of an atom in a cell with periodic boundary 
C conditions. Alternatively, it finds the atoms within a sphere 
C centered at an arbitrary point. It also allows to update the atomic
C positions one at a time, what is useful in Montecarlo simulations.
C Written by J.M.Soler. Nov'96.
C *********** INPUT **************************************************
C CHARACTER*4 MODE       : MODE='CELL' => Initialize or reshape cell
C                          MODE='MOVE' => Move atom(s)
C                          MODE='FIND' => Find neighbours
C INTEGER NX             : Space dimension
C REAL*8  CELL(NX,NX)    : Unit cell vectors CELL(IXYZ,IVECT)
C REAL*8  RANGE          : Maximum distance of neighbours required
C INTEGER NA             : Number of atoms
C REAL*8  XA(NX,NA)      : Atomic positions in cartesian coordinates
C INTEGER NAMOVE         : Number of atoms to be moved
C                          (not used unless MODE='MOVE')
C INTEGER IAMOVE(NAMOVE) : Index(es) of atom(s) moved (not used
C                          unless MODE='MOVE' and 0<NAMOVE<NA)
C INTEGER IA0            : Atom whose neighbours are needed.
C                          If IA0=0, point X0 is used as origin instead
C                          (not used unless MODE='FIND')
C INTEGER ISC            : Single-counting switch (0=No, 1=Yes).
C                          If ISC=1, only neighbours with JA.LE.IA0
C                          are included in JAN
C                          (not used unless MODE='FIND' and IA0.NE.0)
C REAL*8  X0(NX)         : Origin from which atoms are to be found,
C                          in cartesian coordinates.
C                          (not used unless MODE='FIND' and IA0=0)
C INTEGER NNA            : Size of arrays JAN, XIJ and R2IJ
C                          (not used unless MODE='FIND')
C *********** OUTPUT *************************************************
C REAL*8  CELL(NX,NX)  : Unit cell vectors CELL(IXYZ,IVECT)
C                        The output cell is generated only when input
C                        CELL has zero volume
C INTEGER NNA          : Number of 'neighbour' atoms within RANGE of
C                        atom IA0 or position X0 (only for MODE='FIND')
C INTEGER JAN(NNAin)   : Atom-index of neighbours (only for MODE='FIND')
C REAL*8  XIJ(NX,NNAin): Vectors from atom IA0 or point X0 to neighbours
C                        in cartesian coordinates (only for MODE='FIND')
C REAL*8  R2IJ(NNAin)  : Squared distances to neighbours
C                        (only for MODE='FIND')
C *********** UNITS **************************************************
C Units of CELL, RANGE, XA and X0 are arbitrary but must be equal.
C All vectors in cartesian coordinates.
C *********** SUBROUTINES USED ***************************************
C DISMIN, DOT, RECCEL, VOLCEL
C *********** BEHAVIOUR **********************************************
C This is a 'remembering' routine, that saves a single copy of required 
C   information on the system. Therefore, it cannot be used
C   simultaneously for different cells or sets of atoms.
C A call with MODE='CELL' is required if the cell is changed. This also
C   updates all atomic positions with no need of a MODE='MOVE' call.
C A call with MODE='MOVE' is required if any atoms are moved without
C   reshaping the cell, before any subsequent 'FIND' calls.
C The routine knows if it has not been ever called, so that initial 
C   calls with MODE='CELL' and MODE='MOVE' are implicit but not required
C If MODE='MOVE' and NAMOVE=NA, all the atomic positions are
C   reinitialized, and the list IAMOVE is not used. This is also
C   true if MODE='CELL', irrespective of the value of NAMOVE
C This routine works always with periodic boundary conditions.
C   If periodic boundary conditions are not desired, you must either
C   define a CELL large enough to contain all the atoms without
C   'interaction' (i.e. distances smaller than RANGE) between 
C   different cells, or make CELL vectors identical to zero, in which
C   case an appropiate CELL is generated automatically by RANGER.
C   The size of this cell is determined by parameters DXMARG and DXRANG
C   below, and should be safe for small atomic displacements, but notice
C   that, if atoms move too much after the cell is generated (i.e. after
C   the last MODE='CELL' call), spureous 'interactions' may occur.
C   Do not make CELL extremely large to avoid this, since internal
C   array memory increases with cell volume.
C If the number of neighbour atoms found is larger than the size of
C   the arrays JAN, XIJ and R2IJ, i.e. if NNAout > NNAin, these arrays
C   are filled only up to their size NNAin. With dynamic memory
C   allocation, this allows to find first the required array sizes
C   and then find the neighbours. Notice however that no warning is
C   given, so that you should always check that NNAout.LE.NNAin.
C Different ranges can be used for different atoms or origins, but for
C   good performance, the largest range should be used in the initial
C   call (with MODE='CELL').
C There are no limitations regarding cell shape or size. The range may
C   be larger than the cell size, in which case many 'images' of the
C   same atom will be included in the neighbour list, with different
C   interatomic vectors XIJ and distances R2IJ.
C The atom IA0 itself is included in the neighbour list, with zero
C   distance. You have to discard it if you want so.
C CPU time and memory scale linearly with the number of atoms, for
C   sufficiently large numbers.
C This is not an extremely optimized routine. The emphasis has been
C   rather put on functionality. It is not intended to vectorize or
C   parallelize well either. However, it is expected to be very
C   reasonably efficient on scalar machines, since most of the time
C   will be spent on the relatively simple and optimized 'search 
C   section' at the end.
C Code from beginning of loop 120 to label 130 was designed to exclude
C   neighbour cells which are inside a parallelepiped containing the
C   range sphere, but outside this sphere. This takes a nonneglegible
C   CPU time compared to one 'FIND' call for each atom (i.e. a MD step),
C   but may be worth if CELL is fixed or changes rarely compared to the
C   atomic displacements. But this is only true when CELL is strongly
C   nonorthorhombic (like fcc, bcc or hex) or when parameter NCR.GE.3
C   If you have doubts, don't worry: all this is normally irrelevant
C   unless you do Parrinello-Rahman dynamics.
C *********** USAGE **************************************************
C Example of usage for a molecular dynamics simulation:
C    DIMENSION CELL(3,3),JAN(MAXNNA),R2IJ(MAXNNA),XIJ(3,MAXNNA),XA(3,NA)
C    Define CELL and initial positions XA
C    DO ITER = 1,NITER                (Molecular dynamics iteration)
C      CALL RANGER('MOVE',3,CELL,RANGE,NA,XA,NA,0,0,1,X0,
C                  NNA,JAN,XIJ,R2IJ)
C      Initialize to zero all atomic forces FA(IX,IA)
C      DO IA = 1,NA                   (Loop on atoms)
C        NNA = MAXNNA
C        CALL RANGER('FIND',3,CELL,RANGE,NA,XA,0,0,IA,1,X0,
C                    NNA,JAN,XIJ,R2IJ)
C        IF (NNA.GT.MAXNNA) STOP 'Parameter MAXNNA too small'
C        DO IN = 1,NNA                (Loop on neighbours of IA)
C          JA = JAN(IN)               (Atomic index of neighbour)
C          RIJ = SQRT(R2IJ(IN))       (Interatomic distance)
C          IF (RIJ .GT. 1.D-12) THEN  (Discard atom itself)
C            Find interatomic force FIJ( RIJ )
C            DO IX = 1,3
C              FA(IX,IA) = FA(IX,IA) - FIJ * XIJ(IX,IN) / RIJ
C              FA(IX,JA) = FA(IX,JA) + FIJ * XIJ(IX,IN) / RIJ
C            ENDDO
C          ENDIF
C        ENDDO
C      ENDDO
C      Move atomic positions XA       (Molecular dynamics step)
C    ENDDO
C
C Example of usage for a Montecarlo simulation:
C    DIMENSION CELL(3,3),JAN(MAXNNA),R2IJ(MAXNNA),XIJ(3,MAXNNA),
C              XA(3,NA),XNEW(3)
C    Define CELL and initial positions XA
C    CALL RANGER('CELL',3,CELL,RANGE,NA,XA,0,0,0,0,XNEW,
C                 NNA,JAN,XIJ,R2IJ)
C    Find here each atom's interaction energy EA(IA)
C    DO ITER = 1,NITER                (Montecarlo iteration)
C      Choose atom moved IA and its trial new position XNEW
C      NNA = MAXNNA
C      CALL RANGER('FIND',3,CELL,RANGE,NA,XA,0,0,0,0,XNEW,
C                   NNA,JAN,XIJ,R2IJ)
C      IF (NNA.GT.MAXNNA) STOP 'Parameter MAXNNA too small'
C      EANEW = 0.D0
C      DO IN = 1,NNA                  (Loop on neighbours of IA)
C        JA = JAN(IN)                 (Atomic index of neighbour)
C        IF (JA.NE.IA)                (Discard atom itself)
C          RIJ = SQRT(R2IJ(IN))       (Interatomic distance)
C          EANEW = EANEW + VIJ(RIJ)   (Add interaction energy with JA)
C        ENDIF
C      ENDDO
C      IF (EXP(-(EANEW-EA(IA))/TEMP).GT.RAND()) THEN   (Move atom)
C        EA(IA) = EANEW
C        DO IX = 1,3
C          XA(IX,IA) = XNEW(IX)
C        ENDDO
C        CALL RANGER('MOVE',3,CELL,RANGE,NA,XA,1,IA,0,0,XNEW,
C                    NNA,JAN,XIJ,R2IJ)
C      ENDIF
C    ENDDO
C *********** ALGORITHM **********************************************
C The unit cell is divided in 'mesh-cells', and a list of atoms in
C   each cell is stored. To look for the neighbours of an atom, the
C   distances to the atoms in the neighbour cells are calculated and
C   compared to the range.
C The list of atoms within a cell is stored as an 'ordered list', of
C   the kind so popular in C, using pointers from an atom to the next.
C To deal with periodic boundary conditions, the mesh is 'extended'
C   on each side of the unit cell, and an index from the extended to
C   the unextended meshes is stored. In the extended mesh, the 
C   'index-shifts' and the vector distances from a mesh point (within
C   the unit cell) to its neighbour mesh points are independent of 
C   the point, and they can thus be stored only once.
C To calculate the vector from an atom to its neighbours, the 
C   position of each atom relative to its mesh cell is also stored.
C   Thus, the vector between two atoms is the vector between their
C   cells plus the difference between their positions relative to
C   their cells.
C *********** LANGUAGE ***********************************************
C Illegal ANSI Fortran77 features used: IMPLICIT NONE, INCLUDE
C Illegal ANSI Fortran90 features used: none.
C *********** HISTORY AND CHANGES ************************************
C This is an improved version of routine NEIGHB, first written in C by
C   J.M.Soler in March 1995 and later translated to Fortran in Nov'96
C The added capabilities of RANGER over NEIGHB are:
C   - Using a variable space dimension
C   - Finding the atoms within a sphere centered at an arbitrary point
C   - Updating the position of a single atom
C   - Finding the number of neighbours without any arrays
C   - Automatic cell generation for nonperiodic boundary conditions
C All the algorithms and coding developed by J.M.Soler. November 1996.
C ********************************************************************

C
C  Modules
C
      use precision

      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER         MODE*4
      INTEGER           NA, NAMOVE, NNA, NX
      INTEGER           IA0, IAMOVE(*), ISC, JAN(*)
      real(dp)  CELL(*), RANGE, R2IJ(*),
     .                  X0(NX), XA(NX,NA), XIJ(NX,*)

C NCR is the ratio between range radius and mesh-planes distance.
C It fixes the size (and number) of mesh cells.
C Recommended values are between 1 and 3
      INTEGER NCR
      PARAMETER ( NCR = 2 )

C DXMARG and DXRANG are used for automatic CELL generation
C DXMARG is the minimum margin relative to coordinate range
C DXRANG is the minimum margin relative to RANGE
C EPS is a small number to be subtracted from 1
      real(dp) DXMARG, DXRANG, EPS 
      PARAMETER ( DXMARG = 0.1D0  )
      PARAMETER ( DXRANG = 1.0D0  )
      PARAMETER ( EPS    = 1.D-14 )

C Variable-naming hints:
C   Character(s) indicates
C          A     Atom
C          D     Difference
C          E     Extended (mesh)
C          I,J   Index
C          M     Mesh cell or point (lower vertex of mesh cell)
C          N     Neighbour or Number (if first character)
C          R     distance or vector modulus
C          V     Vertex
C          X     coordinate or basis vector
C          1     lower bound
C          2     upper bound or square

C Internal functions, variables and arrays
C REAL*8  CELMSH(MX*MX) Mesh-cell vectors
C REAL*8  DMX(MX)       In-cell atomic position in mesh coordinates 
C REAL*8  DOT()         Finds the scalar product of two vectors
C REAL*8  DPLANE        Distance between lattice or mesh planes
C REAL*8  DRM           Minimum distance between two mesh cells
C REAL*8  DX(3)         Vector between two atoms
C REAL*8  DX0M(MX)      Origin position within mesh cell
C REAL*8  DXAM(MX)      Atom position within mesh cell
C REAL*8  DXM(MX)       Minimum vector between two mesh cells
C REAL*8  DXMARG        Parameter defined above
C REAL*8  DXNM(MX,MAXNM) Cartesian vector between neighbour mesh points
C REAL*8  DXRANG        Parameter defined above
C REAL*8  EPS           Parameter defined above
C LOGICAL FRSTME        First call to subroutine? 
C INTEGER IA            Atom index
C INTEGER IAM           Atom-to-move index
C INTEGER IA1M(NM)      Pointer to first atom in mesh cell
C INTEGER IANEXT(NA)    Pointer to next atom in mesh cell
C INTEGER IAPREV(NA)    Pointer to previous atom in mesh cell
C INTEGER IDNM(MAXNM)   Index-distance between neighbour mesh points 
C INTEGER IM            Mesh index  
C INTEGER IEM           Extended-mesh index  
C INTEGER IEMA(NA)      Extended-mesh index of atoms
C INTEGER I1EMX(MX)     Minimum vaue of extended-mesh-coordinate indices
C INTEGER I2EMX(MX)     Maximum vaue of extended-mesh-coordinate indices
C INTEGER I1MX(MX)      Minimum vaue of mesh-coordinate indices
C INTEGER I2MX(MX)      Maximum vaue of mesh-coordinate indices
C INTEGER IMX(MX)       Mesh-cell index for each mesh vector
C INTEGER IN            Neighbour-mesh-cell index
C INTEGER I1NX(MX)      Minimum neighbour-cell-coordinate indices
C INTEGER I2NX(MX)      Maximum neighbour-cell-coordinate indices
C INTEGER INX(MX)       Neighbour-cell-coordinate indices
C LOGICAL INSIDE        Are two mesh cells within each other's range?
C INTEGER IMESH(NEM)    Correspondence between extended and normal mesh
C INTEGER IV            Vertex index
C INTEGER IVX           Vertex coordinate index
C INTEGER IX            Cartesian coordinte index
C INTEGER IXX           Double cartesian coordinte index
C INTEGER JA            Atom index
C INTEGER JEM           Extended-mesh index
C INTEGER J1NX(MX)      Minimum vertex-coordinate indices
C INTEGER J2NX(MX)      Maximum vertex-coordinate indices
C INTEGER JX            Cartesian coordinte index
C LOGICAL MOVALL        Move all atoms?
C INTEGER NAM           Number of atoms to move
C INTEGER NCR           Parameter defined above
C INTEGER NEM           Number of extended-mesh cells
C INTEGER NEMX(MX)      Extended-mesh cells in each mesh direction
C INTEGER NM            Number of mesh cells
C INTEGER NMX(MX)       Mesh cells in each mesh direction
C INTEGER NNAMAX        Input value of NNA (size of JAN, R2IJ and XIJ)
C INTEGER NNM           Number of neighbour mesh cells
C INTEGER NNMMAX        Maximum number of neighbour mesh cells
C INTEGER NNX(MX)       Neighbour-cell ranges
C LOGICAL NULCEL        Null cell?
C REAL*8  R2            Squared distance between two atoms
C REAL*8  RANGE2        Square of range
C REAL*8  RCELL(MX*MX)  Reciprocal cell vectors
C         RECCEL()      Finds reciprocal lattice vectors
C REAL*8  RNGMAX        Maximum range
C REAL*8  RMCELL(MX*MX) Reciprocal mesh-cell vectors
C REAL*8  RRANGE        Slightly reduced range
C REAL*8  XDIFF         Range of atom coordinates
C REAL*8  XMAX          Maximum atom coordinate
C REAL*8  XMIN          Minimum atom coordinate

      INTEGER
     .  IA, IAM, IEM, IM, 
     .  IN, IX, IXX, JA, JEM, JM, JX,
     .  NAM, NM, NEM, NNAMAX, NNM, NNMMAX, MAXNA
c
c     Auxiliary variable to avoid compiler warnings
c
      integer j_aux   

      integer, dimension(:), allocatable, save ::
     .  INX, I1NX, I2NX, J1NX, J2NX

      integer, dimension(:), allocatable, save ::
     .  IANEXT, IAPREV, IEMA, I1EMX, I2EMX, IMX, I1MX, I2MX,
     .  NEMX, NMX, NNX, IA1M, IMESH, IDNM

      real(dp)
     .  DISMIN, DOT, DPLANE, 
     .  R2, RANGE2, RNGMAX, RRANGE,
     .  XDIFF, XMARG, XMAX, XMIN

      double precision, dimension(:), allocatable, save ::
     .  DMX, DX, DX0M

      double precision, dimension(:), allocatable, save ::
     .  CELMSH, RCELL, RMCELL

      double precision, dimension(:,:), allocatable, save ::
     .  DXAM, DXNM

      LOGICAL
     .  FRSTME, INSIDE, MOVALL, NULCEL

      EXTERNAL
     .  DISMIN, DOT, RECCEL, memory

      SAVE
     .  FRSTME, IAM, IEM, IM, 
     .  NEM, NM, NNM, RANGE2, RNGMAX, RRANGE, MAXNA

      DATA FRSTME /.TRUE./
      DATA MAXNA  / 0 /

C Allocate local memory - check for change in number of atoms
C and if there has been one then re-initialise
      if (NA.gt.MAXNA) then
        if (allocated(IANEXT)) then
          call memory('D','I',size(IANEXT),'ranger')
          deallocate(IANEXT)
        endif
        if (allocated(IAPREV)) then
          call memory('D','I',size(IAPREV),'ranger')
          deallocate(IAPREV)
        endif
        if (allocated(IEMA)) then
          call memory('D','I',size(IEMA),'ranger')
          deallocate(IEMA)
        endif
        if (allocated(DXAM)) then
          call memory('D','D',size(DXAM),'ranger')
          deallocate(DXAM)
        endif
        allocate(IANEXT(NA))
        call memory('A','I',NA,'ranger')
        allocate(IAPREV(NA))
        call memory('A','I',NA,'ranger')
        allocate(IEMA(NA))
        call memory('A','I',NA,'ranger')
        allocate(DXAM(NX,NA))
        call memory('A','D',NX*NA,'ranger')
        MAXNA = NA
        FRSTME = .FALSE.
      endif
      if (.not.allocated(I1EMX)) then
        allocate(I1EMX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(I2EMX)) then
        allocate(I2EMX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(IMX)) then
        allocate(IMX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(I1MX)) then
        allocate(I1MX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(I2MX)) then
        allocate(I2MX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(NEMX)) then
        allocate(NEMX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(NMX)) then
        allocate(NMX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(NNX)) then
        allocate(NNX(NX))
        call memory('A','I',NX,'ranger')
      endif
      if (.not.allocated(CELMSH)) then
        allocate(CELMSH(NX*NX))
        call memory('A','D',NX*NX,'ranger')
      endif
      if (.not.allocated(RCELL)) then
        allocate(RCELL(NX*NX))
        call memory('A','D',NX*NX,'ranger')
      endif
      if (.not.allocated(RMCELL)) then
        allocate(RMCELL(NX*NX))
        call memory('A','D',NX*NX,'ranger')
      endif
      allocate(INX(NX))
      call memory('A','I',NX,'ranger')
      allocate(I1NX(NX))
      call memory('A','I',NX,'ranger')
      allocate(I2NX(NX))
      call memory('A','I',NX,'ranger')
      allocate(J1NX(NX))
      call memory('A','I',NX,'ranger')
      allocate(J2NX(NX))
      call memory('A','I',NX,'ranger')
      allocate(DMX(NX))
      call memory('A','D',NX,'ranger')
      allocate(DX(NX))
      call memory('A','D',NX,'ranger')
      allocate(DX0M(NX))
      call memory('A','D',NX,'ranger')

C Cell-mesh initialization section
      IF (MODE.EQ.'CELL' .OR. MODE.EQ.'cell' .OR.
     .    FRSTME .OR. RANGE.GT.RNGMAX) THEN

C       Start time counter (this is for debugging)
*       CALL TIMER( 'RANGER1', 1 )

C       Store range for comparison in subsequent calls
        RNGMAX = RANGE

C       Reduce the range slitghtly to avoid numerical-roundoff
C       ambiguities
        RRANGE = RANGE * (1.D0 - EPS)
        RANGE2 = RRANGE**2

C       Check if CELL must be generated automatically
        NULCEL = .TRUE.
        DO 20 IXX = 1,NX*NX
          IF (CELL(IXX) .NE. 0.D0) NULCEL = .FALSE.
   20   CONTINUE
        IF (NULCEL) THEN
          DO 40 IX = 1,NX
C           Find atom position bounds
            XMIN =  1.D30
            XMAX = -1.D30
            DO 30 IA = 1,NA
              XMIN = MIN( XMIN, XA(IX,IA) )
              XMAX = MAX( XMAX, XA(IX,IA) )
   30       CONTINUE
C           Determine 'cell margins' to prevent intercell interactions
            XDIFF = XMAX - XMIN
            XMARG = MAX( RANGE*DXRANG, XDIFF*DXMARG )
C           Define orthorrombic cell
            IXX = IX + NX * (IX-1)
            CELL(IXX) = XDIFF + 2.D0 * XMARG
   40     CONTINUE
        ENDIF

C       Find reciprocal cell vectors (not multiplied by 2*pi)
        CALL RECCEL( NX, CELL, RCELL, 0 )

C       Find number of mesh divisions
        NM = 1
        DO 50 IX = 1,NX
          IXX = 1 + NX * (IX-1)
          DPLANE = 1.D0 / SQRT( DOT( RCELL(IXX), RCELL(IXX), NX ) )
          NMX(IX) = 0.999D0 * DPLANE / (RRANGE / NCR)
          IF (NMX(IX) .LE. 0) NMX(IX) = 1
          NM = NM * NMX(IX)
   50   CONTINUE

C       Find mesh-cell vectors
        IXX = 0
        DO 70 IX = 1,NX
          DO 60 JX = 1,NX
            IXX = IXX + 1
            CELMSH(IXX) = CELL(IXX)  / NMX(IX)
            RMCELL(IXX) = RCELL(IXX) * NMX(IX)
   60     CONTINUE
   70   CONTINUE

C       Find index-range of neighbour mesh cells and of extended mesh
        NNM = 1
        NEM = 1
        DO 80 IX = 1,NX
          IXX = 1 + NX * (IX-1)
          DPLANE = 1.D0 / SQRT( DOT(RMCELL(IXX),RMCELL(IXX),NX) )
          NNX(IX) = RRANGE / DPLANE + 1
          J1NX(IX) = 0
          J2NX(IX) = 1
          I1NX(IX) = - NNX(IX)
          I2NX(IX) = + NNX(IX)
          I1MX(IX) = 0
          I2MX(IX) = NMX(IX) - 1
          I1EMX(IX) = - NNX(IX)
          I2EMX(IX) = NMX(IX) + NNX(IX) - 1
          NEMX(IX) = NMX(IX) + 2 * NNX(IX)
          NNM = NNM * (1+2*NNX(IX))
          NEM = NEM * NEMX(IX)
   80   CONTINUE

C  Allocate arrays whose dimensions are now known
        if (allocated(IA1M)) then
          call memory('D','I',size(IA1M),'ranger')
          deallocate(IA1M)
        endif
        if (allocated(IDNM)) then
          call memory('D','I',size(IDNM),'ranger')
          deallocate(IDNM)
        endif
        if (allocated(DXNM)) then
          call memory('D','D',size(DXNM),'ranger')
          deallocate(DXNM)
        endif
        if (allocated(IMESH)) then
          call memory('D','I',size(IMESH),'ranger')
          deallocate(IMESH)
        endif
        allocate(IA1M(NM))
        call memory('A','I',NM,'ranger')
        allocate(IDNM(NNM))
        call memory('A','I',NNM,'ranger')
        allocate(DXNM(NX,NNM))
        call memory('A','D',NX*NNM,'ranger')
        allocate(IMESH(NEM))
        call memory('A','I',NEM,'ranger')

C       Find which mesh cells are actually within range
        NNMMAX = NNM
        NNM = 0
        DO 170 IN = 1,NNMMAX
          j_aux = in
          CALL INDARR( -1, NX, I1NX, I2NX, INX, 1, j_aux )
          INSIDE = .TRUE.
C         From here to label 130 is generally not worth unless CELL is
C         very nonorthorrombic (like fcc, bcc or hex) and changes rarely
*         DO 120 IV = 1,2**NX
*           j_aux = iv
*           CALL INDARR( -1, NX, J1NX, J2NX, IVX, 1, j_aux )
*           DO 100 IX = 1,NX
*             DXM(IX) = 0.D0
*             DO 90 JX = 1,NX
*               IXX = IX + NX * (JX-1)
*               DXM(IX) = DXM(IX) + CELMSH(IXX) * (INX(JX)+IVX(JX))
*  90         CONTINUE
* 100       CONTINUE
*           DRM = DISMIN( NX, CELMSH, DXM )
*           IF (DRM .LT. RRANGE) THEN
*             INSIDE = .TRUE.
*             GOTO 130
*           ENDIF
* 120     CONTINUE
*           INSIDE = .FALSE.
* 130     CONTINUE
          IF (INSIDE) THEN
            NNM = NNM + 1
C           IDNM is the extended-mesh-index distance between
C           neighbour mesh cells
            IDNM(NNM) = INX(NX)
            DO 140 IX = NX-1,1,-1
              IDNM(NNM) = INX(IX) + NEMX(IX) * IDNM(NNM)
  140       CONTINUE
C           DXNM is the vector distance between neighbour mesh cells
            DO 160 IX = 1,NX
              DXNM(IX,NNM) = 0.D0
              DO 150 JX = 1,NX
                IXX = IX + NX * (JX-1)
                DXNM(IX,NNM) = DXNM(IX,NNM) + CELMSH(IXX) * INX(JX)
  150         CONTINUE
  160       CONTINUE
          ENDIF
  170   CONTINUE

C       Find correspondence between extended and reduced (normal) meshes
        DO 180 IEM = 1,NEM
          j_aux = iem
          CALL INDARR( -1, NX, I1EMX, I2EMX, IMX, 1, j_aux )
          CALL INDARR( +1, NX, I1MX,  I2MX,  IMX, 1, IM  )
          IMESH(IEM) = IM
  180   CONTINUE

C       Stop time counter
*       CALL TIMER( 'RANGER1', 2 )

C       Set 'move all atoms' switch
        MOVALL = .TRUE.
      ELSE
        MOVALL = .FALSE.
      ENDIF
C End of cell initialization section

C Atom-positions (relative to mesh) initialization section
      IF (MODE.EQ.'MOVE' .OR. MODE.EQ.'move' .OR. MOVALL) THEN
        IF (NAMOVE .EQ. NA) MOVALL = .TRUE.

C       Start time counter
*       CALL TIMER( 'RANGER2', 1 )

        IF (MOVALL) THEN
          NAM = NA
C         Initialize 'atoms in mesh-cell' lists
          DO 200 IA = 1,NA
            IANEXT(IA) = 0
            IAPREV(IA) = 0
  200     CONTINUE
          DO 210 IM = 1,NM
            IA1M(IM) = 0
  210     CONTINUE
        ELSE
          NAM = NAMOVE
        ENDIF

C       Loop on moved atoms
        DO 240 IAM = 1,NAM

C         Select atom to move
          IF (MOVALL) THEN
            IA = IAM
          ELSE
            IA = IAMOVE(IAM)
C           Supress atom from its previous mesh-cell
            JA = IAPREV(IA)
            IF (JA.NE.0) IANEXT(JA) = IANEXT(IA)
            JA = IANEXT(IA)
            IF (JA.NE.0) IAPREV(JA) = IAPREV(IA)
            IEM = IEMA(IA)
            IM = IMESH(IEM)
            IF (IA1M(IM) .EQ. IA) IA1M(IM) = JA
          ENDIF
          
C         Find mesh-cell in which atom is
          DO 220 IX = 1,NX
            IXX = 1 + NX * (IX-1)
            DMX(IX) = DOT( RMCELL(IXX), XA(1,IA), NX )
            IMX(IX) = INT( DMX(IX) + 1000.D0 ) - 1000
            DMX(IX) = DMX(IX) - IMX(IX)
            IMX(IX) = MOD( IMX(IX) + 1000 * NMX(IX), NMX(IX) )
  220     CONTINUE
          CALL INDARR( +1, NX, I1EMX, I2EMX, IMX, 1, IEM )
          CALL INDARR( +1, NX, I1MX,  I2MX,  IMX, 1, IM )
          IEMA(IA) = IEM

C         Put atom first in its new mesh-cell
          JA = IA1M(IM)
          IF (JA .NE. 0) IAPREV(JA) = IA
          IANEXT(IA) = JA
          IA1M(IM) = IA

C         Find atomic position relative to mesh
          DO 230 IX = 1,NX
            DXAM(IX,IA) = 0.D0
            DO 225 JX = 1,NX
              IXX = IX + NX * (JX-1)
              DXAM(IX,IA) = DXAM(IX,IA) + CELMSH(IXX) * DMX(JX)
  225       CONTINUE
  230     CONTINUE

  240   CONTINUE

C       Stop time counter
*       CALL TIMER( 'RANGER2', 2 )
      ENDIF
C End of atom-positions initialization section

C Search section
      IF (MODE.EQ.'FIND' .OR. MODE.EQ.'find') THEN
        RRANGE = RANGE * (1.D0 - EPS)
        RANGE2 = RANGE**2

C       Find the mesh cell of the center of the sphere
        IF (IA0.LE.0) THEN
C         Find mesh cell of position X0
          DO 250 IX = 1,NX
            IXX = 1 + NX * (IX-1)
            DMX(IX) = DOT( RMCELL(IXX), X0, NX )
            IMX(IX) = INT( DMX(IX) + 1000.D0 ) - 1000
            DMX(IX) = DMX(IX) - IMX(IX)
            IMX(IX) = MOD( IMX(IX) + 1000 * NMX(IX), NMX(IX) )
  250     CONTINUE
          CALL INDARR( +1, NX, I1EMX, I2EMX, IMX, 1, IEM )
          DO 270 IX = 1,NX
            DX0M(IX) = 0.D0
            DO 260 JX = 1,NX
              IXX = IX + NX * (JX-1)
              DX0M(IX) = DX0M(IX) + CELMSH(IXX) * DMX(JX)
  260       CONTINUE
  270     CONTINUE
        ELSE
C         Find mesh cell of atom IA0
          IEM = IEMA(IA0)
          DO 280 IX = 1,NX
            DX0M(IX) = DXAM(IX,IA0)
  280     CONTINUE
        ENDIF

C       Loop on neighbour mesh cells and on the atoms within them
C       This is usually the only time-consuming loop
        NNAMAX = NNA
        NNA = 0
        DO 330 IN = 1,NNM
          JEM = IEM + IDNM(IN)
          JM = IMESH(JEM)
C         Loop on atoms of neighbour cell.
C         Try first atom in this mesh-cell
          JA = IA1M(JM)
  300     CONTINUE
          IF (JA .NE. 0) THEN
C           Check that single-counting exclusion does not apply
            IF (IA0.LE.0 .OR. ISC.EQ.0 .OR. JA.LE.IA0) THEN
C             Find vector and distance to atom JA
              R2 = 0.D0
              DO 310 IX = 1,NX
                DX(IX) = DXNM(IX,IN) + DXAM(IX,JA) - DX0M(IX)
                R2 = R2 + DX(IX)**2
  310         CONTINUE
C             Check if atom JA is within range
              IF (R2 .LE. RANGE2) THEN
                NNA = NNA + 1
C               Check that array arguments are not overflooded
                IF (NNA .LE. NNAMAX) THEN
                  JAN(NNA) = JA
                  DO 320 IX = 1,NX
                    XIJ(IX,NNA) = DX(IX)
  320             CONTINUE
                  R2IJ(NNA) = R2
                ENDIF
              ENDIF
            ENDIF
C           Take next atom in this mesh-cell and go to begining of loop
            JA = IANEXT(JA)
            GOTO 300
          ENDIF
  330   CONTINUE
      ENDIF
C End of search section

C Deallocate local memory
      call memory('D','I',size(INX),'ranger')
      deallocate(INX)
      call memory('D','I',size(I1NX),'ranger')
      deallocate(I1NX)
      call memory('D','I',size(I2NX),'ranger')
      deallocate(I2NX)
      call memory('D','I',size(J1NX),'ranger')
      deallocate(J1NX)
      call memory('D','I',size(J2NX),'ranger')
      deallocate(J2NX)
      call memory('D','D',size(DMX),'ranger')
      deallocate(DMX)
      call memory('D','D',size(DX),'ranger')
      deallocate(DX)
      call memory('D','D',size(DX0M),'ranger')
      deallocate(DX0M)

C This is the unique return point
      FRSTME = .FALSE.

      RETURN
      END 

      SUBROUTINE RECCEL( N, A, B, IOPT )

C  CALCULATES RECIPROCAL LATTICE VECTORS B.
C  THEIR PRODUCT WITH DIRECT LATTICE VECTORS A IS 1 (IF IOPT=0) OR
C  2*PI (IF IOPT=1). N IS THE SPACE DIMENSION.
C  WRITTEN BY J.M.SOLER.

      use precision
      implicit none

      integer, intent(in) :: n, iopt
      real(dp) A(N,N),B(N,N)

      real(dp) :: c, ci
      integer i

      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*ACOS(-1.D0)

      IF (N .EQ. 1) THEN
        B(1,1) = C / A(1,1)
      ELSEIF (N .EQ. 2) THEN
        C = C / (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        B(1,1) =  A(2,2)*C
        B(1,2) = (-A(2,1))*C
        B(2,1) = (-A(1,2))*C
        B(2,2) =  A(1,1)*C
      ELSEIF (N .EQ. 3) THEN
        B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
        B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
        B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
        B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
        B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
        B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
        B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
        B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
        B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
        DO 20 I=1,3
          CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
          B(1,I)=B(1,I)*CI
          B(2,I)=B(2,I)*CI
          B(3,I)=B(3,I)*CI
  20    CONTINUE
      ELSE
        WRITE(6,*) 'RECCEL: NOT PREPARED FOR N =', N
        STOP
      ENDIF
      END




      SUBROUTINE INDARR( IOPT, ND, I1, I2, I, J1, J )

C ********************************************************************
C Finds the global index in a multidimensional array from the idexes
C in each dimension, or viceversa (the first is an explicit solution
C of the standard index-resolution problem that the compiler solves
C each time an array element is referenced).
C Written by J.M.Soler. Nov'96.
C *********** INPUT **************************************************
C INTEGER IOPT   : IOPT>0 => from I to J. IOPT<0 => from J to I.
C INTEGER ND     : Number of array dimensions (indexes)
C INTEGER I1(ND) : Minimum value of array indexes
C INTEGER I2(ND) : Maximum value of array indexes (i.e. ARRAY(I1:I2))
C INTEGER J1     : Minimum value of global index (typically 1)
C *********** INPUT OR OUTPUT (DEPENDING OF IOPT) ********************
C INTEGER I(ND)  : Array indexes in each dimension
C INTEGER J      : Global index (first dimension increasing fastest)
C *********** BEHAVIOUR **********************************************
C Indexes I() are taken as periodic, i.e. their modulus I2(ID)-I1(ID)+1
C   is taken before using them. This simplifies its use as indexes of a
C   mesh with periodic boundary conditions. This modulus operation is
C   also done with J, so that the output I() are always within range.
C If IOPT=0, nothing is done.
C *********** USAGE **************************************************
C 	Sample usage to find the Laplacian of a function defined in a mesh 
C with periodic boundary conditions in a space of variable dimension
C    SUBROUTINE LAPLACIAN( ND, N, DX, F, FLAPL )
C    PARAMETER (MAXD = 3)
C    DIMENSION N(ND),DX(ND),F(*),FLAPL(*),I1(MAXD),I2(MAXD),I(MAXD)
C    NMESH = 1
C    DO ID = 1,ND
C      I1(ID) = 1
C      I2(ID) = N(ID)
C      NMESH = NMESH * N(ID)
C    ENDDO
C    DO IMESH = 1, NMESH
C      CALL INDARR( -1, ND, I1, I2, I, 1, IMESH )
C      FLAPL(IMESH) = 0.
C      DO ID = 1,ND
C        DO K = -1,1,2
C          I(ID) = I(ID) + K
C          CALL INDARR( +1, ND, I1, I2, I, 1, JMESH )
C          FLAPL(IMESH) = FLAPL(IMESH) + F(JMESH) / DX(ID)**2
C          I(ID) = I(ID) - K
C        ENDDO
C        FLAPL(IMESH) = FLAPL(IMESH) - 2. * F(IMESH) / DX(ID)**2
C      ENDDO
C    ENDDO
C    END
C *********** LANGUAGE ***********************************************
C Illegal ANSI Fortran77 features used: IMPLICIT NONE
C Illegal ANSI Fortran90 features used: none.
C ********************************************************************

C Next line is non-standard but may be supressed
      use precision

      IMPLICIT NONE
      INTEGER  ND
      INTEGER  I(ND), I1(ND), I2(ND), ID, IOPT, J, J1, K, N

      IF (IOPT .GT. 0) THEN
        J = 0
        DO 10 ID = ND,1,-1
          N = I2(ID) - I1(ID) + 1
          K = I(ID) - I1(ID)
          K = MOD( K + 1000 * N, N )
          J = K + N * J
   10   CONTINUE
        J = J + J1
      ELSEIF (IOPT .LT. 0) THEN
        K = J - J1
        DO 20 ID = 1,ND
          N = I2(ID) - I1(ID) + 1
          I(ID) = I1(ID) + MOD( K, N )
          K = K / N
   20   CONTINUE
      ENDIF

      END



