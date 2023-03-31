! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
C**************************************************************************
C  This module contains subroutines required for Z-matrix based geometry  *
C  optimisation in SIESTA.                                                *
C                                                                         *
C  Written by Rainer Hoft (UTS) & Julian Gale (Curtin), March - May 2005  *
C**************************************************************************
      module zmatrix
        use alloc,     only : re_alloc, de_alloc
        use precision, only : dp
        use units,     only : Ang, eV, deg, pi
        use sys,       only : die
        implicit none
C  Variables

C  iZmat holds the pointers to the atoms to which the components
C  of Zmat relate
        integer,  pointer, save :: iZmat(:)
C  VaryZmat(3*na) indicates whether the component of Zmat is fixed,
        logical,  pointer, save :: VaryZmat(:)
C  Zmat(3*na) holds the Z-matrix components for each atom
        real(dp), pointer, save :: Zmat(:)
C  ZmatForce(3*na) holds the forces on each Zmat coordinate
        real(dp), pointer, save :: ZmatForce(:)
C  ZmatForceVar(nVar) holds the forces on additional constrained coordinates
        real(dp), pointer, save :: ZmatForceVar(:)
C  lUseZmatrix indicates whether Zmatrix approach is to be used
        logical,           save :: lUseZmatrix = .false.

C  nZmol is the number of molecules with Z matrices defined
        integer,           save :: nZmol
C  nZmolAtoms(nZmol) indicates the number of atoms in each molecule
        integer,  pointer, save :: nZmolAtoms(:)
C  nZmolStartAtom(nZmol) is a pointer to the first atom of each molecule
        integer,  pointer, save :: nZmolStartAtom(:)
C  nZcart is the number of cartesian blocks defined in the Zmatrix
        integer,           save :: nZcart
C  nZcartAtoms(nZcart) indicates the number of atoms in each cartesian block
        integer,  pointer, save :: nZcartAtoms(:)
C  nZcartStartAtom(nZcart) is a pointer to the first atom of each cartesian block
        integer,  pointer, save :: nZcartStartAtom(:)
C  Variables to specify the input (and output) units of the Zmatrix
        integer,           save :: ZmatUnitsLen
        integer,           save :: ZmatUnitsAng
C  Variables to specify the maximum force tolerance for CG optimization
        real(dp),          save :: ZmatForceTolLen
        real(dp),          save :: ZmatForceTolAng
C  Variables to specify the maximum displacement during CG optimization
        real(dp),          save :: ZmatMaxDisplLen
        real(dp),          save :: ZmatMaxDisplAng
C  ZmatType(3*na) defines the component type:
C    1 - angle
C    2 - bond length
C    3 - pure cartesian (molecule)
C    4 - scaled cartesian (molecule)
C    5 - fractional cartesian (molecule)
C    6 - pure cartesian
C    7 - scaled cartesian (scaled by lattice constant)
C    8 - fractional cartesian (scaled by lattice vectors)
        integer,  pointer, save :: ZmatType(:)
C  nVars specifies the number of constrained variables
        integer,           save :: nVars
C  ZmatVarNames(nVars) gives the names of the constrained variables
        character(len=10), pointer, save :: ZmatVarNames(:)
C  iZmattoVars(3*na) is the index in the variable array of a 
C     symbolic coordinate
        integer, pointer,  save :: iZmattoVars(:) 
C  iVarstoZmat(nVars) is the Zmat index of the first symbolic 
C    coordinate corresponding to the variable
        integer, pointer,  save :: iVarstoZmat(:)
C  lZmatVarsDef(nVars) tells us whether the variable name has been 
C    given a value
        logical, pointer,  save :: lZmatVarsDef(:)
C  lCalcAllForces specifies whether forces for fixed coorinates should
C    be calculated
        logical,           save :: lCalcAllForces
C  iNextDept(3*na) gives the linked list of dependencies on symbolic coordinates
C    (note it is set up so that each coordinate can only depend on one other)
        integer, pointer,  save :: iNextDept(:)
C  coeffA and coeffB contains the linear relationships between variables
C       v2 = a*v1 + B
        real(dp), pointer, save :: coeffA(:)
        real(dp), pointer, save :: coeffB(:)


      ! AG: convenience variables
        integer, pointer,  save :: species(:) 
        real(dp),  save :: scale_length, scale_angle 
        real(dp),  save :: zmatrix_alat
        logical,   save :: unknown_cell

      public :: read_Zmatrix, lUseZmatrix, iofaZmat
      public :: CartesianForce_to_ZmatForce, write_Zmatrix
      public :: write_canonical_Zmatrix
      public :: write_canonical_ucell_and_Zmatrix
      public :: VaryZmat, Zmat, ZmatForce, ZmatForceVar
      public :: iZmattoVars, ZmatType, Zmat_to_Cartesian
      public :: coeffA, coeffB, iNextDept
      public :: ZmatForceTolLen, ZmatForceTolAng
      public :: ZmatMaxDisplLen, ZmatMaxDisplAng
      public :: nZmol, nZmolStartAtom, nZmolAtoms
      public :: nZcart, nZcartStartAtom, nZcartAtoms

      private 

      contains

C  Subroutines

        subroutine read_Zmatrix(na,nSpecies,alat,unit_cell,
     $                          lOrigin,origin)

        use fdf
        use sys,         only : die
        use parallel,    only : Node

C Passed variables
        integer,  intent(in)                :: na
        integer,  intent(out)               :: nSpecies(na)
        real(dp), intent(in)                :: alat
        real(dp), intent(in)                :: unit_cell(3,3)
        logical,  intent(in)                :: lOrigin
        real(dp), intent(in)                :: origin(3)

C Local variables
        character(len= 10)            :: angleunits
        character(len= 10)            :: lengthunits
        integer                       :: i
        integer                       :: j
        integer                       :: m
        integer                       :: k
        integer                       :: nStart
        integer                       :: nAtoms
!        logical                       :: eof
        logical,                 save :: firsttime = .true.
        integer                       :: input_type
        integer                       :: units_type
        character(len=10), parameter  :: unitslen_default = 'Bohr'
        character(len=10), parameter  :: unitsang_default = 'rad' 
        real(dp),          parameter  :: ftollen_default = 1.55574d-3
        real(dp),          parameter  :: ftolang_default = 3.56549d-3
        real(dp),          parameter  :: dxmaxlen_default = 0.2d0
        real(dp),          parameter  :: dxmaxang_default = 0.003d0
        logical,           parameter  :: lCalcAllForces_default=.false.
        logical                       :: found
        integer                       :: nnames
        integer                       :: nreals
        real(dp)                      :: ZmatVal(3)
        real(dp)                      :: coeffBVal(3)
        integer                       :: indi
        integer                       :: indi1
        integer                       :: flagnr
        real(dp)                      :: A
        real(dp)                      :: B

        type(block_fdf)            :: bfdf
        type(parsed_line), pointer :: pline

        real(dp), external    :: volcel
        real(dp)              :: dummy_volume

       zmatrix_alat = alat     ! AG: For future processing

       dummy_volume = volcel(unit_cell)
       unknown_cell = (abs(dummy_volume) < 1.0e-6_dp)

C Nullify pointers
        if (firsttime) then
          nullify(iZmat,nZmolAtoms,nZmolStartAtom)
          nullify(nZcartAtoms,nZcartStartAtom)
          nullify(Zmat,VaryZmat,ZmatForce,ZmatForceVar)
          nullify(ZmatType,iZmattoVars,iVarstoZmat)
          nullify(lZmatVarsDef,iNextDept)
          nullify(coeffA,coeffB)
          nullify(species)

C Allocate Zmatrix arrays
          call re_alloc( Zmat,       1, 3*na, 'Zmat',      'zmatrix' )
          call re_alloc( species,    1, 3*na, 'species',   'zmatrix' )
          call re_alloc( ZmatForce,  1, 3*na, 'ZmatForce', 'zmatrix' )
          call re_alloc( iZmat,      1, 3*na, 'iZmat',     'zmatrix' )
          call re_alloc( VaryZmat,   1, 3*na, 'VaryZmat',  'zmatrix' )
          call re_alloc( nZmolAtoms, 1, na, 'nZmolAtoms',  'zmatrix' )
          call re_alloc( nZmolStartAtom, 1, na, 'nZmolStartAtom',
     &                   'zmatrix' )
          call re_alloc( nZcartAtoms, 1, na, 'nZcartAtoms', 'zmatrix' )
          call re_alloc( nZcartStartAtom, 1, na, 'nZcartStartAtom',
     &                   'zmatrix' )
          call re_alloc( ZmatType, 1, 3*na, 'ZmatType', 'zmatrix' )
          call re_alloc( iZmattoVars, 1, 3*na, 'iZmattoVars',
     &                   'zmatrix' )
          call re_alloc( iNextDept, 1, 3*na, 'iNextDept', 'zmatrix' )
          call re_alloc( coeffA, 1, 3*na, 'coeffA', 'zmatrix' )
          call re_alloc( coeffB, 1, 3*na, 'coeffB', 'zmatrix' )
          call re_alloc( ZmatForceVar, 1, 3*na, 'ZmatForceVar',
     &                   'zmatrix' )
          call re_alloc( iVarstoZmat, 1, 3*na, 'iVarstoZmat',
     &                   'zmatrix' )
          call re_alloc( lZmatVarsDef, 1, 3*na, 'lZmatVarsDef',
     &                   'zmatrix' )
          nullify(ZmatVarNames)
          allocate(ZmatVarNames(3*na))
!          call re_alloc( ZmatVarNames, 1, 3*na, 'ZmatVarNames',
!     &                   'zmatrix' )

C Read units for lengths and angles
          lengthunits = fdf_string('ZM.UnitsLength',unitslen_default)
          if (leqi(lengthunits,'Bohr')) then
            ZmatUnitsLen = 0
          elseif (leqi(lengthunits,'Ang').or.
     .            leqi(lengthunits,'Angstrom')) then
            ZmatUnitsLen = 1
          else
            call die('Invalid Zmatrix length units')
          endif

          angleunits = fdf_string('ZM.UnitsAngle',unitsang_default)
          if (leqi(angleunits,'rad').or.leqi(angleunits,'radians')) 
     .        then
            ZmatUnitsAng = 0
          elseif (leqi(angleunits,'deg').or.
     .            leqi(angleunits,'degrees')) then
            ZmatUnitsAng = 1
          else
            call die('Invalid Zmatrix angular units')
          endif

C Read maximum force tolerance for lengths and angles        
          ZmatForceTolLen = fdf_physical('ZM.ForceTolLength',
     .                                   ftollen_default,'Ry/Bohr')
          ZmatForceTolAng = fdf_physical('ZM.ForceTolAngle',
     .                                   ftolang_default,'Ry/rad')

C Read maximum displacement per CG step for lengths and angles
          ZmatMaxDisplLen = fdf_physical('ZM.MaxDisplLength',
     .                                   dxmaxlen_default,'Bohr')
          ZmatMaxDisplAng = fdf_physical('ZM.MaxDisplAngle',
     .                                   dxmaxang_default,'rad')


C Check if we should calculate all the forces
          lCalcAllForces = fdf_boolean('ZM.CalcAllForces',
     .                                  lCalcAllForces_default)
        endif

C Check whether a Z-matrix has been input
        lUseZmatrix = fdf_block('Zmatrix',bfdf)

C If not Z matrix return
        if (.not.lUseZmatrix) goto 999

C Output information about Z matrix units only for Z matrix case
        if (Node.eq.0) then
          write (6,'(''read_Zmatrix: Length units: '',a10)') lengthunits
          write (6,'(''read_Zmatrix: Angle  units: '',a10)') angleunits
          write (6,'(''read_Zmatrix: Force tolerances:'')')
          write (6,'(''read_Zmatrix:    for lengths = '',
     .                f12.6,'' Ry/Bohr'')') ZmatForceTolLen
          write (6,'(''read_Zmatrix:    for angles  = '',
     .                f12.6,'' Ry/rad'',/)') ZmatForceTolAng
          write (6,'(''read_Zmatrix: Maximum displacements:'')')
          write (6,'(''read_Zmatrix:    for lengths = '',
     .               f12.6,'' Bohr'')') ZmatMaxDisplLen
          write (6,'(''read_Zmatrix:    for angles  = '',
     .                 f12.6,'' rad'')') ZmatMaxDisplAng
        endif

C Initialise number of molecules and Cartesian blocks
        nZmol = 0
        nZcart = 0
        nVars = 0

C Read data block up to endblock statement
!
!   NOTE: We use a function to search for substrings
!         in the line, NOT the exact token match
!         implied by "fdf_bsearch".
!
!        eof = .false.
        nAtoms = 0
        if (Node.eq.0) write(*,"(a)") "%block Zmatrix"

        do while (fdf_bline(bfdf,pline))
          if (Node.eq.0) write(*,"(a)") trim(pline%line)

          if (fdf_substring_search(pline, 'molecule')) then
C Found new molecule          
            nZmol = nZmol + 1
            nZmolAtoms(nZmol) = 0
            nZmolStartAtom(nZmol) = nAtoms + 1
            input_type = 0
            if (fdf_substring_search(pline, 'scale')) then
              units_type = 1
            elseif (fdf_substring_search(pline, 'frac')) then
              if (unknown_cell) then
                call die("Cannot use fractional coordinates " //
     $                   " if unit cell is not specified")
              endif
              units_type = 2
            else
              units_type = 0
            endif
          elseif (fdf_substring_search(pline, 'cart')) then
C Found new cartesian block
            nZcart = nZcart + 1
            nZcartAtoms(nZcart) = 0
            nZcartStartAtom(nZcart) = nAtoms + 1
            input_type = 1
            units_type = 0
          elseif (fdf_substring_search(pline, 'scale')) then
C Found new cartesian(scaled) block
            nZcart = nZcart + 1
            nZcartAtoms(nZcart) = 0
            nZcartStartAtom(nZcart) = nAtoms + 1
            input_type = 1
            units_type = 1
          elseif (fdf_substring_search(pline, 'frac')) then
C Found new cartesian(fractional) block
!AG
            if (unknown_cell) then
              call die("Cannot use fractional coordinates " //
     $                 " if unit cell is not specified")
            endif
            nZcart = nZcart + 1
            nZcartAtoms(nZcart) = 0
            nZcartStartAtom(nZcart) = nAtoms + 1
            input_type = 1
            units_type = 2
          elseif (fdf_substring_search(pline, 'constant')) then
            input_type = 2
          elseif (fdf_substring_search(pline, 'variable')) then
            input_type = 3
          elseif (fdf_substring_search(pline, 'constraint')) then
            input_type = 4

          else 
C Check that the lines contains the correct type of data
            if (input_type.eq.0) then
C Molecule (1-4: integers, 5-7: reals or symbolic names)
              if (.not. (fdf_bmatch(pline, 'iiiisss'))) then
                call die('read_Zmatrix: (iiiisss) Error in syntax') 
              endif
            elseif (input_type.eq.1) then
C Cartesian (1: integer, 2-4: reals or symbolic names)
              if (.not. (fdf_bmatch(pline, 'isss') )) then
                call die ('read_Zmatrix: (isss) Error in syntax')
              endif
            elseif (input_type.eq.2.or.input_type.eq.3) then
C Constant/variable definition (1: string, 2: real)
              ! Should this be allowed to be 'nv'?
              if (.not. (fdf_bmatch(pline, 'nr'))) then
                call die ('read_Zmatrix: (nr) Error in syntax')
              endif
            elseif (input_type.eq.4) then
C Constraint definition (1-2: string, 3-4: reals)
              ! Should this be allowed to be 'nnvv'?
              if (.not. (fdf_bmatch(pline, 'nnrr'))) then
                call die ('read_Zmatrix: (nrrr) Error in syntax')
              endif
            endif

C Increment molecule/cartesian block count
            if (input_type.eq.0) then
              nZmolAtoms(nZmol) = nZmolAtoms(nZmol) + 1
            else if (input_type.eq.1) then
              nZcartAtoms(nZcart) = nZcartAtoms(nZcart) + 1
            endif

C Molecule/cartesian block
            if (input_type.eq.0.or.input_type.eq.1) then

C Increment atom number 
              nAtoms = nAtoms + 1
C Check number of atoms against input number
              if (nAtoms.gt.na) then
                call die('read_Zmatrix: Too many atoms in Z-matrix')
              endif

C Assign values for species and Z matrix
              nSpecies(nAtoms) = fdf_bintegers(pline,1)
            endif

C Read dependencies if its a molecule
            if (input_type.eq.0) then
              nStart = nZmolStartAtom(nZmol)
              iZmat(3*nAtoms-2) = fdf_bintegers(pline,2) + nStart-1
              iZmat(3*nAtoms-1) = fdf_bintegers(pline,3) + nStart-1
              iZmat(3*nAtoms)   = fdf_bintegers(pline,4) + nStart-1
            elseif (input_type.eq.1) then
              iZmat(3*nAtoms-2) = 1
              iZmat(3*nAtoms-1) = 1
              iZmat(3*nAtoms) = 1
            endif

C Assign type to coordinates
            if (input_type.eq.0) then
              if (nZmolAtoms(nZmol).eq.1) then
C Molecule, first atom
                do k= 1, 3
                  ZmatType(3*(nAtoms-1)+k) = 3+units_type
                enddo
              else
C Molecule, not first atom
                ZmatType(3*nAtoms-2) = 2
                ZmatType(3*nAtoms-1) = 1
                ZmatType(3*nAtoms)   = 1
              endif
            elseif (input_type.eq.1) then
C Cartesian block              
              do k= 1, 3 
                ZmatType(3*(nAtoms-1)+k) = 6+units_type
              enddo
            endif

C Read reals/symbols for values of coordinates
            if (input_type.eq.0.or.input_type.eq.1) then
              nreals = 0
              nnames = 0
              do k= 1, 3
                flagnr = 4-3*input_type+k
                indi = 3*(nAtoms-1)+k

                if (fdf_bmatch(pline, 'n', flagnr-1)) then
C String/symbol
                  nnames = nnames + 1
                  i = 1
                  found = .false.
                  do while (i.le.nVars.and..not.found)
                    found = leqi(fdf_bnames(pline, nnames),
     .                           ZmatVarNames(i))
                    i = i + 1
                  enddo
                  if (found) then
                    iZmattoVars(indi) = i-1
                    VaryZmat(indi) = .false.
                    j = iVarstoZmat(i-1) 
C Check length/angle dependency incompatibility
                    if (ZmatType(indi).eq.1.and.
     .                  ZmatType(j).ne.1) then
                      call die('read_Zmatrix: error - angle/length
     .                          dependency')                       
                    endif
                    do while (iNextDept(j).ne.0) 
                      j = iNextDept(j)
                    enddo
                    iNextDept(j) = indi
                  else
                    nVars = nVars + 1
                    lZmatVarsDef(nVars) = .false.
                    ZmatVarNames(nVars) = fdf_bnames(pline, nnames)
                    iZmattoVars(indi) = nVars
                    iVarstoZmat(nVars) = indi
                  endif
                  coeffA(indi) = 1.0d0
                  coeffB(indi) = 0.0d0
                else
C Explicit values                
                  nreals = nreals + 1
                  Zmat(indi) = fdf_breals(pline, nreals)
                  iZmattoVars(indi) = 0
                  coeffA(indi) = 1.0d0
                  coeffB(indi) = 0.0d0
C Optional flags for control of optimisation
                  if (fdf_bnintegers(pline).ge.flagnr) then
                    VaryZmat(indi) = (fdf_bintegers(pline,flagnr).ne.0)
                  else
                    VaryZmat(indi) = .true.
                  endif
                endif
                iNextDept(indi) = 0
              enddo
            endif

C Read variable/constant definitions
            if (input_type.eq.2.or.input_type.eq.3) then
              i = 1
              found = .false.
              do while (i.le.nVars.and..not.found) 
                found = leqi(ZmatVarNames(i),fdf_bnames(pline,1))
                i = i+1
              enddo
              if (found) then
                if (lZmatVarsDef(i-1)) then
                  call die('read_Zmatrix: Multiple definition of
     .                      Z-matrix symbol')
                else
                  indi = iVarstoZmat(i-1)
                  do while (indi.ne.0) 
                    Zmat(indi) = fdf_breals(pline,1)
                    indi = iNextDept(indi)
                  enddo
                  lZmatVarsDef(i-1) = .true.
                endif
                if (input_type.eq.2) then
C Fixed, symbolic                  
                  VaryZmat(iVarstoZmat(i-1)) = .false.
                else
C Vary, symbolic                  
                  VaryZmat(iVarstoZmat(i-1)) = .true.
                endif
              else
                call die('read_Zmatrix: Invalid symbol in Z-matrix')
              endif
            endif

C Read constraint definitions              
            if (input_type.eq.4) then
              i = 1
              found = .false.
              do while (i.le.nVars.and..not.found)
                found = leqi(ZmatVarNames(i),fdf_bnames(pline,1))
                i = i + 1
              enddo
              if (found) then
                if (lZmatVarsDef(i-1)) then
                  call die('read_Zmatrix: Multiple definition of
     .                      Z-matrix symbol')
                else
                  indi = iVarstoZmat(i-1)
                  lZmatVarsDef(i-1) = .true.
                endif
              else
                call die('read_Zmatrix: Invalid symbol in Z-matrix')
              endif

              i = 1
              found = .false.
              do while (i.le.nVars.and..not.found)
                found = leqi(ZmatVarNames(i),fdf_bnames(pline,2))
                i = i + 1
              enddo
              if (found) then
                if (.not.lZmatVarsDef(i-1)) then
                  call die('read_Zmatrix: Invalid dependency for
     .                      Z-matrix symbol')
                else
                  indi1 = iVarstoZmat(i-1)
                endif
              else
                call die('read_Zmatrix: Invalid symbol in Z-matrix')
              endif
C If coefficient A=0 then its really a constant definition
              if (fdf_breals(pline,1).eq.0.0d0) then
                call die('read_Zmatrix: erroneous constraints
     .                      definition - rather put in constants
     .                      definition ')
              endif

C Put values into the newly defined Zmat components
              j = indi
              do while (j.ne.0) 
                Zmat(j) = fdf_breals(pline,1)*Zmat(indi1) +
     .                    fdf_breals(pline,2)
                j = iNextDept(j)
              enddo
C Decide which is the independent variable
              if (indi.gt.indi1) then
                A = fdf_breals(pline,1)
                B = fdf_breals(pline,2)
              else
                A = 1.0d0/fdf_breals(pline,1)
                B = -fdf_breals(pline,2)/fdf_breals(pline,1)
                j = indi
                indi = indi1
                indi1 = j
                VaryZmat(indi1) = VaryZmat(indi)
              endif
              VaryZmat(indi) = .false.
C Now indi1 is the independent var with indi = A*indi1+B
C Link the dependency lists of the two variables and update
C   coeffiecients and iZmattoVars array
              j = indi1
              do while (iNextDept(j).ne.0) 
                j = iNextDept(j)
              enddo
              iNextDept(j) = indi
              j = indi
              do while (j.ne.0)
                iZmattoVars(j) = iZmattoVars(indi1)
                coeffA(j) = A*coeffA(indi1)
                coeffB(j) = B + A*coeffB(indi1)
                j = iNextDept(j)
              enddo

            endif

          endif
C End of parsing input line            

        enddo
  10    continue
C Done reading Z-matrix
        if (Node.eq.0) write(*,"(a)") "%endblock Zmatrix"

        species(1:nAtoms) = nspecies(1:nAtoms)

C Check that all variables have been defined
        do i = 1,nVars
          if (.not.lZmatVarsDef(i)) then
            call die('read_Zmatrix: Symbol has not been given a 
     .                value')
          endif
        enddo

C Scale all coordinates according to their type            
        do i = 1,3*nAtoms
C Angle
          if (ZmatType(i).eq.1) then
            Zmat(i) = Zmat(i) * (1.0d0+(deg-1.0d0)*ZmatUnitsAng)
            coeffB(i) = coeffB(i) * (1.0d0+(deg-1.0d0)*ZmatUnitsAng)
          endif

C Bond length/pure cartesian
          if (ZmatType(i).eq.2.or.ZmatType(i).eq.3.or.
     .        ZmatType(i).eq.6) then
            Zmat(i) = Zmat(i) * (1.0d0+(Ang-1.0d0)*ZmatUnitsLen)
            coeffB(i) = coeffB(i) * (1.0d0+(Ang-1.0d0)*ZmatUnitsLen)
          endif
             
C Scaled cartesian
          if (ZmatType(i).eq.4.or.ZmatType(i).eq.7) then
            Zmat(i) = Zmat(i) * alat
            coeffB(i) = coeffB(i) * alat
          endif
                    
        enddo

C Translate pure cartesians by the origin
        do i=1,nAtoms
          if ((ZmatType(3*i).eq.3.or.ZmatType(3*i).eq.6)
     .                  .and.lOrigin) then
            Zmat(3*i-2) = Zmat(3*i-2) + origin(1)
            Zmat(3*i-1) = Zmat(3*i-1) + origin(2)
            Zmat(3*i) = Zmat(3*i) + origin(3)
          endif
        enddo

C Fractional cartesian 
        do i=1,nAtoms
          if (ZmatType(3*i).eq.5.or.ZmatType(3*i).eq.8) then
            do k=1,3
              ZmatVal(k) = Zmat(3*(i-1)+k)
              coeffBVal(k) = coeffB(3*(i-1)+k)
            enddo
            do k=1,3
              Zmat(3*(i-1)+k) = unit_cell(k,1)*ZmatVal(1) +
     .                          unit_cell(k,2)*ZmatVal(2) +
     .                          unit_cell(k,3)*ZmatVal(3)
              coeffB(3*(i-1)+k) = unit_cell(k,1)*coeffBVal(1) +
     .                          unit_cell(k,2)*coeffBVal(2) +
     .                          unit_cell(k,3)*coeffBVal(3)
            enddo
          endif
        enddo

C Set remaining undefined dependencies to 1 (they are not used but if they 
C   are 0 the program might crash later)
        do m=1,nZmol
          nStart = nZmolStartAtom(m)
          iZmat(3*nStart-2) = 1
          iZmat(3*nStart-1) = 1
          iZmat(3*nStart) = 1
          iZmat(3*nStart+2) = 1
          iZmat(3*nStart+3) = 1
          iZmat(3*nStart+6) = 1
        enddo



C Check Z matrix data for consistency
C   - are all atoms defined in terms of previous atoms?
C   - are the three dependencies for every atom distinct?

C Loop over molecules        
        do m = 1,nZmol
          nStart = nZmolStartAtom(m)
          nAtoms = nZmolAtoms(m)

C PB BEGIN
C         Check whether any atom info has been read.
          if (nAtoms .eq. 0) then
             if (Node .eq. 0) then
                write (6,'(A)') 'read_Zmatrix: Ill defined Zmatrix:'
                write (6,'(A,i3)')
     .               'read_Zmatrix: no atoms defined at molecule', m
             endif
             call die()
          endif
C PB END

C Second atom
          if (nAtoms.gt.1) then
            if (iZmat(3*nStart+1).ne.nStart) then
              if (Node.eq.0) then
                write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
                write(6,'(''read_Zmatrix: molecule nr '',i7,
     .                    ''; atom nr 2'')') m
              endif
              call die()
            endif
          endif

C Third atom
          if (nAtoms.gt.2) then
            if (.not.(iZmat(3*nStart+4).eq.nStart.and.iZmat(3*nStart+5)
     .          .eq.nStart+1.or.iZmat(3*nStart+5).eq.nStart.and.
     .          iZmat(3*nStart+4).eq.nStart+1)) then
              if (Node.eq.0) then
                write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
                write(6,'(''read_Zmatrix: molecule nr '',i7,
     .                    ''; atom nr 3'')') m
              endif
              call die ()
            endif
          endif

C Fourth atom and up -> general case          
          do i = nStart+3,nStart+nAtoms-1
            if (iZmat(3*i-2).gt.i-1.or.
     .          iZmat(3*i-2).eq.iZmat(3*i-1)) then
              if (Node.eq.0) then
                write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
                write(6,'(''read_Zmatrix: molecule nr '',i7,
     .                    ''; atom nr '',i7)') m,i-nStart
              endif
              call die()
            endif
            if (iZmat(3*i-1).gt.i-1.or.
     .          iZmat(3*i-1).eq.iZmat(3*i)) then
              if (Node.eq.0) then
                write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
                write(6,'(''read_Zmatrix: molecule nr '',i7,
     .                    ''; atom nr '',i7)') m,i-nStart
              endif
              call die()
            endif
            if (iZmat(3*i).gt.i-1.or.
     .          iZmat(3*i).eq.iZmat(3*i-2)) then
              if (Node.eq.0) then
                write(6,'(''read_Zmatrix: Ill defined Zmatrix:'')')
                write(6,'(''read_Zmatrix: molecule nr '',i7,
     .                    ''; atom nr '',i7)') m,i-nStart
              endif
              call die()
            endif
          enddo
C End looping over molecules          
        enddo

        if (ZmatUnitsAng == 0) then
          scale_angle = 1.0_dp
        else
          scale_angle = deg
        endif

        if (ZmatUnitsLen == 0) then
          scale_length = 1.0_dp
        else
          scale_length = Ang
        endif

C Output Z-matrix coordinates
        if (Node.eq.0) then
           print *, "nZmol: ", nZmol
           print *, "nSpecies: ", nSpecies
          call write_Zmatrix()
        endif

C Return point
  999   continue
        firsttime = .false.

        end subroutine read_Zmatrix

        subroutine Zmat_to_Cartesian(Cartesian)
C
C  Passed variables
C
        real(dp), intent(inout) ::  Cartesian(:,:)
C        
C  Local Variables
C
        integer               :: i
        integer               :: m
        integer               :: nStart
        integer               :: nAtoms
        real(dp)              :: RelatedC(3,3)
        real(dp)              :: x
        real(dp)              :: y
        real(dp)              :: z

        integer               :: index_i, index_j
        real(dp)              :: phi_ref

C  Loop over molecules

        do m = 1,nZmol

          nStart = nZmolStartAtom(m)
          nAtoms = nZmolAtoms(m)


C  Loop over atoms within each molecule generating Cartesian coordinates
C  Use general Z2CGen subroutine which takes special first 3 atoms into account
          do i = nStart,nStart+nAtoms-1

            ! reference atoms in the order i,j,k

            RelatedC(1:3,1) = Cartesian(1:3,iZmat(3*i-2))
            RelatedC(1:3,2) = Cartesian(1:3,iZmat(3*i-1))
            RelatedC(1:3,3) = Cartesian(1:3,iZmat(3*i))

            ! Pass explicitly the azimuth angle
            ! of the atom "i" to which this one 'binds' with respect
            ! to j
            index_i = iZmat(3*i-2)
            if (ZmatType(3*index_i) == 1) then  ! an angle
               phi_ref = Zmat(3*index_i)
            else
               ! i is the first atom in the molecule...
               ! We pass then phi_ji - pi, where phi_ji
               ! is the azimuth of j with respect to i.
               !
               index_j = iZmat(3*i-1)
               phi_ref = Zmat(3*index_j) - pi
            endif

            call Z2CGen(i-nStart,Zmat(3*i-2),Zmat(3*i-1),Zmat(3*i),
     .                      RelatedC,x,y,z, phi_ref)
            Cartesian(1,i) = x
            Cartesian(2,i) = y
            Cartesian(3,i) = z
          enddo

        enddo

C  Loop over cartesian blocks
        do m = 1,nZcart

          nStart = nZcartStartAtom(m)
          nAtoms = nZcartAtoms(m)

          do i = nStart,nStart+nAtoms-1
            Cartesian(1,i) = Zmat(3*i-2)
            Cartesian(2,i) = Zmat(3*i-1)
            Cartesian(3,i) = Zmat(3*i)
          enddo

        enddo

        end subroutine Zmat_to_Cartesian

        subroutine iofaZmat()
        
        implicit none
C
C  Local variables
C
        integer               :: m
        integer               :: i
        integer               :: k
        integer               :: nStart
        integer               :: nAtoms
        integer               :: jindex
        character(len=4)      :: lenstr
        character(len=4)      :: angstr
C
        if (ZmatUnitsLen.eq.0) then
          lenstr = 'Bohr'
        elseif (ZmatUnitsLen.eq.1) then
          lenstr = 'Ang'
        endif
        if (ZmatUnitsAng.eq.0) then
          angstr = 'rad'
        elseif (ZmatUnitsAng.eq.1) then
          angstr = 'deg'
        endif

        write(6,'(/,''zmatrix: Atomic forces (eV/'',a4,
     .            ''; eV/'',a4,'')'')') lenstr,angstr
        write(6,'(''zmatrix: (No information if symbols are used)'')')
C  Print molecule forces in the user's units
        do m = 1,nZmol
          nStart = nZmolStartAtom(m)
          nAtoms = nZmolAtoms(m)
          write(6,'(''molecule'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
          write(6,'(i6,3f12.6)') nStart,
     .          (ZmatForce(3*(nStart-1)+k)*(1+(Ang-1)*ZmatUnitsLen)/eV
     .          ,k=1,3)
          write(6,'(i6,3f12.6)') (i,ZmatForce(3*i-2)*
     .                            (1+(Ang-1)*ZmatUnitsLen)/eV,
     .             ZmatForce(3*i-1)*(1+(deg-1)*ZmatUnitsAng)/eV,
     .             ZmatForce(3*i)*(1+(deg-1)*ZmatUnitsAng)/eV,
     .             i=nStart+1,nStart+nAtoms-1)
        enddo

C  Print cartesian forces in the user's units
        do m = 1,nZcart
          nStart = nZcartStartAtom(m)
          nAtoms = nZcartAtoms(m)
          write(6,'(''cartesian'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
          write(6,'(i6,3f12.6)') (i,
     .          (ZmatForce(3*(i-1)+k)*(1+(Ang-1)*ZmatUnitsLen)/eV
     .          ,k=1,3),i=nStart,nStart+nAtoms-1)
        enddo
        write (6,'(/)') 

C  Print variables' forces in user's units
        write(6,'(''Variable forces (eV/'',a4,
     .            ''; eV/'',a4,'')'')') lenstr,angstr
        do i = 1,nVars
           jindex = iVarstoZmat(i)
           if (.not.VaryZmat(jindex)) cycle

          if (ZmatType(iVarstoZmat(i)).eq.1) then
            write(6,*) ZmatVarNames(i), ZmatForceVar(i)*
     .                      (1+(deg-1)*ZmatUnitsAng)/eV
          else 
            write(6,*) ZmatVarNames(i), ZmatForceVar(i)*
     .                      (1+(Ang-1)*ZmatUnitsLen)/eV
          endif
        enddo

        return

        end subroutine iofaZmat

        subroutine CartesianForce_to_ZmatForce(na,Cartesian,
     .                                         CartesianForce)

        use parallel,    only : IONode
#ifdef MPI
        use mpi_siesta
#endif

        implicit none

C
C  Passed variables
C
        integer,  intent(in)  :: na
        real(dp), intent(in)  :: Cartesian(3,na)
        real(dp), intent(in)  :: CartesianForce(3,na)
C
C  Local variables
C
        integer                     :: i
        integer                     :: j
        integer                     :: k
        integer                     :: ia
        integer                     :: vi
        real(dp), parameter         :: hfrac = 1e-6
        real(dp), pointer           :: CartesianB(:,:)=>null(),
     $                                 CartesianF(:,:)=>null()
        real(dp), pointer           :: ZmatF(:)=>null(),
     $                                 ZmatB(:)=>null()
        real(dp)                    :: h
#ifdef MPI
        integer                     :: MPIerror
#endif        

      if (IOnode) then
C       Allocate local workspace arrays
        call re_alloc( CartesianF, 1, 3, 1, na, 'CartesianF',
     &                 'zmatrix' )
        call re_alloc( CartesianB, 1, 3, 1, na, 'CartesianB',
     &                 'zmatrix' )
        call re_alloc( ZmatF, 1, 3*na, 'ZmatF', 'zmatrix' )
        call re_alloc( ZmatB, 1, 3*na, 'ZmatB', 'zmatrix' )

C  Initialise forces & cartesian transformation matrices
        ZmatForce(1:3*na) = 0.0_dp
        ZmatForceVar(1:nVars) = 0.0_dp
        CartesianF(1:3,1:na) = 0.0_dp
        CartesianB(1:3,1:na) = 0.0_dp

C  Reset the FD and BD Z-Zmatrices
        ZmatF(1:3*na) = Zmat(1:3*na)
        ZmatB(1:3*na) = Zmat(1:3*na)

C Zero the Forward & Backward Cartesian matrices:
        CartesianF = 0.0_dp
        CartesianB = 0.0_dp

C  Loop over all atoms
        do ia = 1,na
C  Loop over all coordinates of the atom
          do k = 1,3
C  Index for this coordinate          
            i = 3*(ia-1) + k
C  Only proceed if this coordinate is allowed to vary or CalcAllForces is true            
            if (VaryZmat(i).or.lCalcAllForces) then
C  Set the h-increment for the coordinate
              if (Zmat(i).eq.0.0_dp) then 
                h = hfrac
              else 
                h = hfrac*Zmat(i)
              endif
C  Update the FD and BD Zmatrices for this coordinate
              ZmatF(i) = Zmat(i) + h
              ZmatB(i) = Zmat(i) - h
C  Calculate the force on this coordinate
              if (iZmattoVars(i).eq.0.or.lCalcAllForces) then
                ZmatForce(i) = GetForce(na,Cartesian,CartesianForce,
     .                                  CartesianF,CartesianB,
     .                                  ZmatF,ZmatB,
     .                                  i,h,.false.)
              endif
C  Check if variable is independent              
              vi = iZmattoVars(i)
              j = 0
              if (vi.ne.0) then
                j = iVarstoZmat(vi)
              endif
C  Calculate the force on this variable
              if (j.eq.i) then
C  Update the FD and BD Z-matrices for all dependents
                do while (iNextDept(j).ne.0) 
                  j = iNextDept(j)
                  ZmatF(j) = Zmat(j) + h*coeffA(j)
                  ZmatB(j) = Zmat(j) - h*coeffA(j)
                enddo
                ZmatForceVar(vi) = GetForce(na,Cartesian,CartesianForce,
     .                                  CartesianF,CartesianB,
     .                                  ZmatF,ZmatB,
     .                                  i,h,.true.)
C  reset the FD and BD Z-matrices
                j = i
                do while (iNextDept(j).ne.0)
                  j = iNextDept(j)
                  ZmatF(j) = Zmat(j)
                  ZmatB(j) = Zmat(j)
                enddo
              endif
              ZmatF(i) = Zmat(i)
              ZmatB(i) = Zmat(i)
            endif
          enddo
C  Done looping over 3 coordinate of this atom          
C  From now on the previous atoms will always be unchanged
          CartesianF(1:3,ia) = Cartesian(1:3,ia)
          CartesianB(1:3,ia) = Cartesian(1:3,ia)
        enddo
C  Done looping over all atoms

C  Deallocate workspace arrays
        call de_alloc( ZmatB, 'ZmatB', 'zmatrix' )
        call de_alloc( ZmatF, 'ZmatF', 'zmatrix' )
        call de_alloc( CartesianB, 'CartesianB', 'zmatrix' )
        call de_alloc( CartesianF, 'CartesianF', 'zmatrix' )
      endif !IOnode

#ifdef MPI
        call MPI_Bcast(ZmatForce,3*na,MPI_double_precision,0,
     .                 MPI_Comm_World,MPIerror)
        call MPI_Bcast(ZmatForceVar,nVars,MPI_double_precision,0,
     .                 MPI_Comm_World,MPIerror)
#endif        

      end subroutine CartesianForce_to_ZmatForce

      function GetForce(na,Cartesian,CartesianForce,
     .                  CartesianF,CartesianB,
     .                  ZmatF,ZmatB,i,h,doAll)

      implicit none

C  Passed variables
      integer, intent(in)           :: i
      logical, intent(in)           :: doAll
      integer,  intent(in)          :: na
      real(dp), intent(in)          :: Cartesian(3,na)
      real(dp), intent(in)          :: CartesianForce(3,na)
      real(dp), intent(inout)       :: CartesianF(3,na)
      real(dp), intent(inout)       :: CartesianB(3,na)
      real(dp), intent(inout)       :: ZmatF(3*na)
      real(dp), intent(inout)       :: ZmatB(3*na)
      real(dp), intent(in)          :: h


C  Local variables
      real(dp)        :: GetForce
      integer         :: m
      integer         :: prevm
      integer         :: j
      integer         :: l
      integer         :: ca
      real(dp)        :: RelatedC(3,3)
      real(dp)        :: X,Y,Z

      integer               :: index_i, index_j
      real(dp)              :: phi_ref_F, phi_ref_B

C  Initialize force
      GetForce = 0.0d0
C  Add force contributions from all the affected cartesian coordinates 
      m = 1
      prevm = 0
C  Loop over affected Zmat Coords            
      j = i
      do while (j.ne.0) 
C  l is the atom number for index j            
        l = (j+2)/3
        if (ZmatType(j).lt.6) then
C  Molecule              
C  Find correct molecule
          do while (l.gt.nZmolStartAtom(m)+nZmolAtoms(m)-1)
            m = m + 1
          enddo
C  Only add force contributions if this is the first affected 
C     coordinate for a specific molecule (otherwise it has been done already)
          if (m.ne.prevm) then
            prevm = m
C  Initialise FD and BD cartesians that come before changed atom
            do ca = nZmolStartAtom(m),l-1
              CartesianF(1:3,ca) = Cartesian(1:3,ca)
              CartesianB(1:3,ca) = Cartesian(1:3,ca)
            enddo
C  Add force contributions from all atoms later than and
C     including the changed one in this molecule
            do ca = l,nZmolStartAtom(m)+nZmolAtoms(m)-1
C  Forward cartesian     
           
              RelatedC(1:3,1) = CartesianF(1:3,iZmat(3*ca-2))
              RelatedC(1:3,2) = CartesianF(1:3,iZmat(3*ca-1))
              RelatedC(1:3,3) = CartesianF(1:3,iZmat(3*ca))

            ! Pass explicitly the azimuth angle
            ! of the atom "i" to which this one 'binds' with respect
            ! to j
              index_i = iZmat(3*ca-2)
              if (ZmatType(3*index_i) == 1) then ! an angle
                 phi_ref_F = ZmatF(3*index_i)
                 phi_ref_B = ZmatB(3*index_i)
              else
               ! i is the first atom in the molecule...
               ! We pass then phi_ji - pi, where phi_ji
               ! is the azimuth of j with respect to i.
               !
               index_j = iZmat(3*ca-1)
               phi_ref_F = ZmatF(3*index_j) - pi
               phi_ref_B = ZmatB(3*index_j) - pi
            endif

              call Z2CGen(ca-nZmolStartAtom(m),ZmatF(3*ca-2),
     .                    ZmatF(3*ca-1),ZmatF(3*ca),RelatedC,
     .                    CartesianF(1,ca),CartesianF(2,ca),
     .                    CartesianF(3,ca), phi_ref_F)

C  Backward cartesian               

              RelatedC(1:3,1) = CartesianB(1:3,iZmat(3*ca-2))
              RelatedC(1:3,2) = CartesianB(1:3,iZmat(3*ca-1))
              RelatedC(1:3,3) = CartesianB(1:3,iZmat(3*ca))

              call Z2CGen(ca-nZmolStartAtom(m),ZmatB(3*ca-2),
     .                    ZmatB(3*ca-1),ZmatB(3*ca),RelatedC,
     .                    CartesianB(1,ca),CartesianB(2,ca),
     .                    CartesianB(3,ca), phi_ref_B)

C  Add in force constribution     
              X = CartesianForce(1,ca)*
     .            (CartesianF(1,ca)-CartesianB(1,ca))/(2*h)
              Y = CartesianForce(2,ca)*
     .            (CartesianF(2,ca)-CartesianB(2,ca))/(2*h)
              Z = CartesianForce(3,ca)*
     .            (CartesianF(3,ca)-CartesianB(3,ca))/(2*h)
              GetForce = GetForce + X+Y+Z
C testing
C             write(6,*) 'test:',i,ca,X,Y,Z
            enddo
            endif
        else  
C  Cartesian              
          GetForce = GetForce+CartesianForce(j-3*(l-1),l)*coeffA(j)
        endif  
        if (doAll) then
          j = iNextDept(j)
        else 
          j=0
        endif
      enddo

      end function GetForce


      subroutine write_canonical_Zmatrix(filename,iunit)
!
!     Write a "canonical" Zmatrix with the following characteristics:
!     1. No symbolic variables or constants are used.
!     2. The position coordinates of the first atom in each molecule
!        are absolute cartesian coordinates.
!     3. Any "crystal" (non-molecule) coordinates are also absolute cartesians.
!     4. There is no provision for output of constraints.
!     5. The units used are those initially specified by the user, and are
!        specified in fdf form.
!
!     The output is to file "OUT.CANON.ZMATRIX" by default, but a new filename
!     might be specified as an optional argument. If a unit number is used
!     output is directly appended to it.
!     
        implicit none
        
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: iunit
C
C Local variables
C
        integer             :: nStart
        integer             :: nAtoms
        integer             :: k
        integer             :: i, iu, js
        integer             :: m
        character(len=4)    :: lenstr
        character(len=4)    :: angstr
        character(len=90)    :: fname

        integer, dimension(3)             :: mask

C  Output Z matrix information
        if (ZmatUnitsLen.eq.0) then
          lenstr = 'Bohr'
        elseif (ZmatUnitsLen.eq.1) then
          lenstr = 'Ang'
        endif
        if (ZmatUnitsAng.eq.0) then
          angstr = 'rad'
        elseif (ZmatUnitsAng.eq.1) then
          angstr = 'deg'
        endif

      if (present(iunit)) then
         iu = iunit
      else
         if (present(filename)) then
            fname = filename
         else
            fname = "OUT.CANON.ZMATRIX"
         endif
         call io_assign( iu )
         open(iu, file=trim(fname), form='formatted',
     $        position='rewind', status='unknown')
      endif

C  Write molecule coordinates in user's units
      write(iu,"(2a)") "ZM.UnitsLength ", trim(lenstr)
      write(iu,"(2a)") "ZM.UnitsAngle ", trim(angstr)
      write(iu,"(a)") "%block Zmatrix"
        do m = 1,nZmol
          write(iu,"(a)") "molecule_cartesian"
          nStart = nZmolStartAtom(m)
          nAtoms = nZmolAtoms(m)
          js = 3*(nStart-1)
          ! First atom
          write(iu,'(i2,3i4,3f16.8,3i4)') species(nStart),
     $          0, 0, 0,
     $         (Zmat(js+k)/scale_length,k=1,3),
     $         (to_int(VaryZmat(js+k)),k=1,3)

          ! Start this loop with the second atom in the molecule
          ! Mask the references for the second and third atom
          do i=nStart+1,nStart+nAtoms-1
             js = 3*i-3
             select case (i-nStart+1) 
             case (2) 
                mask = (/ 1, 0, 0 /)
             case (3)
                mask = (/ 1, 1, 0 /)
             case (4:)
                mask = (/ 1, 1, 1 /)
             end select
             ! Use an "intra-molecule" numbering scheme for
             ! the reference atoms
             write(iu,'(i2,3i4,3f16.8,3i4)') species(i),
     $         (iZmat(js+1:js+3)-nStart+1)*mask,
     .            Zmat(3*i-2)/scale_length,
     .            Zmat(3*i-1)/scale_angle,
     .            Zmat(3*i)/scale_angle,
     $            (to_int(VaryZmat(js+k)),k=1,3)

          enddo
        enddo

C  Write cartesian coordinates in user's units        
        do m = 1,nZcart
          nStart = nZcartStartAtom(m)
          nAtoms = nZcartAtoms(m)
          write(iu,"(a)") "cartesian"
          do i=nStart,nStart+nAtoms-1
             js = 3*(i-1)
             write(iu,'(i2,3f16.8,3i4)') species(i),
     $                            (Zmat(js+k)/scale_length,k=1,3),
     $                            (to_int(VaryZmat(js+k)),k=1,3)
          enddo

        enddo
        write(iu,"(a)") "%endblock Zmatrix"
        if (.not. present(iunit)) then
           call io_close(iu)
        endif

      end subroutine write_canonical_Zmatrix

      subroutine write_Zmatrix

        implicit none
C
C Local variables
C
        integer             :: nStart
        integer             :: nAtoms
        integer             :: k
        integer             :: i
        integer             :: m
        character(len=4)    :: lenstr
        character(len=4)    :: angstr

C  Output Z matrix information
        if (ZmatUnitsLen.eq.0) then
          lenstr = 'Bohr'
        elseif (ZmatUnitsLen.eq.1) then
          lenstr = 'Ang'
        endif
        if (ZmatUnitsAng.eq.0) then
          angstr = 'rad'
        elseif (ZmatUnitsAng.eq.1) then
          angstr = 'deg'
        endif


C  Write molecule coordinates in user's units
        write(6,'(/,''zmatrix: Z-matrix coordinates: ('',a4,
     .            ''; '',a4,'')'')') lenstr,angstr
        write(6,'(''zmatrix: '',
     $              ''(Fractional coordinates have been converted '',
     $              ''to cartesian)'')')
        do m = 1,nZmol
          nStart = nZmolStartAtom(m)
          nAtoms = nZmolAtoms(m)
          write(6,'(''molecule'',i5,'' ('',i6,'' atoms)'')') m,nAtoms
          write(6,'(3f16.8)') (Zmat(3*(nStart-1)+k)/(1.0d0+(Ang-1.0d0)*
     .                         ZmatUnitsLen),k=1,3)
          do i=nStart+1,nStart+nAtoms-1
            write(6,'(3f16.8)') 
     .            Zmat(3*i-2)/(1.0d0+(Ang-1.0d0)*ZmatUnitsLen),
     .            Zmat(3*i-1)/(1.0d0+(deg-1.0d0)*ZmatUnitsAng),
     .            Zmat(3*i)/(1.0d0+(deg-1.0d0)*ZmatUnitsAng)
          enddo
        enddo

C  Write cartesian coordinates in user's units        
        do m = 1,nZcart
          nStart = nZcartStartAtom(m)
          nAtoms = nZcartAtoms(m)
          write(6,'(''cartesian block'',i5,'' ('',i6,'' atoms)'')') 
     .              m,nAtoms
          write(6,'(3f16.8)')((Zmat(3*(i-1)+k)/(1+(Ang-1)*ZmatUnitsLen),
     .                         k=1,3),i=nStart,nStart+nAtoms-1)
        enddo
        write(6,'(/)')

        call print_variables()

      end subroutine write_Zmatrix

      subroutine Z2CGen(atomNr,r,theta,phi,RelatedC,x,y,z,phi_ref)
      use parallel, only: IONode
C
C  Subroutine for generation of cartesian coordinates from Z-matrix
C  given the following special cases:
C    - the first entry is pure cartesian coordinates
C    - the second entry is spherical coordinates relative to the first atom
C    - the third entry is general Z-matrix coordinates, but with the torsion atom 
C      a dummy atom 1 unit in the z direction above the second
C    - the fourth entry onwards are general Z-matrix coordinates
C
        implicit none
C
C  Passed variables
C
        integer,  intent(in)     :: atomNr
        real(dp), intent(in)     :: r
        real(dp), intent(in)     :: theta
        real(dp), intent(in)     :: phi
        real(dp), intent(inout)  :: RelatedC(3,3)
        real(dp), intent(out)    :: x
        real(dp), intent(out)    :: y
        real(dp), intent(out)    :: z
!
        real(dp), intent(in)     :: phi_ref
C        
        real(dp), dimension(3), parameter ::
     $        ex = (/ 1.0_dp, 0.0_dp, 0.0_dp /),
     $        ez = (/ 0.0_dp, 0.0_dp, 1.0_dp /)

        real(dp), dimension(3)   :: third_atom_coords, tac2
        real(dp), dimension(3)   :: unit_vector, r2
        real(dp) :: r2_mod,  angle_2

        real(dp) :: r_1, phi_1, theta_1
        real(dp) :: x_0, y_0, z_0

        if (atomNr.eq.0) then
          x = r
          y = theta
          z = phi
        elseif (atomNr.eq.1) then

          x = r*sin(theta)*cos(phi) + RelatedC(1,1)
          y = r*sin(theta)*sin(phi) + RelatedC(2,1)
          z = r*cos(theta) + RelatedC(3,1)

        elseif (atomNr.eq.2) then

                x_0 = RelatedC(1,2)
                y_0 = RelatedC(2,2)
                z_0 = RelatedC(3,2)
!               Note, the atom this atom binds to and the previous one
                r2 = RelatedC(:,1) - RelatedC(:,2)

                r2_mod = sqrt(dot_product(r2,r2))
                theta_1 = acos(r2(3)/r2_mod)        ! Could this be not general enough?
                ! Theta should be always between 0 and pi.
                ! If it is not, it should be forced to be, changing its sign and
                ! changing phi to phi+/- pi
                ! This is not currently enforced!

                if (
     $               (abs(r2(2)) .gt. 1.0e-8_dp)
     $               .and.
     $               (abs(r2(1)) .gt. 1.0e-8_dp)
     $             ) then
                   ! Compute explicitly the phi angle

                   phi_1 = atan2(r2(2),r2(1))

                else
                   ! Take it from the Zmatrix itself
                   phi_1 = phi_ref
                endif

                r_1 = r2_mod

          ! Change temporarily the  coordinates of the reference
          ! atoms for the third atom to:

          ! Atom2 along the x axis, at a distance r_1 from the origin.
          ! Remember that RelatedC(:,1) refers to the atom to which this
          ! one "binds" (the "i" in the "i j k" in the Zmatrix spec)
          RelatedC(:,1) = r_1 * ex

          ! Atom 1 at the origin
          RelatedC(:,2) = 0.0_dp  ! Note the index: atom1 is now j

          ! Fake reference atom at (0,0,1) over the second 
          ! Atom k in the Zmatrix spec
          RelatedC(:,3) = RelatedC(:,1)  + ez

          call Z2C(r,theta,phi,RelatedC,x,y,z)
!
! We need to rotate and translate the third atom

           ! Original (reduced) positions

            third_atom_coords(1:3) = (/ x, y, z /)  !just computed above

          ! First, rotate by phi_1 around the z axis.
          ! Use Rodrigues formula (see 'rotation matrix' in Wikipedia)
            third_atom_coords(1:3) =
     $           rotate(third_atom_coords,ez,phi_1)

            ! New rotation axis defined by a vector lying on the xy
            ! plane, forming a phi_1-90 degree angle with the x axis:

            unit_vector = rotate(ex,ez, phi_1-pi/2.0_dp)

            ! Rotate 
            ! with the unit_vector as axis and an angle of 90-theta_1

            angle_2 = pi/2.0_dp - theta_1
            tac2(:) = rotate(third_atom_coords,unit_vector,angle_2)
!
!           Now translate

            third_atom_coords = tac2 +  (/ x_0, y_0, z_0 /)
!
            x = third_atom_coords(1)
            y = third_atom_coords(2)
            z = third_atom_coords(3)

            
        elseif(atomNr.gt.2) then
          call Z2C(r,theta,phi,RelatedC,x,y,z)
        endif

      end subroutine Z2CGen

      subroutine Z2C(r,theta,phi,RelatedC,x,y,z)
C
C  Subroutine for generation of Cartesian coordinates from Z-Matrix
C
C  Julian Gale, NRI, Curtin University, March 2004
C
        implicit none
C
C  Passed variables
C
        real(dp), intent(in)  :: r
        real(dp), intent(in)  :: theta
        real(dp), intent(in)  :: phi
        real(dp), intent(in)  :: RelatedC(3,3)
        real(dp), intent(out) :: x
        real(dp), intent(out) :: y
        real(dp), intent(out) :: z
C
C  Local variables
C
        real(dp)          :: rji
        real(dp)          :: rn
        real(dp)          :: rp
        real(dp)          :: xi
        real(dp)          :: yi
        real(dp)          :: zi
        real(dp)          :: xj
        real(dp)          :: yj
        real(dp)          :: zj
        real(dp)          :: xji
        real(dp)          :: yji
        real(dp)          :: zji
        real(dp)          :: xk
        real(dp)          :: yk
        real(dp)          :: zk
        real(dp)          :: xki
        real(dp)          :: yki
        real(dp)          :: zki
        real(dp)          :: xn
        real(dp)          :: yn
        real(dp)          :: zn
        real(dp)          :: xp
        real(dp)          :: yp
        real(dp)          :: zp
C
C  Find coordinates for related atoms
C
        xi = RelatedC(1,1)
        yi = RelatedC(2,1)
        zi = RelatedC(3,1)
!       vec_i = RelatedC(1:3,1)
        xj = RelatedC(1,2)
        yj = RelatedC(2,2)
        zj = RelatedC(3,2)
!       vec_j = RelatedC(1:3,2)
        xk = RelatedC(1,3)
        yk = RelatedC(2,3)
        zk = RelatedC(3,3)
!       vec_k = RelatedC(1:3,3)
C
C  Find unit vector along j->i vector
C
        xji = xi - xj
        yji = yi - yj
        zji = zi - zj
!       vec_ji = vec_i - vec_j
        rji = xji*xji + yji*yji + zji*zji
        rji = sqrt(rji)
        xji = xji/rji
        yji = yji/rji
        zji = zji/rji
!       call normalize(vec_ji)
C
C  Find j->k vector
C
        xki = xk - xj
        yki = yk - yj
        zki = zk - zj
!       vec_kj = vec_k - vec_j
C
C  Find unit vector normal to the i-j-k plane
C
        xn = yji*zki - yki*zji
        yn = zji*xki - zki*xji
        zn = xji*yki - xki*yji
!       vec_n = cross_product(vec_ji,vec_kj)
        rn = xn*xn + yn*yn + zn*zn
        rn = sqrt(rn)
        xn = xn/rn
        yn = yn/rn
        zn = zn/rn
!       call normalize(vec_n)
C
C  Find unit vector normal to the other 2 directions already found
C
C  Since original vectors are normalised the result should be likewise
C
        xp = yn*zji - yji*zn
        yp = zn*xji - zji*xn
        zp = xn*yji - xji*yn
!       vec_p = cross_product(vec_n,vec_ji)
C
C  Find distances along each unit vector
C
        rji = r*cos(theta)
        rn  = r*sin(theta)*sin(phi)
        rp  = r*sin(theta)*cos(phi)
C
C  Multiply unit vector by distances and add to origin to get position
C
!       rji enters with a minus sign as the convention for the angle is
!       that i is at the vertex, so an angle of less than pi/2 has the
!       atom "looking back" at j, and so opposite to vec_ji.
!
        x = xi - rji*xji + rn*xn + rp*xp
        y = yi - rji*yji + rn*yn + rp*yp
        z = zi - rji*zji + rn*zn + rp*zp
!       vec_pos = vec_i - rji*vec_ji + rn*vec_n + rp*vec_p
!       x = vec_pos(1)
!       y = vec_pos(2)
!       z = vec_pos(3)
C
      end subroutine z2c

!-------------------------------------------------------------

      function user_print(i,offset) result(x)
!
!     Converts variables and constants to the
!     format originally supplied by the user.
!
      use m_cell, only  : cart2frac

      ! Index of Zmatrix variable
      integer, intent(in)  :: i
      real(dp), intent(in), optional :: offset
      real(dp)             :: x  ! final value

      real(dp) :: xin
      real(dp) :: r(3), rfrac(3)
      integer  :: natom, k
!
!     Offset: prepared to deal with origin shift... *****
!
C  ZmatType(3*na) defines the component type:
C    1 - angle
C    2 - bond length
C    3 - pure cartesian (molecule)
C    4 - scaled cartesian (molecule)
C    5 - fractional cartesian (molecule)
C    6 - pure cartesian
C    7 - scaled cartesian (scaled by lattice constant)
C    8 - fractional cartesian (scaled by lattice vectors)

      xin = Zmat(i)
      if (present(offset)) then
         xin = xin + offset
      endif

      select case (ZmatType(i))

      case(1)
         x =  xin/scale_angle
      case(2,3,6)
         x=   xin/scale_length
      case(4,7)
         x =  xin / zmatrix_alat
      case(5,8)
   
         ! Get other cart coords and compute fractional coords
         ! we assume that any use of fractional coordinates is
         ! on a whole-atom basis, i.e., there are no bond lenghts
         ! or angles involved, so Zmat(base+1:base+3) is homogeneous
   
         natom = 1 + (i-1)/3
         r(1:3) = Zmat(3*(natom-1)+1:3*(natom-1)+3)
         call cart2frac(r,rfrac)
         k = i - 3*(natom-1)
         x = rfrac(k)

      case default
         call die("Wrong type for Zmatrix coordinate")
      end select

      end function user_print

!----------------------------------------------------
      subroutine print_variables()
      use siesta_geom, only: ucell
      use m_cell, only  : celli

      integer  :: i, jindex
      integer  :: ivars

      external          io_assign, io_close
      real(dp), external    :: volcel
      real(dp)              :: dummy_volume

         !     Update celli, just in case
         !     Wasteful for now, but we cannot assume
         !     that ucell is known in any other case,
         !     and celli is not universally updated yet.


       dummy_volume = volcel(ucell)
       unknown_cell = (abs(dummy_volume) < 1.0e-6_dp)
       if (.not. unknown_cell) then
          call reclat(ucell,celli,0)
       endif

      call io_assign( ivars )
      open(ivars, file="ZMATRIX_VARS", form='formatted',
     $     position='rewind', status='unknown')

        write(6,"(a)") "Z-matrix Symbol Section -------"
        write(6,"(a)") "Variables"
        do i = 1, nVars
           jindex = iVarstoZmat(i)
           if (.not.VaryZmat(jindex)) cycle
           print *, ZmatVarNames(i), user_print(jindex)
           write(unit=ivars, fmt="(a20,g25.15)")
     $          ZmatVarNames(i), user_print(jindex)
        enddo
        call io_close(ivars)

        write(6,"(a)") "Constants"
        do i = 1, nVars
           jindex = iVarstoZmat(i)
           if (VaryZmat(jindex)) cycle
           print *, ZmatVarNames(i), user_print(jindex)
        enddo
        write(6,"(a,/)") "------------ End of Z-matrix Information"

        end subroutine print_variables

      subroutine upper2lower(string,nchar)
C
C  upper2lower accepts a string of nchar characters and replaces
C  any lowercase letters by uppercase ones.
C
      character string*(*)
      integer nchar, ncopy, i, itemp
      intrinsic :: char
C
      ncopy = nchar
      if (ncopy.le.0) ncopy = len(string)
      do i = 1,ncopy
        if (lge(string(I:I),'A').and.lle(string(I:I),'Z')) then
          itemp = ichar(string(i:i))-ichar('A')+ichar('a')
          string(I:I) = char(itemp)
        endif
      enddo
      end subroutine upper2lower

      function cross_product(a,b) result(c)
      real(dp), dimension(3), intent(in) :: a, b
      real(dp), dimension(3)             :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      end function cross_product

      function rotate(x,v,ang) result(y)
      real(dp), dimension(3), intent(in) :: x, v
      real(dp), intent(in)               :: ang
      real(dp), dimension(3)             :: y

      ! Use Rodrigues' formula to rotate a vector x around
      ! an axis defined by the vector v, by an angle 'ang'
      ! The rotation is counterclockwise looking down from
      ! the head to the foot of v.
      ! v does not need to be a unit vector
      ! See, for example,  Wikipedia 

      real(dp)                :: cos_ang, sin_ang, dot_wx, v2
      real(dp), dimension(3)  :: cross_wx, w

      ! Construct a unit vector from v

      v2 = dot_product(v,v)
      if ( v2 == 0.0_dp) RETURN
      w  = v/sqrt(v2)

      cos_ang = cos(ang)
      sin_ang = sin(ang)
      cross_wx = cross_product(w,x)
      dot_wx = dot_product(w,x)
      y = x*cos_ang + cross_wx * sin_ang +
     $                dot_wx * w * ( 1.0_dp - cos_ang)

      end function rotate

      subroutine normalize(v)
      real(dp), dimension(3), intent(inout) :: v

      real(dp)                :: v2

      ! Construct a unit vector from v

      v2 = dot_product(v,v)
      if ( v2 == 0.0_dp) call die("Normalize: zero vector input")
      v  = v/sqrt(v2)

      end subroutine normalize

      function to_int(p) result (i)
      logical, intent(in)   :: p
      integer               :: i

      if (p) then
         i = 1
      else
         i = 0
      endif
      end function to_int

      subroutine write_canonical_ucell_and_Zmatrix(filename,iunit)
      use m_cell, only: write_canonical_ucell
!
!     Write a "canonical" Zmatrix with the following characteristics:
!     1. No symbolic variables or constants are used.
!     2. The position coordinates of the first atom in each molecule
!        are absolute cartesian coordinates.
!     3. Any "crystal" (non-molecule) coordinates are also absolute cartesians.
!     4. There is no provision for output of constraints.
!     5. The units used are those initially specified by the user, and are
!        specified in fdf form.
!     6. Unit cell information in fdf-compatible format is also included.
!
!     The output is to file "OUT.UCELL.ZMATRIX" by default, but a new filename
!     might be specified as an optional argument. If a unit number is used
!     output is directly appended to it.
!     
!
        implicit none
        
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: iunit
C
C Local variables
C
        integer             :: iu
        character(len=90)    :: fname

      if (present(iunit)) then
         iu = iunit
      else
         if (present(filename)) then
            fname = filename
         else
            fname = "OUT.UCELL.ZMATRIX"
         endif
         call io_assign( iu )
         open(iu, file=trim(fname), form='formatted',
     $        position='rewind', status='unknown')
      endif

      if (.not. unknown_cell) then
         call write_canonical_ucell(iunit=iu)
      endif
      call write_canonical_Zmatrix(iunit=iu)

      if (.not. present(iunit)) then
         call io_close(iu)
      endif

      end subroutine write_canonical_ucell_and_Zmatrix

      end module zmatrix
