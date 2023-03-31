! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!-----------------------------------------------------------------------------
!
! module fsiesta
!
! Support routines for siesta-as-a-subroutine.
! This version requires SIESTA to be compiled together with its calling program. 
! It allows multiple siesta processes to run and communicate through MPI.
!
! Public procedures provided by this module:
!   siesta_launch   : Starts a siesta process
!   siesta_units    : Sets the physical units for calls to siesta_forces
!   siesta_forces   : Calculates atomic forces, energy, and stress
!   siesta_get      : Calculates other properties
!   siesta_quit     : Finishes a siesta process
!
! Public parameters, variables and arrays provided by this module:
!   none
!
! Interfaces of public procedures:
!
!   subroutine siesta_launch( label, nnodes, mpi_comm, launcher, localhost )
!     character(len=*),intent(in) :: label    : Name of siesta process
!                                               (prefix of its .fdf file)
!     integer,optional,intent(in) :: nnodes   : Number of MPI processes
!                                               reserved for each siesta process
!     integer,optional,intent(in) :: mpi_comm : MPI communicator defined by the
!                                               calling program for siesta use
!     character(len=*),optional,intent(in):: launcher (not used in this version)
!     logical,optional,intent(in) :: localhost : will siesta run at localhost?
!                                                (not used in this version)
!   end subroutine siesta_launch
!
!   subroutine siesta_units( length, energy )
!     character(len=*),intent(in) :: length : Physical unit of length
!     character(len=*),intent(in) :: energy : Physical unit of energy
!   end subroutine siesta_units
!
!   subroutine siesta_forces( label, na, xa, cell, energy, fa, stress )
!     character(len=*), intent(in) :: label      : Name of siesta process
!     integer,          intent(in) :: na         : Number of atoms
!     real(dp),         intent(in) :: xa(3,na)   : Cartesian coords
!     real(dp),optional,intent(in) :: cell(3,3)  : Unit cell vectors
!     real(dp),optional,intent(out):: energy     : Total energy
!     real(dp),optional,intent(out):: fa(3,na)   : Atomic forces
!     real(dp),optional,intent(out):: stress(3,3): Stress tensor
!   end subroutine siesta_forces
!
!   subroutine siesta_get( label, property, value, units )
!     character(len=*), intent(in) :: label      : Name of siesta process
!     character(len=*), intent(in) :: property   : Name of required magnitude
!     real(dp),         intent(out):: value      : Value of the magnitude
!                                                (various dimensions overloaded)
!     character(len=*), intent(out):: units      : Name of physical units
!   end subroutine siesta_get
!     
!   subroutine siesta_quit( label )
!     character(len=*),intent(in) :: label       : Name of siesta process
!   end subroutine siesta_quit
!
! Properties available through siesta_get:
!
!   property name:  'atomic_numbers'
!            size:  na (number of atoms per unit cell)
!            units: ' '
!
! Usage:
! - The typical expected profiles for calling siesta_launch are
!   - call siesta_launch(myLabel) => a siesta process is launched for each
!     label, with as many MPI processes as they use the same label. Each
!     siesta process will read a different myLabel.fdf data file.
!   - call siesta_launch(singleLabel,siestaNodes) => a number nSiestaProc =
!     (totalNodes/siestaNodes) of siesta processes are launched, each with 
!     siestaNodes MPI processes. All will read the same file singleLabel.fdf
!     The MPI processes of any given siesta process have conscutive ranks
!     so that siestaNode = int(MPInode/siestaNodes)
!   - call siesta_launch(singleLabel,mpi_comm=siestaComm) => the same as 
!     previous one, but with more control and flexibility in distributing 
!     MPI processes among siesta processes (notice: 'mpi_comm=' is mandatory)
! - Using nnodes/mpi_comm AND different labels is possible but care must be
!   paid to their consistency, i.e. different labels cannot occur within the
!   same siesta process.
! - A data file named label.fdf (e.g. H2O.fdf) must be present in the working
!   directory for each label used in siesta_launch. Parameter NumberOfAtoms
!   in label.fdf must coincide with input argument na in siesta_forces.
!   See siesta manual for format and required parameters in label.fdf
! - A statement 'MD.TypeOfRun forces' must be present in label.fdf, to istruct
!   siesta not to perform any dynamics or relaxation of its own, and to accept
!   instead the geometries from the master program
! - A pseudopotential data file for each atomic species must also be present
! - The call to siesta_launch can be omited in serial execution, and when
!   neither nnodes nor mpi_comm are present. If the call is present, it must
!   be before the first call to siesta_forces.
! - siesta_units may be called either before or after siesta_launch
! - The stress is defined as dE/d(strain)/Volume, with a positive sign
!   when the system tends to contract (negative pressure)
! - Output unit 6 should be open explicitly, assigning it to a file, before
!   calling siesta_launch. Otherwise, the output of print (or write(unit=6)),
!   after calling siesta_launch, may be written on a file named MAIN_OUTPUT.
!
! Sample usage for serial MD simulation:
!   use fsiesta, only: siesta_forces, siesta_get
!   character(len=32):: label='H2O', units
!   integer:: na
!   real*8 :: energy, dipole(3), cell(3,3), stress(3,3)
!   real*8 :: fa(3,maxAtoms), xa(3,maxAtoms)
!   do mdStep = 1,nSteps
!     ... find the new geometry. Angstroms & eV are assumed.
!     call siesta_forces( label, na, xa, cell, energy, fa, stress )
!   end do
!   call siesta_get( label, 'dipole', dipole, units )
!
! Sample usage for a nudged elastic band calculation:
!   use mpi,     only: MPI_Comm_World
!   use fsiesta, only: siesta_launch, siesta_units, siesta_forces
!   character(len=20):: label='Si_vacancy'
!   integer:: error, MPI_Comm_Siesta, myNode, myNudge
!   integer:: na, nNodes, nNudges, siestaNodes
!   real*8 :: energy, cell(3,3), stress(3,3), fa(3,maxAtoms), xa(3,maxAtoms)
!   close(6)
!   open(6,file='neb.out')
!   call MPI_Comm_Rank( MPI_Comm_World, myNode, error )
!   call MPI_Comm_Size( MPI_Comm_World, nNodes, error )
!   ! Set distribution of MPI processes by
!     siestaNodes = nNodes / nNudges
!     call siesta_launch( label, siestaNodes )
!   ! Or equivalently (but not simultaneously)
!     myNudge = myNode / nNudges
!     call MPI_Comm_Split( MPI_Comm_World, myNudge, 0, MPI_Comm_Siesta, error )
!     call siesta_launch( label, mpi_comm=MPI_Comm_Siesta )
!   call siesta_units( 'bohr', 'Ryd' )
!   do iStep = relaxSteps
!     ... find my nudge's new geometry
!     call siesta_forces( label, na, xa, cell, energy, fa, stress )
!   end do
!
! Behaviour:
! - Different siesta processes are distiguished because they either:
!   - have different labels in siesta_launch
!   - belong to different MPI groups, as set by nnodes or mpi_comm
! - If mpi_comm is present, nnodes is ignored
! - If neither nnodes nor mpi_comm are present, the MPI processes are
!   assigned to the siesta processes according to the input labels
! - If siesta_units is not called, length='Ang', energy='eV' are
!   used by default. If it is called more than once, the units in the
!   last call become in effect.
! - The physical units set by siesta_units are used for all the siesta
!   processes launched
! - If siesta_forces is called without a previous call to siesta_launch
!   for that label, an internal call siesta_launch(label) is generated
! - If argument cell is not present in the call to siesta_forces, or if
!   the cell has zero volume, it is assumed that the system is a molecule,
!   and a supercell is generated automatically by siesta so that the 
!   different images do not overlap. In this case the stress returned
!   has no physical meaning.
! - The following events result in a stopping error message:
!   - mpi_comm is not a valid MPI communicator
!   - two MPI processes in the same siesta process have a different label
!   - two different labels are used in any calls within the same MPI process
!   - siesta_launch is called twice
!   - siesta_quit is called before a call to siesta_launch or siesta_forces
!   - file label.fdf is not found
!   - NumberOfAtoms in label.fdf is different from na input in siesta_forces
!   - An input argument of siesta_forces differs in two MPI processes within
!     the same siesta process
!
! Written by J.M.Soler. Oct.2010
!-----------------------------------------------------------------------------

MODULE fsiesta

! Used module parameters and procedures
  use precision,        only: dp              ! Double precision real kind
  use sys,              only: die             ! Termination routine
  use m_siesta_init,    only: siesta_init     ! Siesta initialization
  use m_siesta_analysis,only: siesta_analysis ! Post processing calculations
  use m_siesta_move,    only: siesta_move     ! Move atoms
  use m_siesta_forces,  only: &
    external_siesta_forces => siesta_forces   ! Atomic force calculation
  use m_siesta_end,     only: siesta_end      ! Finish siesta
  use siesta_master,    only: siesta_server        ! Is siesta a server?
  use siesta_master,    only: siesta_subroutine    ! Is siesta a subroutine?
  use siesta_master,    only: setMasterUnits       ! Set physical units
  use siesta_master,    only: setCoordsFromMaster  ! Set atomic positions
  use siesta_master,    only: getForcesForMaster   ! Get atomic forces
  use siesta_master,    only: getPropertyForMaster ! Get other properties
  use siesta_master,    only: siesta_server        ! Is siesta a server?
  use siesta_master,    only: siesta_subroutine    ! Is siesta a subroutine?
  use siesta_master,    only: input_file           ! fdf data file
#ifdef MPI
  use mpi_siesta, only: MPI_Comm_World => true_MPI_Comm_World ! The true
                                                              ! MPI_COMM_WORLD
  use mpi_siesta, only: MPI_Comm_Siesta => MPI_Comm_World ! What siesta uses
                                                          ! as MPI_Comm_World
  use mpi_siesta, only: MPI_Integer           ! Integer data type
  use mpi_siesta, only: MPI_Character         ! Character data type
  use mpi_siesta, only: MPI_Double_Precision  ! Real double precision type
  use mpi_siesta, only: MPI_Max               ! Maximum-option switch
#endif

  implicit none

PUBLIC :: &
  siesta_launch, &! Start a siesta process
  siesta_units,  &! Set physical units
  siesta_forces, &! Calculate atomic forces, energy, and stress
  siesta_get,    &! Calculate other magnitudes
  siesta_quit     ! Finish siesta process

PRIVATE ! Nothing is declared public beyond this point

! Global module variables
  integer,      parameter :: maxLenLabel = 32
  integer,      parameter :: maxLenFile = 300
  logical,           save :: siesta_launched = .false.
  logical,           save :: siesta_quitted  = .false.
  logical,           save :: analysed = .false.
  logical,           save :: relaxed = .false.
  integer,           save :: step  = 0
  character(len=32), save :: xunit = 'Ang'
  character(len=32), save :: eunit = 'eV'
  character(len=32), save :: funit = 'eV/Ang'
  character(len=32), save :: sunit = 'eV/Ang**3'
  character(len=maxLenLabel):: myLabel = ' '
  character(len=*),parameter:: mainOutFileDef = 'MAIN_OUTPUT'
  character(len=maxLenFile) :: mainOutFile = mainOutFileDef

interface siesta_get
  module procedure siesta_get_rank0, siesta_get_rank1, siesta_get_rank2
end interface 

CONTAINS

!---------------------------------------------------

subroutine siesta_launch( label, nNodes, mpi_comm, launcher, localhost )
  implicit none
  character(len=*),  intent(in) :: label    ! Name of the siesta process
  integer, optional, intent(in) :: nNodes   ! Number of MPI processes to be used
  integer, optional, intent(in) :: mpi_comm ! MPI communicator to be used
  character(len=*),optional,intent(in):: launcher  ! Not used in this version
  logical,         optional,intent(in):: localhost ! Not used in this version

#ifdef MPI
  logical:: initialized, labelFound, mainOutFileOpened
  integer:: error, iColor, iNode, lenLabel, &
            myColor, myLenLabel, myNode, mySiestaNode, &
            nColors, siestaNodes, totNodes
  integer,  allocatable:: color(:)
  character(len=maxLenLabel):: rootLabel
  character(len=maxLenLabel),allocatable:: colorLabel(:), nodeLabel(:)
#endif

! DEBUG
!  character(len=32):: output_file
!  inquire( unit=6, name=output_file )
!  print*,'siesta_launch: output file = ', trim(output_file)
! END DEBUG

! Check that siesta has not been launched yet
  if (siesta_launched .or. siesta_quitted) &
    print*, 'siesta_launch: ERROR: siesta process already launched'

! Declare siesta as a server subroutine and set input file
  siesta_server     = .true.
  siesta_subroutine = .true.
  input_file = trim(label)//'.fdf'

! Store label
  myLabel = label

#ifdef MPI

! Initialise MPI unless siesta is running as a subroutine 
! of a driver program that may have initialized MPI already
  call MPI_Initialized( initialized, error )
  if (.not.initialized) call MPI_Init( error )

! Get total number of MPI processes (nodes) and my index among them
  call MPI_Comm_Size( MPI_Comm_World, totNodes, error )
  call MPI_Comm_Rank( MPI_Comm_World, myNode, error )

! Find maximum label length
  myLenLabel = len(trim(label))
  call MPI_AllReduce( myLenLabel, lenLabel, 1, MPI_Integer, &
                      MPI_Max, MPI_Comm_World, error )
  if (lenLabel>maxLenLabel) &
    call die('siesta_launch: ERROR: parameter maxLenLabel too small')

! Start parallel siesta process
  if (present(mpi_comm)) then                   ! Use input MPI communicator

    ! Check that mpi_comm is a valid MPI communicator
    call MPI_Comm_Size( mpi_comm, siestaNodes, error )
    if (error/=0) &
      call die('siesta_launch: ERROR: mpi_comm not a valid MPI communicator')

    ! Assign communicator to siesta
    MPI_Comm_Siesta = mpi_comm

  elseif (present(nNodes) .and. nNodes>0) then  ! Create new MPI communicator

    call MPI_Comm_Rank( MPI_Comm_World, myNode, error )
    myColor = myNode/nNodes
    call MPI_Comm_Split( MPI_Comm_World, myColor, 0, MPI_Comm_Siesta, error )

  else   ! Create new MPI communicator according to label(s)

    ! Collect all labels
    allocate( nodeLabel(totNodes) )
    call MPI_AllGather( myLabel,   maxLenLabel, MPI_Character, &
                        nodeLabel, maxLenLabel, MPI_Character, &
                        MPI_Comm_World, error )

    ! Assing processor colors according to labels
    allocate( color(totNodes), colorLabel(totNodes) )
    nColors = 0
    do iNode = 1,totNodes
      labelFound = .false.
      do iColor = 1,nColors
        if (nodeLabel(iNode)==colorLabel(iColor)) then
          labelFound = .true.
          color(iNode) = iColor
          exit ! iColor loop
        end if
      end do ! iColor
      if (.not.labelFound) then
        nColors = nColors + 1
        colorLabel(nColors) = nodeLabel(iNode)
        color(iNode) = nColors
      end if
    end do ! iNode

    ! Create new communicator
    myColor = color(myNode+1)
    call MPI_Comm_Split( MPI_Comm_World, myColor, 0, MPI_Comm_Siesta, error )
    deallocate( color, colorLabel, nodeLabel )

  end if ! (present(mpi_comm))

! Check that my siesta process has a unique label
  rootLabel = myLabel
  call MPI_Bcast( rootLabel, maxLenLabel, MPI_Character, &
                  0, MPI_Comm_Siesta, error )
  if (myLabel/=rootLabel) &
    call die('siesta_launch: ERROR: label mismatch in siesta process')

! Store name of output file of the calling program
  inquire( unit=6, name=mainOutFile )
  inquire( file=mainOutFile, opened=mainOutFileOpened )
  if (.not.mainOutFileOpened) mainOutFile = mainOutFileDef
! DEBUG
!  if (myNode==0) print*,'siesta_launch: mainOutFile= ',trim(mainOutFile)
! END DEBUG

! Create a directory for each siesta process, e.g.
!   mkdir -p H2O_proc07
!   cp -n H2O_proc07.* *.fdf *.psf *.vps *.ion H2O_proc07 2> /dev/null
  call MPI_Comm_Rank( MPI_Comm_Siesta, mySiestaNode, error )
  if (mySiestaNode==0) then
    call system('mkdir -p '//trim(label))
    call system('cp -n ' // trim(label) // '.* ' // trim(label))
    call system('cp -n *.fdf *.vps *.psf *.ion ' &
                // trim(label) // ' 2> /dev/null')
  endif

#endif

! Initialize flags and counters
  siesta_launched = .true.
  siesta_quitted  = .false.
  analysed        = .false.
  relaxed         = .false.
  step = 0

end subroutine siesta_launch

!---------------------------------------------------

subroutine siesta_units( length, energy )
  implicit none
  character(len=*), intent(in) :: length, energy
  xunit = length
  eunit = energy
  funit = trim(eunit)//'/'//trim(xunit)
  sunit = trim(eunit)//'/'//trim(xunit)//'**3'
  call setMasterUnits( xunit, eunit )
end subroutine siesta_units

!---------------------------------------------------

subroutine siesta_forces( label, na, xa, cell, energy, fa, stress )
  implicit none
  character(len=*),   intent(in) :: label
  integer,            intent(in) :: na
  real(dp),           intent(in) :: xa(3,na)
  real(dp), optional, intent(in) :: cell(3,3)
  real(dp), optional, intent(out):: energy
  real(dp), optional, intent(out):: fa(3,na)
  real(dp), optional, intent(out):: stress(3,3)

  real(dp):: e, f(3,na), myCell(3,3), s(3,3)
  integer :: myNode=0

#ifdef MPI
  integer :: error, n
  real(dp):: c(3,3), x(3,na)
#endif

! Launch siesta, if not yet done
  if (.not.siesta_launched) call siesta_launch( label )

! Check label
  if (label/=myLabel) call die('siesta_forces: ERROR: label mismatch')

! Set unit cell vectors
  if (present(cell)) then
    myCell = cell
  else
    myCell = 0
  end if

#ifdef MPI
! Check that input arguments are equal in all MPI processes
  call MPI_Comm_Rank( MPI_Comm_Siesta, myNode, error )
  n = na
  x = xa
  c = myCell
  call MPI_Bcast( n, 1, MPI_Integer, 0, MPI_Comm_Siesta, error )
  call MPI_Bcast( x, 3*na, MPI_Double_Precision, 0, MPI_Comm_Siesta, error )
  call MPI_Bcast( c, 3*3, MPI_Double_Precision, 0, MPI_Comm_Siesta, error )
  if (n/=na .or. any(x/=xa) .or. any(c/=myCell)) &
    call die('siesta_forces: ERROR: input mismatch among MPI processes')
#endif

#ifdef MPI
! Change directory and set output file
  call chdir(trim(label))
  close(unit=6)
  open(unit=6, file=trim(label)//'.out', position='append')
#endif

! Copy master's coordinates to master repository
  call setCoordsFromMaster( na, xa, myCell )

! BEGIN DEBUG: Print coords
  if (myNode==0) then
    write(6,'(/,2a)')          'siesta_forces: label = ', trim(label)
    write(6,'(3a,/,(3f12.6))') 'siesta_forces: cell (',trim(xunit),') =',myCell
    write(6,'(3a,/,(3f12.6))') 'siesta_forces: xa (',trim(xunit),') =', xa
    write(6,*) ' '
  end if
! END DEBUG

! Move atoms (positions will be read from master repository)
  if (step==0) then
    call siesta_init()
  else
    call siesta_move( step, relaxed )
  end if
  step = step + 1

! Calculate forces
  call external_siesta_forces( step )

! Get from master repository the forces for master
  call getForcesForMaster( na, e, f, s )

! BEGIN DEBUG: Print forces
  if (myNode==0) then
    write(6,'(/,3a,f12.6)')    'siesta_forces: energy (',trim(eunit),') =', e
    write(6,'(3a,/,(3f12.6))') 'siesta_forces: stress (',trim(sunit),') =', s
    write(6,'(3a,/,(3f12.6))') 'siesta_forces: forces (',trim(funit),') =', f
    write(6,*) ' '
  end if
! END DEBUG

! Copy results to output arguments
  if (present(energy)) energy = e
  if (present(fa))     fa     = f
  if (present(stress)) stress = s

! Flag that post processing analysis has not been done for this geometry
  analysed = .false.

#ifdef MPI
! Go back to parent directory and reset output file
  call chdir('..')
  close(unit=6)
  open(unit=6, file=trim(mainOutFile), position='append')
#endif

end subroutine siesta_forces

!---------------------------------------------------

subroutine siesta_get_rank0( label, property, value, units )
  character(len=*), intent(in) :: label      ! Name of siesta process
  character(len=*), intent(in) :: property   ! Name of required magnitude
  real(dp),         intent(out):: value      ! Value of the magnitude
  character(len=*), intent(out):: units      ! Name of physical units
  real(dp):: v(1)
  call siesta_get_value( label, property, 1, v, units )
  value = v(1)
end subroutine siesta_get_rank0

!---------------------------------------------------

subroutine siesta_get_rank1( label, property, value, units )
  character(len=*), intent(in) :: label      ! Name of siesta process
  character(len=*), intent(in) :: property   ! Name of required magnitude
  real(dp),         intent(out):: value(:)   ! Value(s) of the magnitude
  character(len=*), intent(out):: units      ! Name of physical units
  integer:: vsize
  vsize = size(value)
  call siesta_get_value( label, property, vsize, value, units )
end subroutine siesta_get_rank1

!---------------------------------------------------

subroutine siesta_get_rank2( label, property, value, units )
  character(len=*), intent(in) :: label      ! Name of siesta process
  character(len=*), intent(in) :: property   ! Name of required magnitude
  real(dp),         intent(out):: value(:,:) ! Value(s) of the magnitude
  character(len=*), intent(out):: units      ! Name of physical units
  integer:: vsize
  vsize = size(value)
  call siesta_get_value( label, property, vsize, value, units )
end subroutine siesta_get_rank2

!---------------------------------------------------

recursive subroutine siesta_get_value( label, property, vsize, value, units )
  character(len=*), intent(in) :: label         ! Name of siesta process
  character(len=*), intent(in) :: property      ! Name of required magnitude
  integer,          intent(in) :: vsize         ! Size of value array
  real(dp),         intent(out):: value(vsize)  ! Value(s) of the magnitude
  character(len=*), intent(out):: units         ! Name of physical units

  character(len=132):: error, message

#ifdef MPI
! Change directory and set output file
  call chdir(trim(label))
  close(unit=6)
  open(unit=6, file=trim(label)//'.out', position='append')
#endif

! Check that siesta has been launched
  if (.not.siesta_launched) call die('siesta_get: ERROR: siesta not launched')

! Check label
  if (label/=myLabel) call die('siesta_get: ERROR: label mismatch')

! Get property from master repository
  call getPropertyForMaster( property, vsize, value, units, error )

! Check that this property is returned with the right size
  if (error=='unknown_property') then
    if (analysed) then
      write(message,*) 'siesta_get: ERROR: not prepared for property = ', &
                        trim(property)
      call die( trim(message) )
    else ! (.not.analysed) => try again after performing post-processing
      call siesta_analysis( relaxed )
      analysed = .true.
      call siesta_get_value( label, property, vsize, value, units )
    end if ! (analysed)
  elseif (error=='wrong_size') then
    write(message,*) 'siesta_get: ERROR: wrong array size for property = ', &
                      trim(property)
    call die( trim(message) )
  end if ! (error=='unknown_property')

#ifdef MPI
! Go back to parent directory and reset output file
  call chdir('..')
  close(unit=6)
  open(unit=6, file=trim(mainOutFile), position='append')
#endif

end subroutine siesta_get_value

!---------------------------------------------------

subroutine siesta_quit( label )
  implicit none
  character(len=*), intent(in) :: label

#ifdef MPI
! Change directory and set output file
  call chdir(trim(label))
  close(unit=6)
  open(unit=6, file=trim(label)//'.out', position='append')
#endif

  if (.not.siesta_launched) then
    call die('siesta_quit: ERROR: no siesta process launched')
  elseif (label /= myLabel) then
    call die('siesta_quit: ERROR: label mismatch')
  endif

! Finish siesta process
  call siesta_end()

  siesta_launched = .false.
  siesta_quitted  = .true.

#ifdef MPI
! Go back to parent directory and reset output file
  call chdir('..')
  close(unit=6)
  open(unit=6, file=trim(mainOutFile), position='append')
#endif

end subroutine siesta_quit

END MODULE fsiesta

