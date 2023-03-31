!     -------------------------------------------------------------
!     Copyright (c) Universite de Montreal
!     M-A.Malouin and F.El-Mellouhi, 2005
!     
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!     A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
!     OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!     SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!     THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!     -------------------------------------------------------------
!     Program use to write OpenDx file in .dx format from files
!     produce by siesta runs (use in VEP OpenDx program 'Dx View.net')
!
!     Write data for: 1) Electronic density 
!                     2) Atoms 
!                     3) Lattice Vectors
!                     4) Unit Cell Frame 
!                   
!     Written in fortran 90/95 ... use the same compiler as with siesta
!     for correct reading of unformatted binary file
!
!     The following compilation has been working fine...on linux and mac_os_x
!     e.g: xlf95 (or xlf90) -o DxFormat DxFormat.f 
!          ifort -o DxFormat DxFormat.f90 (you must change the file extension for this one)
!
!     -------------------------------------------------------------
!     This program convert siesta output from 
!     .RHO, .DRHO, .VH, .VT, .IOCH, .TOCH, .LDOS and .XV files
!     in .dx format and combine them in one big file for use 
!     with OpenDx visualization program DxView.net
!    
!     Notes on options:
!           1) Siesta uses periodic cell representation and some 
!              shifting corrections must in general be done for
!              isosurfaces to ensure centered visualization 
!           2) Atoms and isosurface are not writen in the same
!              coordinates system by Siesta and translation must
!              be done to ensure consistency with isosurface
!              (Unfortunately this work only for orthogonal cell
!               right now...but this will eventually be fixed) 
!           3) When comparing your results with your initial 
!              configuration, viewing only change data can be 
!              useful, so you can substract from your data any 
!              configuration that you want to use as a base 
!              reference to view only the changes
!              (Note: The difference is applied before conversion
!                     to .dx file format for OpenDx)
!
!     Special note: If the file systemLabel.xyz is found, the
!                   legend labels will be replaced by their 
!                   real labels instead of their default id code
!
!     Usage: Dxformat system_label [-s, -t, -d diff_label]
!       system_label: the one defined in the .fdf input file used 
!                     by Siesta for your simulation
!       Options: -s --> shifting of coordinates for isosurfaces (1)
!                -t --> translation of coordinates for atoms (2)
!                -d --> use files from second system_label 'diff_label'
!                       as Reference base data (3)
!
!     PLEASE READ THE README FILE FOR MORE INFORMATIONS ABOUT USING THIS
!     PROGRAM TOGETHER WITH SIESTA AND OPENDX
!
!     Version 1.0
!     @author Marc-Andre Malouin
!     @date July 12, 2005
!     -------------------------------------------------------------
      Program DxFormat
      Implicit none
    
      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

      Integer, Parameter :: MAX_LENGTH = 80, EXT_LENGTH = 5, OP_LENGTH = 2 ! For character array variables
      Integer, Parameter :: READ_UNIT = 1, WRITE_UNIT = 2, DIFF_UNIT = 3, XYZ_UNIT = 4 ! IO/Units
      
!     Index for loops
      Integer :: i

!     To check for i/o and allocation error
      Integer :: status = 0 ! If not zero after use: an error occured !
   
!     Various file names and extensions
      Character(MAX_LENGTH) :: system_label, diff_label, filename
!     The different type of isosurface process by this program
      Character(EXT_LENGTH), Parameter :: FILE_FORMAT(8)=(/".RHO ",".DRHO", ".VH  ", ".VT  ",".IOCH",".TOCH",".LDOS",".XV  "/)
      Character(EXT_LENGTH), Parameter :: ISO = ".ISO ", DX = ".dx  ", XYZ = ".xyz "
     
      Character(6), Parameter :: NO_ISO = '"None"'
      Integer :: nb_iso = 0, xv = 0
      Logical :: is_fileFound(7)= (/.false.,.false.,.false.,.false., .false.,.false.,.false./)

!     To write Unit Cell and Lattice Vectors info only once for iso files
      Logical :: add_uclv = .true.

!     Command line option input buffer
      Character(MAX_LENGTH) :: buffer
    
!     Command line arguments
      Integer :: nargs ! Number of arguments
      Integer, External :: iargc ! Function to retreive nargs
      
!     Valid command line options
      Character(OP_LENGTH), Parameter :: OPT_SHIFT_ISO = "-s"
      Character(OP_LENGTH), Parameter :: OPT_TRANSLATE_ATOMS = "-t"
      Character(OP_LENGTH), Parameter :: OPT_DIFF_ISO = "-d"
      
!     Default options behavior
      Logical :: shift_iso = .false.
      Logical :: translate_atoms = .false.
      Logical :: diff = .false. 
      
!     --- Io/Formats---
 10   Format ("Option '",A, "' use : ",A) ! For translation option
 20   Format ('Bad option argument: ',A) ! For bad argument

!     -------------------------------------------------------------
!     ---Reading command line arguments---
!     -------------------------------------------------------------

!     Number of arguments
      nargs = iargc()
      If(nargs == 0) Then
         Call printInfo()
         Call printUsage()
         Stop
      End if
      
!     First argument require: the system_label
      Call getarg(1, system_label)
      Write(*,*) "Processing for system_label: ", system_label    

!     No difference by default (difference only if diff_label != system_label)
      diff_label = system_label 

!     Reading optional arguments
      i=2
      Do while (i <= nargs)
         Call getarg(i, buffer)
         
         Select case (buffer)
         Case (OPT_SHIFT_ISO)
            Write(*,10) trim(buffer), 'shifting of coordinates for isosurfaces'
            shift_iso = .true.
         Case (OPT_TRANSLATE_ATOMS)
            Write(*,10) trim(buffer), 'translation of coordinates for atoms'
            translate_atoms = .true.
         Case (OPT_DIFF_ISO)
            i=i+1
            Call getarg(i, diff_label)

            If(trim(diff_label) == trim(system_label)) Then
               Write(*,10) trim(buffer), 'BUT diff_label = system_label...option skip (not used)'
            Else
               Write(*,10) trim(buffer), "using files from system_label '" // trim(diff_label) // "' as reference for difference"
            End if
         Case DEFAULT
            Write(*,20) trim(buffer)
            Call printUsage()
            Stop
         End select
         i=i+1
      End do
           
!     -------------------------------------------------------------
!     ---Processing files---
!     -------------------------------------------------------------

      filename = trim(system_label)//trim(DX)
      Open(UNIT=WRITE_UNIT, FILE=trim(filename), STATUS='Unknown', ACTION='write', IOSTAT=status, FORM='formatted')
      
!     IOError checking                                                       
      If(status /= 0) Then        
         Write(*,*) 'Error while creating or replacing the file: ' // trim(filename)
         Stop ! End of program                                
      End if   
      
      Do i=1, size(FILE_FORMAT)
        
!     Opening file
         filename = trim(system_label)//trim(FILE_FORMAT(i))
         Write(*,*)
         Write(*,*) 'Reading data from file: ', trim(filename)
         
         If(i == size(FILE_FORMAT)) Then
            Open(UNIT=READ_UNIT,FILE=filename,STATUS='old',ACTION='READ',IOSTAT=status,FORM='formatted')
         Else
            Open(UNIT=READ_UNIT,FILE=filename,STATUS='old',ACTION='READ', IOSTAT=status,FORM='unformatted')
         End if

         If(status == 0) Then
            is_fileFound(i) = .true.
         
            If(trim(diff_label) /= trim(system_label)) Then
               diff = .true.
!     Opening file for difference option
               filename = trim(diff_label)//trim(FILE_FORMAT(i))
              
               If(i == size(FILE_FORMAT)) Then
                  Open(UNIT=DIFF_UNIT,FILE=filename,STATUS='old',ACTION='READ',IOSTAT=status,FORM='formatted')
               Else
                  Open(UNIT=DIFF_UNIT,FILE=filename,STATUS='old',ACTION='READ',IOSTAT=status,FORM='unformatted')
               End if
              
               If(status /= 0) Then
                  Write(*,*) 'File "' // trim(filename) // '" for difference not found...option not used'
                  diff = .false. ! Option not used
               End if
            End if

!      Processing files
            If(i == size(FILE_FORMAT)) Then
               
!      Is the systemLabel.xyz present ?
               filename = trim(system_label)//trim(XYZ)
               Open(UNIT=XYZ_UNIT,FILE=filename,STATUS='old',ACTION='READ',IOSTAT=status,FORM='formatted')

               If(status == 0) Then
                  Write(*,*) '"' // trim(filename) // '" file found: changing legend labels to strings instead of numbers'
                  Call processAtoms(FILE_FORMAT(i), translate_atoms, diff,.true.) ! Yes
               Else
                  Call processAtoms(FILE_FORMAT(i), translate_atoms, diff,.false.) ! No
               End if

               xv = 1
               Close(XYZ_UNIT)

            Else
               Call processIso(FILE_FORMAT(i), shift_iso, diff, add_uclv)
               add_uclv = .false.
               nb_iso = nb_iso + 1
            End if
            
            If(diff) Close(DIFF_UNIT)
            
         Else ! An error occured or the file was not found
            Write(*,*) 'Error while opening the file: ' // trim(filename) // '  --> Skipping file...'
         End if
      
         Close(READ_UNIT)

      End do
      
!     Grouping members and data together in .dx file
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) '#'
      Write(WRITE_UNIT,*) '#    GROUPING MEMBERS '
      Write(WRITE_UNIT,*) '#'
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) 'object "iso" array type string rank 1 shape ',(EXT_LENGTH+1),' items ', nb_iso+1, ' data follows'
      Write(WRITE_UNIT,*) NO_ISO
      Do i=1, size(is_fileFound)
         If(is_fileFound(i)) Write(WRITE_UNIT,*) '"',trim(FILE_FORMAT(i)),'"'
      End do
      Write(WRITE_UNIT,*) 'attribute "XV" value number ', xv
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) 'object "all" class group'
      Write(WRITE_UNIT,*) 'member 0 value "iso"'
      Do i=1, size(is_fileFound)
         If(is_fileFound(i)) Write(WRITE_UNIT,*) 'member "' // trim(FILE_FORMAT(i)) // '" value "iso' // trim(FILE_FORMAT(i)) // '"'
      End do
      If(xv == 1) Write(WRITE_UNIT,*) 'member "' // trim(FILE_FORMAT(size(FILE_FORMAT))) // '" value "molecule"'
      If(nb_iso /= 0) Then
         Write(WRITE_UNIT,*) 'member "vectors" value "vectors' // trim(ISO) // '"'
         Write(WRITE_UNIT,*) 'member "cell" value "cell' // trim(ISO) // '"'
      End if
      
      close(WRITE_UNIT)
      
      End program

!     ------------------------------------------------------------------
!     Subroutines for processing files
!     ------------------------------------------------------------------
      Subroutine processIso (type, translate, use_diff, add_uclv)
      Implicit none

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

      Integer, Parameter :: EXT_LENGTH = 5 ! For type length
      Integer, Parameter :: READ_UNIT = 1, WRITE_UNIT = 2, DIFF_UNIT = 3 ! IO/Units
      Character(EXT_LENGTH), Parameter :: ISO = ".ISO "

!     Parameters
      Character(EXT_LENGTH), Intent(IN) :: type ! type = isosurface type (RHO, DRHO, VH, VT, IOCH or TOCH)
      Logical, Intent(IN) :: translate ! Translation correction ?
      Logical, Intent(IN) :: add_uclv ! Write Cell And Vector info ?
      Logical, Intent(IN) :: use_diff ! Use difference option or not

      Logical :: diff

!     For storing isosurfaces data read from input_file
!     (N.B: We don't use nsd...just use to read an extra token in file...number spin density information)

      Integer :: nsd, ntp ! ntp = ntm(1) * ntm(2) * ntm(3) (mesh volume)
      Integer :: ntm(3), tntm(3) ! The mesh dimension in x, y, z 
                                 ! (tntm to test system compatibility when use_diff is used)
      real(dp) :: ucell(3,3), tcell(3,3) ! Lattice unit cell 
                                         ! (tcell to test system compatibility when use_diff is used)


!     3D Electronic density projected in one dimension (at each x,y,z of mesh)
!     [dt is after translation (shift of origin for correct visualization)]
      real(sp), Allocatable :: ed(:), edt(:), ted(:) ! ted for difference option (when use_diff is specified)
     
      Integer :: delta(3) ! To shift the origin by an amount delta in x, y, z

!     For allocation errors testing
      Integer :: status

!     Dummy variable
      Integer :: temp 

!     Index for various loops  (ix,iy,iz,iix,iiy,iiz are use for triple loop for translation)
      Integer :: i, j, ix, iy, iz, iix, iiy, iiz
     

!     --- Io/Formats---
 200  Format (3F30.25) ! For printing matrix rows (lattice vectors) of ucell
 210  Format (A, ' ---> ',3F30.25)
 220  Format (3F15.10) ! A (x,y,z) point

!     --Reading--
      Read(READ_UNIT) ucell
      Read(READ_UNIT) ntm, nsd
      
      ntp = ntm(1) * ntm(2) * ntm(3)
      
!     Difference (Testing system compatibility)
      If(use_diff) Then
         diff = .true.

         Read(DIFF_UNIT) tcell
         Read(DIFF_UNIT) tntm, temp

         If(tntm(1) /= ntm(1) .OR. tntm(2) /= ntm(2) .OR. tntm(3) /= ntm(3) .OR. temp /= nsd) Then
            diff = .false.
         Else 
            Do i=1, 3
               If(tcell(i,1) /= ucell(i,1) .OR. tcell(i,2) /= ucell(i,2) .OR. tcell(i,3) /= ucell(i,3)) Then
                  diff = .false.
                  Exit
               End if
            End do
         End if
                 
         If(.NOT. diff) Write(*,*) 'Incompatible system specified for difference...option not used'
      Else
         diff = .false.
      End if
   
!     Printing values in terminal   
      If(add_uclv) Then 
         Write(*,*)
         Write(*,*) 'Common data for all isosurface files:'
         Write(*,*) 'Ntm size = ', ntp, ' = ', ntm(1), ' * ', ntm(2), ' * ', ntm(3), '     ( nsd = ', nsd,')'
         Write(*,*) 'Cell Vectors:'
         Write(*,210) "x'",(ucell(1,j), j=1, 3)
         Write(*,210) "y'",(ucell(2,j), j=1, 3)
         Write(*,210) "z'",(ucell(3,j), j=1, 3)
      End if

!     Array allocation       
      status = 0
      Allocate(ed(ntp),edt(ntp), STAT=status)
      If(diff) Allocate(ted(ntm(1)), STAT=status)
      If(status /= 0) Then ! Allocation error checking         
         Write(*,*) 'Array allocation error (overflow ???)...ending program !'
         Stop
      End if
      
      Write(*,*)
      Write(*,*) 'Reading charge density ... '
      Do iz=0, ntm(3)-1
         Write(*,*) 'Reading plane ', iz+1, ' ...'
         Do iy=0, ntm(2)-1
            Read(READ_UNIT) (ed(ix+iy*ntm(1)+iz*ntm(1)*ntm(2)), ix=1,ntm(1))
!     Difference
            If(diff) Then
               Read(DIFF_UNIT) (ted(i),i=1,ntm(1))
               Do i = 1, ntm(1)
                  ed(i+iy*ntm(1)+iz*ntm(1)*ntm(2)) = ed(i+iy*ntm(1)+iz*ntm(1)*ntm(2)) - ted(i)
               End do
            End if
         End do
      End do

!     Shifting of the origin for centered visualization of isosurface
      If(translate) Then
           
         Write(*,*)
         Write(*,*) 'Shifting of coordinates for visualization ...'
         Do i=1, 3
            delta(i) = ntm(i)/2.0 ! ucell(i,i)/2.0 *  ntm(i)/ucell(i,i) 
         End do
         Write(*,*) 'Delta shifting coefficients in x, y, z:', delta
            
         Do iz=1, ntm(3)
            Do iy=1, ntm(2)
               Do ix=1, ntm(1)
                  
!     Wise translation for not going out of bounds
                  If(ix > ntm(1)/2) Then
                     iix = ix - delta(1)
                  Else if(ix < ntm(1)/2) Then
                     iix = ix + delta(1)
                  End if
                  If(iy > ntm(2)/2) Then
                     iiy = iy - delta(2)
                  Else if(iy < ntm(2)/2) Then
                     iiy = iy + delta(2)
                  End if
                  If(iz > ntm(3)/2) Then
                     iiz = iz - delta(3)
                  Else if(iz < ntm(3)/2) Then
                     iiz = iz + delta(3)
                  End if
                  
!     Bounds testing (Just in case...should not happen !)
                  If(iix < 1) Stop 'iix < 0'
                  If(iiy < 1) Stop 'iiy < 0'
                  If(iiz < 1) Stop 'iiz < 0'
                  If(iix > ntm(1)) Stop 'iix > cell'
                  If(iiy > ntm(2)) Stop 'iiy > cell'
                  If(iiz > ntm(3)) Stop 'iiz > cell'
         
!     Translation of values around the new shifted origin
                  i = ix + (iy-1)*ntm(1) + (iz-1)*ntm(1)*ntm(2)
                  j = iix + (iiy-1)*ntm(1) + (iiz-1)*ntm(1)*ntm(2)
                  edt(j) = ed(i)
               End do
            End do
         End do
         
!     Final step (copy of translated coordinates)
         ed = edt

      End if

!     Writing data in .dx format 
      
!     Unit cell and lattice vectors info
      If(add_uclv)  Call write_UCLV(ucell, ISO)     
      
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) '#' 
      Write(WRITE_UNIT,*) '#    ' // trim(type) // 'FILE INFO'
      Write(WRITE_UNIT,*) '#'      
      Write(WRITE_UNIT,*) '############################################' 

      Write(WRITE_UNIT,*) 'object "positions' // trim(type) // '" class gridpositions counts', ntm(3), ntm(2), ntm(1)
                
      Write(WRITE_UNIT,*) 'origin 0 0 0'
      Do i=3, 1, -1
         Write(WRITE_UNIT,*) 'delta',((ucell(i,j)/ntm(i)),j=1,3)
      End do

      Write(WRITE_UNIT,*) 'object "connections' // trim(type) // '" class gridconnections counts', ntm(3), ntm(2), ntm(1)
      Write(WRITE_UNIT,*) 'object "data' // trim(type) // '" class array type float rank 0 items', ntp, 'data follows'
      
      Write(WRITE_UNIT,*) (ed(i), i=1, ntp)

      Write(WRITE_UNIT,*) 'object "iso' // trim(type) // '" class field'
      Write(WRITE_UNIT,*) 'component "positions" value "positions' // trim(type) // '"'
      Write(WRITE_UNIT,*) 'component "connections" value "connections' // trim(type) // '"'
      Write(WRITE_UNIT,*) 'component "data" value "data' // trim(type) // '"'
      Write(WRITE_UNIT,*)

      End subroutine
!     ------------------------------------------------------------------ 
      Subroutine processAtoms (type, translate, Use_diff, add_xyz)
      Implicit none

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

      Integer, Parameter :: MAX_LENGTH = 80, EXT_LENGTH = 5 ! For type length
      Integer, Parameter :: READ_UNIT = 1, WRITE_UNIT = 2, DIFF_UNIT = 3, XYZ_UNIT = 4 ! IO/Units

!     Parameters
      Character(EXT_LENGTH), Intent(IN) :: type ! type = atom type (.XV here)
      Logical, Intent(IN) :: translate
      Logical, Intent(IN) :: use_diff
      Logical, Intent(IN) :: add_xyz

      Logical :: diff, xyz

!     For storing atoms data from input file
      real(dp) :: ucell(3,3), tcell(3,3) ! Lattice unit cell (tcell to test system compatibility when use_diff is used)
      real(dp) :: N(3) !To translate the origin by an amount dx, dy and dz
      Character(MAX_LENGTH), allocatable :: atom_label(:) ! Atom labels (code number or string if add_xyz specified) 
      Integer, allocatable :: atomic_number(:) ! Atoms atomic number
      real(dp), allocatable :: x(:), y(:), z(:) ! Coordinates 
      real(dp) :: tx, ty, tz ! For difference opion (when use_diff is specified)
      Integer :: nb_atoms ! Total number of atoms (future size of all allocatable variables above)
   
      Integer :: ndiff_atoms ! Total number of different atoms
      Logical, allocatable :: diff_atom(:) ! If the atom moved or not (Use only if use_diff is specified)

!     For allocation error testing
      Integer :: status

!     Dummy variable
      Integer :: temp      

!     Index for various loops
      Integer :: i, j
      
!     ---Io/Formats----
 200  Format (3F30.25) ! For printing matrix rows (lattice vectors) of ucell
 210  Format (A, ' ---> ',3F30.25)
 220  Format (3F15.10) ! A (x,y,z) point

!     --Reading--
      Do i=1, 3
         Read(READ_UNIT,*) (ucell(i,j), j=1,3)
      End do
      Read(READ_UNIT,*) nb_atoms ! Number of atoms
   
!     Get .xyz file ready to read
      If(add_xyz) Then
         xyz = .true.
         Read(XYZ_UNIT,*) temp ! Skipping first line (the numbers of atoms)
         If(temp /= nb_atoms) Then
            xyz = .false.
            Write(*,*) 'Warning ! Number of atoms in file '//trim(type)//' and file .xyz do not match...legend labels not read...'
         End if
      Else
         xyz = .false.
      End if
      
!     Difference (Testing for system compatibility)
      If(use_diff) Then
         diff = .true.
         Do i=1, 3
            Read(DIFF_UNIT,*) (tcell(i,j), j=1,3)
         End do
         Read(DIFF_UNIT,*) temp ! Number of atoms
        
         If(temp /= nb_atoms) Then
            diff = .false.
         Else
            Do i=1, 3
               If(tcell(i,1) /= ucell(i,1) .OR. tcell(i,2) /= ucell(i,2) .OR. tcell(i,3) /= ucell(i,3)) Then
                  diff = .false.
                  Exit
               End if
            End do
         End If
                  
         If(.NOT. diff) Write(*,*) 'Incompatible system specified for difference...option not used'
      Else
         diff = .false.
      End if
   
!     Printing of values in terminal
      Write(*,*) 'Number of atoms = ', nb_atoms
      Write(*,*) 'Cell Vectors:'
      Write(*,210) "x'",(ucell(1,j), j=1, 3)
      Write(*,210) "y'",(ucell(2,j), j=1, 3)
      Write(*,210) "z'",(ucell(3,j), j=1, 3)

!     Array allocation
      status = 0
      Allocate(x(nb_atoms), y(nb_atoms), z(nb_atoms),STAT=status)
      If(status /= 0) Then ! Allocation error checking
         Write(*,*) 'Array allocation error (overflow ???)...ending program !'
         Stop
      End if
      Allocate (atom_label(nb_atoms), atomic_number(nb_atoms), diff_atom(nb_atoms), STAT=status)
      If(status /= 0) Then ! Allocation error checking
         Write(*,*) 'Array allocation error (overflow ???)...ending program !'
         Stop
      End if
 
      ndiff_atoms = 0
      Do i=1,nb_atoms       
         
         Read(READ_UNIT,*)atom_label(i),atomic_number(i),x(i),y(i),z(i) 
         
!        Replacing labels with text if .xyz file specified
         If(xyz) Read(XYZ_UNIT,*) atom_label(i), tx, ty, tz ! Only first data will be use..skipping already read coordinates x, y, z

         diff_atom(i) = .true. ! All atoms are different from nothing by default...unless otherwise specified

!     Difference
         If(diff) Then 
            Read(DIFF_UNIT,*)temp,temp,tx,ty,tz
!     If the atom moved, we keep it (it is different)...otherwise we remove it
            If(x(i) == tx .AND. y(i) == ty .AND. z(i) == tz) Then
               diff_atom(i) = .false.
               ndiff_atoms = ndiff_atoms + 1
            End if
         End if
      End do
      
!     Last consistency test for difference option
      If(nb_atoms == ndiff_atoms) Then
         Write(*,*) 'File use for difference is identical to original file...option not used'
         ndiff_atoms = 0
         Do i=1, nb_atoms
            diff_atom(i) = .true.
         End do
      End if
     
!     Translation of coordinates for consistency with isosurface
      If(translate) Then
         
         Do i=1, 3
!     Norm of lattice vectors
            N(i) = sqrt(ucell(i,1)**2 + ucell(i,2)**2 + ucell(i,3)**2)
         End do

         Write(*,*)
         Write(*,*) "Delta translation coefficients for atoms in each lattice vectors direction x', y', z': "
         Write(*,200) N/2

!     --Translating atom coordinates--
         Do i=1, nb_atoms
            x(i) = x(i) + (ucell(1,1) + ucell(2,1) + ucell(3,1))/2
            y(i) = y(i) + (ucell(1,2) + ucell(2,2) + ucell(3,2))/2
            z(i) = z(i) + (ucell(1,3) + ucell(2,3) + ucell(3,3))/2
         End do
      End if
      
      Write(*,*)
      Write(*,*) 'Writing atoms data...'

!     Writing to file
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) '#' 
      Write(WRITE_UNIT,*) '#    ATOMS COORDINATES AND COLOR INFO'
      Write(WRITE_UNIT,*) '#'      
      Write(WRITE_UNIT,*) '############################################'
        
!     Atoms positions
      Write(WRITE_UNIT,*) 'object "atomcoord" array type float rank 1 shape 3 items   ',  (nb_atoms-ndiff_atoms), ' data follows' 
      Do i=1, nb_atoms
         If(diff_atom(i)) Write(WRITE_UNIT,220) x(i), y(i), z(i)
      End do

!     Atoms size
      Write(*,*) '  Default atoms size: Base on unique atomic number...'

      Write(WRITE_UNIT,*) 'object "atomsize" array type float rank 0 items ',  (nb_atoms-ndiff_atoms), ' data follows'   
      Do i=1, nb_atoms
         If(diff_atom(i)) Write(WRITE_UNIT,*) (atomic_number(i) + 0.0)
      End do
!     Atoms code number id
      Write(*,*) '  Default atoms color: Base on unique code number...'

      Write(WRITE_UNIT,*)  &
        'object "label" array type string rank 1 shape ', &
         MAX_LENGTH,' items ', (nb_atoms-ndiff_atoms), ' data follows'   
      Do i=1, nb_atoms
         If(diff_atom(i)) Write(WRITE_UNIT,*) ('"'//trim(atom_label(i))//'"')
      End do
      Write(WRITE_UNIT,*) 'attribute "dep" value string "positions"'

      Write(WRITE_UNIT,*) 'object "atoms' // trim(type) // '" class field'
      Write(WRITE_UNIT,*) 'component "positions" value "atomcoord"'  
      Write(WRITE_UNIT,*) 'component "data" value "atomsize"'
      Write(WRITE_UNIT,*) 'component "id" value "label"'
      Write(WRITE_UNIT,*) 

!     Unit cell and lattice vectors info
      Call write_UCLV(ucell, type)
      
      Write(WRITE_UNIT,*) '############################################'
      Write(WRITE_UNIT,*) 'object "molecule" class group'
      Write(WRITE_UNIT,*) 'member "atoms" value "atoms' // trim(type) // '"'
      Write(WRITE_UNIT,*) 'member "vectors" value "vectors' // trim(type) // '"'
      Write(WRITE_UNIT,*) 'member "cell" value "cell' // trim(type) // '"'
      Write(WRITE_UNIT,*)
      
      End subroutine
!     ------------------------------------------------------------------
      Subroutine write_UCLV (ucell, label)
      Implicit none

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)
      
      Integer, Parameter :: WRITE_UNIT = 2 ! IO/Units
      
!     Parameters
      real(dp), Intent(IN) :: ucell(3,3) ! Lattice unit cell
      Character(5), Intent(IN) :: label ! For objects name

!     Index for various loops 
      Integer :: i, j
      
!     ---Io/Formats----
 200  Format (3F30.25) ! For printing matrix rows (lattice vectors) of ucell
 210  Format (A, ' ---> ',3F30.25)
 220  Format (3F15.10) ! A (x,y,z) point
      
!     LATTICE_VECTORS
      Write(WRITE_UNIT,*) '#    ' // trim(label) // ' LATTICE VECTOR INFO'
      Write(WRITE_UNIT,*) 'object "lattice_vectors" ' // trim(label) //  &
              '" class array type float rank 1 shape 3 items 3 data follows'   
      
!     Writing data (cell vectors)
      Do i=1, 3
         Write(WRITE_UNIT,200) (ucell(i,j), j=1, 3)
      End do
      Write(WRITE_UNIT,*) 'object "lattice_origin' // trim(label) // &
                           '" class array type float rank 1 shape 3 items 3 data follows'
      
!     Writing data (origin of cell)
      Do i=1, 3
         Write(WRITE_UNIT,220) 0.0,0.0,0.0
      End do
      
      Write(WRITE_UNIT,*) 'object "vectors' // trim(label) // '" class field'
      Write(WRITE_UNIT,*) 'component "data" value "lattice_vectors' // trim(label) // '"' 
      Write(WRITE_UNIT,*) 'component "positions" value "lattice_origin' // trim(label) // '"'
      
!     UNIT CELL FRAME
      Write(WRITE_UNIT,*) '#    ' // trim(label) // ' UNIT CELL FRAME INFO:'
      
      Write(WRITE_UNIT,*) 'object "cell_bonds' // trim(label) // '" class array type int rank 1 shape 2 items 12 data follows'
      
!     Writing data [One connection between #a and #b]
!     (N.B: See positions of each corner point below...)
      Write(WRITE_UNIT,*) ' 0  1' 
      Write(WRITE_UNIT,*) ' 0  2' 
      Write(WRITE_UNIT,*) ' 0  3' 
      Write(WRITE_UNIT,*) ' 1  4' 
      Write(WRITE_UNIT,*) ' 1  5' 
      Write(WRITE_UNIT,*) ' 3  5' 
      Write(WRITE_UNIT,*) ' 3  6' 
      Write(WRITE_UNIT,*) ' 2  6' 
      Write(WRITE_UNIT,*) ' 2  4' 
      Write(WRITE_UNIT,*) ' 7  5' 
      Write(WRITE_UNIT,*) ' 7  6' 
      Write(WRITE_UNIT,*) ' 7  4' 
      
      Write(WRITE_UNIT,*) 'attribute "element type" string "lines"'
      Write(WRITE_UNIT,*) 'object "cell_corners' // trim(label) // '" class array type float rank 1 shape 3 items 8 data follows'
      
!     Writing data
      Write(WRITE_UNIT,220) 0.0       , 0.0        , 0.0 ! #0 Origin
      Write(WRITE_UNIT,220) ucell(1,1), ucell(1,2) , ucell(1,3) ! #1 x
      Write(WRITE_UNIT,220) ucell(2,1), ucell(2,2) , ucell(2,3) ! #2 y
      Write(WRITE_UNIT,220) ucell(3,1), ucell(3,2) , ucell(3,3) ! #3 z
      Write(WRITE_UNIT,220) (ucell(1,1)+ucell(2,1)), &
                            (ucell(1,2)+ucell(2,2)), (ucell(1,3)+ucell(2,3))   ! #4 x+y
      Write(WRITE_UNIT,220) (ucell(1,1)+ucell(3,1)), &
                            (ucell(1,2)+ucell(3,2)), (ucell(1,3)+ucell(3,3))   ! #5 x+z
      Write(WRITE_UNIT,220) (ucell(2,1)+ucell(3,1)), &
                            (ucell(2,2)+ucell(3,2)), (ucell(2,3)+ucell(3,3))   ! #6 y+z
      Write(WRITE_UNIT,220)  &
           (ucell(1,1)+ucell(2,1)+ucell(3,1)), &
           (ucell(1,2)+ucell(2,2)+ucell(3,2)), &
           (ucell(1,3)+ucell(2,3)+ucell(3,3)) ! #7 x+y+z
      
      Write(WRITE_UNIT,*) 'object "bonds' // trim(label) // &
                           '" array type float rank 0 items 12 data follows'
      
!     Writing data
      Do i=1, 12
         Write(WRITE_UNIT,*) '1.0'
      End do
      Write(WRITE_UNIT,*) 'attribute "dep" string "connections"'
      
      Write(WRITE_UNIT,*) 'object "cell' // trim(label) // '" class field'
      Write(WRITE_UNIT,*) 'component "data" value "bonds' // trim(label) // '"'
      Write(WRITE_UNIT,*) 'component "positions" value "cell_corners' // trim(label) // '"'
      Write(WRITE_UNIT,*) 'component "connections" value "cell_bonds' // trim(label) // '"'
      Write(WRITE_UNIT,*)
      
      End subroutine

!     ------------------------------------------------------------------
!     Subroutines for printing usage and program description in terminal
!     ------------------------------------------------------------------
      Subroutine printUsage ()
      Implicit none
      
      Write(*,*) 'USAGE: Dxformat system_label [-s, -t, -d diff_label]'
      Write(*,*) '  system_label: the one defined in the .fdf input ' 
      Write(*,*) '                used by Siesta for your simulation'
      Write(*,*) '  Options: -s --> shifting realspace grid points for isosurfaces (1)'
      Write(*,*) '           -t --> translation of atoms coordinates (2)'
      Write(*,*) '           -d --> Use files from second system_label'
      Write(*,*) "                 'diff_label' as reference base data (3)"
      Write(*,*) '(Note: Use DxFormat alone to see program description)'
      
      End subroutine
!     ------------------------------------------------------------------
      Subroutine printInfo ()
      Implicit none

      Write(*,*) 'This program converts Siesta output from .RHO, .DRHO, '
      Write(*,*) '.VH, .VT, .IOCH, .TOCH, .LDOS and .XV files in .dx '
      Write(*,*) '.dx format and combine them in one big file for use '
      Write(*,*) 'with OpenDx visualization program DxView.net '
      Write(*,*)
      Write(*,*) 'Notes on options:'
      Write(*,*) ' 1) Siesta uses periodic cell representation and some '
      Write(*,*) '    shifting corrections must in general be done for'
      Write(*,*) '    isosurfaces to ensure centered visualization '
      Write(*,*) ' 2) Atoms and isosurface are not writen in the same '
      Write(*,*) '    coordinates system by Siesta and translation must'
      Write(*,*) '    be done to ensure consistency with isosurface'
      Write(*,*) '    (Unfortunately this works only for orthogonal cell'
      Write(*,*) '     for now...but this will eventually be fixed) '
      Write(*,*) ' 3) When comparing your results with your initial '
      Write(*,*) '    configuration, viewing only change data can be '
      Write(*,*) '    useful, so you can substract from your data any '
      Write(*,*) '    configuration that you want to use as a base '
      Write(*,*) '    reference to view only the changes'
      Write(*,*)'    (Note: The difference is applied before conversion'
      Write(*,*) '    to .dx file format for OpenDx)'
      Write(*,*)
      Write(*,*) 'Special note: If the file systemLabel.xyz is found, the'
      Write(*,*) '              legend labels will be replaced by their '
      Write(*,*) '              real labels instead of their default id code'
      Write(*,*)
      Write(*,*) 'For further info please read the README file'
      Write(*,*)
     
      End subroutine
!     ------------------------------------------------------------------
      
  

