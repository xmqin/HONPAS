! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module m_struct_init
      private
      public :: struct_init
      
      CONTAINS

      subroutine struct_init()
      USE siesta_options
      use units, only: Ang
      use m_ioxv
      USE alloc,          only: re_alloc, de_alloc

      use siesta_geom
      use atmfuncs,       only : izofis
      use atomlist,       only : elem, iza
      use fdf
      use m_mpi_utils,    only : broadcast
      use periodic_table, only : symbol
      use m_iostruct,     only : read_struct
      use siesta_cml
      use files,          only : slabel
      use zmatrix,        only : lUseZmatrix
      use parallel,       only : IOnode
#ifdef NCDF_4
      use m_steps, only : inicoor, fincoor
      use m_exp_coord, only : read_exp_coord_options
      use m_exp_coord, only : exp_coord_init
      use m_exp_coord, only : exp_coord_next
#endif
      use siesta_master,  only: coordsFromMaster  ! Get coords from master prog
      use siesta_master,  only: setPropertyForMaster  ! Set a property to be
                                                      ! sent to master program
      use siesta_master,  only: siesta_server         ! Is siesta a server?

      implicit none

      real(dp), external :: volcel
      external :: automatic_cell, bonds, iozm, shaper
      
      character(len=80) :: dummy_str
      logical           :: foundxv, foundzm
      real(dp)          :: local_charnet

      real(dp):: bcell(3,3), tmp_cell(3,3)
      integer :: i, ia, ix
      integer:: nbcell  ! Auxiliary to call shaper routine
      real(dp), pointer :: xang(:,:)
      character(len=22) :: dyntyp
      character(len=60) :: restart_file
      logical :: found_restart, step_back

!------------------------------------------------------------------------- BEGIN
#ifdef DEBUG
      call write_debug( '  PRE struct_init' )
#endif

      !! Read number of atoms and coordinates, and unit cell
      use_struct_file = fdf_get("UseStructFile",.false.)
      use_struct_file = fdf_get("MD.UseStructFile",
     $     use_struct_file)    ! For legacy symbol
      if (use_struct_file) then
         call read_struct( na_u, tmp_cell) ! Sets na_u, xa, and isa
      else
         call coor(na_u,tmp_cell)  ! Sets na_u, xa, and isa
      endif
      ucell = tmp_cell  ! initialize module variable

!     Prepare iza here: it might be needed by ioxv      
      nullify( iza, va )
      call re_alloc( iza, 1, na_u, 'iza', 'struct_init' )
      call re_alloc( va, 1, 3, 1, na_u, 'va', 'struct_init' )

      ! Options read here instead of in siesta_options
      usesavexv = fdf_get('MD.UseSaveXV', usesaveddata)
      usesavezm = fdf_get('MD.UseSaveZM', usesaveddata)
      writic    = fdf_get('WriteCoorInitial', .true.)

      ! Read Z-matrix coordinates and forces from file
      if (lUseZmatrix) then
        foundzm = .false.
        if (usesavezm) then
            call iozm('read',ucell,vcell,xa,foundzm)
            if (IOnode) then
                if (.not.foundzm) then
                   write(6,'(/,a)') 'siesta: WARNING: ZM file not found'
                else
                    write(6,'(/,a)') 
     .        "! Info in XV file prevails over previous structure input"
                endif
            endif
        endif
      endif

! Read cell shape and atomic positions from a former run ..............
#ifdef NCDF_4
      dummy_str = fdf_get('MD.TypeOfRun','none')
      if ( leqi(dummy_str,'explicit') ) then
         ! Read in file name
         call read_exp_coord_options( )
         call exp_coord_init(slabel,na_u,inicoor,fincoor)
         if ( inicoor > 0 ) then
            call exp_coord_next(inicoor,na_u,xa)
         end if
      end if
#endif
      foundxv = .false.
      if (usesavexv) then
        call ioxv('read', ucell, vcell, na_u, isa, iza, xa, va, foundxv)
        if (IOnode) then
           if (.not.foundxv) then
              write(6,'(/,a)') 'siesta: WARNING: XV file not found'
           else
              write(6,"(a)")
     $       "! Info in XV file prevails over previous structure input"
           endif
        endif
!       For a Verlet/Nose/PR/NPR run with more than one time step, if
!       the RESTART file is not found, backward-propagate the atomic
!       positions by one time step using the Euler method

        istart = fdf_get('MD.InitialTimeStep',1)
        ifinal = fdf_get('MD.FinalTimeStep',1)
        if (foundxv .and. (ifinal-istart>0)) then
           dyntyp = fdf_get('MD.TypeOfRun','CG')
           step_back=.true.
           if (leqi(dyntyp,'verlet')) then
             restart_file = trim(slabel) // '.VERLET_RESTART'
           else if (leqi(dyntyp,'nose')) then
             restart_file = trim(slabel) // '.NOSE_RESTART'
           else if (leqi(dyntyp,'parrinellorahman')) then
             restart_file = trim(slabel) // '.PR_RESTART'
           else if (leqi(dyntyp,'noseparrinellorahman')) then
             restart_file = trim(slabel) // '.NPR_RESTART'
           else
             step_back=.false.
           endif
           if (step_back) then
              if (IOnode) then
                 inquire( file=restart_file, exist=found_restart )
              endif
              call broadcast(found_restart)
              if (.not. found_restart) then
                 ! dt_default = 1.0 (fs)
                 ! The user should not rely on the default here...
                 ! The parameter dt_default is private to read_options 
                 dt = fdf_get('MD.LengthTimeStep',1.0_dp,'fs')
                 xa=xa-dt*va
                 if (IOnode) then
                    write(6,'(5a)') 'WARNING: ', trim(restart_file),
     $                           ' not found--reading only XV file',
     $                           ' and moving back 1 time step using',
     $                           ' Euler'
                 endif
              endif
           endif
        endif
      endif
! ..................

! Find if siesta is a server of forces, unless we already know it .....
      if (.not.siesta_server) then
        dummy_str = fdf_get('MD.TypeOfRun',"none")
        siesta_server = (leqi(dummy_str,'forces') .or.
     .                    leqi(dummy_str,'master'))
      end if
! .....................

! Read cell shape and atomic positions from driver program
      if (siesta_server) then
        call setPropertyForMaster('atomic_numbers', 
     .                             na_u, real(iza(:),dp), ' ')
        call coordsFromMaster( na_u, xa, ucell )
      end if
! .....................

! Dump initial coordinates to output ..................................
      if ( writic.and.(IOnode) ) then
        write(6,'(/a)') 'siesta: Atomic coordinates (Bohr) and species'
        write(6,"('siesta: ',2x,3f10.5,i3,3x,i6)")
     .           ( (xa(ix,ia), ix=1,3), isa(ia), ia, ia=1, na_u)
      endif
! ..................

! Automatic cell generation ...........................................
      if (volcel(ucell) .lt. 1.0d-8) then
         local_charnet = fdf_get('NetCharge',0.0_dp)
         call automatic_cell(ucell,scell,na_u,xa,isa,local_charnet)
      endif

! Initialize atomic velocities to zero ................................
      if (.not. foundxv) then
        ! AG ** What happens with iozm call?
        va(1:3,1:na_u) = 0.0_dp
        vcell(1:3,1:3) = 0.0_dp
      endif
! ..................


! Find system shape ...................................................
      call shaper( ucell, na_u, isa, xa, shape, nbcell, bcell )
      if (IOnode) then
        write(6,'(/,2a)') 'siesta: System type = ', shape
      endif

! Find interatomic distances (output in file BONDS)
      if (IOnode) then
        rmax_bonds = fdf_physical("MaxBondDistance", 6.0_dp, "Bohr")
        call bonds( ucell, na_u, isa, xa,
     $       rmax_bonds, trim(slabel)// ".BONDS" )
      endif

! Output of initial system details:

      if (cml_p) then
! We need the names of the elements on node 0
        nullify(elem)
        allocate( elem(na_u) )
!        call re_alloc( elem, 1, na_u, 'elem', 'struct_init' )
        do i = 1, na_u
           elem(i) = symbol(izofis(isa(i)))
        enddo

        call cmlStartModule(xf=mainXML, title='Initial System')
        nullify(xang)
        call re_alloc( xang, 1, 3, 1, na_u, 'xang', 'struct_init' )
        xang = xa(1:3,1:na_u)/Ang
        call cmlAddMolecule(xf=mainXML, natoms=na_u,
     .       coords=xa/Ang, elements=elem, atomRefs=cisa)
        call de_alloc( xang, 'xang', 'struct_init' )

        call cmlAddLattice(xf=mainXML, cell=ucell, 
     .       units='siestaUnits:angstrom', dictref='siesta:ucell')
        call cmlAddProperty(xf=mainXML, value=trim(shape), 
     .       dictref='siesta:shape')
        call cmlEndModule(xf=mainXML)
      endif

#ifdef DEBUG
      call write_debug( '  POS struct_init' )
#endif
      END subroutine struct_init
      end module m_struct_init
