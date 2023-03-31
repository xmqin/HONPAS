! This file has been contributed by Roberto Robles.
!
! It is included to provide some perspective on ways to
! plot the spin-texture information.
! It might need some re-formatting of the default spin-texture
! data provided by the program 'spin_texture' in this directory.
! An alternative to the procedure indicated is to select a single
! band with 'spin_texture', and do some semi-automatic processing
! of the resulting file for input to gnuplot or similar tools.
!
! One needs to plot a vector field where each point k is assigned a
! vector (Sx, Sy, Sz). To achieve this satisfactorily is not
! completely trivial, and different techniques can be used.
!
! Program to extract the spin texture from file spin_texture.dat as
!  printed by SIESTA. It requires the number of bands and their index
!  provided interactively.
!
!  If executed as "read_spin_texture_xsf" it gives the output as the
!  head of a xsf file, elsewhere it produces a .xyz file.
!
      program read_spin_texture

      double precision, allocatable :: kpoints(:,:), st(:,:,:), en(:,:)

      integer, allocatable :: bands(:)

      integer :: nkp, nen, nb, band
      double precision:: ef

      logical :: xyz=1

      character (len = 100) :: file, ab, name

      open(unit=9,file='spin_texture.dat',status='old')

      read(9,'(11x,i4,16x,i4,17x,f12.6)') nkp, nen, ef

      write(*,*) nkp, nen, ef

      allocate( kpoints(3,nkp) )
      allocate( st(3,nkp,nen) )
      allocate( en(nkp,nen) )
      allocate( bands(nen) )

      call get_command_argument(0,name)
      if(name.eq."read_spin_texture_xyz" &
         .OR.name.eq."./read_spin_texture_xyz") then
        xyz = .true.
      elseif(name.eq."read_spin_texture_xsf" &
         .OR.name.eq."./read_spin_texture_xsf") then 
        xyz = .false.
      endif

      if (xyz) then
        write(*,*) "Output in .xyz format"
      else
        write(*,*) "Output for .xsf format"
      endif

      do ik=1,nkp
        read(9,'(/14x,3f12.6/)') (kpoints(j,ik),j=1,3)
!        write(*,'(i4,3f12.6)') ik,(kpoints(j,ik),j=1,3)
          do ie=1,nen
            read(9,'(7x,f12.5,3f8.4)') en(ik,ie),(st(j,ik,ie),j=1,3)
!            write(*,'(i4,f12.5,3f8.4)') ie,en(ik,ie),(st(j,ik,ie),j=1,3)
          enddo
      enddo

      write(*,'(a)') "Number of bands: "
      read(*,'(i4)') nb
      do ib=1,nb
        write(*,'(a,i4)') "Index of band # ", ib
        read(*,'(i4)') bands(ib)
        if (bands(ib) .gt. nen) then
          write(*,*) "Band index can not be bigger than ", nen
          stop
        endif

        write(ab,'(i0)') bands(ib)
        if (xyz) file = "st_band_"//trim(ab)//".xyz"
        if (.not.xyz) file = "st_band_"//trim(ab)//".xsf"

        ifile = bands(ib)+100

        open(ifile,file=file)

        if (xyz) write(ifile,'(i4/)') nkp
        if (.not.xyz) write(ifile,'(a/a)') "set xsfStructure {","ATOMS"

        do ik=1,nkp
          if(xyz) write(ifile,'(a,3f12.6,3f8.4)') "X ",&
               (kpoints(j,ik)*10,j=1,2), (en(ik,bands(ib))-ef), &
               (st(j,ik,bands(ib)),j=1,3)
          if(.not.xyz) write(ifile,'(i3,3f12.6,3f8.4)') 0,&
               (kpoints(j,ik)*10,j=1,2), (en(ik,bands(ib))-ef), &
               (st(j,ik,bands(ib)),j=1,3)
        enddo
      enddo

      if (.not.xyz) write(ifile,'(a)') "}"

      end program read_spin_texture
