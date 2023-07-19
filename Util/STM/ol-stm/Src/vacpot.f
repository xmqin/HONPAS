
      subroutine vacpot(V0)

C Routine to calculate the value of the potential in the vacuum.
C Calculates the average of the potential on a plane at a given
C hight 'zvac' defined by the user, which should be in the middle
C of the vacuum region of the slab
C Returns the potential in eV!!!
C
C P. Ordejon, November 2004

      use fdf
      use precision

      implicit none

      double precision
     .  V0


C Internal variables

      integer 
     .  Ind, ic, ix, iy, iz, ip, mesh(3), np, nspin, npt,
     .  k1, k2, k3, unit1

      real, dimension(:), allocatable ::
     .  VH

      double precision
     .  cell(3,3), dxdm(3,3), z, zvac

      character
     .  sname*75, fname*80, paste*80

      logical
     .  found

      external 
     .  io_assign, io_close, paste


C Initialize variables

      V0 = 0.0d0
      npt = 0

C Read value of the vacuum Z position
      zvac = fdf_physical('STM.VacZ',1.0d40,'Bohr')
      if (zvac .gt. 0.99999d40) then
        write(6,*) 'ERROR: You must specify STM.VacZ in input'
        stop
      endif

C Assign file name and open file
      sname = fdf_string('SystemLabel','siesta')
      fname = paste( sname, '.VH' )
      call io_assign(unit1)
      inquire( file=fname, exist=found)
      if (.not. found) then
        write(6,*) 'ERROR: File ', fname, ' is not present'
        stop
      endif
      open(unit=unit1, file=fname, status='old', form='unformatted')

C Read info from VH file...
      read(unit1) cell
      read(unit1) mesh, nspin
      np = mesh(1) * mesh(2) * mesh(3)
      if (nspin .ne. 1) stop 'error: wrong spin in VH file'
      allocate (VH(np))
      call memory('A','S',np,'readvh')
      Ind=0
      do iz=1,mesh(3)
        do iy=1,mesh(2)
C read in potential
          read(unit1) (VH(Ind+ip),ip=1,mesh(1))
          Ind=Ind+mesh(1)
        enddo
      enddo
      call io_close(unit1)
C 

      IF (.not. monoclinic(cell)) then
        WRITE(6,*) 'error: the code only accepts monoclinic cells'
        WRITE(6,*) '       with Z as the vertical axis'
        STOP
      ENDIF

      do ic = 1,3
        do ix = 1,3
          dxdm(ix,ic) = cell(ix,ic) / mesh(ic)
        enddo
      enddo

C Loop on mesh points
      do k1 = 0,mesh(1)-1
      do k2 = 0,mesh(2)-1
C z-direction is scanned from top to bottom
      do k3 = mesh(3)-1,0,-1

C        Find mesh index of this point
        ip = 1 + k1 + mesh(1) * k2 + mesh(1) * mesh(2) * k3

C Calculate Z coordinate of this point:
        z = dxdm(3,3) * k3

        if (dabs(z-zvac) .lt. dxdm(3,3)) then
          V0 = V0 + VH(ip)
          npt = npt+1
        endif
      enddo
      enddo
      enddo

      V0 = V0/npt

C Potentials in Siesta are dumped to VH file in Ry. 
C Convert to eV

      V0 = V0 * 13.60581454d0

      write(6,*) 'Vacuum potential V0= ',V0,' eV'
      write(6,*)

      CONTAINS
      function monoclinic(cell)
      real(dp), intent(in) :: cell(3,3)
      logical monoclinic

      real(dp), parameter :: tol = 1.0e-8_dp

      monoclinic =  (abs(CELL(3,1)) < tol
     $         .and. abs(CELL(3,2)) < tol
     $         .and. abs(CELL(1,3)) < tol
     $         .and. abs(CELL(2,3)) < tol )

      print *, "monoclinic: ", monoclinic
      print *, cell

      end function monoclinic

      end





