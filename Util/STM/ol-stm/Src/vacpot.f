! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

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
      use m_gridfunc
      use units, only: eV, Ang

      implicit none

      type(gridfunc_t) :: gf
      real(dp) :: V0

C Internal variables

      integer nspin, lun, n3, i, npt
      real, dimension(:,:), allocatable ::  average
      real(dp) :: z_coord, zvac, delta

      character  sname*75, fname*80
      external :: io_assign, io_close

! Read value of the vacuum Z position
      zvac = fdf_physical('STM.VacZ',huge(1.0_dp),'Bohr')
      if (zvac .gt. 1.0e10_dp) then
        write(6,*) 'ERROR: You must specify STM.VacZ in input'
        stop
      endif

! Assign file name and open file
      sname = fdf_string('SystemLabel','siesta')
      fname = trim(sname)//'.VH'
      call read_gridfunc(fname,gf)
      IF (.not. monoclinic_z(gf%cell)) then
        WRITE(6,*) 'error: the code only accepts monoclinic cells'
        WRITE(6,*) '       with Z as the vertical axis'
        STOP
      ENDIF
      if (gf%nspin .ne. 1) stop 'error: wrong spin in VH file'
      call get_planar_average(gf,3,average)

! Print V_ave(z) file for plotting
      call io_assign(lun)
      open(unit=lun,file=trim(sname)//'.v_ave_z',form="formatted",
     $     status="unknown", position="rewind")
      write(lun,"(a)") "#    z (Ang)     V_ave(z) (eV) "
      delta = gf%cell(3,3) / gf%n(3)
      n3 = size(average,dim=1)
      do i = 1, n3
         z_coord = gf%origin(3) + (i-1) * delta
         write(lun,"(f10.4,1x,f10.4)") z_coord/Ang, average(i,1)/eV
      enddo
      call io_close(lun)

      ! Find V0 by rough interpolation, as it was done in the
      ! old version.
      ! Since we should be in a flattish region, just average over
      ! close points. If there are no close points, we are out of range.

      V0 = 0.0_dp
      npt = 0
      do i = 1, n3
         z_coord = gf%origin(3) + (i-1) * delta
         if (abs(z_coord-zvac) < delta) then
            npt = npt + 1
            V0 = V0 + average(i,1)
         endif
      enddo
      if (npt == 0) then
         call die("Zvac out of range")
      endif

!     Potentials in Siesta are dumped to VH file in Ry. 
!     Convert to eV

      V0 = V0/eV/npt

      write(6,*) 'Vacuum potential V0= ',V0,' eV'

      end





