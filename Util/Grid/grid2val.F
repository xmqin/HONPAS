! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      program grid2val
      use f2kcli

c****************************************************************************
c GRID2VAL
c
c This program reads files with info on the grid from SIESTA
c and evaluates the corresponding function on given points
c
c Compilation:
c    
c        You need to use the f2kcli.F90 file located in the main Src
c        directory:
c
c        $(FC) $(FFLAGS) -o grid2val ../Src/f2kcli.F90 grid2val.F
c
c      If you use GFORTRAN, use instead:
c
c
c        $(FC) $(FFLAGS) -DGFORTRAN -o grid2val ../Src/f2kcli.F90 grid2val.F
c
c      FC and FFLAGS should be compatible with those you used to
c      compile Siesta, as grid2val reads a binary file whose data
c      layout can change. This includes 32bit vs 64 bit flags, if any.
c
c Usage:
c     
c        grid2val gridfile 
c
c   Points where the function is to be evaluated are read from
c   standard input, in the following (free) format
c
c   int1, int2, xf1, xf2, xf3,
c
c   where int1 and int2 are two arbitrary integers (for example,
c   the atom index and atomic species number, and xf1, xf2, xf3
c   are *fractional* coordinates.
c   This format is meant to be compatible with the .STRUCT files.
c
c   The values of the function are written to standard output.
c   Other information is written to standard error.
c
c   CAVEATS: 
c
c     -  For Files holding information for several spin components,
c        only the first is processed.
c
c     -  The function is computed using linear interpolation.
c
c****************************************************************************

      implicit none

      integer           ipt, isp, ix, iy, iz, i, j,
     .                  mesh(3), nspin, nargs, iostat

      character         fnamein*75

      integer, parameter :: dp = selected_real_kind(14,100)

      real, dimension(:,:,:), allocatable  :: rho
      real(dp)  :: cell(3,3)
      real  :: xfrac(3), val, fmin, fmax

c ---------------------------------------------------------------------------


      nargs = command_argument_count()
      if (nargs /= 1)  Stop "Usage: grid2val filename"

      call get_command_argument(1,value=fnamein,status=iostat)

c read function from the 3D grid --------------------------------------------

      open( unit=1, file=fnamein, form="unformatted", status='old' )

      read(1) cell
  
      write(0,*) 
      write(0,*) 'Cell vectors'
      write(0,*) 
      write(0,*) cell(1,1),cell(2,1),cell(3,1)
      write(0,*) cell(1,2),cell(2,2),cell(3,2)
      write(0,*) cell(1,3),cell(2,3),cell(3,3)

      read(1) mesh, nspin

      write(0,*) 
      write(0,*) 'Grid mesh: ',mesh(1),'x',mesh(2),'x',mesh(3)
      write(0,*) 
      write(0,*) 'nspin = ',nspin
      if (nspin > 1) write(0,*)
     $     "** Only 1st spin info can be read at this point"

      allocate(rho(0:mesh(1)-1, 0:mesh(2)-1, 0:mesh(3)-1))

          do iz=0,mesh(3)-1
             do iy=0,mesh(2)-1
                read(1) (rho(ix,iy,iz),ix=0,mesh(1)-1)
             enddo
          enddo

      close(1)

      fmin = minval(rho)
      fmax = maxval(rho)
!      print *, "size of rho:", shape(rho), size(rho,dim=1)
       write(0,*) "minval, maxval:", fmin, fmax

      write(0,*) "(Ready to read points: i, j, xf1,xf2,xf3)"
      do
!         print *, "Enter fractional coordinates"
         read(5,iostat=iostat,fmt=*) i, j, xfrac(:)
         if (iostat /= 0) exit
         call evaluate(rho(:,:,:),xfrac,val)
         write(6,"(i3,i3,3f10.6,3x,f12.5)")  i, j, xfrac, val
      enddo

      CONTAINS

      subroutine evaluate(d,xfrac,val)
c
      implicit none

      real, intent(in) ::  d(0:,0:,0:)
      real, intent(in) ::  xfrac(3)           ! Reduced coordinates of point
      real, intent(out) :: val

      integer n(3), lo(3), hi(3)
      real    r(3), x(3), y(3)
      integer i, j, k
      real ::  nk

      n(1) = size(d,dim=1)
      n(2) = size(d,dim=2)
      n(3) = size(d,dim=3)

c           Find the right 3D "grid cube" and the reduced coordinates
c           of the point in it. The double mod assures that negative
c           numbers are well treated (the idea is to bring the coordinates
c           to the [0,n(k)) interval)
c 
            
            do k = 1, 3
               nk = real(n(k))
               r(k) =  modulo(n(k)*xfrac(k),nk)
               lo(k) = int(r(k))
               hi(k) = mod ( lo(k)+1, n(k) )
               x(k) = r(k) - lo(k)
               y(k) = 1 - x(k)
!               print *, "rk, lok, hik, xk, yk:",
!     $              r(k), lo(k), hi(k), x(k), y(k)

            enddo

c      compute charge density by linear interpolation

            val     = d(lo(1),lo(2),lo(3)) * y(1) * y(2) * y(3) +
     &                d(lo(1),lo(2),hi(3)) * y(1) * y(2) * x(3) +
     &                d(lo(1),hi(2),lo(3)) * y(1) * x(2) * y(3) +
     &                d(lo(1),hi(2),hi(3)) * y(1) * x(2) * x(3) +
     &                d(hi(1),lo(2),lo(3)) * x(1) * y(2) * y(3) +
     &                d(hi(1),lo(2),hi(3)) * x(1) * y(2) * x(3) +
     &                d(hi(1),hi(2),lo(3)) * x(1) * x(2) * y(3) +
     &                d(hi(1),hi(2),hi(3)) * x(1) * x(2) * x(3) 


            end subroutine evaluate


      end program grid2val


