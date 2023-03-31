! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      program grid1d

c****************************************************************************
c Calculates the value of a function (rho, drho ..) on a 1D grid in a given
c line, from the values of that function in a 3D grid (like out of SIESTA).
c It uses (tri)linear 3D interpolation (Numerical Recipes, 2nd ed. p. 117)
c
c Input.  
c  - fname : name of file with 3D information,
c  - nx    : number of points along the 1D
c  - ucell : .false. (default) plots using peridic bound. cond.
c            .true.  plots only within the unit cell, zero outside
c  - rawinp: controls type of input 
c   * .true. (raw-input) based on two 3D vectors: 
c        - ori(3): origin (point in 3D space where the segment starts)
c        - av(3) : vector that defines the 1D plot.
c        x-scale in Ang starting at origin.
c   * .false. (default): based on 2 points on the line (2 atomic positions)
c        - Two 3D vectors with the two positions, r0 and r1
c        - A 1D vector (scalar) a2 defining the length of the segment
c          in the basis of v1=r1-r0
c        - A 1D vector (scalar) o2 defining the origin of the plot, 
c          from r0 and again in the basis of v1=r1-r0
c        x-scale in Ang starting at r0
c
c Al input coordinates in Angstrom, internally too.
c
c Written by Emilio Artacho, December 1999
c****************************************************************************

      implicit none

      logical           outsid, ucell, rawinp

      integer           nx, ix, i, n, ip, np, is, 
     .                  mesh(3), m(3), mp1(3), nspin

      character         fname*25, fform*12, ffform*12
      parameter         (fform = 'unformatted')

      integer           maxp
      parameter         (maxp = 10000000)
      integer           ind,i1,i2

      real              rho(maxp,2)

      double precision  r0(3), r1(3), v1(3), av(3), ori(3), o2, a2, x,
     .                  vperp(3), r(3), rn(3), alp(3), sid(3), xfinal,
     .                  cell(3,3), celli(3,3), celmod(3), avmod, v1mod,
     .                  rhox, rho000, rho100, rho010, rho001, rho110, 
     .                  rho101, rho011, rho111, xincr, Ang, eps, orides

      data              ucell  /.false./
     .                  rawinp /.false./
     .                  Ang    /0.529177d0/
     .                  eps    /1.d-6/

c ---------------------------------------------------------------------------


c read data for defining the segment and 1D grid on it ------------------------

      read(5,*) fname
      write(6,"(a,16x,a/a/a,4x,a/a)") 
     .         '#             - G R I D  1 D -',
     .         'Emilio Artacho, December 1999', '#',
     .         '# 1D-grid values for 3D function in  : ', 
     .         fname, '#'
      read(5,*) nx
      read(5,*) ucell
      read(5,*) rawinp

      if ( rawinp ) then

c input via 1 vector and an origin defining a segment in space --------------

         read(5,*) (ori(i), i=1,3)
         read(5,*) (av(i), i=1,3)

         orides = 0.d0

         avmod = 0.d0
         do i=1, 3
            avmod = avmod + av(i)*av(i)
         enddo
         avmod = dsqrt(avmod)
         if (avmod .lt. eps) stop 'Vector has null size'

      else

c input via 2 points in the plane: 2 atom positions -------------------------

         read(5,*) (r0(i), i=1,3)
         read(5,*) (r1(i), i=1,3)
         read(5,*) o2
         read(5,*) a2

         write(6,"(a,3f12.6)") 
     .      '# Two points within the segment (Ang): ', (r0(i), i=1,3)
         write(6,"(a,3f12.6)") 
     .      '#                                      ', (r1(i), i=1,3)
         write(6,"(a)") '#'
         
         do i = 1, 3
           v1(i) = r1(i) - r0(i)
           ori(i) = r0(i) + o2*v1(i)
           av(i) = a2*v1(i)
         enddo

c module of av and v1 ------------------------------------------------------

         avmod = 0.d0
         v1mod = 0.d0
         do i = 1, 3
            avmod = avmod + av(i)*av(i)
            v1mod = v1mod + v1(i)*v1(i)
         enddo
         avmod = dsqrt(avmod)
         v1mod = dsqrt(v1mod)

         orides = o2*v1mod

      endif

      write(6,"(a,3f12.6)")
     .      '# Position of the origin in 1D plot  : ', (ori(i), i=1,3)
      write(6,"(a,3f12.6)")
     .      '# Vector defining the 1D window      : ', (av(i), i=1,3)
      write(6,"(a)") '#'

c ---------------------------------------------------------------------------


c read function from the 3D grid --------------------------------------------

      open( unit=1, file=fname, status='old', form=fform )

      ffform = fform
      if (ffform .eq. 'formatted') then
         read(1,*) cell
         read(1,*) mesh, nspin
         np = mesh(1) * mesh(2) * mesh(3)
         if (np .gt. maxp) stop 'grid1d: Parameter MAXP too small'
         ind = 0
         do is = 1,nspin
            do i1 = 1,mesh(3)
               do i2 = 1,mesh(2)
                  read(1,*) (rho(ind+ip,is),ip=1,mesh(1))
                  ind = ind + mesh(1)
               enddo
            enddo
         enddo
!!! OLD         read(1,*) ( (rho(ip,is), ip = 1, np), is = 1, nspin)
      else
         read(1) cell
         read(1) mesh, nspin
         np = mesh(1) * mesh(2) * mesh(3)
         if (np .gt. maxp) stop 'grid1d: Parameter MAXP too small'
        ind = 0
        do is = 1,nspin
          do i1 = 1,mesh(3)
            do i2 = 1,mesh(2)
              read(1) (rho(ind+ip,is),ip=1,mesh(1))
              ind = ind + mesh(1)
            enddo
          enddo
        enddo
!!! OLD         read(1) ( (rho(ip,is), ip = 1, np), is = 1,nspin)
      endif
      close(1)

      do i = 1, 3
         do n = 1, 3
            cell(i,n) = cell(i,n) * Ang
         enddo
      enddo

      write(6,"(a,3f12.6)")
     .   '# Cell vectors (Ang)                 : ', (cell(i,1), i=1,3)
      write(6,"(a,3f12.6)")
     .   '#                                      ', (cell(i,2), i=1,3)
      write(6,"(a,3f12.6)")
     .   '#                                      ', (cell(i,3), i=1,3)
      write(6,"(a)") '#'
      write(6,'(a,3i6)') 
     .   '# mesh                               : ', mesh
      write(6,'(a,3i6)') 
     .   '# nspin                              : ', nspin
      write(6,'(a)') '#'

c inverse of matrix of cell vectors -----------------------------------------

      call inver( cell, celli )

c lattice-vector moduli and size of sides of grid cell (small 'cube') -------

      do n = 1, 3
         celmod(n) = 0.d0
         do i = 1, 3
            celmod(n) = celmod(n) + cell(i,n)*cell(i,n)
         enddo
         sid(n) = celmod(n)/dfloat( mesh(n) )
      enddo


c scan the segment in fractional coordinates of av --------------------------

      xincr = 0.d0
      if (nx .gt. 1) xincr = 1.d0 / dfloat(nx -1)

c loop in 1D grid -----------------------------------------------------------

      write(6,"(a,8x,a,15x,a/a)") '#', 'x', 'F','#'

      x = -xincr
      do ix = 1, nx
         x = x + xincr

c find 3D position ----------------------------------------------------------

         do i = 1, 3
            r(i) = ori(i) + x*av(i) 
         enddo

c find 'cube' around it via components of r in cell vectors -----------------

         outsid = .false.

         do n = 1, 3
            rn(n) = 0.d0
            do i = 1, 3
               rn(n) = rn(n) + celli(n,i)*r(i)
            enddo

c m(i) stores de integers that define the lower corner in the grid: ---------
c the cube is between m1, m1+1, m2, m2+1, m3, and m+1, unless m.lt.1 --------

            if ( rn(n) .lt. 0.d0 ) then
               m(n) = int( rn(n)*mesh(n) )
            else
               m(n) = int( rn(n)*mesh(n) ) + 1
            endif
            mp1(n) = m(n) + 1

c fractional coordinates of the point within the small 'cube' ---------------

            alp(n) =  rn(n)*mesh(n) - dfloat(m(n)-1) 

c bring the points to the unit cell (using the traslational symmetry) -------

            if (m(n).lt.1) then
c              m(n) = m(n) + (1 - (m(n)-1)/mesh(n))*mesh(n)
               m(n) = mod(m(n)-1, mesh(n)) + mesh(n) + 1
               outsid = .true.
            endif
            if (m(n).gt.mesh(n)) then
               m(n) = mod(m(n)-1, mesh(n)) + 1
               outsid = .true.
            endif

            if (mp1(n).lt.1) then
c              mp1(n) = mp1(n) + (1 - (mp1(n)-1)/mesh(n))*mesh(n)
               mp1(n) = mod(mp1(n)-1, mesh(n)) + mesh(n) + 1
               outsid = .true.
            endif
            if (mp1(n).gt.mesh(n)) then
               mp1(n) = mod(mp1(n)-1, mesh(n)) + 1
               outsid = .true.
            endif

         enddo

c only do the interpolation if the 2D-plot point is within the 3D cell ------

         if (outsid .and. ucell) then
            rhox = 0.d0
         else
              
c values for the function at the 'cube' points ------------------------------
         
            rho000 = rho(  m(1)                +
     .                    (m(2)  -1) * mesh(1) +
     .                    (m(3)  -1) * mesh(1)*mesh(2) , 1)
            rho100 = rho(  mp1(1)              +
     .                    (m(2)  -1) * mesh(1) +
     .                    (m(3)  -1) * mesh(1)*mesh(2) , 1)
            rho010 = rho(  m(1)                +
     .                    (mp1(2)-1) * mesh(1) +
     .                    (m(3)  -1) * mesh(1)*mesh(2) , 1)
            rho001 = rho(  m(1)                +
     .                    (m(2)  -1) * mesh(1) +
     .                    (mp1(3)-1) * mesh(1)*mesh(2) , 1)
            rho110 = rho(  mp1(1)              +
     .                    (mp1(2)-1) * mesh(1) +
     .                    (m(3)  -1) * mesh(1)*mesh(2) , 1)
            rho101 = rho(  mp1(1)              +
     .                    (m(2)  -1) * mesh(1) +
     .                    (mp1(3)-1) * mesh(1)*mesh(2) , 1)
            rho011 = rho(  m(1)                +
     .                    (mp1(2)-1) * mesh(1) +
     .                    (mp1(3)-1) * mesh(1)*mesh(2) , 1)
            rho111 = rho(  mp1(1)              +
     .                    (mp1(2)-1) * mesh(1) +
     .                    (mp1(3)-1) * mesh(1)*mesh(2) , 1)
             
c for debugging -------------------------------------------------------------

c           write(6,*) m
c           write(6,*) mp1
c           write(6,*) alp 
c           write(6,"(8f9.5)")rho000,rho100,rho010,rho001,
c    .                        rho110,rho101,rho011,rho111

c interpolate ---------------------------------------------------------------

            rhox = 
     .           (1.d0-alp(1))*(1.d0-alp(2))*(1.d0-alp(3))*rho000 +
     .                 alp(1) *(1.d0-alp(2))*(1.d0-alp(3))*rho100 +
     .           (1.d0-alp(1))*      alp(2) *(1.d0-alp(3))*rho010 +
     .           (1.d0-alp(1))*(1.d0-alp(2))*      alp(3) *rho001 +
     .                 alp(1) *      alp(2) *(1.d0-alp(3))*rho110 +
     .                 alp(1) *(1.d0-alp(2))*      alp(3) *rho101 +
     .           (1.d0-alp(1))*      alp(2) *      alp(3) *rho011 +
     .                 alp(1) *      alp(2) *      alp(3) *rho111 

         endif

c write down ----------------------------------------------------------------

         xfinal = avmod*x + orides
         write(6,"(2d15.6)") xfinal, rhox

      enddo

      stop
      end


c****************************************************************************

c Inverts a 3x3 matrix

      subroutine inver(a,b)

      implicit none

      integer              i
      double precision     a(3,3), b(3,3), c

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(2,1) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(2,3) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

      do i = 1, 3
         c=1.d0/(a(1,i)*b(i,1) + a(2,i)*b(i,2) + a(3,i)*b(i,3) )
         b(i,1)=b(i,1)*c
         b(i,2)=b(i,2)*c
         b(i,3)=b(i,3)*c
      enddo

      end


c****************************************************************************

      subroutine inver2(a,b,n)

      implicit none

      integer               n, i, j, k
      double precision      a(n,n), b(n,n), x

      do i = 1, n
         do j = 1, n
           b(i,j) = a(i,j)
         enddo
      enddo

      do 4 i = 1, n
         x = b(i,i)
         b(i,i) = 1.d0
         do j = 1, n
            b(j,i) = b(j,i)/x
         enddo
         do 4 k = 1, n
            if(k-i) 2,4,2
 2          x = b(i,k)
            b(i,k) = 0.d0
            do j = 1, n
               b(j,k) = b(j,k) - b(j,i)*x
            enddo
 4    continue

      return
      end

