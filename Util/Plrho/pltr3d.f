! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine pltr3d( v, gv, nc, color )

c **********************************************************************
c Plots a colored 3-dimensional triangle.
c Written by J.M.Soler. Nov. 1997.
c ********* Input ******************************************************
c real   v(3,3) : Vertex coordinates v(ixyz,ivertex)
c real   gv(3,3) : Gradient of the function whose isosurface determines
c                 the triangle as part of a solid surface.
c                 The brightness of illumination is determined by the
c                 angle between gv and the light source direction, and
c                 it is interpolated between the three vertices.
c integer nc         : Number of color functions which determine the
c                      actual color plotted (normally nc=0 or 1).
c real   color(3,nc) : Values of the color functions at the vertices
c **********************************************************************
C    6  10        20        30        40        50        60        7072

      implicit none
      integer  icback, icolor, nc
      real     color(3,nc), gv(3,3), v(3,3)
      external icback, icolor, pgqwin, pgqvsz

c Common arrays
      include 'plrho.h'

c Internal parameters
c maxc   : Maximun number of color functions
c eps    : Small number to avoid divisions by zero
      integer maxc
      real eps
      parameter ( maxc   = 3   )
      parameter ( eps    = 1.d-15 )

c Internal variables
      integer
     .  i, ic, icb, iv, ix, ix1, ix2, iy, iy1, iy2, jv
      real
     .  dx, dy, b, b1, b2, bv(3), c(maxc), c1(maxc), c2(maxc),
     .  g(3), g1(3), g2(3), r, ray, wxmax, wxmin, wymax, wymin,
     .  x, x1, x2, xmin, xmax,
     .  y, y1, y2, ymin, ymax, 
     .  z, z1, z2
      logical
     .  frstme
      save
     .  dx, dy, frstme, xmin, xmax, ymin, ymax
      data
     .  frstme /.true./

c Find window size in user units and pixels
      if (frstme) then
        call pgqwin( xmin, xmax, ymin, ymax )
        call pgqvsz( 3, wxmin, wxmax, wymin, wymax )
        ixmin = max( 0, nint(wxmin) )
        iymin = max( 0, nint(wymin) )
        ixmax = min( maxx-1, nint(wxmax) )
        iymax = min( maxy-1, nint(wymax) )
        dx = (xmax - xmin) / (ixmax - ixmin + 1)
        dy = (ymax - ymin) / (iymax - iymin + 1)
*       write(6,*) 'pltr3d: xmin, xmax =', xmin, xmax
*       write(6,*) 'pltr3d: ymin, ymax =', ymin, ymax
*       write(6,*) 'pltr3d: ixmin, ixmax =', ixmin, ixmax
*       write(6,*) 'pltr3d: iymin, iymax =', iymin, iymax
      endif

c Initialize z-buffer and pixmap
      if (frstme) then
        icb = icback( nc )
        do iy = iymin,iymax
          do ix = ixmin,ixmax
            zbuff(ix,iy) = -1.e30
            pixmap(ix,iy) = icb
          enddo
        enddo
      endif

c Find y-range of triangle
      y1 = min( v(2,1), v(2,2), v(2,3) )
      y2 = max( v(2,1), v(2,2), v(2,3) )
      iy1 = iymin + int( (y1-ymin)/dy - 0.5 ) + 1
      iy2 = iymin + int( (y2-ymin)/dy - 0.5 )
      iy1 = min( iy1, iymax )
      iy1 = max( iy1, iymin )
      iy2 = min( iy2, iymax )
      iy2 = max( iy2, iymin )

c Fill triangle, line by line, on pixmap
      do iy = iy1,iy2
        y = ymin + (iy+0.5) * dy
        x1 = 1.e30
        x2 = -1.e30

c       Loop on triangle edges
        do iv = 1,3
          jv = iv + 1
          if (jv .eq. 4) jv = 1

c         Find if the triangle edge intersects line iy
          r = (y - v(2,iv)) / (v(2,jv) - v(2,iv) + eps)
          if (r.gt.0. .and. r.lt.1.) then

c           Find intersection values
            x = v(1,iv) + (v(1,jv)-v(1,iv)) * r
            z = v(3,iv) + (v(3,jv)-v(3,iv)) * r
            do i = 1,3
              g(i) = gv(i,iv) + (gv(i,jv)-gv(i,iv)) * r
            enddo
*           b = bv(iv) + (bv(jv) - bv(iv)) * r
            do ic = 1,nc
              c(ic) = color(iv,ic) + (color(jv,ic) - color(iv,ic)) * r
            enddo

c           Find leftmost and rightmost intersections
            if (x .lt. x1) then
              x1 = x
              z1 = z
*             b1 = b
              do i = 1,3
                g1(i) = g(i)
              enddo
              do ic = 1,nc
                c1(ic) = c(ic)
              enddo
            endif
            if (x .gt. x2) then
              x2 = x
              z2 = z
*             b2 = b
              do i = 1,3
                g2(i) = g(i)
              enddo
              do ic = 1,nc
                c2(ic) = c(ic)
              enddo
            endif

          endif
        enddo

c       Draw line
        ix1 = ixmin + int( (x1-xmin)/dx - 0.5 ) + 1
        ix2 = ixmin + int( (x2-xmin)/dx - 0.5 )
        do ix = ix1,ix2
          x = xmin + (ix+0.5) * dx
          r = (x - x1) / (x2 - x1 + eps)
          z = z1 + (z2 - z1) * r
          if (z .gt. zbuff(ix,iy)) then
*           b = b1 + (b2 - b1) * r
*           b = max( b, 0.0 )
            do i = 1,3
              g(i) = g1(i) + (g2(i) -g1(i)) * r
            enddo
            if (g(3) .lt. 0.) then
              b = ray( g )
              do ic = 1,nc
                c(ic) = c1(ic) + (c2(ic) - c1(ic)) * r
              enddo
              pixmap(ix,iy) = icolor( b, nc, c )
*             write(6,*) 'pltr3d: ix,iy,pixmap =',ix,iy,pixmap(ix,iy)
              zbuff(ix,iy) = z
            endif
          endif
        enddo
      enddo

      frstme = .false.
      end




