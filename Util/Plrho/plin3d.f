! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine plin3d( v1, v2, nc, color1, color2 )

c **********************************************************************
c Plots a colored 3-dimensional line.
c Written by J.M.Soler. Nov. 1997.
c ********* Input ******************************************************
c real    v1(3)      : First vertex coordinates
c real    v2(3)      : Second vertex coordinates
c integer nc         : Number of color functions which determine the
c                      actual color plotted (normally nc=0 or 1).
c real    color1(nc) : Value of the color functions at 1st vertex
c real    color2(nc) : Value of the color functions at 2nd vertex
c **********************************************************************
C    6  10        20        30        40        50        60        7072

      implicit none
      integer  icback, icolor, nc
      real     color1(nc), color2(nc), v1(3), v2(3)
      external icback, icolor, pgqwin, pgqvsz

c Common arrays
      include 'plrho.h'

c Internal parameters
c maxc   : Maximun number of color functions
c eps    : Small number to avoid divisions by zero
      integer maxc
      real eps
      parameter ( maxc   = 3      )
      parameter ( eps    = 1.d-12 )

c Internal variables
      integer
     .  i, ic, icb, idy, iv, ix, ix1, ix2, iy, iy1, iy2, jv
      real
     .  dx, dy, b, b1, b2, bv(3), c(maxc), c1(maxc), c2(maxc),
     .  g(3), g1(3), g2(3), r, ray, rr, vv, wxmax, wxmin, wymax, wymin,
     .  x, x1, x2, xmin, xmax,
     .  y, y1, y2, ymin, ymax, 
     .  z, z1, z2
      logical
     .  frstme
      save
     .  dx, dy, frstme, xmin, xmax, ymin, ymax
      data
     .  frstme /.true./

*     write(6,*) 'plin3d: point 1'

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
*       write(6,*) 'plin3d: xmin, xmax =', xmin, xmax
*       write(6,*) 'plin3d: ymin, ymax =', ymin, ymax
*       write(6,*) 'plin3d: ixmin, ixmax =', ixmin, ixmax
*       write(6,*) 'plin3d: iymin, iymax =', iymin, iymax
      endif

*     write(6,*) 'plin3d: point 2'

c Initialize z-buffer and pixmap
*     if (frstme) then
*       icb = icback( nc )
*       do iy = iymin,iymax
*         do ix = ixmin,ixmax
*           zbuff(ix,iy) = -1.e30
*           pixmap(ix,iy) = icb
*         enddo
*       enddo
*     endif

*     write(6,*) 'plin3d: point 3'

c Check that both vertices are not equal
      vv = (v2(1)-v1(1))**2 + (v2(2)-v1(2))**2 + (v2(3)-v1(3))**2
      if (sqrt(vv) .lt. eps) goto 999

c Find reflection bright from the normal to v(2)-v(1)
      g(1) =   v2(2) - v1(2)
      g(2) = - v2(1) + v1(1)
      g(3) = - ( (v2(1)-v1(1))**2 + (v2(2)-v1(2))**2 ) /
     .         ( abs(v2(3)-v1(3)) + eps )
      b = ray( g )

*     write(6,*) 'plin3d: point 4'

c Find range of line
      if (v1(1) .lt. v2(1)) then
        x1 = v1(1)
        x2 = v2(1)
        y1 = v1(2)
        y2 = v2(2)
        z1 = v1(3)
        z2 = v2(3)
      else
        x1 = v2(1)
        x2 = v1(1)
        y1 = v2(2)
        y2 = v1(2)
        z1 = v2(3)
        z2 = v1(3)
      endif
      ix1 = ixmin + nint( (x1-xmin)/dx )
      ix2 = ixmin + nint( (x2-xmin)/dx )
      iy1 = iymin + nint( (y1-ymin)/dy )
      idy = nint( sign(1.,y2-y1) )
      rr = (x2-x1)**2 + (y2-y1)**2 + eps

      if (sqrt(rr).lt.dx .and. sqrt(rr).lt.dy) goto 999

*     write(6,*) 'plin3d: point 5'

c Draw line on pixmap
      ix = ix1
      iy = iy1
   10 continue

        x = xmin + (ix+0.5) * dx
        y = ymin + (iy+0.5) * dy

        r = ( (x-x1)*(x2-x1) + (y-y1)*(y2-y1) ) / rr
        z = z1 + (z2-z1) * r

        if ( ix.ge.ixmin .and. ix.le.ixmax .and.
     .       iy.ge.iymin .and. iy.le.iymax .and.
     .       z .gt. zbuff(ix,iy)) then
          do ic = 1,nc
            c(ic) = color1(ic) + (color2(ic) - color1(ic)) * r
          enddo
          pixmap(ix,iy) = icolor( b, nc, c )
*         write(6,*) 'plin3d: ix,iy,pixmap =',ix,iy,pixmap(ix,iy)
          zbuff(ix,iy) = z
        endif

        if (abs(y-y1)/(x-x1+eps) .gt. abs(y2-y1)/(x2-x1+eps)) then
          ix = ix + 1
        else
          iy = iy + idy
        endif
          
      if (r .lt. 1.) goto 10

*     write(6,*) 'plin3d: point 6'

  999 frstme = .false.
      end




