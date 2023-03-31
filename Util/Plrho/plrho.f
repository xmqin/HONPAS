! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      program plrho

c Plots the electron (spin) density, potential or LDOS as 
c a colored solid surface.
c Written by J.M.Soler. Nov/1997. Last modified Jan/2001.

      implicit none

c Common arrays
      include 'plrho.h'

c Internal variables
      character
     .  name*75, fform*12, fname*80, output*20, task*15
      logical
     .  found
      integer
     .  i, i1, i2, i3, icmin, icmax, ip, is, j, ind,
     .  mesh(3), ncolor, np, nsm, nspin, nt, pgopen
      real
     .  alpha, beta, cell(3,3), colavg, colmax, colmin,
     .  fvalue, f1max, f1min, f2max, f2min, f2zero,
     .  gamma, origin(3), rho, rx, ry,
     .  vpxmax, vpxmin, vsxmax, vsxmin,
     .  vpymax, vpymin, vsymax, vsymin,
     .  x, xmax, xmean, xmin, xrange,
     .  y, ymax, ymean, ymin, yrange
      real, allocatable ::
     .  f(:,:)
      double precision
     .  dcell(3,3)
      external
     .  grays, pgclos, pgenv, pgopen, pgpixl, plsurf,
     .  rotate
      data
     .  ncolor /0/
     .  origin /3*0.0/

c Read plot data
      open( unit=1, file='plrho.dat',
     .      status='old', form='formatted' )
      read(1,*) name
      read(1,*) task
      read(1,*) alpha, beta, gamma
      read(1,*) fvalue
      read(1,*) f2min, f2zero, f2max
      read(1,*) fform
      read(1,*) output
      close(1)

c Read density
      if (task .eq. 'ldos') then
        fname = trim(name)//'.LDOS'
      else
        fname = trim(name)//'.RHO'
      endif
      nsm = 1
      np = 0
      nspin = 0
      call iorho( 'read', fname, dcell, mesh, nsm, np, nspin, 
     .            f, found )
      if (found) then
        allocate( f(np,2) )
        call iorho( 'read', fname, dcell, mesh, nsm, np, nspin, 
     .              f, found )
        do i = 1,3
          do j= 1,3
            cell(j,i) = dcell(j,i)
          enddo
        enddo
      else
        write(6,*) 'plrho: ERROR: file not found: ', fname
        stop
      endif

      write(6,'(/,a,/,(3f12.6))') 'plrho: cell =', cell
      write(6,'(a,3i6)') 'plrho: mesh  =', mesh
      write(6,'(a,3i6)') 'plrho: nspin =', nspin

c Reorder spin density as sum and difference
      if (nspin .eq. 2) then
        ncolor = 1
        do ip = 1,np
          rho = f(ip,1) + f(ip,2)
          f(ip,2) = f(ip,1) - f(ip,2)
          f(ip,1) = rho
        enddo
      elseif (task .eq. 'spin') then
        stop 'plrho: data not spin polarized'
      endif

c Read potential
      if (task.eq.'vt' .or. task.eq.'vh') then
        ncolor = 1
        if (task .eq. 'vt') fname = trim(name)//'.VT'
        if (task .eq. 'vh') fname = trim(name)//'.VH'
        call iorho( 'read', fname, dcell, mesh, nsm, np, nspin, 
     .              f(1,2), found )
        if (.not.found) then
          write(6,*) 'plrho: ERROR: file not found: ', fname
          stop
        endif
      endif
        
c Rotate cell vectors
      call rotate( 3, cell, alpha, beta, gamma )

c Find the required window limits
      xmin =  1.e30
      xmax = -1.e30
      ymin =  1.e30
      ymax = -1.e30
      do i3 = 0,1
      do i2 = 0,1
      do i1 = 0,1
        x = cell(1,1) * i1 + cell(1,2) * i2 + cell(1,3) * i3
        y = cell(2,1) * i1 + cell(2,2) * i2 + cell(2,3) * i3     
        xmin = min( x, xmin )
        xmax = max( x, xmax )
        ymin = min( y, ymin )
        ymax = max( y, ymax )
      enddo
      enddo
      enddo
      xmin = xmin - 1.
      xmax = xmax + 1.
      ymin = ymin - 1.
      ymax = ymax + 1.

c Find range of shape and color functions
      f1min = 1.e30
      f1max = -1.e30
      do ip = 1,np
        f1min = min( f1min, f(ip,1) )
        f1max = max( f1max, f(ip,1) )
      enddo
      write(6,'(a,2e15.6)') 'plrho: f1min, f1max =', f1min, f1max
*     if (ncolor .eq. 1) then
*       f2min = 1.e30
*       f2max = -1.e30
*       do ip = 1,np
*         f2min = min( f2min, f(ip,2) )
*         f2max = max( f2max, f(ip,2) )
*       enddo
*       write(6,*) 'plrho: f2min, f2max =', f2min, f2max
*     endif

c Invert potential
      if (task.eq.'vh' .or. task.eq.'vt') then
        do ip = 1,np
          f(ip,2) = -f(ip,2)
        enddo
      endif

c Open pgplot window
      if (pgopen(output) .le. 0) stop 'pgtest: not able to open device'
      call pgqcol( icmin, icmax )
      call pgscir( icmin, icmax )
      call pgenv( xmin, xmax, ymin, ymax, 1, -2 )
      write(6,'(a,2i6)') 'plrho: color index range =', icmin, icmax

c Resize viewport to occupy full window
      call pgqvsz( 3, vsxmin, vsxmax, vsymin, vsymax )
      call pgqvp( 3, vpxmin, vpxmax, vpymin, vpymax )
      rx = (vsxmax-vsxmin) / (vpxmax-vpxmin)
      ry = (vsymax-vsymin) / (vpymax-vpymin)
      if (rx .lt. ry) then
        ry = ry / rx
        rx = 1.
      else
        rx = rx / ry
        ry = 1.
      endif
      xmean = (xmax + xmin) / 2.
      ymean = (ymax + ymin) / 2.
      xrange = xmax - xmin
      yrange = ymax - ymin
      xmin = xmean - xrange * rx / 2.
      xmax = xmean + xrange * rx / 2.
      ymin = ymean - yrange * ry / 2.
      ymax = ymean + yrange * ry / 2.
*     call pgsvp( 0.25, 0.75, 0.25, 0.75 )
      call pgsvp( 0., 1., 0., 1. )
      call pgswin( xmin, xmax, ymin, ymax )

c Create color table
      if (ncolor .eq. 0) then
        call grays()
      else
        call redblu( f2min, f2zero, f2max )
      endif

c Generate solid surface (ploting on pixmap array)
      call plsurf( origin, cell, mesh, mesh,
     .              f(1,1), fvalue, ncolor, f(1,2),
     .              nt, colmin, colavg, colmax )
      write(6,'(a,i8)') 'plrho: number of triangles =', nt
      if (ncolor .gt. 0) write(6,'(a,3f12.6)')
     . 'plrho: colmin, colavg, colmax =', colmin, colavg, colmax

c Draw bonding framework
*     call platom( name, alpha, beta, gamma )

c Copy pixmap to window and show it
      call pgpixl( pixmap, maxx, maxy,
     .             ixmin+1, ixmax+1, iymin+1, iymax+1,
     .             xmin, xmax, ymin, ymax )
      call pgclos
      end



      include 'icolor.f'
      include 'iorho.f'
      include 'plsurf.f'
      include 'pltr3d.f'
      include 'ray.f'
      include 'rotate.f'

*     include 'chkdim.f'
*     include 'dot.f'
*     include 'neighb.f'
*     include 'platom.f'
*     include 'plin3d.f'
*     include 'ranger.f'

