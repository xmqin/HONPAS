! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine platom( sname, alpha, beta, gamma )

c ******************************************************
c Plots a frame of atom bonds. J.M.Soler. June'99
c ************ Input ***********************************
c character*(*) sname : System name
c ******************************************************

      implicit      none
      character*(*) sname
      real          alpha, beta, gamma

c Internal parameters and arrays
c   maxa   : Maximun number of atoms
c   maxna  : Maximun number of neigbors
c   rtol   : Maximun accepted value of rij/(ri+rj)
c   rijmin : Minimun bond length
c   rijmax : Maximun bond length
      integer maxa, maxna
      double precision rijmax, rijmin, rtol
      parameter ( maxa   =  1000 )
      parameter ( maxna  =   100 )
      parameter ( rtol   = 1.2d0 )
      parameter ( rijmin = 1.d-6 )
      parameter ( rijmax = 10.d0 )

      logical           found
      character         fname*30
      integer           i1, i2, i3, ia, isa(maxa), iv, ix, iza(maxa),
     .                  j, j1, j2, j3, ja, jn, jna(maxna), na, nna
      real              scell(3,3), xi(3), xj(3)
      double precision  cell(3,3), r2ij(maxna), range, ratom, rij,
     .                  xa(3,maxa), xij(3,maxna)
      external          ratom

c Look for coordinates file
      fname = trim(sname)//'.XV'
      inquire( file=fname, exist=found )
      if (.not.found) return

c Read coordinates
      open( unit=1, file=fname, status='old' )
      do iv = 1,3
        read(1,*) (cell(ix,iv),ix=1,3)
      enddo
      read(1,*) na
      call chkdim( 'platom', 'maxa', maxa, na, 1 )
      do ia = 1,na
        read(1,*) isa(ia),iza(ia),(xa(ix,ia),ix=1,3)
      enddo
      close( unit=1 )

C Initialize neighb subroutine 
      nna = maxna
      call neighb( cell, rijmax, na, xa, 0, 0,
     .             nna, jna, xij, r2ij )

c Loop on atoms
      do ia = 1,na

c       Look for neighbors
        nna = maxna
        call neighb( cell, rijmax, na, xa, ia, 0,
     .               nna, jna, xij, r2ij )
        call chkdim( 'platom', 'maxna', maxna, nna, 1 )

c       Draw bonds
        do jn = 1,nna
          ja = jna(jn)
          if (ja .le. ia) goto 10
          rij = sqrt( r2ij(jn) )
          range = (ratom(iza(ia))+ratom(iza(ja))) *
     .            rtol / 0.5292d0
          if ( rij .gt. rijmin .and.
     .         rij .lt. range ) then
            do ix = 1,3
              xi(ix) = xa(ix,ia)
              xj(ix) = xa(ix,ia) + xij(ix,jn)
            enddo
            call rotate( 1, xi, alpha, beta, gamma )
            call rotate( 1, xj, alpha, beta, gamma )
            call plin3d( xi, xj, 0, 0., 0. )
          endif
   10   enddo
      enddo

c Draw cell frame
      do i1 = 0,1
      do i2 = 0,1
      do i3 = 0,1
        do j = 1,3
          j1 = i1
          j2 = i2
          j3 = i3
          if (j.eq.1) j1 = mod(i1+1,2)
          if (j.eq.2) j2 = mod(i2+1,2)
          if (j.eq.3) j3 = mod(i3+1,2)
          do ix = 1,3
            xi(ix) = cell(ix,1)*i1 + cell(ix,2)*i2 + cell(ix,3)*i3
            xj(ix) = cell(ix,1)*j1 + cell(ix,2)*j2 + cell(ix,3)*j3
*           xi(ix) = xi(ix) * 0.999999 + 0.000001
*           xj(ix) = xj(ix) * 0.999999 + 0.000001
          enddo
          call rotate( 1, xi, alpha, beta, gamma )
          call rotate( 1, xj, alpha, beta, gamma )
*         call plin3d( xi, xj, 0, 0., 0. )
        enddo
      enddo
      enddo
      enddo

      end


      double precision function ratom( z )

c Returns the atomic radius. J.M.Soler June'99

      implicit none
      integer  z

      integer n
      parameter ( n = 5 )

      character*2      si(n)
      integer          i, zi(n)
      double precision ri(n)

      data (zi(i),si(i),ri(i),i=1,n) /
     .  1,'H', 0.373,   6,'C', 0.772,   7,'N', 0.549, 
     .  8,'O', 0.604,  15,'P', 1.105 / 

      do i = 1,n
        if (zi(i) .eq. z) then
          ratom = ri(i)
          return
        endif
      enddo
      write(6,*) 'ratom: no data for Z =', z
      ratom = 0.d0
      end


