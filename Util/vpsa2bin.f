! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
c$Id: vpsa2bin.f,v 1.5 2001/07/22 14:32:10 wdpgaara Exp $
c
      program vpsa2bin
c
      implicit none
c
c     This program converts pseudopotential "VPS" files (created by 
c     the ATOM program) from ASCII to Binary.
c
c     Use 12 significant digits (G-format) for the double precision
c     quantities.
c
      integer nrmax, lmax
      parameter (nrmax=2000, lmax=4)
c
      double precision a, b, zion
      character*2 nameat, corr
      character*3 rel
      character*4 core
      character*10 ray(6), title(7)
c
      integer i, j, lo, nrp, npotd, npotu
      double precision  r(nrmax), viod(lmax,nrmax), viou(lmax,nrmax),
     &                  cdc(nrmax), cdv(nrmax)
c
      character*70 message
c
      integer iargc, nargs
      character*70 ascii_file, binary_file
c
c     Let's get the files from the command line:
c
      nargs = iargc()
      if (nargs .ne. 2) then
         write(0,*) 'Usage: vpsa2bin ascii_file binary_file'
         stop 
      endif
c
      call getarg(1,ascii_file)
      call getarg(2,binary_file)
c
c      open files
c
      open(unit=1,file=ascii_file,form='formatted',status='old')
      rewind(1)
      open(unit=2,file=binary_file,form='unformatted',status='unknown')
      rewind(2)
c
      read(1,9000) nameat, corr, rel, core
      read(1,9010) (ray(j),j=1,6), (title(j),j=1,7)
      read(1,9015) npotd, npotu, nrp, a, b, zion
      if (nrp.gt.nrmax) stop 'NRMAX'

      write(2) nameat, corr, rel, core, (ray(j),j=1,6), 
     &         (title(j),j=1,7), npotd, npotu, nrp, a, b, zion
c
 9000 format(1x,a2,1x,a2,1x,a3,1x,a4)
 9010 format(1x,6a10,/,1x,7a10)
 9015 format(1x,2i3,i5,3g20.12)
c      
c     Note the format. Change if needed.
c
 8000 format(1x,i2)
 9030 format(4(g20.12))
 9040 format(1x,a)
c
c     Radial grid
c
      read(1,9040) message
      write(6,9040) message
      read(1,9030) (r(j),j=1,nrp)
      write(2) (r(j),j=1,nrp)
c 
c     "Down" potentials
c
      do 30 i = 1, npotd
         read(1,9040) message
         write(6,9040) message
         read(1,8000) lo
         if (lo .lt. 0 .or. lo .gt. (lmax-1)) stop 'LMAX'
         read(1,9030) (viod(lo+1,j),j=1,nrp)
         write(2) lo, (viod(lo+1,j),j=1,nrp)
   30 continue
c
c     "Up" potentials
c
      do 35 i = 1, npotu
         read(1,9040) message
         write(6,9040) message
         read(1,8000) lo
         if (lo .lt. 0 .or. lo .gt. (lmax-1)) stop 'LMAX'
         read(1,9030) (viou(lo+1,j),j=1,nrp)
         write(2) lo, (viou(lo+1,j),j=1,nrp)
   35 continue
c
c     Core and valence charge
c
      read(1,9040) message
      write(6,9040) message
      read(1,9030) (cdc(j),j=1,nrp)
      read(1,9040) message
      write(6,9040) message
      read(1,9030) (cdv(j),j=1,nrp)
c
      write(2) (cdc(j),j=1,nrp)
      write(2) (cdv(j),j=1,nrp)
c
      stop
c
      end

