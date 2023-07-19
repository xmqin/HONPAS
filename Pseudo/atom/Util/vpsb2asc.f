c$Id: vpsb2asc.f,v 1.3 2001/07/22 14:29:43 wdpgaara Exp $
c
      program vpsb2asc
c
      implicit none
c
c     This program converts pseudopotential "VPS" files (created by 
c     the ATOM program) from Binary to ASCII.
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
         write(0,*) 'Usage: vpsb2asc binary_file ascii_file'
         stop 
      endif
c
      call getarg(1,binary_file)
      call getarg(2,ascii_file)
c
c      open files
c
      open(unit=2,file=binary_file,form='unformatted',status='old')
      rewind(2)
      open(unit=6,file=ascii_file,form='formatted',status='unknown')
      rewind(6)
c
      read(2) nameat, corr, rel, core, (ray(j),j=1,6), 
     &         (title(j),j=1,7), npotd, npotu, nrp, a, b, zion
      if (nrp .gt. nrmax) stop 'NRMAX'

      write(6,9000) nameat, corr, rel, core
      write(6,9010) (ray(j),j=1,6), (title(j),j=1,7)
      write(6,9015) npotd, npotu, nrp, a, b, zion
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
      read(2) (r(j),j=1,nrp)
      write(6,9040) 'Radial grid follows' 
      write(6,9030) (r(j),j=1,nrp)
c 
c     "Down" potentials
c
      do 30 i = 1, npotd
         read(2) lo, (viod(lo+1,j),j=1,nrp)
         write(6,9040) 'Down Pseudopotential follows (l on next line)' 
         if (lo .gt. (lmax-1)) stop 'LMAX'
         write(6,8000) lo
         write(6,9030) (viod(lo+1,j),j=1,nrp)
   30 continue
c
c     "Up" potentials
c
      do 35 i = 1, npotu
         read(2) lo, (viou(lo+1,j),j=1,nrp)
         write(6,9040) 'Up Pseudopotential follows (l on next line)' 
         if (lo .gt. (lmax-1)) stop 'LMAX'
         write(6,8000) lo
         write(6,9030) (viou(lo+1,j),j=1,nrp)
   35 continue
c
c     Core and valence charge
c
      read(2) (cdc(j),j=1,nrp)
      read(2) (cdv(j),j=1,nrp)
c
      write(6,9040) 'Core charge follows' 
      write(6,9030) (cdc(j),j=1,nrp)
      write(6,9040) 'Valence charge follows' 
      write(6,9030) (cdv(j),j=1,nrp)
c
      stop
c
      end
