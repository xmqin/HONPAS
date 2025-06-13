
      program gnubands

***********************************************************************
* Utility for transorming SIESTA bands into gnuplot format
* Written by Emilio Artacho, February 1999.
***********************************************************************
* Read "systemlabel.bands" of SIESTA from standard input and use
* standard output fot gnuplot plotting
***********************************************************************

      implicit none

      integer           maxk, maxb, maxs
      parameter         (maxk=1000, maxb=100000, maxs=2)

      integer           nk, nspin, nband, ik, is, ib
      double precision  e(maxb, maxs, maxk), k(maxk)
      double precision  ef, kmin, kmax, emin, emax
      logical           overflow
c ---------------------------------------------------------------------

      read(5,*) ef
      read(5,*) kmin, kmax
      read(5,*) emin, emax
      read(5,*) nband, nspin, nk

      overflow = (nband.gt.maxb) .or. (nk.gt.maxk) .or. (nspin.gt.maxs)
      if (overflow) stop 'Dimensions in gnubands too small'

!      write(6,"(2a)") '# GNUBANDS: Utility for SIESTA to transform ',
!     .                'bands output into Gnuplot format'
!      write(6,"(a)") '#'
!      write(6,"(2a)") '#                                           ',
!     .                '       Emilio Artacho, Feb. 1999'
!      write(6,"(2a)") '# ------------------------------------------',
!     .                '--------------------------------'
!      write(6,"(a,f10.4)")  '# E_F               = ', ef
!      write(6,"(a,2f10.4)") '# k_min, k_max      = ', kmin, kmax
!      write(6,"(a,2f10.4)") '# E_min, E_max      = ', emin, emax
!      write(6,"(a,3i6)")    '# Nbands, Nspin, Nk = ', nband, nspin, nk
!      write(6,"(a)") '#'
!      write(6,"(a)") '#        k            E'
!      write(6,"(2a)") '# ------------------------------------------',
!     .                '--------------------------------'

      read(5,*) (k(ik),((e(ib,is,ik),ib=1,nband), is=1,nspin), ik=1,nk)

!      do ik = 1, nk
      do ib = 1, nband
         do ik = 1, nk
            write(6,"(3f14.6)") k(ik), e(ib,1,ik)-ef, e(ib,2,ik)-ef
        enddo
          write(6,*) 
!         write(6,*) "k2"
      enddo

      stop
      end
