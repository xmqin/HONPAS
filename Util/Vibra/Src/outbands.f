! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine outbands(iunits, nspin, maxuo, nuo, maxk, nk, 
     .                   nlines, lastk, label, kpoint, ek, ef)

      use fdf
C *********************************************************************
C Writes output with energies at k-point lines.
C Extracted from bands.f (written by J.Soler, August 1997),
C P.Ordejon, August 1998
C **************************** INPUT **********************************
C integer iunits              : = 0 leave unchanged
C                               = 1 convert from Ry to eV
C integer nspin               : Number of 'spin' components
C                               (=1 for phonons)
C integer maxuo               : maximum number of eigen. per k-point
C                               (first dimension of ek)
C integer nuo                 : number of eigenvalues per k-point
C integer maxk                : Last dimension of kpoint and ek
C integer nk                  : Number of band k-points
C integer nlines              : Number of band k-lines
C integer lastk(nlines)       : Index of last k-point in each k-line
C character*8 label(nlines)   : Label of each k-line
C real*8 kpoint(3,maxk)       : k point vectors
C real*8 ek(maxuo,maxk,nspin) : Eigenvalues
C real*8 ef                   : Fermi energy (=0.0  for phonons)
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C ***************** BEHAVIOUR *****************************************
C - Dispersion bands are saved in file sys_name.bands, where
C   sys_name is the value of fdf label SystemLabel, or 'dispersion'
C   by default.
C *********************************************************************
      implicit          none
      integer           iunits, maxk, maxuo, nk, nlines, nuo, nspin
      integer           lastk(nlines)
      double precision  ef, ek(maxuo,maxk,nspin), kpoint(3,maxk)
      character         label(nlines)*8
C *********************************************************************

C  Internal variables ...
      integer
     .  ik, il, io, ispin, iu

      double precision
     .  emax, emin, eV, path

      character(len=150):: 
     .  fname, sname, string

      external          io_assign, io_close
C ...

      if (iunits .eq. 0) then
        eV = 1.0d0
      else if (iunits .eq. 1) then
        eV = 1.d0 / 13.60580d0 
      else
        stop 'outdispl:  ERROR, wrong iunits'
      endif
C       Find name of output file and open it
        sname = fdf_string('SystemLabel','dispersion')
        fname = trim(sname)//'.bands'
        call io_assign(iu)
*       write(6,*) 'bands: iu,fname=', iu, fname
        open( iu, file=fname, status='unknown')

C       Write Reference Energy (for compatibility with band.f 
C       of Siesta)

        write(iu,*) ef/eV

C       Find and write the ranges of k and ek
        path = 0.d0
        emax = ek(1,1,1)
        emin = ek(1,1,1)
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          do ispin = 1,nspin
            do io = 1, nuo
              emax = max( emax, ek(io,ik,ispin) )
              emin = min( emin, ek(io,ik,ispin) )
            enddo
          enddo
        enddo
        write(iu,*) 0.d0, path
        write(iu,*) emin/eV, emax/eV

C       Write eigenvalues
C       (the 1 is to provide compatibility with bands.f of Siesta)

        write(iu,*) nuo, nspin, nk
        path = 0.d0
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          write(iu,'(f10.6,10f12.4,/,(10x,10f12.4))')
     .      path, ((ek(io,ik,ispin)/eV,io=1,nuo),ispin=1,nspin)
        enddo

C       Write abscisas of line ends and their labels
        write(iu,*) nlines
        il = 1
        path = 0.d0
        do ik = 1,nk
          if (ik .gt. 1)
     .      path = path + sqrt( (kpoint(1,ik)-kpoint(1,ik-1))**2 +
     .                          (kpoint(2,ik)-kpoint(2,ik-1))**2 +
     .                          (kpoint(3,ik)-kpoint(3,ik-1))**2 )
          if (ik .eq. lastk(il)) then
C           Put label between quotes
            if (label(il) .eq. ' ') then
              string = "' '"
            else
              string = "'"//trim(label(il))//"'"
            endif
            write(iu,'(f12.6,3x,a)') path, string
            il = il + 1
          endif
        enddo

C       Close output file
        call io_close(iu)


      end


