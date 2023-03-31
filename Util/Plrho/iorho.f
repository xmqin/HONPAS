! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iorho( task, fname, cell, mesh, nsm, maxp, nspin,
     .                  f, found )

C *********************************************************************
C Saves/recovers the electron density at the mesh points.
C This simplified version reads only files written in serial mode.
C Writen by J.Soler July 1997.
C *************************** INPUT **********************************
C character*(*) task      : 'read'/'READ' or 'write'/'WRITE'
C character*(*) fname     : File name for input or output
C integer nsm             : Number of sub-mesh points per mesh point
C                           (not used in this version)
C integer maxp            : First dimension of array rho
C integer nspin           : Second dimension of array rho
C ************************** OUTPUT **********************************
C integer maxp            : Required first dimension of array rho,
C                           equal to mesh(1)*mesh(2)*mesh(3)
C                           Set only when task='read' and required
C                           value is larger than input value
C integer nspin           : Number of spin polarizations (1 or 2)
C logical found           : Were data found? (only when task='read')
C ******************** INPUT or OUTPUT (depending on task) ***********
C real*8  cell(3,3)       : Lattice vectors
C integer mesh(3)         : Number of mesh divisions of each
C                           lattice vector
C real    f(maxp,nspin)   : Electron density
C                           Notice single precision in this version
C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C ******************** BEHAVIOUR **************************************
C If task='read', and the values of maxp or nspin on input are less than
C those required to copy the array f from the file, then the required
C values of maxp and nspin are returned on output, but f is not read.
C *********************************************************************

      implicit          none

C Arguments
      character*(*)     fname, task
      integer           maxp, mesh(3), nspin, nsm
      real              f(maxp,nspin)
      double precision  cell(3,3)
      logical           found

c Internal variables and arrays
      character fform*11
      integer   i2, i3, ind, ip, is, np, ns

c Fix whether formatted or unformatted files will be used
      fform = 'unformatted'

c Look for data file
      inquire( file=fname, exist=found )
      if (.not.found) return

c Read unit cell vectors, number of mesh points and spin components
      open( unit=1, file=fname, status='old', form=fform )
      if (fform .eq. 'formatted') then
        read(1,*) cell
        read(1,*) mesh, ns
      else
        read(1) cell
        read(1) mesh, ns
      endif

c Read density (only if array f is large enough)
      np = mesh(1) * mesh(2) * mesh(3)
      if (ns.gt.nspin .or. np.gt.maxp) then
        maxp = np
      else
        if (fform .eq. 'formatted') then
          do is = 1,ns
            ind = 0
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
                read(1,*) (f(ind+ip,is),ip=1,mesh(1))
                ind = ind + mesh(1)
              enddo
            enddo
          enddo
        else
          do is = 1,ns
            ind = 0
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
                read(1) (f(ind+ip,is),ip=1,mesh(1))
                ind = ind + mesh(1)
              enddo
            enddo
          enddo
        endif
      endif
      close(1)
      nspin = ns
      end
