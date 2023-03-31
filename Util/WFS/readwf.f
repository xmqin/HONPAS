! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!


        program readwf

c****************************************************************************
c READWF  version 0.0.1
c
c This program READWF reads a the coefficients of wavefunctions
c on the expansion of the atomic orbitals basis set, as written
c by SIESTA (unformatted) and writes them in ascii, user-friendly
c form.
c
c Written by P. Ordejon, June 2003
c****************************************************************************
c USAGE:
c
c This program reads files generated from SIESTA, with information of
c the wavefunctions coefficients (SystemLabel.WFS) and writes
c this information in user-friendly form to an output file.
c
c The program needs two input files:
c
c 1) Main input file, read by standard input. A sample of input file is:
c
c    --- begin input file ---
c        h2o.WFS
c        h2o.ascii
c        0.00001
c    --- end input file ---
c
c    where:
c    - The first line is the name of the wavefunctions file genterated
c      by Siesta.
c    - The second line is the name of the output file
c    - The third line is a real number, specifying a thresshold to plot
c      the coefficients of the wavefunctions. If the norm of the weight
c      of a wavefunction on a given orbital is smaller than this thresshold
c      then it is not printed. This is useful for large systems, with
c      a very large number of basis orbitals.
c
c 2) The file with the information genetared by SIESTA, containing the
c    information on the wavefunctions. In example above, h2o.WFS
c****************************************************************************



        implicit none

        integer, parameter :: dp = selected_real_kind(10,100)

        integer nkmax
        parameter (nkmax=100000)

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,
     .          nwflist(nkmax), iw, indwf, iaorb, cnfigfio,
     .          j, iphorb, nuotot, jj

        character(len=20) labelfis, symfio
        character fname*33, oname*33

        real(dp) k(3), repsi,impsi,energy,thress
 
        read(5,*) fname
        read(5,*) oname
        read(5,*) thress

        iu = 10
        io = 11

        open(iu, file=fname, form='unformatted', status='old' )
        open(io, file=oname, form='formatted', status='new' )

        rewind (iu)

        read(iu) nk
        if (nk .gt. nkmax) then
          write(6,*) 'nkmax too small; increase it to ',nk
          stop
        endif
        read(iu) nspin
        read(iu) nuotot

        write(io,*)
        write(io,'(a22,2x,i6)') 'Nr of k-points = ',nk
        write(io,'(a22,2x,i6)') 'Nr of Spins = ',nspin
        write(io,'(a22,2x,i6)') 'Nr of basis orbs = ',nuotot
        write(io,*)


        do iik = 1,nk
          do iispin = 1,nspin

          read(iu) ik,k(1),k(2),k(3)
          if (ik .ne. iik) stop 'error in index of k-point'
          read(iu) ispin
          if (ispin .ne. iispin) stop 'error in index of spin'
          read(iu) nwflist(ik)

          write(io,*)
          write(io,'(a72)')    ' ***************************************
     .********************************'
          write(io,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik,
     .                                         k(1),k(2),k(3)
          write(io,'(a22,2x,i6)') 'Spin component = ',ispin
          write(io,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist(ik)


C Loop over wavefunctions 

          do iw = 1,nwflist(ik)

            read(iu) indwf
            read(iu) energy

            write(io,*)
            write(io,'(a22,2x,i6)') 'Wavefunction = ', indwf
            write(io,'(a22,2x,f10.6)') 'Energy (eV) = ', energy
            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'
            write(io,'(a72)')  '  Atom  Species Orb-global  Orb-in-atom
     . Orb-type      Re(psi)   Im(psi)'

            do jj = 1,nuotot
              read(iu) 
     .          iaorb,labelfis,j,
     .          iphorb, cnfigfio,
     .          symfio,
     .          repsi, impsi
              if (j .ne. jj) stop 'error in index of basis orbital'
              if (sqrt(repsi**2 + impsi**2) .gt. thress) then
              write(io,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,2(f10.6))') 
     .          iaorb,labelfis,j,
     .          iphorb, cnfigfio,
     .          symfio,
     .          repsi, impsi
              endif
            enddo


            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'

          enddo
        enddo
      enddo


      close (iu)
      close (iu)

      stop
      end
