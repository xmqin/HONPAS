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


        program readwfx

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
c the wavefunctions coefficients (NEW STYLE, SystemLabel.WFSX) and writes
c this information in user-friendly form to an output file.
c
c The program needs two input files:
c
c 1) Main input file, read by standard input. A sample of input file is:
c
c    --- begin input file ---
c        h2o.WFSX
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
c    information on the wavefunctions. In example above, h2o.WFSX
c****************************************************************************



        implicit none

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,
     .          nwflist, iw, indwf, j, nuotot, jj

        character fname*33, oname*33
        integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
        character(len=20), allocatable, dimension(:) :: symfio,labelfis

        real*4, allocatable, dimension(:,:) :: psi
        logical gamma

        real*8 k(3), energy,thress
 
        read(5,*) fname
        read(5,*) oname
        read(5,*) thress

        thress = thress**2  ! We use the square later on

        iu = 10
        io = 11

        open(iu, file=fname, form='unformatted', status='old' )
        open(io, file=oname, form='formatted', status='new' )

        rewind (iu)

        read(iu) nk, gamma

        read(iu) nspin
        read(iu) nuotot

        if (gamma) then
         allocate(psi(1,nuotot))
        else
         allocate(psi(2,nuotot))
        endif
        allocate(iaorb(nuotot), labelfis(nuotot),
     $           iphorb(nuotot), cnfigfio(nuotot),
     $           symfio(nuotot))

         read(iu) (iaorb(j),labelfis(j),
     .            iphorb(j), cnfigfio(j),
     .            symfio(j), j=1,nuotot)

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
          read(iu) nwflist

          write(io,*)
          write(io,'(a72)')    ' ***************************************
     .********************************'
          write(io,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik,
     .                                         k(1),k(2),k(3)
          write(io,'(a22,2x,i6)') 'Spin component = ',ispin
          write(io,'(a22,2x,i6)') 'Num. wavefunctions = ',nwflist



C Loop over wavefunctions 

          do iw = 1,nwflist

            read(iu) indwf
            read(iu) energy

            write(io,*)
            write(io,'(a22,2x,i6)') 'Wavefunction = ', indwf
            write(io,'(a22,2x,f10.6)') 'Energy (eV) = ', energy
            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'
            write(io,'(a72)')  '  Atom  Species Orb-global  Orb-in-atom
     . Orb-type      Re(psi)   Im(psi)'

            read(iu) (psi(1:,j), j=1,nuotot)
            do jj = 1,nuotot
              if (dot_product(psi(1:,jj),psi(1:,jj)) .gt. thress)
     $               then
              write(io,'(i6,5x,a20,1x,i10,8x,i3,7x,i1,a20,2(f10.6))') 
     .          iaorb(jj),labelfis(jj),jj,
     .          iphorb(jj), cnfigfio(jj),
     .          symfio(jj), psi(1:,jj)
              endif
            enddo

            write(io,'(a72)')  ' ---------------------------------------
     .--------------------------------'

          enddo
        enddo
      enddo


      close (iu)
      end
