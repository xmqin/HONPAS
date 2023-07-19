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


        program info_wfsx

c****************************************************************************
c USAGE:
c
c This program reads files generated from SIESTA, with information of
c the wavefunctions coefficients (NEW STYLE, SystemLabel.WFSX) and writes
c summary information in user-friendly form to standard output
c
c The input file must be called WFSX for now.
c****************************************************************************

        implicit none

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,
     .          nwflist, iw, indwf, j, nuotot, jj

        logical gamma

        real*8 k(3), energy

        iu = 10
        io = 6

        open(iu, file="WFSX", form='unformatted', status='old' )
        rewind (iu)

        read(iu) nk, gamma
        read(iu) nspin
        read(iu) nuotot

        read(iu)    

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

          write(io,'(a,2x,i6,2x,3f10.6,2x,a,i1,2x,a,i6)')
     $             "k-point: ",ik, k(1:3),
     $             " Spin: ", ispin, "# of wfs: ", nwflist

C Loop over wavefunctions 

          do iw = 1,nwflist

            read(iu) indwf
            read(iu) energy

            write(io,'(a,2x,i6,3x,a,f14.8)')
     $           '  ==> Wfn: ', indwf, "Energy (eV): ", energy
            read(iu) 

          enddo
        enddo
      enddo

      close (iu)
      end program info_wfsx
