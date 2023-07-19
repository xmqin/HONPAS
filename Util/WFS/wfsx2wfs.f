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


        program wfsx2wfs

!
!       Converts a WFSX file (well-packed, single precision), to
!       a WFS file (old style wf file).
!
        implicit none

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,
     .          nwflist, iw, indwf, j, nuotot, jj

        character fname*33, oname*33
        integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
        character(len=20), allocatable, dimension(:) :: symfio,labelfis

        real*4, allocatable, dimension(:,:) :: psi
        logical gamma

        real*8 k(3), energy
 
        iu = 10
        io = 11

        open(iu, file="WFSX", form='unformatted', status='old' )
        open(io, file="WFS", form='unformatted', status='new' )

        rewind (iu)
        read(iu) nk, gamma
        read(iu) nspin
        read(iu) nuotot

        rewind (io)
        write(io) nk
        write(io) nspin
        write(io) nuotot

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

        do iik = 1,nk
          do iispin = 1,nspin

          read(iu) ik,k(1),k(2),k(3)
          write(io) ik,k(1),k(2),k(3)
          if (ik .ne. iik) stop 'error in index of k-point'
          read(iu) ispin
          write(io) ispin
          if (ispin .ne. iispin) stop 'error in index of spin'
          read(iu) nwflist
          write(io) nwflist

C Loop over wavefunctions 

          do iw = 1,nwflist

            read(iu) indwf
            read(iu) energy
            write(io) indwf
            write(io) energy

            read(iu) (psi(1:,j), j=1,nuotot)

            do jj = 1,nuotot
               if (gamma) then
                  write(io)
     .                 iaorb(jj),labelfis(jj),jj,
     .                 iphorb(jj), cnfigfio(jj),
     .                 symfio(jj), dble(psi(1,jj)), 0.d0
               else
                  write(io)
     .                 iaorb(jj),labelfis(jj),jj,
     .                 iphorb(jj), cnfigfio(jj),
     .                 symfio(jj), dble(psi(1,jj)), dble(psi(2,jj))
               endif
            enddo

          enddo
        enddo
      enddo

      close (iu)
      close (io)
      end
