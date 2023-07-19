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
!       Converts  a WFS file (old style wf file) to a WFSX file 
!       (well-packed, single precision)
!       Since the WFSX file contains also the k-point weights, it
!       is necessary to have a KP file, with the same format as
!       in SIESTA, to hold that information.
!
!       NOTE: The heuristic nk=1 --> gamma calculation is used
!       to determine whether to use a complex or a real array
!       for the wavefunctions in the WFSX file. Use with care
!       if this is not true...  a crude detection is used.
!
        implicit none

        integer, parameter :: dp = selected_real_kind(14,100)
        integer, parameter :: sp = selected_real_kind(6,30)

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,
     .          nwflist, iw, indwf, j, nuotot, jj

        character fname*33, oname*33
        integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
        character(len=20), allocatable, dimension(:) :: symfio,labelfis

        real(dp), allocatable, dimension(:) :: psir, psii
        real(dp), allocatable, dimension(:) :: wk

        logical gamma, orb_info_written

        real(dp) k(3), energy
        integer :: jj_dummy
 
        iu = 10
        io = 11

        open(iu, file="WFS", form='unformatted', status='old' )
        open(io, file="WFSX", form='unformatted', status='new' )

        rewind (iu)
        read(iu) nk
        read(iu) nspin
        read(iu) nuotot

        allocate(wk(nk))
        gamma = (nk .eq. 1)   ! Suspect
        open(unit=4,file="KP",form="formatted")
        read(4,*) 
        do ik=1, nk
           read(4,*) iik, k(1:3), wk(ik)
        enddo
        close(4)

        rewind (io)
        write(io) nk, gamma
        write(io) nspin
        write(io) nuotot

        allocate(psir(nuotot),psii(nuotot))

        allocate(iaorb(nuotot), labelfis(nuotot),
     $           iphorb(nuotot), cnfigfio(nuotot),
     $           symfio(nuotot))

        orb_info_written = .false.

        do iik = 1,nk
          do iispin = 1,nspin

          read(iu) ik,k(1),k(2),k(3)
          if (ik .ne. iik) stop 'error in index of k-point'
          read(iu) ispin
          if (ispin .ne. iispin) stop 'error in index of spin'
          read(iu) nwflist

          if (orb_info_written) then
             write(io) ik,k(1),k(2),k(3), wk(iik)
             write(io) ispin
             write(io) nwflist
          endif

C Loop over wavefunctions 

          do iw = 1,nwflist

            read(iu) indwf
            read(iu) energy

            do jj = 1,nuotot
               read(iu)
     .                 iaorb(jj),labelfis(jj),jj_dummy,
     .                 iphorb(jj), cnfigfio(jj),
     .                 symfio(jj), psir(jj), psii(jj)
            enddo

            if (.not. orb_info_written) then
               write(io) (iaorb(j),labelfis(j),
     .            iphorb(j), cnfigfio(j),
     .            symfio(j), j=1,nuotot)

               if ( gamma .and. (sum(k(1:3)) /= 0.0)) then
                  STOP "nk=1 but not gamma!"
               endif
               write(io) ik,k(1),k(2),k(3), wk(ik)
               write(io) ispin
               write(io) nwflist
               orb_info_written = .true.
            endif

            write(io) indwf
            write(io) energy
            if (gamma) then
               write(io) (real(psir(j),kind=sp), j=1,nuotot)
            else
               write(io) (real(psir(j),kind=sp),
     $                    real(psii(j),kind=sp), j=1,nuotot)
            endif

          enddo
        enddo
      enddo

      close (iu)
      close (io)
      end
