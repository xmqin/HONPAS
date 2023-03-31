! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!


        program wfsnc2wfsx
#ifndef CDF
        print *, "netCDF not enabled"
#else
    
    use netcdf

!
!       Converts  a WFS.nc file (new_diagk intermediate file) to a WFSX file 
!       (well-packed, single precision)
!       Since the WFSX file contains also the k-point weights, it
!       is necessary to have a KP file, with the same format as
!       in SIESTA, to hold that information.
!
!       Orbital info is meaningless in this version
!       Eigenvalue information newly added to WFS.nc
!
!       Alberto Garcia, March 3, 2009

        implicit none

        integer, parameter :: dp = selected_real_kind(14,100)
        integer, parameter :: sp = selected_real_kind(6,30)

        real(dp), parameter :: eV = 13.60580

        integer  :: ncid 
        integer  :: norbs_id, nspin_id, nk_id, nvectors_id
        integer  :: ncomplex_id, wf_id, eigval_id

        integer io,iu, nk, nspin, ik, iik, ispin, iispin,   &
               nwflist, iw, indwf, j, nuotot, jj, iret, nvectors, n_dummy, &
               ncomplex

        character fname*33, oname*33
        integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
        character(len=20), allocatable, dimension(:) :: symfio,labelfis

        real(sp), allocatable, dimension(:) :: psir, psii
        real(sp), allocatable, dimension(:) :: eigval
        real(dp), allocatable, dimension(:) :: wk
        real(dp), allocatable, dimension(:,:) :: kp

        logical gamma, orb_info_written

        integer :: jj_dummy
 

        call check( nf90_open('WFS.nc',NF90_NOWRITE,ncid))
        call check( nf90_inq_dimid(ncid,'nspin',nspin_id) )
        call check( nf90_inquire_dimension(ncid, dimid=nspin_id, len=nspin) )
        ! The (original)WFS.nc file predates the new convention to support NC/SOC files
        ! Remove this check if its format is updated.
        if (nspin > 2) STOP "This utility does not work for nspin>2"

        call check( nf90_inq_dimid(ncid,'norbs',norbs_id) )
        call check( nf90_inquire_dimension(ncid, dimid=norbs_id, len=nuotot) )
        call check( nf90_inq_dimid(ncid,'nk',nk_id) )
        call check( nf90_inquire_dimension(ncid, dimid=nk_id, len=nk) )
        call check( nf90_inq_dimid(ncid,'ncomplex',ncomplex_id) )
        call check( nf90_inquire_dimension(ncid, dimid=ncomplex_id, len=ncomplex) )
        call check( nf90_inq_dimid(ncid,'nvectors',nvectors_id) )
        call check( nf90_inquire_dimension(ncid, dimid=nvectors_id, len=nvectors) )
        call check( nf90_inq_varid(ncid, "eigval", eigval_id) )
        call check( nf90_inq_varid(ncid, "wf", wf_id) )

        allocate(eigval(nvectors))
        allocate(psir(nuotot),psii(nuotot))

        allocate(wk(nk),kp(3,nk))
        gamma = (ncomplex == 1)
        open(unit=4,file="KP",form="formatted")
        read(4,*)  n_dummy
        if (n_dummy /= nk) STOP "bad nk"
        do ik=1, nk
           read(4,*) iik, kp(1:3,ik), wk(ik)
        enddo
        close(4)

        io = 11
        open(io, file="WFSX", form='unformatted', status='new' )
        rewind (io)
        write(io) nk, gamma
        write(io) nspin
        write(io) nuotot


        allocate(iaorb(nuotot), labelfis(nuotot), &
                iphorb(nuotot), cnfigfio(nuotot), &
                symfio(nuotot))
        iaorb = 1
        labelfis = "pepe"
        iphorb = 1
        cnfigfio = 1
        symfio = "s"
        
        write(io) (iaorb(j),labelfis(j), iphorb(j), cnfigfio(j), &
                 symfio(j), j=1,nuotot)

        do iik = 1,nk
          do iispin = 1,min(4,nspin)

             write(io) iik,kp(1:3,iik), wk(iik)
             write(io) iispin
             write(io) nvectors

!      Loop over wavefunctions 

             iret = nf90_get_var(ncid, eigval_id, eigval,  &
                        start = (/1, 1, 1 /), count = (/nvectors, 1, 1/) )
             call check(iret)

          do iw = 1,nvectors
            write(io) iw
            write(io) real(eigval(iw),kind=dp) * eV

            if (gamma) then
               iret = nf90_get_var(ncid, wf_id, psir,   &
                        start = (/1, 1, iw, iik, iispin /),   &
                        count = (/1, nuotot, 1, 1, 1/) )
               call check(iret)
               write(io) (psir(j), j=1,nuotot)
            else
               iret = nf90_get_var(ncid, wf_id, psir,    &
                        start = (/1, 1, iw, iik, iispin /), &
                        count = (/1, nuotot, 1, 1, 1/) )
               iret = nf90_get_var(ncid, wf_id, psii,    &
                        start = (/2, 1, iw, iik, iispin /),   &
                        count = (/1, nuotot, 1, 1, 1/) )
               call check(iret)
               write(io) (psir(j), psii(j), j=1,nuotot)
            endif

          enddo
        enddo
      enddo

      close (io)
      call check( nf90_close(ncid) )
CONTAINS

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check
#endif
      end program wfsnc2wfsx

