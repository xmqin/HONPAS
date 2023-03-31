! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program lwf2cdf
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

!
! Converts a LWF (binary) file to netCDF format
!
! Usage: lwf2cdf
!
!        The input file must be named "LWF" (a symbolic link would do)
!        The output file is called "WannierFunctions.nc"
!
use netcdf

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)

integer ncid, norbs_id, nspin_id, nnzs_id, nbands_id
integer numc_id, listc_id, c_id

integer   ::    norbs  ! Number of atomic orbitals
integer   ::    nnzs   ! Maximum number of WF (bands) in which a given orb
                       ! has a non-zero coefficient
integer   ::    nspin  ! Number of spins 
integer   ::    nbands ! Number of bands

integer   ::    iog, ios, ispin, j, iostat

integer, dimension(:), allocatable  :: numc
integer, dimension(:,:), allocatable  :: listc
real(dp), dimension(:,:,:), allocatable :: c

!-----------------------------------------------------

open(unit=1,file="LWF",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File LWF cannot be opened"
  STOP
endif

read(1) norbs, nspin
print *, "Norbs, nspin: ", norbs, nspin


allocate(numc(1:norbs))
read(1) (numc(iog),iog=1,norbs)
!
nnzs = maxval(numc(1:norbs))
allocate(listc(1:nnzs,1:norbs))

!
!   Index information
!
listc(1:nnzs,1:norbs) = 0
do iog = 1, norbs
  read(1) (listc(j,iog),j=1,numc(iog))
enddo
nbands = maxval(listc)

call check( nf90_create('WannierFunctions.nc',NF90_CLOBBER,ncid))
!
!     Dimensions
!
      call check( nf90_def_dim(ncid,'norbs',norbs,norbs_id))  !"Number of basis orbitals"
      call check( nf90_def_dim(ncid,'nspin',nspin,nspin_id))   !"Number of spin components"
      call check( nf90_def_dim(ncid,'nnzs',nnzs,nnzs_id))     !"Max. number of bands touching an orb"
      call check( nf90_def_dim(ncid,'nbands',nbands,nbands_id))  !"Number of bands"
!
!     Variables

      call check( nf90_def_var(ncid,'numc',nf90_int,(/norbs_id/),numc_id))
      call check( nf90_put_att(ncid,numc_id,'Description',"Number of bands touching a given orb"))
      call check( nf90_def_var(ncid,'listc',nf90_int,(/nnzs_id,norbs_id/),listc_id))
      call check( nf90_put_att(ncid,listc_id,'Description',"List of bands touching a given orb"))
      call check( nf90_def_var(ncid,'c',nf90_float,(/nnzs_id,norbs_id,nspin_id/),c_id))
      call check( nf90_put_att(ncid,c_id,'Description',"Wannier Function coefficients"))

      call check( nf90_enddef(ncid))
!
!     Fill-in  values
!
call check( nf90_put_var(ncid,numc_id,numc(1:norbs),count=(/norbs/)))
call check( nf90_put_var(ncid,listc_id,listc,count=(/nnzs,norbs/)))

deallocate(listc)
!
!  C
!
   allocate(c(1:nnzs,1:norbs,1:nspin))
   c(:,:,:) = 0.0_dp
   do ispin = 1, nspin        
      do iog = 1, norbs
            read(1) (c(j,iog,ispin),j=1,numc(iog))
      enddo
   enddo

   close(1)

   call check( nf90_put_var(ncid, c_id, c(1:nnzs,1:norbs,1:nspin),  & 
                              start = (/1, 1, 1 /), &
                              count = (/nnzs, norbs, nspin/) ))

   deallocate(numc,c)

   call check( nf90_close(ncid))
   print *, "Output in file WannierFunctions.nc"

CONTAINS

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check

#endif /* CDF */
end program lwf2cdf
