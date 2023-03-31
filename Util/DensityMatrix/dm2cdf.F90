! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program dm2cdf
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

!
! Converts a DM (binary) file to netCDF format
!
! Usage: dm2cdf
!
!        The input file must be named "DM" (a symbolic link would do)
!        The output file is called "DensityMatrix.nc"
!
!        The program will request the number of orbitals in the auxiliary
!        "interaction supercell". This can be obtained by grepping the OUT
!        file for "orbitals".
!
use netcdf

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)

integer ncid, norbs_id, nspin_id, nnzs_id, scf_step_id
integer no_s_id, indxuo_id
integer numd_id, row_pointer_id, column_id, dm_id

integer   ::    norbs  ! Number of atomic orbitals
integer   ::    no_s   ! Number of orbitals in interaction supercell
integer   ::    nnzs   ! Total number of orbital interactions
integer   ::    nspin  ! Number of spins 

integer   ::    iog, ios, ispin, j, iostat

integer, dimension(:), allocatable  :: numd_global 
integer, dimension(:), allocatable  :: global_row_pointer 
integer, dimension(:), allocatable  :: column
integer, dimension(:), allocatable  :: indxuo
real(dp), dimension(:,:), allocatable :: dm

!-----------------------------------------------------

open(unit=1,file="DM",form="unformatted",status="old",action="read", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File DM cannot be opened"
  STOP
endif

read(1) norbs, nspin
print *, "Norbs, nspin: ", norbs, nspin, ". Enter norbs in supercell: "
read *, no_s

allocate(numd_global(1:norbs),global_row_pointer(1:norbs))

read(1) (numd_global(iog),iog=1,norbs)
nnzs = sum(numd_global(1:norbs))
!
allocate(column(1:nnzs),dm(1:nnzs,1:nspin))

call check( nf90_create('DensityMatrix.nc',NF90_CLOBBER,ncid))
!
!     Dimensions
!
      call check( nf90_def_dim(ncid,'norbs',norbs,norbs_id))  !"Number of basis orbitals"
      call check( nf90_def_dim(ncid,'no_s',no_s,no_s_id))      !"Number of orbitals in supercell"
      call check( nf90_def_dim(ncid,'nspin',nspin,nspin_id))   !"Number of spin components"
      call check( nf90_def_dim(ncid,'nnzs',nnzs,nnzs_id))     !"Number of non-zero interactions"
      call check( nf90_def_dim(ncid,'scf_step',NF90_UNLIMITED,scf_step_id)) !"Index of SCF step"

!
!     Variables
!
      call check( nf90_def_var(ncid,'numd',nf90_int,(/norbs_id/),numd_id))
      call check( nf90_put_att(ncid,numd_id,'Description',"Number of interactions of a given orbital"))
      call check( nf90_def_var(ncid,'row_pointer',nf90_int,(/norbs_id/),row_pointer_id))
      call check( nf90_put_att(ncid,row_pointer_id,'Description',"Index (minus 1) of the start of a given row"))
      call check( nf90_def_var(ncid,'column',nf90_int,(/nnzs_id/),column_id))
      call check( nf90_put_att(ncid,column_id,'Description',"Column index of a given element"))
      call check( nf90_def_var(ncid,'dm',nf90_float,(/nnzs_id,nspin_id,scf_step_id/),dm_id))
      call check( nf90_put_att(ncid,dm_id,'Description',"Density matrix"))

      if (norbs /= no_s) then
         call check( nf90_def_var(ncid,'indxuo',nf90_int,(/no_s_id/),indxuo_id))
         call check( nf90_put_att(ncid,indxuo_id,'Description',"Index of equivalent orb in unit cell"))
      endif

      call check( nf90_enddef(ncid))

!
!     Fill-in  values
!

call check( nf90_put_var(ncid,numd_id,numd_global(1:norbs),count=(/norbs/)))

 ! Compute global row pointer
   global_row_pointer(1) = 0
   do iog=2,norbs
      global_row_pointer(iog) = global_row_pointer(iog-1) + numd_global(iog-1)
   enddo
call check(nf90_put_var(ncid,row_pointer_id,global_row_pointer,count=(/norbs/)))

!
!   Column information
!
do iog = 1, norbs
  read(1) (column(global_row_pointer(iog)+j),j=1,numd_global(iog))
enddo
call check( nf90_put_var(ncid,column_id,column,count=(/nnzs/)))

if (norbs /= no_s) then
   allocate(indxuo(1:no_s))
   do ios = 1, no_s
      indxuo(ios) = mod(ios,norbs)
      if (indxuo(ios) == 0) indxuo(ios) = norbs
   enddo
   call check( nf90_put_var(ncid,indxuo_id,indxuo,count=(/no_s/)))
   deallocate(indxuo)
endif

!
!  DM
!
   do ispin = 1, nspin        
      do iog = 1, norbs
            read(1) (dm(global_row_pointer(iog)+j,ispin),j=1,numd_global(iog))
      enddo
   enddo

   close(1)

   call check( nf90_put_var(ncid, dm_id, dm(1:nnzs,1:nspin),  & 
                              start = (/1, 1, 1 /), &
                              count = (/nnzs, nspin, 1 /) ))

   deallocate(numd_global,global_row_pointer,column,dm)

   call check( nf90_close(ncid))
   print *, "Output in file DensityMatrix.nc"

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
end program dm2cdf
