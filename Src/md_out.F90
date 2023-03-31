! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module md_out

use units, only: Ang
use precision, only: dp
use sys,      only: die
use files,    only: slabel
use m_energies, only: etot

implicit none

public :: md_v_format
#ifdef CDF
public :: md_netcdf
#endif

CONTAINS

subroutine md_v_format(na,isa,xa,cell)

!
! For V-compatible "movie" output
!
integer, intent(in)                 :: na
integer, dimension(na), intent(in)   :: isa 
real(dp), dimension (3,na), intent(in) :: xa
real(dp), dimension (3,3), intent(in)  :: cell

real(dp), dimension (3,3)  :: celli
real(dp), dimension (3)  :: frac
integer, dimension(:), allocatable, save   :: natoms

integer, save :: iomd
logical, save :: first = .true.

integer  :: nspecies, i
external :: reclat, io_assign

if (first) then

    call io_assign(iomd)
    open(unit=iomd, file=trim(slabel)//".MD_CAR", &
         form='formatted', position='append', &
         action="write", status='unknown')

       nspecies = maxval(isa(1:na))
       allocate(natoms(nspecies))
       natoms(:) = 0
       do i = 1,na
          natoms(isa(i)) = natoms(isa(i)) + 1
       enddo

     first = .false.
endif

write(iomd,"(a)") "---" // trim(slabel) //"---"
write(iomd,"(f10.1)") 1.0
do i=1, 3
   write(iomd,"(3f16.9)") cell(:,i)/Ang
enddo
write(iomd,"(30i6)") natoms(:)
write(iomd,"(a)") "Direct"

call reclat(cell,celli,0)
       do i=1, na
          frac(:) =  matmul(transpose(celli),xa(:,i))
          write(iomd,"(3f16.9)") frac(:)
       enddo
call pxfflush(iomd)

end subroutine md_v_format

#ifdef CDF
!--------------------------------------------------------

subroutine md_netcdf(na,isa,iza,xa,va,cell,vcell,varcel, &
                     temp, eks, tot_energy, volume, Psol)

use netcdf

integer, intent(in)                 :: na
integer, dimension(na), intent(in)   :: isa, iza 
real(dp), dimension (3,na), intent(in) :: xa, va
real(dp), dimension (3,3), intent(in)  :: cell, vcell
logical, intent(in)                    :: varcel
real(dp), intent(in)               :: temp, eks, tot_energy
real(dp), intent(in)               :: volume, Psol


integer iret
integer  :: ncid 
integer  :: xyz_id, atom_id, step_id, abc_id
integer  :: xa_id, va_id, cell_id, vcell_id
integer  :: eks_id, etot_id, temp_id, psol_id
integer  :: isa_id, iza_id
integer  :: volume_id

integer, save  :: step_no

integer        :: atom_no

!
!  Attempt to open. If the file does not exist, create and initialize it
!
    iret = nf90_open(trim(slabel)//".MD.nc",NF90_SHARE+NF90_WRITE,ncid)

    if (iret /= nf90_noerr) then

       ! Need to create the dataset first...

       iret = nf90_create(trim(slabel)//".MD.nc",NF90_SHARE+NF90_WRITE,ncid)
       call check(iret)
       iret = nf90_def_dim(ncid,'xyz',3,xyz_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'abc',3,abc_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'atom',na,atom_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'step',NF90_UNLIMITED,step_id)
       call check(iret)
       step_no = 0

       iret = nf90_def_var(ncid,'isa',nf90_int,(/atom_id/),isa_id)
       call check(iret)
       iret = nf90_put_att(ncid,isa_id,'Description',"Species index")
       call check(iret)
       iret = nf90_def_var(ncid,'iza',nf90_int,(/atom_id/),iza_id)
       call check(iret)
       iret = nf90_put_att(ncid,iza_id,'Description',"Atomic number")
       call check(iret)

       iret = nf90_def_var(ncid,'temp',nf90_double,(/step_id/),temp_id)
       call check(iret)
       iret = nf90_put_att(ncid,temp_id,'Description',"Temperature in K")
       call check(iret)

       iret = nf90_def_var(ncid,'psol',nf90_double,(/step_id/),psol_id)
       call check(iret)
       iret = nf90_put_att(ncid,psol_id,'Description',"Pressure in Kbar")
       call check(iret)

       iret = nf90_def_var(ncid,'eks',nf90_double,(/step_id/),eks_id)
       call check(iret)
       iret = nf90_put_att(ncid,eks_id,'Description',"Kohn-Sham energy in ")
       call check(iret)

       iret = nf90_def_var(ncid,'etot',nf90_double,(/step_id/),etot_id)
       call check(iret)
       iret = nf90_put_att(ncid,etot_id,'Description',"Total energy in")
       call check(iret)

       iret = nf90_def_var(ncid,'xa',nf90_double,(/xyz_id,atom_id,step_id/),xa_id)
       call check(iret)
       iret = nf90_put_att(ncid,xa_id,'Description', &
             "Atomic coordinates in cartesian Angstrom: xyz, ia, step")
       call check(iret)

       iret = nf90_def_var(ncid,'va',nf90_double,(/xyz_id,atom_id,step_id/),va_id)
       call check(iret)
       iret = nf90_put_att(ncid,va_id,'Description', &
             "Atomic velocities in cartesian Angstrom/time?: xyz, ia, step")
       call check(iret)

          iret = nf90_def_var(ncid,'cell',nf90_double,(/xyz_id,abc_id,step_id/),cell_id)
       call check(iret)
          iret = nf90_put_att(ncid,cell_id,'Description', &
               "Variable cell vectors in Ang: xyz, abc, step")
       call check(iret)

          iret = nf90_def_var(ncid,'vcell',nf90_double,(/xyz_id,abc_id,step_id/),vcell_id)
       call check(iret)
          iret = nf90_put_att(ncid,vcell_id,'Description', &
             "Cell vectors' velocities in Ang/time?: xyz, abc, step")
       call check(iret)

          iret = nf90_def_var(ncid,'volume',nf90_double,(/step_id/),volume_id)
       call check(iret)
          iret = nf90_put_att(ncid,volume_id,'Description', &
               "Cell volume in Ang**3")
       call check(iret)

       iret = nf90_enddef(ncid)
       call check(iret)
!
!      Put unchanging stuff
!
       iret = nf90_put_var(ncid, isa_id, isa(1:na), start = (/ 1 /) )
       call check(iret)
       iret = nf90_put_var(ncid, iza_id, iza(1:na), start = (/ 1 /) )
       call check(iret)

    endif  ! File creation

! Now append the interesting information

       iret = nf90_inq_dimid(ncid, "step", step_id)
       iret = nf90_inq_dimid(ncid, "xyz", xyz_id)
       iret = nf90_inq_dimid(ncid, "abc", abc_id)
       iret = nf90_inq_dimid(ncid, "atom", atom_id)
       iret = nf90_inquire_dimension(ncid, dimid=step_id, len=step_no)
       call check(iret)
       iret = nf90_inquire_dimension(ncid, dimid=atom_id, len=atom_no)
       call check(iret)

       if (atom_no /= na) call die("Wrong number of atoms in existing NC set!!")

       iret = nf90_inq_varid(ncid, "temp", temp_id)
       call check(iret)
       iret = nf90_inq_varid(ncid, "psol", psol_id)
       call check(iret)
       iret = nf90_inq_varid(ncid, "eks", eks_id)
       call check(iret)
       iret = nf90_inq_varid(ncid, "etot", etot_id)
       call check(iret)

       iret = nf90_inq_varid(ncid, "volume", volume_id)
       call check(iret)

       iret = nf90_inq_varid(ncid, "xa", xa_id)
       call check(iret)
       iret = nf90_inq_varid(ncid, "va", va_id)
       call check(iret)

       iret = nf90_inq_varid(ncid, "cell", cell_id)
       call check(iret)
       iret = nf90_inq_varid(ncid, "vcell", vcell_id)
       call check(iret)

    step_no = step_no + 1

! put values
    iret = nf90_put_var(ncid, temp_id, temp, start = (/ step_no /) )
       call check(iret)

    iret = nf90_put_var(ncid, psol_id, psol, start = (/ step_no /) )
       call check(iret)

    iret = nf90_put_var(ncid, eks_id, eks, start = (/ step_no /) )
       call check(iret)

    iret = nf90_put_var(ncid, etot_id, etot, start = (/ step_no /) )
       call check(iret)

    iret = nf90_put_var(ncid, volume_id, volume, start = (/ step_no /) )
       call check(iret)

    iret = nf90_put_var(ncid, xa_id, xa, start = (/1, 1, step_no /), &
                        count = (/3, na, 1 /) )
       call check(iret)

    iret = nf90_put_var(ncid, va_id, va, start = (/1, 1, step_no /), &
                        count = (/3, na, 1 /) )
       call check(iret)

    iret = nf90_put_var(ncid, cell_id, cell, start = (/1, 1, step_no /), &
                        count = (/3, 3, 1 /) )
       call check(iret)

    iret = nf90_put_var(ncid, vcell_id, vcell, start = (/1, 1, step_no /), &
                        count = (/3, 3, 1 /) )
       call check(iret)
!
!   Close file at every step, to avoid leaving it in an undefined
!   state in the event of a crash
!
    iret = nf90_close(ncid)
    call check(iret)

end subroutine md_netcdf

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) call die("netCDF error: " // NF90_STRERROR(code))
end subroutine check

#endif

end module md_out

       
