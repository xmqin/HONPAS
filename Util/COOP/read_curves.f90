! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module read_curves

implicit none

public :: read_curve_information, mask_to_arrays

CONTAINS

subroutine read_curve_information(dos,coop,  &
                                     mpr_u,no_u,ncbmx,ncb,tit,orb_mask,dtc)
use precision, only: dp
use orbital_set, only: get_orbital_set

  implicit none

  logical, intent(in)           :: dos, coop
  integer, intent(in)           :: mpr_u, no_u, ncbmx
  integer, intent(out)          :: ncb            ! Total number of curves
  character(len=*), intent(out) :: tit(:)
  logical, intent(out)          :: orb_mask(no_u,2,ncbmx)  ! Orbital sets
  real(dp), intent(out)         :: dtc(ncbmx,2)   ! Distances

  integer :: iostat
  character(len=512) :: line
  real(dp) :: dm

  ncb=0
  do
     read(mpr_u,"(a)",iostat=iostat,end=1300) line
     if (iostat /= 0) then
        print *, "|" // trim(line) // "|"
        stop "Wrong title input"
     endif

     ! Exit if emptly line or if explicit marker is used
     if ((len_trim(line) == 0) .or. (line(1:4) == "----")) EXIT

     ncb=ncb+1
     if (ncb.gt.ncbmx) then
        stop "* ERROR * DOS or COOP curves limit exceeded."
     endif
     tit(ncb)=trim(line)

     ! Get first orbital set for this curve
     read(mpr_u,fmt="(a)",iostat=iostat) line
     if (iostat /= 0) then
        print *, "|" // trim(line) // "|"
        stop "Wrong line input"
     endif
     call get_orbital_set(line,orb_mask(:,1,ncb))

     if (dos) then
        CYCLE
     elseif (coop) then

        ! Get distance range
        read(mpr_u,*,iostat=iostat) dtc(ncb,1:2)
        if (iostat /= 0) then
           print *, dtc(ncb,:)
           stop "Wrong distance input"
        endif
        dm=minval(dtc(ncb,1:2))
        if (dm < 0.0_dp .or. dm /= dtc(ncb,1)) then
           STOP "Bad distance range"
        endif

        ! Get second orbital set for this curve
        read(mpr_u,fmt="(a)",iostat=iostat) line
        if (iostat /= 0) then
           print *, "|" // trim(line) // "|"
           stop "Wrong line input"
        endif
        call get_orbital_set(line,orb_mask(:,2,ncb))

     endif

  enddo

1300 continue

  print *, "Total number of curves processed: ", ncb

end subroutine read_curve_information

subroutine mask_to_arrays(ncb,orb_mask,noc,koc)

  implicit none

  integer, intent(in)          :: ncb            ! Total number of curves
  logical, intent(in)          :: orb_mask(:,:)  ! Orbital sets
  integer, intent(out)         :: noc(:)   ! Number of orbitals in sets
  integer, intent(out)         :: koc(:,:)   ! Distances

  integer i, ic, norbs, no_u

  no_u = size(orb_mask,dim=1)
  noc = 0
  koc = 0
  do ic=1,ncb
     norbs = 0
     do i = 1, no_u
        if (orb_mask(i,ic)) then
           norbs = norbs + 1
           noc(ic) = noc(ic) + 1
           koc(ic,norbs) = i
        endif
     enddo
  enddo
     
end subroutine mask_to_arrays
  

end module read_curves
