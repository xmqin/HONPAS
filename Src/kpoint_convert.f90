!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!
subroutine kpoint_convert(ucell,kin,kout,iopt)
! **********************************************************************
! Enables the conversion between Fourier space k-points into reciprocal
! space k-points.
! Created by Nick Papior Andersen, Aug. 2012
! Modified by Nick Papior Andersen, Jan. 2015
! ***************** INPUT **********************************************
! real*8  ucell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! real*8  kin(3)      : k-point in units of [b] or [1/Bohr]
! integer iopt        : -2 => From units of [b] to [1/Bohr], 
!                             Here 'ucell' is the reciprocal cell with 2Pi
!                             This can be obtained by:
!                                'call reclat(cell,rcell,1)'
!                     : -1 => From units of [b] to [1/Bohr]
!                     :  1 => From units of [1/Bohr] to [b]
! ***************** OUTPUT *********************************************
! real*8  kout(3)     : k-point in units of [b] or [1/Bohr]
!
! Allows for conversion between units of reciprocal k-points.
! **********************************************************************
  use precision, only : dp
  use units    , only : Pi
  use sys      , only : die
  
  real(dp), dimension(3,3), intent(in)  :: ucell
  real(dp), dimension(3)  , intent(in)  :: kin
  real(dp), dimension(3)  , intent(out) :: kout
  integer                 , intent(in)  :: iopt
  
! ***********************
! * LOCAL variables     *
! ***********************
  real(dp), dimension(3,3) :: rcell
  
  if ( iopt == 1 ) then
     kout(1) = sum(ucell(:,1) * kin(:)) * 0.5_dp / Pi
     kout(2) = sum(ucell(:,2) * kin(:)) * 0.5_dp / Pi
     kout(3) = sum(ucell(:,3) * kin(:)) * 0.5_dp / Pi
  else if ( iopt == -1 ) then
     call reclat(ucell,rcell,1)
     kout(1) = sum(rcell(1,:) * kin(:))
     kout(2) = sum(rcell(2,:) * kin(:))
     kout(3) = sum(rcell(3,:) * kin(:))
  else if ( iopt == -2 ) then
     kout(1) = sum(ucell(1,:) * kin(:))
     kout(2) = sum(ucell(2,:) * kin(:))
     kout(3) = sum(ucell(3,:) * kin(:))
  else
     call die("Wrong option for kpoint_convert! Only 1, -1 or -2 allowed.")
  end if

end subroutine kpoint_convert
