! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! 
! Reprogrammed tbtrans
! This code has been fully created by;
! Nick Papior Andersen, nickpapior @ gmail.com
! The code has been constructed in 2014 and is
! meant to superseede any previous tbtrans versions.

! One major difference in tbtrans is that we never deal with
! xij arrays.
! We ONLY deal with unitcell offsets
! Hence just after reading in TSHS we convert it
! to the correct format (list_col contains the supercell
! index AND the column in the unitcell)
program tbtrans

  use m_tbt_options, only : kT
  use m_tbt_hs, only : TSHS

  use m_tbtrans

  implicit none

  ! Initialize everything.
  call tbt_init()

  ! Call tbtrans
  call tbt(TSHS, kT)

  call tbt_end()

end program tbtrans

