! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
!  Nick Papior, 2015
module m_ts_io_version

  use precision, only : dp
  implicit none

contains

  subroutine write_TSHS_0(filename, &
       onlyS, Gamma, TSGamma, &
       ucell, na_u, no_l, no_u, no_s, maxnh, nspin,  &
       kscell, kdispl, &
       xa, iza, lasto, &
       numh, listhptr, listh, xij, &
       H, S, Ef, &
       Qtot, Temp, &
       istep, ia1)

    use geom_helper, only : ucorb

! **********************
! * INPUT variables    *
! **********************
    character(len=*), intent(in) :: filename
    logical, intent(in) :: onlyS
    logical, intent(in) :: Gamma, TSGamma
    real(dp), intent(in) :: ucell(3,3)
    integer, intent(in) :: na_u, no_l, no_u, no_s, maxnh, nspin
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: iza(na_u)
    integer, intent(in) :: numh(no_l), listhptr(no_l)
    integer, intent(in) :: listh(maxnh)
    real(dp), intent(in) :: xij(3,maxnh)
    integer, intent(in) :: lasto(0:na_u)
    real(dp), intent(in) :: H(maxnh,nspin), S(maxnh)
    real(dp), intent(in) :: Ef
    integer, intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)
    real(dp), intent(in) :: Qtot,Temp
    integer, intent(in) :: istep, ia1
    
! ************************
! * LOCAL variables      *
! ************************
    integer :: indxuo(no_s)
    integer :: iu
    integer :: ispin, i,j

    external :: io_assign, io_close

    ! Open file
    call io_assign( iu )
    open( iu, file=filename, form='unformatted', status='unknown' )

    ! Write Dimensions Information
    write(iu) na_u, no_u, no_s, nspin, maxnh
       
    ! Write Geometry information
    write(iu) xa
    write(iu) iza
    write(iu) ucell
    
    ! Write k-point samplung information
    write(iu) Gamma
    ! Is this only an S containing object?
    write(iu) onlyS
    write(iu) TSGamma
    write(iu) kscell
    write(iu) kdispl
    write(iu) istep, ia1
    
    write(iu) lasto
    
    if ( .not. Gamma ) then
       do i = 1 , no_s
          indxuo(i) = ucorb(i,no_u)
       end do
       write(iu) indxuo
    endif

    write(iu) numh

    ! Write Electronic Structure Information
    write(iu) Qtot,Temp
    write(iu) Ef
    
    ! Write listh
    do i = 1 , no_u
       write(iu) listh(listhptr(i)+1:listhptr(i)+numh(i))
    end do

    ! Write Overlap matrix
    do i = 1 , no_u
       write(iu) S(listhptr(i)+1:listhptr(i)+numh(i))
    end do
    
    if ( .not. onlyS ) then
       ! Write Hamiltonian 
       do ispin = 1 , nspin 
          do i = 1 , no_u
             write(iu) H(listhptr(i)+1:listhptr(i)+numh(i),ispin)
          end do
       end do
    end if

    if ( .not. Gamma ) then
       do i = 1 , no_u
          write(iu) (xij(j,listhptr(i)+1:listhptr(i)+numh(i)),j=1,3)
       end do
    end if

    ! Close file
    call io_close( iu )

  end subroutine write_TSHS_0
 
end module m_ts_io_version
    

