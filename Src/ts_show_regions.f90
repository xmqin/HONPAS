! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine ts_show_regions(ucell,na_u,xa,N_Elec,Elecs)
  use precision, only : dp
  use units, only : Ang
  use parallel, only : IONode
  use fdf, only : fdf_get

  use m_ts_electype
  use m_ts_method, only : r_aBuf, r_aC, na_Buf, atom_type
  use m_ts_method, only : TYP_BUFFER, TYP_DEVICE
  use m_region
  
  implicit none
! ********************
! * INPUT variables  *
! ********************
  real(dp), intent(in) :: ucell(3,3)
  integer, intent(in)  :: na_u
  real(dp), intent(in) :: xa(3,na_u)
  integer, intent(in) :: N_Elec
  type(Elec), intent(in) :: Elecs(N_Elec)

! ********************
! * LOCAL variables  *
! ********************
  integer :: ia, i, ia_mid
  type(tRgn) :: rgn, rtmp

  if ( .not. IONode ) return

  if ( .not. fdf_get('TS.Atoms.Print',.false.) ) then
     
    write(*,'(/,a)') 'transiesta: Regions of atoms:'
     
    if ( r_aBuf%n > 0 ) then
      call rgn_copy(r_aBuf, rgn)
      call rgn_sort(rgn)
      rgn%name = 'Buffer'
      call rgn_print(rgn, name='//', seq_max=12, indent=3)
    end if

    ! Prepare the device region
    call rgn_copy(r_aC, rtmp)
    do i = 1 , N_Elec
      ia = Elecs(i)%idx_a
      ia_mid = ia + TotUsedAtoms(Elecs(i)) - 1
      call rgn_range(rgn, ia, ia_mid)
      rgn%name = 'Elec.'//trim(Elecs(i)%name)
      call rgn_print(rgn, name='##', seq_max=12, indent=3)
      ! Remove the electrode from the device
      call rgn_remove(rtmp, rgn, rtmp)
    end do
    ! Sort the region
    call rgn_sort(rtmp)
    rtmp%name = 'Device'
    call rgn_print(rtmp, name='--', seq_max=12, indent=3)

    call rgn_delete(rgn, rtmp)
     
    write(*,*) ! new-line

    return
     
  end if

  ! mid-point of device
  ia_mid = (na_u - na_Buf-sum(TotUsedAtoms(Elecs))+1) / 2

  write(*,'(/,a)') 'transiesta: Atomic coordinates and regions (Ang):'
  ia = 1
  do while ( ia <= na_u )
     select case ( atom_type(ia) )
     case ( TYP_BUFFER )
        i = ia
        do while ( i <= na_u )
           ! this assures that we do not have a memory leak
           if ( atom_type(i) /= TYP_BUFFER ) exit
           i = i + 1
        end do
        i = i - ia
        ! steps ia counter
        call out_REGION(ia,i,'Buffer','/')
     case ( TYP_DEVICE )
        ia_mid = ia_mid - 1
        if ( ia_mid == 0 ) then
           write(*,'(tr1,i5,tr2,3(tr2,f12.7),tr8,a)') ia,xa(:,ia)/Ang,'Device'
        else
           write(*,'(tr1,i5,tr2,3(tr2,f12.7))') ia,xa(:,ia)/Ang
        end if
        ia = ia + 1
     case default
        ! electrode position
        do i = 1 , N_Elec
           if ( ia == Elecs(i)%idx_a ) then
              call out_REGION(ia,TotUsedAtoms(Elecs(i)), &
                   trim(name(Elecs(i)))//' electrode','#')
              exit
           end if
        end do
     end select
  end do

  write(*,*)

contains 

  subroutine out_REGION(ia,NA,name,marker,first,last)
    integer, intent(inout) :: ia
    integer, intent(in)    :: NA
    character(len=*), intent(in) :: name
    character(len=1), intent(in) :: marker
    character(len=1), intent(in), optional :: first, last
    integer :: i, mid
    logical :: lfirst, llast
    lfirst = present(first)
    llast = present(last)

    if ( NA < 1 ) return
    mid = (NA+1) / 2
    write(*,'(tr7,a)') repeat(marker,46)
    do i = 1, NA
       if ( i == 1 .and. lfirst ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr1,a1)') &
               ia,marker,xa(:,ia)/Ang,marker,first
       else if ( i == mid ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr5,a)') &
               ia,marker,xa(:,ia)/Ang,marker,trim(name)
       else if ( i == NA .and. llast ) then
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1,tr1,a1)') &
               ia,marker,xa(:,ia)/Ang,marker,last
       else
          write(*,'(tr1,i5,tr1,a1,3(tr2,f12.7),tr2,a1)') &
               ia,marker,xa(:,ia)/Ang,marker
       end if
       ia = ia + 1
    end do
    write(*,'(tr7,a)') repeat(marker,46)

  end subroutine out_REGION

end subroutine ts_show_regions


