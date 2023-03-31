! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module orbital_set
CONTAINS
subroutine get_orbital_set(line,set_mask)
  use main_vars
  use subs, only : txt2wrd, orbital

  implicit none

  character(len=*), intent(in) :: line
  logical, intent(out)         :: set_mask(:)

  print *, "Size of set_mask: ", size(set_mask)
  set_mask(:) = .false.

  call txt2wrd (line, wrd, nwd, nlwmx)
  if (nwd.gt.nlwmx) stop "* Groups per subset limit exceeded."

  if (trim(wrd(1)).eq.'+') then
     do iw=2,nwd
        k=ival(wrd(iw))
        if (k.le.0.or.k.gt.no_u) then
           print *, "Wrong orbital number: ", k
           STOP
        endif
        set_mask(k) = .true.
     enddo

  else

     do iw=1,nwd
        call orbital (wrd(iw), ia, cx, n, l, k)

        if (ia.lt.0) then
           print *, "Error in orb spec: ", trim(wrd(iw))
           STOP
        endif

        if (ia.eq.0) then
           it=0
           do i=1,nspecies
              if (trim(label(i)) .eq. trim(cx)) it=i
           enddo
           if (it.eq.0) then
              print *, "Wrong species: ", trim(cx)
              STOP
           endif
        endif
        if (ia > na_u) then
           print *, "Atom index too big: ", ia
           STOP
        endif
        !!! Could go on checking whether a given
        !!! atom has the orbitals specified, instead
        !!! or giving an empty result
        !
        !             See which orbitals match
        !
        do io=1,no_u
           if ((za(io).eq.ia).or.(ia.eq.0.and.zc(io).eq.it)) then
              if ((zn(io).eq.n.or.n.eq.-1).and. &
                   (zl(io).eq.l.or.l.eq.-1).and. &
                   (zx(io).eq.k.or.k.eq.-1)) then
                 set_mask(io) = .true.
              endif
           endif
        enddo
     enddo

  endif

end subroutine get_orbital_set

end module orbital_set

