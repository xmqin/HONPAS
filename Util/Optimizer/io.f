      subroutine io_assign(lun)
      integer, intent(out) :: lun
      logical used
      integer iostat
c
c     Looks for a free unit and assigns it to lun
c
      do lun= 10, 99
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            if (.not. used) return
      enddo
      end subroutine io_assign







