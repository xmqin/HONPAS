      subroutine get_unit(lun)

C     Get an available Fortran unit number

      integer lun

      integer i
      logical unit_used

      do i = 10, 99
         lun = i
         inquire(lun,opened=unit_used)
         if (.not. unit_used) return
      enddo
      stop 'NO LUNS'
      end
