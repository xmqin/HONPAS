      module sys
      public :: die
      CONTAINS
      subroutine die(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      stop
      end subroutine die
      end module sys

