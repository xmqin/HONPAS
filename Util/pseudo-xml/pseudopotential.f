! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module pseudopotential

      implicit none

      external :: io_assign, io_close

      private

      integer, parameter  :: dp = selected_real_kind(14)
      
      public :: pseudopotential_t
      public :: pseudo_write_formatted

      type pseudopotential_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
        real(dp)                :: gen_zval  ! Generation valence charge
        logical                 :: relativistic
        character(len=10)       :: correlation
        character(len=2)        :: icorr
        character(len=3)        :: irel
        character(len=4)        :: nicore
        real(dp)                :: a
        real(dp)                :: b
        character(len=10)       :: method(6)
        character(len=70)       :: text
        integer                 :: npotu
        integer                 :: npotd
        real(dp), pointer       :: r(:)
        real(dp), pointer       :: chcore(:)
        real(dp), pointer       :: chval(:)
        real(dp), pointer       :: vdown(:,:)
        real(dp), pointer       :: vup(:,:)
        integer, pointer        :: ldown(:)
        integer, pointer        :: lup(:)
      end type pseudopotential_t

        CONTAINS

!----
        subroutine pseudo_write_formatted(fname,p,print_gen_zval)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t), intent(in)     :: p
        logical, intent(in), optional :: print_gen_zval
        integer io_ps, i, j

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown',
     $       action="write",position="rewind")
        write(6,'(3a)') 'Writing pseudopotential information ',
     $       'in formatted form to ', trim(fname)

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,4g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

        write(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
        write(io_ps,8010) (p%method(i),i=1,6), p%text
        if (present(print_gen_zval)) then
           if (print_gen_zval) then
              write(io_ps,8015) p%npotd, p%npotu, p%nr,
     $             p%b, p%a, p%zval, p%gen_zval
           else
              write(io_ps,8015) p%npotd, p%npotu,
     $                          p%nr, p%b, p%a, p%zval
           endif
        else
           write(io_ps,8015) p%npotd, p%npotu, p%nr,
     $                       p%b, p%a, p%zval
        endif

        write(io_ps,8040) "Radial grid follows"
        write(io_ps,8030) (p%r(j),j=2,p%nrval)

        do i=1,p%npotd
           write(io_ps,8040)
     $          "Down Pseudopotential follows (l on next line)"
           write(io_ps,8000) p%ldown(i)
           write(io_ps,8030)
     $          (force_underflow(p%vdown(i,j)), j=2,p%nrval)
        enddo

        do i=1,p%npotu
           write(io_ps,8040)
     $          "Up Pseudopotential follows (l on next line)"
           write(io_ps,8000) p%lup(i)
           write(io_ps,8030)
     $          (force_underflow(p%vup(i,j)), j=2,p%nrval)
        enddo

        write(io_ps,8040) "Core charge follows"
        write(io_ps,8030) (force_underflow(p%chcore(j)),j=2,p%nrval)
        write(io_ps,8040) "Valence charge follows"
        write(io_ps,8030) (force_underflow(p%chval(j)),j=2,p%nrval)

        call io_close(io_ps)
        end subroutine pseudo_write_formatted
!--------
!
      function force_underflow(x) result(res)
      real(dp), intent(in) ::  x
      real(dp)             ::  res
C
C     Avoid very small numbers that might need a three-character
C     exponent field in formatted output
C      
      if (abs(x) .lt. 1.0e-99_dp) then
         res = 0.0_dp
      else
         res = x
      endif

      end function force_underflow

      end module pseudopotential



