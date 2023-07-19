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

      use sys, only: die
      use precision, only: dp
      use flib_spline, only: generate_spline, evaluate_spline
      use atom_options, only: write_ion_plot_files
      
      implicit none

      external :: io_assign, io_close

      private

      public :: pseudopotential_t, pseudo_read, pseudo_header_print
      public :: pseudo_write_formatted, pseudo_reparametrize
      public :: read_ps_conf, pseudo_dump

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

        subroutine pseudo_read(label,p)
        character(len=*), intent(in)   :: label
        type(pseudopotential_t)                    :: p

!       PS information can be in a .vps file (unformatted)
!       or in a .psf file (formatted)

        character(len=40) fname
        logical found

        call vps_init(p)

        fname  = trim(label) // '.vps'
        inquire(file=fname, exist=found)
        if (found) then
           call pseudo_read_unformatted(fname,p)
        else
           fname = trim(label) // '.psf'
           inquire(file=fname, exist=found)
           if (found) then
              call pseudo_read_formatted(fname,p)
           else
              write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',
     .             'Pseudopotential file not found: ', fname
              call die
           endif
        endif
        if (write_ion_plot_files)
     $       call pseudo_dump(trim(label) // ".psdump",p)
        end subroutine pseudo_read
!
        subroutine pseudo_read_unformatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j
        real(dp) :: r2

        call io_assign(io_ps)
        open(io_ps,file=fname,form='unformatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in unformatted form from ', trim(fname)

        read(io_ps) p%name, p%icorr, p%irel, p%nicore,
     .       (p%method(i),i=1,6), p%text,
     .       p%npotd, p%npotu, p%nr, p%b, p%a, p%zval

!
!       Old style vps files should have the right info in text.
!
        call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)

        p%nrval = p%nr + 1
        allocate(p%r(1:p%nrval))
        read(io_ps) (p%r(j),j=2,p%nrval)
        p%r(1) = 0.d0

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif
        do i=1,p%npotd
           read(io_ps) p%ldown(i), (p%vdown(i,j), j=2,p%nrval)
           p%vdown(i,1) = p%vdown(i,2)
        enddo

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
        do i=1,p%npotu
           read(io_ps) p%lup(i), (p%vup(i,j), j=2,p%nrval)
           p%vup(i,1) = p%vup(i,2)
        enddo

        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        read(io_ps) (p%chcore(j),j=2,p%nrval)
        read(io_ps) (p%chval(j),j=2,p%nrval)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))

        call io_close(io_ps)
        end subroutine pseudo_read_unformatted
!----
        subroutine pseudo_read_formatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p

        integer io_ps, i, j, ios
        character(len=70) dummy
        real(dp) :: r2, gen_zval_inline

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in formatted form from ', trim(fname)

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,4g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

        read(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
        read(io_ps,8010) (p%method(i),i=1,6), p%text
        read(io_ps,8015,iostat=ios)
     $       p%npotd, p%npotu, p%nr, p%b, p%a, p%zval,
     $       gen_zval_inline
        if (ios < 0) gen_zval_inline = 0.0_dp
        call read_ps_conf(p%irel,p%npotd-1,p%text,p%gen_zval)
!
!       (Some .psf files contain an extra field corresponding
!       to the ps valence charge at generation time. If that
!       field is not present, the information has to be decoded
!       from the "text" variable.
!
!       "Zero" pseudos have gen_zval = 0, so they need a special case.

        if (p%gen_zval == 0.0_dp) then
           if (gen_zval_inline == 0.0_dp) then
              if (p%method(1) /= "ZEROPSEUDO")
     $             call die("Cannot get gen_zval")
           else
              p%gen_zval = gen_zval_inline
           endif
        endif

        p%nrval = p%nr + 1
        allocate(p%r(1:p%nrval))
        read(io_ps,8040) dummy
        read(io_ps,8030) (p%r(j),j=2,p%nrval)
        p%r(1) = 0.d0

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif
        do i=1,p%npotd
           read(io_ps,8040) dummy 
           read(io_ps,8000) p%ldown(i)
           read(io_ps,8030) (p%vdown(i,j), j=2,p%nrval)
           p%vdown(i,1) = p%vdown(i,2)
        enddo

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
        do i=1,p%npotu
           read(io_ps,8040) dummy 
           read(io_ps,8000) p%lup(i)
           read(io_ps,8030) (p%vup(i,j), j=2,p%nrval)
           p%vup(i,1) = p%vup(i,2)
        enddo

        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        read(io_ps,8040) dummy
        read(io_ps,8030) (p%chcore(j),j=2,p%nrval)
        read(io_ps,8040) dummy
        read(io_ps,8030) (p%chval(j),j=2,p%nrval)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))

        call io_close(io_ps)
        end subroutine pseudo_read_formatted
!------
!----
        subroutine pseudo_write_formatted(fname,p)
        character(len=*), intent(in) :: fname
        type(pseudopotential_t), intent(in)     :: p

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
        write(io_ps,8015) p%npotd, p%npotu, p%nr, p%b, p%a, p%zval,
     $                    p%gen_zval

        write(io_ps,8040) "Radial grid follows"
        write(io_ps,8030) (p%r(j),j=2,p%nrval)

        do i=1,p%npotd
           write(io_ps,8040) "Down potential follows (l on next line)"
           write(io_ps,8000) p%ldown(i)
           write(io_ps,8030) (p%vdown(i,j), j=2,p%nrval)
        enddo

        do i=1,p%npotu
           write(io_ps,8040) "Up potential follows (l on next line)"
           write(io_ps,8000) p%lup(i)
           write(io_ps,8030) (p%vup(i,j), j=2,p%nrval)
        enddo

        write(io_ps,8040) "Core charge follows"
        write(io_ps,8030) (p%chcore(j),j=2,p%nrval)
        write(io_ps,8040) "Valence charge follows"
        write(io_ps,8030) (p%chval(j),j=2,p%nrval)

        call io_close(io_ps)
        end subroutine pseudo_write_formatted
!--------
        subroutine pseudo_dump(fname,p)
!
!       Column-oriented output
!
        character(len=*), intent(in) :: fname
        type(pseudopotential_t), intent(in)     :: p

        integer io_ps, i, j

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown',
     $       action="write",position="rewind")
        write(6,'(3a)') 'Dumping pseudopotential information ',
     $       'in formatted form in ', trim(fname)

 9040    format(i6,1x,7f12.6)
         do j = 1, p%nrval
            write(io_ps,9040) j, p%r(j), (p%vdown(i,j),i=1,p%npotd),
     $                        p%chval(j), p%chcore(j)
         enddo
         call io_close(io_ps)
         end subroutine pseudo_dump

!------

        subroutine vps_init(p)
        type(pseudopotential_t)  :: p
        nullify(p%lup,p%ldown,p%r,p%chcore,p%chval,p%vdown,p%vup)
        end subroutine vps_init

!-------
        subroutine pseudo_header_print(lun,p)
        integer, intent(in) :: lun
        type(pseudopotential_t)  :: p

        integer :: i

 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
        
        write(lun,'(a)') '<pseudopotential_header>'
        write(lun,8005) p%name, p%icorr, p%irel, p%nicore
        write(lun,8010) (p%method(i),i=1,6), p%text
        write(lun,'(a)') '</pseudopotential_header>'

        end subroutine pseudo_header_print

c$$$        subroutine pseudo_header_string(p,s)
c$$$        type(pseudopotential_t)  :: p
c$$$        character(len=*), intent(inout) :: s
c$$$
c$$$        integer :: n, i
c$$$
c$$$ 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
c$$$ 8010   format(1x,6a10,/,1x,a70)
c$$$ 8015   format(1x,2i3,i5,3g20.12)
c$$$        
c$$$        write(s,8005) p%name, p%icorr, p%irel, p%nicore
c$$$        n = len_trim(s) + 1
c$$$        write(s(n:),fmt="(a1)") char(10)
c$$$        n = len_trim(s) + 1
c$$$        write(s(n:),8010) (p%method(i),i=1,6), p%text
c$$$
c$$$        end subroutine pseudo_header_string
!--------
!
      subroutine read_ps_conf(irel,lmax,text,chgvps)
!
!     Attempt to decode the valence configuration used for
!     the generation of the pseudopotential
!     (At least, the valence charge)

      character(len=3), intent(in)  :: irel
      integer, intent(in)           :: lmax
      character(len=70), intent(in) :: text
      real(dp), intent(out)         :: chgvps

      integer  :: l, itext
      real(dp) :: ztot, zup, zdown, rc_read
      character(len=2) :: orb

      chgvps=0.0_dp

            if(irel.eq.'isp') then
               write(6,'(/,2a)')
     .          'Pseudopotential generated from an ',
     .          'atomic spin-polarized calculation'

               write(6,'(/,a)') 'Valence configuration '//
     .                 'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8080)
     $                 orb, zdown, zup, rc_read
 8080             format(a2,f4.2,1x,f4.2,1x,f4.2)
                  chgvps = chgvps + zdown + zup
                  write(6,8085) orb, zdown, zup, rc_read
 8085             format(a2,'(',f4.2,',',f4.2,') rc: ',f4.2)
               enddo

            else
               if(irel.eq.'rel') then
                  write(6,'(/,2a)')
     .          'Pseudopotential generated from a ',
     .                 'relativistic atomic calculation'
                  write(6,'(2a)')
     .          'There are spin-orbit pseudopotentials ',
     .                 'available'
                  write(6,'(2a)')
     .          'Spin-orbit interaction is not included in ',
     .                 'this calculation'
               endif

               write(6,'(/,a)') 'Valence configuration '//
     .                 'for pseudopotential generation:'

               do l=0,min(lmax,3)
                  itext=l*17
                  read(text(itext+1:),err=5000,fmt=8090)
     $                 orb, ztot, rc_read
 8090             format(a2,f5.2,4x,f5.2)
                  chgvps = chgvps + ztot
                  write(6,8095) orb, ztot, rc_read
 8095             format(a2,'(',f5.2,') rc: ',f4.2)
               enddo

           endif
           return

 5000    continue       ! Error return: set chgvps to zero

         end subroutine read_ps_conf
!
!***********---------------------------------------------
!
         subroutine pseudo_reparametrize(p,a,b,label)
         use alloc, only: re_alloc, de_alloc
!
!        Interpolate values into new grid, given by a and b
!
!        Typical new values:  a = 5x10-4, b=10

         type(pseudopotential_t)          :: p
         real(dp), intent(in)             :: a, b
         character(len=*)                 :: label

         real(dp)  :: rmax, rpb, ea, ea2, rr
         integer   :: ir, new_nrval, i, j
         real(dp), dimension(:), pointer   :: func, tmp, new_r
         real(dp), dimension(:,:), pointer :: tmp2

         real(dp), dimension(:), pointer   :: y2 => null()

         rmax = p%r(p%nrval)
         print *, "Reparametrization. rmax: ", rmax
         rpb=b
         ea=exp(a)
         ea2=1.0d0
         ir = 0
         do 
            rr = b*(ea2-1.0d0)
            if (rr > rmax) exit
            ir = ir + 1
            rpb=rpb*ea
            ea2=ea2*ea
         enddo
         new_nrval = ir
         allocate(new_r(new_nrval))   ! Will go in derived type
         print *, "Reparametrization. New nrval: ", new_nrval

         rpb=b
         ea=exp(a)
         ea2=1.0d0
         do ir = 1, new_nrval
            new_r(ir) = b*(ea2-1.0d0)
            rpb=rpb*ea
            ea2=ea2*ea
         enddo
         
        call re_alloc( y2, 1, p%nrval, name='y2', 
     &                 routine='pseudo_reparametrize' )
!-----------------------------------------------------------------------
!       Basic idiom to reparametrize
!       (The alloc module is not used here, as storage is
!        linked to derived type)
!       Use natural spline (zero second derivative)
!
        func => p%chcore
        call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
        allocate(tmp(new_nrval))
        do j = 1, new_nrval
           call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp(j))
        enddo
        nullify(func)
        deallocate(p%chcore)    ! Old data
        p%chcore => tmp         ! Point to new memory area
        nullify(tmp)            ! To re-use tmp
!--------------------------------------------------------------------
        func => p%chval
        call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
        allocate(tmp(new_nrval))
        do j = 1, new_nrval
           call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp(j))
        enddo
        nullify(func)
        deallocate(p%chval)    ! Old data
        p%chval => tmp         ! Point to new memory area
        nullify(tmp)           ! To re-use tmp
        
!
!       Careful with 2D arrays...
!
        allocate(tmp2(p%npotd,new_nrval))
        do i=1,p%npotd
           func => p%vdown(i,:)
           call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
           do j = 1, new_nrval
            call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp2(i,j))
           enddo
           nullify(func)
        enddo
        deallocate(p%vdown)      ! Old data
        p%vdown => tmp2          ! Point to new memory area
        nullify(tmp2)            ! To re-use tmp

        if (p%npotu > 0) allocate(tmp2(p%npotu,new_nrval))
        do i=1,p%npotu         ! Only executed if npotu > 0 ...
           func => p%vup(i,:)
           call generate_spline(p%r,func,p%nrval,y2,0.0_dp,0.0_dp)
           do j = 1, new_nrval
            call evaluate_spline(p%r,func,y2,p%nrval,new_r(j),tmp2(i,j))
           enddo
           nullify(func)
        enddo
        if (p%npotu > 0) then
           deallocate(p%vup  )  ! Old data
           p%vup => tmp2        ! Point to new memory area
           nullify(tmp2)        ! To re-use tmp
        endif

!
!       Now re-set the values
!
        deallocate(p%r)
        p%r => new_r
        p%nrval = new_nrval
        p%nr    = p%nrval - 1
        p%a     = a
        p%b     = b

        call de_alloc( y2, name='y2', routine='pseudo_reparametrize' )

        call pseudo_write_formatted(trim(label)// ".Reparam.psf",p)
        if (write_ion_plot_files) then
           call pseudo_dump(trim(label) // ".Reparam.psdump",p)
        endif

      end subroutine pseudo_reparametrize

      end module pseudopotential



