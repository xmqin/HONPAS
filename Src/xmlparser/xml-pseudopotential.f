      module pseudopotential

      use sys
      use precision
      use ionew
      use m_pseudo_types
      use m_pseudo
      use flib_sax

      integer      :: iostat
      type(xml_t)  :: fxml
      
      private

      public :: pseudopotential_t, pseudo_read, pseudo_header_print

      integer, parameter        :: nrmax = 1500


      type pseudopotential_t
        character(len=2)        :: name
        integer                 :: nr
        integer                 :: nrval
        real(dp)                :: zval
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

!       PS information can be in a :
!               .vps file (unformatted), 
!               .psf file (formatted), 
!       or in a .xml file (xml-type formatted)

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
             fname = trim(label) // '.xml'
             inquire(file=fname, exist=found)
             if (found) then
              call pseudo_read_xml(fname,p)
             else
               write(6,'(/,2a,a20,/)') 'read_pseudo: ERROR: ',
     .              'Pseudopotential file not found: ', fname
               call die
             endif
           endif
        endif
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

        integer io_ps, i, j
        character(len=70) dummy
        real(dp) :: r2

        call io_assign(io_ps)
        open(io_ps,file=fname,form='formatted',status='unknown')
        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in formatted form from ', trim(fname)

 8000   format(1x,i2)
 8005   format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010   format(1x,6a10,/,1x,a70)
 8015   format(1x,2i3,i5,3g20.12)
 8030   format(4(g20.12))
 8040   format(1x,a)

        read(io_ps,8005) p%name, p%icorr, p%irel, p%nicore
        read(io_ps,8010) (p%method(i),i=1,6), p%text
        read(io_ps,8015) p%npotd, p%npotu, p%nr, p%b, p%a, p%zval

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

! For debugging -------------------------------
!        write(6,*)p%name
!        write(6,*)p%nr
!        write(6,*)p%nrval
!        write(6,*)p%zval
!        write(6,*)p%relativistic
!        write(6,*)p%correlation
!        write(6,*)p%icorr
!        write(6,*)p%irel
!        write(6,*)p%nicore
!        write(6,*)p%a
!        write(6,*)p%b
!        write(6,*)p%method(1)
!        write(6,*)p%method(2)
!        write(6,*)p%method(3)
!        write(6,*)p%method(4)
!        write(6,*)p%method(5)
!        write(6,*)p%method(6)
!        write(6,*)p%npotu
!        write(6,*)p%npotd
!        write(6,'(4f20.12)')p%r(1:4)
!        write(6,'(4f20.12)')p%r(100:104)
!        write(6,'(4f20.12)')p%chcore(1:4)
!        write(6,'(4f20.12)')p%chval(1:4)
!        write(6,*)p%ldown(:)
!        write(6,*)p%lup(:)
! ---                         
        call io_close(io_ps)
        end subroutine pseudo_read_formatted
!----
        subroutine pseudo_read_xml(fname,p)
        implicit none 
        character(len=*), intent(in) :: fname
        type(pseudopotential_t)                    :: p
        type(pseudo_t), pointer                    :: psxml

        write(6,'(3a)') 'Reading pseudopotential information ',
     $       'in XML form from ', trim(fname)

        call open_xmlfile(fname,fxml,iostat)
        if (iostat /=0) stop "Cannot open file"
                                                                                
        call xml_parse(fxml, 
     .          begin_element,end_element,pcdata_chunk,verbose=.false.)

        psxml => pseudo

        call xml2psf( psxml, p )

        end subroutine pseudo_read_xml

!----
        subroutine xml2psf( psxml, p )
! Translate the more complete xml format and data structure 
! into the old Siesta type for the pseudopotential
        implicit none 
        type(pseudo_t), intent(in)                 :: psxml
        type(pseudopotential_t), intent(out)       :: p

        integer          :: position, i, il, ir, lmax, lshellint
        character(len=1) :: ispp, lshell
        logical          :: polarized
        real(dp)         :: zeld(0:4), zelu(0:4)
        real(dp)         :: r2


        p%name         = psxml%header%symbol
        p%nr           = psxml%pot(1)%V%grid%npts
        p%nrval        = p%nr + 1
        p%zval         = psxml%header%zval
! relativistic and correlation are not needed in Siesta, so they are
! not translated
        select case(psxml%header%xcfunctionalparametrization)
          case('Ceperley-Alder')
             p%icorr = 'ca'
          case('Wigner')
             p%icorr = 'wi'
          case('Hedin-Lundqvist')
             p%icorr = 'hl'
          case('Gunnarson-Lundqvist')
             p%icorr = 'gl'
          case('von Barth-Hedin')
             p%icorr = 'bh'
          case('Perdew-Burke-Ernzerhof')
             p%icorr = 'pb'
          case('Becke-Lee-Yang-Parr')
             p%icorr = 'bl'
        end select

        select case(psxml%header%relativistic)
          case(.true.)
            p%irel    = 'rel'
            ispp      = 'r'
            polarized = .false.
          case(.false.)
            select case(psxml%header%polarized)
              case(.true.)
                p%irel    = 'isp'
                ispp      = 's'
                polarized = .true.
              case(.false.)
                p%irel    = 'nrl'
                ispp      = ' '
                polarized = .false.
            end select
        end select

        select case(psxml%header%core_corrections)
          case("yes")
            p%nicore = 'pcec'
          case("no")
            p%nicore = 'nc'
        end select

        p%a            = psxml%pot(1)%V%grid%step
        p%b            = psxml%pot(1)%V%grid%scale

        p%method(1)    = psxml%header%creator
        p%method(2)    = psxml%header%date
        read(psxml%header%flavor,'(4a10)') (p%method(i),i=3,6) 

        p%npotu        = psxml%npots_up
        p%npotd        = psxml%npots_down

! Allocate the radial variables and semilocal potentials
        allocate(p%r(1:p%nrval))
        allocate(p%chcore(1:p%nrval))
        allocate(p%chval(1:p%nrval))

        if (p%npotd.gt.0) then
           allocate(p%vdown(1:p%npotd,1:p%nrval))
           allocate(p%ldown(1:p%npotd))
        endif

        if (p%npotu.gt.0) then
           allocate(p%vup(1:p%npotu,1:p%nrval))
           allocate(p%lup(1:p%npotu))
        endif
! ---

! Calculate the points of the logarithmic radial grid 
        do 30 ir = 1, p%nrval
          p%r(ir) = p%b * (exp(p%a*(ir-1))-1)
 30     enddo
! ---

! Translate the valence charge density and the pseudo-core charge density,
! and define the value at the first point of the logarithmic grid
        p%chcore(2:p%nrval) = psxml%core_charge%data(1:p%nr)
        p%chval(2:p%nrval)  = psxml%valence_charge%data(1:p%nr)
        r2=p%r(2)/(p%r(3)-p%r(2))
        p%chcore(1) = p%chcore(2) - r2*(p%chcore(3)-p%chcore(2))
        p%chval(1) = p%chval(2) - r2*(p%chval(3)-p%chval(2))
! ---

        zeld(:) = 0.0d0
        zelu(:) = 0.0d0

        do il = 1, p%npotd
          lshell = psxml%pot(il)%l
          select case(lshell)
            case('s') 
              lshellint = 0
            case('p') 
              lshellint = 1
            case('d') 
              lshellint = 2
            case('f') 
              lshellint = 3
          end select
          p%ldown(il) = lshellint
          p%vdown(il,2:p%nrval) = psxml%pot(il)%V%data(1:p%nr)
          zeld(lshellint) = psxml%pot(il)%occupation
        enddo

        do il = 1, p%npotu
          lshell = psxml%pot(p%npotd+il)%l
          select case(lshell)
            case('s') 
              lshellint = 0
            case('p') 
              lshellint = 1
            case('d') 
              lshellint = 2
            case('f') 
              lshellint = 3
          end select
          p%lup(il) = lshellint
          p%vup(il,2:p%nrval) = psxml%pot(p%npotd+il)%V%data(1:p%nr)
          zelu(lshellint) = psxml%pot(p%npotd+il)%occupation
        enddo

        p%text = ' '
        position = 1
        lmax = max(p%npotd, p%npotu)
        do 240 il = 1, lmax
           if ( .not. polarized) then
              write(p%text(position:),9070) 
     .                                      psxml%pot(il)%n,
     .                                      psxml%pot(il)%l,
     .                                      zeld(il-1)+zelu(il-1), 
     .                                      ispp, 
     .                                      psxml%pot(il)%cutoff
 9070         format(i1,a1,f5.2,a1,' r=',f5.2,'/')
              position = position + 17
           else
              write(p%text(position:),9090)
     .                                      psxml%pot(il)%n,
     .                                      psxml%pot(il)%l,
     .                                      zeld(il-1), 
     .                                      zelu(il-1), 
     .                                      ispp, 
     .                                      psxml%pot(il)%cutoff
 9090         format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
              position = position + 17
           end if
  240   enddo

       
! For debugging -------------------------------------------------------
!        write(6,*)p%name
!        write(6,*)p%nr
!        write(6,*)p%nrval
!        write(6,*)p%zval
!        write(6,*)p%icorr
!        write(6,*)p%irel
!        write(6,*)p%nicore
!        write(6,*)p%a
!        write(6,*)p%b
!        write(6,*)p%method(1)
!        write(6,*)p%method(2)
!        write(6,*)p%method(3)
!        write(6,*)p%method(4)
!        write(6,*)p%method(5)
!        write(6,*)p%method(6)
!        write(6,*)p%text
!        write(6,*)p%npotu
!        write(6,*)p%npotd
!        write(6,'(4f20.12)')p%r(1:4)
!        write(6,'(4f20.12)')p%r(101:104)
!        write(6,'(4f20.12)')p%vdown(1,1:4)
!        write(6,'(4f20.12)')p%vdown(2,1:4)
!        write(6,'(4f20.12)')p%vdown(3,1:4)
!        write(6,'(4f20.12)')p%vdown(4,1:4)
!        write(6,'(4f20.12)')p%vup(1,1:4)
!        write(6,'(4f20.12)')p%vup(2,1:4)
!        write(6,'(4f20.12)')p%vup(3,1:4)
!        write(6,'(4f20.12)')p%chcore(1:4)
!        write(6,'(4f20.12)')p%chval(1:4)
!        write(6,*)p%ldown
!        write(6,*)p%lup
!         do ir = 1, p%nrval
!           write(6,*)p%r(ir), pseudo%pswf(1)%V%data(ir)
!         enddo
! -----
        end subroutine xml2psf
!------
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
        end module pseudopotential



