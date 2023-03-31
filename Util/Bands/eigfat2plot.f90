! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

program eigfat2plot

  use m_getopts
  use f2kcli

  implicit none

  integer           nk, nspin, nband, ik, is, ib
  real              ::  kl, delta(3)
  real, allocatable ::  k(:), pk(:,:)
  real, allocatable ::  fat(:,:,:), eig(:,:,:)


  integer  :: min_band = 1
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: band_interval_set = .false.


  integer  :: fat_u = 2
  integer  :: n_opts, iostat, nargs, nlabels
  character(len=132) :: opt_name, opt_arg
  character(len=132) :: fatfile, string_id
  integer  :: ika

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('hb:B:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('h')
        !        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: eigfat2plot [ -b -B ] <EIGFAT file>"
        write(0,*) "       -b and -B relative to eigfat file! "
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: eigfat2plot [ -b -B ] <EIGFAT file>"
     write(0,*) "       -b and -B relative to eigfat file! "
     STOP
  endif

  call get_command_argument(n_opts,value=fatfile,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get EIGFAT file"
  endif

  open(fat_u,file=trim(fatfile), status='old',action="read",position="rewind")
  read(fat_u,"(a132)") string_id
  read(fat_u,*) nband, nspin, nk

  print *, "# " // trim(string_id)
  print *, "# nband, nspin, nk:", nband, nspin, nk

  if (.not. max_band_set) max_band = nband

  if (min_band_set .and. (min_band < 1)) then
     print "(a)", " ** Min_band implicitly reset to 1..."
     min_band = 1
  endif
  if (min_band_set .and. (min_band > nband)) then
     print "(a,2i5)", " ** Min_band is too large for some k-points: (min_band, nband):", min_band, nband
     STOP
  endif
  if (max_band_set .and. (max_band > nband)) then
     print "(a,2i5)", " ** Max_band is too large for some k-points: (max_band, nband):", max_band, nband
     print "(a)", " ** Max_band will be effectively reset to its maximum allowed value"
     max_band = nband
  endif
  if (max_band_set .and. (max_band < max(1,min_band))) then
     print "(a,2i5)", " ** Max_band is less than the effective min_band: (max_band, eff min_band):", max_band, max(1,min_band)
     STOP
  endif


  print *, "# min_band, max_band: ", min_band, max_band

  allocate(k(nk),pk(3,0:nk))
  allocate(fat(nband,nspin,nk),eig(nband,nspin,nk))

  kl = 0.0
  do ik = 1, nk
     read(fat_u,*) ika, pk(1:3,ika)
     if (ik == 1) then
        pk(1:3,0) = pk(1:3,1)
     endif
     delta(1:3) = pk(1:3,ika) - pk(1:3,ika-1)
     kl = kl + sqrt(dot_product(delta,delta))
     k(ik) = kl
     do is = 1, nspin
        read(fat_u,*) (eig(ib,is,ik),fat(ib,is,ik), ib = 1, nband)
     enddo
  enddo
  close(fat_u)


  do is = 1, nspin
     do ib = min_band, max_band
        do ik = 1 , nk
           ! write spin variable to differentiate
           write(6,"(3f14.6,i3)") k(ik), eig(ib,is,ik), fat(ib,is,ik), is
        enddo
        write(6,"(/)")
     enddo
  enddo

  deallocate(eig,fat,pk,k)

end program eigfat2plot
