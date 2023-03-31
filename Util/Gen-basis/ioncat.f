! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      program ioncat
!
!     Looks inside .ion files (which hold the PAO, KB, and Vna tables)
!     Generates data files for plotting, optionally with zoom feature
!     near rc.

      use precision
      use atm_types, only: species_info, nspecies, species
      use atmfuncs
      use m_getopts

      use basis_io, only: read_ion_ascii

      implicit none

C Internal variables ...................................................
      integer is, i, j
      integer, parameter :: n=1000
      real(dp) :: delta, rc, rmin, rmax, range

      real(dp) :: x, val, grad, rlog
      real(dp) :: gradv(3)
      type(species_info), pointer   :: spp

      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
      integer :: no, nkb

      logical :: show_all    = .false.,
     $           show_shells = .false.,
     $           show_kbshells = .false.,
     $           zoom        = .false.,
     $           near_origin = .false.,
     $           process_orb = .false.,
     $           process_kbp = .false.,
     $           process_vna = .false.,
     $           process_chlocal = .false.,
     $           process_core = .false.
      character(len=1000) nameat
!
!     Process options
!
      n_opts = 0
      do
         call getopts('sijo:k:vclZO:h',opt_name,opt_arg,n_opts,iostat)
         if (iostat /= 0) exit
         select case(opt_name)
         case ('s', '+s')
            show_all = .true.
            write(0,*) "Will show everything ", trim(opt_arg)
         case ('i', '+i')
            show_shells = .true.
         case ('j', '+j')
            show_kbshells = .true.
         case ('o', '+o')
            process_orb = .true.
            write(0,*) "Will process orb ", trim(opt_arg)
            read(opt_arg,*) iorb
         case ('k', '+k')
            process_kbp = .true.
            read(opt_arg,*) ikb
            write(0,*) "Will process KB proj ", ikb
         case ('c', '+c')
            process_core = .true.
            write(0,*) "Will process pseudo core charge "
         case ('l', '+l')
            process_chlocal = .true.
            write(0,*) "Will process Vlocal charge density "
         case ('v', '+v')
            process_vna = .true.
            write(0,*) "Will process Vna "
         case ('Z', '+Z')
            zoom = .true.
            write(0,*) "Will use zoom around rc"
         case ('O', '+O')
            near_origin = .true.
            read(opt_arg,*) rlog
            write(0,*) "Will zoom near zero, up to ", 10.0_dp**(rlog)
         case ('h', '+h')
            call print_help()
            STOP
         case ('?',':')
            write(0,*) "Invalid option: ", opt_arg(1:1)
            write(0,*) "Use -h option for help"
            STOP
         end select
      enddo

      nargs = command_argument_count()
      nlabels = nargs - n_opts + 1
      if (nlabels /= 1)  then
         call print_help()
         STOP
      endif

      call get_command_argument(n_opts,value=nameat,status=iostat)
      if (iostat /= 0) then
          STOP "Cannot get Species_Label"
      endif

      nspecies = 1
      allocate(species(1))
      spp => species(1)
      spp%label = trim(nameat)
      spp%read_from_file = .true.
      call read_ion_ascii(spp)

      no = nofis(1)
      nkb = nkbfis(1)

      if (show_all) then
         write(0,*) "Atomic number, norbs, nkbs: ", izofis(1), no, nkb
         write(0,*) "Orbitals (#, l, z, m, rc):"
         do i= 1, no
            write(6,*) i,
     $               lofio(1,i), zetafio(1,i), mofio(1,i), rcut(1,i)
         enddo
         write(6,*) "KB projs (#, l, m, rc):"
         do i= 1, nkb
            write(6,*) i, lofio(1,-i), mofio(1,-i), rcut(1,-i)
         enddo
         write(6,*) "Vna rcut: ", rcut(1,0)
      endif
!
      if (show_shells) then
         i = 0
         do 
            if (i == no) exit
            i = i + 1
            write(6,fmt="(i3)",advance="no") i
            i = i +  2* lofio(1,i)           ! Skip other m copies
         enddo
         write(6,*)
      endif
!
      if (show_kbshells) then
         i = 0
         do 
            if (i == nkb) exit
            i = i + 1
            write(6,fmt="(i3)",advance="no") i
            i = i +  2* lofio(1,-i)           ! Skip other m copies
         enddo
         write(6,*)
      endif
!
      if (process_orb) then
         if (iorb > no) STOP "no such orbital"
         write(6,*) "# Orbital (#, l, z, m, rc):",
     $    lofio(1,iorb), zetafio(1,iorb), mofio(1,iorb), rcut(1,iorb)
         rc = rcut(1,iorb)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         if (near_origin) then
            rmin = 0.0_dp 
            rmax = 10.0_dp**(rlog)
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,iorb,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      if (process_kbp) then
         if (ikb > nkb) STOP "no such KB projector"
         write(6,*) "# KB proj (#, l, m, rc):",
     $        ikb, lofio(1,-ikb), mofio(1,-ikb), rcut(1,-ikb)
         rc = rcut(1,-ikb)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         if (near_origin) then
            rmin = 0.0_dp 
            rmax = 10.0_dp**(rlog)
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,-ikb,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      if (process_vna) then
         write(6,*) "# Vna, rcut: ", rcut(1,0)
         rc = rcut(1,0)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         if (near_origin) then
            rmin = 0.0_dp 
            rmax = 10.0_dp**(rlog)
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call rphiatm(1,0,x,val,grad)
            write(*,"(3f14.8)") x, val, grad
         enddo
      endif

      if (process_core) then
         write(6,*) "# Core charge, rcut: ", rcore(1)
         rc = rcore(1)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         if (near_origin) then
            rmin = 0.0_dp 
            rmax = 10.0_dp**(rlog)
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call chcore_sub(1,(/x,0.0_dp,0.0_dp/),val,gradv)
            write(*,"(3f14.8)") x, val, gradv(1)
         enddo
      endif

      if (process_chlocal) then
         write(6,*) "# Local ps charge, rcut: ", rchlocal(1)
         rc = rchlocal(1)
         rmin = 0.0_dp
         rmax = 1.05_dp * rc
         if (zoom) then
            rmin = 0.95_dp * rc
            rmax = 1.05_dp * rc
         endif
         if (near_origin) then
            rmin = 0.0_dp 
            rmax = 10.0_dp**(rlog)
         endif
         range = rmax - rmin
         delta = range/n
         do i=0,n
            x = rmin + delta*i
            call psch(1,(/x,0.0_dp,0.0_dp/),val,gradv)
            write(*,"(3f14.8)") x, val, gradv(1)
         enddo
      endif

      end program ioncat

      subroutine print_help()
         print *, "Usage: ioncat [options] Species_Label"
         print *, "Options:"
         print *, " "
         print *, " -s          : Show header information"
         print *, " -i          : Print indexes of unique orbitals"
         print *, " -j          : Print indexes of unique KB projectors"
         print *, " -o ORBINDEX : Generate table for orbital ORBINDEX"
         print *, " -k KBINDEX  : Generate table for KB proj KBINDEX"
         print *, " -v          : Generate table for Vna potential"
         print *, " -c          : Generate table for pseudocore charge"
         print *, " -l          : Generate table for Vlocal charge"
         print *, " -Z          : Zoom in near rc for table generation"
         print *, " -O rlog     : Zoom near zero (up to 10^rlog)"
         print *, " -h          : Print this help message"
         print *, " "
      end subroutine print_help



