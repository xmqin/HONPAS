! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

program readwfx

  use m_getopts
  use f2kcli

! This program READWF reads a the coefficients of wavefunctions
! on the expansion of the atomic orbitals basis set, as written
! by SIESTA (unformatted) and writes them in ascii, user-friendly
! form.

! Written by P. Ordejon, June 2003
! Updated to WFSX format and given more features by A. Garcia.

!  threshold to plot the coefficients of the wavefunctions. If the
!  norm of the weight of a wavefunction on a given orbital is smaller
!  than this threshold then it is not printed. This is useful for
!  large systems, with a very large number of basis orbitals.

implicit none
        
integer, parameter :: dp = selected_real_kind(10,100)
integer, parameter :: sp = selected_real_kind(5,10)

  character(len=200) :: opt_arg
  character(len=10)  :: opt_name 
  integer :: nargs, iostat, n_opts, nlabels

  integer :: io, wfs_u, nk, nspin_blocks, ik, ispin, is, idummy, &
       number_of_wfns, iw, indwf, j, nuotot, jj, nspin_flag, is0

  character(len=256) :: wfsx_file, outfile = ""

  integer, allocatable, dimension(:) :: iaorb,iphorb,cnfigfio
  character(len=20), allocatable, dimension(:) :: symfio,labelfis

  real(sp), allocatable, dimension(:,:) :: psi
  logical  :: gamma, non_coll
  real(dp) :: k(3), eigval, norm, wk
        
  real(dp) :: threshold = 0.0_dp

  real(dp) :: emin        = -huge(1.0_dp)
  real(dp) :: emax        =  huge(1.0_dp)
  logical  :: emin_given  = .false.
  logical  :: emax_given  = .false.
  integer  :: min_band = -huge(1)
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: energies_only = .false.

  real(dp) :: min_eigval, max_eigval
  real(dp) :: min_eigval_in_band_set, max_eigval_in_band_set
  integer  :: nwfmin, nwfmax, wfs_spin_flag

!     Process options
!
  n_opts = 0
  do
     call getopts('hlt:e:m:E:M:b:B:o:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('o')
        read(opt_arg,*) outfile
     case ('m', 'e')
        emin_given = .true.
        read(opt_arg,*) emin
     case ('M', 'E')
        emax_given = .true.
        read(opt_arg,*) emax
     case ('b')
        read(opt_arg,*) min_band
     case ('B')
        read(opt_arg,*) max_band
     case ('t')
        read(opt_arg,*) threshold
     case ('l')
        energies_only = .true.
     case ('h')
        call manual()
        STOP
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Use -h option for manual"
        write(0,*) ""
        call manual()
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Use -h option for manual"
     write(0,*) ""
     call manual()
     STOP
  endif

  call get_command_argument(n_opts,value=wfsx_file,status=iostat)
  if ( iostat /= 0 ) then
     stop "Cannot get WFSX file"
  end if

  threshold = threshold**2  ! We use the square later on

  wfs_u = 10
  io = 6

  open(wfs_u, file=wfsx_file, form='unformatted', status='old' )

  rewind(wfs_u)
  read(wfs_u) nk, gamma

  read(wfs_u) wfs_spin_flag   !  1, 2, or 4 
  non_coll = (wfs_spin_flag >= 4)
  read(wfs_u) nuotot
  read(wfs_u)        !! Symbols, etc


  if (non_coll) then
     nspin_blocks = 1
  else
     nspin_blocks = wfs_spin_flag
  endif

  !-------------------------------------
  
  nwfmax = -huge(1)
  nwfmin = huge(1)
  min_eigval = huge(1.0_dp)
  max_eigval = -huge(1.0_dp)
  min_eigval_in_band_set = huge(1.0_dp)
  max_eigval_in_band_set = -huge(1.0_dp)

  do ik=1,nk
     do is=1,nspin_blocks

        read(wfs_u) idummy, k(1:3), wk
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_u) is0
        read(wfs_u) number_of_wfns
        nwfmax = max(nwfmax,number_of_wfns)
        nwfmin = min(nwfmin,number_of_wfns)

        do iw=1,number_of_wfns
           read(wfs_u) indwf
           read(wfs_u) eigval
           min_eigval = min(min_eigval,eigval)
           max_eigval = max(max_eigval,eigval)
           ! 
           !
           if ((iw>=min_band).and.(iw<=max_band)) then
              min_eigval_in_band_set = min(min_eigval_in_band_set,eigval)
              max_eigval_in_band_set = max(max_eigval_in_band_set,eigval)
           endif
           read(wfs_u)
        enddo
     enddo
  enddo

  print "(a,2i5)", "Minimum/Maximum number of wfs per k-point: ", nwfmin, nwfmax
  print "(a,2f12.4)", "Min_eigval, max_eigval on WFS file: ",  &
       min_eigval, max_eigval

  print "(a,2f12.4)", "Min_eigval, max_eigval in band set : ",  &
          min_eigval_in_band_set, max_eigval_in_band_set

  if (energies_only) STOP   ! energies only

  if (outfile /= "") then
     open(io, file=outfile, form='formatted', status='unknown' )
  endif

  rewind (wfs_u)

  read(wfs_u) nk, gamma
  read(wfs_u) nspin_flag
  read(wfs_u) nuotot

  if (nspin_flag == 4) then
     non_coll = .true.
     nspin_blocks = 1
  else
     non_coll = .false.
     nspin_blocks = nspin_flag
  endif
        
  if (gamma) then
     if (non_coll) then
        allocate(psi(4,nuotot))
     else
        allocate(psi(1,nuotot))
     endif
  else
     if (non_coll) then
        allocate(psi(4,nuotot))
     else
        allocate(psi(2,nuotot))
     endif
  endif

  allocate(iaorb(nuotot), labelfis(nuotot), &
       iphorb(nuotot), cnfigfio(nuotot),    &
       symfio(nuotot))

  read(wfs_u) (iaorb(j),labelfis(j),  &
       iphorb(j), cnfigfio(j), symfio(j), j=1,nuotot)

  write(io,*)
  write(io,'(a22,2x,i6)') 'Nr of k-points = ',nk
  if (non_coll) then
     write(io,'(a,2x,i6)') 'Nr of Spins blocks ' // &
                                  '(non-collinear)  = ', nspin_blocks
  else
     write(io,'(a,2x,i6)') 'Nr of Spins blocks = ',nspin_blocks
  endif
  write(io,'(a,2x,i6)') 'Nr of basis orbs = ',nuotot
  write(io,*)

 kpoints:  do ik = 1,nk
    spin_blocks: do ispin = 1,nspin_blocks

        read(wfs_u) idummy, k(1),k(2),k(3)
        if (idummy .ne. ik) stop 'error in index of k-point'
        ! if ik not in k-range CYCLE kpoints
        read(wfs_u) idummy
        if (.not. non_coll) then
           if (idummy .ne. ispin) stop 'error in index of spin'
           ! if ispin not in requested spin-channel CYCLE spin_blocks
        endif
        read(wfs_u) number_of_wfns

        write(io,*)
        write(io,fmt='(a,2x,i6,2x,3f10.6)',advance="no") 'k-point = ',ik, k(1),k(2),k(3)
        if (.not. non_coll) then
           write(io,fmt='(2x,a,2x,i1)',advance="no") 'Spin component = ',ispin
        endif
        write(io,fmt='(2x,a,1x,i6)') 'Num. wavefunctions = ',number_of_wfns

! Loop over wavefunctions 

        do iw = 1,number_of_wfns

           read(wfs_u) indwf
           read(wfs_u) eigval
           if ( (eigval < emin) .or. (eigval > emax) .or. &
                (iw < min_band) .or. (iw > max_band)) then
              read(wfs_u)  ! skip over wfn data
              CYCLE
           endif

           write(io,'(a22,2x,i6,2x,a,f10.6)') &
                'Wavefunction = ', indwf, 'Eigval (eV) = ', eigval

           if (threshold >= 100*100) then
              read(wfs_u)   ! Skip wf data
              CYCLE         ! Early exit to print only energies
           endif

           write(io,"(60('-'))")
           if (non_coll) then
              write(io,'(a72)') &
                           ' Atom  Species Orb-global Orb-in-atom '//  &
                           '  .Orb-type    [Re(psi)  Im(psi)]Up ' //   &
                           ' [Re(psi)  Im(psi)]Down '
           else
              write(io,'(a72)') &
                    ' Atom  Species Orb-global  Orb-in-atom ' //  &
                    ' .Orb-type      Re(psi)   Im(psi)'
           endif
            
           read(wfs_u) (psi(1:,j), j=1,nuotot)
           do jj = 1,nuotot
              norm = dot_product(psi(1:,jj),psi(1:,jj))  ! valid for all
              if (norm .lt. threshold) CYCLE
              write(io,  &
                   '(i6,5x,a10,1x,i10,2x,i3,2x,i1,a20,2(2(f10.6),2x))')  &
                   iaorb(jj),labelfis(jj),jj, iphorb(jj), cnfigfio(jj),  &
                   symfio(jj), psi(1:,jj)

           enddo

           write(io,"(60('-'))")

          enddo
        enddo spin_blocks
      enddo kpoints

      close (wfs_u)

    CONTAINS
            subroutine manual

      write(6,*) "Usage: readwfx [ options ] WFSX_FILE"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "   -o File   :  set output file (default: use stdout)       "
      write(6,*) "    "
      write(6,*) "   Selection of states to be printed: by eigenvalue range or band index:  "
      write(6,*) "    "
      write(6,*) "   -m Min_e  :  set lower bound of eigenvalue range                    "
      write(6,*) "   -M Max_e  :  set upper bound of eigenvalue range                    "
      write(6,*) "   -b Min_band  :  set minimum band index to be used               "
      write(6,*) "   -B Max_band  :  set maximum band index to be used               "
      write(6,*) "    "
      write(6,*) "   -t threshold :  set threshold for orbital contributions to be printed "
      write(6,*) "                   (use large value to get only energy information)      "
      end subroutine manual

    end program readwfx
