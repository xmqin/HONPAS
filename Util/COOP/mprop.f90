! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!==================================================================
program mprop

  use main_vars
  use subs, only: manual
  use orbital_set, only: get_orbital_set
  use io_hs, only: read_hs_file
  use read_curves, only: read_curve_information, mask_to_arrays

  implicit none

  logical :: gamma_wfsx, got_qcos, non_coll
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  integer :: imin, imax, wfs_spin_flag, nspin_blocks
  
  complex(DP), dimension(2)    :: spinor_1, spinor_2, H_c2
  complex(DP), dimension(2,2)  :: H
  complex(DP) :: c_c1_c2, c_c1_H_c2
  real(dp) :: factor_S, factor_H, dos_value

  ! We use a smearing function of the form f(x) = exp(-(x/smear)**2) / (smear*sqrt(pi))
  ! A weight tolerance of 1.0e-4 corresponds to going about 3*smear on either
  ! side of the eigenvalue. We should set nsigma appropriately

  real(dp), parameter  :: tol_weight = 1.0e-4_dp
  integer, parameter   :: nsigma = 3


  real(dp), parameter  :: tol_overlap = 1.0e-10_dp

  logical, allocatable   :: mask2(:)
  integer, allocatable   :: num_red(:), ptr(:), list_io2(:), list_ind(:)

  logical  :: enough_electrons 

  integer  :: nwfmx, nwfmin
  integer  :: min_band = -huge(1)
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: band_interval_set = .false.
  real(dp) :: min_eigval
  real(dp) :: max_eigval
  real(dp) :: min_eigval_in_file
  real(dp) :: min_eigval_in_band_set
  real(dp) :: max_eigval_in_band_set
  real(dp) :: minimum_spec_eigval = -huge(1.0_dp)
  real(dp) :: maximum_spec_eigval = huge(1.0_dp)
  real(dp) :: ewindow_low = -huge(1.0_dp)
  real(dp) :: ewindow_high = huge(1.0_dp)

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhls:n:m:M:R:b:B:w:W:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('+d')
        debug = .false.
     case ('l')
        energies_only = .true.
     case ('s')
        read(opt_arg,*) smear
     case ('n')
        read(opt_arg,*) npts_energy
     case ('m')
        read(opt_arg,*) minimum_spec_eigval
     case ('M')
        read(opt_arg,*) maximum_spec_eigval
     case ('R')
        ref_line_given = .true.
        ref_line = opt_arg
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('w')
        read(opt_arg,*) ewindow_low
     case ('W')
        read(opt_arg,*) ewindow_high
     case ('h')
        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: mprop [ -d ] [ -h ] MPROP_FILE_ROOT"
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: mprop [ -d ] [ -h ] MPROP_FILE_ROOT"
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=mflnm,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get .mprop file root"
  endif

  band_interval_set = (min_band_set .or. max_band_set)

  !==================================================

  ierr=0

  ! Read type of job

  open(mpr_u,file=trim(mflnm) // ".mpr", status='old')
  read(mpr_u,*) sflnm
  read(mpr_u,*) what

  dos=(trim(what).eq.'DOS')
  coop=(trim(what).eq.'COOP')
  
  !==================================================
  ! Read WFSX file

  write(6,"(a)") "Reading wave-function file: " // trim(sflnm) // ".WFSX..."
  write(6,"(a)") "Energy units are eV"

  open(wfs_u,file=trim(sflnm)//'.WFSX',status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) wfs_spin_flag   !  1, 2, or 4 
  non_coll = (wfs_spin_flag >= 4)
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "WFSX read: nkp, nspin_flag, nnao: ", nkp, wfs_spin_flag, nao

  if (non_coll) then
     nspin_blocks = 1
  else
     nspin_blocks = wfs_spin_flag
  endif
  print *, "NSPIN BLOCKS: ", nspin_blocks

  allocate (ados(npts_energy,nspin_blocks), ww(npts_energy))

  nwfmx = -huge(1)
  nwfmin = huge(1)
  min_eigval = huge(1.0_dp)
  max_eigval = -huge(1.0_dp)
  min_eigval_in_band_set = huge(1.0_dp)
  max_eigval_in_band_set = -huge(1.0_dp)


  do ik=1,nkp
     do is=1,nspin_blocks

        read(wfs_u) idummy, pk(1:3,ik), wk(ik)
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_u) is0
        read(wfs_u) number_of_wfns
        nwfmx = max(nwfmx,number_of_wfns)
        nwfmin = min(nwfmin,number_of_wfns)

        do iw=1,number_of_wfns
           read(wfs_u) iw0
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

  print "(a,2i5)", "Minimum/Maximum number of wfs per k-point: ", nwfmin, nwfmx
  print "(a,2f12.4)", "Min_eigval, max_eigval on WFS file: ",  &
       min_eigval, max_eigval
  print "(a,2f12.4)", "Min_eigval, max_eigval in band set : ",  &
          min_eigval_in_band_set, max_eigval_in_band_set

  min_eigval_in_file = min_eigval    ! Saved for E_Fermi logic

  if (band_interval_set) then

     if (min_band_set .and. (min_band < 1)) then
        print "(a)", " ** Min_band implicitly reset to 1..."
     endif
     if (min_band_set .and. (min_band > nwfmin)) then
        print "(a,2i5)", " ** Min_band is too large for some k-points: (min_band, nwfmin):", min_band, nwfmin
        STOP
     endif
     if (max_band_set .and. (max_band > nwfmin)) then
        print "(a,2i5)", " ** Max_band is too large for some k-points: (max_band, nwfmin):", max_band, nwfmin
        print "(a)", " ** Max_band will be effectively reset to its maximum allowed value"
     endif
     if (max_band_set .and. (max_band < max(1,min_band))) then
        print "(a,2i5)", " ** Max_band is less than the effective min_band: (max_band, eff min_band):", max_band, max(1,min_band)
        STOP
     endif

     min_eigval = min_eigval_in_band_set
     max_eigval = max_eigval_in_band_set
  endif
  print "(a,3i4)", "Implicit band set used: (min, max_min, max_max):",  &
                   max(1,min_band), min(nwfmin,max_band), min(nwfmx,max_band)

  ! min_eigval, max_eigval: Determine which eigenstates are used to compute the curves
  print "(a,2f12.4)", "Minimum and maximum eigenvalues (based on file data and band selection): ",  min_eigval, max_eigval

  if (minimum_spec_eigval > min_eigval) then
     min_eigval = minimum_spec_eigval
     print "(a,f12.4)", "* Minimum eigenvalue changed as per user range request: ",  min_eigval
  endif
  if (maximum_spec_eigval < max_eigval) then
     max_eigval = maximum_spec_eigval
     print "(a,f12.4)", "* Maximum eigenvalue changed as per user range request: ",  max_eigval
  endif

  ! Sanity checks
  print "(a,2f12.4)", "Minimum and maximum eigenvalues to be processed: ",  min_eigval, max_eigval
  if (min_eigval > max_eigval) STOP "Meaningless range. Check -b/-B and -m/-M options"

  ! Here low_e and high_e represent a window for the plot.
  ! Avoid cut tails by extending the eigenvalue range on both sides

  low_e = min_eigval - nsigma*smear
  high_e = max_eigval + nsigma*smear

  print "(a,2f12.4)", "Plotting range adequate for eigenvalue range selected: ",  &
       low_e, high_e

  if (ewindow_low /= -huge(1.0_dp))  then
     low_e = ewindow_low
     print "(a,f12.4)", "* Lower bound of plotting range changed as per user request: ", low_e
  endif
  if (ewindow_high /= huge(1.0_dp))  then
     high_e = ewindow_high
     print "(a,f12.4)", "* Upper bound of plotting range changed as per user request: ", high_e
  endif

  ! Sanity checks
  if (ewindow_low >= ewindow_high) STOP "Meaningless plotting window. Check -w/-W options"
  if (ewindow_low >= max_eigval) STOP "Ewindow_low > max eigenvalue used. Check -b/-B, -m/-M, -w/-W options"
  if (ewindow_high <= min_eigval) STOP "Ewindow_high < min eigenvalue used. Check -b/-B, -m/-M, -w/-W options"

  print "(a,f7.3)", "Using smearing parameter: ", smear
  print "(a,i6,a)", "Using ", npts_energy, " points in energy range"


  e_step = (high_e-low_e)/(npts_energy-1)
  ados(:,1:nspin_blocks) = 0.0_dp

  ! skip four records

  rewind(wfs_u)

  read(wfs_u) 
  read(wfs_u) 
  read(wfs_u) 
  read(wfs_u) 

  do ik=1,nkp
     do is=1,nspin_blocks
        read(wfs_u)
        read(wfs_u)
        read(wfs_u)  number_of_wfns
        do iw=1,number_of_wfns
           read(wfs_u) 
           read(wfs_u) eigval
           if ( (iw>=min_band) .and. (iw<=max_band)) then
              do i = 1, npts_energy
                 energy = low_e + e_step*(i-1)
                 weight = delta(energy-eigval)
                 if (weight < tol_weight) CYCLE           ! Will not contribute
                 ados(i,is) = ados(i,is) + wk(ik) * weight
              enddo
           endif
           read(wfs_u)       ! Skip wfn info
        enddo
     enddo
  enddo

  ! Write "LARGE-SCALE DOS"  
  call write_curve_to_file(trim(sflnm)//".alldos",ados)

  ! Read HSX file
  ! Will pick up atoms, zval, and thus the nominal number of electrons,
  ! but the total charge is read as qtot.

  call read_hs_file(trim(sflnm)//".HSX")

  ztot = 0.0_dp
  do ia = 1, na_u
     ztot = ztot + zval(isa(ia))
  enddo
  if (abs(qtot-ztot) > 1.0e-8_dp) then
     print "(a,f14.4)", "Note: The system is charged: ", qtot-ztot
  endif

  !
  !       Compute integrated total DOS

  allocate(intdos(npts_energy), intebs(npts_energy))
  call io_assign(intdos_u)
  open(intdos_u,file=trim(sflnm)//".intdos",form="formatted", &
       status="unknown",action="write",position="rewind")
  intdos(1) = 0.0_dp
  intebs(1) = 0.0_dp
  write(intdos_u,*) low_e, intdos(1), intebs(1)
  do i = 2, npts_energy
     energy = low_e + e_step*(i-1)

     ! For the spinless case, we need a factor of two for proper normalization here
     dos_value = sum(ados(i,:))
     if (wfs_spin_flag == 1) dos_value = 2*dos_value

     intdos(i) = intdos(i-1) + dos_value * e_step 
     intebs(i) = intebs(i-1) + energy*dos_value * e_step 
     write(intdos_u,*) energy, intdos(i), intebs(i)
  enddo
  call io_close(intdos_u)

  if ( (max(min_band,1) > 1) .OR. (min_eigval > min_eigval_in_file)) then
     write(6,"(a,f10.5,a)") "Not meaningful to compute Fermi energy, as min_band>1 or restricted eigvals"
  else
     enough_electrons = .false.
     ! Look for Fermi Energy
     do i = 2, npts_energy
        if (intdos(i) > qtot) then
           enough_electrons = .true.
           ! Found fermi energy
           energy = low_e + e_step*(i-1)
           ! Correct overshoot, interpolating linearly
           efermi = energy - (intdos(i)-ztot)*e_step/(intdos(i)-intdos(i-1))
           exit
        endif
     enddo
     if (enough_electrons) then
        write(6,"(a,f10.5,a)") "Fermi energy: ", efermi, " (depends on smearing)"
     else
        write(6,"(a,f10.5,a)") "The band set does not contain enough electrons to compute E_Fermi"
     endif
  endif



  if (energies_only) STOP


  !====================

  ! * Orbital list

  allocate(za(no_u), zc(no_u), zn(no_u), zl(no_u), zx(no_u), zz(no_u))
  nao = 0
  do ia=1,na_u
     it = isa(ia)
     io = 0
     do 
        io = io + 1
        if (io > no(it)) exit
        lorb = lquant(it,io)
        do ko = 1, 2*lorb + 1
           nao = nao + 1
           za(nao)=ia
           zc(nao)=it
           zn(nao)=nquant(it,io)
           zl(nao)=lorb
           zx(nao)=ko
           zz(nao)=zeta(it,io)
        enddo
        io = io + 2*lorb
     enddo
  enddo
  if (nao /= no_u) STOP "nao /= no_u"

  ! ==================================

!
! Process orbital sets
!
  allocate(orb_mask(no_u,2,ncbmx))
  allocate (koc(ncbmx,2,no_u))
  call read_curve_information(dos,coop,  &
                                    mpr_u,no_u,ncbmx,ncb,tit,orb_mask,dtc)
  if (dos) then
     orb_mask(:,2,1:ncb) = .true.       ! All orbitals considered
  endif
  call mask_to_arrays(ncb,orb_mask(:,1,:),noc(:,1),koc(:,1,:))
  call mask_to_arrays(ncb,orb_mask(:,2,:),noc(:,2),koc(:,2,:))

!!
  write(6,"('Writing files: ',a,'.stt ...')") trim(mflnm)
  open(stt_u,file=trim(mflnm)//'.stt')
  write(stt_u,"(/'UNIT CELL ATOMS:')")
  write(stt_u,"(3x,i4,2x,i3,2x,a20)") (i, isa(i), label(isa(i)), i=1,na_u)
  write(stt_u,"(/'BASIS SET:')")
  write(stt_u,"(5x,a20,3(3x,a1))") 'spec', 'n', 'l', 'z'
  do it=1,nspecies
     write(stt_u,"(5x,a20)") trim(label(it))
     io = 0
     do 
        io = io + 1
        if (io > no(it)) exit
        write(stt_u,"(3(2x,i2))") nquant(it,io), lquant(it,io), zeta(it,io)
        io = io + 2*lquant(it,io)
     enddo
  enddo

  write(stt_u,"(/'SPIN: ',i2)") nsp
  write(stt_u,"(/'FERMI ENERGY: ',f18.6)") efermi

  write(stt_u,"(/'AO LIST:')")
  taux=repeat(' ',len(taux))
  do io=1,no_u
     taux(1:30)=repeat(' ',30)
     ik=1
     if (zl(io).eq.0) then
        taux(ik:ik)='s'
        ik=0
     elseif (zl(io).eq.1) then
        if (zx(io).eq.1) taux(ik:ik+1)='py'
        if (zx(io).eq.2) taux(ik:ik+1)='pz'
        if (zx(io).eq.3) taux(ik:ik+1)='px'
        ik=0
     elseif (zl(io).eq.2) then
        if (zx(io).eq.1) taux(ik:ik+2)='dxy'
        if (zx(io).eq.2) taux(ik:ik+2)='dyz'
        if (zx(io).eq.3) taux(ik:ik+2)='dz2'
        if (zx(io).eq.4) taux(ik:ik+2)='dxz'
        if (zx(io).eq.5) taux(ik:ik+5)='dx2-y2'
        ik=0
     elseif (zl(io).eq.3) then
        taux(ik:ik)='f'
     elseif (zl(io).eq.4) then
        taux(ik:ik)='g'
     elseif (zl(io).eq.5) then
        taux(ik:ik)='h'
     endif
     write(stt_u,"(3x,i5,2x,i3,2x,a20)",advance='no')  &
                          io, za(io), trim(label(zc(io)))
     if (ik.eq.0) then
        write(stt_u,"(3x,i2,a)") zn(io), trim(taux)
     else
        write(stt_u,"(3x,i2,a,i2.2)") zn(io), trim(taux), zx(io)
     endif
  enddo

     write(stt_u,"(/'KPOINTS:',i7)") nkp
     do ik=1,nkp
        write(stt_u,"(3x,3f9.6)") pk(:,ik)
     enddo

  if (dos) then
     write(stt_u,"(/'PDOS CURVES:')")
     do ic=1,ncb
        write(stt_u,"(3x,a)") trim(tit(ic))
        write(stt_u,"(3x,'AO set I: ',/,15x,12i5)") (koc(ic,1,j),j=1,noc(ic,1))
        write(stt_u,"(3x,'Number of set II orbs:',i8)") noc(ic,2)
     enddo
  endif  ! dos

  if (coop) then
     write(stt_u,"(/'COOP CURVES:')")
     do ic=1,ncb
        write(stt_u,"(3x,a)") trim(tit(ic))
        write(stt_u,"(3x,'Distance range:',2f12.4)") dtc(ic,:)
        write(stt_u,"(3x,'AO set I: ',/,15x,12i5)") (koc(ic,1,j),j=1,noc(ic,1))
        write(stt_u,"(3x,'AO set II:',/,15x,12i5)") (koc(ic,2,j),j=1,noc(ic,2))
     enddo
  endif

  close(stt_u)

  !==================================

  if (ref_line_given) then
     allocate(ref_mask(no_u))
     print *, "Orbital set spec: ", trim(ref_line)
     call get_orbital_set(ref_line,ref_mask)
     do io=1, no_u
        if (ref_mask(io)) write(6,fmt="(i5)",advance="no") io
     enddo
     deallocate(ref_mask)
     write(6,*)
     STOP "bye from ref_line processing"
  endif


  !=====================
 
 if (coop) then
    allocate (coop_vals(npts_energy,nspin_blocks,ncb))
    allocate (cohp_vals(npts_energy,nspin_blocks,ncb))
    coop_vals(:,:,:)=0.d0
    cohp_vals(:,:,:)=0.d0
 endif
 if (dos) then
    allocate (pdos_vals(npts_energy,nspin_blocks,ncb))
    pdos_vals(:,:,:) = 0.0_dp
 endif

  ados(:,:) = 0.0_dp

  !================================================================

  ! * Curves

     ! The first dimension is the number of real numbers per orbital
     ! 1 for real wfs, 2 for complex, and four for the two spinor components
  
     if (non_coll) then
        allocate(wf_single(4,1:no_u))
        allocate(wf(4,1:no_u))
     else
        if (gamma_wfsx) then
           allocate(wf_single(1,1:no_u))
           allocate(wf(1,1:no_u))
        else
           allocate(wf_single(2,1:no_u))
           allocate(wf(2,1:no_u))
        endif
     endif
     allocate (mask2(1:no_u))

     do ic=1,ncb

        no1 = noc(ic,1)
        no2 = noc(ic,2)

        mask2(1:no_u) = .false.
        do i2=1,no2
           io2=koc(ic,2,i2)              ! AO Set II
           mask2(io2) = .true.
        enddo

        ! Create reduced pattern
        ! First pass for checking dimensions

        allocate (num_red(no1))
        do i1=1,no1
           num_red(i1) = 0
           io1=koc(ic,1,i1)              ! AO Set I
           do ii1 = 1,numh(io1)
              ind = listhptr(io1)+ii1
              ii2 = indxuo(listh(ind))      ! Equiv orb in unit cell

              if ( .not. mask2(ii2)) Cycle    ! Is not one of Set II

              ! No distance restrictions for PDOS, just an overlap check
              if (dos) then
                 if ( abs(Sover(ind)) < tol_overlap) CYCLE  ! Not overlapping
              endif
              if (coop) then
                if  ( (dij(ind) < dtc(ic,1)) .or. (dij(ind) > dtc(ic,2)) ) CYCLE
              endif

              num_red(i1) = num_red(i1) + 1
           enddo
        enddo
        allocate (ptr(no1))
        ptr(1)=0
        do i1=2,no1
           ptr(i1)=ptr(i1-1)+num_red(i1-1)
        enddo
        nnz = sum(num_red(1:no1))  

        write(*,"(a,3x,a,2x,a,i6,1x,i12)") 'Curve ', trim(tit(ic)),  &
                                      'Base orbitals and interactions: ', &
                                       no1, nnz

        allocate (list_io2(nnz))
        allocate (list_ind(nnz))

        n_int = 0
        do i1=1,no1
           io1=koc(ic,1,i1)              ! AO Set I
           do ii1 = 1,numh(io1)
              ind = listhptr(io1)+ii1
              ii2 = indxuo(listh(ind))

              if ( .not. mask2(ii2)) Cycle

              ! No distance restrictions for PDOS, just an overlap check
              if (dos) then
                 if ( abs(Sover(ind)) < tol_overlap) CYCLE  ! Not overlapping
              endif
              if (coop) then
                if  ( (dij(ind) > dtc(ic,2)) .or. (dij(ind) < dtc(ic,1)) ) CYCLE
              endif

              n_int = n_int + 1
              list_io2(n_int) = ii2
              list_ind(n_int) = ind
           enddo
        enddo
        if (n_int .ne. nnz) then
           print *, "n_int, nnz:", n_int, nnz
           STOP "mismatch"
        endif


     !Stream over file, without using too much memory

        rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        if (debug) print *, "Number of k-points, spins: ", nkp, nsp
        do ik=1,nkp
           if (debug) print *, "k-point: ", ik
           do is=1,nspin_blocks
              read(wfs_u)
              read(wfs_u)
              read(wfs_u)  number_of_wfns
              if (debug) print *, "  Number of wfns: ", number_of_wfns
              do iw=1,number_of_wfns
                 if (debug) print *, "     wfn: ", iw
                 read(wfs_u) 
                 read(wfs_u) eigval

                 ! Early termination of iteration if outside range
                 if (eigval < min_eigval .or. eigval > max_eigval) then
                    read(wfs_u)   ! Still need to read this
                    CYCLE
                 endif

                 ! Use only the specified band set
                 if ( (iw<min_band) .or. (iw>max_band)) then
                    read(wfs_u)   ! Still need to read this
                    CYCLE
                 endif

                 read(wfs_u) (wf_single(:,io), io=1,no_u)
                 ! Use a double precision form in what follows
                 wf(:,:) = real(wf_single(:,:), kind=dp)

                 ! This block will be repeated for every curve,
                 ! but we will divide by the number of curves before writing out
                    do i = 1, npts_energy
                       energy = low_e + e_step*(i-1)
                       ww(i) = delta(energy-eigval)
                       ados(i,is) = ados(i,is) + wk(ik) * ww(i)
                    enddo
                    ! 
                    ! Find the interesting energy region for this state
                    !
                    imin = npts_energy
                    do i = 1, npts_energy
                       if (ww(i) > tol_weight) then
                          imin = i
                          exit
                       endif
                    enddo
                    imax = 1
                    do i =  npts_energy, 1, -1
                       if (ww(i) > tol_weight) then
                          imax = i
                          exit
                       endif
                    enddo

                 do i1 = 1, no1
                    io1=koc(ic,1,i1)              ! AO Set I

                      do i2 = 1,num_red(i1)
                         ind_red = ptr(i1)+i2
                         io2 = list_io2(ind_red)
                         ind = list_ind(ind_red)

                                ! (qcos, qsin) = conjg(C_1) * [S or H] * C_2
                                ! We might want to avoid recomputing the "S" pure wf parts

                                if (non_coll) then

                                   ! Use 'dp' to keep the wfs in double precision
                                   ! (Recall that wf() is now dp; converted right after reading)
                                   spinor_1 = [ cmplx(wf(1,io1),wf(2,io1), dp), &
                                                      cmplx(wf(3,io1),wf(4,io1), dp) ]
                                   spinor_2 = [ cmplx(wf(1,io2),wf(2,io2), dp), &
                                                cmplx(wf(3,io2),wf(4,io2), dp) ]

                                   ! We take the signs for 1,2 and 2,1
                                   ! from the construction of Ebs_Haux in 'compute_energies'
                                   if (h_spin_dim == 8) then
                                      H(1,1) = cmplx(Hamilt(ind,1), Hamilt(ind,5), dp)
                                      H(1,2) = cmplx(Hamilt(ind,3), -Hamilt(ind,4), dp)
                                      H(2,1) = cmplx(Hamilt(ind,7), Hamilt(ind,8), dp)
                                      H(2,2) = cmplx(Hamilt(ind,2), Hamilt(ind,6), dp)
                                   else   ! nsp=4; just non-collinear; no SOC
                                      H(1,1) = cmplx(Hamilt(ind,1), 0.0_dp, dp)
                                      H(1,2) = cmplx(Hamilt(ind,3), -Hamilt(ind,4), dp)
                                      H(2,1) = cmplx(Hamilt(ind,3), Hamilt(ind,4), dp)
                                      H(2,2) = cmplx(Hamilt(ind,2), 0.0_dp, dp)
                                   endif

                                   ! For DOS/COOP, take as weight the "complete spinor" product
                                   ! The dot_product is directly sum(conjg(a1)*a2) for complex arrays
                                   ! The kind returned is "the highest" of the arguments. In this case,
                                   ! it will be 'dp' if the 'spinor_X' variables are 'dp'.
                                   c_c1_c2 = dot_product(spinor_1,spinor_2)
                                   qcos= real(c_c1_c2, dp) 
                                   qsin= aimag(c_c1_c2)   
                                   
                                   ! For COHP, insert non-trivial H matrix
                                   H_c2 = matmul( H, spinor_2 )
                                   c_c1_H_c2  = dot_product( spinor_1, H_c2 )
                                   qcos_H= real(c_c1_H_c2, dp) 
                                   qsin_H= aimag(c_c1_H_c2)

                                else
                                   ! These have explicit spin quantum numbers (is)
                                   if (gamma_wfsx) then
                                      qcos = wf(1,io1)*wf(1,io2) 
                                      qsin = 0.0_dp
                                      qcos_H = qcos * Hamilt(ind,is)
                                      qsin_H = 0.0_dp
                                   else
                                      qcos = (wf(1,io1)*wf(1,io2) + &
                                           wf(2,io1)*wf(2,io2))
                                      qcos_H = qcos * Hamilt(ind,is)
                                      qsin = (wf(1,io1)*wf(2,io2) - &
                                           wf(2,io1)*wf(1,io2))
                                      qsin_H = qsin * Hamilt(ind,is)
                                   endif
                                endif

                                ! k*R_12    (r_2-r_1)
                                alfa=dot_product(pk(1:3,ik),xij(1:3,ind))

                                ! Crb = Real(C_1*conjg(C_2)*exp(-i*alfa)) * S_12
                                ! or    Real(conjg(C_1)*C_2)*exp(+i*alfa)) * S_12
                                ! Common factor computed here
                                factor_S =  Sover(ind) * (qcos*cos(alfa)-qsin*sin(alfa)) * wk(ik)
                                factor_H =  (qcos_H*cos(alfa)-qsin_H*sin(alfa)) * wk(ik)

                                if (dos) then
                                   ! Note that there is no factor of 2
                                   do i = imin, imax
                                      pdos_vals(i,is,ic)=pdos_vals(i,is,ic)  +  factor_S * ww(i)
                                   enddo
                                endif

                                if (coop) then
                                   !! COOP is basically the off-diagonal DOS, chosen on the basis of particular bonds
                                   do i = imin, imax
                                      coop_vals(i,is,ic)=coop_vals(i,is,ic)  +  factor_S * ww(i)
                                      cohp_vals(i,is,ic)=cohp_vals(i,is,ic)  +  factor_H * ww(i)
                                   enddo
                                endif  ! coop

                        enddo   ! i2
                    enddo  ! i1

                 enddo   ! iwf
              enddo      ! is
           enddo         ! ik
           
           deallocate (num_red)
           deallocate (ptr)
           deallocate (list_io2)
           deallocate (list_ind)

     !     Output curves
     !
           if (dos) then
              call write_curve_to_file(trim(mflnm)// "." // trim(tit(ic)) // '.pdos', &
                                       pdos_vals(:,:,ic))
           endif
           if (coop) then
              call write_curve_to_file(trim(mflnm)// "." // trim(tit(ic)) // '.coop', &
                                       coop_vals(:,:,ic))
              call write_curve_to_file(trim(mflnm)// "." // trim(tit(ic)) // '.cohp', &
                                       cohp_vals(:,:,ic))
           endif
           
        enddo    ! ic

!--------------------------------------------------------

     !
     !      Simple DOS output
     !
     ! Divide by the number of curves        
     call write_curve_to_file(trim(sflnm)//".ados",ados/ncb)


CONTAINS

  subroutine write_curve_to_file(filename,array)
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: array(:,:)

    integer  :: i
    real(dp) :: energy
    
    ! tab_u, low_e, e_step, npts_energy by host association
    
    open(tab_u,file=trim(filename))
    !
    ! Header
    !
    select case ( wfs_spin_flag)
    case ( 1 )
       ! Two identical 'spin' columns, and the sum (complete value) in the third column
       write(tab_u,"(a1,12x,'ENERGY',3(14x,a))") '#', 's1', 's2=s1', 'total'
    case ( 2 )
       ! Two 'spin' columns, and the sum (complete value) in the third column
       write(tab_u,"(a1,12x,'ENERGY',3(14x,a))") '#', 's1', 's2', 'total'
    case ( 4 )
       ! A single column with  the complete value
       write(tab_u,"(a1,12x,'ENERGY',14x,a)") '#', 'total'
    end select
    !
    ! Values   
    !
    do i=1, npts_energy
       energy = low_e + e_step*(i-1)
       select case ( wfs_spin_flag)
       case ( 1 )
          ! Two identical 'spin' columns, and the sum (complete value) in the third column
          write(tab_u,"(f20.8,3(5x,f13.8))")  &   
               energy, array(i,1), array(i,1), 2*array(i,1)
       case ( 2 )
          ! Two 'spin' columns, and the sum (complete value) in the third column
          write(tab_u,"(f20.8,3(5x,f13.8))")  &      
               energy, (array(i,is),is=1,nspin_blocks), sum(array(i,:))
       case ( 4 )
          ! A single column with  the complete value
          write(tab_u,"(f20.8,5x,f13.8)") energy, array(i,1)
       end select
    enddo
    close(tab_u)

  end subroutine write_curve_to_file
  
  function delta(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res

    if ((abs(x) > 8*smear)) then
       res = 0.0_dp
       RETURN
    endif

    res = exp(-(x/smear)**2) / (smear*sqrt(pi))

  end function delta
end program mprop


