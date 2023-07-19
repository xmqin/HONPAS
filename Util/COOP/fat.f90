!==================================================================
!
program fatband
!
!
! Computes the generalized projection of each eigenvector on a given
! orbital set.
! The syntax for orbital sets is the same as in mprop, so this
! program can re-use as description file the same DOS-type file
! as mprop.
!
!

  use main_vars
  use orbital_set, only: get_orbital_set
  use io_hs, only: read_hs_file
  use read_curves, only: read_curve_information, mask_to_arrays

  implicit none

  logical :: gamma_wfsx, got_qcos
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  real(dp) :: factor

  real(dp), parameter  :: tol_overlap = 1.0e-10_dp

  logical, allocatable   :: mask2(:)
  integer, allocatable   :: num_red(:), ptr(:), list_io2(:), list_ind(:)
  real(dp), allocatable  :: eig(:,:), fat(:,:)


  integer  :: nwfmx, nwfmin
  integer  :: min_band = 1
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

  integer  :: fat_u, ib, nbands

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhls:n:m:M:R:b:B:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('l')
        energies_only = .true.
     case ('R')
        ref_line_given = .true.
        ref_line = opt_arg
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('h')
        call manual()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Usage: fat [ options ] DescriptionFile "
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: fat [ options ] DescriptionFile "
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=mflnm,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get .mpr file root"
  endif

  band_interval_set = (min_band_set .or. max_band_set)

  !==================================================

  ierr=0

  ! Read type of job

  open(mpr_u,file=trim(mflnm) // ".mpr", status='old')
  read(mpr_u,*) sflnm
  read(mpr_u,*) what
  if (trim(what) /= "DOS") STOP "Fatbands needs DOS-style jobfile"
  
  !==================================================
  ! Read WFSX file

  write(6,"(a)") "Reading wf file: " // trim(sflnm) // ".WFSX"

  open(wfs_u,file=trim(sflnm)//'.WFSX',status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) nsp
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "WFSX read: nkp, nsp, nnao: ", nkp, nsp, nao

  nwfmx = -huge(1)
  nwfmin = huge(1)
  min_eigval = huge(1.0_dp)
  max_eigval = -huge(1.0_dp)
  min_eigval_in_band_set = huge(1.0_dp)
  max_eigval_in_band_set = -huge(1.0_dp)

  do ik=1,nkp
     do is=1,nsp

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

  print "(a,2i5)", " Minimum/Maximum number of wfs per k-point: ", nwfmin, nwfmx
  print "(a,2f12.4)", "Min_eigval, max_eigval on WFS file: ",  &
       min_eigval, max_eigval
  print "(a,2f12.4)", "Min_eigval, max_eigval in band set : ",  &
          min_eigval_in_band_set, max_eigval_in_band_set

  if (.not. band_interval_set) then

     min_band = 1      ! Already set by default
     max_band = nwfmin

  else

     if (min_band_set .and. (min_band < 1)) then
        print "(a)", " ** Min_band implicitly reset to 1..."
        min_band = 1
     endif
     if (min_band_set .and. (min_band > nwfmin)) then
        print "(a,2i5)", " ** Min_band is too large for some k-points: (min_band, nwfmin):", min_band, nwfmin
        STOP
     endif
     if (max_band_set .and. (max_band > nwfmin)) then
        print "(a,2i5)", " ** Max_band is too large for some k-points: (max_band, nwfmin):", max_band, nwfmin
        print "(a)", " ** Max_band will be effectively reset to its maximum allowed value"
        max_band = nwfmin
     endif
     if (max_band_set .and. (max_band < min_band)) then
        print "(a,2i5)", " ** Max_band is less than min_band: (max_band, min_band):", max_band, min_band
        STOP
     endif

     min_eigval = min_eigval_in_band_set
     max_eigval = max_eigval_in_band_set
  endif
  print "(a,3i4)", "Band set used: (min, max):",  min_band, max_band

  ! Read HSX file
  ! Will pick up atoms, zval, and thus the nominal number of electrons,
  ! but the total charge is read as qtot.

  call read_hs_file(trim(sflnm)//".HSX")
  if (gamma_wfsx .neqv. gamma) STOP "Gamma mismatch"

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

  ! Give values to flags for reuse of the reading routine
  dos = .true.
  coop = .false.
  call read_curve_information(.true.,.false.,  &
                                    mpr_u,no_u,ncbmx,ncb,tit,orb_mask,dtc)

  orb_mask(:,2,1:ncb) = .true.       ! All orbitals considered

  call mask_to_arrays(ncb,orb_mask(:,1,:),noc(:,1),koc(:,1,:))
  call mask_to_arrays(ncb,orb_mask(:,2,:),noc(:,2),koc(:,2,:))

!!
  write(6,"('Writing files: ',a,'.stt ...')") trim(mflnm)
  open(stt_u,file=trim(mflnm)//'.info')
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

     write(stt_u,"(/'FATBAND ORBITAL SETS:')")
     do ic=1,ncb
        write(stt_u,"(3x,a)") trim(tit(ic))
        write(stt_u,"(3x,'AO set I: ',/,15x,12i5)") (koc(ic,1,j),j=1,noc(ic,1))
        write(stt_u,"(3x,'Number of set II orbs:',i8)") noc(ic,2)
     enddo

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
 
  ! * Fatband weights

     nbands = max_band - min_band + 1
     allocate(eig(nbands,nspin), fat(nbands,nspin))

     if (gamma) then
        allocate(wf(1,1:no_u))
     else
        allocate(wf(2,1:no_u))
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

              ! No distance restrictions for PDOS or FATBANDS, just an overlap check

                 if ( abs(Sover(ind)) < tol_overlap) CYCLE  ! Not overlapping

              num_red(i1) = num_red(i1) + 1
           enddo
        enddo
        allocate (ptr(no1))
        ptr(1)=0
        do i1=2,no1
           ptr(i1)=ptr(i1-1)+num_red(i1-1)
        enddo
        nnz = sum(num_red(1:no1))  

        write(*,"(a,3x,a,2x,a,i6,1x,i12)") 'Fatband coeffs set: ', trim(tit(ic)),  &
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

                 if ( abs(Sover(ind)) < tol_overlap) CYCLE  ! Not overlapping

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
     
        open(fat_u,file=trim(mflnm)// "." // trim(tit(ic)) // '.EIGFAT')
        write(fat_u,"(a,2i5)") "# " // trim(sflnm) // " min_band, max_band: ", min_band, max_band
        write(fat_u,"(3i6)")   nbands, min(nsp,2), nkp

        do ik=1,nkp

           write(fat_u,"(i4,3(1x,f10.5))")  ik, pk(1:3,ik)

           do is=1,nsp

              ib = 0
              fat(:,is) = 0.0_dp

              read(wfs_u) 
              read(wfs_u) 
              read(wfs_u)  number_of_wfns
              do iw=1,number_of_wfns
                 read(wfs_u) 
                 read(wfs_u) eigval

                 ! Use only the specified band set
                 if ( (iw<min_band) .or. (iw>max_band)) then
                    read(wfs_u)   ! Still need to read this
                    CYCLE
                 endif

                 ib = ib + 1
                 eig(ib,is) = eigval    ! This will be done for every curve... harmless
                 read(wfs_u) (wf(:,io), io=1,nao)


                 do i1 = 1, no1
                    io1=koc(ic,1,i1)              ! AO Set I

                      do i2 = 1,num_red(i1)
                         ind_red = ptr(i1)+i2
                         io2 = list_io2(ind_red)
                         ind = list_ind(ind_red)
                             
                                ! (qcos, qsin) = C_1*conjg(C_2)
                                !AG: Corrected:  (qcos, qsin) = conjg(C_1)*(C_2)
                                ! We might want to avoid recomputing this

                                if (gamma) then
                                   qcos = wf(1,io1)*wf(1,io2) 
                                   qsin = 0.0_dp
                                else
                                   qcos= (wf(1,io1)*wf(1,io2) + &
                                        wf(2,io1)*wf(2,io2))
                                   qsin= (wf(1,io1)*wf(2,io2) - &
                                        wf(2,io1)*wf(1,io2))
                                endif

                             ! k*R_12    (r_2-r_1)
                             alfa=dot_product(pk(1:3,ik),xij(1:3,ind))

                             ! Crb = Real(C_1*conjg(C_2)*exp(-i*alfa)) * S_12
                             !AG: This one better --  or Real(conjg(C_1)*C_2)*exp(+i*alfa)) * S_12
                             ! Common factor computed here
                             factor =  (qcos*cos(alfa)-qsin*sin(alfa))

                             fat(iw,is) = fat(iw,is) + Sover(ind)*factor

                        enddo   ! i2
                    enddo  ! i1

                 enddo   ! iwf
                 write(fat_u,"(4(4x,f10.4,f9.5))")   &
                      (eig(ib,is),fat(ib,is),ib=1,nbands)
              enddo      ! is

           enddo         ! ik
           
           deallocate (num_red)
           deallocate (ptr)
           deallocate (list_io2)
           deallocate (list_ind)

              close(fat_u)

        enddo    ! ic

 CONTAINS

      subroutine manual()

      write(6,"('* FAT(BANDS) PROGRAM')")
      write(6,"('  Alberto Garcia, ICMAB-CSIC, 2012 ')")
      write(6,*)
      write(6,"('    FAT calculates eigenvector projections ')")
      write(6,"('    using output files obtained with SIESTA. The atomic orbital (AO)')")
      write(6,"('    sets are defined in an input file (MLabel.mpr).')")
      write(6,"('  ')")
      write(6,*) "Usage: fat [ options ] MPROP_FILE_BASENAME"
      write(6,*) "Options:"
      write(6,*) "           -h:  print manual                    "
      write(6,*) "           -d:  debug                    "
      write(6,*) "           -l:  print summary of energy information         "
      write(6,*) "    "
      write(6,*) "   Selection of eigenstates to be used: "
      write(6,*) "    "
      write(6,*) "   -b Min_band  :  set minimum band index to be used               "
      write(6,*) "   -B Max_band  :  set maximum band index to be used               "
      write(6,*) "    "
      write(6,*)
      write(6,"('* .mpr FILE STRUCTURE')")
      write(6,"('         SLabel                   # Name of the siesta output files')")
      write(6,"('         DOS                      # DOS option of mprop is mandatory')")
      write(6,"('    /-  As many blocks as projections wanted ]')")
      write(6,"('    |    projection_name         # DOS projection name')")
      write(6,"('    \-   Subset of AO (*)        # Subset of orbitals included')")
      write(6,"('     (*) See below how to define subsets of AO')")
      write(6,"('     A final line with leading chars  ----  can signal the end of the input')")
      write(6,*)
      write(6,"('* INPUT FILES')")
      write(6,"('    [output files from SIESTA >=  2.4.1]')")
      write(6,"('    SLabel.WFSX and SLabel.HSX (new format)')")
      write(6,*)
      write(6,"('* OUTPUT FORMAT')")
      write(6,*) 
      write(6,*) " MLabel.CurveName.EIGFAT    :  File with eigenvalue and projection info "
      write(6,"('    [An information file with extension .info will always be generated]')")
      write(6,*)
      write(6,"('* PROJECTION AND CURVES NAMES')")
      write(6,"('    Alphanumerical string up to 30 char. with no spaces')")
      write(6,"('* SUBSET OF AO USING ORDER NUMBERS')")
      write(6,"('    List of integer numbers preceeded by a + symbol')")
      write(6,"('    Each number refers to one AO in the final list of AO of SIESTA')")
      write(6,"('    Example: + 23 65 78')")
      write(6,"('* SUBSET OF AO USING ATOM_SHELL NOTATION')")
      write(6,"('    List of atoms and shell groups of AO')")
      write(6,"('    General notation: ATOM_SHELL')")
      write(6,"('     > ATOM:  Atomic symbol refers to all the atoms of that type')")
      write(6,"('              Integer number refers to the N-th atom in unit cell')")
      write(6,"('     > SHELL: Integer1+Letter+Integer2')")
      write(6,"('               > Integer1 refers to the n quantum number')")
      write(6,"('               > Letter   refers to the l quantum number (s,p,d,f,g,h)')")
      write(6,"('               > Integer2 refers to a single AO into the n-l shell')")
      write(6,"('                   Alternatively, alphanumerical strings can be used')")
      write(6,"('                     p-shells   1  y    d-shells   1  xy   4  xz')")
      write(6,"('                                2  z               2  yz   5  x2-y2')")
      write(6,"('                                3  x               3  z2')")
      write(6,"('    Particular cases:')")
      write(6,"('     > Just ATOM is indicated: all the AO of the atom will be included')")
      write(6,"('     > No value for Integer2:  all the AO of the shell will be included')")
      write(6,"('    Example: Ca_3p Al 4_4d3 5 O_2py')")
      stop

      end subroutine manual


end program fatband


