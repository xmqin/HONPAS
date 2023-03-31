! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!==================================================================
program spin_texture

  ! Computes the spin texture from a set of spinor wavefunctions.
  ! Alberto Garcia, February 2019
  ! Based on 'dm_creator' and code from Roberto Robles.
  
  use main_vars
  use subs, only: manual_spin_texture
  use io_hs, only: read_hs_file

  implicit none

  logical :: gamma_wfsx, got_qcos, logical_dummy, non_coll
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  integer :: nwfmx, nwfmin, nspin_blocks
  real(dp) :: factor, qsol
  character(len=256) :: file_prefix

  integer  :: min_band = -huge(1)
  integer  :: max_band =  huge(1)
  logical  :: min_band_set = .false.
  logical  :: max_band_set = .false.
  logical  :: band_interval_set = .false.
  real(dp) :: min_eigval_in_band_set
  real(dp) :: max_eigval_in_band_set
  real(dp) :: maximum_spec_eigval = huge(1)
  real(dp) :: minimum_spec_eigval = -huge(1)
  real(dp) :: min_eigval
  real(dp) :: max_eigval

  real(dp) :: st(0:3)
  complex(DP), dimension(2)    :: spinor_1, spinor_2
  complex(DP) :: d11, d12, d21, d22, kphs
  
  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhm:M:b:B:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('+d')
        debug = .false.
     case ('b')
        read(opt_arg,*) min_band
        min_band_set = .true.
     case ('B')
        read(opt_arg,*) max_band
        max_band_set = .true.
     case ('m')
        read(opt_arg,*) minimum_spec_eigval
     case ('M')
        read(opt_arg,*) maximum_spec_eigval
     case ('h')
        call manual_spin_texture()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: spin_texture [ -d ] [ -h ] " // &
                "                    [-m MIN_EIGVAL] [-M MAX_EIGVAL]" // &
                "                    [-b MIN_BAND] [-M MAX_BAND]   SystemLabel"
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=file_prefix,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get SystemLabel"
  endif

  band_interval_set = (min_band_set .or. max_band_set)

  !==================================================

  ierr=0

  !==================================================
  ! Read WFSX file

  write(0,"(1x,a,'.WFSX ...')") trim(file_prefix)

  open(wfs_u,file=trim(file_prefix)//'.WFSX',status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) nsp
  non_coll = (nsp >= 4)
  if (.not. (non_coll)) then
     STOP "Spin texture not available for collinear-spin"
  endif
  
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "WFSX read: nkp, nsp, nnao: ", nkp, nsp, nao

  if (non_coll) then
     nspin_blocks = 1
  else
     nspin_blocks = nsp
  endif
  if (debug) print *, "NSPIN BLOCKS: ", nspin_blocks

  nwfmx = 0
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
           if ((iw>=min_band).and.(iw<=max_band)) then
              min_eigval_in_band_set = min(min_eigval_in_band_set,eigval)
              max_eigval_in_band_set = max(max_eigval_in_band_set,eigval)
           endif
           read(wfs_u)
        enddo
     enddo
  enddo

  write(0,*) "Minimum/Maximum number of wfs per k-point: ", nwfmin, nwfmx
  write(0,"(a,2f12.4)") "min_eigval, max_eigval on WFS file: ",  &
       min_eigval, max_eigval
  write(0,"(a,2f12.4)") "Min_eigval, max_eigval in band set : ",  &
          min_eigval_in_band_set, max_eigval_in_band_set

  if (band_interval_set) then

     if (min_band_set .and. (min_band < 1)) then
        write(0,"(a)") " ** Min_band implicitly reset to 1..."
     endif
     if (min_band_set .and. (min_band > nwfmin)) then
        write(0,"(a,2i5)") " ** Min_band is too large for some k-points: (min_band, nwfmin):", min_band, nwfmin
        STOP
     endif
     if (max_band_set .and. (max_band > nwfmin)) then
        write(0,"(a,2i5)") " ** Max_band is too large for some k-points: (max_band, nwfmin):", max_band, nwfmin
        write(0,"(a)") " ** Max_band will be effectively reset to its maximum allowed value"
     endif
     if (max_band_set .and. (max_band < max(1,min_band))) then
        write(0,"(a,2i5)") " ** Max_band is less than the effective min_band: "// &
             "(max_band, eff min_band):", max_band, max(1,min_band)
        STOP
     endif

     min_eigval = min_eigval_in_band_set
     max_eigval = max_eigval_in_band_set
  endif
  
  write(0,"(a,3i4)") "Implicit band set used: (min, max_min, max_max):",  &
                   max(1,min_band), min(nwfmin,max_band), min(nwfmx,max_band)

  ! min_eigval, max_eigval: Determine which eigenstates are used to compute the curves
  write(0,"(a,2f12.4)") "Minimum and maximum eigenvalues " // &
       "(based on file data and band selection): ",  min_eigval, max_eigval

  if (minimum_spec_eigval > min_eigval) then
     min_eigval = minimum_spec_eigval
     write(0,"(a,f12.4)") "* Minimum eigenvalue changed as per user range request: ",  min_eigval
  endif
  if (maximum_spec_eigval < max_eigval) then
     max_eigval = maximum_spec_eigval
     write(0,"(a,f12.4)") "* Maximum eigenvalue changed as per user range request: ",  max_eigval
  endif

  ! Sanity checks
  write(0,"(a,2f12.4)") "Minimum and maximum eigenvalues to be processed: ",  min_eigval, max_eigval
  if (min_eigval > max_eigval) STOP "Meaningless range. Check -b/-B and -m/-M options"


  ! Read HSX file
  ! Will pick up atoms, zval, and thus N_electrons

  call read_hs_file(trim(file_prefix)//".HSX")

  ztot = 0.0_dp
  do ia = 1, na_u
     ztot = ztot + zval(isa(ia))
  enddo

  !-------------------------------------------------------------------

  if (minimum_spec_eigval > min_eigval)   min_eigval = minimum_spec_eigval
  if (maximum_spec_eigval < max_eigval)  max_eigval = maximum_spec_eigval

  write(6,"(a,2f12.4)") "# min_eigval, max_eigval considered: ",  &
       min_eigval, max_eigval

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


!!
  write(0,"('Writing files: ',a,'.stt ...')") trim(file_prefix)
  open(stt_u,file=trim(file_prefix)//'.stt')
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

  close(stt_u)


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
  
     !===============================================
     !Stream over file, without using too much memory

     rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        if (debug) write(0,*) "Number of k-points, spins: ", nkp, nsp
        do ik=1,nkp

           write(6,'(/A,i4,A,3f12.6,A)')   &
                 "k-point ", ik, " (",(pk(j,ik),j=1,3)," ) (Bohr^-1)"
           write(6,'(A)') "     ie     e(eV)      Sx      Sy      Sz  "

           do is=1,nspin_blocks
              read(wfs_u)
              read(wfs_u)
              read(wfs_u)  number_of_wfns
              if (debug) write(0,*) "  Number of wfns: ", number_of_wfns
              do iw=1,number_of_wfns
                 if (debug) write(0,*) "     wfn: ", iw
                 read(wfs_u) 
                 read(wfs_u) eigval
                 if (debug) write(0,*) "   eigval: ", eigval
                 ! Early termination of iteration
                 if (eigval < min_eigval .or. eigval > max_eigval) then
                    read(wfs_u)   ! Still need to read this
                    if (debug) write(0,*) "---- skipped: not in e range"
                    CYCLE
                 endif
                 if (iw < min_band .or. iw > max_band) then
                    read(wfs_u)   ! Still need to read this
                    if (debug) write(0,*) "---- skipped: not in band range"
                    CYCLE
                 endif

                 read(wfs_u) (wf_single(:,io), io=1,nao)
                 wf(:,:) = real(wf_single(:,:), kind=dp)

                 st(:) = 0.0_dp
                 do io1 = 1, nao
                    do ii1 = 1,numh(io1)
                       ind = listhptr(io1)+ii1
                       ii2 = listh(ind)
                       io2 = indxuo(ii2)

                       ! k*R_12    (r_2-r_1)
                       alfa=dot_product(pk(1:3,ik),xij(1:3,ind))
                       kphs = exp(cmplx(0.0_dp,-1.0_dp, dp) * alfa )

                          ! Use 'dp' to keep the wfs in double precision                        
                          ! (Recall that wf() is now dp; converted right after reading)         

                       spinor_1 = [ cmplx(wf(1,io1),wf(2,io1), dp), &
                            cmplx(wf(3,io1),wf(4,io1), dp) ]
                       spinor_2 = [ cmplx(wf(1,io2),wf(2,io2), dp), &
                            cmplx(wf(3,io2),wf(4,io2), dp) ]

                       d11 = spinor_1(1) * conjg(spinor_2(1)) * kphs
                       d12 = spinor_1(1) * conjg(spinor_2(2)) * kphs
                       d21 = spinor_1(2) * conjg(spinor_2(1)) * kphs
                       d22 = spinor_1(2) * conjg(spinor_2(2)) * kphs

                       ! Recall: dm(:,3) = real(d12);  dm(:,4) = -aimag(d12),
                       ! so this matches with the convention in 'spnvec'
                       ! st(0) is the "charge";  x:1, y:2, z:3: spin components
                       ! (this is shifted with respect to the original code by R. Robles)
                       ! Units: Bohr magnetons
                       D12 = 0.5_dp * (D12 + dconjg(D21))
                       st(0) = st(0) + Sover(ind) * (real(D11,dp) + real(D22,dp))
                       st(1) = st(1) + Sover(ind) * 2.0_dp * real(D12,dp)
                       st(2) = st(2) - Sover(ind) * 2.0_dp * aimag(D12)
                       st(3) = st(3) + Sover(ind) * (real(D11,dp) - real(D22,dp))
                          
                    enddo   ! i2
                 enddo  ! i1
                 !
                 ! Write spin texture to standard output
                 !
                 write(6,'(i7,f12.5,3f8.4)') iw, eigval, (st(j),j=1,3)

              enddo   ! iwf
           enddo      ! is
        enddo         ! ik

      end program spin_texture




