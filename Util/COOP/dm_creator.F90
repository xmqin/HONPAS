! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!==================================================================
program dm_creator

  use main_vars
  use subs, only: manual_dm_creator
  use io_hs, only: read_hs_file
#ifdef CDF
  use iodm_netcdf
#endif

  implicit none

  logical :: gamma_wfsx, got_qcos, logical_dummy, non_coll
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  integer :: nwfmx, nspin_blocks
  real(dp) :: factor, qsol
  real(dp), dimension(:,:), allocatable :: DMout

  complex(DP), dimension(2)    :: spinor_1, spinor_2
  complex(DP) :: d11, d12, d21, d22, kphs
  
  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhm:M:',opt_name,opt_arg,n_opts,iostat)
     if (iostat /= 0) exit
     select case(opt_name)
     case ('d')
        debug = .true.
     case ('+d')
        debug = .false.
     case ('m')
        read(opt_arg,*) minimum_spec_energy
     case ('M')
        read(opt_arg,*) maximum_spec_energy
     case ('h')
        call manual_dm_creator()
     case ('?',':')
        write(0,*) "Invalid option: ", opt_arg(1:1)
        write(0,*) "Use -h option for manual"
        STOP
     end select
  enddo

  nargs = command_argument_count()
  nlabels = nargs - n_opts + 1
  if (nlabels /= 1)  then
     write(0,*) "Usage: dm_creator [ -d ] [ -h ] [-m MIN_ENERGY] [-M MAX_ENERGY] SystemLabel"
     write(0,*) "Use -h option for manual"
     STOP
  endif

  call get_command_argument(n_opts,value=sflnm,status=iostat)
  if (iostat /= 0) then
     STOP "Cannot get SystemLabel"
  endif

  !==================================================

  ierr=0

  !==================================================
  ! Read WFSX file

  write(6,"(1x,a,'.WFSX ...')") trim(sflnm)

  open(wfs_u,file=trim(sflnm)//'.WFSX',status='old',form='unformatted')
  read(wfs_u) nkp, gamma_wfsx
  allocate (wk(nkp), pk(3,nkp))

  read(wfs_u) nsp
  read(wfs_u) nao
  read(wfs_u)        !! Symbols, etc
  if (debug) print *, "WFSX read: nkp, nsp, nnao: ", nkp, nsp, nao

  if ( nsp == 8 ) then
    if (debug) print *, "WFSX from spin-orbit calculation"
    nspin_blocks = 1
    non_coll = .true.
  else if ( nsp == 4 ) then
    if (debug) print *, "WFSX from non-collinear calculation"
    nspin_blocks = 1
    non_coll = .true.
  else
     nspin_blocks = nsp
  endif
  print *, "NSPIN BLOCKS: ", nspin_blocks

  nwfmx = 0
  min_energy = huge(1.0_dp)
  max_energy = -huge(1.0_dp)

  do ik=1,nkp
     do is=1,nspin_blocks

        read(wfs_u) idummy, pk(1:3,ik), wk(ik)
        if (idummy /= ik) stop "ik index mismatch in WFS file"
        read(wfs_u) is0
        read(wfs_u) number_of_wfns
        nwfmx = max(nwfmx,number_of_wfns)

        do iw=1,number_of_wfns
           read(wfs_u) iw0
           read(wfs_u) eigval
           min_energy = min(min_energy,eigval)
           max_energy = max(max_energy,eigval)
           read(wfs_u)
        enddo
     enddo
  enddo

  print *, " Maximum number of wfs per k-point: ", nwfmx
  print "(a,2f12.4)", "Min_energy, max_energy on WFS file: ",  &
       min_energy, max_energy


  ! Read HSX file
  ! Will pick up atoms, zval, and thus N_electrons

  call read_hs_file(trim(sflnm)//".HSX")

  ztot = 0.0_dp
  do ia = 1, na_u
     ztot = ztot + zval(isa(ia))
  enddo

  !-------------------------------------------------------------------

  if (minimum_spec_energy > min_energy)   min_energy = minimum_spec_energy
  if (maximum_spec_energy < max_energy)  max_energy = maximum_spec_energy

  print "(a,2f12.4)", "Min_energy, max_energy used for DM: ",  &
       min_energy, max_energy

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
  write(6,"('Writing files: ',a,'.stt ...')") trim(sflnm)
  open(stt_u,file=trim(sflnm)//'.stt')
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
  
     nnz = sum(numh(1:nao))
     allocate(DMout(nnz,nspin))
     DMout(:,:) = 0.0_dp

     !===============================================
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
                 if (debug) print *, "   eigval: ", eigval
                 ! Early termination of iteration
                 if (eigval < min_energy .or. eigval > max_energy) then
                    read(wfs_u)   ! Still need to read this
                    if (debug) print *, "-------------- skipped"
                    CYCLE
                 endif

                 read(wfs_u) (wf_single(:,io), io=1,nao)
                 wf(:,:) = real(wf_single(:,:), kind=dp)
                 
                 do io1 = 1, nao
                    do ii1 = 1,numh(io1)
                       ind = listhptr(io1)+ii1
                       ii2 = listh(ind)
                       io2 = indxuo(ii2)

                       ! k*R_12    (r_2-r_1)
                       alfa=dot_product(pk(1:3,ik),xij(1:3,ind))
                       kphs = exp(cmplx(0.0_dp,-1.0_dp, dp) * alfa )

                       if (non_coll) then

                          ! Use 'dp' to keep the wfs in double precision                        
                          ! (Recall that wf() is now dp; converted right after reading)         

                          spinor_1 = [ cmplx(wf(1,io1),wf(2,io1), dp), &
                               cmplx(wf(3,io1),wf(4,io1), dp) ]
                          spinor_2 = [ cmplx(wf(1,io2),wf(2,io2), dp), &
                               cmplx(wf(3,io2),wf(4,io2), dp) ]

                          d11 = spinor_1(1) * conjg(spinor_2(1)) * kphs * wk(ik)
                          d12 = spinor_1(1) * conjg(spinor_2(2)) * kphs * wk(ik)
                          d21 = spinor_1(2) * conjg(spinor_2(1)) * kphs * wk(ik)
                          d22 = spinor_1(2) * conjg(spinor_2(2)) * kphs * wk(ik)

                          if (nspin == 8) then
                             ! taken from diag3k
                             dmout(ind,1) = dmout(ind,1) + real(d11, dp)
                             dmout(ind,2) = dmout(ind,2) + real(d22, dp)
                             dmout(ind,3) = dmout(ind,3) + real(d12, dp)
                             dmout(ind,4) = dmout(ind,4) - aimag(d12)
                             dmout(ind,5) = dmout(ind,5) + aimag(d11)
                             dmout(ind,6) = dmout(ind,6) + aimag(d22)
                             dmout(ind,7) = dmout(ind,7) + real(d21, dp)
                             dmout(ind,8) = dmout(ind,8) + aimag(d21)

                          else ! nspin = 4
                             ! taken from diag2k
                             D12 = 0.5_dp * (D12 + dconjg(D21))
                             dmout(ind,1) = dmout(ind,1) + real(d11, dp)
                             dmout(ind,2) = dmout(ind,2) + real(d22, dp)
                             dmout(ind,3) = dmout(ind,3) + real(d12, dp)
                             dmout(ind,4) = dmout(ind,4) - aimag(d12)
                          endif
                          
                       else
                                !AG: Corrected:  (qcos, qsin) = conjg(C_1)*(C_2)

                                if (gamma_wfsx) then
                                   qcos = wf(1,io1)*wf(1,io2) 
                                   qsin = 0.0_dp
                                else
                                   qcos= (wf(1,io1)*wf(1,io2) + &
                                        wf(2,io1)*wf(2,io2))
                                   qsin= (wf(1,io1)*wf(2,io2) - &
                                        wf(2,io1)*wf(1,io2))
                                endif

                             ! Real(conjg(C_1)*C_2)*exp(+i*alfa)) * S_12
                             ! Common factor computed here
                             factor =  (qcos*cos(alfa)-qsin*sin(alfa)) * wk(ik)
                             DMout(ind,is) = DMout(ind,is) + factor
                             
                          endif  ! coll or non-coll
                          
                        enddo   ! i2
                    enddo  ! i1

                 enddo   ! iwf
              enddo      ! is
           enddo         ! ik

!
!          Factor of two for non-spin-polarized cases
!
           if (nsp == 1)  DMout(:,:) = DMout(:,:) * 2.0_dp 

!
! Compute total number of electrons
           !
           ! ** UPDATE THIS FOR SPINORS

           qsol = 0.0_dp
           do is = 1,min(nspin,2)
              do io = 1,nnz
                 qsol = qsol + DMout(io,is) * Sover(io)
              enddo
           enddo
           print *, "Total number of electrons in DM: ", qsol
           print *, "Total number of electrons in system: ", ztot

           
#ifdef CDF
           call setup_dm_netcdf_file( nnz, nao, nspin,    &
                                 no_s, indxuo,            &
                                 numh,  listhptr, listh)
           call write_dm_netcdf( nao, nnz, nspin, DMout, overwrite = .true. )
#endif
           call iodm( "write", nnz, nao, nspin, numh, listhptr, listh, DMout, logical_dummy)

end program dm_creator




