!==================================================================
program dm_creator

  use main_vars
  use subs, only: manual_dm_creator
  use io_hs, only: read_hs_file
#ifdef CDF
  use iodm_netcdf
#endif

  implicit none

  logical :: gamma_wfsx, got_qcos, logical_dummy
  integer :: ii1, ii2, ind, ind_red, no1, no2, n_int, nnz
  integer :: nwfmx
  real(dp) :: factor, qsol
  real(dp), dimension(:,:), allocatable :: DMout

  ! We use a smearing function of the form f(x) = exp(-(x/smear)**2) / (smear*sqrt(pi))
  ! A weight tolerance of 1.0e-4 corresponds to going about 3*smear on either
  ! side of the eigenvalue. We should set nsigma appropriately

  real(dp), parameter  :: tol_weight = 1.0e-4_dp
  integer, parameter   :: nsigma = 3

  real(dp), parameter  :: tol_overlap = 1.0e-10_dp

  !
  !     Process options
  !
  n_opts = 0
  do
     call getopts('dhls:n:m:M:',opt_name,opt_arg,n_opts,iostat)
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

  allocate (ados(npts_energy,nsp), ww(npts_energy))

  nwfmx = 0
  min_energy = huge(1.0_dp)
  max_energy = -huge(1.0_dp)

  do ik=1,nkp
     do is=1,nsp

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


  ! Here low_e and high_e represent a window for the plot, to
  ! avoid cut tails

  low_e = min_energy - nsigma*smear
  high_e = max_energy + nsigma*smear

  e_step = (high_e-low_e)/(npts_energy-1)
  ados(:,1:nsp) = 0.0_dp

  ! skip four records

  rewind(wfs_u)

  read(wfs_u) 
  read(wfs_u) 
  read(wfs_u) 
  read(wfs_u) 

  do ik=1,nkp
     do is=1,nsp
        read(wfs_u)
        read(wfs_u)
        read(wfs_u)  number_of_wfns
        do iw=1,number_of_wfns
           read(wfs_u) 
           read(wfs_u) eigval
           do i = 1, npts_energy
              energy = low_e + e_step*(i-1)
              weight = delta(energy-eigval)
              if (weight < tol_weight) CYCLE           ! Will not contribute
              ados(i,is) = ados(i,is) + wk(ik) * weight
           enddo
           read(wfs_u)       ! Skip wfn info
        enddo
     enddo
  enddo

  call io_assign(idos)
  open(idos,file=trim(sflnm)//".alldos",form="formatted", &
       status="unknown",action="write",position="rewind")
  write(idos,*) "#  Energy   LARGE-SCALE DOS"

  do i = 1, npts_energy
     energy = low_e + e_step*(i-1)
     write(idos,*) energy, (ados(i,is), is=1,nsp)
  enddo

  call io_close(idos)


  ! Read HSX file
  ! Will pick up atoms, zval, and thus N_electrons

  call read_hs_file(trim(sflnm)//".HSX")
  if (gamma_wfsx .neqv. gamma) STOP "Gamma mismatch"

  ztot = 0.0_dp
  do ia = 1, na_u
     ztot = ztot + zval(isa(ia))
  enddo

  !
  !       Compute integrated total DOS
  !       Here we double the DOS for the case of nspin=1,
  !       to get the correct number of states.

  allocate(intdos(npts_energy))
  call io_assign(intdos_u)
  open(intdos_u,file=trim(sflnm)//".intdos",form="formatted", &
       status="unknown",action="write",position="rewind")
  intdos(1) = 0.0_dp
  write(intdos_u,*) low_e, intdos(1)
  do i = 2, npts_energy
     energy = low_e + e_step*(i-1)
     intdos(i) = intdos(i-1) + sum(ados(i,:)) * e_step * 2.0_dp /nsp
     write(intdos_u,*) energy, intdos(i)
  enddo
  call io_close(intdos_u)

  ! Look for Fermi Energy
  do i = 2, npts_energy
     if (intdos(i) > ztot) then
        ! Found fermi energy
        energy = low_e + e_step*(i-1)
        ! Correct overshoot, interpolating linearly
        efermi = energy - (intdos(i)-ztot)*e_step/(intdos(i)-intdos(i-1))
        exit
     endif
  enddo

  write(6,"(a,f10.5,a,f10.5)") "Fermi energy: ", efermi, " for a smearing of: ", smear

  if (energies_only) STOP
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

  close(stt_u)

  !==================================

  !================================================================

  ! * Curves

     if (gamma) then
        allocate(wf(1,1:no_u))
     else
        allocate(wf(2,1:no_u))
     endif

     nnz = sum(numh(1:nao))
     allocate(DMout(nnz,nspin))
     DMout(:,:) = 0.0_dp

     !Stream over file, without using too much memory

     rewind(wfs_u)

        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
        read(wfs_u) 
     
        if (debug) print *, "Number of k-points, spins: ", nkp, nsp
        do ik=1,nkp
           if (debug) print *, "k-point: ", ik
           do is=1,nsp
              read(wfs_u)
              read(wfs_u)
              read(wfs_u)  number_of_wfns
              if (debug) print *, "  Number of wfns: ", number_of_wfns
              do iw=1,number_of_wfns
                 if (debug) print *, "     wfn: ", iw
                 read(wfs_u) 
                 read(wfs_u) eigval
                 ! Early termination of iteration
                 ! Note that we keep a few more states on the sides, due to
                 ! the smearing
                 if (eigval < min_energy .or. eigval > max_energy) then
                    read(wfs_u)   ! Still need to read this
                    CYCLE
                 endif

                 read(wfs_u) (wf(:,io), io=1,nao)

                 do io1 = 1, nao
                    do ii1 = 1,numh(io1)
                       ind = listhptr(io1)+ii1
                       ii2 = listh(ind)
                       io2 = indxuo(ii2)

                             
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
                             factor =  (qcos*cos(alfa)-qsin*sin(alfa)) * wk(ik)

                             DMout(ind,is) = DMout(ind,is) + factor

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

           qsol = 0.0_dp
           do is = 1,nspin
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

!--------------------------------------------------------


CONTAINS

  function delta(x) result(res)
    real(dp), intent(in) :: x
    real(dp)             :: res

    if ((abs(x) > 8*smear)) then
       res = 0.0_dp
       RETURN
    endif

    res = exp(-(x/smear)**2) / (smear*sqrt(pi))

  end function delta
end program dm_creator




