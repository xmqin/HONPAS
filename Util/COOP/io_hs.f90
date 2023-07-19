module io_hs

implicit none
public :: read_hs_file

CONTAINS
subroutine read_hs_file(fname)
  use main_vars
  use precision, only: sp

  character(len=*), intent(in) :: fname

  integer, allocatable  :: ibuff(:)
  real(sp), allocatable  :: hbuff(:)
  real(sp), allocatable  :: buff3(:,:)

  integer numx, ind

  write(6,"(1x,a)",advance='no') trim(fname)
  open(hs_u,file=trim(fname),status='old',form='unformatted')

  read(hs_u,iostat=iostat) nnao, no_s, nspin, nh
  print *, "nnao, no_s, nspin, nh:",  nnao, no_s, nspin, nh
  if (iostat /= 0) STOP "nnao, no_s..."
  if (nnao /= nao) STOP "norbs inconsistency"
  no_u = nao

  read(hs_u,iostat=iostat) gamma
  if (iostat /= 0) STOP "gamma"
  IF (DEBUG) PRINT *, "GAMMA=", gamma
  if (.not. gamma) then
     allocate(indxuo(no_s))
     read(hs_u) (indxuo(i),i=1,no_s)
  else
     allocate(indxuo(no_u))
     do i=1,no_u
        indxuo(i) = i
     enddo
  endif

  if (debug) print *, "HS read: nh, nsp, nnao: ", nh, nspin, nnao
  if (nnao.ne.nao) STOP " nnao .ne. nao in HS"
  if (wfs_x.and.(nspin.ne.nsp)) STOP " nspin .ne. nsp in HS"
  nsp=nspin
  allocate (numh(nao), listhptr(nao), listh(nh))

       allocate (hamilt(nh,nspin))
       allocate (Sover(nh))
       allocate (xij(3,nh),dij(nh))

  read(hs_u,iostat=iostat) (numh(io), io=1,no_u)         ! numhg
  if (iostat /= 0) STOP "numh(io)"
  do io=1,no_u
     if (debug) print *, "numhg ", io, numh(io)
  enddo

  numx = maxval(numh(:))
  allocate(ibuff(numx), hbuff(numx), buff3(3,numx))


  ! Create listhptr 
  listhptr(1)=0
  do io=2,no_u
     listhptr(io)=listhptr(io-1)+numh(io-1)
  enddo
  if (listhptr(no_u)+numh(no_u).gt.nh) STOP "nh overflow in HS"

  do io=1,no_u
     read(hs_u,iostat=iostat) (ibuff(im), im=1,numh(io))
     if (iostat /= 0) STOP "listh"
     do im=1,numh(io)
        listh(listhptr(io)+im) = ibuff(im)
        if (debug) print *, "listh ", io, im, listh(listhptr(io)+im)
     enddo
  enddo

  do is=1,nspin
     do io=1,no_u
        read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,numh(io)
           hamilt(listhptr(io)+im,is) = hbuff(im)
           if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
     enddo
  enddo
  !
  !       Read overlap matrix
  !
  do io=1,no_u
     read(hs_u,iostat=iostat) (hbuff(im), im=1,numh(io))
     if (iostat /= 0) STOP "Overlap matrix read error"
     do im=1,numh(io)
        Sover(listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
     enddo
  enddo

  read(hs_u,iostat=iostat) qtot, temp_in_file 
  if (debug) print *, "QTOT, Temp in file: ", qtot, temp_in_file
  if (iostat /= 0) then
     if (debug) print *, "iostat:", iostat
     STOP "qtot, temp in file"
  endif

  !
  !        Always read xijk
  !
  do io=1,no_u
     read(hs_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,numh(io))
     if (iostat /= 0) STOP "xij(k)"
     do im=1,numh(io)
        ind = listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        xij(1:3,ind) = buff3(1:3,im)
        dij(ind) = sqrt(dot_product(buff3(:,im),buff3(:,im))) / Ang
     enddo
  enddo
  !
  !        Read auxiliary info
  !
  read(hs_u) nspecies
  allocate(label(nspecies), zval(nspecies), no(nspecies))
  read(hs_u) (label(is),zval(is),no(is), is=1,nspecies)
  naoatx = maxval(no(1:nspecies))
  allocate (nquant(nspecies,naoatx), lquant(nspecies,naoatx), &
       zeta(nspecies,naoatx))
  do is=1, nspecies
     do io=1, no(is)
        read(hs_u) nquant(is,io), lquant(is,io), zeta(is,io)
     enddo
  enddo
  read(hs_u) na_u
  allocate(isa(na_u))
  allocate(iaorb(no_u), iphorb(no_u))
  read(hs_u) (isa(ia), ia=1,na_u)
  read(hs_u) (iaorb(io), iphorb(io), io=1,no_u)

  close(hs_u)
  deallocate(ibuff, hbuff, buff3)

end subroutine read_hs_file
end module io_hs
