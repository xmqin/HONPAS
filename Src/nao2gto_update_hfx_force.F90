! This file is part of the ONPAS package.
!
!
! Module coded by xmqin, Oct.12, 2013
! Edited by  xmqin, Nov. 27, 2016
!

      subroutine update_hfx_force( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,  &
                                   na, isa, xa, indxua, cell, nsc, maxnd, numd,   &
                                   listdptr, listd, Dscf, maxnh, numh, listhptr,  &
                                   listh, Fal )

      use precision
      use parallel,         only : Node, Nodes
      use parallelsubs,     only : GetNodeOrbs, GlobalToLocalOrb, &
                                   LocalToGlobalOrb, WhichNodeOrb
#ifdef MPI
      use mpi_siesta
#endif
      use alloc,            only : re_alloc, de_alloc
      use hfx_gradient,     only : build_hfx_gradient
      use hfx_types,  only: D2Sindx
      use get_DM
  
      implicit none

! ----------------------- Input and Output -------------------------
      integer, intent(in) :: &
        maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,          &
        iaorb(norb), indxua(na), iphorb(norb), isa(na),           &
        listd(*), listh(maxnh), listdptr(nuo), listhptr(nuo),         &
        numd(nuo), numh(nuo),  nsc(3)

      real(dp), intent(in) ::  cell(3,3), xa(3,na)
      real(dp), intent(in) ::  Dscf(maxnd,nspin)  !  Sparse DM for nuo
      real(dp), intent(inout) :: Fal(3,nua)      

! -------------------TEMPOS, INTERNAL VARIABLES--------------------
      integer &
         i, ia, io, iio, is, ispin, j, n
!
! Global buffers for the storage of sparse matrix !
#ifdef MPI
      integer &
        BNode, nuog, maxnumh, MPIerror, indg
      integer, save :: maxnhg, maxndg
      integer, dimension(:), allocatable, save      ::   &
        numhg, listhptrg, listhg
      integer, dimension(:), allocatable, save      ::   &
        numdg, listdptrg, listdg
      real(dp), dimension(:,:), allocatable, save   ::  Dscfg
#endif

      integer num_u,num_v,num_m,num_n
      integer iu,ju,jo,ind
      real(dp) time_start, time_end
!      real(dp), dimension(:,:,:), allocatable, save  ::  DM_tmp
!      real(dp), dimension(:,:), allocatable, save  :: P_max
      real(dp), dimension(:), allocatable, save  :: Dscf_max
      external memory
! ----------------------------------------------------------------------

     call cpu_time(time_start)
#ifdef MPI
     
      allocate(numhg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(listhptrg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(numdg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
      allocate(listdptrg(nuotot))
      call memory('A','I',nuotot,'update_hfx')
     
!      call cpu_time(time_start) 
! Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
          numdg(io) = numd(iio)
        endif
        call MPI_Bcast(numhg(io),1,MPI_integer,BNode, &
                       MPI_Comm_World,MPIerror)
        call MPI_Bcast(numdg(io),1,MPI_integer,BNode, &
                       MPI_Comm_World,MPIerror)
      enddo

! Build global listhptr
      listhptrg(1) = 0
      listdptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
        listdptrg(io) = listdptrg(io-1) + numdg(io-1)
      enddo

! Globalse listh
      maxnhg = listhptrg(nuotot) + numhg(nuotot)
      maxndg = listdptrg(nuotot) + numdg(nuotot)
      allocate(listhg(maxnhg))
      allocate(listdg(maxndg))
      call memory('A','I',maxnhg,'update_hfx')
      call memory('A','I',maxndg,'update_hfx')

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = &
                 listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
            listdg(listdptrg(io)+1:listdptrg(io)+numdg(io)) = &
                 listd(listdptr(iio)+1:listdptr(iio)+numd(iio))
          enddo
        endif

        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,&
                       BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(listdg(listdptrg(io)+1),numdg(io),MPI_integer,&
                       BNode,MPI_Comm_World,MPIerror)
      enddo
! We trans the sparse matrix to full matrix. What's its nuo type?

      allocate(Dscfg(maxndg,nspin))
      call memory('A','D',maxndg*nspin,'update_hfx')
      Dscfg(1:maxndg,1:nspin)=0.0_dp

!   Globalise Dscf
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ispin = 1,nspin
            do jo = 1,numd(iio)
           Dscfg(listdptrg(io)+jo,ispin) = Dscf(listdptr(iio)+jo,ispin)
            enddo
          enddo
        endif
        do ispin = 1,nspin
          call MPI_Bcast(Dscfg(listdptrg(io)+1,ispin),numdg(io), &
                         MPI_double_precision,BNode, &
                         MPI_Comm_World,MPIerror)
        enddo
      enddo
!      call cpu_time(time_end1)
!      if (node.eq.0) then
!      write(6,'(a,f16.6,a)') "global index",time_end1-time_start, "secs"
!      endif

#endif

!--------------------end for global Dscfg and Hmatg--------------------

#ifdef MPI
!---------------tranf sparse Dm to full matrix--------------------------
!      call sparse_dense( nspin, nuotot, nuotot, norb, maxndg, numdg, &
!                         listdptrg, listdg, Dscfg, DM_tmp )
      allocate(Dscf_max(maxndg))
      Dscf_max(1:maxndg)= 0.0_dp

!      write(61,*) Dscfg
      allocate(D2Sindx(norb,norb))
      D2Sindx(:,:) = 0

      call build_indx(nuotot, norb, maxndg, numdg, listdptrg, listdg, D2Sindx)

      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuotot, na, &
                           isa, maxndg, numdg, listdptrg, listdg, &
                           Dscfg, Dscf_max )
!      call cpu_time(time_end2)
!      if (node.eq.0) then
!      write(6,'(a,f16.6,a)')'sparse to dense',time_end2-time_end1,'secs'
!      endif

      call build_hfx_gradient( nspin, norb, iaorb, iphorb, nuotot, nuotot, nua, na, isa, xa, &
                               indxua, cell, nsc, maxnhg, numhg, listhptrg, listhg,  &
                               Dscfg, Dscf_max, Fal )
      

!      call sparse_dense( nspin, nuotot, nuotot, norb, maxnd, numd, &
!                         listdptr, listd, Dscf, DM_tmp )

!      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuotot, na, &
!                           isa, maxnd, numd, listdptr, listd, &
!                           DM_tmp, P_max )
      
!      call build_hfx_gradient( nspin, norb, iaorb, iphorb, nuotot, nuotot, nua, na, isa, xa,  &
!                          indxua, cell, nsc, maxnh, numh, listhptr, listh,       &
!                          samexa, Dscfg, Dscf_max, Fal )
#endif

      call cpu_time(time_end)
      if(Node.eq.0) write(6,*) " All Force time :: ",  time_end-time_start, " s"

      deallocate(Dscfg)
      deallocate(Dscf_max)

#ifdef MPI
      deallocate(listhg)
      deallocate(listdg)
      deallocate(listhptrg)
      deallocate(listdptrg)
      deallocate(numhg)
      deallocate(numdg)
      deallocate(D2Sindx)
#endif

      end subroutine update_hfx_force
