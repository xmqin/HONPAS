! This file is part of the ONPAS package.
!
!
! Module coded by xmqin, Oct.12, 2013
! Edited by  xmqin, Nov. 27, 2016
!
! Parallel update HFX(u,v) to Hamiltonian, 
! u : Local PAOs of unitcell per node,
! v : Global PAOs of supercell    
!
! HFX(u,v) = (um|vn) * Dm (m,n)
! HFX(u,v) 2D structure ---> Hmat(ind) 1D sparse structure
!
! Hmat(pure DFT)<--- HFX (this subroutine)----> HFX (HSE06)
!
!
! Possible errors from sparse Hamiltonian and density matrix,!!!!!!!!!!!!!
! In siesta Method, density matrix has been stored using the same sparse
! pattern of Hamiltonian. This treatment is robust in pure DFT where the effective
! pure-DFT potential is only determined by the electron density and its gradients.
!
! rho(r) = sum_uv [Duv*psi_u(r-Ru)*psi_v(r-Rv)],  Veff[rho(r)] /= 0 only when
! psi_u(r-R_u) overlaps with psi_v(r-R_v).
! To calculate rho and Veff,  we only need to store a sparse subset of Duv that
! u overlaps with v even though the actual Duv is unknown before SCF. 

! However, HFX is dependent with the density matrix of real space:
! rho(r,r') = sum_uv [Duv*psi_u(r-Ru)*psi_v(r'-Rv)] 
! rho(r,r') ~ exp[-a(r-r'-Ru-Rv)] , a prop to \sqrt(Eg) of semiconductors
! 

      subroutine update_hfx( nspin, norb, iaorb, iphorb, nuo, nuotot, nua,  &
                             na, isa, xa, indxua, cell, nsc, maxnd, numd,   &
                             listdptr, listd, Dscf, maxnh, numh, listhptr,  &
                             listh, samexa, Hmat, Exc)

      use precision
      use parallel,         only : Node, Nodes
      use parallelsubs,     only : GetNodeOrbs, GlobalToLocalOrb, &
                                   LocalToGlobalOrb, WhichNodeOrb
#ifdef MPI
      use mpi_siesta
#endif
!      use atm_types,        only : species, species_info
      use alloc,            only : re_alloc, de_alloc
      use hfx_potential,    only : build_hfx_potential
      use hfx_types,  only: D2Sindx
      use get_DM
  
      implicit none

! ----------------------- Input and Output -------------------------
      integer, intent(in) :: &
        maxnd, maxnh, nua, na, norb, nspin, nuo, nuotot,          &
        iaorb(norb), indxua(na), iphorb(norb), isa(na),           &
        listd(*), listh(maxnh), listdptr(nuo), listhptr(nuo),         &
        numd(nuo), numh(nuo),  nsc(3)

      logical,  intent(in) :: samexa  ! For structure opt

      real(dp), intent(in) ::  cell(3,3), xa(3,na)
      real(dp), intent(in) ::  Dscf(maxnd,nspin)  !  Sparse DM for nuo
      
      real(dp), intent(inout) :: Hmat(maxnh,nspin) ! Sparse Ham for nuo
      real(dp), intent(inout) :: Exc        ! Total Exc for nuotot

! -------------------TEMPOS, INTERNAL VARIABLES--------------------
      integer &
         i, ia, io, iio, is, ispin, j, n
!
! Global buffers for the storage of sparse matrix !
#ifdef MPI
      integer &
        BNode, nuog, indg, maxnumh, MPIerror
      integer, save :: maxnhg, maxndg
      integer, dimension(:), allocatable, save      ::   &
        numhg, listhptrg, listhg
      integer, dimension(:), allocatable, save      ::   &
        numdg, listdptrg, listdg
      real(dp), dimension(:,:), allocatable, save   ::  Dscfg, Hfxg, Hfxg_tmp

!      real(dp), dimension(:,:,:), allocatable, save ::  Hg_EXX
#endif

      integer num_u,num_v,num_m,num_n
      integer iu,ju,jo,ind

      real(dp) E_HFX, time_start, time_end, time_start2, time_end2
!      real(dp), dimension(:,:,:), allocatable, save  ::  &
!                                                      H_EXX, DM_tmp
!      real(dp), dimension(:,:), allocatable, save  :: P_max
      real(dp), dimension(:), allocatable, save  :: Dscf_max

      real(dp) spin_factor

      external memory
! ----------------------------------------------------------------------

!      call cpu_time(time_start)
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

!      call MPI_Bcast(numhg(1:nuotot),nuotot,MPI_integer,BNode, &
!                       MPI_Comm_World,MPIerror)

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

      allocate(Hfxg(maxnhg,nspin))
      call memory('A','D',maxnhg*nspin,'update_hfx')
      Hfxg(1:maxnhg,1:nspin)=0.0_dp


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

      allocate(Dscf_max(maxndg))
!      call memory('A','D',norb*norb,'update_hfx')
      Dscf_max(1:maxndg)=0.0_dp
!      write(6,*) "node, maxn", node, maxndg, maxnhg, maxnh, maxnd

!------------build HFX -----------------

#ifdef MPI
!---------------tranf sparse Dm to full matrix--------------------------
!      call sparse_dense( nspin, nuotot, nuotot, norb, maxndg, numdg, &
!                         listdptrg, listdg, Dscfg, DM_tmp )
!      nullify(D2Sindx)
!      call re_alloc(subshell,1,norb)
    
!      if(node.eq.0) write(6,*) maxnd, numdg
      !nullify(D2Sindx)    
      call re_alloc( D2Sindx, 1, norb, 1, norb,  &
               name='D2Sindx', routine='hfx_update' )
      D2Sindx(:,:) = 0

      call build_indx(nuotot, norb, maxndg, numdg, listdptrg, listdg, D2Sindx)

      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuotot, na, &
                           isa, maxndg, numdg, listdptrg, listdg, &
                           Dscfg, Dscf_max )

!      if (node.eq.0) then
!      write(6,'(a,f16.6,a)')'sparse to dense',time_end2-time_end1,'secs'
!      endif

!      call cpu_time(time_start2)

      call build_hfx_potential( nspin, norb, iaorb, iphorb, nuotot, nuotot, na, isa, xa, &
                                indxua, cell, nsc, maxnhg, numhg, listhptrg, listhg,  &
                                samexa, Dscfg, Dscf_max, Hfxg )
!
! part H_EXX(norb,norb,nspin) per node : parallel over list_mn and list_uv !
! u and m are local in this node, to calc H_EXX(u,m,nspin),otherwise it
! will be zero ! 
!
!      if (allocated(Hg_EXX)) then
!         deallocate(Hg_EXX)
!      endif

      allocate(Hfxg_tmp(maxnhg,nspin))
      Hfxg_tmp(1:maxnhg,1:nspin)=0.0d0

! Get all H_EXX
!      call cpu_time(time_start2)
      call MPI_AllReduce( Hfxg(1,1), Hfxg_tmp(1,1),                &
                          maxnhg*nspin, MPI_double_precision,    &
                          MPI_Sum, MPI_Comm_World, MPIerror )
      Hfxg=Hfxg_tmp
      deallocate(Hfxg_tmp)

!      call cpu_time(time_end2)
!      if(Node.eq.0) WRITE(6,*) " Allreduce time :: " , time_end2 - time_start2, " s "
!#else

!      call sparse_dense( nspin, nuotot, nuotot, norb, maxnd, numd, &
!                         listdptr, listd, Dscf, DM_tmp )

!      call get_pmax_shell( nspin, norb, iaorb, iphorb, nuotot, na, &
!                           isa, maxnd, numd, listdptr, listd, &
!                           Dscf, Dscf_max )
     
!      call build_hfx_potential( nspin, norb, iaorb, iphorb, nuotot, nuotot, na, isa, xa,  &
!                          indxua, cell, nsc, maxnh, numh, listhptr, listh,       &
!                          samexa, DM_tmp, P_max, H_EXX )
#endif
  
!      call cpu_time(time_end)
     
!      if(Node.eq.0) WRITE(6,*) " All HFX time :: " , time_end - time_start, " s "

      if(nspin.eq.1) then
         spin_factor = 1.0d0
      else
         spin_factor = 2.0d0
      endif



! Here H_EXX and DM is global
! Calculate EXX from dense HFX and DM.
      E_HFX = 0.0_dp
      do ispin=1,nspin
         do io=1,nuotot
           do jo=1,numhg(io)
              ind = listhptrg(io)+jo
              E_HFX=E_HFX-0.25d0*spin_factor*Hfxg(ind,ispin) &
                    *Dscfg(ind,ispin)*0.25d0
           enddo
         enddo
      enddo
       
      Exc = Exc + E_HFX
       
      if(node.eq.0)then
         write(6,'(a12,1x,f18.6,1x,a3)') 'HFX energy:',E_HFX*13.60580_dp, "eV"
      endif

! Global H_EXX to local Hmat, io: global orbital, iio: local in this
! node !
#ifdef MPI 
      do ispin=1,nspin
        do iio=1,nuo
           call LocalToGlobalOrb(iio,Node,Nodes,io)
           do j = 1,numh(iio)
              ind = listhptr(iio) + j
              indg = listhptrg(io) + j
              jo = listh(ind)
              Hmat(ind,ispin)= -Hfxg(indg,ispin)*spin_factor*0.5d0*0.25d0 &
                               + Hmat(ind,ispin)
          enddo
        enddo
      enddo

!     WRITE(6,*)">>>>>>>>>>>>>>>>>"
!     write(6,*) Hmat
!#else
!      do ispin=1,nspin
!        do io=1,nuo
!           do j = 1,numh(io)
!              ind = listhptr(io) + j
!              jo = listh(ind)
!              !you need to a combine code  for ns=1, and ns=2
!              Hmat(ind,ispin)=-HFX(ind,ispin)*spin_factor*0.5d0*0.25d0 &
!                              + Hmat(ind,ispin)
!          enddo
!        enddo
!      enddo
#endif
!     call cpu_time(time_end)
     
      deallocate(Hfxg)
      deallocate(Dscfg)
      deallocate(Dscf_max)
      deallocate(listhg)
      deallocate(listdg)
      deallocate(listhptrg)
      deallocate(listdptrg)
      deallocate(numhg)
      deallocate(numdg)
      deallocate(D2Sindx)


      end subroutine update_hfx
