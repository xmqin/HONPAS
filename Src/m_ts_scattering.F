      MODULE m_ts_scattering
!
! Routines that are used for the scattering system calculation
!
!=============================================================================
! CONTAINS:
!          1) consHS
!          2) consHSp
!          3) getGF_Tri
!          4) getsfe
!          5) GFGammaGF
!          6) weightDM
!          7) weightDMC

      implicit none

#ifdef MPI 
      public :: consHSp, getsfe, getGF_Tri, GFGammaGF,
     .          weightDM, weightDMC
#else
      public :: consHS, getsfe, getGF_Tri, GFGammaGF, 
     .          weightDM, weightDMC
#endif

      private

      CONTAINS

      subroutine consHS( nuo, no,  maxnh,  
     .                  gamma, numh, listhptr, listh, 
     .                  H, S, xij, indxuo,  kpoint,
     .                  Ha,Sa,ef,nou,NGL,NGR,NBUFL)

C *********************************************************************
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxuo               : Last dimension of xij
C                               Must be at least max(indxuo)
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix locally
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix locally
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  kpoint(3)           : k point vectors
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C

      use precision, only: dp

      implicit          none

      integer           maxnh, no, nuo, nou,nbufl,ngl,ngr

      integer           indxuo(no), listh(maxnh), numh(nuo), 
     .                  listhptr(nuo)
      real(dp)  
     .                  H(maxnh), kpoint(3),  
     .                  S(maxnh),xij(3,maxnh)

      real(dp)  
     .                  ef
      complex(dp)  
     .                  Ha(nou,nou), Sa(nou,nou)
      logical gamma

C  Internal variables .............................................
      integer
     .  io, ind, iuo,
     .  j, jo, juo, io1,jo1
      real(dp)   kxij
      complex(dp) cphase
      complex(dp), dimension (:,:), allocatable :: Hw,Sw

C  ....................


      allocate(Hw(nuo,nuo))
      allocate(Sw(nuo,nuo))

      Hw = (0.d0,0.d0)
      Sw = (0.d0,0.d0)

      Ha = (0.d0,0.d0)
      Sa = (0.d0,0.d0)

      do io = 1,nuo
        do j = 1,numh(io)
          ind = listhptr(io) + j
          jo = listh(ind)
          iuo = indxuo(io)
          juo = indxuo(jo)
          if(.not. Gamma) then
            kxij = (kpoint(1) * xij(1,ind) +
     .              kpoint(2) * xij(2,ind) +
     .              kpoint(3) * xij(3,ind))
            cphase = exp(dcmplx(0d0,1d0)*kxij)
c             write(6,*)kxij,iuo,juo
C FDN Changed index ordering
            Sw(iuo,juo) = Sw(iuo,juo) + S(ind)*cphase
            Hw(iuo,juo) = Hw(iuo,juo) + H(ind)*cphase
          else
            Sw(iuo,juo) =  S(ind)
            Hw(iuo,juo) =  H(ind)
C FDN
          endif
        enddo
      enddo
C ....................


      do io1=1,nou
        io=io1+NBUFL
        do jo1=1,nou
          jo=jo1+NBUFL
          Ha(jo1,io1) = Hw(jo,io)
          Sa(jo1,io1) = Sw(jo,io)
        end do                 
      end do

      deallocate(Hw)
      deallocate(Sw)

c Symmetrize:


      do io = 1,nou
        do jo = 1,io-1
C FDN
          Sa(jo,io) = 0.5d0*(Sa(jo,io) + dconjg(Sa(io,jo)) )
          Ha(jo,io) = 0.5d0*(Ha(jo,io) + dconjg(Ha(io,jo)) )
          Sa(io,jo) =  dconjg(Sa(jo,io))
          Ha(io,jo) =  dconjg(Ha(jo,io))
C FDN
        enddo                 
      enddo                  

c
c     Make sure Left-Right part is zero:
c
      do io = nou-NGR+1,nou
        do jo = 1,NGL
          Sa(io,jo) = dcmplx(0d0,0d0)
          Ha(io,jo) = dcmplx(0d0,0d0)
          Sa(jo,io) = dcmplx(0d0,0d0)
          Ha(jo,io) = dcmplx(0d0,0d0)
        end do
      end do

c
c     Shift Energy-zero to ef (initially from SIESTA) :
c
      call zaxpy(nou*nou,dcmplx(-ef,0.d0),Sa,1,Ha,1)


      return
      end subroutine consHS

!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


c
c     Put dummy subroutine here so that compilers do not complain
c     about empty file in serial compilation.
c
C      subroutine dummy_sub2
C      end
c
#ifdef MPI
      subroutine consHSp(nuo, no,   maxnh,  maxnhg,
     .                  gamma, numh, listhptr, 
     .                  numhg, listhptrg,listhg,
     .                  H, S, xijg, indxuo,  kpoint, nuotot,
     .                  Ha,Sa,
     .                  Node, Nodes,ef,nou,NGL,NGR,NBUFL)
C *********************************************************************
C integer nuo                 : Number of basis orbitals in unit cell
C                               local to this processor
C integer no                  : Number of basis orbitals in supercell
C integer maxnh               : Maximum number of orbitals interacting  
C integer maxnd               : First dimension of listd / DM
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix locally
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               hamiltonian matrix locally
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not used if only gamma point)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  kpoint(3)           : k point vectors
C integer nuotot              : total number of orbitals per unit cell
C                               over all processors
C integer Node                : local processor
C integer Nodes               : number of processors
C
C  Modules
C
C FDN
C      use parallel
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
C FDN
      use mpi_siesta
      use precision, only : dp

      implicit          none

      integer 
     .  MPIerror

      integer           maxnh, maxnhg, no,
     .                  nuo, nuotot, Node, Nodes,
     .                  nou,nbufl,ngl,ngr

      integer         
     .                  numh(nuo), listhptr(nuo)
      integer           indxuo(no), listhg(maxnhg), 
     .                  numhg(nuotot), listhptrg(nuotot)
      real(dp)  
     .                  H(maxnh), kpoint(3),  
     .                  S(maxnh),xijg(3,maxnhg)

      real(dp)  
     .                  ef
      complex(dp)  
     .                  Ha(nou,nou), Sa(nou,nou)
      logical gamma

C  Internal variables .............................................
      integer
     .  BNode, io, iio, ind, iuo,
     .  j, jo, juo,  nuog , io1,jo1
      real(dp)
     .  kxij
      real(dp), dimension(:), allocatable, save ::
     .  Snew,Hnew
      complex(dp) cphase
      complex(dp), dimension (:,:), allocatable :: Hw,Sw

C  ....................

C Globalise list arrays - assumes listh and listd are the same


C Create new distribution of H and S
      nuog = nuotot
      allocate(Snew(maxnhg))
      allocate(Hnew(maxnhg))

      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
            do jo = 1,numh(iio)
              Hnew(listhptrg(io)+jo) = H(listhptr(iio)+jo)
            enddo
          do jo = 1,numh(iio)
            Snew(listhptrg(io)+jo) = S(listhptr(iio)+jo)
          enddo
        endif
          call MPI_Bcast(Hnew(listhptrg(io)+1),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
        call MPI_Bcast(Snew(listhptrg(io)+1),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
      enddo

C Find eigenvalues .........................................

        allocate(Hw(nuotot,nuotot))
        allocate(Sw(nuotot,nuotot))


             Hw = (0.d0,0.d0)
             Sw = (0.d0,0.d0)

             Ha = (0.d0,0.d0)
             Sa = (0.d0,0.d0)

          do io = 1,nuog
            do j = 1,numhg(io)
              ind = listhptrg(io) + j
              jo = listhg(ind)
              iuo = indxuo(io)
              juo = indxuo(jo)
               if(.not. Gamma) then
                   kxij = (kpoint(1) * xijg(1,ind) +
     .             kpoint(2) * xijg(2,ind) +
     .             kpoint(3) * xijg(3,ind))
                  cphase = exp(dcmplx(0d0,1d0)*kxij)
C FDN Changed index ordering
                  Sw(iuo,juo) = Sw(iuo,juo) + Snew(ind)*cphase
                  Hw(iuo,juo) = Hw(iuo,juo) + Hnew(ind)*cphase
               else
                  Sw(juo,iuo) =  Snew(ind)
                  Hw(juo,iuo) =  Hnew(ind)
               endif

            enddo
          enddo
C ....................

C Free local memory
      deallocate(Snew)
      deallocate(Hnew)


      do io1=1,nou
         io=io1+NBUFL
         do jo1=1,nou
            jo=jo1+NBUFL
            Ha(jo1,io1) = Hw(jo,io)
            Sa(jo1,io1) = Sw(jo,io)
         end do                 !jo1
      end do

      deallocate(Hw)
      deallocate(Sw)

c Symmetrize:


      do io = 1,nou
         do jo = 1,io-1
C FDN
            Sa(jo,io) = 0.5d0*(Sa(jo,io) + dconjg(Sa(io,jo)) )
            Ha(jo,io) = 0.5d0*(Ha(jo,io) + dconjg(Ha(io,jo)) )
            Sa(io,jo) =  dconjg(Sa(jo,io))
            Ha(io,jo) =  dconjg(Ha(jo,io))
C FDN
         enddo                  !io
      enddo                     !jo


c
c     Make sure Left-Right part is zero:
c
      do io = nou-NGR+1,nou
         do jo = 1,NGL
            Sa(io,jo) = dcmplx(0d0,0d0)
            Ha(io,jo) = dcmplx(0d0,0d0)
            Sa(jo,io) = dcmplx(0d0,0d0)
            Ha(jo,io) = dcmplx(0d0,0d0)
         end do
      end do

c
c     Shift Energy-zero to ef (initially from SIESTA) :
c
      call zaxpy(nou*nou,dcmplx(-ef,0.d0),Sa,1,Ha,1)

      end subroutine consHSp
#endif



!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


C ##################################################################
C ## Calculating Full Greens functions of                         ## 
C ## contact coupled to electrodes using block tri-diagonal       ##
C ##                                                              ##          
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ##  Modified by F. D. Novaes                                    ##
C ##################################################################

      subroutine getGF_Tri(USEBULK,no,NGL,NGR,Hi,Si,
     &     HLL,HRR,ZEnergy,GF,ierr,useTri)


C FDN
      use sys, only: die
      use m_ts_aux_rout, only : splittri, collectri, invmat, transpose,
     .                          csolveg, ts_alloc, ts_dealloc
      use m_ts_options, only: updatedmcr, voltfdf
      use precision, only: dp
C FDN
      implicit none                              


C ======================================================================
c     Note: Indices are starting at 0 in this routine
C ======================================================================
C INPUT:
      logical USEBULK           ! if true self-energy only is input else
c                                 z*S-H-Sigma for bulk is in sfe

      integer no                !no : no. states
      complex(dp) Hi(0:(no*no)-1),Si(0:(no*no)-1)
      complex(dp) ZEnergy 
      integer NGL,NGR
      complex(dp) HLL(NGL,NGL)   ! Left electrode GF
      complex(dp) HRR(NGR,NGR)   ! Right electrode GF
C FDN 
C Use the tridiagonality ?
      logical :: useTri
C FDN 
      

C ======================================================================
C INPUT/OUTPUT:


C ======================================================================
C OUTPUT:
C      complex*16 GF(0:(no*no)-1)  ! 1/(Sc*Energy - H)
      complex(dp), dimension(:), pointer :: GF
      integer ierr              !inversion err
C ======================================================================
c Helpers, tempos ...
      complex(dp), dimension(:), allocatable :: H
      integer, dimension (:), allocatable:: ipvt

      integer i,j,ii,jj
      complex(dp) a,b
C --------------------------------------------------------------------

      integer n1,n2,n3
C FDN 
      integer :: allocStat
      complex(dp), dimension(:,:), pointer ::  a11,a12,
     .                                           a21,a22,a23,
     .                                               a32,a33
      complex(dp), dimension(:,:), pointer ::  b11,b12,
     .                                           b21,b22,b23,
     .                                               b32,b33
      complex(dp), dimension(:,:), pointer ::  x12,x32
C FDN
      complex(dp), dimension(:,:), pointer ::  x21,x23   
      logical :: contactGFonly   

c----------------------------------------------------------------------
      
      external timer, memory

      
    
c=======================================================================
C BEGIN
c=======================================================================

c      call timer1('GFT',1) 
      
c ================================================================
      nullify(a11,a12,a21,a22,a23,a32,a33)
      nullify(b11,b12,b21,b22,b23,b32,b33)
      nullify(x12,x32,x21,x23)
      ierr=0
! FDN
! Compute only GF(2,2) ?
! A difference could be made when updatedmcr(.false.) and volt=0.0 ...
! b11 and b33 are needed only when there is a bias, and because of the
! integration scheme, where both L and R electrodes are used as references.
      contactGFonly = updatedmcr .and. (voltfdf == 0.d0)      

      if(.not.USEBULK) write(6,*) 
     &  "WARNING: Tridiag only works for USEBULK for finite bias!"

c Setup full H (this part will be changed later..)

      allocate(H(no*no),stat=allocStat)
      if(allocStat /=0) call die("Could not allocate in getGF_Tri")
      call memory('A','Z',size(H),"getGF_Tri")
 
c
c     Contact z*S-H without self-energies:
c
      if(useTri) then
         do j = 0, no-1
            do i = 0, no-1
               H(i + 1 + j*no) = 
     .              ZEnergy*Si(i+j*no)-Hi(i+j*no) 
            end do                 !i
         end do                    !j
       else
C GF not yet in memory control !!!
         allocate(GF(no*no),stat=allocStat)
         if(allocStat /=0) call die("Could not allocate in getGF_Tri")
         do j = 0, no-1
            do i = 0, no-1
               H(i + 1 + j*no) = 
     .              ZEnergy*Si(i+j*no)-Hi(i+j*no) 
                GF(i + 1 + no*j) = dcmplx(0d0,0d0)
            end do                 !i
            GF(j + 1 + no*j) = dcmplx(1d0,0d0)
         end do                    !j
       end if 
     
c
c     Add Left Self-energy
c
      do j = 1, NGL
         do i = 1, NGL
            if(USEBULK) then 
               H(i + (j-1)*no) = HLL(i,j)
            else
               H(i + (j-1)*no) = H(i + (j-1)*no) - HLL(i,j)
            end if              !USEBULK
         end do                 !i
      end do                    !j

c
c     Add Right Self-energy
c

      do j = no-NGR+1,no
         jj=j-(no-NGR)
         do i = no-NGR+1,no
            ii=i-(no-NGR)
            if(USEBULK) then 
               H(i + (j-1)*no) = HRR(ii,jj)
            else
               H(i + (j-1)*no) = H(i + (j-1)*no) - HRR(ii,jj)
            end if              !USEBULK
         end do                 !i
      end do                    !j

      if ( .not. useTri) then
! Invert directly
         allocate(ipvt(no),stat=allocStat)
         if(allocStat /=0) call die("Could not allocate in getGF_Tri") 
         call memory('A','I',size(ipvt),"getGF_Tri")      

         call csolveg(no,no,H,GF,ipvt,ierr)
   
         call memory('D','Z',size(H),"getGF_Tri")
         deallocate(H,stat=allocStat)
         if(allocStat /=0) call die("Could not deallocate in getGF_Tri")
 
         call memory('D','I',size(ipvt),"getGF_Tri")
         deallocate(ipvt,stat=allocStat)
         if(allocStat /=0) call die("Could not deallocate in getGF_Tri")

      else

! Use the fact that H is tridiagonal       
c ================================================================
  
c allocate:
      n1 = NGL
      n2 = no - (NGL+NGR)
      n3 = NGR
              
      call ts_alloc(a11,n1,n1,"getGF_Tri")
      call ts_alloc(a12,n1,n2,"getGF_Tri")
      call ts_alloc(a21,n2,n1,"getGF_Tri")      
      call ts_alloc(a22,n2,n2,"getGF_Tri")
      call ts_alloc(a23,n2,n3,"getGF_Tri")      
      call ts_alloc(a32,n3,n2,"getGF_Tri")
      call ts_alloc(a33,n3,n3,"getGF_Tri")
      
      call splittri(n1,n2,n3,H,a11,a12,a21,a22,a23,a32,a33)
      call memory('D','Z',size(H),"getGF_Tri")
      deallocate(H,stat=allocStat)
      if(allocStat /=0) call die("Could not deallocate in getGF_Tri")
      
      

     
c====================================================================
       
      

c     inv(a11) -> store in b11
      call ts_alloc(b11,n1,n1,"getGF_Tri")
      call invmat(n1,a11,b11,ierr)
      call ts_dealloc(a11,"getGF_Tri")
      

c     inv(a33) -> store in b33
      call ts_alloc(b33,n3,n3,"getGF_Tri")
      call invmat(n3,a33,b33,ierr)
      call ts_dealloc(a33,"getGF_Tri")
      

c     x12 = inv(a11).a12 
      call ts_alloc(x12,n1,n2,"getGF_Tri")
      a=dcmplx(1.d0,0.d0)
      b=dcmplx(0.d0,0.d0) 
      call zgemm('N','N', n1, n2, n1 ,a ,b11, n1, a12,
     &     n1, b, x12, n1)     
      call ts_dealloc(a12,"getGF_Tri")
      if(contactGFonly) call ts_dealloc(b11,"getGF_Tri")
      

c     a22 - a21.x12 --> a22
      a=dcmplx(-1.d0,0.d0)
      b=dcmplx(1.d0,0.d0) 
      call zgemm('N','N', n2, n2, n1 ,a ,a21, n2, x12,
     &  n1, b, a22, n2)      
      if(contactGFonly) then
         call ts_dealloc(x12,"getGF_Tri")
         call ts_dealloc(a21,"getGF_Tri")
      end if
      
c     x32 = inv(a33).a32
      call ts_alloc(x32,n3,n2,"getGF_Tri")
      a=dcmplx(1.d0,0.d0)
      b=dcmplx(0.d0,0.d0) 
      call zgemm('N','N', n3, n2, n3 ,a ,b33, n3, a32,
     &     n3, b, x32, n3)   
      call ts_dealloc(a32,"getGF_Tri")
      if(contactGFonly) call ts_dealloc(b33,"getGF_Tri")
       

c     a22 - a23.x32 --> a22
      a=dcmplx(-1.d0,0.d0)
      b=dcmplx(1.d0,0.d0) 
      call zgemm('N','N', n2, n2, n3 ,a ,a23, n2, x32,
     &  n3, b, a22, n2)
      if(contactGFonly) then
         call ts_dealloc(x32,"getGF_Tri")
         call ts_dealloc(a23,"getGF_Tri")
      end if      
      
c now a22 contains (a22 - a21.x12 - a23.x32) 

c b22 = inv(a22 - a21.x12 - a23.x32) 
      call ts_alloc(b22,n2,n2,"getGF_Tri")
      call invmat(n2,a22,b22,ierr)   ! now we've got b22
      call ts_dealloc(a22,"getGF_Tri")

      
      if(.not.contactGFonly) then
         
c  b12 = -x12.b22, b21 = b12^T:
         call ts_alloc(b12,n1,n2,"getGF_Tri")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n1 , n2 , n2 ,a ,x12, n1, b22,
     &     n2, b, b12, n1)      ! got b12
         call ts_dealloc(x12,"getGF_Tri")
         
C FDN
C      call transpose(b12,n1,n2,b21) ! got b21
C b21 = -b22.a21.inv(a11)
C  a21.inv(a11) --> x21
         call ts_alloc(x21,n2,n1,"getGF_Tri")
         a=dcmplx(1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2, n1, n1 ,a ,a21, n2, b11, 
     &     n1, b, x21, n2)
         call ts_dealloc(a21,"getGF_Tri")
         
! b21 = -b22.x21
         call ts_alloc(b21,n2,n1,"getGF_Tri")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2 , n1 , n2 ,a ,b22, n2, x21, 
     &     n2, b, b21, n2)      ! got b21


! b11 is needed only when there is a bias.
! b11 = inv(a11) - b12.x21
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(1.d0,0.d0) 
         call zgemm('N','N', n1 , n1 , n2 ,a ,b12, n1, x21, 
     &     n2, b, b11, n1)
         call ts_dealloc(x21,"getGF_Tri")
C FDN 
 

c  b32 = -x32.b22, b23 = b32^T:
         call ts_alloc(b32,n3,n2,"getGF_Tri")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n3 , n2 , n2 ,a ,x32, n3, b22,
     &     n2, b, b32, n3)      ! got b32  
         call ts_dealloc(x32,"getGF_Tri")
     
C FDN
C      call transpose(b32,n3,n2,b23) ! got b23
C b23 = -b22.a23.inv(a33)
C   a23.inv(a33) --> x23
         call ts_alloc(x23,n2,n3,"getGF_Tri")
         a=dcmplx(1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2, n3, n3 ,a ,a23, n2, b33, 
     &     n3, b, x23, n2)
         call ts_dealloc(a23,"getGF_Tri")
         
C b23 = -b22.x23
         call ts_alloc(b23,n2,n3,"getGF_Tri")
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(0.d0,0.d0) 
         call zgemm('N','N', n2 , n3 , n2 ,a ,b22, n2, x23, 
     &     n2, b, b23, n2)      ! got b23


! b33 is needed only when there is a bias.
! b33 = inv(a33) - b32.x23
         a=dcmplx(-1.d0,0.d0)
         b=dcmplx(1.d0,0.d0) 
         call zgemm('N','N', n3 , n3 , n2 ,a ,b32, n3, x23, 
     &     n2, b, b33, n3)
         call ts_dealloc(x23,"getGF_Tri")

      end if ! .not.contactGFonly
         

c
c     now collect into GF
c     

 
C GF not yet in memory control !!!
      allocate(GF(no*no),stat=allocStat)
      if(allocStat /=0) call die("Could not allocate in getGF_Tri") 
    
      call collectri(n1,n2,n3,b11,b12,b21,b22,b23,b32,b33,GF)
 
      call ts_dealloc(b22,"getGF_Tri")
      
      if(.not.contactGFonly) then
         call ts_dealloc(b11,"getGF_Tri")
         call ts_dealloc(b12,"getGF_Tri")
         call ts_dealloc(b21,"getGF_Tri")
         call ts_dealloc(b23,"getGF_Tri")
         call ts_dealloc(b32,"getGF_Tri")
         call ts_dealloc(b33,"getGF_Tri")
      end if
      
      end if ! .not. useTri
 
     
C ====================================================================
      RETURN
      END subroutine getGF_Tri
C ====================================================================



!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------



C ##################################################################
C ## Subroutine which read-in GF for the 1x1 surface cell         ##
C ## and returns the GF for the NA1xNA2 surface cell              ##
C ## by summing over the appropriate q-points                     ## 
C ##                                                              ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Changed to F90 by Jose-Luis Mozos, jlm@kanigo.icmab.es       ##
C ##################################################################
C version 5. Dec. 2000

C FDN  nq,q,wq no longer dummys
C     The orihinal interface was like this
C
C      subroutine getsfe(USEBULK,jgfu,HAA,SAA,ZEnergy,nq,q,wq,
C     &     NG2,nua,lasto,NA1,NA2,NG,SFE,GAMMA,Nodes,points2read,
C     &     errorgs)

      subroutine getsfe(USEBULK,jgfu,HAA,SAA,ZEnergy,
     &     NG2,NG,SFE,GAMMA,Nodes,points2read,
     &     errorgs,ikinp)

C FDN

! FDN
!      use parallel
      use m_ts_aux_rout, only : csolve 
      use precision, only : dp
! FDN

#ifdef MPI
      use mpi_siesta, only : MPI_Bcast, MPI_Comm_Rank, MPI_AllReduce, 
     . MPI_Sum, DAT_dcomplex, MPI_Integer, MPI_Status_Size, 
     . MPI_Comm_World
               
#endif



c-----------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------
      real(dp) EPS
      parameter(EPS=1d-7)

      real(dp) Pi
      parameter(Pi=3.14159265358979323846264338327950288419717d0)
      
 
c ==================================================================
c INPUT

      logical frstme
      integer icall
      data      frstme /.true./
      save frstme,icall

#ifdef MPI
      integer  MPIerror, Request, Status(MPI_Status_Size)

#endif


      logical USEBULK           ! if true return self-energy only else
c                                 return z*S-H-Sigma for bulk

      integer jgfu              !fileunit from which we read GF

      
      integer NG                ! size of GF in NA1xNA2 supercell

     
      

      complex(dp) ZEnergy 

      integer nq                !no. q-points <= NA1*NA2
C FDN No longer dummys
C      real*8 q(2,nq),wq(nq)     !q-points and weights
      real(dp) q(2,1),wq(1)     !q-points and weights
C FDN


c ----------- 
c 1x1 surface cell
c
      integer NG2               ! size of used part of GAA in UC.

C FDN nq no longer a dummy
C      complex*16 HAA(NG2,NG2,nq),SAA(NG2,NG2,nq)
      complex(dp) HAA(NG2,NG2,1),SAA(NG2,NG2,1)
C FDN
c -----------
      integer Nodes             ! number of processors
      integer points2read


C FDN
      integer ikpar,ikinp
C FDN

c======================================================================
c OUTPUT: 
      complex(dp) SFE(NG,NG) !Electrode atoms self-energy for USEBULK true
c                           !else z*S - H - Sigma for ideal surface
      complex(dp)  GAMMA(NG,NG)  ! 0.5*(Sigma - Sigma^dagger) 
      logical errorgs

c=====================================================================
c Helpers, tempos etc...
      integer Node              ! local processor
      integer ierr,itmp              ! inversion error      
      complex(dp) ZEi, wgfi
      integer iEni
      integer iq,inode
      integer io
      integer jo
      complex(dp) zenloc

      complex(dp), allocatable :: WRK1(:,:),WRK2(:,:),GAA(:,:,:)
      integer,dimension(:), allocatable :: ipvt
#ifdef MPI
      complex(dp), allocatable :: cbuffer(:,:,:)
#endif
      integer  ntot

      external memory

c=======================================================================
c     BEGIN:
c=======================================================================

C FDN Used just to calculate each SE for each kparr
      nq=1
      q=0.0
      wq=1.0
C FDN

      
      if(frstme) icall=1
      if(.not.frstme) icall=icall+1
      frstme=.false.
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

       errorgs =.false.

       ntot =  NG2*NG2*nq

      if(points2read .gt. Nodes) then
         if(Node.eq.0)
     .   write(6,*) 
     .    'ERROR GetSFE: points2read=',points2read,' > Nodes=',Nodes
          errorgs = .true.
          return
      end if
      if(points2read .le. 0) then
         if(Node.eq.0)
     .   write(6,*) 'ERROR GetSFE: points2read=',points2read,'< 0'  
         errorgs = .true.
         return
      end if

c
c     Init:
c
       allocate(WRK1(NG,NG))
       allocate(WRK2(NG,NG))
       allocate(ipvt(NG))
c       call memory('A','Z',2*size(WRK1),'getsfe')
       allocate(GAA(NG2,NG2,nq))


#ifdef MPI
c In order for the inversion to work when inode > points2read
       GAA=dcmplx(0d0,0d0)
       do iq=1,nq
          do io=1,NG2
             GAA(io,io,iq)=dcmplx(1d0,1d0)
          end do
       end do
#endif



c       call memory('A','Z',size(GAA),'getsfe')


      do io=1,NG
         do jo=1,NG
            WRK1(jo,io) = dcmplx(0d0,0d0)
            WRK2(jo,io) = dcmplx(0d0,0d0)
            SFE(jo,io)  = dcmplx(0d0,0d0)
         end do                 !jo
         SFE(io,io) = dcmplx(1d0,0d0) !to be used in inv.
      end do                    !io

c Loop over nodes. If root send stuff out to relevant nodes, if
c not root receive information from root.
c If there is more nodes than points left to read in file (=points2read)
c then just send out the data of the last point to the remaining nodes.


#ifdef MPI
       allocate(cbuffer(NG2,NG2,nq))
#endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      do inode = 0, Nodes-1
cccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ------ root read -------

         if(Node .eq. 0) then
            if(inode .lt. points2read) 
C FDN Now also reads ikpar
     .           read(jgfu) iEni,ZEi,wgfi,ikpar !still enough points to feed to processors
C FDN
         endif                  !node0
#ifdef MPI
      if(Node.eq.inode) then
         zenloc = zenergy
      endif
        
      call MPI_Bcast(zenloc,1,DAT_dcomplex,inode,
     .     MPI_Comm_World,MPIerror)  
      
      call MPI_Bcast(iEni,1,MPI_integer,0,
     .     MPI_Comm_World,MPIerror)  
#else
           zenloc = zenergy

#endif
           if(Node.eq.0) then
             if(cdabs(ZEi-ZEnloc) .gt. EPS) then
              write(6,*)
     .        'ERROR: getsfe: ZEnergy=',ZEi,' expected:',ZEnloc,
     .        ' in processor ', inode
                errorgs = .true.
                return
             end if
C FDN
             if(.not. (ikinp.eq.ikpar)) then
              write(6,*)
     .        'ERROR: getsfe: ikpar=',ikpar,' expected:',ikinp,
     .        ' in processor ', inode
               errorgs = .true.
               return
             end if
C FDN
           end if



         if(iEni .eq. 1) then
           if(Node .eq. 0) then
             read(jgfu) HAA
             read(jgfu) SAA
           endif
           
                     
#ifdef MPI
               call MPI_Bcast(HAA(1,1,1),NG2*NG2*nq,DAT_dcomplex,
     .              0,MPI_Comm_World,MPIerror)  
               call MPI_Bcast(SAA(1,1,1),NG2*NG2*nq,DAT_dcomplex,
     .              0,MPI_Comm_World,MPIerror)           
#endif
          end if              !iEni==1
       


          if(inode .lt. points2read)  then

#ifdef MPI
               if(inode.eq.0.and.Node.eq.0)then
#endif

                   read(jgfu) GAA            !still enough points to feed to processors

#ifdef MPI
                elseif (Node.eq.0) then

                 read(jgfu) cbuffer
               
                call MPI_ISend(cbuffer(1,1,1),ntot,DAT_dcomplex,
     .             inode,1,MPI_Comm_World,Request,MPIerror) 
                call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.inode) then
                call MPI_IRecv(GAA(1,1,1),ntot,DAT_dcomplex,
     .            0,1,MPI_Comm_World,Request,MPIerror)

                call MPI_Wait(Request,Status,MPIerror)

               endif !inode=0 & Node=0

               if (inode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
               endif
#endif
           end if  ! inode.lt.points2read                    
c     -----------------------

c     -----------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do                    !loop over nodes to distribute GF's etc ..
cccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
       deallocate(cbuffer)
#endif

C FDN : CUrrently not used
Cc
Cc     Expansion of unit-cell: (1 x 1) -> (NA1 x NA2) 
Cc     Layer-wise ordering!
Cc
Cc----------------------------------------------
C      do iq=1,nq                !q-points
Cc----------------------------------------------
C
C         io=0
C         do iau=1,nua
C            do ia2=0,NA2-1
C               do ia1=0,NA1-1
C                  do iuo=1+lasto(iau-1),lasto(iau)  
C                     
C                     io=io+1
C
C                     jo=0
C                     do jau=1,nua
C                        do ja2=0,NA2-1
C                           do ja1=0,NA1-1
C                              do juo=1+lasto(jau-1),lasto(jau)  
C                                 
C                                 jo=jo+1
C
C                     
C          cphase  =  exp(dcmplx(0d0,-2d0*Pi)  *  
C     .      ( (ia1-ja1)*q(1,iq) + (ia2-ja2)*q(2,iq)) )    
C
C          ccphase  =  exp(dcmplx(0d0,2d0*Pi)  *  
C     .      ( (ia1-ja1)*q(1,iq) + (ia2-ja2)*q(2,iq)) )    
C
C
Cc
Cc     Here we employ the time-reversal symmetry: 
Cc     GF(1x1)(q)(i,j) = GF(1x1)(-q)(j,i)  
Cc     Note: wq*0.5 since wq already sums up to 1!
C
C
C          WRK1(jo,io) = WRK1(jo,io) + dcmplx(0.5d0*wq(iq),0d0)*
C     .       (cphase*GAA(juo,iuo,iq)+ccphase*GAA(iuo,juo,iq))
C
C          WRK2(jo,io) = WRK2(jo,io) + dcmplx(0.5d0*wq(iq),0d0)*
C     .       (cphase*(ZEnergy*SAA(juo,iuo,iq)-HAA(juo,iuo,iq))
C     .     + ccphase*(ZEnergy*SAA(iuo,juo,iq)-HAA(iuo,juo,iq)))
C             
Cc debug
Cc          WRK2(jo,io) = WRK2(jo,io) + 
Cc     .         dcmplx(0.5d0*wq(iq),0d0)*
Cc     .         (cphase*SAA(juo,iuo,iq)+ccphase*SAA(iuo,juo,iq))
C
C                              
C                              enddo !iuo
C                           enddo !ia1
C                        enddo    !ia2
C                     enddo    !iau
C                     
C                  enddo !juo
C               enddo !ja1
C            enddo    !ja2
C         enddo    !jau
Cc----------------------------------------------
C      enddo                     !q-points
Cc----------------------------------------------
C FDN

C FDN  
       do io = 1,NG
          do jo = 1,NG
             WRK1(io,jo)=GAA(io,jo,1)
             WRK2(io,jo)=ZEnergy*SAA(io,jo,1)-HAA(io,jo,1)
          enddo
       enddo
C FDN

         
c     Find inverse of GAA --> SFE
      ierr=0
      itmp=0

      


      call csolve(NG,WRK1,SFE,ipvt,ierr)



      if(ierr.ne.0) write(90+Node,*) 
     .     "err hola: ,Node,icall,ierr,zenloc,NG" 
      if(ierr.ne.0) write(90+Node,*) 
     .     "err hola: ",Node,icall,ierr,zenloc,NG 


#ifdef MPI
      if (Nodes.gt.1) 
     . call MPI_AllReduce(ierr,itmp,1,MPI_integer,MPI_sum,
     . MPI_Comm_World,MPIerror)
      
      ierr=itmp      
#endif

      IF(IERR.NE.0) THEN
c         errorgs=.true.
         if(Node.eq.0) then
            write(6,*) 'ERROR: getsfe 1 MATRIX INVERSION FAILED'
            write(6,*) 'ERROR: ierr=',ierr
         endif 
      END IF




c     SFE = Z*S - HAA - Sigma ---> WRK2 - SFE = Sigma -> WRK1:
      do io=1,NG
         do jo=1,NG
            WRK1(jo,io) = WRK2(jo,io)-SFE(jo,io)
         end do                 !jo
      end do                    !io


c
c Gamma 
c
      do io=1,NG
         do jo=1,NG
            GAMMA(jo,io) = 
     .           dcmplx(0d0,0.5d0)*(WRK1(jo,io) - 
     .           dconjg(WRK1(io,jo)))
         end do                 !jo
      end do                    !io
      
      if(.not.USEBULK) then
c      
c    Sigma=(Z*S - H) - inv[GF] ---> SFE
c
         do io=1,NG
            do jo=1,NG
               SFE(jo,io) = WRK1(jo,io)
            end do              !jo
         end do                 !io

      end if                    !USEBULK


c     debug
c      do io=1,NG
c         do jo=1,NG
c            SFE(jo,io) = WRK2(jo,io)
c         end do                 !jo
c      end do                    !io


       call memory('D','Z',2*size(WRK1),'getsfe')
       call memory('D','Z',size(GAA),'getsfe')
       deallocate(WRK1)
       deallocate(WRK2)
       deallocate(GAA)
       deallocate(ipvt)


c======================================================================
      return
      end subroutine getsfe
c======================================================================


!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------



      subroutine GFGammaGF(LEFTFLAG,nc,ng,GF,Gamma,GGG)

C ======================================================================
c  This routine returns GGG=(-i)*GF.Gamma.GF, where GF is a (nc)x(nc)
c  matrix and the first[LEFTFLAG true]/last[LEFTFLAG false] states
c  corresponds to the (ng) Left/Right electrode states.
c  Gamma is a (ng)x(ng) matrix.
C ======================================================================
      use precision, only : dp
      implicit none
C ======================================================================

      logical PRINTALOT
c      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)

c
C ======================================================================
C INPUT:

      
      logical LEFTFLAG          !True if LEFT
      integer nc                !no. states in contact region
      integer ng                !if LeftFLAG the ng first states of nc
c                               !are left electrode states
c                               !else it is the ng last right electrode
c                               !states.
      complex(dp) Gamma(ng,ng)   !i(Sigma-Sigma^dagger)/2
      complex(dp) GF(nc,nc)

C ======================================================================
C OUTPUT:
C      
      complex(dp) GGG(nc,nc)     !GF.GAMMA.GF

C ======================================================================
c Helpers, tempos ...
      real(dp), allocatable ::  swrk1(:,:),swrk2(:,:)
      complex(dp), allocatable ::  swrk3(:,:),swrk4(:,:)
      integer i,j,ie,je

c=======================================================================
C BEGIN
c=======================================================================

      if(nc.lt.ng) then
         write(6,*) 'ERROR: GFGammaGF: nc=',nc,'<ng=',ng
         stop 'ERROR: GFGammaGF: nc<ng'
      end if
      
      call timer('TRANS2',1)
      
      allocate(swrk2(nc,nc))
      allocate(swrk4(nc,nc))
      
      
      
c     do i=1,nc
c     do j=1, nc 
c     GGG(i,j) =  dcmplx(0d0,0.d0)
c          swrk2(i,j) =  0.d0
c     swrk4(i,j) =  (0.d0,0.d0)
c     enddo
c     enddo
      
      swrk2 = 0.0d0
      swrk4 = dcmplx(0d0,0.d0)
      GGG = dcmplx(0d0,0.d0)


      if(LEFTFLAG) then
         
         allocate(swrk1(ng,ng))
         allocate(swrk3(ng,nc))

         do je=1,ng
            do ie=1,ng
               swrk1(ie,je) = DREAL(Gamma(ie,je))
            end do                
         end do                
         
c     do je=1,ng
c     do i=1,nc
c     swrk3(je,i) =(0.0d0,0.0d0)
c     end do                
c     end do                

         swrk3 = dcmplx(0d0,0.d0)

         do i=1,nc
            do j=1, nc 
               swrk4(i,j) = GF(j,i)
            enddo
         enddo
         
         do i=1,nc
            do ie=1, ng
               do je=1,ng
             swrk3(ie,i) =  swrk3(ie,i) +
     .                 swrk1(je,ie) *  swrk4(je,i)
          enddo
       enddo
      enddo
      
 
      do i=1,nc
         do j=1,nc
           do ie=1,ng
              swrk2(j,i) =  swrk2(j,i) +
     .             dreal (swrk3(ie,i) *  dconjg(swrk4(ie,j)))
           enddo
        enddo
      enddo
 
 
      do i=1,nc
         do j=1,i
            GGG(j,i) = swrk2(j,i) *  dcmplx(0d0,-1.d0)
            GGG(i,j) = GGG(j,i)
         end do                 !j
      end do                    !i
      
      else
 
         do i=1,nc
            do j=1, nc
               swrk4(i,j) = GF(j,i)
            enddo
         enddo
 
         allocate(swrk1(nc-ng+1:nc,nc-ng+1:nc))
         allocate(swrk3(nc-ng+1:nc,nc))
 
 
         do je=1,ng
            do ie=1,ng
               swrk1(ie+nc-ng,je+nc-ng) = DREAL(Gamma(ie,je))
            end do                
         end do                
 
         do je=nc-ng+1,nc
            do i=1,nc
               swrk3(je,i) =(0.0d0,0.0d0)
            end do                
         end do                
         
         do i=1,nc
           do j=1, nc
              swrk4(i,j) = GF(j,i)
           enddo
        enddo
        
        do i=1,nc
           do ie= nc-ng+1,nc
              do je= nc-ng+1,nc
                 swrk3(ie,i) =  swrk3(ie,i) +
     .                swrk1(je,ie) *  swrk4(je,i)
              enddo
           enddo
        enddo
 
 
        do i=1,nc
           do j=1,nc
              do ie= nc-ng+1,nc
             swrk2(j,i) =  swrk2(j,i) +
     .                dreal (swrk3(ie,i) *  dconjg(swrk4(ie,j)))
          enddo
       enddo
      enddo
 
 
      do i=1,nc
         do j=1,i
            GGG(j,i) = swrk2(j,i) *  dcmplx(0d0,-1.d0)
            GGG(i,j) = GGG(j,i)
         end do                 !j
      end do                    !i
 
      endif
 

      deallocate(swrk1)
      deallocate(swrk2)
      deallocate(swrk3)
      deallocate(swrk4)

      call timer('TRANS2',2)

C ====================================================================
      RETURN
      END subroutine GFGammaGF
C ====================================================================


!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


C ##################################################################
C ##   Mixing the Density matrixes according to the smallest      ##
C ##    realspace integral                                        ##
C ##                                                              ##
C ##  Version 011200  Kurt Stokbro, ks@mic.dtu.dk                 ##
C ##################################################################


      subroutine weightDM(nc,NGL,NGR,DML,DMR,DMneqL,DMneqR,
     .     EDML,EDMR)
C ======================================================================
c  This routine find weight for the DM integral. On output
C  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
C  EDML:= w EDML +(1-w) EDMR
C  In left electrode w=1 and in right electrode w=0
C ======================================================================

#ifdef MPI
      use mpi_siesta, only: MPI_Comm_Rank, MPI_Comm_World
#endif
      use sys, only : die
      use precision, only : dp

      implicit none

      logical PRINTALOT
      parameter(PRINTALOT=.FALSE.)
c      parameter(PRINTALOT=.TRUE.)

c
C ======================================================================
C INPUT:

      
      integer nc                !no. states in density matrix
      integer NGL               !no. states in left electrode
      integer NGR               !no. states in right electrode
      real(dp) DML(nc,nc),DMR(nc,nc)      !contour part of DM integration
      real(dp) DMneqL(nc,nc),DMneqR(nc,nc)!real axis part of DM integration
      real(dp) EDML(nc,nc),EDMR(nc,nc)    !L,R estimates of Edm

C ======================================================================
c Helpers, tempos ...
      real(dp) w,wl,wr,wsum
      real(dp) maxdif2            !largest integration error
      real(dp) maxndif2           !largest normalized integration error
      
      integer i,j

      integer Node

#ifdef MPI
      integer
     .  MPIerror
#endif


c=======================================================================
c BEGIN:
c=======================================================================


C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif

      if (nc.lt. (NGL+NGR))  then
         if(Node.eq.0)
     .   write(6,*) 'ERROR: weightDM: nc=',nc,'<=',NGL+NGR
         call die ('ERROR: weightDM: nc<,NGL+NGR')
      end if
     
      maxdif2=0.0d0
      maxndif2=0.0d0
      do j=1,nc
         do i=1,nc
            wr=DMneqR(i,j)*DMneqR(i,j)
            wl=DMneqL(i,j)*DMneqL(i,j)
            wsum=wr+wl
            if(wsum .gt. 0d0) then
                w = wl/wsum
             else
                w=0.5d0
             end if

c            write(*,*) i,j,DMneqR(i,j),DMneqL(i,j),DML(i,j),DMR(i,j)
c  In the leads we always use either left or right integration scheme.

            if (i .le. NGL .or. j .le. NGL) then
! FDN With TriDiag, it must be the opposite !! 
! GF(1,3) and GF(3,1) are not zero !, only in TriDiag ...
!               w = 1.0d0
               w = 0.d0
            endif
            if (i .gt. nc-NGR .or. j .gt. nc-NGR) then
! FDN With TriDiag, it must be the opposite !!
!               w = 0.0d0
               w = 1.d0
            endif
c weight DM and DML
            DML(i,j) = w   *(DML(i,j)+DMneqR(i,j))
     .           +(1.0-w)*(DMR(i,j)+DMneqL(i,j))
            EDML(i,j) = w*EDML(i,j)+(1.0-w)*EDMR(i,j)
C make an integration error estimate
! FDN
!            dif=DML(i,j)+DMneqR(i,j)-(DMR(i,j)+DMneqL(i,j))
!            dif2=dif*dif
!            ndif2=dif2*wl*wr/(wsum*wsum)
!            if (dif2 .gt. maxdif2) then
!               maxdif2=dif2
!               idif2=i
!               jdif2=j
!            endif
!            if (ndif2 .gt. maxndif2) then
!               maxndif2=ndif2
!               indif2=i
!               jndif2=j
!            endif
         enddo                  !i
      enddo                     !j

! FDN
!      if(Node.eq.0) then
!
!      if (PRINTALOT) then
!         write(joutfile,*)
!     .    'Integration error estimate',indif2,jndif2, sqrt(maxndif2)
!         write(joutfile,*)
!     .    'unnormalized error ',idif2,jdif2, sqrt(maxdif2)
!      endif
!
!      endif

C ====================================================================
      RETURN
      END subroutine weightDM
C ====================================================================


C FDN Routine for complex matrices
      subroutine weightDMC(nc,NGL,NGR,DML,DMR,DMneqL,DMneqR,
     .     EDML,EDMR)
C ======================================================================
c  This routine find weight for the DM integral. On output
C  DML := w (DML+DMneqR) + (1-w) (DMR+DMneqL)
C  EDML:= w EDML +(1-w) EDMR
C  In left electrode w=1 and in right electrode w=0
C ======================================================================

#ifdef MPI
      use mpi_siesta, only: MPI_Comm_Rank, MPI_Comm_World
#endif
      use sys, only : die
      use precision, only : dp

      implicit none

      logical PRINTALOT
      parameter(PRINTALOT=.FALSE.)
c      parameter(PRINTALOT=.TRUE.)

c
C ======================================================================
C INPUT:

      
      integer nc                !no. states in density matrix
      integer NGL               !no. states in left electrode
      integer NGR               !no. states in right electrode
      complex(dp) DML(nc,nc),DMR(nc,nc)      !contour part of DM integration
      complex(dp) DMneqL(nc,nc),DMneqR(nc,nc)!real axis part of DM integration
      complex(dp) EDML(nc,nc),EDMR(nc,nc)    !L,R estimates of Edm

C ======================================================================
c Helpers, tempos ...
      real(dp) w,wl,wr,wsum
      real(dp) maxdif2            !largest integration error
      real(dp) maxndif2           !largest normalized integration error
      
      integer i,j

      integer Node

#ifdef MPI
      integer
     .  MPIerror
#endif


c=======================================================================
c BEGIN:
c=======================================================================


C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
#else
      Node = 0
#endif


      if (nc.lt. (NGL+NGR))  then
         if(Node.eq.0)
     .   write(6,*) 'ERROR: weightDM: nc=',nc,'<=',NGL+NGR
         call die ('ERROR: weightDM: nc<,NGL+NGR')
      end if
     
      maxdif2=0.0d0
      maxndif2=0.0d0
      do j=1,nc
         do i=1,nc
            wr=cdabs(DMneqR(i,j))**2
            wl=cdabs(DMneqL(i,j))**2
            wsum=wr+wl
            if(wsum .gt. 0d0) then
                w = wl/wsum
             else
                w=0.5d0
             end if

c            write(*,*) i,j,DMneqR(i,j),DMneqL(i,j),DML(i,j),DMR(i,j)
c  In the leads we always use either left or right integration scheme.

            if (i .le. NGL .or. j .le. NGL) then
! FDN With TriDiag, it must be the opposite !!
!               w = 1.0d0
               w = 0.d0
            endif
            if (i .gt. nc-NGR .or. j .gt. nc-NGR) then
! FDN With TriDiag, it must be the opposite !!
!               w = 0.0d0
               w = 1.d0
            endif
c weight DM and DML
            DML(i,j) = w   *(DML(i,j)+DMneqR(i,j))
     .           +(1.0-w)*(DMR(i,j)+DMneqL(i,j))
            EDML(i,j) = w*EDML(i,j)+(1.0-w)*EDMR(i,j)
C make an integration error estimate
! FDN
!            dif=cdabs(DML(i,j)+DMneqR(i,j)-(DMR(i,j)+DMneqL(i,j)))
!            dif2=dif*dif
!            ndif2=dif2*wl*wr/(wsum*wsum)
!            if (dif2 .gt. maxdif2) then
!               maxdif2=dif2
!               idif2=i
!               jdif2=j
!            endif
!            if (ndif2 .gt. maxndif2) then
!               maxndif2=ndif2
!               indif2=i
!               jndif2=j
!            endif
         enddo                  !i
      enddo                     !j

! FDN
!      if(Node.eq.0) then

!      if (PRINTALOT) then
!         write(joutfile,*)
!     .    'Integration error estimate',indif2,jndif2, sqrt(maxndif2)
!         write(joutfile,*)
!     .    'unnormalized error ',idif2,jdif2, sqrt(maxdif2)
!      endif
!
!      endif

C ====================================================================
      RETURN
      END subroutine weightDMC
C ====================================================================









      END MODULE m_ts_scattering
