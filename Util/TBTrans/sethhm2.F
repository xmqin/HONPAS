C ##################################################################
C ##               Setup Hamiltonian in k-space                   ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################
c

C FDN cell,kscell and kdispl added as dummy
C      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,
C     &     hsfile,nua,lasto,nuo,nspin,Hk,Sk,Hk2,Sk2)    ! ->
      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,
     &     hsfile,nua,lasto,nuo,nspin,cell,kscell,kdispl,
     &     Hk,Sk,Hk2,Sk2)    ! ->
   

      use tsread2
C FDN
      use m_tbt_gf
      use parallel, only : IOnode
C FDN
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------

c=======================================================================
      real*8 EPS
      parameter(EPS=1.0d-8)
c=======================================================================


c INPUT

      integer joutfile          !info-out file unit
      logical tkham             ! true if full H_k hamiltonian is generated

c      integer nuo               !No. states in unitcell (expected in read-in)
      real*8 kpoint(3) 
      integer ispin             !set up H for ispin
      character*33 hsfile       !H,S parameter file


c-----------------------------------------------------------------------
c READ-IN from HS-file/OUTPUT


      integer nua               ! No. atoms in unitcell
      integer mo             ! Number of orbitals in supercell
      integer nuo            ! Number of basis  orbitals
      integer mno             ! Number of orbitals interacting
      integer nspin          ! Spin polarization (1 or 2)

c      integer numh(maxuo)        ! Number of nonzero elements of each row
C                               ! of hamiltonian matrix
c      integer listh(maxno,maxuo) ! Nonzero hamiltonian-matrix element column
C                               ! indexes for each matrix row
c      integer indxuo(maxo)      ! Index of equivalent orbital in unit cell
C                               !  Unit cell orbitals must be the first in
C                               ! orbital lists, i.e. indxuo.le.nuo, with
C                               ! nuo the number of orbitals in unit cell  

c      real*8  H(maxno,maxuo,nspin) ! Hamiltonian in sparse form
c      real*8  S(maxno,maxuo)     ! Overlap in sparse form
      real*8  qtot              ! Total number of electrons
      real*8  temp              ! Electronic temperature for Fermi smearing

      logical Gamma             ! true if Gamma

c      real*8  xij(3,maxno,maxuo) ! Vectors between orbital centers (sparse)
c                               (not read/written if only gamma point ..
c                                but must be written here!!)

      real*8 cell(3,3)          ! unit cell
C FDN
      integer kscell(3,3)
      real*8 kdispl(3)
      
C FDN

c      integer nua               ! No. atoms in unitcell
c      integer isa(maxua)         ! Species index of each atom
c      real*8 xa(3,maxua)         ! Atomic coordinates (Bohr)

    

      integer, dimension (:), pointer:: listh, listhptr,
     .                        numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs


c-----------------------------------------------------------------------
c OUTPUT
c      complex*16 Hk(nuo*nuo), Sk(nuo*nuo)
c      complex*16 Hk2(nuo*nuo), Sk2(nuo*nuo)

      complex*16, dimension (:), pointer:: Hk, Sk, Hk2, Sk2
c      integer, dimension (:), pointer:: lasto
c      integer lasto(0:maxua)   ! Index of last orbital of each atom
c-----------------------------------------------------------------------
c Helpers
      integer nb           ! No. basis orbitals on atoms (the same on all..)
c      real*8 xo(3,maxuo)        ! Atomic coordinates (Bohr)
      real*8, allocatable ::   xo(:,:) 

      integer nuotot,notot,maxnh
      integer ia
      integer i,j,io,jo,iuo,juo,j2
      real*8 k(3),kxij,rcell(3,3)
      real*8 recell(3,3)
      complex*16 cphase

c      integer ix(maxnh)
      integer, dimension (:), pointer :: ix
      integer icoi,icoa
      integer iprop,inn,it,in,ind
      real*8   xc      
      logical tinit,tlast
c-----------------------------------------------------------------------
c SAVED:
C      save numh,listh,indxuo,Gamma
C      save H,S,efs,rcell,xij,listhptr,nuotot
C      save ix

c=======================================================================
c BEGIN
      
       iprop = 3
       icoa = 0
       icoi = 0

       if (.not. LFrstTime .and. LJob) then

C        nullify(H)
C        nullify(S)
C        nullify(xij)
C        nullify(indxuo)
C        nullify(listh)
C        nullify(listhptr)
C        nullify(numh)
C        nullify(efs)
C        nullify(ix)
 
        H=>LH
        if (.not. Lgamma) xij=>Lxij
        S=>LS
        efs=>Lefs
        listh=>Llisth
        listhptr=>Llisthptr
        numh=>Lnumh
        indxuo=>Lindxuo
        allocate(ix(size(Lix)))
        ix=>Lix
        nuotot=Lnuotot
        rcell=Lrcell
        gamma=Lgamma


       end if
       if (.not. RFrstTime .and. RJob) then

C        nullify(H)
C        nullify(S)
C        nullify(xij)
C        nullify(indxuo)
C        nullify(listh)
C        nullify(listhptr)
C        nullify(numh)
C        nullify(efs)
C        nullify(ix)
 
        H=>RH
        if (.not. Rgamma) xij=>Rxij
        S=>RS
        efs=>Refs
        listh=>Rlisth
        listhptr=>Rlisthptr
        numh=>Rnumh
        indxuo=>Rindxuo
        allocate(ix(size(Rix)))
        ix=>Rix
        nuotot=Rnuotot
        rcell=Rrcell
        gamma=Rgamma

       end if

c-----------------------------------------------------------------
c     Read-in Hamiltonian/Overlap parameters from HS-lattice-file
c     Only read-in firsttime:
c
            
       if(tinit) then
c-----------------------------------------------------------------         
              nullify(H)
              nullify(S)
              nullify(xij)
              nullify(indxuo)
              nullify(listh)
              nullify(listhptr)
              nullify(numh)
              nullify(efs)
              nullify(xa)

         if (LFrstTime .and. LJob) then
              nullify(LH)
              nullify(LS)
              nullify(Lxij)
              nullify(Lindxuo)
              nullify(Llisth)
              nullify(Llisthptr)
              nullify(Lnumh)
              nullify(Lefs)      
         end if

         if (RFrstTime .and. RJob) then
              nullify(RH)
              nullify(RS)
              nullify(Rxij)
              nullify(Rindxuo)
              nullify(Rlisth)
              nullify(Rlisthptr)
              nullify(Rnumh)
              nullify(Refs)      
         end if

C FDN kscell and kdispl added as dummys

      call TSiohs('read',
     . hsfile, gamma, nua, nuotot,notot,nspin,
     .              maxnh,numh, listhptr, listh, H, S, qtot, temp,
     .              xij, indxuo, efs, cell, isa, lasto, xa,
     .              kscell, kdispl)

C FDN

           nuo = nuotot
       
              allocate(Hk(nuo*nuo))
              allocate(Sk(nuo*nuo))
              allocate(Hk2(nuo*nuo))
              allocate(Sk2(nuo*nuo))
       
         if (LFrstTime .and. LJob) then
            allocate(LH00(nuo*nuo))
            allocate(LS00(nuo*nuo))
            allocate(LH01(nuo*nuo))
            allocate(LS01(nuo*nuo)) 
C Alocar matrizes que sao salvas
            call alloc_gf_vars(LH,LS,Lxij,Lindxuo,Llisth,Llisthptr,Lnumh
     . ,Lefs,Lix,notot,nuotot,maxnh,nspin,gamma)
         end if         

         if (RFrstTime .and. RJob) then
            allocate(RH00(nuo*nuo))
            allocate(RS00(nuo*nuo))
            allocate(RH01(nuo*nuo))
            allocate(RS01(nuo*nuo))
C Alocar matrizes que sao salvas 
            call alloc_gf_vars(RH,RS,Rxij,Rindxuo,Rlisth,Rlisthptr,Rnumh
     . ,Refs,Rix,notot,nuotot,maxnh,nspin,gamma)     
         end if
        

         if (IOnode) then
            write(joutfile,*) 'unit cell:'

            do j=1,3
               write(joutfile,'(3F8.4)') (cell(i,j),i=1,3)
            end do
         end if ! IOnode 
         
         call reclat(cell,rcell,1) !reciprocal of cell incl. 2Pi!

         call reclat(cell,recell,0)

         if(.not. Gamma) then

            
          allocate(xo(3,nuo))
          
          if (associated(ix)) nullify(ix)
          
        
          allocate (ix(maxnh))
        
c
c     Transform xij so there is no k-dep. phase within uc.
c

c     ... but first find orbital coordinate

            do ia=1,nua
               do iuo=lasto(ia-1)+1,lasto(ia)
                  xo(1,iuo)=xa(1,ia)
                  xo(2,iuo)=xa(2,ia)
                  xo(3,iuo)=xa(3,ia)
               end do           !iuo
            end do              !ia in uc


           do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              xij(1,ind)=xij(1,ind)-(xo(1,juo)-xo(1,iuo))
              xij(2,ind)=xij(2,ind)-(xo(2,juo)-xo(2,iuo))
              xij(3,ind)=xij(3,ind)-(xo(3,juo)-xo(3,iuo))
              xc =0.0
              do  i = 1,3
                xc = xc + xij(i,ind)*recell(i,iprop)
              end do
              ix(ind) = nint(xc)
              icoa = max(icoa,nint(xc))
              icoi = min(icoi,nint(xc))
            enddo
            enddo

             deallocate( xo )
         end if                 ! not Gamma         

         
         deallocate( isa )
         deallocate( xa )
            
        if (LFrstTime .and. LJob ) then
            call cp_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs,ix,
     . rcell,nuotot, gamma)
C          LFrstTime=.false.
        end if         

        if (RFrstTime .and. RJob ) then
            call cp_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs,ix,
     . rcell,nuotot, gamma)
C         RFrstTime=.false.          
        end if


         tinit = .false.
         return
c-----------------------------------------------------------------
         end if                    !tinit
c-----------------------------------------------------------------
         if(tlast) then
!            deallocate(H)
!            deallocate(S)
!            deallocate(xij)
!            deallocate(indxuo)
!            deallocate(listh)
!            deallocate(numh)
            return
         endif
         
   


      k=kpoint
c-----------------------------------------------------------------
      if (tkham)then
c-----------------------------------------------------------------

c
c Setup H,S for this k-point:
c
      do inn = 1,nuo*nuo
        Hk(inn) = dcmplx(0.d0,0.d0)
        Sk(inn) = dcmplx(0.d0,0.d0)
      enddo

      if(.not.Gamma) then



          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
              kxij = (k(1) * xij(1,ind) +
     .                k(2) * xij(2,ind) +
     .                k(3) * xij(3,ind) )
              cphase = cdexp(dcmplx(0d0,1d0)*kxij)
              inn = iuo+(juo-1)*nuo
              Hk(inn) = Hk(inn)+H(ind,ispin)*cphase
              Sk(inn) = Sk(inn)+S(ind)*cphase
            enddo
          enddo


      else !Gamma!

        do io = 1,nuo
          do j = 1,numh(io)
            ind = listhptr(io) + j
            jo = listh(ind)
            inn = io+(jo-1)*nuo
            Hk(inn) = Hk(inn)+H(ind,ispin)
            Sk(inn) = Sk(inn)+S(ind)
          enddo
        enddo

      end if                    !Gamma or not

c
c     Symmetrize and *make EF the energy-zero*!!!
c
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))

           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) )
     .           - efs(ispin)*Sk(it)
           Hk(in) =  dconjg(Hk(it))

        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in) = Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
         Hk(in) = Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in))
     .           - efs(ispin)*Sk(in)
      enddo


c-----------------------------------------------------------------
      else                      !tkham
c-----------------------------------------------------------------

c Setup Transfer H,S for this k||-point:
c
          if(gamma)then
           write(6,*) 'Transfer matrix not possible with gamma'
           stop 'Transfer matrix not possible with gamma'
          endif


         do inn = 1,nuo*nuo
           Hk(inn) = dcmplx(0.d0,0.d0)
           Sk(inn) = dcmplx(0.d0,0.d0)
           Hk2(inn) = dcmplx(0.d0,0.d0)
           Sk2(inn) = dcmplx(0.d0,0.d0)
         enddo

          do iuo = 1,nuo
            do j = 1,numh(iuo)
              ind = listhptr(iuo) + j
              jo = listh(ind)
              juo = indxuo(jo)
               kxij = (k(1) * xij(1,ind) +
     .                 k(2) * xij(2,ind) +
     .                 k(3) * xij(3,ind) -
     .                 k(iprop) *xij(iprop,ind) )
              cphase = cdexp(dcmplx(0d0,1d0)*kxij)

              inn = iuo+(juo-1)*nuo

              if(ix(ind).eq.0) then
                 Hk(inn) = Hk(inn)+H(ind,ispin)*cphase
                 Sk(inn) = Sk(inn)+S(ind)*cphase    
               else if(ix(ind).eq.1) then
                 Hk2(inn) = Hk2(inn)+H(ind,ispin)*cphase
                 Sk2(inn) = Sk2(inn)+S(ind)*cphase
               endif

            enddo
          enddo

c
c     Symmetrize and *make EF the energy-zero*!!!
c
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))

           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) )
     .           - efs(ispin)*Sk(it)
           Hk(in) =  dconjg(Hk(it))

        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in)=Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
         Hk(in)=Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in))
     .           - efs(ispin)*Sk(in)

      enddo
      

      do iuo = 1,nuo
         do juo = 1,nuo
          in = iuo+(juo-1)*nuo
             Hk2(in)=Hk2(in) - efs(ispin)*Sk2(in)
         enddo
      enddo    
        
c-----------------------------------------------------------------
       endif                      !tkham
c-----------------------------------------------------------------

       
c ===============================================================
       return
       END
c ===============================================================











