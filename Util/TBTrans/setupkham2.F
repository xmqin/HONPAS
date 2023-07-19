C ##################################################################
C ##               Setup Hamiltonian in k-space                   ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Modified by F.D.Novaes, fdnovaes@icmab.es                    ##
C ##################################################################
c

        subroutine setupkham2(joutfile,tinit,kpoint,ispin,hsfile,
     &  lasto,Hk,Sk,nua,nuo,nspin)

      use tsread2
C FDN 
      use m_tbt_kpts, only : reclat
C FDN
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------

c=======================================================================
      real*8 EPS
      parameter(EPS=1.0d-8)
      logical PRINTALOT
      parameter(PRINTALOT=.false.)

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
c      real*8 efs(nspin)         ! Fermi energy for each spin

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

      complex*16, dimension (:), pointer:: Hk, Sk
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

      integer inn,it,in,ind
      logical tinit,tlast

C FDN
      integer kscell(3,3)
      real*8 kdispl(3)
C FDN
c-----------------------------------------------------------------------
c SAVED:
      save numh,listh,indxuo,Gamma
      save H,S,efs,rcell,xij,listhptr,nuotot

c=======================================================================
c BEGIN

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
C FDN
      call TSiohs( 'read',hsfile, gamma, nua, nuotot,notot,nspin,
     .              maxnh,numh, listhptr, listh, H, S, qtot, temp,
     .              xij, indxuo, efs, cell, isa, lasto, xa,
     .              kscell, kdispl  )
C FDN
           nuo = nuotot
         


         allocate(Hk(nuo*nuo))
         allocate(Sk(nuo*nuo))

         if(PRINTALOT) then
            write(joutfile,*) 'unit cell:'
            do j=1,3
               write(joutfile,'(3F8.4)') (cell(i,j),i=1,3)
            end do
         end if
         
         call reclat(cell,rcell,1) !reciprocal of cell incl. 2Pi!

         call reclat(cell,recell,0)

         if(.not. Gamma) then

          allocate(xo(3,nuo))
      
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
            enddo
            enddo

c             call memory('D','D',size(xo),'setuphk')
             deallocate( xo )
         end if                 ! not Gamma         

c         call memory('D','I',size(isa),'setuphk')
c         call memory('D','D',size(xa),'setuphk')
         deallocate( isa )
         deallocate( xa )
            

         tinit = .false.
         return
c-----------------------------------------------------------------
         end if                    !tinit
c-----------------------------------------------------------------






       k=kpoint


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
         Sk(in)=Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
         Hk(in)=Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in))
     .           - efs(ispin)*Sk(in)
      enddo



c ===============================================================
       return
       END
c ===============================================================











