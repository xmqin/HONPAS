      program ReadTSHSnetcdf

C integer nuotot              : Number of basis orbitals per unit cell
C integer nspin               : Spin polarization (1 or 2)
C integer maxnh               : First dimension of listh, H, S and
C                               second of xij
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row (-1)
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  qtot                : Total number of electrons
C real*8  temp                : Electronic temperature for Fermi smearing
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not read/written if only gamma point)
C integer nuo              : Number of basis orbitals per supercell
C integer nua                 : No. atoms in unitcell
C integer isa(nua)            : Species index of each atom
C integer lasto(0:nua)        : Index of last orbital of each atom
C real*8 xa(3,nua)            : Atomic coordinates (Bohr)
C real*8 ef                   : Fermi energy
C real*8 efs(nspin)           : Fermi energy for each spin
c real*8  cell(3,3)           : Unit cell vectors in real space cell(ixyz,ivec)
C integer indxuo(notot)       : Index of equivalent orbital in unit cell
C                             Unit cell orbitals must be the first in
C                             orbital lists, i.e. indxuo.le.nuo, with
C                             nuo the number of orbitals in unit cell


C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

      use netcdf

      implicit          none



      character         finname*33
      logical           gamma
      integer           maxnh, nuo, nspin

      double precision  qtot, temp


      integer, dimension (:), pointer:: listh, listhptr,
     .                        numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs

c Internal variables and arrays
      integer    im, is, iu, j, ju, k, mnh, ns, i,ic
      integer    ih,hl,nut
      logical    found
cccc
      integer nua                
      double precision ef   
      double precision  cell(3,3),rcell(3,3)       


      integer notot


c================================================================
c NetCDF: Begin
      integer  iret             !error status return
      integer  ncid             ! NetCDF id
* dimension ids
      integer nR_dim, nuo_dim, no_dim 
      integer maxnh_dim, nspin_dim, nua_dim,nuaa_dim
* variable ids
      integer  listh_id,numh_id
      integer  lasto_id
      integer  H_id, S_id, efs_id
      integer  xa_id, xij_id
      integer  indxuo_id,rcell_id
* variable shapes
      integer  listh_dims(1),numh_dims(1)
      integer  lasto_dims(1)
      integer  H_dims(2),S_dims(1),efs_dims(1)
      integer  xa_dims(2),xij_dims(2)
      integer  indxuo_dims(1),rcell_dims(2)

c NetCDF: End
c================================================================
c
c write to mathematica:
      integer iuo,juo,jo,inn,ind,ispin
      real*8 kp(3),kxij
      complex*16 cphase
      complex*16, dimension (:), pointer:: Hk, Sk
c================================================================


cccc
      write(*,*) "input file:"
      read(5,*) finname

      
C Check if input file exists
      inquire( file=finname, exist=found )

      if (found) then
C Open file
         iu=10
         open( iu, file=finname, status='old' ,
     .        form = 'unformatted')
C     Read dimensions
         read(iu) nua, nuo, nspin, notot, maxnh
         write(6,*) 'dim:', nua, nuo, nspin, notot, maxnh
c  allocations 

                        
         allocate(isa(nua))            
         allocate(lasto(0:nua))
         allocate(xa(3,nua))
         allocate(indxuo(notot))
         
         read(iu)(isa(i),i=1,nua)
         read(iu)(lasto(i),i=0,nua)
         read(iu)((xa(j,i),j=1,3),i=1,nua)
         read(iu)(indxuo(i),i=1,notot)
         

C Allocate local array for global numh
         allocate(numh(nuo))

         nut = 0
         do ih = 1,nuo
            read(iu) numh(ih)
            nut = nut + numh(ih)
         enddo
         

         allocate(listhptr(nuo))
         
C     Create listhptr
         listhptr(1) = 0
         do hl = 2,nuo
            listhptr(hl) = listhptr(hl-1) + numh(hl-1)
         enddo
         
          
         allocate(listh(maxnh))

         ic = 0
         do ih = 1,nuo
            hl = ih
            do im = 1,numh(hl)
               ic = ic +1
               read(iu) listh(listhptr(hl)+im)
            enddo
         enddo



C Read Hamiltonian

         allocate(H(maxnh,nspin))
         
         do is = 1,nspin
            do ih = 1,nuo
               hl = ih
               do im = 1,numh(hl)
                  read(iu) H(listhptr(hl)+im,is)
               enddo
            enddo
         enddo



C Read Overlap matrix
         
         allocate(S(maxnh))
         
         do ih = 1,nuo
            hl = ih
            do im = 1,numh(hl)
               read(iu) S(listhptr(hl)+im)
            enddo
         enddo
c         write(*,*) S
         
         allocate(efs(nspin)) 
           
         read(iu) qtot,temp
         if(nspin.gt.1) read(iu) efs
         if(nspin.eq.1) then 
            read(iu) ef
            efs(1)=ef
         end if

         read(iu) cell
         
c reciprocal lattice cell:
         call reclat(cell,rcell,1)
c         write(*,*) cell,rcell

         read(iu) gamma  
         if(gamma) write(*,*) 'This is a gamma file: may not work'

         if (.not.gamma) then
C Read interorbital vectors for K point phasing

            allocate(xij(3,maxnh)) 
           
            
            
c            do ih = 1,maxnh
c               do i=1,3
c                  xij(i,ih)=0d0
c               enddo
c            enddo


            do ih = 1,nuo
               hl = ih
               do im = 1,numh(hl)
                  read(iu) (xij(k,listhptr(hl)+im),k=1,3)
              enddo
            enddo

         endif
c         write(*,*) xij
            
c     Close file
      else
         write(6,*) 'iohs: ERROR: file not found: ', finname
         stop 'iohs: ERROR: file not found'
      endif
      

c debug 
      write(77,*) "isa: ",isa
      write(77,*) "lasto: ",lasto
      write(77,*) "xa: ",xa
      write(77,*) "indxuo: ",indxuo
      write(77,*) "listh: ",listh
      write(77,*) "numh: ",numh
      write(77,*) "qtot: ",qtot
      write(77,*) "efs,ef: ",efs,ef
      write(77,*) "gamma: ",gamma
      write(77,*) "cell: ",cell
      write(77,*) "H:",H
      write(77,*) "S:",S
      write(77,*) "xij: ",xij



ccccccccccccccccccccccccc
c netCDF Begin
ccccccccccccccccccccccccc

      iret = nf90_create('HS.nc', NF90_CLOBBER, ncid)
      call check_err(iret)


*     define dimensions
      iret = nf90_def_dim(ncid, 'nR', 3, nR_dim)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'nuo', nuo, nuo_dim)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'no', notot, no_dim)
      call check_err(iret)

      iret = nf90_def_dim(ncid, 'maxnh', maxnh, maxnh_dim)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'ns', nspin, nspin_dim)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'nua', nua, nua_dim)
      call check_err(iret)
      iret = nf90_def_dim(ncid, 'nuaa', nua+1, nuaa_dim)
      call check_err(iret)


*     define variables
      listh_dims(1)      = maxnh_dim
      numh_dims(1)       = nuo_dim
      lasto_dims(1)      = nuaa_dim

      H_dims(1)          = maxnh_dim
      H_dims(2)          = nspin_dim
      S_dims(1)          = maxnh_dim

      efs_dims(1)        = nspin_dim
      xa_dims(1)         = nR_dim
      xa_dims(2)         = nua_dim

      xij_dims(1)        = nR_dim      
      xij_dims(2)        = maxnh_dim
      
      indxuo_dims(1)     = no_dim  
      rcell_dims(1)      = nR_dim
      rcell_dims(2)      = nR_dim


      iret = nf90_def_var(ncid, 'listh',
     .     NF90_INT,listh_dims,listh_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'numh', 
     .     NF90_INT,numh_dims, numh_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'lasto', 
     .     NF90_INT,lasto_dims, lasto_id)
      call check_err(iret)
         
      iret = nf90_def_var(ncid, 'H', NF90_DOUBLE,H_dims, H_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'S', NF90_DOUBLE, S_dims, S_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'efs', NF90_DOUBLE, efs_dims, efs_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'xa', NF90_DOUBLE, xa_dims, xa_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'rcell',
     .     NF90_DOUBLE,rcell_dims,rcell_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'indxuo', 
     .     NF90_INT, indxuo_dims, indxuo_id)
      call check_err(iret)

      iret = nf90_def_var(ncid, 'xij',NF90_DOUBLE,xij_dims,xij_id)
      call check_err(iret)


*     assign attributes

      iret = nf90_put_att(ncid, listh_id, 'long_name',
     .     'listh')
      call check_err(iret)

      iret = nf90_put_att(ncid, numh_id, 'long_name',
     .     'Number of non-zero elements in row')
      call check_err(iret)

      iret = nf90_put_att(ncid, lasto_id, 'long_name',
     .     'index of last orbital of atom')
      call check_err(iret)

      call check_err(iret)
      iret = nf90_put_att(ncid, H_id, 'long_name',
     .     'Hamiltonian, sparse')
      call check_err(iret)

      iret = nf90_put_att(ncid, S_id, 'long_name',
     .     'Overlap matrix, sparse')
      call check_err(iret)

      iret = nf90_put_att(ncid, efs_id, 'long_name',
     .     'Fermi energy, spin dependent')
      call check_err(iret)

      iret = nf90_put_att(ncid, xa_id, 'long_name',
     .     'Atomic coordinates')
      call check_err(iret)


      iret = nf90_put_att(ncid, rcell_id, 'long_name',
     .     'Reciprocal lattice')
      call check_err(iret)

      iret = nf90_put_att(ncid, indxuo_id, 'long_name',
     .     'Index, sparse')
      call check_err(iret)

      iret = nf90_put_att(ncid, xij_id, 'long_name',
     .     'Vector between atoms, sparse')
      call check_err(iret)


      iret = nf90_put_att(ncid, NF90_GLOBAL, 
     .     'source', 'readTSHScdf')

      call check_err(iret)

*     leave define mode
      iret = nf90_enddef(ncid)
      call check_err(iret)


* write to netCDF:

      iret = nf90_put_var(ncid, listh_id, listh(1:),start=(/1/))
      call check_err(iret)
      
      iret = nf90_put_var(ncid, numh_id, numh(1:),start=(/1/))
      call check_err(iret)

      iret = nf90_put_var(ncid, lasto_id, lasto(0:),start=(/1/))
      call check_err(iret)

c      write(*,*) "lasto:",lasto
c      do i =0,nua
c         write(*,*) i,lasto(i)
c      end do
 
      iret = nf90_put_var(ncid, H_id, H(1:,1:),start=(/1,1/))
      call check_err(iret)

c     write S to netCDF:
      iret = nf90_put_var(ncid, S_id, S(1:),start=(/1/))
      call check_err(iret)

      iret = nf90_put_var(ncid, efs_id, efs(1:),start=(/1/),
     .     count=(/nspin/))
      call check_err(iret)
     
     
      iret = nf90_put_var(ncid, xa_id, xa(1:,1:),start=(/1,1/),
     .     count=(/3,nua/))
      call check_err(iret)

 
      iret = nf90_put_var(ncid, rcell_id, rcell(1:,1:),
     .     start=(/1,1/),count=(/3,3/))
      call check_err(iret)
 
      iret = nf90_put_var(ncid, indxuo_id, indxuo(1:),
     .     start=(/1/))
      call check_err(iret)
      
      iret = nf90_put_var(ncid, xij_id, xij)
      call check_err(iret)


c I dont know why I have to do like this to make it work:
      do k=1,notot
         iret = nf90_put_var(ncid, indxuo_id, indxuo(k),
     .        start=(/ k /))
         call check_err(iret)      
      end do
    
c I dont know why I have to do like this to make it work:
      do i=1,maxnh
         do j=1,nspin
            iret = nf90_put_var(ncid, H_id, H(i,j),start=(/i,j/))
            call check_err(iret)      
         end do
      end do
      
c I dont know why I have to do like this to make it work:
      do i=1,maxnh
         iret = nf90_put_var(ncid, S_id, S(i),start=(/i/))
         call check_err(iret)      
      end do
       
c I dont know why I have to do like this to make it work:
      do i=1,3
         do k=1,maxnh
            iret = nf90_put_var(ncid, xij_id, xij(i,k),start=(/i,k/))
            call check_err(iret)
         end do
      end do
      
 

ccccccccccccccccccccccccc
c netCDF End
ccccccccccccccccccccccccc


      if(.true.) then
CCCCCCCCCCCCCCCCC Write to mathematica for the gamma-point as check:
         kp(1)=0d0
         kp(2)=0d0
         kp(3)=0d0
         ispin=1
         allocate(Hk(nuo*nuo))
         allocate(Sk(nuo*nuo))
         
         do i=1,nuo*nuo
            Sk(i)=dcmplx(0d0,0d0)
            Hk(i)=dcmplx(0d0,0d0)
         end do

         do iuo = 1,nuo
            do j = 1,numh(iuo)
               ind = listhptr(iuo) + j
               jo = listh(ind)
               juo = indxuo(jo)
               kxij = (kp(1) * xij(1,ind) +
     .                 kp(2) * xij(2,ind) +
     .                 kp(3) * xij(3,ind) )
               cphase = cdexp(dcmplx(0d0,1d0)*kxij)
               inn = iuo+(juo-1)*nuo
               Hk(inn) = Hk(inn) + H(ind,ispin)*cphase
               Sk(inn) = Sk(inn) + S(ind)*cphase
               write(99,*) ind,S(ind),iuo,juo
            enddo
         enddo

         call writezmat(nuo,"HAM0",Hk)
         call writezmat(nuo,"OVL0",Sk)

         deallocate(Hk)
         deallocate(Sk)
      end if






c=========================

      CONTAINS

c=========================


      subroutine check_err(status)
      integer, intent(in):: status
      if (status .ne. nf90_noerr) then
         print  *, trim(nf90_strerror(status))
         stop 'netcdf Error: Stopped'
      endif
      end subroutine check_err



      subroutine reclat(A,B,IOPT)

C  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
C  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(3,3),B(3,3)
      integer iopt
      PI=ACOS(-1.D0)
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      IF (IOPT.EQ.1) C=2.D0*PI
      DO 20 I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
  20  CONTINUE

      return
      END subroutine reclat


C ====================================================================
C Write out a complex matrix in a format so that it can be plotted
C directly with Mathematica using e.g. ListDensityPlot
C ====================================================================


      subroutine writezmat(N,name,Z)
      implicit none
      
      integer N
      complex*16 Z(N,N)
      character*4 name ! Use 4 chars as a name

      integer i,j,unit1
      external io_assign,io_close


c
c     BEGIN
c
      unit1=10

      open (unit1,
     &     file=name,
     &     status='unknown')
      

      write(unit1,*) '   '
      write(unit1,*) name,' ={'
      do j=1,N-1
         write(unit1,*) '{'
         write(unit1,'(F12.8,a3,F12.8,a1)') 
     .        (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),',',i=1,N-1)
         write(unit1,'(F12.8,a3,F12.8,a2)') (dreal(Z(i,j)),'+I*',
     .        dimag(Z(i,j)),'},',
     .        i=N,N)
      end do  
      j=N
      write(unit1,*) '{'
      write(unit1,'(F12.8,a3,F12.8,a1)') 
     .     (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),',',i=1,N-1)
      write(unit1,'(F12.8,a3,F12.8,a1)') 
     .     (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),'}',
     .     i=N,N)
      write(unit1,*) '};'
         


      return
      end subroutine writezmat



c======================================================================




      END

