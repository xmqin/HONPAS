C ##################################################################
C ##            read-in header of Greens function file            ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

      subroutine read_green(joutfile,jgfu,
     .     EFermi,NEn,contour,wgf,
     .     nua,NA1,NA2,nq,ng,wq,q,kxy,iq)


      implicit none

      real*8 EPS
      parameter(EPS=1d-7)

      logical PRINTALOT
      parameter(PRINTALOT=.FALSE.)
c      parameter(PRINTALOT=.TRUE.)

c INPUT
      integer joutfile          ! out file
      integer jgfu               !unit of gf-file
C FDN kxy(2) ==> kxy(3)
      real*8 kxy(3)             !surface kxy point used for checking
C FDN


c these are the values expected:
c they will be compared to the read-in and ERROR messages will appear
      real*8 efermi             !The required Fermi energy
      integer NEn,ng
      complex*16 contour(NEn),wgf(NEn)
      integer nua,NA1,NA2       !no. atoms in uc and no. repetitions in A1,A2
      integer nq                !no. q-points
      logical tleft



c READ-IN values
      real*8 efermii            !The required Fermi energy
      integer NEni,ngi
      complex*16, dimension(:), allocatable ::
     . contouri,wgfi

      integer nqi
c      real*8 q(2,nq),wq(nq)  !q-point=k_|| point and weigth

      real*8, dimension (:,:), pointer:: q
      real*8, dimension (:), pointer:: wq

      character*70 gftitle

c Helpers..

      complex*16 ctmp
      integer iEn,iq,i,j
      real*8 qtmp(2),wqtmp
      logical err

c
c     BEGIN
c      

      read(jgfu) gftitle
      read(jgfu) EFermii,NEni
      read(jgfu) nua,NA1,NA2,nqi

      if(PRINTALOT) write(joutfile,*) 'read GF: ',gftitle

      if(NEni .ne. NEn)  then
         write(joutfile,*) 'read GF: ERROR: NEn=',NEni,' expected:', NEn
         stop 'read GF: read-in ERROR: No. energypoints not OK'
      end if


      allocate(contouri(nen))
      allocate(wgfi(nen))

C FDN Surface q points =1 (no repeated cells) 
      nq=1
C FDN  
C FDN q(2,nqi)--> q(3,nqi)
      allocate(q(3,nqi))
C FDN
      allocate(wq(nqi))


       
      read(jgfu) contouri,wgfi,q,wq

      read(jgfu) ng

     
      err=.false.
      if (dabs((q(1,iq)-kxy(1)) ) > EPS ) err=.true.
      if (dabs((q(2,iq)-kxy(2)) ) > EPS ) err=.true.
      if (dabs((q(3,iq)-kxy(3)) ) > EPS ) err=.true. 
      if(err) then
         write(joutfile,*) 'read_green: Unexpected q-points'
         write(joutfile,'(5I,3F12.7)') iq,q(1,iq),q(2,iq),q(3,iq)  
         write(joutfile,'(5I,3F12.7)') iq,kxy(1),kxy(2),kxy(3) 
         stop 'read_green: Unexpected q-points'
      end if
C FDN

 
c     check on contours.


      do iEn=1,NEn
         if(PRINTALOT) then
            write(joutfile,'(2F10.6,2x,2F10.6)') 
     .           contouri(iEn), contour(iEn)
         end if
         ctmp=contouri(iEn)-contour(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(joutfile,*) ' Warning: contours differ by >', EPS
         end if
         if(cdabs(ctmp).GT.10d0*EPS) then 
            write(joutfile,*) 
     &           ' ERROR: contours differ by >', 10.d0*EPS
            stop ' ERROR: read_green: contours differ'
         end if
      end do

      do iEn=1,NEn
         ctmp=wgfi(iEn)-wgf(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(joutfile,*) 
     &           ' ERROR: contour weights differ by >',EPS
            stop ' ERROR: read_green: contour weights differ'
         end if
      end do

      


      deallocate(contouri)
      deallocate(wgfi)

      return
      end





