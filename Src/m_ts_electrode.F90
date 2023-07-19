MODULE m_ts_electrode
!
! Routines that are used for Electrodes GFs calculations
!
!=============================================================================
! CONTAINS:
!          1) calc_green
!          2) green
!          3) sethhm2


  implicit none

  public :: calc_green, green, sethhm2

  private

  CONTAINS


       subroutine calc_green(nv,zenergy,h00,s00,h01,s01, &
       gs,zdos,joutfile,tleft)

       use m_ts_aux_rout, only : csolveg
       use precision, only : dp

       implicit none
       integer nv

       integer nv2,iter
       integer ierr             !error in inversion
       integer i,j,ic,ic2
       integer joutfile
       logical tleft

       complex(dp) a,b,zdos
       real(dp) ro
       real(dp), parameter :: accur=1.d-15
       complex(dp) ZEnergy 
       complex(dp) h00(0:nv*nv-1),s00(0:nv*nv-1)
       complex(dp) h01(0:nv*nv-1),s01(0:nv*nv-1)
       complex(dp) gs(0:nv*nv-1)


       integer, dimension (:), allocatable:: ipvt
       complex(dp), dimension (:), allocatable:: &
         rh,rh1,rh3,alpha,beta,ab,ba,gb,gs2

       real(dp) Pi
       parameter(Pi=3.14159265358979323846264338327950288419717d0)

       allocate(ipvt(nv))
       allocate(rh(0:2*nv*nv))
       allocate(rh1(0:2*nv*nv))
       allocate(rh3(0:4*nv*nv))
       allocate(alpha(0:nv*nv-1))
       allocate(beta(0:nv*nv-1))
       allocate(ba(0:nv*nv-1))
       allocate(ab(0:nv*nv-1))
       allocate(gb(0:nv*nv-1))
       allocate(gs2(0:nv*nv-1))

       a=(1.d0,0.d0)
       b=(0.d0,0.d0)
       nv2 =2*nv

! FDN
! gb = Z*S00-H00
! alpha = -(Z*S01-H01)        
      do i=0,nv*nv-1
        gb(i) = zenergy*s00(i)-h00(i)
        alpha(i) = h01(i)-zenergy*s01(i)
      end do

! FDN
! gs = Z*S00-H00
! gs2 = Z*S00-H00
      do i=0,nv*nv-1
        gs(i) = gb(i)
        gs2(i) = gb(i)
      end do                

! FDN
! beta = -(Z*S10-H10)
      do j=0,nv-1
       do i=0,nv-1
        ic = i + nv*j
        ic2 = j + nv*i
        beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
      end do


        iter=0
1000    continue
        iter=iter+1


! FDN
! nv2=2*nv
! rh = -(Z*S01-H01) ,j<nv
! rh = -(Z*S10-H10) ,j>nv
      do j=0,nv2-1
         do i=0,nv-1
              ic =i + j*nv
              ic2=i + (j - nv)*nv
              if(j.lt.nv)then
               rh(ic) = alpha(ic)
              else
               rh(ic) = beta(ic2)
              endif                
         end do                 
      end do                  

! FDN
! rh3 = Z*S00-H00
      do i=0,nv*nv-1
        rh3(i) = gb(i)
      end do                

! FDN
! rh =  rh3^(-1)*rh
! rh = t0
      call csolveg(nv,nv2,rh3,rh,ipvt,ierr) 

      IF(IERR.NE.0) THEN
         write(joutfile,*) 'ERROR: calc_green 1 MATRIX INVERSION FAILED'
         write(joutfile,*) 'ERROR: LAPACK INFO = ',IERR
      END IF


! FDN
! nv2=2*nv
! rh1 = -(Z*S01-H01) ,j<nv
! rh1 = -(Z*S10-H10) ,j>nv
! a = 1
! b = 0
      do j=0,nv-1
         do i=0,nv2-1
             ic =i + j*nv
             ic2 =i-nv + j*nv
              if(i.lt.nv)then
               rh1(i + nv2*j) = alpha(ic)
              else
               rh1(i + nv2*j) = beta(ic2)
             end if
         end do                
      end do
! FDN
! rh3 = 1*rh1*rh + 0*rh3                   
! rh3 = -(Z*S01-H01)*t0
      call zgemm('N','N',nv2,nv2,nv,a,rh1,nv2,rh,nv,b,rh3,nv2)

! FDN
! aplha = -(Z*S01-H01)*t0
! ba = -(Z*S10-H10)*t0b
      do j=0,nv-1
         do i=0,nv2-1
             ic =i + j*nv
             ic2 =i-nv + j*nv
              if(i.lt.nv)then
                alpha(ic) = rh3(i + nv2*j) 
              else
                ba(ic2) = -rh3(i + nv2*j) 
             end if
         end do                 
      end do                   
      do j=nv,nv2-1
         do i=0,nv2-1
              ic=i + (j - nv)*nv
              ic2=i - nv + (j - nv)*nv
              if(i.lt.nv)then
                ab(ic)= -rh3(i + nv2*j) 
              else
                beta(ic2)= rh3(i + nv2*j) 
             end if
         end do                 
      end do                   
 
      do i=0,nv*nv-1
       gb(i) =  gb(i) + ba(i) + ab(i)
       gs(i) =  gs(i) + ab(i) 
       gs2(i) =  gs2(i) + ba(i) 
      end do                    

      ro =-1.0
      do j =0,nv*nv-1
        ro =max(ro,dreal(ab(j))**2+dimag(ab(j))**2)
      end do                   
      ro =dsqrt(ro)

      if(ro.gt.accur) go to 1000


      do i=0,nv*nv-1
        rh3(i) = gs(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,ierr)

      IF(IERR.NE.0) THEN
         write(joutfile,*) 'ERROR: calc_green 2 MATRIX INVERSION FAILED'
         write(joutfile,*) 'ERROR: LAPACK INFO = ',IERR
      END IF

      do i=0,nv*nv-1
        gs(i) = rh(i)
      end do



      do i=0,nv*nv-1
        rh3(i) = gs2(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,ierr)


      IF(IERR.NE.0) THEN
         write(joutfile,*) 'ERROR: calc_green 3 MATRIX INVERSION FAILED'
         write(joutfile,*) 'ERROR: LAPACK INFO = ',IERR
      END IF

      do i=0,nv*nv-1
        gs2(i) = rh(i)
      end do


      do i=0,nv*nv-1
        rh3(i) = gb(i)
        rh(i) = 0.0
      end do

      do j=0,nv-1
        rh(j*(nv+1)) = 1.d0
      end do

      call csolveg(nv,nv,rh3,rh,ipvt,ierr)

      IF(IERR.NE.0) THEN
         write(joutfile,*) 'ERROR: calc_green 4 MATRIX INVERSION FAILED'
         write(joutfile,*) 'ERROR: LAPACK INFO = ',IERR
      END IF

      do i=0,nv*nv-1
        gb(i) = rh(i)
      end do


!      ----      DOS     -----

      do i=0,nv*nv-1
        alpha(i) = h01(i)-zenergy*s01(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,gs2,nv,alpha,nv,b,ab,nv)
      do i=0,nv*nv-1
        alpha(i) = ab(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,alpha,nv,gb,nv,b,ab,nv)


      do j=0,nv-1
       do i=0,nv-1
        ic = i + nv*j
        ic2 = j + nv*i
        beta(ic) = dconjg(h01(ic2))-zenergy*dconjg(s01(ic2))
       end do
      end do
      call zgemm('N','N',nv,nv,nv,a,gs,nv,beta,nv,b,ba,nv)
      do i=0,nv*nv-1
        beta(i) = ba(i)
      end do
      call zgemm('N','N',nv,nv,nv,a,beta,nv,gb,nv,b,ba,nv)


      do i=0,nv*nv-1
        rh3(i) = 0.0
      end do

      call zgemm('N','N',nv,nv,nv,a,gb,nv,s00,nv,b,rh3,nv)
      call zgemm('N','C',nv,nv,nv,a,ab,nv,s01,nv,a,rh3,nv)
      call zgemm('N','N',nv,nv,nv,a,ba,nv,s01,nv,a,rh3,nv)


      zdos =0.0

      do j=0,nv-1
        zdos = zdos + (rh3(j*(nv+1)))
      end do



      if(tleft) then

       do i=0,nv*nv-1
        gs(i) =  gs2(i) 
       end do


      endif

       deallocate(ipvt)
       deallocate(rh)
       deallocate(rh1)
       deallocate(rh3)
       deallocate(alpha)
       deallocate(beta)
       deallocate(ba)
       deallocate(ab)
       deallocate(gb)
       deallocate(gs2)

      return
      end subroutine calc_green

!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------



! ##################################################################
! ## Driver subroutine for calculating the (ideal)                ##
! ## Left surface Greens function.                                ##          
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
! FDN nspin added as dummy
      subroutine green(joutfile, &
         NEn,contour,wgf,efermi,zbulkdos,tjob,nspin)
! FDN      

      use fdf, only : fdf_convfac, fdf_integer, fdf_string
      use precision, only : dp
! FDN
      use m_ts_kpoints, only : ts_nkpnt, ts_kpoint, ts_kweight
      use m_ts_options, only : GFFileL, GFFileR, calcGF
      use files, only: label_length
! FDN

!======================================================================
      implicit none
!======================================================================

! Dimensions for input to green

      logical PRINTALOT
!      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)
      
!=======================================================================
! INPUT:      
      character*20 slabel       ! System Label (to name output files)
      

      

      logical tjob !True if Left, False if Right
      
      integer joutfile          !unit-number of out-file

      character(len=label_length+5) hsfile        !name of HS-input file
      integer NEn ! no. contour points
      complex(dp) contour(NEn),wgf(NEn) !energy contours and weights for GF
      real(dp) efermi             ! the Fermi-energy we REQUIRE the electrode
                                ! to have (e.g. when a voltage is applied)
      integer ispin             !spin number for which to generate gf
!=======================================================================

!   k_|| == q-points:
      real(dp), pointer:: q(:,:),wq(:) ! k_|| and their weights 
      integer na1,na2           ! Replication of unitcell33
      integer nq                ! no. q-points (<= na1*na2 for gamma) 
      real(dp) kpoint(3)          !3D k-point (q,kz)
      data kpoint /0d0, 0d0, 0d0/
!------------------------------------------------------------------     
! Hamiltonian/Overlap for given k

      integer nspin
      integer nua !no atom in uc
!      integer lasto(0:NG1) ! Index of last orbital of each atom
      integer, dimension (:), pointer:: lasto
! (we use that the number of atoms < NG1)


! ==================================================================
!     OUTPUT GAA:            


      character*33 gffile        ! name of gffile
      integer jgfu               ! output gf-file unit
!      complex*16 HAA(NGAA,mqpt)  ! Real-space bulk Hamiltonian
!      complex*16 SAA(NGAA,mqpt)  ! Real-space bulk Overlap
!      complex*16 GAAq(NGAA,mqpt) ! Inverse ideal GF (z*SAA - HAA - SigmaAA)

!      complex*16 GS(ngaa1)
      complex(dp), allocatable :: GS(:),HAA(:,:),SAA(:,:), GAAq(:,:)

      integer nucuse             !Truncate GF to *last/first* nucuse atoms
!                                !         for   *left/right*

      integer, dimension (:), allocatable :: lastou
!      integer lastou(0:NG1) ! Index of last orbital of each *used* atom
! (we use that the number of atoms < NG2)
      integer nou               !no. used orbitals     

      complex(dp) zbulkdos(NEn) 

!=======================================================================
!     Helpers, workspace, tempos etc...

      character*70 gftitle      !title of gf
      

      logical tinit,tlast

!     LEFT/RIGHT sign in GF calc.: Exp(i*LRsign*Kz*Z):       
      integer LRSign

      character*5 tag
      character stag
      character*6  gfjob

      logical exist1
      logical tdos,tkham
      
      complex(dp), dimension (:), pointer:: H00, S00, H01, s01

!      complex*16 h00(ngaa1),s00(ngaa1)
!      complex*16 h01(ngaa1),s01(ngaa1)
!      complex*16 ab(ngaa1),ba(ngaa1)
       complex(dp), allocatable :: ab(:),ba(:)


      integer i,l1,l2,ia,ia2
      integer iqpt,iq
      integer ng1tmp
      

      integer ngaa,ngaa1
      integer NG1      ! Number of basis  orbitals
      integer NG2      ! Number of orbitals used

      real(dp)  factor
      real(dp), allocatable:: eig(:)
      

      complex(dp) ZEnergy
      complex(dp) ZSEnergy
      complex(dp) zdos
      integer iEn

! FDN Celula unitaria, kscell and kdispl
      real(dp) cell(3,3)
      integer kscell(3,3)
      real(dp)  kdispl(3)
      integer allocstat
! FDN

      external io_assign,io_close

! FDF-stuff:
      character*33 paste,itemfdf     
!      real*8 fdf_convfac
      external paste

!=======================================================================
! BEGIN:
!=======================================================================

       tdos =.true.

! FDN Put ispin=1 for the initialization procedures
       ispin=1
! FDN

       factor =fdf_convfac('Ry','eV')

      if(tjob)then
!!   left
       gfjob='LEFT  '
       LRSign=1
       tag='Left'
       stag='L'
       gffile = GFFileL
      else
!!   left
       gfjob='RIGHT '
       LRSign=-1
       tag='Right'
       stag='R'
       gffile = GFFileR
      endif
!
      write(joutfile,*) 'Begin ',gfjob
      
! --------------

! output .GF files to this name


      itemfdf =paste('TS.ReplicateA1',tag)
      NA1=fdf_integer(itemfdf,1)
      itemfdf =paste('TS.ReplicateA2',tag)
      NA2=fdf_integer(itemfdf,1)

      itemfdf =paste('TS.HSFile',tag)
      hsfile = fdf_string(itemfdf,paste(slabel,'.TSHS'))

! --------------
!
! on first call the H,S read-in is initialized:
!
      tinit=.true.
      tlast=.false.
      tkham=.false.

      nullify(h00)
      nullify(s00)
      nullify(h01)
      nullify(s01)
      nullify(lasto)
! FDN initialize cell variable
      cell=0.d0
      kscell=0.d0
      kdispl=0.d0
! FDN

! FDN cell,kscell, kdispl variables added
      call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,hsfile, &
          nua,lasto,ng1,nspin,cell,kscell,kdispl, &
          H00,s00,h01,s01)    ! ->
! FDN


      ng1tmp = ng1
      ngaa1=ng1*ng1
      allocate(GS(ngaa1))
!
!     check dimensions:
!

      itemfdf = paste('TS.NumUsedAtoms',tag) 
      nucuse = fdf_integer(itemfdf,nua) !default use all atoms in uc.


!     Truncate lasto to used atoms: lastou:

      allocate(lastou(0:ng1))

      lastou(0)=0
      ia2=0
         if(LRsign .EQ. 1) then !Left
!     use nucuse *last* atoms           
            do ia=nua-nucuse ,nua
               lastou(ia2)=lasto(ia)
               ia2=ia2+1
            end do
         else                   !Right
!     use nucuse *first* atoms            
            do ia=0,nucuse
               lastou(ia2)=lasto(ia)
               ia2=ia2+1               
            end do
         end if                 !L or R            
         
! No. used orbitals in uc:
         nou=0
         do ia=1,nucuse
            nou=nou + (lastou(ia)-lastou(ia-1))
         end do

         ng2 = nou
         ngaa=ng2*ng2

!
!     Get  k_|| - points == q-points:
!
! FDN cell,kscell,kdispl added as dummys
!         call mkqgrid(joutfile,NA1,NA2,nq,q,wq,cell,kscell,kdispl)
        nullify(q, wq)
        nq=ts_nkpnt
        allocate(q(3,ts_nkpnt))
        q=ts_kpoint
        allocate(wq(ts_nkpnt))
        wq=ts_kweight
! FDN
      

!--------------------------------------------------------------
         
      write(joutfile,*) 'Size: ',ngaa

!=======================================================================
!
! Read-in bulk parameters and do calculation:
!
!=======================================================================

!cccccccccccccccccccccccccccccccxccccccccccccccccccccccccccccccc
! in future versions green will be called from gengf with
! ispin as argument and generate the GF for spin=ispin
! FDN
!      ispin = 1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! FDN
!      write(joutfile,*) 'GF: Spin number -> ',ispin
! FDN
    
      inquire(file=gffile,exist=exist1)

! If it does not find the file, calculate the GF
      if(.not.calcGF .and. .not.exist1) calcGF = .true.

      if (calcGF) then
      if(exist1) then
       write(joutfile,*) gffile,' already exist, will be overwritten!'
      end if

      call io_assign(jgfu)
      OPEN(FILE=gffile,UNIT=jgfu, &
          FORM='UNFORMATTED')

       
      gftitle = fdf_string('TS.GFTitle','genGF')
!

!     Write initial info to gffile:
!
         write(jgfu) gftitle
         write(jgfu) EFermi,NEn
         write(joutfile,*)'Efermi: ',EFermi
         write(jgfu) nucuse,NA1,NA2,nq
         write(jgfu) contour,wgf,q,wq
         write(jgfu) NG2



!-------------------
         tkham=.true.

         kpoint = 0.0

! FDN cell variable added
         call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,hsfile, &
             nua,lasto,NG1tmp,nspin,cell,kscell,kdispl, &
             H00,S00,H01,S01)
! FDN


! FDN Checking allocation ...
!      write(joutfile,*) allocated(eig),allocated(ab),allocated(ba)
      allocate(eig(ng1), STAT=allocstat)
      if ( allocstat .ne. 0 ) then
       write(joutfile,*) 'Could not allocate eig(ng1) in routine green'
       stop
      end if         
!      write(joutfile,*) allocated(eig),allocated(ab),allocated(ba)
      allocate(ab(ngaa1), STAT=allocstat)
      if ( allocstat .ne. 0 ) then
       write(joutfile,*) 'Could not allocate ab(ngaa1) in routine green'
       stop
      end if
!      write(joutfile,*) allocated(eig),allocated(ab),allocated(ba) 
      allocate(ba(ngaa1), STAT=allocstat)
      if ( allocstat .ne. 0 ) then
       write(joutfile,*) 'Could not allocate ba(ngaa1) in routine green'
       stop
      end if
! FDN



! FDN
      deallocate(ab)
      deallocate(ba)
! FDN


      deallocate(eig)

!------------------

!     initialize ideal GF for all Energies:
!     ^^^^^^^^^^

! FDN Dimension changed (NGAA,nq) ---> (NGAA,1)
       allocate(HAA(NGAA,1))
       allocate(SAA(NGAA,1))
       allocate(GAAq(NGAA,1))
! FDN

!=========================================================
! FDN Spin Loop
!=========================================================

      do ispin=1,nspin

! FDN
      write(joutfile,*) 'GF: Spin number -> ',ispin
! FDN



!     initialize GAAq,HAA,SAA
! FDN Loop over nq to 1
          do iq=1,1
             do i=1,NGAA
               GAAq(i,iq)=dcmplx(0d0,0d0)
               HAA(i,iq)=dcmplx(0d0,0d0)
               SAA(i,iq)=dcmplx(0d0,0d0)
             end do           !i
          end do              !iq
! FDN



!     and bulk density of states

! FDN Loop placed here
! ==========================================================================
!     Loop over q=k_|| - points:
      do iqpt=1,nq


        kpoint(1) = q(1,iqpt)    !in (xy)
        kpoint(2) = q(2,iqpt)    !in (xy)
! FDN In general k// may have z component
        kpoint(3) = q(3,iqpt)
! FDN

! ==========================================================================


      do iEn=1,NEn
         zbulkdos(iEn)=dcmplx(0d0,0d0)
      end do


! ======================================================================
      do iEn = 1, NEn  ! loop over ALL Contour points
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
        ZEnergy = contour(iEn)

! FDN Information now written here
        write(jgfu) iEn,contour(iEn),wgf(iEn),iqpt
! FDN



            h00 =0.0
            s00 =0.0
            h01 =0.0
            s01 =0.0


      tinit=.false.
      tkham=.false.

! FDN cell variable added
      call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,hsfile, &
          nua,lasto,NG1tmp,nspin,cell,kscell,kdispl, &
          H00,s00,h01,s01)    ! ->
! FDN




       if(iEn.eq.1) then
        if(LRsign .EQ. -1) then !Right
            i=0
            do l2 = 1,NG2
               do l1 = 1,NG2
                  i=i+1
! FDN Second variable put as one
                  HAA(i,1) = H00(l1+NG1*(l2-1))
                  SAA(i,1) = S00(l1+NG1*(l2-1))
! FDN
               end do           ! l1
            end do              ! l2
        else                   !left
            i=0
           do l2 = 1,NG2
              do l1 = 1,NG2
                  i=i+1

! FDN Second variable put as one
                  HAA(i,1) = H00(l1+(NG1-NG2)+NG1*(l2+(NG1-NG2)-1))
                  SAA(i,1) = S00(l1+(NG1-NG2)+NG1*(l2+(NG1-NG2)-1))
! FDN
              end do           ! l1
           end do              ! l2
        end if                 !L or R

!
!     Shift so efermi is energy zero
!
! FDN Second variable put as one
! Commented since it does not shift, as it is ... shifted in sethhm2
! efermi is already  zero here !
!            call zaxpy(NGAA,dcmplx(efermi,0.d0),
!     .           SAA(1,1),1,HAA(1,1),1)
! FDN

       end if             

        zsenergy = zenergy-efermi



         call calc_green(ng1,zsenergy,h00,s00,h01,s01, &
       gs,zdos,joutfile,tjob)


! FDN In case of spin, this sum must be checked
        if(tdos) zbulkdos(iEn) = zbulkdos(iEn) + wq(iqpt)*zdos

        if(LRsign .EQ. -1) then !Right
            i=0
            do l2 = 1,NG2
               do l1 = 1,NG2
                  i=i+1
! FDN Second variable put as one
                  GAAq(i,1) = gs(l1+NG1*(l2-1))
! FDN
               end do           ! l1
            end do              ! l2
         else                   !left
            i=0
            do l2 = 1,NG2
               do l1 = 1,NG2
                  i=i+1
! FDN Second variable put as one
                  GAAq(i,1) = gs(l1+(NG1-NG2)+ NG1*(l2+(NG1-NG2)-1) )
! FDN
               end do           ! l1
            end do              ! l2
         end if                 !L or R

! FDN outside energy loop
!C ======================================================================
!      end do                    !loop over k_|| - points (iqpt)
!C ======================================================================
! FDN

! FDN
!         write(jgfu) iEn,contour(iEn),wgf(iEn)
! FDN
         if(iEn .eq. 1) then
            write(jgfu) HAA
            write(jgfu) SAA
         end if

           write(jgfu) GAAq

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=
               end do           ! the loop over ALL Energies
!=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-

! FDN kpoint loop
! ======================================================================
      end do                    !loop over k_|| - points (iqpt)
! ======================================================================
! FDN


! FDN spin loop                                                            |C     dimension changed from 2 to 3
! ======================================================================     |      allocate(q(3,nqi))
      end do                    !loop over spin
! ======================================================================     |      allocate(wq(nqi))
! FDN

      write(joutfile,*) 'Got ',gfjob,' Electrode GF'  

      call io_close(jgfu)

      write(joutfile,*) 'Done '

      tlast=.true.

! FDN cell variable added
      call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,hsfile, &
          nua,lasto,NG1tmp,nspin,cell,kscell,kdispl, &
          H00,s00,h01,s01)    ! ->
! FDN

      deallocate(HAA)
      deallocate(SAA)
      deallocate(GAAq)
      end if ! calcGF

      
      deallocate(q)
      deallocate(wq)
      deallocate(GS)
      deallocate(lastou)
      deallocate(H00)
      deallocate(s00)
      deallocate(H01)
      deallocate(s01)
      deallocate(lasto)  
      
      
! ======================================================================
      return
      end subroutine green
! ======================================================================





!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


! ##################################################################
! ##               Setup Hamiltonian in k-space                   ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
!
! This routine use TSiohs1.f which can both read and write out the
! TSHS1 format. The same routine is used by SIESTA to write out.
! The label "1" after setupkham leaves room for other formats, like
! "2" etc..

! FDN Modified to include the unit cell of the electrodes as  dummy,
! also kscell and kdispl
! The original version was like this:
!      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,
!     &     hsfile,nua,lasto,nuo,nspin,Hk,Sk,Hk2,Sk2)    ! ->
! FDN

      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin, &
          hsfile,nua,lasto,nuo,nspin,cell,kscell,kdispl, &
          Hk,Sk,Hk2,Sk2)    ! ->


      use m_ts_io, only: ts_iohs
      use m_ts_kpoints, only: ts_gamma_scf
      use files, only: label_length
      use precision, only: dp
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

!=======================================================================
      real(dp) EPS
      parameter(EPS=1.0d-8)
!=======================================================================


! INPUT

      integer joutfile          !info-out file unit
      logical tkham             ! true if full H_k hamiltonian is generated

!      integer nuo               !No. states in unitcell (expected in read-in)
      real(dp) kpoint(3) 
      integer ispin             !set up H for ispin
      character(len=label_length+5) hsfile       !H,S parameter file


!-----------------------------------------------------------------------
! READ-IN from HS-file/OUTPUT


      integer nua               ! No. atoms in unitcell
      integer nuo            ! Number of basis  orbitals
      integer nspin          ! Spin polarization (1 or 2)
! FDN Electrode spin
      integer Enspin
! FDN

!      integer numh(maxuo)        ! Number of nonzero elements of each row
!                               ! of hamiltonian matrix
!      integer listh(maxno,maxuo) ! Nonzero hamiltonian-matrix element column
!                               ! indexes for each matrix row
!      integer indxuo(maxo)      ! Index of equivalent orbital in unit cell
!                               !  Unit cell orbitals must be the first in
!                               ! orbital lists, i.e. indxuo.le.nuo, with
!                               ! nuo the number of orbitals in unit cell  

!      real*8  H(maxno,maxuo,nspin) ! Hamiltonian in sparse form
!      real*8  S(maxno,maxuo)     ! Overlap in sparse form
      real(dp)  qtot              ! Total number of electrons
      real(dp)  temp              ! Electronic temperature for Fermi smearing

      logical Gamma             ! true if Gamma

!      real*8  xij(3,maxno,maxuo) ! Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point ..
!                                but must be written here!!)

      real(dp) cell(3,3)          ! unit cell
      real(dp) ef                 ! Fermi energy from parameterfile

!      integer nua               ! No. atoms in unitcell
!      integer isa(maxua)         ! Species index of each atom
!      real*8 xa(3,maxua)         ! Atomic coordinates (Bohr)


!      integer, pointer ::   isa(:),numh(:),indxuo(:),listh(:,:)
!      real*8, pointer ::   H(:,:,:),S(:,:),xij(:,:,:),xa(:,:),efs(:)

      integer, dimension (:), pointer:: listh, listhptr, &
                             numh,indxuo,isa,lasto
      real(dp), dimension (:,:), pointer:: H,xij,xa
      real(dp), dimension (:), pointer:: S

! FDN Temporary array xijtemp so that there will be a complex
!     phase dependence on unit cell
      real(dp), dimension (:,:), pointer:: xijtemp
! FDN

! FDN
      integer kscell(3,3)
      real(dp)  kdispl(3)
      character(8) :: iotask
! FDN


!-----------------------------------------------------------------------
! OUTPUT
!      complex*16 Hk(nuo*nuo), Sk(nuo*nuo)
!      complex*16 Hk2(nuo*nuo), Sk2(nuo*nuo)

      complex(dp), dimension (:), pointer:: Hk, Sk, Hk2, Sk2
!      integer, dimension (:), pointer:: lasto
!      integer lasto(0:maxua)   ! Index of last orbital of each atom
!-----------------------------------------------------------------------
! Helpers
!      real*8 xo(3,maxuo)        ! Atomic coordinates (Bohr)
      real(dp), allocatable ::   xo(:,:) 

      integer nuotot,notot,maxnh
      integer ia
      integer i,j,io,jo,iuo,juo
      real(dp) k(3),kxij,rcell(3,3)
      real(dp) recell(3,3)
      complex(dp) cphase

!      integer ix(maxnh)
      integer, dimension (:),allocatable:: ix
      integer icoi,icoa
      integer iprop,inn,it,in,ind
      real(dp)   xc      
      logical tinit,tlast
      integer :: icrap1, icrap2 ! dummy variables
! FDN
      logical ts_gamma
! FDN
!-----------------------------------------------------------------------
! SAVED:
      save numh,listh,indxuo,Gamma
      save H,S,ef,rcell,xij,listhptr,nuotot
      save ix

!=======================================================================
! BEGIN
      
       iprop = 3
       icoa = 0
       icoi = 0

!-----------------------------------------------------------------
!     Read-in Hamiltonian/Overlap parameters from HS-lattice-file
!     Only read-in firsttime:
!    
        
       if(tinit) then
!-----------------------------------------------------------------         
              nullify(H)
              nullify(S)
              nullify(xij)
              nullify(indxuo)
              nullify(listh)
              nullify(listhptr)
              nullify(numh)
              nullify(xa)

! FDN xijtemp var
              nullify(xijtemp)
! FDN
       


! FDN kscell, kdispl added;nspin --> Espin
!      call TSiohs2( 'read',hsfile, gamma, nua, nuotot,notot,Enspin, &
!                   maxnh,numh, listhptr, listh, H, S, qtot, temp, &
!                   xij, indxuo, ef, efs, cell, isa, lasto, xa, &
!                   kscell, kdispl)

       iotask='read'
       call ts_iohs(iotask,gamma, .false., nuotot, notot, Enspin, indxuo, &
                    maxnh, numh, listhptr, listh, H, S, qtot, temp, xij, &
                    label_length+5, hsfile, nua, lasto, isa, ef, cell, &
                    kscell, kdispl, ts_gamma, xa, icrap1, icrap2)  
! FDN
    
! FDN Check if electrode has the same spin
      if( nspin.ne.Enspin ) then
       write(*,*) 'Differente spin in Electrode !!'
       write(*,*) 'Stoping code'
       stop
      end if
  

           nuo = nuotot
              allocate(Hk(nuo*nuo))
              allocate(Sk(nuo*nuo))
              allocate(Hk2(nuo*nuo))
              allocate(Sk2(nuo*nuo))
     
! FDN  Alocate xijtemp
              allocate(xijtemp(3,size(xij)/3))
! FDN
         
         write(joutfile,*) 'unit cell:'

         do j=1,3
            write(joutfile,'(3F8.4)') (cell(i,j),i=1,3)
         end do
         
         call reclat(cell,rcell,1) !reciprocal of cell incl. 2Pi!

         call reclat(cell,recell,0)


         if(.not. Gamma) then

            
          allocate(xo(3,nuo))
          
          if (allocated(ix))deallocate(ix)
          
        
          allocate (ix(maxnh))
        
!
!     Transform xij so there is no k-dep. phase within uc.
!

!     ... but first find orbital coordinate

            do ia=1,nua
               do iuo=lasto(ia-1)+1,lasto(ia)
                  xo(1,iuo)=xa(1,ia)
                  xo(2,iuo)=xa(2,ia)
                  xo(3,iuo)=xa(3,ia)
               end do           !iuo
            end do              !ia in uc


! FDN copy xij to xijtemp
           xijtemp=xij
! FDN


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

! FDN copy xijtemp to xij
         xij=xijtemp
         deallocate(xijtemp)
! FDN       

         
         deallocate( isa )
         deallocate( xa )
     

         tinit = .false.
         return
!-----------------------------------------------------------------
         end if                    !tinit
!-----------------------------------------------------------------
         if(tlast) then
            deallocate(H)
            deallocate(S)
            deallocate(xij)
            deallocate(indxuo)
            deallocate(listh)
            deallocate(numh)                     
            return
         endif
               

         k=kpoint

!-----------------------------------------------------------------
      if (tkham)then
!-----------------------------------------------------------------

!
! Setup H,S for this k-point:
!
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
              kxij = (k(1) * xij(1,ind) + &
                     k(2) * xij(2,ind) + &
                     k(3) * xij(3,ind) )
              cphase = exp(dcmplx(0d0,1d0)*kxij)
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
!
!     Symmetrize and *make EF the energy-zero*!!!
!
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))
! FDN 
           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) ) &
                - ef*Sk(it)
           Hk(in) =  dconjg(Hk(it))

       
        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in)=Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
! FDN
         Hk(in)=Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in)) &
                - ef*Sk(in) 
      enddo

!-----------------------------------------------------------------
      else                      !tkham
!-----------------------------------------------------------------

! Setup Transfer H,S for this k||-point:
!
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
               kxij = (k(1) * xij(1,ind) + &
                      k(2) * xij(2,ind) + &
                      k(3) * xij(3,ind) - &
                      k(iprop) *xij(iprop,ind) )
              cphase = exp(dcmplx(0d0,1d0)*kxij)

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

!
!     Symmetrize and *make EF the energy-zero*!!!
!
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))
! FDN
           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) ) &
                - ef*Sk(it)

           Hk(in) =  dconjg(Hk(it))

        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in)=Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
! FDN
         Hk(in)=Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in)) &
                - ef*Sk(in)

      enddo


      do iuo = 1,nuo
         do juo = 1,nuo
          in = iuo+(juo-1)*nuo
! FDN
          Hk2(in)=Hk2(in) - ef*Sk2(in) 
         enddo
      enddo

!-----------------------------------------------------------------
       endif                      !tkham
!-----------------------------------------------------------------
       
! ===============================================================
       return  
       END subroutine sethhm2
! ===============================================================






END MODULE m_ts_electrode
