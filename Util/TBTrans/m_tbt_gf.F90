module m_tbt_gf

! ----------------------------------------------------------------------
! Module that contains the routines and variables used for to obtain 
! the electrodes surface green's functions and self energies
!
! F.D.Novaes, fdnovaes@icmab.es
! ----------------------------------------------------------------------
! CONTAINS:
! 1) alloc_gf_vars
! 2) cp_gf_vars
! 3) green
! 4) sethhm2


implicit none

save

logical :: LFrstTime=.True., RFrstTime=.True., LJob, RJob, WGFFiles=.false.

complex*16, allocatable, dimension (:) :: LH00, LS00, LH01, LS01
complex*16, allocatable, dimension (:) :: RH00, RS00, RH01, RS01
complex*16, allocatable, dimension(:) :: LGS,RGS
complex*16, allocatable, dimension(:,:) :: LHAA,LSAA,LGAAq,RHAA,RSAA,RGAAq
integer, pointer, dimension (:) :: Llasto,Rlasto
!complex*16, allocatable, dimension(:) ::  Lzbulkdos,Rzbulkdos

! Saved Variables at green4.F

integer :: Lng1, Rng1,NGL2,Lngaa,NGR2,Rngaa,nuaL,nuaR

! Saved Left Variables at sethhm2
integer :: Lnuotot
integer, allocatable, dimension (:), target :: Llisth, Llisthptr, Lnumh, Lindxuo
integer, allocatable, dimension (:), target :: Lix
double precision, allocatable, dimension (:,:), target :: LH,Lxij
double precision, allocatable, dimension (:), target :: LS,Lefs
logical LGamma             ! true if Gamma
real(8), dimension(3,3) ::  Lrcell

! Saved Right Variables at sethhm2
integer :: Rnuotot
integer, allocatable, dimension (:), target :: Rlisth, Rlisthptr, Rnumh, Rindxuo
integer, allocatable, dimension (:), target :: Rix
double precision, allocatable, dimension (:,:), target :: RH, Rxij
double precision, allocatable, dimension (:), target :: RS,Refs
logical RGamma             ! true if Gamma
real(8), dimension(3,3) ::  Rrcell


contains

subroutine alloc_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs, ix, &
        notot, nuotot, maxnh, nspin, gamma)

double precision, allocatable, dimension (:,:) :: H, xij
integer, allocatable, dimension (:) :: ix
double precision, allocatable, dimension (:) :: S,efs
integer, allocatable, dimension (:) :: listh, listhptr, numh, indxuo

! Dimesnions

integer :: notot, nuotot, maxnh, nspin
logical :: gamma

allocate(H(maxnh,nspin))
allocate(S(maxnh))
allocate(indxuo(notot))
allocate(listh(maxnh))
allocate(listhptr(nuotot))
allocate(numh(nuotot))
allocate(efs(nspin))
allocate(ix(maxnh)) 
if (.not.gamma) allocate(xij(3,maxnh))


end subroutine alloc_gf_vars

subroutine cp_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs,ix,rcell,&
         nuotot, gamma)

use sys, only : die

double precision, dimension (:,:), pointer :: H, xij
double precision, dimension (:), pointer :: S,efs
integer, dimension (:), pointer :: listh, listhptr, numh, indxuo
integer, dimension (:), pointer :: ix
integer :: nuotot
real(8), dimension(3,3) ::  rcell
logical :: gamma


if (LJob) then

   LH=H
   if (.not. gamma) Lxij=xij
   LS=S
   Lefs=efs
   Llisth=listh
   Llisthptr=listhptr
   Lnumh=numh
   Lindxuo=indxuo
   Lix=ix
   Lnuotot=nuotot
   Lrcell=rcell
   Lgamma=gamma


else if (RJob) then
   
   RH=H
   if (.not. gamma) Rxij=xij
   RS=S
   Refs=efs
   Rlisth=listh
   Rlisthptr=listhptr
   Rnumh=numh
   Rindxuo=indxuo
   Rix=ix
   Rnuotot=nuotot
   Rrcell=rcell
   Rgamma=gamma

else

call die("In routine cp_gf_vars: It should be either LJob or RJob !!")

end if

end subroutine cp_gf_vars


! ##################################################################
! ## Driver subroutine for calculating the (ideal)                ##
! ## Left surface Greens function.                                ##          
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################

! FDN Added wq as dummy
      subroutine green(joutfile,nq,q,wq, &
          NEn,contour,wgf,efermi,zbulk,tjob,tbtjob, &
          ispin_fdf, k_fdf, zenergy_fdf,iEn_fdf)
      

! FDN
      use m_tbt_options, only : sppol, Lhsfile, Rhsfile, Lnucuse, &
                  Rnucuse
      use sys, only : die
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
      integer nq                ! Number of surface q-points
! FDN q may have a z component
      real*8 q(3,nq)         ! Surface q-point
! FDN
      character*20 slabel       ! System Label (to name output files)
      character*20 llabel       ! Lead Label (to name output files)
      character*59 sname        ! System Name

      logical tjob !True if Left, False if Right
      logical tbtjob 
      integer joutfile          !unit-number of out-file

      character*33 hsfile        !name of HS-input file
      integer NEn ! no. contour points
      complex*16 contour(NEn),wgf(NEn) !energy contours and weights for GF
      real*8 efermi             ! the Fermi-energy we REQUIRE the electrode
                                ! to have (e.g. when a voltage is applied)
      integer ispin             !spin number for which to generate gf
    
!=======================================================================

!   k_|| == q-points in simple (non-repeated) cell:
      real*8, pointer:: q1(:,:),wq1(:) ! k_|| and their weights 
! FDN
      real*8, dimension(nq) :: wq
! FDN
      real*8, pointer:: q2(:,:) ! shifted q1 points
      integer na1,na2           ! Replication of unitcell
      integer nq1                ! no. q1-points (<= na1*na2 for gamma) 
      real*8 kpoint(3)          !3D k-point (q,kz)
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
!      complex*16 HAA(NGAA,nq1)  ! Real-space bulk Hamiltonian
!      complex*16 SAA(NGAA,nq1)  ! Real-space bulk Overlap
!      complex*16 GAAq(NGAA,nq1) ! Inverse ideal GF (z*SAA - HAA - SigmaAA)

!      complex*16 GS(ngaa1)
      complex*16, allocatable :: GS(:),HAA(:,:),SAA(:,:), GAAq(:,:)

      integer nucuse             !Truncate GF to *last/first* nucuse atoms
!                                !         for   *left/right*

      integer, dimension (:), allocatable :: lastou
!      integer lastou(0:NG1) ! Index of last orbital of each *used* atom
! (we use that the number of atoms < NG2)
      integer nou               !no. used orbitals     

!      complex*16 zbulkdos(NEn)
      complex*16 zbulk 

!=======================================================================
!     Helpers, workspace, tempos etc...

      character*70 gftitle      !title of gf
      character*33 paramfile    !parameter file with this "tag"

      logical tinit,tlast

!     LEFT/RIGHT sign in GF calc.: Exp(i*LRsign*Kz*Z):       
      integer LRSign

      character*5 tag
      character stag
      character*6  gfjob
      character*33 gffile_default

      logical exist1
      logical tdos,tkham
      
      complex*16, dimension (:), pointer:: H00, S00, H01, s01

!      complex*16 h00(ngaa1),s00(ngaa1)
!      complex*16 h01(ngaa1),s01(ngaa1)
!      complex*16 ab(ngaa1),ba(ngaa1)
       complex*16, allocatable :: ab(:),ba(:)


      integer i,l1,l2,ia,ia2
      integer iq1,iq
      integer ju,ng1tmp
      integer ierror

      integer ngaa,ngaa1
      integer NG1      ! Number of basis  orbitals
      integer NG2      ! Number of orbitals used

      real*8  factor
      real*8, allocatable:: eig(:)
      

      complex*16 ZEnergy
      complex*16 ZSEnergy
      complex*16 zdos
      integer iEn

! FDN 
      real*8 cell(3,3)
      integer kscell(3,3),j
      real*8 kdispl(3),tol
      complex*16 zenergy_fdf
      integer ispin_fdf,iEn_fdf
      real*8 k_fdf(3)
! FDN

      external io_assign,io_close

! FDF-stuff:
      character*33 paste,header,itemfdf     
      real*8 fdf_convfac
      external paste


!=======================================================================
! BEGIN:
!=======================================================================
       tdos =.true.


       factor =fdf_convfac('Ry','eV')


      if(tjob)then
!!   left
       gfjob='LEFT  '
       LRSign=1
       tag='Left'
       stag='L'
      else
!!   right
       gfjob='RIGHT '
       LRSign=-1
       tag='Right'
       stag='R'
      endif

! output .GF files to this name

! FDN
      if ( LFrstTime .or. RFrstTime ) then 
! FDN

      if (LJob) hsfile = Lhsfile
      if (RJob) hsfile = Rhsfile 

!=======================================================================
!
! Read-in bulk parameters and do calculation:
!
!=======================================================================
        nullify(h00)
        nullify(s00)
        nullify(h01)
        nullify(s01)
        nullify(lasto)

! FDN To be changed !!!
        allocate(h00(1))
        allocate(s00(1))
        allocate(h01(1))
        allocate(s01(1))

! initialize and get nspinin         
        tinit=.true.
        tlast=.false.
        tkham=.false.
! Before read, dimension of H00,s00,h01 and s01 = 0
        ng1=0
        call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin, &
          hsfile, nua,lasto,ng1,nspin,cell,kscell,kdispl, &
          H00,s00,h01,s01) ! ->

! FDN To be changed !!!
        deallocate(h00)
        deallocate(s00)
        deallocate(h01)
        deallocate(s01)


        if (LJob) then
          Lng1=ng1

          ngaa1=ng1*ng1
          allocate(LGS(ngaa1))

          itemfdf = paste('TS.NumUsedAtoms',tag) 
!          nucuse = fdf_integer(itemfdf,nua) !default use all atoms in uc.
          if (Lnucuse == 0) then
             nucuse=nua
          else 
             nucuse=Lnucuse
          end if

          nuaL=Lnucuse
     
!            write(*,*) "Used no. atoms:",nucuse,nua
!     Truncate lasto to used atoms: lastou:

          allocate(lastou(0:ng1))
          lastou=0

          lastou(0)=0
          ia2=0   
          if(LRsign .EQ. 1) then !Left
!     use nucuse *last* atoms            
              do ia=nua-nucuse ,nua
                 lastou(ia2)=lasto(ia)
                 ia2=ia2+1
              end do
           else                !Right
!     use nucuse *first* atoms            
              do ia=0,nucuse
                 lastou(ia2)=lasto(ia)
                 ia2=ia2+1               
              end do
           end if              !L or R            
         
!     No. used orbitals in uc:
           nou=0
           do ia=1,nucuse
              nou=nou + (lastou(ia)-lastou(ia-1))
           end do

           NGL2 = nou
           Lngaa=NGL2*NGL2
           deallocate(lastou)

!           do iEn=1,NEn
!              zbulkdos(iEn)=dcmplx(0d0,0d0)
!           end do
           

! FDN Changing the dimensions
           allocate(LHAA(LNGAA,1))  
           allocate(LSAA(LNGAA,1))
           allocate(LGAAq(LNGAA,1))

        end if
        if (RJob) then
          Rng1=ng1


          ngaa1=ng1*ng1
          allocate(RGS(ngaa1))
          itemfdf = paste('TS.NumUsedAtoms',tag) 
!          nucuse = fdf_integer(itemfdf,nua) !default use all atoms in uc.
          if (Rnucuse == 0) then
             nucuse = nua 
          else
             nucuse=Rnucuse
          end if

          nuaR=Rnucuse

!            write(*,*) "Used no. atoms:",nucuse,nua
!     Truncate lasto to used atoms: lastou:

          allocate(lastou(0:ng1))
          lastou=0

          lastou(0)=0
          ia2=0
          if(LRsign .EQ. 1) then !Left
!     use nucuse *last* atoms            
             do ia=nua-nucuse ,nua
                lastou(ia2)=lasto(ia)
                ia2=ia2+1
             end do
          else                !Right
!     use nucuse *first* atoms            
             do ia=0,nucuse
                lastou(ia2)=lasto(ia)
                ia2=ia2+1               
             end do
          end if              !L or R            
         
!     No. used orbitals in uc:
          nou=0
          do ia=1,nucuse
             nou=nou + (lastou(ia)-lastou(ia-1))
          end do

          NGR2 = nou
          Rngaa=NGR2*NGR2       
          deallocate(lastou)

!          do iEn=1,NEn
!             zbulkdos(iEn)=dcmplx(0d0,0d0)
!          end do
! FDN Changing the dimensions
          allocate(RHAA(RNGAA,1))  
          allocate(RSAA(RNGAA,1))
          allocate(RGAAq(RNGAA,1)) 

        end if 


        if(sppol.AND.nspin<2) write(joutfile,*) & 
          "WARNING/ERROR: Spin-polarized electrodes was expected"
        if(.not.sppol .AND. nspin>1)  write(joutfile,*) &
          "WARNING/ERROR: Spin-polarized electrodes was NOT expected"

      end if ! LFirstTime .or. RFirstTime

! FDN Added ....

!      if (.not. WGFFiles ) then
         iEn=iEn_fdf

         if (LJob .and. (.not. LFrstTime)) then

            if ( iEn == 1 ) then


              tinit=.false.
              tkham=.false.
              tlast=.false.
              call sethhm2(joutfile,tinit,tkham,tlast,k_fdf, &
                   ispin_fdf,hsfile,nuaL,Llasto,Lng1,nspin, &
                   cell,kscell,kdispl,LH00,Ls00,Lh01,Ls01) 
            end if ! iEn == 1
            if(iEn.eq.1) then
                  
               i=0
               do l2 = 1,NGL2
                  do l1 = 1,NGL2
                     i=i+1
! FDN iq1 ==> 1
          LHAA(i,1) = LH00(l1+(LNG1-NGL2)+LNG1*(l2+(LNG1-NGL2)-1))
          LSAA(i,1) = LS00(l1+(LNG1-NGL2)+LNG1*(l2+(LNG1-NGL2)-1))
! FDN
                  end do  ! l1
               end do     ! l2
              
               NGAA=NGL2*NGL2
               call zaxpy(NGAA,dcmplx(efermi,0.d0), &
                      LSAA(1,1),1,LHAA(1,1),1)

            end if        ! iEn.eq.1

            zsenergy = zenergy_fdf-efermi
            call calc_green(Lng1,zsenergy,Lh00,Ls00,Lh01,Ls01, &
                   LGS,zdos,joutfile,tjob,tdos)
! A Fazer .... colocar wk ...           
            if(tdos) zbulk =  zdos


            i=0
                  do l2 = 1,NGL2
                     do l1 = 1,NGL2
                        i=i+1
! FDN iq1 ==> 1
           LGAAq(i,1) = Lgs(l1+(LNG1-NGL2)+ LNG1*(l2+(LNG1-NGL2)-1) )
! FDN
                     end do     ! l1
                  end do        ! l2

         end if ! LJob


         if (RJob .and. .not.RFrstTime) then

            if ( iEn == 1 ) then



              tinit=.false.
              tkham=.false.
              tlast=.false.
              call sethhm2(joutfile,tinit,tkham,tlast,k_fdf, &
                   ispin_fdf,hsfile,nuaR,Rlasto,Rng1,nspin, &
                   cell,kscell,kdispl,RH00,Rs00,Rh01,Rs01) 
            end if ! iEn
 
            if(iEn.eq.1) then
                 
                     i=0
                     do l2 = 1,NGR2
                        do l1 = 1,NGR2
                           i=i+1
! FDN iq1 ==> 1
               RHAA(i,1) = RH00(l1+RNG1*(l2-1))
               RSAA(i,1) = RS00(l1+RNG1*(l2-1))
! FDN
                        end do  ! l1
                     end do     ! l2                  
         
!
!     Shift so efermi is energy zero
!
! FDN iq1 ==> 1
                  NGAA=NGR2*NGR2
                  call zaxpy(NGAA,dcmplx(efermi,0.d0), &
                      RSAA(1,1),1,RHAA(1,1),1)
! FDN

             end if ! iEn.eq.1 

             zsenergy = zenergy_fdf-efermi


             call calc_green(Rng1,zsenergy,Rh00,Rs00,Rh01,Rs01, &
                   RGS,zdos,joutfile,tjob,tdos) 
! A Fazer .... wk
            if(tdos) zbulk =   zdos
             i=0
                  do l2 = 1,NGR2
                     do l1 = 1,NGR2
                        i=i+1
! FDN iq1 ==> 1
                        RGAAq(i,1) = Rgs(l1+RNG1*(l2-1))
! FDN
                     end do     ! l1
                  end do        ! l2           

         end if !RJob

!      end if ! .not. WGFFiles

! FDN
      if (LFrstTime .and. tjob ) then
        LFrstTime=.false.     
! Alcocar LH00, LS00, etc ..  
      end if
      if (RFrstTime .and. (.not.tjob) ) then
        RFrstTime=.false.
! Alocar RH00, RS00, etc ...
      end if 
! FDN     


! ======================================================================
      return
      end subroutine green
! ======================================================================




      subroutine sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin, &
          hsfile,nua,lasto,nuo,nspin,cell,kscell,kdispl, &
          Hk,Sk,Hk2,Sk2)    ! ->

! ##################################################################
! ##               Setup Hamiltonian in k-space                   ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##################################################################
!
   

      use tsread2
! FDN
      use parallel, only : IOnode
      use m_tbt_kpts, only : reclat
! FDN
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

!=======================================================================
      real*8 EPS
      parameter(EPS=1.0d-8)
!=======================================================================


! INPUT

      integer joutfile          !info-out file unit
      logical tkham             ! true if full H_k hamiltonian is generated

!      integer nuo               !No. states in unitcell (expected in read-in)
      real*8 kpoint(3) 
      integer ispin             !set up H for ispin
      character*33 hsfile       !H,S parameter file


!-----------------------------------------------------------------------
! READ-IN from HS-file/OUTPUT


      integer nua               ! No. atoms in unitcell
      integer mo             ! Number of orbitals in supercell
      integer nuo            ! Number of basis  orbitals
      integer mno             ! Number of orbitals interacting
      integer nspin          ! Spin polarization (1 or 2)

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
      real*8  qtot              ! Total number of electrons
      real*8  temp              ! Electronic temperature for Fermi smearing

      logical Gamma             ! true if Gamma

!      real*8  xij(3,maxno,maxuo) ! Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point ..
!                                but must be written here!!)

      real*8 cell(3,3)          ! unit cell
! FDN
      integer kscell(3,3)
      real*8 kdispl(3)
      
! FDN

!      integer nua               ! No. atoms in unitcell
!      integer isa(maxua)         ! Species index of each atom
!      real*8 xa(3,maxua)         ! Atomic coordinates (Bohr)

    

      integer, dimension (:), pointer:: listh, listhptr, &
                             numh,indxuo,isa,lasto
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs


!-----------------------------------------------------------------------
! OUTPUT
      complex*16 Hk(:), Sk(:)
      complex*16 Hk2(:), Sk2(:)


!      integer, dimension (:), pointer:: lasto
!      integer lasto(0:maxua)   ! Index of last orbital of each atom
!-----------------------------------------------------------------------
! Helpers
      integer nb           ! No. basis orbitals on atoms (the same on all..)
!      real*8 xo(3,maxuo)        ! Atomic coordinates (Bohr)
      real*8, allocatable ::   xo(:,:) 

      integer nuotot,notot,maxnh
      integer ia
      integer i,j,io,jo,iuo,juo,j2
      real*8 k(3),kxij,rcell(3,3)
      real*8 recell(3,3)
      complex*16 cphase

!      integer ix(maxnh)
      integer, dimension (:), pointer :: ix
      integer icoi,icoa
      integer iprop,inn,it,in,ind
      real*8   xc      
      logical tinit,tlast
!-----------------------------------------------------------------------
! SAVED:
!      save numh,listh,indxuo,Gamma
!      save H,S,efs,rcell,xij,listhptr,nuotot
!      save ix

!=======================================================================
! BEGIN
      
       iprop = 3
       icoa = 0
       icoi = 0

       if (.not. LFrstTime .and. LJob) then


 
        H=>LH
        if (.not. Lgamma) xij=>Lxij
        S=>LS
        efs=>Lefs
        listh=>Llisth
        listhptr=>Llisthptr
        numh=>Lnumh
        indxuo=>Lindxuo
!        allocate(ix(size(Lix)))
        ix=>Lix
        nuotot=Lnuotot
        rcell=Lrcell
        gamma=Lgamma


       end if
       if (.not. RFrstTime .and. RJob) then


 
        H=>RH
        if (.not. Rgamma) xij=>Rxij
        S=>RS
        efs=>Refs
        listh=>Rlisth
        listhptr=>Rlisthptr
        numh=>Rnumh
        indxuo=>Rindxuo
!        allocate(ix(size(Rix)))
        ix=>Rix
        nuotot=Rnuotot
        rcell=Rrcell
        gamma=Rgamma

       end if

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
              nullify(efs)
              nullify(xa)


! FDN kscell and kdispl added as dummys

      call TSiohs('read', &
      hsfile, gamma, nua, nuotot,notot,nspin, &
                   maxnh,numh, listhptr, listh, H, S, qtot, temp, &
                   xij, indxuo, efs, cell, isa, lasto, xa, &
                   kscell, kdispl)

! FDN

           nuo = nuotot


       
         if (LFrstTime .and. LJob) then
            allocate(LH00(nuo*nuo))
            allocate(LS00(nuo*nuo))
            allocate(LH01(nuo*nuo))
            allocate(LS01(nuo*nuo)) 
! Alocar matrizes que sao salvas
            call alloc_gf_vars(LH,LS,Lxij,Lindxuo,Llisth,Llisthptr,Lnumh &
      ,Lefs,Lix,notot,nuotot,maxnh,nspin,gamma)
         end if         

         if (RFrstTime .and. RJob) then
            allocate(RH00(nuo*nuo))
            allocate(RS00(nuo*nuo))
            allocate(RH01(nuo*nuo))
            allocate(RS01(nuo*nuo))
! Alocar matrizes que sao salvas 
            call alloc_gf_vars(RH,RS,Rxij,Rindxuo,Rlisth,Rlisthptr,Rnumh &
      ,Refs,Rix,notot,nuotot,maxnh,nspin,gamma)     
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
            call cp_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs,ix, &
      rcell,nuotot, gamma)
!          LFrstTime=.false.
        end if         

        if (RFrstTime .and. RJob ) then
            call cp_gf_vars(H,S,xij,indxuo,listh,listhptr,numh,efs,ix, &
      rcell,nuotot, gamma)
!         RFrstTime=.false.          
        end if


         tinit = .false.
         return
!-----------------------------------------------------------------
         end if                    !tinit
!-----------------------------------------------------------------
         if(tlast) then

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

!
!     Symmetrize and *make EF the energy-zero*!!!
!
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))

           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) ) &
                - efs(ispin)*Sk(it)
           Hk(in) =  dconjg(Hk(it))

        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in) = Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
         Hk(in) = Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in)) &
                - efs(ispin)*Sk(in)
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

!
!     Symmetrize and *make EF the energy-zero*!!!
!
      do iuo = 1,nuo
         do juo = 1,iuo-1
           it = juo+(iuo-1)*nuo
           in = iuo+(juo-1)*nuo

           Sk(it) = 0.5d0*( Sk(it) + dconjg(Sk(in)) )
           Sk(in) =  dconjg(Sk(it))

           Hk(it) = 0.5d0*( Hk(it) + dconjg(Hk(in)) ) &
                - efs(ispin)*Sk(it)
           Hk(in) =  dconjg(Hk(it))

        enddo

         in = iuo+(iuo-1)*nuo
         Sk(in)=Sk(in) - dcmplx(0d0,1d0)*dimag(Sk(in))
         Hk(in)=Hk(in) - dcmplx(0d0,1d0)*dimag(Hk(in)) &
                - efs(ispin)*Sk(in)

      enddo
      

      do iuo = 1,nuo
         do juo = 1,nuo
          in = iuo+(juo-1)*nuo
             Hk2(in)=Hk2(in) - efs(ispin)*Sk2(in)
         enddo
      enddo    
        
!-----------------------------------------------------------------
       endif                      !tkham
!-----------------------------------------------------------------
       
! ===============================================================
       return
       END subroutine sethhm2
! ===============================================================














end module m_tbt_gf
