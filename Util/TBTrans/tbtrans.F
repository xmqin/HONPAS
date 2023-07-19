      PROGRAM tbtrans 
C
C ##################################################################
C ##                           The                                ##
C ##     Localized Basis-set Nonequilibrium Greens function       ## 
C ##                      Transport Program                       ##
C ##                                                              ## 
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##           MIC, Technical University of Denmark (DTU)         ##
C ##                                                              ##
C ##################################################################
C
C
C Tight-binding transport program
C Copyright by Mads Brandbyge, 1999-2005.
C The use of this program is allowed for non-for-profit research only.
C Copy or disemination of all or part of this package is not
C permitted without prior and explicit authorization by the author.
C

C------------------------------------------------------------------
C
C    Modified by F. D. Novaes 
C
C------------------------------------------------------------------
 
      use hamil2
      use fdf
C Module Routines 
      use m_tbt_options, only : tbt_read_options
C Module Variables
      use m_tbt_options, only : Emin, Emax, Volt, ncontour, GFeta, 
     .     USEBULK, neigch, sppol,
     .     RNode, NBUFATL, NBUFATR, hsfile,isoat1,isoat2,
     .     NA1L,NA2L,nqL,NA1R,NA2R,nqR,CalcIeig  
      use files, only : slabel
      use parallel, only : Node, Nodes, IOnode
      use sys, only : die
      use m_tbt_kpts, only : get_kp_on_node
      use m_fdf_global, only : fdf_global_get
C Module Routines and Variables
      use m_tbt_gf
#ifdef MPI
      use mpi_siesta
      use m_mpi_utils
#endif MPI

c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------

      logical PRINTALOT
      parameter(PRINTALOT=.false.)

c======================================================================

      real*8 eV
      parameter ( eV = 1.d0 / 13.60580d0 )
      real*8 Pi
      parameter(Pi=3.141592653589793238462643383279503d0)

      integer ns,maxuo
c======================================================================
c READ-IN: 
c   Hamiltonian and Overlap matrix elements for given k:

c     contact part of Ham./Ovl.             
       

c      complex*16, dimension(:,:), allocatable ::
c     . Hk,Sk  
c     we will only accept maxuo = nou presently

      integer nua               !no. atoms in unit cell 
c     we will only accept maxua = nua presently

      integer nuo               !no. orbs. in unit cell
      integer nou               !no. of orbs. from readin (nuo) 
c                               !minus buffer-orbitals
      integer nspin             !nspin=1 or nspin=2(spin polarized)
c======================================================================
c FROM FDF
      real*8, pointer:: kxy(:,:),wkxy(:) ! k_|| points and their weights 
#ifdef MPI
      real*8, pointer :: kxybuff(:,:),wkxybuff(:)
#endif MPI
      integer nkxy                ! Number of surface kxy-points
#ifdef MPI
      integer nkxybuff
      integer iNode
#endif


      integer NBUFL,NBUFR       ! No. buffer states, L/R



      integer kscell(3,3)
      real*8  kdispl(3)

      logical genlogic
      integer  iukp,lkp,jkp




c ----------------------------------------------------------------
c     Greens functions:
c ----------------------------------------------------------------
c     Energy Contour:
c

      complex*16, dimension(:), allocatable ::
     . contour, WGF ! Equilibrium contour

      real*8, dimension(:), allocatable ::  bulkdosL, bulkdosR ! bulk DOS
#ifdef MPI
      real*8, dimension(:), allocatable ::  bulkdosbuffL,bulkdosbuffR
#endif

      complex*16 :: zbulk
      
       complex*16, dimension(:,:), allocatable ::
     . SFEL,GAMMAL,SFER,GAMMAR
      complex*16, dimension(:,:,:), allocatable ::
     . HAAL,SAAL, HAAR,SAAR  


      complex*16, dimension(:), allocatable ::  H,S

      real*8 , dimension(:), allocatable :: teig !eigentransmissions
      real*8 , dimension(:), allocatable :: tav,pdosav,tdosav !kxy-point averages
#ifdef MPI
      real*8 , dimension(:,:,:), allocatable :: teigbuff
      real*8 , dimension(:), allocatable :: tavbuff,pdosavbuff,
     .                                      tdosavbuff 
#endif MPI

c
c     Read-in GF(q)'s
c
C      real*8, dimension (:,:), pointer:: qL,qR
C      real*8, dimension (:), pointer:: wqL,wqR

       integer,dimension(:), allocatable ::
     . lastoL,! L uc. atoms (from read-in superc)
     . lastoR! R uc. atoms (from read-in superc)
      
      integer, dimension (:), pointer:: lasto
      complex*16, dimension (:), pointer:: Hk, Sk
      complex*16, dimension (:), allocatable:: Hiso, Siso, aux
      real*8, dimension (:), allocatable:: aux2



c ----------------------------------------------------------------
      integer joutfile          ! unit-number of out-file
      parameter(joutfile=6)
c======================================================================
c OUTPUT:
      real*8 tottrans           !total transmission
      complex*16, dimension(:), allocatable :: tt !transmission matrix
      real*8 totdos             ! DOS
      real*8 pdos               !PDOS
#ifdef MPI
      real*8, allocatable, dimension(:,:,:) :: prtbuff
#endif MPI

c ==================================================================
c     Left/Right inverse Greens functions:
      integer NGL,NGR           !no. orbs. in L/R electrode parts
c ==================================================================

c Fermi levels
      real*8 EFermiL,EFermiR,EFermi0

      real*8 Current,f_L,f_R,fermi,Temp,dE
#ifdef MPI
      real*8 Currentbuff
#endif MPI

c No. contour points (=NCONTOUR)
      integer NEn

c ==================================================================
c Helpers, workspace, tempos etc...

      integer jgfL,jgfR         !File units for Left,Right GF

      integer isobegin,isoend,ncp,isoTDOSBegin,isoTDOSEnd
      
      integer mc                !max no. contact states 
      integer nc                !no. contact states
      
      complex*16, dimension(:), allocatable :: Sc
      complex*16, dimension(:), allocatable :: GF, GFRGF
      complex*16 ZEnergy
      
      real*8  kpoint(3)

      integer i,j,io,jo,io1,jo1,ia,ia2,ispin,ikxy
      integer iEn
      integer unit1,unit2,unit3,unit4,unit5,unit6
      logical tinit,FirstTime
      
      
c Externals:
C      External io_assign,io_close
      complex*16 zdotc
      External zdotc            ! BLAS function

c Debug:
      
#ifdef MPI
      real*8, allocatable :: eigbuff(:,:)
      integer, allocatable :: ncpbuff(:) 
#endif
      real*8, allocatable :: eig(:)        ! Eigenvalues of cluster




c FDF-stuff:
      logical overflow          ! True = Some of the dimensions is too small
      character*75 fname
     
     
 
      character paste*33

      external io_assign, io_close, paste
  


#ifdef MPI
      integer :: MPIerror, Status(MPI_Status_Size),SendNode
#endif MPI




c=======================================================================
c BEGIN:
c=======================================================================

#ifdef MPI
      call MPI_Init( MPIerror )
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror ) 
#endif MPI
      IOnode = (Node == 0) 
      RNode = 0

      overflow = .false.
      FirstTime=.true.

C Initialize fdf   
      call reinit(slabel)

C Read Input values
      call tbt_read_options()


c        
c     We take the mean Fermi energy of the electrodes to be zero!
c
      EFermiL =  Volt*eV/2.0d0 !Ry
      EFermiR = -Volt*eV/2.0d0 !Ry
      EFermi0 =       0.0d0
      Temp = 0.0d0
c
c     Setup contour for Transmission 
c      
      allocate(contour(ncontour))
      allocate(wgf(ncontour))
      call mkRealContour(joutfile,Emin,Emax,GFeta,NCONTOUR,contour,WGF)
      allocate(bulkdosL(ncontour))
      allocate(bulkdosR(ncontour))


      NEn=ncontour
      dE = (Emax-Emin)/NEn
      
      nspin=1
      if(sppol) then 
         nspin=2
      end if      


      if (IOnode) then
        genlogic = fdf_block('kgrid_Monkhorst_Pack',iukp)
        if ( genlogic ) then
         do lkp = 1,3
           read(iukp,*) (kscell(jkp,lkp),jkp=1,3), kdispl(lkp)
         enddo
C       close(iukp)
        else
          write(*,*) 'Must specify k-points with Monkhorst_Pack'
          write(*,*) 'Stoping code'  
          stop
        end if
      end if ! IOnode

#ifdef MPI
      call broadcast(kscell(:,:))
      call broadcast(kdispl(:))
#endif MPI
     
      fname=paste(slabel,'.TSHS') 

      call get_kp_on_node(fname,nkxy,kxy,wkxy,kscell,kdispl)
     


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ispin=1,nspin ! Loop over spins,
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         if (IOnode) then
c
c     Open output files
c

            call io_assign(unit1)
            fname=paste(slabel,'.TRANS' )
            if(nspin.gt.1) then
               if(ispin.eq.1) fname=paste(slabel,'.UP.TRANS' )
               if(ispin.eq.2) fname=paste(slabel,'.DOWN.TRANS' )
            end if
            open (unit1, file=fname,
     &        status='unknown', blank='zero', err=1314)

            call io_assign(unit4)
            fname=paste(slabel,'.AVTRANS' )
            if(nspin.gt.1) then
               if(ispin.eq.1) fname=paste(slabel,'.UP.AVTRANS' )
               if(ispin.eq.2) fname=paste(slabel,'.DOWN.AVTRANS' )
            end if
            open (unit4, file=fname,
     &        status='unknown', blank='zero', err=1314)


            call io_assign(unit5)
            fname=paste(slabel,'.LDOS' )
            if(nspin.gt.1) then
               if(ispin.eq.1) fname=paste(slabel,'.UP.LDOS' )
               if(ispin.eq.2) fname=paste(slabel,'.DOWN.LDOS' )
            end if
            open (unit5, file=fname,
     &        status='unknown', blank='zero', err=1314)


            call io_assign(unit6)
            fname=paste(slabel,'.RDOS' )
            if(nspin.gt.1) then
               if(ispin.eq.1) fname=paste(slabel,'.UP.RDOS' )
               if(ispin.eq.2) fname=paste(slabel,'.DOWN.RDOS' )
            end if
            open (unit6, file=fname,
     &        status='unknown', blank='zero', err=1314)


            if(neigch.gt.0) then
               call io_assign(unit2)
               fname=paste(slabel,'.TEIG' )
               if(nspin.gt.1) then
                  if(ispin.eq.1) fname=paste(slabel,'.UP.TEIG' )
                  if(ispin.eq.2) fname=paste(slabel,'.DOWN.TEIG')
               end if
               open (unit2, file=fname,
     &           status='unknown', blank='zero', err=1314)
            end if


            if ( CalcIeig ) then
               call io_assign(unit3)
               fname=paste(slabel,'.IEIG' )
               if(nspin.gt.1) then
                  if(ispin.eq.1) fname=paste(slabel,'.UP.IEIG' )
                  if(ispin.eq.2) fname=paste(slabel,'.DOWN.IEIG')
               end if
               open (unit3,
     &         file=fname,
     &         status='unknown', blank='zero', err=1314)
            end if ! CalcIeig

         end if ! IOnode 

 
c kxy-averaged transmission
         allocate(tav(NEn))
         allocate(pdosav(NEn))
         allocate(tdosav(NEn))
         tav=0d0
         pdos=0d0
         tdosav=0d0

c==========================================================================
c     H,S INITIALIZATION BEGIN
c
c     Read-in Hamiltonian/Overlap parameters from HS-file: 
c     Firsttime!
c==========================================================================
         if (FirstTime) then

c     
c Left GAAL:
c 
    
            LJob=.true.
            RJob=.false.
            call green(joutfile,nkxy,kxy,wkxy,
     &        NEn,contour,WGF,EFermiL,zbulk,.true.,.true.,
     & ispin,kpoint,zenergy,iEn)

            LJob=.false.
            RJob=.true.
            call green(joutfile,nkxy,kxy,wkxy,
     &        NEn,contour,WGF,EFermiR,zbulk,.false.,.true.,
     &  ispin, kpoint, zenergy,iEn)
           

            allocate(HAAL(NGL2,NGL2,nqL))
            allocate(SAAL(NGL2,NGL2,nqL))
         
            NGL=NGL2*NA1L*NA2L
         
            allocate(SFEL(NGL,NGL))
            allocate(lastoL(0:nuaL)) 


            nqR=1
            allocate(HAAR(NGR2,NGR2,nqR))
            allocate(SAAR(NGR2,NGR2,nqR))
            NGR=NGR2*NA1R*NA2R
            allocate(SFER(NGR,NGR))
            allocate(lastoR(0:nuaR))


            tinit = .true.        
            call setupkham2(joutfile,tinit,kpoint,ispin,hsfile,
     .       lasto,Hk,Sk,nua,nuo,ns) ! ->
       
            maxuo = nuo

       
c
c Deal with buffer atoms 
c
      
c the first NBUFATL atoms will be removed
            NBUFL=0
            do ia=1,NBUFATL
               NBUFL=NBUFL+(lasto(ia)-lasto(ia-1))
            end do                    !ia
c the last NBUFATR atoms will be removed
            NBUFR=0
            do ia=nua - NBUFATR+1,nua
               NBUFR=NBUFR+(lasto(ia)-lasto(ia-1))
            end do                    !ia

c
c No. states minus buffers
c     
            nou = nuo - (NBUFL+NBUFR)


c
c expected no. states on Electrode atoms:
c
            lastoL(0)=0
            ia2=0
            do ia=NBUFATL+1,NBUFATL + nuaL*NA1L*NA2L, NA1L*NA2L
               ia2=ia2+1
               lastoL(ia2)=lastoL(ia2-1) + (lasto(ia) - lasto(ia-1))
            end do                 !ia
            if(lastoL(nuaL) .ne. NGL2) then
#ifdef MPI
              call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI
              if (IOnode) then
                write(joutfile,*) 'ERROR: lastoL,NGL2',lastoL,NGL2
                call die('ERROR: Unexpected no. orbs. in L elec.')
              end if ! IOnode
            end if  ! lastoL(nuaL) .ne. NGL2

            lastoR(0)=0
            ia2=0
            do ia=nua-(nuaR*NA1R*NA2R + NBUFATR)+1,nua-NBUFATR, 
     . NA1R*NA2R
               ia2=ia2+1
               lastoR(ia2)=lastoR(ia2-1) + (lasto(ia) - lasto(ia-1))
            end do                 !ia
            if(lastoR(nuaR) .ne. NGR2) then
#ifdef MPI
              call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI
              if (IOnode) then
                 write(joutfile,*) 'ERROR: lastoR,NGR2',lastoR,NGR2
                 call die('ERROR: Unexpected no. orbs. in R elec.')
              end if ! IOnode
            end if

c
c No. contact states
c
            nc = nou - (NGL+NGR)
            if(nc.le.0) then
#ifdef MPI
              call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI
               write(joutfile,*) 
     &           'Transport : ERROR : nc = ',nc,'!!'
               call die('Transport : ERROR : nc <=0')
            end if


       
            Current = 0.d0 
            bulkdosL=0.d0
            bulkdosR=0.d0
 
c==========================================================================
c     H,S INITIALIZATION END
c==========================================================================
            mc  = maxuo -NGR -NGL


            isoTDOSBegin = lasto(NBUFATL+nuaL*NA1L*NA2L)+1-NBUFL-NGL
            isoTDOSEnd = lasto(nua-NBUFATR-nuaR*NA1R*NA2R)-NBUFL-NGL
c CHOOSE ATOM IN TERMS OF GLOBAL ATOM NUMERBS (i.e. INCLUDE ELECTRODES HERE)
            if ( isoat1 == 0 )
     .         isoat1 = 1 + NBUFATL + nuaL*NA1L*NA2L !default
            if ( isoat2 == 0 )
     .         isoat2 = nua - NBUFATR - nuaR*NA1R*NA2R !default
            
            isobegin = lasto(isoat1-1)+1
            isoend = lasto(isoat2)
c     write(joutfile,*) 'iso=',isoat1,isoat2,isobegin,isoend,ncp
            isobegin = isobegin - NBUFL - NGL
            isoend = isoend - NBUFL - NGL
c     write(joutfile,*) 'iso=',isoat1,isoat2,isobegin,isoend,ncp
            isobegin = max(1,isobegin)
            isoend = min(nuo,isoend)
c     Default: IN CASE OF ERRORS, USE CONTACT REGION
            if(isobegin.gt.isoend) then
               isobegin = NBUFL+NGL+1
               isoend = isobegin+nc-1
            end if
            ncp = isoend - isobegin + 1
            if (IOnode) then
               write(joutfile,*) repeat('=',62) 
               write(joutfile,'(a30,I10,a1,I10,a1)') 
     .           'Projection Region: atoms : [',isoat1,';',isoat2,']'
               write(joutfile,'(a30,I10,a1,I10,a6,I10)') 
     .    'Projection Region: states: [',isobegin,';',isoend,'] Tot:'
     .    ,ncp
            end if ! IOnode

 
       
            if(isobegin.gt.isoend) then
              isobegin = 1
              isoend = nc
            end if
            ncp = isoend - isobegin + 1

            FirstTime=.false.
         end if ! FirstTime
         
         if ( IOnode ) then
            write(*,*) 
            write(*,'(a,I4)') '  Starting K-Point loop for spin:',ispin
            write(*,*) 
         end if ! IOnode
c============================================================
         do ikxy = 1,nkxy           ! kxy-point grid
c============================================================

         kpoint(1)=kxy(1,ikxy)
         kpoint(2)=kxy(2,ikxy)
         kpoint(3)=kxy(3,ikxy)


        call setupkham2(joutfile,tinit,kpoint,ispin,hsfile,
     .  lasto,Hk,Sk,nua,nuo,ns)    ! ->



c
c     Remove buffer states
c
        allocate(H(maxuo*maxuo))
        allocate(S(maxuo*maxuo))


        do io1=1,nou
            io=io1+NBUFL
            do jo1=1,nou
              jo=jo1+NBUFL
              H(jo1 + nou*(io1-1))=Hk(jo+nuo*(io-1))
              S(jo1 + nou*(io1-1))=Sk(jo+nuo*(io-1))
           end do                 !jo1
        end do                 !io1




c
c     Make sure Left-Right part is zero:
c
        do io = nou-NGR+1,nou
           do jo = 1,NGL
              H(jo + nou*(io-1))=dcmplx(0d0,0d0)
              S(jo + nou*(io-1))=dcmplx(0d0,0d0)            
              H(io + nou*(jo-1))=dcmplx(0d0,0d0)
              S(io + nou*(jo-1))=dcmplx(0d0,0d0)
           end do              !io
        end do                 !jo

c
c     Contact Region Overlap
c
       allocate(Sc(mc*mc))
         do io=1,nc
            io1 = io + NGL
            do jo=1,nc
               jo1 = jo + NGL
               Sc(jo+nc*(io-1)) = S(jo1+nou*(io1-1))
            end do              !jo
         end do                 !io


       allocate(GF(mc*mc))
       allocate(GFRGF(mc*mc))
       allocate(tt(nc*nc))


c ============================================================================
c Eigenstates of isolated region:

c
c     Eigenvalues of the isolated PDOS region:
c
       if ( CalcIEig ) then
         allocate(Hiso(ncp*ncp))
         allocate(Siso(ncp*ncp))
         do j=1,ncp
            jo=j+isobegin-1+NGL
            do i=1,ncp
               io=i+isobegin-1+NGL
               Hiso(i+ncp*(j-1)) = H(io+nou*(jo-1))
               Siso(i+ncp*(j-1)) = S(io+nou*(jo-1))
            end do                !i
         end do                   !j


         allocate(eig(maxuo))
         allocate(aux(ncp*ncp))
         allocate(aux2(4*ncp))

                
          call cdiag(Hiso,ncp,Siso,ncp,ncp,eig,aux,ncp,aux2)


#ifdef MPI
          if ( .not. allocated(eigbuff) ) allocate(eigbuff(nkxy,maxuo))
          if ( .not. allocated(ncpbuff) ) allocate(ncpbuff(nkxy))
          if (size(eigbuff) /= nkxy*maxuo ) then
             deallocate(eigbuff)
             if ( ikxy > 1 ) then
                call die('Problems with eigbuff dimension !!')
             end if
             allocate(eigbuff(nkxy,maxuo))
          end if
          eigbuff(ikxy,1:maxuo)=eig(1:maxuo)
          ncpbuff(ikxy)=ncp
#else 
         write(unit3,*) '# Kxy-point: ',ikxy
         write(unit3,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
         do i=1,ncp
            write(unit3,*) i,(eig(i))/eV
         end do
         if(ikxy.eq.nkxy) call io_close(unit3)
#endif MPI
     
       
         deallocate(eig)
         deallocate(Hiso) 
         deallocate(Siso)
         deallocate(aux)
         deallocate(aux2)
       end if ! CalcIEig

       
c ============================================================================

#ifndef MPI
       write(joutfile,'(a22,1X,i2,1X,a16,i3)') 
     .      'Transport for spin ->',ispin,', kxy-point -> ', ikxy

       write(joutfile,'(a16,f5.3,1X,f5.3,1X,f5.3,a6,f6.4)') 
     . '         kxy=(',
     .      kxy(1,ikxy),kxy(2,ikxy),kxy(3,ikxy), '),w= ' ,wkxy(ikxy)
       write(joutfile,*) 
     .   '   E           T(E)            TotDOS         PDOS'
#endif MPI

C Headers for output file; for each k-point
       if ( IOnode ) then
          write(unit1,*) '# Kxy-point: ',ikxy
          write(unit1,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)

          if(neigch.gt.0 .and. neigch.le.nc) then
             write(unit2,*) '# Kxy-point: ',ikxy
             write(unit2,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
          end if
       endif  ! IOnode

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
      do 7913 iEn = 1, NCONTOUR      ! loop over energy points
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
    
         ZEnergy=contour(iEn)

c----------------------------------------------------------------------
c At this point we need the transformation from
c G(q) in (1x1) to G(0) in (NA1xNA2)
c
c SFEL/SFER  <------ GAAL/GAAR
c----------------------------------------------------------------------
c LEFT:  
         allocate(GAMMAL(NGL,NGL))        


 
         LJob=.true.
         RJob=.false. 
         call green(joutfile,nkxy,kxy,wkxy,
     &        NEn,contour,WGF,EFermiL,zbulk,.true.,.true.,
     & ispin,kpoint,zenergy,iEn)     
    
         call getsfe(USEBULK,joutfile,jgfL,
     &        HAAL,SAAL,ZEnergy,nqL,
     &        ngL2,nuaL,lastoL,NA1L,NA2L,
     &        NGL,SFEL,GAMMAL)

         deallocate(GAMMAL) 


         bulkdosL(iEn)=bulkdosL(iEn)-(dimag(zbulk)/Pi)*wkxy(ikxy)


c debug
c         if(iEn .eq.1) write(*,*) 
c     .        'hola: tbt H(1,28)=',H(1+nou*(28-1))


c RIGHT:  
         allocate(GAMMAR(NGR,NGR))

         LJob=.false.
         RJob=.true.
         call green(joutfile,nkxy,kxy,wkxy,
     &        NEn,contour,WGF,EFermiR,zbulk,.false.,.true., 
     &  ispin, kpoint, zenergy,iEn)


         bulkdosR(iEn)=bulkdosR(iEn)-(dimag(zbulk)/Pi)*wkxy(ikxy)



         call getsfe(USEBULK,joutfile,jgfR,
     &        HAAR,SAAR,ZEnergy,nqR,
     &        ngR2,nuaR,lastoR,NA1R,NA2R,
     &        NGR,SFER,GAMMAR)


         deallocate(GAMMAR) 


c==========================================================================


         CALL transmission(USEBULK,joutfile,nou,H,S,
     &        nc,NGL,SFEL,NGR,SFER,ZEnergy,
     &        GF,GFRGF,tottrans,tt)
c on return GFRGF contains GF^dagger.GammaR.GF                 
         

! DOS
         totdos = (-1.0d0/Pi)*DIMAG(zdotc(nc*nc,Sc,1,GF,1))

! Calculate the "excluded" DOS and then subtract from totdos
         pdos=0d0
         do j=isoTDOSBegin, isobegin-1
            do i=1,nc
               pdos = pdos -
     .         (1.0d0/Pi)*DIMAG(Sc(i+nc*(j-1))*GF(i+nc*(j-1)))
            end do              !i
         end do                 !j

         do j=isoend+1, isoTDOSEnd
            do i=1,nc
               pdos = pdos -
     .         (1.0d0/Pi)*DIMAG(Sc(i+nc*(j-1))*GF(i+nc*(j-1)))
            end do              !i
         end do
         pdos = totdos - pdos ! Subtract

         if(nspin.eq.1) then !incl. both spins
            tottrans = tottrans*2.0d0
            pdos=pdos*2.0d0 
            totdos=totdos*2.0d0
         end if


      tav(iEn)    = tav(iEn) + tottrans*wkxy(ikxy)
      tdosav(iEn) = tdosav(iEn) + totdos*wkxy(ikxy)
      pdosav(iEn) = pdosav(iEn) + pdos*wkxy(ikxy)

#ifdef MPI
        if ( .not. allocated(prtbuff) ) 
     . allocate(prtbuff(NCONTOUR,nkxy,3))
       
        if (size(prtbuff) /= NCONTOUR*nkxy*3 ) then
           deallocate(prtbuff)
           if ( ikxy > 1 ) then
              call die('Problems with prtbuff dimension !!')
           end if
           allocate(prtbuff(NCONTOUR,nkxy,3))
        end if
        prtbuff(iEN,ikxy,1)=tottrans
        prtbuff(iEN,ikxy,2)=totdos
        prtbuff(iEN,ikxy,3)=pdos 
#else
         write(joutfile,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,tottrans,totdos,pdos


         write(unit1,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,tottrans,totdos,pdos

#endif MPI

c
c     Call Eigenchannel Transmission
c
         if(neigch.gt.0 .and. neigch.le.nc) then
              allocate(teig(nc))
              teig = 0.0
            call teigchana(nc,tt,neigch,teig)
#ifdef MPI
            if (.not. allocated(teigbuff)) 
     .    allocate(teigbuff(nkxy,NCONTOUR,neigch))
           teigbuff(ikxy,iEn,1:neigch)=teig(1:neigch)
#else
           write(unit2,'(F9.5,1X,30(1X,f8.5))')
     &           DREAL(ZEnergy)/eV,(teig(i),i=1,neigch)
#endif MPI


               deallocate(teig)
         endif                  !eigenchannel analysis

         f_L =  fermi(DREAL(ZEnergy),EFermiL,Temp)
         f_R =  fermi(DREAL(ZEnergy),EFermiR,Temp)
    
         Current = Current + (f_L-f_R)*tottrans*dE*wkxy(ikxy)
         

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
 7913    continue               ! the loop over ALL Energies
C =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=

   
      deallocate(H)
      deallocate(S)
    
      deallocate(Sc)
      deallocate(GF)
      deallocate(GFRGF)
      deallocate(tt)
    


#ifndef MPI
      write(unit1,*) " "
      if(neigch.gt.0)      
     . write(unit2,*) " "
#endif MPI

#ifndef MPI
! Leave a blank line inbetween two kpoints
       write(unit1,*)
       if(neigch.gt.0) write(unit2,*)
       if ( CalcIEig ) write(unit3,*)
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end do                    !loop over surface kxy-points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      do iEn = 1,NCONTOUR

      enddo

      
#ifdef MPI
      if (IOnode) then
         if (nkxy > 0) then
            do ikxy = 1,nkxy
               if ( CalcIEig ) then
                 write(unit3,*) '# Kxy-point: ',ikxy, ' of Node', Node
                 write(unit3,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
                 do i = 1,ncpbuff(ikxy)
                    write(unit3,*) i,(eigbuff(ikxy,i))/eV
                 enddo
               end if ! CalcIeig

               write(joutfile,'(a22,1X,i2,1X,a16,i3,a9,i6)') 
     .      'Transport for spin ->',ispin,', kxy-point -> ', ikxy,
     .      ', of Node', Node
               write(joutfile,'(a16,f5.3,1X,f5.3,1X,f5.3,a6,f6.4)') 
     . '         kxy=(',
     .      kxy(1,ikxy),kxy(2,ikxy),kxy(3,ikxy), '),w= ' ,wkxy(ikxy)
               write(joutfile,*) 
     .   '   E           T(E)            TotDOS         PDOS'

               write(unit1,*) '# Kxy-point: ',ikxy, ' of Node', Node
               write(unit1,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)



               if(neigch.gt.0) then
                  write(unit2,*) '# Kxy-point: ',ikxy, ' of Node', Node
                  write(unit2,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
               end if 
               do iEn = 1,NCONTOUR
                  ZEnergy=contour(iEn)
                  write(joutfile,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,prtbuff(iEN,ikxy,1),prtbuff(iEN,ikxy,2),
     .        prtbuff(iEN,ikxy,3)

         write(unit1,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,prtbuff(iEN,ikxy,1),prtbuff(iEN,ikxy,2),
     .        prtbuff(iEN,ikxy,3)

                  if(neigch.gt.0)
     &             write(unit2,'(F9.5,1X,30(1X,f8.5))')
     &           DREAL(ZEnergy)/eV,(teigbuff(ikxy,iEn,i),i=1,neigch)

               enddo

! Leave a blank line inbetween two kpoints
               write(unit1,*)
               if(neigch.gt.0) write(unit2,*)
               if ( CalcIEig ) write(unit3,*)

            enddo
         end if ! nkxy > 0 
      end if  ! IOnode

#ifdef MPI
      call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI

      if (.not. IOnode) then
         call MPI_Send(nkxy,1,MPI_integer,
     .          RNode,1,MPI_Comm_World,MPIerror)
         if ( nkxy > 0 ) then
            call MPI_Send(kxy(1,1),nkxy*3,MPI_double_precision,
     .          RNode,1,MPI_Comm_World,MPIerror)
            call MPI_Send(wkxy(1),nkxy,MPI_double_precision,
     .          RNode,2,MPI_Comm_World,MPIerror)
            if ( CalcIeig ) then
              call MPI_Send(ncpbuff,nkxy,MPI_integer,
     .          RNode,3,MPI_Comm_World,MPIerror)
              call MPI_Send(eigbuff(1,1),nkxy*maxuo,MPI_double_precision
     .          ,RNode,4,MPI_Comm_World,MPIerror)
            end if ! CalcIeig

            call MPI_Send(prtbuff(1,1,1),NCONTOUR*nkxy*3,
     .  MPI_double_precision,RNode,5,MPI_Comm_World,MPIerror)
            call MPI_Send(teigbuff(1,1,1),NCONTOUR*nkxy*neigch,
     .  MPI_double_precision,RNode,6,MPI_Comm_World,MPIerror)
         end if ! nkxy > 0
      end if

      if (IOnode) then
         do iNode = 1,Nodes-1
            call MPI_Recv(nkxybuff,1,
     .          MPI_integer,iNode,1,MPI_Comm_World,Status,MPIerror)
            
            nullify(kxybuff,wkxybuff)
            if ( nkxybuff > 0 ) then
               allocate(kxybuff(3,nkxybuff))
               allocate(wkxybuff(nkxybuff))
               call MPI_Recv(kxybuff(1,1),nkxybuff*3,
     .MPI_double_precision,iNode,1,MPI_Comm_World,Status,MPIerror)
               call MPI_Recv(wkxybuff(1),nkxybuff,MPI_double_precision,
     .          iNode,2,MPI_Comm_World,Status,MPIerror)
          
               if ( CalcIeig ) then
                 if(allocated(ncpbuff)) deallocate(ncpbuff)
                 allocate(ncpbuff(nkxybuff))
                 call MPI_Recv(ncpbuff,nkxybuff,
     .          MPI_integer,iNode,3,MPI_Comm_World,Status,MPIerror)

                 if(allocated(eigbuff)) deallocate(eigbuff)
                 allocate(eigbuff(nkxybuff,maxuo))
                 call MPI_Recv(eigbuff(1,1),nkxybuff*maxuo,
     .  MPI_double_precision,iNode,4,MPI_Comm_World,Status,MPIerror)
               end if ! CalcIeig

               if(allocated(prtbuff)) deallocate(prtbuff)
               allocate(prtbuff(NCONTOUR,nkxybuff,3))
               call MPI_Recv(prtbuff(1,1,1),NCONTOUR*nkxybuff*3,
     .  MPI_double_precision,iNode,5,MPI_Comm_World,Status,MPIerror)

               if(allocated(teigbuff)) deallocate(teigbuff)
               allocate(teigbuff(nkxybuff,NCONTOUR,neigch))
               call MPI_Recv(teigbuff(1,1,1),NCONTOUR*nkxybuff*neigch,
     .  MPI_double_precision,iNode,6,MPI_Comm_World,Status,MPIerror)

               do ikxy = 1,nkxybuff
                  if ( CalcIEig ) then
                   write(unit3,*) '# Kxy-point: ',ikxy, ' of Node', Node
                    write(unit3,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
                    do i = 1,ncpbuff(ikxy)
                       write(unit3,*) i,(eigbuff(ikxy,i))/eV
                    enddo
                  end if ! CalcIeig


                  write(joutfile,'(a22,1X,i2,1X,a16,i3,a9,i6)') 
     .      'Transport for spin ->',ispin,', kxy-point -> ', ikxy,
     .      ', of Node', iNode
               write(joutfile,'(a16,f5.3,1X,f5.3,1X,f5.3,a6,f6.4)') 
     . '         kxy=(',
     .      kxybuff(1,ikxy),kxybuff(2,ikxy),kxybuff(3,ikxy), 
     .  '),w= ' ,wkxybuff(ikxy)
                  write(joutfile,*) 
     .   '   E           T(E)            TotDOS         PDOS'

               write(unit1,*) '# Kxy-point: ',ikxy, ' of Node', iNode
         write(unit1,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxybuff(1,ikxy),kxybuff(2,ikxy),', 
     .   w= ',wkxybuff(ikxy)  


               if(neigch.gt.0) then
                  write(unit2,*) '# Kxy-point: ',ikxy, ' of Node', Node
                  write(unit2,'(a6,f5.3,1X,f5.3,a7,f6.4)') 
     .       '# kxy=',kxy(1,ikxy),kxy(2,ikxy),', w= ',wkxy(ikxy)
               end if

                  do iEn = 1,NCONTOUR
                     ZEnergy=contour(iEn)
                     write(joutfile,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,prtbuff(iEN,ikxy,1),prtbuff(iEN,ikxy,2),
     .        prtbuff(iEN,ikxy,3)  

                     
         write(unit1,'(F9.5,3e16.8)')
     .        DREAL(ZEnergy)/eV,prtbuff(iEN,ikxy,1),prtbuff(iEN,ikxy,2),
     .        prtbuff(iEN,ikxy,3)

                     if(neigch.gt.0)
     &                write(unit2,'(F9.5,1X,30(1X,f8.5))')
     &           DREAL(ZEnergy)/eV,(teigbuff(ikxy,iEn,i),i=1,neigch)                 
                  enddo 

! Leave a blank line inbetween two kpoints
                  write(unit1,*)
                  if(neigch.gt.0) write(unit2,*)
                  if ( CalcIEig ) write(unit3,*) 

               enddo
            end if ! nkxybuff > 0

         enddo

         

      end if ! IOnode


      call MPI_AllReduce(Current,Currentbuff,1,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror)

C Convert to Ampere
      Currentbuff = Currentbuff*13.6058*38.73*1E-06

      if (IOnode) then
         write(joutfile,*) repeat('=',62)
         write(joutfile,*) 
         write(joutfile,*) 'Results:'
         write(joutfile,*) 
         write(joutfile,'(1x,a,2ES16.8)') ' Voltage, Current(A) = ',
     .                                Volt,Currentbuff
         write(joutfile,*)
         write(joutfile,*) repeat('=',62)
 
      end if

#else
C Convert to Ampere
         Current = Current*13.6058*38.73*1E-06

         write(joutfile,*) repeat('=',62)
         write(joutfile,*) 
         write(joutfile,*) 'Results:'
         write(joutfile,*)
         write(joutfile,'(1x,a,2ES16.8)') 'Voltage, Current(A) = ',
     $                                  Volt,Current
         write(joutfile,*)
         write(joutfile,*) repeat('=',62)


      if (nspin==1) bulkdosL=bulkdosL*2.d0
      if (nspin==1) bulkdosR=bulkdosR*2.d0
      do iEn = 1,NCONTOUR
         write(unit5,'(F9.5,e16.8)') DREAL(contour(iEn))/eV,
     .   bulkdosL(iEn)
         write(unit6,'(F9.5,e16.8)') DREAL(contour(iEn))/eV,
     .   bulkdosR(iEn)
      enddo
#endif


#ifdef MPI

      allocate(tavbuff(NEn))
      call MPI_AllReduce(tav,tavbuff,NEn,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror) 
      allocate(tdosavbuff(NEn))
      call MPI_AllReduce(tdosav,tdosavbuff,NEn,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror)
      allocate(pdosavbuff(NEn))
      call MPI_AllReduce(pdosav,pdosavbuff,NEn,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror)

      allocate(bulkdosbuffL(NEn))
      call MPI_AllReduce(bulkdosL,bulkdosbuffL,NEn,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror)

      allocate(bulkdosbuffR(NEn))
      call MPI_AllReduce(bulkdosR,bulkdosbuffR,NEn,MPI_double_precision,
     .            MPI_sum,MPI_Comm_World,MPIerror)
      

      if (IOnode) then

         if (nspin==1) bulkdosbuffL=bulkdosbuffL*2.d0
         if (nspin==1) bulkdosbuffR=bulkdosbuffR*2.d0

         do iEn = 1 ,NEn
            write(unit4,'(F9.5,3e16.8)') DREAL(contour(iEn))/eV,
     .   tavbuff(iEn),tdosavbuff(iEn),pdosavbuff(iEn)

            write(unit5,'(F9.5,e16.8)') DREAL(contour(iEn))/eV,
     .   bulkdosbuffL(iEn)
            write(unit6,'(F9.5,e16.8)') DREAL(contour(iEn))/eV,
     .   bulkdosbuffR(iEn)

         enddo  
      end if

#else
      do iEn=1,NEn
         write(unit4,'(F9.5,3e16.8)')  DREAL(contour(iEn))/eV,
     .        tav(iEn),tdosav(iEn),pdosav(iEn)
      end do
#endif MPI

      
      if (IOnode) call io_close(unit4)     
      deallocate(tav)
      deallocate(tdosav)
      deallocate(pdosav)
#ifdef MPI
      deallocate(tavbuff)
      deallocate(tdosavbuff)
      deallocate(pdosavbuff)
      deallocate(bulkdosbuffL)
      deallocate(bulkdosbuffR)
#endif MPI

      if (IOnode) call io_close(unit1)
      if (IOnode) then
         if(neigch.gt.0) call io_close(unit2)
      end if



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do !loop over spin
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      

      deallocate(bulkdosL)
      deallocate(bulkdosR)



#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif
1314  continue
c======================================================================
      end
c======================================================================




c--------------------------------------------------------
      function fermi(e,ef,T)
      real*8 e,ef,T,tmp,fermi
      if(T.eq.0) then
         if(e>ef) fermi = 0.
         if(e<ef) fermi = 1.
         if(e.eq.ef) fermi = .5
      else
         tmp = (e-ef)/T
         if(abs(tmp).lt.40.) then
            fermi = 1./(1.+exp(tmp)) 
         else
            if(tmp.lt.-40.) fermi = 1.
            if(tmp.gt.40.) fermi = 0.
         endif
      endif
      return
      end
c--------------------------------------------------------
      
      
