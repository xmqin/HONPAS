      MODULE m_ts_in_siesta
!
! Routines that are called from SIESTA routines
!
!=============================================================================
! CONTAINS:
!          1) transiesta
!          2) TSVHFix
!          3) TSVoltage2
!          4) setup_ts_kpoint_grid

  
      implicit none

      public :: transiesta, TSVHFix, TSVoltage2, setup_ts_kpoint_grid

      private

      CONTAINS


C ##################################################################
C ##                                                              ##       
C ##                       "TRANSIESTA"                           ##
C ##                                                              ##       
C ##          Non-equilibrium Density Matrix Subroutine           ##
C ##                   to be called from SIESTA                   ##
C ##                                                              ## 
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##               Kurt Stokbro, ks@mic.dtu.dk                    ## 
C ##               Mikroelektronik Centret (MIC)                  ##
C ##           Technical University of Denmark (DTU)              ##
C ##                                                              ##
C ##################################################################
C
C
C Tight-binding density-matrix/transport program for the SIESTA
C package.
C Copyright by Mads Brandbyge, 1999, 2000, 2001, 2002.
C The use of this program is allowed for not-for-profit research only.
C Copy or disemination of all or part of this package is not
C permitted without prior and explicit authorization by the authors.
C
C ================================================================
C Modified by F. D. Novaes to allow calculations with k-point 
C samplings other than gamma point.
C ================================================================
C
C
C Assumes listh and listd are the same
C
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C UNITS:
C     Ry, and Bohr is used throughout the program unless explicitly 
C     stated.
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine transiesta(gamma,lasto,
     .  nuo, no, nspin, nua, maxnh, numh, listhptr,
     .  listh,  xij, indxuo, Hs, Ss, Dnew, Enew, ef, nuotot,
     .  TSiscf,qtot)

#ifdef MPI
      use mpi_siesta, only: MPI_Bcast, MPI_Comm_Rank, DAT_dcomplex,
     . MPI_Sum, DAT_double, MPI_Comm_World, MPI_Integer, MPI_logical
C FDN
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb, 
     .                         LocalToGlobalOrb
C FDN
#endif
      use m_ts_options, only : voltfdf, USEBULK, TriDiag, updatedmcr, 
     .                         FixQ, NBUFATL, NBUFATR, nvolt, GFFileL, 
     .                         GFFileR
      use m_ts_contour, only : setupcontour, distrcontour
      use m_ts_io, only : read_green
      use m_ts_electrode, only : green
#ifdef MPI
      use m_ts_scattering, only : consHSp, getsfe, getGF_Tri, 
     .                            GFGammaGF, weightDM, weightDMC
#else
      use m_ts_scattering, only : consHS, getsfe, getGF_Tri,
     .                            GFGammaGF, weightDM, weightDMC
#endif
      use sys, only : die
      use precision, only : dp

c======================================================================
      implicit none
c======================================================================

      integer NGL2,NA1L,NA2L,nqL,nuaL
      integer NGR2,NA1R,NA2R,nqR,nuaR
      integer NCONTOUR0,NCONTOUR
c======================================================================
      real(dp) eV
      parameter ( eV = 1.d0 / 13.60580d0 )
      real(dp) Pi
      parameter(Pi=3.141592653589793238462643383279503d0)
      real(dp) toPi
      parameter(toPi=2d0/Pi)
c======================================================================


      character*33 gffileLeft,gffileRight

C From fdf read voltfdf, Volt=voltfdf/eV
      real(dp) Volt               ! the bias voltage

C      integer NBUFATL,NBUFATR   ! No. buffer atoms, L/R
      integer NBUFL,NBUFR       ! No. buffer states, L/R



c =====================================================================
c   FROM SIESTA:
c =====================================================================
      integer nua               ! No. atoms in unit cell
      logical gamma             !True if gamma point calculation
      
     

      
      

      integer lasto(0:nua)      ! last-orbital index on atoms
      integer nspin             ! Number of spin components (1 or 2)
c                               ! Can only be 1 for the time being!

      integer no                ! no. orbitals in supercell
      integer nuo               ! no. orbs. in unit cell
      integer maxnh             ! Maximum number of nonzero elements of 
c                               ! each row of hamiltonian matrix
      integer numh(nuo)         ! Number of nonzero elements of each row
      integer listh(maxnh)      ! Nonzero hamiltonian-matrix element
c                               ! column indexes for each matrix row
c                               ! of hamiltonian matrix
      integer listhptr(nuo)     ! Pointer to each row (-1) of the
c                               ! hamiltonian matrix
c                               ! column indexes for each matrix row
      real*8  Hs(maxnh,nspin)   !Hamiltonian in sparse form
      real*8  Ss(maxnh)         ! Overlap in sparse form
      real*8  xij(3,maxnh)      ! Vectors between orbital centers (sparse)
      integer indxuo(no)        ! Index of equivalent orbital in unit cell
c                               !  Unit cell orbitals must be the first in
c                               ! orbital lists, i.e. indxuo.le.nuo, with
c                               ! nuo the number of orbitals in unit cell  
      

      integer TSiscf

c The fermi-energy of the CONTACT (without electrodes) is initially determined
c by SIESTA. This is taken as the starting energy shift of the contact region:
c The electrodes in equilibrium has per definition Fermi energy at ZERO.
      real(dp) ef                 ! Fermi energy (used in CONTACT) 
      real(dp) qtot
c =====================================================================
c   BACK TO SIESTA:
c =====================================================================

      real(dp) Dnew(maxnh,nspin) ! Output Density Matrix
C FDN
      
      complex(dp) ccphase
C FDN
      real(dp) Enew(maxnh,nspin) ! Output Energy-Density Matrix

c =====================================================================
c     Greens functions:
c =====================================================================
c     Energy Contour:
      complex(dp), dimension (:), pointer:: contour,wgf

      integer NPARACONTOUR   ! No. points on parallel contour
c     Energy points and their weight distributed on nodes
      complex(dp), dimension (:,:), allocatable ::
     .     paracontour,paraWGF
c     (point,node) --> part of contour (L,R,V,N) Left,Right,Voltage, No voltage

      character, dimension(:,:), allocatable ::  contourpart

c     No. actual energy points to read-in for given parallel energy point
      integer, dimension(:), allocatable ::  points2read

c     Read-in bulk H, S, for a given q-points
      complex(dp), dimension(:,:,:), allocatable :: 
     .     HAAL,SAAL,HAAR,SAAR 

c     q-points and their weights:
      real(dp), dimension (:,:), pointer:: qL,qR
      real(dp), dimension (:), pointer:: wqL,wqR

c ==================================================================
c     Left/Right inverse Greens functions:
      integer NGL,NGR           !no. orbs. in L/R electrode parts

c     Electrode GF's
c     Gamma=i*0.5*(Sigma - Sigma^dagger) 
      complex(dp), dimension (:,:), allocatable:: 
     .     SFEL,SFER,GAMMAL,GAMMAR

c ==================================================================
c Fermi levels
      real(dp) EFermiL,EFermiR,EFermi0
c===================================================================

c Internal representation of H,S,DM and other quantities:    

      integer nou               !no. of orbs. from SIESTA minus buffer-orbitals
c                               !i.e. no. used orbitals

c Full Green's function, H and S of contact:
      complex(dp), dimension(:), allocatable :: H,S
      complex(dp), dimension(:), pointer :: GF

c     To save space we will save the temporary DML,EDML in DM,EDM. 
c     contour part of DM integration 
c     real axis part of DM integration
c     L,R estimates of Edm

      real(dp), dimension(:), allocatable :: 
     .     DM,DMR,DMneqL,DMneqR,EDM,EDMR

C FDN
      complex(dp), dimension(:), allocatable ::
     .  DMCplx,DMRCplx,DMneqLCplx,DMneqRCplx,EDMCplx,EDMRCplx,
     .  work1Cplx

      logical :: kpargamma
      real(dp) :: kpartol
C FDN


      real(dp) QC0                ! Init(SIESTA) no. electrons in update region
      real(dp) QC                 ! no. electrons in update region
      real(dp) QCn0

c No. contour points (NCONTOUR if VOLT=0 else NCONTOUR + NVOLT)
      integer NEn

c     No. orbitals on atoms
      integer, allocatable :: norbslist(:) 
      integer maxbaso           !Max. no. of basis orbitals on atoms

c ==================================================================
c Helpers, workspace, tempos, auxilaries etc...
      
      real(dp) k(3),kxij,qcloc
      complex(dp) cphase
      
      integer ind
      integer jgfL,jgfR         !File units for Left,Right GF

c      L uc. atoms (from read-in superc)
c      R uc. atoms (from read-in superc)
c      Index of equiv. orb. in uc 

      integer, dimension(:), allocatable :: indxuo1,lastoL,lastoR
      
      complex(dp), dimension(:), allocatable :: work2
      real(dp), dimension(:), allocatable :: work1

      complex(dp), dimension(:), allocatable :: dosdum

            

      complex(dp) Z,W
      integer iPE,ispin

      integer i,j,ic,jc
      integer ia,iau,ia2
      integer io,jo,iuo,juo,iio,is
      integer ierr,itmp         !matrix inversion failure
      integer Nodes
      integer nuotot

#ifdef MPI
      integer  BNode,  maxnhg
      integer, dimension(:), allocatable,save ::
     .  numhg, listhptrg, listhg
 
      real(dp), dimension(:,:), allocatable, save ::
     .  xijg, Dnewg, Enewg
      logical, dimension (:), allocatable ::  listudg
#endif

      logical existgf,ioinC,joinC,updatedm,isvolt,errorgf
      logical errorgs
      logical, dimension (:), allocatable ::  listud

      real(dp) renorm

      logical firsttime
      data firsttime /.true./

      character  Part

c Files
      external io_assign, io_close
      integer jout
      
c Save 
      save firsttime 
      save gffileLeft,gffileRight
      save listud 

#ifdef MPI
      save listudg
#endif

      save nou,contour,WGF,QC0,QCN0,norbslist,maxbaso

      save NPARACONTOUR,paracontour,paraWGF
      save contourpart,points2read

      save EFermi0,EFermiL,EFermiR,NEn
      save Volt,isvolt
      save NBUFL,NBUFR
      save lastoL,lastoR

      save NGL2,NA1L,NA2L,nqL,nuaL
      save NGR2,NA1R,NA2R,nqR,nuaR

C FDN
      save NGL,NGR,kpargamma
C FDN

      save jout
      save NCONTOUR,NCONTOUR0
      save indxuo1
      
      complex(dp) zdotc
      external zdotc ! BLAS function

      character paste*33
      external paste
      
      logical IONode
      integer Node

#ifdef MPI
      integer MPIerror
#endif

C FDN
c ==================================================================
c           K-points variables
c ==================================================================
      integer, save ::  nkpar
      integer :: ikpar
      real(dp), dimension(:,:), allocatable, save :: kpar
      real(dp), dimension(:), allocatable, save :: wkpar
      logical, save :: diffkpar=.false.

C  Preciso colocar um caso geral de eletrodos com potencialmente
C  um numero diferente de pontos k, como por exemplo um eletrodo sendo
C  um bulk e o outro um nanotubo

C FDN

c=======================================================================
c BEGIN:
c=======================================================================


C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror) 
#else
      Node = 0
      Nodes = 1 
#endif
      IONode=.false.
      if(Node.eq.0) IONode=.true.
      
      nullify(GF)

c ============================================================
c       ----------- BEGIN: first time only ------------
      if(firsttime) then
c ============================================================
         
c     Open output file):
c
         if(Node.eq.0) then
c     call TSprversion
            jout = 6
         end if !node 0

         Volt=voltfdf/eV 
 
         isvolt=.false.
         if(dabs(Volt).GT.0.001) isvolt=.true.
                  

#ifdef MPI
c note: efis not global variables!!
      call MPI_Bcast(ef,1,DAT_double,0, MPI_Comm_World,MPIerror)
#endif

c
c     We take the common Fermi energy of the electrodes to be zero!
c
         EFermiL =  Volt*eV/2.0d0 !Ry
         EFermiR = -Volt*eV/2.0d0 !Ry
         EFermi0 =       0.0d0

c
c     Setup contour
c
         call setupcontour(jout,NEn,EFermi0,EFermiL,EFermiR,
     .                  WGF,contour,NCONTOUR,NCONTOUR0)
      

c     Distribute contour points on nodes
c
         allocate( paracontour(1:NCONTOUR,0:Nodes-1)  )
         allocate( paraWGF(1:NCONTOUR,0:Nodes-1)      )
         allocate( contourpart(1:NCONTOUR,0:Nodes-1)  )
         allocate( points2read(1:NCONTOUR)            )


c (note we're allocating more than we need here if there're more than 1 node)

         call distrcontour(jout,
     .        WGF,contour,NCONTOUR,NCONTOUR0,NVOLT,isvolt,Nodes,
     .        NPARACONTOUR,paracontour,paraWGF,
     .        contourpart,points2read)


c  the first NBUFATL atoms will be removed
         NBUFL=0
         do ia=1,NBUFATL
            NBUFL=NBUFL+(lasto(ia)-lasto(ia-1))
         end do                 !ia
c     the last NBUFATR atoms will be removed
         NBUFR=0
         do ia=nua - NBUFATR+1,nua
            NBUFR=NBUFR+(lasto(ia)-lasto(ia-1))
         end do                 !ia

c
c No. states minus buffers
c     

#ifdef MPI
         nou = nuotot - (NBUFL+NBUFR)
#else
         nou = nuo - (NBUFL+NBUFR)
#endif

      endif ! firsttime



      if(firsttime) then
c -------------------------------------------------------
c
c     GF generation:
c

         if(Node.eq.0) then
            allocate(dosdum(NEn))
C FDN nspin added as variable
            CALL green(jout,
     .           NEn,contour,WGF,EFermiL,dosdum,.true.,nspin )
            CALL green(jout,
     .           NEn,contour,WGF,EFermiR,dosdum,.false.,nspin )
C FDN
            
            gffileLeft=GFFileL
            gffileRight=GFFileR

            inquire(file=gffileLeft,exist=existgf)
            if(.not.existgf) 
     .           write(6,*) 'ERROR: Cannot open ',gffileLeft
            inquire(file=gffileRight,exist=existgf)
            if(.not.existgf) 
     .           write(6,*) 'ERROR: Cannot open ',gffileRight
         endif                  !node=0


#ifdef MPI
      call MPI_Barrier(MPI_Comm_World,MPIerror)
      call MPI_Bcast(existgf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
      if(.not.existgf) call die('ERROR: Cannot open GF file')


c -------------------------------------------------------
c     Read-in header of Green's functions
c
       nullify(wqR); nullify(wqL); nullify(qR); nullify(qL)

c Left


       if(Node.eq.0) then
         call io_assign(jgfL)
         OPEN(FILE=gffileLeft,UNIT=jgfL,FORM='UNFORMATTED')
         call read_green(jgfL,NEn,contour,WGF,
     .        nuaL,NA1L,NA2L,nqL,NGL2,wqL,qL,errorgf)
         call io_close(jgfL)
        endif     

#ifdef MPI 
       call MPI_Bcast(errorgf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

       if(errorgf) call die('ERROR in reading Left GF file') 

#ifdef MPI
      call MPI_Bcast(nuaL,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NA1L,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NA2L,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nqL,1,MPI_integer,0,MPI_Comm_World,MPIerror)


      if(Node.ne.0) then
C FDN in general q// may have z component
       allocate(qL(3,nqL))
C FDN
       allocate(wqL(nqL))
      endif
C FDN 2*nqL -> 3*nqL
      call MPI_Bcast(qL(1,1),3*nqL,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(wqL,nqL,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NGL2,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

c --------------------------
c Right

      if(Node.eq.0) then

c     Open and read header of GF-file

       call io_assign(jgfR)
       OPEN(FILE=gffileRight,UNIT=jgfR,FORM='UNFORMATTED')
       call read_green(jgfR,NEn,contour,WGF,
     .      nuaR,NA1R,NA2R,nqR,NGR2,wqR,qR,errorgf)
       call io_close(jgfR)

C FDN Compare the kpoints of Left and Right
       if (nqR == nqL ) then
        kpartol=0.0001
        do I = 1,nqL
          do J = 1,3
            if ( abs(qL(J,I)-qR(J,I)) > kpartol ) then
              diffkpar=.true.
            end if
          enddo
        enddo 
        if ( diffkpar ) then
         write(*,*) 'Left/Right k-par must be equal !!'
         write(*,*) 'qL'
         do I = 1,nqL
           write(*,*) qL(1:3,I)
         enddo
         write(*,*) 'qR'
         do I = 1,nqR
           write(*,*) qR(1:3,I)
         enddo
         call die('Stopping Code!')
        end if 

        nkpar=nqL
        if(allocated(kpar)) deallocate(kpar)
        allocate(kpar(3,nkpar))
        kpar=qL
        if(allocated(wkpar)) deallocate(wkpar)
        allocate(wkpar(nkpar))
        wkpar=wqL

       else
        write(*,*) 'Left/Right k-par must be equal !!'
        write(*,*) 'nqR=',nqR,'nqL=',nqL
        call die('Stopping Code!')

       endif
       if ( nkpar > 1 ) then
        kpargamma=.false.
       else
        if ( kpar(1,1) /= 0d0 .or. 
     .       kpar(2,1) /= 0d0 .or.
     .       kpar(3,1) /= 0d0 ) then
             
         kpargamma=.false.     
        
        else
       
         kpargamma=.true.         
 
        end if 
          
       end if
C FDN

      endif

C FDN
#ifdef MPI
      call MPI_Bcast(nkpar,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      if (Node /= 0 ) then
        if(allocated(kpar)) deallocate(kpar)
        allocate(kpar(3,nkpar))
        if(allocated(wkpar)) deallocate(wkpar)
        allocate(wkpar(nkpar))
      end if
      call MPI_Bcast(kpar(1,1),3*nkpar,DAT_double,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(wkpar,nkpar,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(kpargamma,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
C FDN

#ifdef MPI
       call MPI_Bcast(errorgf,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

       if(errorgf) call die('ERROR in reading Right GF file') 

#ifdef MPI
      call MPI_Bcast(nuaR,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NA1R,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NA2R,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nqR,1,MPI_integer,0,MPI_Comm_World,MPIerror)

      if(Node.ne.0) then
C FDN
       allocate(qR(3,nqR))
       allocate(wqR(nqR))
      endif

C FDN
      call MPI_Bcast(qR(1,1),3*nqR,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(wqR,nqR,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(NGR2,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

c -------------------------------------------------------    
         NGL=NGL2*NA1L*NA2L
         NGR=NGR2*NA1R*NA2R
c -------------------------------------------------------    
         
c
c expected no. states on Electrode atoms: 
c
         allocate(lastoL(0:nuaL)) 
         lastoL(0)=0
         ia2=0

         do ia=NBUFATL+1,NBUFATL + nuaL*NA1L*NA2L, NA1L*NA2L
            ia2=ia2+1
            lastoL(ia2)=lastoL(ia2-1) + (lasto(ia) - lasto(ia-1))
         end do                 !ia


         if(lastoL(nuaL) .ne. NGL2) then
            if(Node.eq.0)
     .      write(6,*) 'ERROR: lastoL,NGL2',lastoL,NGL2
            call die('ERROR: Unexpected no. orbs. in L elec.')
         end if

         allocate(lastoR(0:nuaR))    
         lastoR(0)=0
         ia2=0

         do ia=nua-(nuaR*NA1R*NA2R + NBUFATR)+1,nua-NBUFATR, NA1R*NA2R  
            ia2=ia2+1
            lastoR(ia2)=lastoR(ia2-1) + (lasto(ia) - lasto(ia-1))
         end do                 !ia



         if(lastoR(nuaR) .ne. NGR2) then
            if(Node.eq.0)
     .      write(6,*) 'ERROR: lastoR,NGR2',lastoR,NGR2

            call die('ERROR: Unexpected no. orbs. in R elec.')
         end if

         if(Node.eq.0) then
         write(jout,*) 'L-ELEC: lastoL: ',(lastoL(ia),ia=0,nuaL)
         write(jout,*) 'R-ELEC: lastoR: ',(lastoR(ia),ia=0,nuaR)
         endif


#ifdef MPI
         if(.not. gamma) then
         allocate(indxuo1(no))          
            indxuo1 = indxuo    !just a copy of indxuo
         else                   !gamma:
         allocate(indxuo1(nuotot))          
            do io=1,nuotot
               indxuo1(io)=io   !just orb. number
            end do
         end if
#else
         if(.not. gamma) then
         allocate(indxuo1(no))
            indxuo1 = indxuo    !just a copy of indxuo
         else                   !gamma:
         allocate(indxuo1(nuo))
            do io=1,nuo
               indxuo1(io)=io   !just orb. number
            end do
         end if
#endif

c
c No. orbitals on each atom:
c
         allocate(norbslist(nua))
         maxbaso=0
         do iau=1,nua
            norbslist(iau)=lasto(iau)-lasto(iau-1)
            if(norbslist(iau).gt. maxbaso) maxbaso = norbslist(iau)
         end do                 !iau


      end if                    !firsttime
c -------------------------------------------------------


#ifdef MPI

C Allocate local memory for global list arrays

      allocate(numhg(nuotot))
      allocate(listhptrg(nuotot))

C Globalise numh
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numhg(io) = numh(iio)
        endif

        call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
      enddo

C Build global listhptr
      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo

C Globalize listh
      maxnhg = listhptrg(nuotot)+numhg(nuotot)
      allocate(listhg(maxnhg))

      do io = 1,nuotot
      
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
           call GlobalToLocalOrb(io,Node,Nodes,iio)
           listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) =
     .          listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
        endif
        
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     .    BNode,MPI_Comm_World,MPIerror)
      enddo

      allocate(Dnewg(maxnhg,nspin))
      allocate(Enewg(maxnhg,nspin))
      if(.not.gamma) allocate(xijg(3,maxnhg))


      do io = 1,nuotot
         call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(iio)
              Dnewg(listhptrg(io)+jo,is) = Dnew(listhptr(iio)+jo,is)
              Enewg(listhptrg(io)+jo,is) = Enew(listhptr(iio)+jo,is)
            enddo
          enddo
          if(.not.gamma) then
           do jo = 1,numh(iio)
             xijg(1:3,listhptrg(io)+jo) = xij(1:3,listhptr(iio)+jo)
           enddo
          endif
        endif

        do is = 1,nspin
          call MPI_Bcast(Enewg(listhptrg(io)+1,is),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
          call MPI_Bcast(Dnewg(listhptrg(io)+1,is),numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)
        enddo

        if(.not.gamma) 
     .  call MPI_Bcast(xijg(1,listhptrg(io)+1),3*numhg(io),DAT_double,
     .      BNode,MPI_Comm_World,MPIerror)

      enddo

#endif

c FIRSTTIME

      if(TSiscf.eq.1) then
c     if(firsttime) then
         if(Node.eq.0) write(6,*)"TRANSIESTA: Initializing lists"

c -------------------------------------------------------
c
c     IMPORTANT: Here we set the region which we do not change in DM
c                and initialize.
c

         if(allocated(listud)) deallocate(listud)
         allocate(listud(maxnh))
         listud = .false.

#ifdef MPI
         if(allocated(listudg)) deallocate(listudg)
         allocate(listudg(maxnhg))
         listudg = .false.
#endif

      do io = 1,nuo
#ifdef MPI
            call  LocalToGlobalOrb(io, Node, Nodes, iio)
            iuo = indxuo1(iio)
#else
            iuo = indxuo1(io)
#endif
            ic = iuo - NBUFL
            do j = 1,numh(io)
               ind = listhptr(io) + j
               jo = listh(ind)
               juo = indxuo1(jo)

               jc = juo - NBUFL

               ioinC=.false.
               joinC=.false.

               if(UseBulk) then !update only where Ham is updated:
                  if(ic.gt.NGL .and. ic.lt.nou-NGR+1) ioinC=.true.
                  if(jc.gt.NGL .and. jc.lt.nou-NGR+1) joinC=.true.

C FDN From fdf Now
                if( .not. updatedmcr ) then                
c
c              update the cross DM terms between R,L and  C
c
                  updatedm = ioinC .or. joinC
                else
c
c              update DM only in C
c
                 updatedm = ioinC .and. joinC
C FDN
                end if
               else
                  if(ic.ge.1 .and. ic.le.nou) ioinC=.true.
                  if(jc.ge.1 .and. jc.le.nou) joinC=.true.
                  updatedm = ioinC .and. joinC
               end if           !UseBulk

c     Always fix bonds between Left/Right electrodes:
               if(jc.le.NGL .and. ic.ge.(nou-NGR+1))
     .              updatedm=.false.
               if(ic.le.NGL .and. jc.ge.(nou-NGR+1))
     .              updatedm=.false.
               listud(ind)=updatedm
            end do              !j
         end do                 !io

#ifdef MPI
         do io = 1,nuotot
            iuo = indxuo1(io)
            ic = iuo - NBUFL

            do j = 1,numhg(io)
               ind = listhptrg(io) + j
               jo = listhg(ind)
               juo = indxuo1(jo)
               jc = juo - NBUFL

               ioinC=.false.
               joinC=.false.

               if(UseBulk) then !update only where Ham is updated:
                  if(ic.gt.NGL .and. ic.lt.nou-NGR+1) ioinC=.true.
                  if(jc.gt.NGL .and. jc.lt.nou-NGR+1) joinC=.true.
C FDN From fdf Now
                  if( .not. updatedmcr ) then                
c
c              update the cross DM terms between R,L and  C
c
                     updatedm = ioinC .or. joinC
                  else
c
c              update DM only in C
c
                     updatedm = ioinC .and. joinC
C FDN
                  end if
               else
                  if(ic.ge.1 .and. ic.le.nou) ioinC=.true.
                  if(jc.ge.1 .and. jc.le.nou) joinC=.true.
                  updatedm = ioinC .and. joinC
               end if           !UseBulk

               if(jc.le.NGL .and. ic.ge.(nou-NGR+1))
     .              updatedm=.false.
               if(ic.le.NGL .and. jc.ge.(nou-NGR+1))
     .              updatedm=.false.
               listudg(ind)=updatedm
            end do              !j
         end do                 !io
#endif


c -------------------------------------------------------
c
c     Calculate initial L-C-R region charge
c

         QCn0 = 0d0
         do ispin = 1,nspin
            do ind = 1,maxnh
                 if(.not.listud(ind)) then
                      QCn0 = QCn0 + Dnew(ind,ispin)*Ss(ind)
                 endif
            end do              !io
         end do                 !ispin

#ifdef MPI
      call MPI_AllReduce(qcn0,qcloc,1,
     .  DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
       qcn0 = qcloc

#endif

        if(Node.eq.0) then 
           write(jout,'(a,f8.3)') 'From SIESTA: Efermi=',ef/eV
           write(jout,'(a,f8.3)')'qc0 :',Qc0
           write(jout,'(a,f8.3)')'qcn0 :',Qcn0
           write(jout,*) 'nou,NGL,NGR:',nou,NGL,NGR
        end if    


        qc0 = qtot -qcn0

c      end if                    !firsttime
      end if                    !TSiscf ==1

   

c----------------------------------------------
c Allocation

        allocate(H(nou*nou) )
        allocate(S(nou*nou))
C FDN Complex for k point sampling
        if ( kpargamma ) then
         allocate(DM(nou*nou))
         allocate(EDM(nou*nou))
        else
         allocate(DMCplx(nou*nou))
         allocate(EDMCplx(nou*nou))
        end if
C FDN

      
        allocate(SFEL(NGL,NGL))
        allocate(SFER(NGR,NGR))
        allocate(GAMMAL(NGL,NGL))
        allocate(GAMMAR(NGR,NGR))
     
C FDN New dimensions
C        allocate(HAAL(NGL2,NGL2,nqL))
C        allocate(SAAL(NGL2,NGL2,nqL))
C        allocate(HAAR(NGR2,NGR2,nqR))
C        allocate(SAAR(NGR2,NGR2,nqR))
        allocate(HAAL(NGL2,NGL2,1))
        allocate(SAAL(NGL2,NGL2,1))
        allocate(HAAR(NGR2,NGR2,1))
        allocate(SAAR(NGR2,NGR2,1))
C FDN


        if(isvolt)then
C FDN Complex for k point sampling
         allocate(work2(nou*nou))
         if ( kpargamma ) then      
          allocate(DMR(nou*nou))
          allocate(DMneqL(nou*nou))
          allocate(DMneqR(nou*nou))
          allocate(EDMR(nou*nou))
         else
          allocate(DMRCplx(nou*nou))
          allocate(DMneqLCplx(nou*nou))
          allocate(DMneqRCplx(nou*nou))
          allocate(EDMRCplx(nou*nou))
         end if
        endif

c----------------------------------------------


C FDN  Put outside kpar loop

c -------------------------------------------------------
c     Read-in header of Green's functions
c
      nullify(wqR); nullify(wqL); nullify(qR); nullify(qL)



         if(Node.eq.0) then
c Left

         call io_assign(jgfL)
          open(file=gffileLeft,unit=jgfL,form='unformatted')

          call read_green(jgfL,NEn,contour,WGF,
     .              nuaL,NA1L,NA2L,nqL,NGL2,wqL,qL,errorgf)

          endif
#ifdef MPI

      call MPI_Bcast(nqL,1,MPI_integer,0,MPI_Comm_World,MPIerror)

      if(Node.ne.0) then
C FDN
         allocate(qL(3,nqL))
         allocate(wqL(nqL))
      endif

C FDN
      call MPI_Bcast(qL(1,1),3*nqL,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(wqL,nqL,DAT_double,0,MPI_Comm_World,MPIerror)
#endif

c --------------------------
c Right
         if(Node.eq.0) then
          call io_assign(jgfR)
          open(file=gffileRight,unit=jgfR,form='unformatted')
          call read_green(jgfR,NEn,contour,WGF,
     .              nuaR,NA1R,NA2R,nqR,NGR2,wqR,qR,errorgf)
         endif

#ifdef MPI
      call MPI_Bcast(nqR,1,MPI_integer,0,MPI_Comm_World,MPIerror)

      if(Node.ne.0) then
C FDN
        allocate(qR(3,nqR))
        allocate(wqR(nqR))
      endif

C FDN
      call MPI_Bcast(qR(1,1),3*nqR,DAT_double,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(wqR,nqR,DAT_double,0,MPI_Comm_World,MPIerror)
#endif

C FDN Spin loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ispin = 1,nspin !not really implemented yet!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C FDN kpar loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do ikpar=1,nkpar           ! kpar loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         k(:)=kpar(:,ikpar)


c     
c     Init DM, EDM
c
C FDN 
        if (kpargamma) then    
         DM    =  0d0
         EDM   =  0d0
        else
C FDN Complex
         DMCplx=dcmplx(0d0,0d0)
         EDMCplx=dcmplx(0d0,0d0)
        end if
C FDN

         if(isvolt) then
C FDN
          if (kpargamma) then            
            DMR    = 0d0
            DMneqL = 0d0
            DMneqR = 0d0
            EDMR   = 0d0
          else
C FDN Complex
            DMRCplx    = dcmplx(0d0,0d0)
            DMneqLCplx = dcmplx(0d0,0d0)
            DMneqRCplx = dcmplx(0d0,0d0)
            EDMRCplx   = dcmplx(0d0,0d0)
          end if
C FDN
 
         end if                

c setup hamiltonian

#ifdef MPI 

         call  consHSp(nuo, no,   maxnh,  maxnhg,
     .                 gamma, numh, listhptr,
     .                 numhg, listhptrg,listhg,
     .                 Hs(1,ispin), Ss, xijg, indxuo1, k, nuotot,
     .                 H,S, Node, Nodes,ef,nou,NGL,NGR,NBUFL)

#else

         call consHS( nuo, no,  maxnh,gamma, numh, listhptr,
     .                listh,Hs(1,ispin), Ss, xij, indxuo1, 
     .                k,H,S,ef,nou, NGL,NGR,NBUFL)
#endif

      
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
      do 7913 iPE = 1, NPARACONTOUR ! loop over parallel contour points
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=

C         allocate(GF(nou*nou))
         

c       if(Node.eq.0.and.iPE.eq.1) write(6,*)iPE

         Z = paracontour(iPE,Node)    ! ZEnergy point for this node
         Part = contourpart(iPE,Node) ! and place on the contour

c     Weight of energy point:
c     This is where we include the factor-of-two for spin and
c     and the (1/Pi) from DM = Im[G]/Pi
         if(nspin.eq.1) then
            W = (2d0/Pi)*paraWGF(iPE,Node)
         else
            W = (1d0/Pi)*paraWGF(iPE,Node)
         end if

c----------------------------------------------------------------------
c At this point we need the transformation from
c G(q) in (1x1) to G(0) in (NA1xNA2)
c
c SFEL/SFER  <------ GAAL/GAAR

c         if(Node.eq.0) write(6,*) "hola getsfe L", ipe
c LEFT:



         call getsfe(USEBULK,jgfL,HAAL,SAAL,Z,
     &        NGL2,NGL,SFEL,GAMMAL,
     &        Nodes,points2read(iPE),errorgs,ikpar)

C FDN

#ifdef MPI
       call MPI_Bcast(errorgs,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
       if(errorgs) call die('ERROR in getsfe') 

c        if(Node.eq.0) write(6,*) "hola getsfe R", ipe
c RIGHT:

         call getsfe(USEBULK,jgfR,HAAR,SAAR,Z,
     &        NGR2,NGR,SFER,GAMMAR,
     &        Nodes,points2read(iPE),errorgs,ikpar)

C FDN
          
#ifdef MPI
       call MPI_Bcast(errorgs,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
       if(errorgs) call die('ERROR in getsfe') 

       
c Calculate Full GF at energy point
       ierr=0

       call  getGF_Tri(USEBULK,nou,NGL,NGR,H,S,SFEL,SFER,Z,GF,ierr,
     . TriDiag)
      
#ifdef MPI
C FDN
      if (Nodes.gt.1)  then
        call MPI_AllReduce(ierr,itmp,1,MPI_integer,MPI_sum,
     .       MPI_Comm_World,MPIerror)      
        ierr=itmp 
      end if
C FDN     
#endif
      IF(IERR.NE.0) THEN
         if(Node.eq.0) then
            write(jout,*) 'ERROR: getGF MATRIX INVERSION FAILED'
         endif 
      END IF


c -------------------------------------------------------------
c     No voltage/Left Equilibrium Contour:
         if(Part .eq. 'N' .or. Part .eq. 'L') then
            do i=1,nou*nou
C FDN
              if (kpargamma) then
               DM(i)  =  DM(i) - DIMAG(GF(i)*W)
               EDM(i) = EDM(i) - DIMAG(Z*GF(i)*W)
              else
C FDN Complex for k sampling
               DMCplx(i) = DMCplx(i) - GF(i)*W
               EDMCplx(i) = EDMCplx(i) - Z*GF(i)*W
              end if
C FDN
            end do              !i
    

c -------------------------------------------------------------
c     Right Equilibrium Contour:
         else if(Part .eq. 'R') then
            do i=1,nou*nou
C FDN
              if (kpargamma) then
               DMR(i)  =  DMR(i) - DIMAG(GF(i)*W)
               EDMR(i) = EDMR(i) - DIMAG(Z*GF(i)*W)
              else 
C FDN Complxe for k sampling
               DMRCplx(i)  =  DMRCplx(i) - GF(i)*W
               EDMRCplx(i) = EDMRCplx(i) - Z*GF(i)*W
              end if
C FDN
            enddo

        
c -------------------------------------------------------------
c     Non-equilibrium part:
         else if(Part .eq. 'V') then


c   Calculate rho^R(e) -> work2 ! NOTE R!!
        call GFGammaGF(.false.,nou,NGR,GF,GAMMAR,work2)
        
c note DM = DML+DMneqR or DM=DMR+DMneqL
c note we use "+" here!

            do i=1,nou*nou
C FDN
              if (kpargamma) then
               DMneqR(i) = DMneqR(i) +DIMAG(work2(i)*W)
               EDM(i)    = EDM(i) + DIMAG(Z*work2(i)*W)
              else
C FDN Complexe for k sampling
               DMneqRCplx(i) = DMneqRCplx(i) +work2(i)*W
               EDMCplx(i)    = EDMCplx(i) + Z*work2(i)*W
              end if
C FDN
            end do
c  Calculate rho^L(e) -> work2 ! NOTE L!!

        call GFGammaGF(.true.,nou,NGL,GF,GAMMAL,work2)
        
c     note we use "-" here!
            do i=1,nou*nou
C FDN
              if (kpargamma) then
               DMneqL(i) = DMneqL(i) - DIMAG(work2(i)*W)
               EDMR(i)   = EDMR(i) + DIMAG(Z*work2(i)*W)
              else
C FDN Complex for k sampling
               DMneqLCplx(i) = DMneqLCplx(i) - work2(i)*W
               EDMRCplx(i)   = EDMRCplx(i) + Z*work2(i)*W
              end if
C FDN
            end do
c -------------------------------------------------------------
         else
            if(Node.eq.0) 
     .      write(6,*) 'TRANSIESTA: ERROR in contour setup!'
           call die ('TRANSIESTA: ERROR in contour setup!')
         end if                 !which part of the contour ??
! As a pointer it MUST be deallocated
         deallocate(GF)
         

C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
 7913    continue               ! the loop over ALL contour points
c                               ! on each processor
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=

c     on return the result is in DM and EDM

#ifdef MPI

C Global reduction of density matrices

C FDN kpargamma
      if (kpargamma) then
       allocate(work1(nou*nou))

       work1 = 0.0d0
       call MPI_AllReduce(DM,work1,nou*nou,
     .     DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
       DM = work1

       work1 = 0.0d0
       call MPI_AllReduce(EDM,work1,nou*nou,
     .     DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
       EDM = work1
      
       if (isvolt) then

          work1 = 0.0d0
          call MPI_AllReduce(DMR,work1,nou*nou,
     .        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
          DMR = work1
          work1 = 0.0d0

          call MPI_AllReduce(DMneqL,work1,nou*nou,
     .        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
          DMneqL = work1
          work1 = 0.0d0
 
          call MPI_AllReduce(DMneqR,work1,nou*nou,
     .        DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         DMneqR = work1
         work1 = 0.0d0
         call MPI_AllReduce(EDMR,work1,nou*nou,
     .       DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
         EDMR = work1
        
       endif

       deallocate(work1)
      else

       allocate(work1Cplx(nou*nou))

       work1Cplx = 0.0d0
       call MPI_AllReduce(DMCplx,work1Cplx,nou*nou,
     .     DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
       DMCplx = work1Cplx

       work1Cplx = 0.0d0

       call MPI_AllReduce(EDMCplx,work1Cplx,nou*nou,
     .     DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
       EDMCplx = work1Cplx

       if (isvolt) then

          work1Cplx = 0.0d0
          call MPI_AllReduce(DMRCplx,work1Cplx,nou*nou,
     .        DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
          DMRCplx = work1Cplx
          work1Cplx = 0.0d0

          call MPI_AllReduce(DMneqLCplx,work1Cplx,nou*nou,
     .        DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
          DMneqLCplx = work1Cplx
          work1Cplx = 0.0d0
 
          call MPI_AllReduce(DMneqRCplx,work1Cplx,nou*nou,
     .        DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
         DMneqRCplx = work1Cplx
         work1Cplx = 0.0d0
         call MPI_AllReduce(EDMRCplx,work1Cplx,nou*nou,
     .       DAT_dcomplex,MPI_sum,MPI_Comm_World,MPIerror)
         EDMRCplx = work1Cplx
        
       endif

       deallocate(work1Cplx)


      end if
C FDN kpragamma
#endif

      if (isvolt  ) then 
C FDN
        if ( ikpar == nkpar .and. ispin == nspin ) then
         deallocate(work2)
        end if
       if (kpargamma) then
         call weightDM(nou,NGL,NGR,DM,DMR,DMneqL,
     .        DMneqR,EDM,EDMR)
       else
         call weightDMC(nou,NGL,NGR,DMCplx,DMRCplx,
     . DMneqLCplx,DMneqRCplx,EDMCplx,EDMRCplx)
       end if
C FDN
      endif
      
c--------------------------------------------------------------------
c
c     The original Hamiltonian from SIESTA was shifted ef: 
c  -- thus we have a shift in EDM:

C FDN
      if (kpargamma) then 
       call daxpy(nou*nou,ef,DM,1,EDM,1)
      else 
C FDN
       call zaxpy(nou*nou,dcmplx(ef,0.d0),DMCplx,1,EDMCplx,1)
      end if
C FDN

c
c     Update Dnew, Enew written in sparse SIESTA-format:
c     

c first set Dnew,Enew = 0 for elements we want to update


#ifdef MPI
C FDN Initialized only once
        if ( ikpar == 1 ) then
        do ind = 1,maxnhg
          if(listudg(ind)) then
            Dnewg(ind,ispin) = 0d0
            Enewg(ind,ispin) = 0d0
          end if
        end do
        endif
C FDN


        do io = 1,nuotot
           iuo = indxuo1(io)
           ic = iuo - NBUFL
           do j = 1,numhg(io)
             ind = listhptrg(io) + j
             jo = listhg(ind)
             juo = indxuo1(jo)
             jc = juo - NBUFL
             if(listudg(ind)) then ! update
C FDN
               if(kpargamma) then
                Dnewg(ind,ispin) = Dnewg(ind,ispin) +
     .               DM(ic+nou*(jc-1))
                Enewg(ind,ispin) = Enewg(ind,ispin) +
     .               EDM(ic+nou*(jc-1))
               else
                kxij = (k(1)*xijg(1,ind) +
     .                 k(2)*xijg(2,ind) +
     .                 k(3)*xijg(3,ind) )
                cphase =  exp( dcmplx(0d0,-1.0d0)*kxij )
                ccphase = exp( dcmplx(0d0,+1.0d0)*kxij )
      
        Dnewg(ind,ispin) = Dnewg(ind,ispin) + Dimag(
     . dcmplx(0.5d0*wkpar(ikpar),0d0)*(
     . cphase*DMCplx(ic+nou*(jc-1))+ccphase*DMCplx(jc+nou*(ic-1))
     .                                )          )

        Enewg(ind,ispin) = Enewg(ind,ispin) + Dimag(
     . dcmplx(0.5d0*wkpar(ikpar),0d0)*(
     . cphase*EDMCplx(ic+nou*(jc-1))+ccphase*EDMCplx(jc+nou*(ic-1))
     .                                )          )

               end if
C FDN
             end if             !update DM,EDM
       
          end do                !j
       end do                   !io
! FDN Teste
!      stop
#else
C FDN Initialized only once
        if ( ikpar == 1 ) then
        do ind = 1,maxnh
          if(listud(ind)) then
            Dnew(ind,ispin) = 0d0
            Enew(ind,ispin) = 0d0
          end if
        end do
        endif
C FDN

        do io = 1,nuo
           iuo = indxuo1(io)
           ic = iuo - NBUFL 
           do j = 1,numh(io)
              ind = listhptr(io) + j
              jo = listh(ind)
              juo = indxuo1(jo)
              jc = juo - NBUFL 
              if(listud(ind)) then ! update

C FDN
                if ( kpargamma ) then

          Dnew(ind,ispin) = Dnew(ind,ispin) +
     .               DM(ic+nou*(jc-1))
          Enew(ind,ispin) = Enew(ind,ispin) +
     .               EDM(ic+nou*(jc-1))
            

                else

                 kxij = (k(1)*xij(1,ind) +
     .                 k(2)*xij(2,ind) +
     .                 k(3)*xij(3,ind) )
                 cphase =  exp( dcmplx(0d0,-1.0d0)*kxij )
                 ccphase = exp( dcmplx(0d0,+1.0d0)*kxij )

        Dnew(ind,ispin) = Dnew(ind,ispin) + Dimag(
     . dcmplx(0.5d0*wkpar(ikpar),0d0)*(
     . cphase*DMCplx(ic+nou*(jc-1))+ccphase*DMCplx(jc+nou*(ic-1))
     .                                )          )

        Enew(ind,ispin) = Enew(ind,ispin) + Dimag(
     . dcmplx(0.5d0*wkpar(ikpar),0d0)*(
     . cphase*EDMCplx(ic+nou*(jc-1))+ccphase*EDMCplx(jc+nou*(ic-1))
     .                                )          )

                end if 

              end if             !update DM,EDM


C FDN

               
          end do                !j
       end do                   !io
#endif


C FDN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        end do               ! ikpar
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C FDN



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do                    !ispin
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


#ifdef MPI

C Distribute Dnew/Enew 

      do io = 1,nuo
        call LocalToGlobalOrb(io,Node,Nodes,iio)
          do is = 1,nspin
            do jo = 1,numh(io)
              Dnew(listhptr(io)+jo,is)= Dnewg(listhptrg(iio)+jo,is)
              Enew(listhptr(io)+jo,is)= Enewg(listhptrg(iio)+jo,is)
            enddo
          enddo
      enddo

       deallocate(Dnewg)
       deallocate(Enewg)

#endif


c 
c     Charge Renormalization -- Fix charge in update region
c
      QC = 0d0
      do ispin = 1,nspin
         do ind = 1,maxnh
            if(listud(ind))
     .           QC = QC + Dnew(ind,ispin)*Ss(ind)
         end do                 !io
      end do                    !ispin

#ifdef MPI
      call MPI_AllReduce(qc,qcloc,1,
     .DAT_double,MPI_sum,MPI_Comm_World,MPIerror)
      qc = qcloc
#endif


      renorm = 1.0d0
      
c      if(FIXQ) then
c         renorm = QC0/QC
c         do ispin = 1,nspin
c            do ind = 1,maxnh
c               if(listud(ind)) then
c                  Dnew(ind,ispin)=Dnew(ind,ispin)*renorm
c                  Enew(ind,ispin)=Enew(ind,ispin)*renorm
c               end if 
c            end do                 
c         end do                 !ispin
c      end if                    !FIXQ

      if(Node.eq.0) then
c     write(jout,'(a25,1X,3(f8.3,1X))')
c     .        'TranSiesta: QC,QC0,renorm:',QC,QC0,renorm
         write(jout,'(a22,1X,3(f8.3,1X))')
     .        'TranSiesta: Qsol,Qtot:',QC+qcn0,qtot,renorm
      endif 
c--------------------------------------------------------------------

c     Close Files

      if(Node.eq.0) then
       call io_close(jgfL)
       call io_close(jgfR)
      endif
      
c--------------------------------------------------------------------

c     Deallocate:

      deallocate(H)
      deallocate(S)
C FDN
      if ( kpargamma ) then
       deallocate(DM)
       deallocate(EDM)
      else
       deallocate(DMCplx)
       deallocate(EDMCplx)
      end if
C FDN
c      deallocate(GF)

      deallocate(SFEL)
      deallocate(SFER)
      deallocate(GAMMAL)
      deallocate(GAMMAR)

      deallocate(HAAL)
      deallocate(SAAL)
      deallocate(HAAR)
      deallocate(SAAR)

      
      if(isvolt) then
C FDN
      if ( kpargamma ) then
         deallocate(DMR)
         deallocate(DMneqL)
         deallocate(DMneqR)
         deallocate(EDMR)
       else
         deallocate(DMRCplx)
         deallocate(DMneqLCplx)
         deallocate(DMneqRCplx)
         deallocate(EDMRCplx)
       end if
C FDN  
      endif


#ifdef MPI
      deallocate(numhg)
      deallocate(listhptrg)
      deallocate(listhg)
      if(.not.gamma) deallocate(xijg)
#endif

c ===================================================================
      firsttime = .false.       !first time -- not any more!
      return                    ! succesfull return
      end subroutine transiesta


!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      subroutine TSVHFix(  mesh, nsm, v )
C
C  Modules
C
      use precision, only : dp, grid_p
C FDN put the only attribute
      use parallel, only : ProcessorY
C FDN
#ifdef MPI
      use mpi_siesta, only : MPI_AllReduce, MPI_Comm_Rank, MPI_Sum,
     . DAT_double, MPI_Integer, MPI_Comm_World
C FDN parallelsubs
      use parallelsubs, only : HowManyMeshPerNode
C FDN
#endif

      implicit          none
      real(grid_p)      v(*)
      integer           mesh(3), nsm


c Internal variables
      integer           i1, i2, i3, imesh, ntemp,
     .                  nlp, Node, Nodes, meshl(3),
     .                  ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset,
     .                  Zoffset, Py, Pz, i30, i20, meshnsm(3),
     .                  NRemY, NRemZ
#ifdef MPI
      integer           MPIerror, npl
#endif
      real(dp)          vav, vtot, temp


C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif


C Find local number of mesh points
        meshnsm(1) = mesh(1)/nsm
        meshnsm(2) = mesh(2)/nsm
        meshnsm(3) = mesh(3)/nsm
#ifdef MPI
        !! print *, "N:",Node, "TSVHFix-- mesh:", mesh
        call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
        meshl(1) = meshl(1)*nsm
        meshl(2) = meshl(2)*nsm
        meshl(3) = meshl(3)*nsm
#else
        meshl(1) = mesh(1)
        meshl(2) = mesh(2)
        meshl(3) = mesh(3)
#endif

C Check that ProcessorY is a factor of the number of processors
        if (mod(Nodes,ProcessorY).gt.0) then
          write(6,'(''ERROR: ProcessorY must be a factor of the'',
     .      '' number of processors!'')')
          stop
        endif
        ProcessorZ = Nodes/ProcessorY

C Calculate blocking sizes
        BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
        NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
        BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
        NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

C Calculate coordinates of current node in processor grid
        Py = (Node/ProcessorZ)+1
        Pz = Node - (Py - 1)*ProcessorZ + 1

C Calculate starting point for grid
        Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
        Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

         vtot =0.d0
         nlp =0

        imesh = 0
        i30 = Zoffset - 1
        do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            do i1 = 0,meshl(1)-1
              imesh = imesh + 1
              if (i30.eq.0) then
                 nlp = nlp +1 
                 vtot = vtot + v(imesh)
              endif
            enddo
          enddo
        enddo

#ifdef MPI
      temp=0.d0
        call MPI_AllReduce(vtot,temp,1,DAT_double,MPI_Sum,
     .  MPI_Comm_World,MPIerror)
       vtot = temp
       ntemp=0
        call MPI_AllReduce(nlp,ntemp,1,MPI_integer,MPI_Sum,
     .    MPI_Comm_World,MPIerror)
       nlp = ntemp
 
#endif
        vav=vtot/float(nlp)

        imesh = 0
        
        do i3 = 0,meshl(3)-1
           do i2 = 0,meshl(2)-1
              do i1 = 0,meshl(1)-1
                 imesh = imesh + 1
                 v(imesh) = v(imesh) - vav
              enddo
           enddo
        enddo


c        if(Node.eq.0) write(6,*)'Av Pot',Node,vav
      end subroutine TSVHFix

!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      subroutine TSvoltage2(cell, mesh, nsm, v, isvolt )


C
C  Modules
C
      use precision, only : dp, grid_p
      use fdf, only : fdf_physical
C FDN put the only attribute
      use parallel, only : ProcessorY
C FDN
#ifdef MPI
C FDN parallelsubs
      use parallelsubs, only : HowManyMeshPerNode
C FDN
      use mpi_siesta, only : MPI_Bcast, MPI_Comm_Rank, DAT_double,
     . MPI_Comm_World
#endif

      implicit          none

      integer           mesh(3), nsm
      real(grid_p)      v(*)
      real(dp)          cell(3,3), dot, field(3)
      logical           isvolt 
      external          dot


c Internal variables
      real(dp) eV
      parameter ( eV = 1.d0 / 13.60580d0 )
 
      logical           frstme
      integer           i1, i2, i3, imesh,  ix, ivc,
     .                  Node, Nodes, meshl(3),
     .                  ProcessorZ, BlockSizeY, BlockSizeZ, Yoffset,
     .                  Zoffset, Py, Pz, i30, i20, meshnsm(3),
     .                  NRemY, NRemZ

#ifdef MPI
      integer           MPIerror, npl
#endif

      real(dp)  f(3),v0,Volt,Lvc,dLvc,vcdir(3),VoltL,VoltR

      parameter(ivc=3)          ! The voltage-direction

      save              f,frstme,Volt,VoltL,VoltR,v0,Lvc


      data frstme /.true./

C Get Node number
#ifdef MPI
      call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
      call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

c Find and store the electric field only the first time
      if (frstme) then
        frstme = .false.
        isvolt = .false.
        do ix = 1,3
          field(ix) = 0.d0
        enddo

 
c     Read the voltage block from the fdf input file
       if (Node.eq.0)Volt = fdf_physical('TS.Voltage',0.0d0,'Ry')

#ifdef MPI
      call MPI_Bcast(Volt,1,DAT_double,0,MPI_Comm_World, MPIerror)
#endif


         if(dabs(Volt).GT.0.001d0*eV) then 
            isvolt=.true.
         else
            isvolt=.false.
            return              !not found and field set to zero in beginning
         endif                  !found
c

         Lvc = sqrt(dot(cell(1,ivc),cell(1,ivc),3))

         do ix=1,3
            vcdir(ix) = cell(ix,ivc)/Lvc
            field(ix) = -vcdir(ix)*Volt/Lvc
         end do                 !i


         if(Node.eq.0) then        
          write(6,*)
          write(6,'(a,f6.3,a)')'TSvoltage: Bias ', Volt/eV,'V'
          write(6,'(a,3(f6.3,a))')
     .         'TSvoltage: In unit cell direction = {',
     .         vcdir(1),',',vcdir(2),',',vcdir(3),'}'
         endif 

      end if                    !frstme

c --------------------- end: first time only -----------------------




c     Find quantities in mesh coordinates

         dLvc = Lvc/max( mesh(ivc), 1 ) !

         VoltL = 0.5d0*Volt     ! Left voltage
         VoltR = -0.5d0*Volt    ! Right voltage

c     field in [0;Lvc]: v = e*x = f*index
         do ix = 1,3
            f(ix) = 0d0         !init
         enddo
         f(ivc) = -Volt*dLvc/Lvc
         v0 = VoltL

      if (isvolt) then

C Find local number of mesh points
        meshnsm(1) = mesh(1)/nsm
        meshnsm(2) = mesh(2)/nsm
        meshnsm(3) = mesh(3)/nsm
#ifdef MPI
        !! print *, "N:",Node, "TSVoltage2-- mesh:", mesh
        call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
        meshl(1) = meshl(1)*nsm
        meshl(2) = meshl(2)*nsm
        meshl(3) = meshl(3)*nsm
#else
        meshl(1) = mesh(1)
        meshl(2) = mesh(2)
        meshl(3) = mesh(3)
#endif

C Check that ProcessorY is a factor of the number of processors
        if (mod(Nodes,ProcessorY).gt.0) then
          write(6,'(''ERROR: ProcessorY must be a factor of the'',
     .      '' number of processors!'')')
          stop
        endif
        ProcessorZ = Nodes/ProcessorY

C Calculate blocking sizes
        BlockSizeY = (meshnsm(2)/ProcessorY)*nsm
        NRemY = (mesh(2) - ProcessorY*BlockSizeY)/nsm
        BlockSizeZ = (meshnsm(3)/ProcessorZ)*nsm
        NRemZ = (mesh(3) - ProcessorZ*BlockSizeZ)/nsm

C Calculate coordinates of current node in processor grid
        Py = (Node/ProcessorZ)+1
        Pz = Node - (Py - 1)*ProcessorZ + 1

C Calculate starting point for grid
        Yoffset = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY)
        Zoffset = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ)

c Add the electric field potential to the input potential
        imesh = 0
        i30 = Zoffset - 1
        do i3 = 0,meshl(3)-1
          i30 = i30 + 1
          i20 = Yoffset - 1
          do i2 = 0,meshl(2)-1
            i20 = i20 + 1
            do i1 = 0,meshl(1)-1
              imesh = imesh + 1
              v(imesh) = v(imesh) + v0 + f(1)*i1 + f(2)*i20 + f(3)*i30
            enddo
          enddo
        enddo

      endif

      end subroutine TSvoltage2



!-------------------------------------------------------------------------
!*************************************************************************
!-------------------------------------------------------------------------


      subroutine setup_ts_kpoint_grid( ucell )

! SIESTA Modules
      USE fdf, only       : fdf_defined
      USE m_find_kgrid, only : find_kgrid
      USE parallel, only  : Node
      USE precision, only : dp       
#ifdef MPI
      USE mpi_siesta, only : MPI_Bcast, MPI_logical, MPI_Comm_World
#endif

! TS Modules
      USE m_ts_kpoints, only : ts_kpoint, ts_nkpnt, ts_eff_kgrid_cutoff,
     . ts_firm_displ, ts_kdispl, ts_gamma_scf, ts_kscell, ts_kweight, 
     . ts_scf_kgrid_first_time, ts_spiral, ts_user_requested_mp,ts_maxk,
     . setup_ts_scf_kscell, setup_ts_scf_kscell, ts_write_k_points

! Local Variables
      implicit none
      real(dp) :: ucell(3,3)

#ifdef MPI
      integer :: MPIerror
#endif

      if (ts_scf_kgrid_first_time) then
       nullify(ts_kweight,ts_kpoint)
       if (Node.eq.0) then
        ts_spiral = fdf_defined('SpinSpiral')
       endif
#ifdef MPI
       call MPI_Bcast(ts_spiral,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

       call setup_ts_scf_kscell(ucell, ts_firm_displ)

       ts_scf_kgrid_first_time = .false.

      else
       if ( ts_user_requested_mp    ) then
          ! no need to set up the kscell again
       else
          ! This was wrong in the old code
          call setup_ts_scf_kscell(ucell, ts_firm_displ)
       endif
      endif

      call find_kgrid(ucell,ts_kscell,ts_kdispl,ts_firm_displ,     
     .             (.not. ts_spiral),                    
     .             ts_nkpnt,ts_kpoint,ts_kweight, ts_eff_kgrid_cutoff)

      ts_maxk = ts_nkpnt
      ts_gamma_scf =  (ts_nkpnt == 1 .and.  
     .         dot_product(ts_kpoint(:,1),ts_kpoint(:,1)) < 1.0e-20_dp)

      if (Node .eq. 0) call ts_write_k_points()

      end subroutine setup_ts_kpoint_grid

      END MODULE m_ts_in_siesta
