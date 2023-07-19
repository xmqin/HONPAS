C ##################################################################
C ##                                                              ## 
C ##   Setup complex energy contour by calling contour routines   ##
C ##                                                              ##   
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Changed to F90 by Jose-Luis Mozos, jlm@kanigo.icmab.es       ##
C ##################################################################

      subroutine setupcontour(joutfile,NEn,
     .     EFermi0,EFermiL,EFermiR,
     .     WGF,contour,ncontour,ncontour0,nvolt)


c======================================================================
      implicit none
c======================================================================
c INPUT
      integer joutfile          !output file unit
      real*8 EFermi0            !Eq. Fermi energy
      real*8 EFermiL,EFermiR    !Left/Right Fermi energies (finite voltage)

c     Read from FDF:
      real*8 Emin,Emax,GFeta,kT
      integer nline,npol,ncircle,ncontour,nvolt,ncontour0

c OUTPUT

      complex*16, dimension (:), pointer:: contour,wgf ! Total "contour"

      integer NEn               ! no. contour points to be used 
c                               ! (may differ from NCONTOUR)

c PARAMS
      real*8 eV
      parameter ( eV = 1.d0/13.60580d0 )


C Dimensions of complex contour (defaults)
      integer NLINED,NPOLD,NCIRCLED,NVOLTD
      parameter (NLINED   =       3 )
      parameter (NPOLD    =       5 )
      parameter (NCIRCLED =      12 )
      parameter (NVOLTD =     5 )
c Helpers etc.
      logical isvolt            !True if finite voltage

      logical baddim,tcmplxcon,trealcon,trealgau

      integer method

      real*8 Volt
      real*8 temp_default       !default temperature
      parameter ( temp_default  = 1.900d-3 ) 


c      complex*16 contmp(NCONTOUR),WGFtmp(NCONTOUR) 
      complex*16, dimension (:), allocatable:: contmp,wgftmp 

      integer i,jc,ju
     
      integer Node

      external io_assign,io_close

      character*20 slabel       ! System Label (to name output files)
      character slabel_default*59
      character paste*33,fname*33
C


c=======================================================================
c BEGIN:
c=======================================================================
C Get Node number
      Node = 0


      isvolt=.false.
      Volt=EFermiL-EFermiR
      if(dabs(Volt).GT.0.001d0) isvolt=.true.

      if(Node.eq.0) then
      slabel_default  = 'gengf'
      slabel = fdf_string('SystemLabel',slabel_default)
      endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   now  the contour info from FDF


      if(Node.eq.0)then


        tcmplxcon = fdf_boolean('TS.ComplexContour',.true.)

        emin=fdf_physical('TS.ComplexContour.Emin',-2.d0,'Ry')
        gfeta=fdf_physical('TS.Eta',0.0005d0,'Ry')

        kT = fdf_physical('ElectronicTemperature',temp_default,'Ry') 

        trealgau =  fdf_boolean('TS.GaussContourReal',.false.)
      endif

        trealcon = .not.tcmplxcon
  
c
c contour for DOS:
c      
      if(trealcon)then


         if(Node.eq.0)then
         emin=fdf_physical('TS.RealContour.Emin',1.d0,'Ry')
         emax=fdf_physical('TS.RealContour.Emax',1.d0,'Ry')
         ncontour=fdf_integer('TS.TBT.NPoints',100)
         endif


         allocate(contour(ncontour))
         allocate(wgf(ncontour))
         call mkRealContour(joutfile,Emin,Emax,
     &         GFeta,NCONTOUR,contour,WGF)
         NEn=NCONTOUR
         if(Node.eq.0) then
          write(joutfile,*) NEn,' energy points'
          call io_assign( ju )
          fname = paste(slabel,'.TSCONTOUR')
          open( ju, file=fname, status='unknown' )
          do i=1,NEn
            write(ju,'(4F8.3,i4)') contour(i),WGF(i),i
          end do
          call io_close( ju )
         endif
         RETURN                 !!!!!!!!!!!! get back
      endif

c
c complex contour for integration of DOS
c
      if(tcmplxcon) then
         if(Node.eq.0) then
          npol=fdf_integer('TS.ComplexContour.NPoles',NPOLD)
          ncircle=fdf_integer('TS.ComplexContour.NCircle',NCIRCLED)
          nline=fdf_integer('TS.ComplexContour.NLine',NLINED)
          nvolt=fdf_integer('TS.Bias.NumPts',NVOLTD)
         endif


         if(isvolt) then
            ncontour=2*(npol+nline+ncircle)+nvolt
         else
            ncontour=  (npol+nline+ncircle)
         end if                 !volt
         ncontour0 = npol+nline+ncircle

         allocate(contour(ncontour))
         allocate(wgf(ncontour))
         

c
c Equilibrium Contour
c
         if(.not.isvolt) then
            allocate(contmp(ncontour0))
            allocate(wgftmp(ncontour0))
            
            call mkCplxContour(
     &           joutfile,1,Emin-Volt*eV,EFermi0,
     &           kT,GFeta,NCONTOUR0,contmp,WGFtmp, !Equilibrium part
     &           ncircle,nline,npol)
            jc=0
            do i=1,NCONTOUR0
               jc=jc+1
               contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
               WGF(jc)=-WGFtmp(i)
               
            end do              !i

            deallocate(contmp)
            deallocate(wgftmp)

         else                   !finite voltage
c
c Non-Equilibrium Contour
c 

c Contour sequence:  
c (i)   equilibrium contour corresponding to LEFT EFermiL 
c (ii)  equilibrium contour corresponding to RIGHT EFermiR 
c (iii) real-axis contour between EFermiL and EfermiR

c (i):
            if(Node.eq.0) write(joutfile,*) 'Left equilibrium:'

            allocate(contmp(ncontour0))
            allocate(wgftmp(ncontour0))

            call mkCplxContour(
     &           joutfile,1,Emin + (EFermiL-EFermi0),EFermiL,
     &           kT,GFeta,NCONTOUR0,contmp,WGFtmp,
     &           ncircle,nline,npol)
            jc=0
            do i=1,NCONTOUR0
               jc=jc+1
               contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
               WGF(jc)=-WGFtmp(i)
            end do              !i

c (ii):
            if(Node.eq.0) write(joutfile,*) 'Right equilibrium:'

            call mkCplxContour(joutfile,1,
     &           Emin + (EFermiR-EFermi0),EFermiR,
     &           kT,GFeta,NCONTOUR0,contmp,WGFtmp,
     &           ncircle,nline,npol)
            do i=1,NCONTOUR0
               jc=jc+1
               contour(jc)=contmp(i)
c     
c     Note we put a minus here because the integral we want is the
c     negative of the pole-sum and C+L integral!!
c     
               WGF(jc)=-WGFtmp(i)
            end do              !i

            deallocate(contmp)
            deallocate(wgftmp)
 
            allocate(contmp(nvolt))
            allocate(wgftmp(nvolt))

c (iii):
            if(Node.eq.0) 
     .      write(joutfile,*) 'Real-axis: From EFermiR to EFermiL'

            if (.not. trealgau) then
               method=2
            else
               method=3
            endif

            call mkCplxContour(joutfile,method,EFermiR,EFermiL,
     &           kT,GFeta,NVOLT,contmp,WGFtmp,
     &           0,0,0)
            do i=1,NVOLT
               jc=jc+1
               contour(jc)=contmp(i)
               WGF(jc)=WGFtmp(i)
            end do              !i

            deallocate(contmp)
            deallocate(wgftmp)

         end if                 !finite voltage


         NEn=jc

         if(Node.eq.0) 
     .   write(joutfile,*) NEn,' energy points'


         if(Node.eq.0) then
         fname = paste(slabel,'.CONTOUR')
         call io_assign( ju )
         open( ju, file=fname, status='unknown' )
         do i=1,NEn
            write(ju,'(4F8.3,i4)') contour(i),WGF(i),i 
         end do
         call io_close( ju )
         end if               
         end if                 !tcmplxcon

c ===================================================================
      return
      end
c ===================================================================
