C ##################################################################
C ## Driver subroutine for calculating the (ideal)                ##
C ## Left surface Greens function.                                ##          
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##                                                              ##
C ## Modified by F. D. Novaes, fdnovaes@icmab.es                  ##
C ##################################################################


      subroutine green(joutfile,nq,q,wq,
     &     NEn,contour,wgf,efermi,zbulk,tjob,tbtjob,
     &     ispin_fdf, k_fdf, zenergy_fdf,iEn_fdf)
      

      use greensub

      use m_tbt_gf
      use m_tbt_options, only : sppol, Lhsfile, Rhsfile, Lnucuse, 
     .             Rnucuse



c======================================================================
      implicit none
c======================================================================

C Dimensions for input to green

      logical PRINTALOT
c      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)
      
c=======================================================================
c INPUT:
      integer nq                ! Number of surface q-points

      real*8 q(3,nq)         ! Surface q-point

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
    
c=======================================================================

c   k_|| == q-points in simple (non-repeated) cell:
      real*8, pointer:: q1(:,:),wq1(:) ! k_|| and their weights 

      real*8, dimension(nq) :: wq

      real*8, pointer:: q2(:,:) ! shifted q1 points
      integer na1,na2           ! Replication of unitcell
      integer nq1                ! no. q1-points (<= na1*na2 for gamma) 
      real*8 kpoint(3)          !3D k-point (q,kz)
      data kpoint /0d0, 0d0, 0d0/
c------------------------------------------------------------------     
c Hamiltonian/Overlap for given k

      integer nspin
      integer nua !no atom in uc
c      integer lasto(0:NG1) ! Index of last orbital of each atom
      integer, dimension (:), pointer:: lasto
c (we use that the number of atoms < NG1)


c ==================================================================
c     OUTPUT GAA:            


      character*33 gffile        ! name of gffile
      integer jgfu               ! output gf-file unit
c      complex*16 HAA(NGAA,nq1)  ! Real-space bulk Hamiltonian
c      complex*16 SAA(NGAA,nq1)  ! Real-space bulk Overlap
c      complex*16 GAAq(NGAA,nq1) ! Inverse ideal GF (z*SAA - HAA - SigmaAA)

c      complex*16 GS(ngaa1)
      complex*16, allocatable :: GS(:),HAA(:,:),SAA(:,:), GAAq(:,:)

      integer nucuse             !Truncate GF to *last/first* nucuse atoms
c                                !         for   *left/right*

      integer, dimension (:), allocatable :: lastou
c      integer lastou(0:NG1) ! Index of last orbital of each *used* atom
c (we use that the number of atoms < NG2)
      integer nou               !no. used orbitals     

C      complex*16 zbulkdos(NEn)
      complex*16 zbulk 

c=======================================================================
c     Helpers, workspace, tempos etc...

      character*70 gftitle      !title of gf
      character*33 paramfile    !parameter file with this "tag"

      logical tinit,tlast

c     LEFT/RIGHT sign in GF calc.: Exp(i*LRsign*Kz*Z):       
      integer LRSign

      character*5 tag
      character stag
      character*6  gfjob
      character*33 gffile_default

      logical exist1
      logical tdos,tkham
      
      complex*16, dimension (:), pointer:: H00, S00, H01, s01

c      complex*16 h00(ngaa1),s00(ngaa1)
c      complex*16 h01(ngaa1),s01(ngaa1)
c      complex*16 ab(ngaa1),ba(ngaa1)
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


      real*8 cell(3,3)
      integer kscell(3,3),j
      real*8 kdispl(3),tol
      complex*16 zenergy_fdf
      integer ispin_fdf,iEn_fdf
      real*8 k_fdf(3)


      external io_assign,io_close

c FDF-stuff:
      character*33 paste,header,itemfdf     
      real*8 fdf_convfac
      external paste

c=======================================================================
c BEGIN:
c=======================================================================
       tdos =.true.


       factor =fdf_convfac('Ry','eV')


      if(tjob)then
c!   left
       gfjob='LEFT  '
       LRSign=1
       tag='Left'
       stag='L'
      else
c!   left
       gfjob='RIGHT '
       LRSign=-1
       tag='Right'
       stag='R'
      endif

c output .GF files to this name


      if ( LFrstTime .or. RFrstTime ) then 


      if (LJob) hsfile = Lhsfile
      if (RJob) hsfile = Rhsfile 

c=======================================================================
c
c Read-in bulk parameters and do calculation:
c
c=======================================================================
        nullify(h00)
        nullify(s00)
        nullify(h01)
        nullify(s01)
        nullify(lasto)

        if (LFrstTime .and. tjob ) then
          nullify(Lh00)
          nullify(Ls00)
          nullify(Lh01)
          nullify(Ls01)
          nullify(Llasto) 
        end if        
        if (RFrstTime .and. (.not.tjob) ) then
          nullify(Rh00)
          nullify(Rs00)
          nullify(Rh01)
          nullify(Rs01)
          nullify(Rlasto)
        end if 


c initialize and get nspinin         
        tinit=.true.
        tlast=.false.
        tkham=.false.
        call sethhm2(joutfile,tinit,tkham,tlast,kpoint,ispin,
     &     hsfile, nua,lasto,ng1,nspin,cell,kscell,kdispl,
     &     H00,s00,h01,s01) ! ->

        if (LJob) then
          Lng1=ng1
          nuaL=nua
          ngaa1=ng1*ng1
          allocate(LGS(ngaa1))

          itemfdf = paste('TS.NumUsedAtoms',tag) 
C          nucuse = fdf_integer(itemfdf,nua) !default use all atoms in uc.
          if (Lnucuse == 0) then
             nucuse=nua
          else 
             nucuse=Lnucuse
          end if
c            write(*,*) "Used no. atoms:",nucuse,nua
c     Truncate lasto to used atoms: lastou:

          allocate(lastou(0:ng1))
          lastou=0

          lastou(0)=0
          ia2=0   
          if(LRsign .EQ. 1) then !Left
c     use nucuse *last* atoms            
              do ia=nua-nucuse ,nua
                 lastou(ia2)=lasto(ia)
                 ia2=ia2+1
              end do
           else                !Right
c     use nucuse *first* atoms            
              do ia=0,nucuse
                 lastou(ia2)=lasto(ia)
                 ia2=ia2+1               
              end do
           end if              !L or R            
         
c     No. used orbitals in uc:
           nou=0
           do ia=1,nucuse
              nou=nou + (lastou(ia)-lastou(ia-1))
           end do

           ngL2 = nou
           Lngaa=ngL2*ngL2
           deallocate(lastou)

C           do iEn=1,NEn
C              zbulkdos(iEn)=dcmplx(0d0,0d0)
C           end do
           


           allocate(LHAA(LNGAA,1))  
           allocate(LSAA(LNGAA,1))
           allocate(LGAAq(LNGAA,1))

        end if
        if (RJob) then
          Rng1=ng1
          nuaR=nua

          ngaa1=ng1*ng1
          allocate(RGS(ngaa1))
          itemfdf = paste('TS.NumUsedAtoms',tag) 
C          nucuse = fdf_integer(itemfdf,nua) !default use all atoms in uc.
          if (Rnucuse == 0) then
             nucuse = nua 
          else
             nucuse=Rnucuse
          end if
c            write(*,*) "Used no. atoms:",nucuse,nua
c     Truncate lasto to used atoms: lastou:

          allocate(lastou(0:ng1))
          lastou=0

          lastou(0)=0
          ia2=0
          if(LRsign .EQ. 1) then !Left
c     use nucuse *last* atoms            
             do ia=nua-nucuse ,nua
                lastou(ia2)=lasto(ia)
                ia2=ia2+1
             end do
          else                !Right
c     use nucuse *first* atoms            
             do ia=0,nucuse
                lastou(ia2)=lasto(ia)
                ia2=ia2+1               
             end do
          end if              !L or R            
         
c     No. used orbitals in uc:
          nou=0
          do ia=1,nucuse
             nou=nou + (lastou(ia)-lastou(ia-1))
          end do

          ngR2 = nou
          Rngaa=ngR2*ngR2       
          deallocate(lastou)


          allocate(RHAA(RNGAA,1))  
          allocate(RSAA(RNGAA,1))
          allocate(RGAAq(RNGAA,1)) 

        end if 


        if(sppol.AND.nspin<2) 
     .     write(joutfile,*) "WARNING/ERROR: 
     .     Spin-polarized electrodes was expected"
        if(.not.sppol .AND. nspin>1)  write(joutfile,*) "WARNING/ERROR: 
     .     Spin-polarized electrodes was NOT expected"

      end if ! LFirstTime .or. RFirstTime



!      if (.not. WGFFiles ) then
         iEn=iEn_fdf

         if (LJob .and. (.not. LFrstTime)) then

            if ( iEn == 1 ) then
              Lh00 =0.0
              Ls00 =0.0
              Lh01 =0.0
              Ls01 =0.0

              tinit=.false.
              tkham=.false.
              tlast=.false.
              call sethhm2(joutfile,tinit,tkham,tlast,k_fdf,
     .              ispin_fdf,hsfile,nuaL,Llasto,Lng1,nspin,
     .              cell,kscell,kdispl,LH00,Ls00,Lh01,Ls01) 
            end if ! iEn == 1
            if(iEn.eq.1) then
                  
               i=0
               do l2 = 1,NGL2
                  do l1 = 1,NGL2
                     i=i+1

          LHAA(i,1) = LH00(l1+(LNG1-NGL2)+LNG1*(l2+(LNG1-NGL2)-1))
          LSAA(i,1) = LS00(l1+(LNG1-NGL2)+LNG1*(l2+(LNG1-NGL2)-1))

                  end do  ! l1
               end do     ! l2
              
               NGAA=ngL2*ngL2
               call zaxpy(NGAA,dcmplx(efermi,0.d0),
     .                 LSAA(1,1),1,LHAA(1,1),1)

            end if        ! iEn.eq.1

            zsenergy = zenergy_fdf-efermi
            call calc_green(Lng1,zsenergy,Lh00,Ls00,Lh01,Ls01,
     &              LGS,zdos,joutfile,tjob,tdos)
C A Fazer .... colocar wk ...           
            if(tdos) zbulk =  zdos


            i=0
                  do l2 = 1,NGL2
                     do l1 = 1,NGL2
                        i=i+1

           LGAAq(i,1) = Lgs(l1+(LNG1-NGL2)+ LNG1*(l2+(LNG1-NGL2)-1) )

                     end do     ! l1
                  end do        ! l2

         end if ! LJob


         if (RJob .and. .not.RFrstTime) then

            if ( iEn == 1 ) then
              Rh00 =0.0
              Rs00 =0.0
              Rh01 =0.0
              Rs01 =0.0

              tinit=.false.
              tkham=.false.
              tlast=.false.
              call sethhm2(joutfile,tinit,tkham,tlast,k_fdf,
     .              ispin_fdf,hsfile,nuaR,Rlasto,Rng1,nspin,
     .              cell,kscell,kdispl,RH00,Rs00,Rh01,Rs01) 
            end if ! iEn
 
            if(iEn.eq.1) then
                 
                     i=0
                     do l2 = 1,NGR2
                        do l1 = 1,NGR2
                           i=i+1

               RHAA(i,1) = RH00(l1+RNG1*(l2-1))
               RSAA(i,1) = RS00(l1+RNG1*(l2-1))

                        end do  ! l1
                     end do     ! l2                  
         
c
c     Shift so efermi is energy zero
c

                  NGAA=ngR2*ngR2
                  call zaxpy(NGAA,dcmplx(efermi,0.d0),
     .                 RSAA(1,1),1,RHAA(1,1),1)


             end if ! iEn.eq.1 

             zsenergy = zenergy_fdf-efermi


             call calc_green(Rng1,zsenergy,Rh00,Rs00,Rh01,Rs01,
     &              RGS,zdos,joutfile,tjob,tdos) 
C A Fazer .... wk
            if(tdos) zbulk =   zdos
             i=0
                  do l2 = 1,NGR2
                     do l1 = 1,NGR2
                        i=i+1

                        RGAAq(i,1) = Rgs(l1+RNG1*(l2-1))

                     end do     ! l1
                  end do        ! l2           

         end if !RJob

!      end if ! .not. WGFFiles


      if (LFrstTime .and. tjob ) then
        LFrstTime=.false.     

      end if
      if (RFrstTime .and. (.not.tjob) ) then
        RFrstTime=.false.

      end if 



C ======================================================================
      return
      end
C ======================================================================
