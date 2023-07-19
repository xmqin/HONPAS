MODULE m_ts_io
!
! Routines that are used for Input and Output of files 
!
!=============================================================================
! CONTAINS:
!          1) read_green
!          2) ts_iohs


  implicit none

  public :: read_green, ts_iohs

  private

  CONTAINS


! ##################################################################
! ##            read-in header of Greens function file            ##
! ##                            By                                ##
! ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
! ##                                                              ##
! ## Changed to F90 by Jose-Luis Mozos, jlm@icmab.es              ##
! ##################################################################

      subroutine read_green(jgfu,NEn,contour,wgf, &
                              nua,NA1,NA2,nq,ng,wq,q,errorgf)


      use precision, only: dp

      implicit none

      real(dp) EPS
      parameter(EPS=1d-7)

      logical PRINTALOT
!      parameter(PRINTALOT=.FALSE.)
      parameter(PRINTALOT=.TRUE.)

! INPUT
      integer jgfu               !unit of gf-file


! these are the values expected:
! they will be compared to the read-in and ERROR messages will appear
      integer NEn,ng
      complex(dp) contour(NEn),wgf(NEn)
      integer nua,NA1,NA2       !no. atoms in uc and no. repetitions in A1,A2
      integer nq                !no. q-points

      logical errorgf


! READ-IN values
      real(dp) efermii            !The required Fermi energy
      integer NEni
      complex(dp), dimension(:), allocatable :: contouri,wgfi

      integer nqi

!     q-point=k_|| point and weigth

      real(dp), dimension (:,:), pointer:: q
      real(dp), dimension (:), pointer:: wq

      character*70 gftitle

! Helpers..

      complex(dp) ctmp
      integer iEn


!=======================================================================
! BEGIN:
!=======================================================================

       errorgf = .false.       


      read(jgfu) gftitle
      read(jgfu) EFermii,NEni

      read(jgfu) nua,NA1,NA2,nqi

!      write(6,*) 'read GF: ',gftitle

      if(NEni .ne. NEn)  then
         write(6,*) 'read GF: ERROR: NEn=',NEni,' expected:', NEn
            errorgf = .true. 
            return
      end if


      allocate(contouri(nen))
      allocate(wgfi(nen))

      nq = nqi

! FDN in general q// may have z component
!     dimension changed from 2 to 3
      allocate(q(3,nqi))
! FDN
      allocate(wq(nqi))

      read(jgfu) contouri,wgfi,q,wq
      read(jgfu) ng


!     check on contours.
      do iEn=1,NEn

         ctmp=contouri(iEn)-contour(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(*,*) ' Warning: contours differ by >', EPS
         end if
         if(cdabs(ctmp).GT.10d0*EPS) then 
            write(*,*) &
                 ' ERROR: contours differ by >', 10.d0*EPS
            errorgf = .true. 
            return
         end if
      end do

      do iEn=1,NEn
         ctmp=wgfi(iEn)-wgf(iEn)
         if(cdabs(ctmp).GT.EPS) then 
            write(*,*)  &
                ' ERROR: contour weights differ by >',EPS
            errorgf = .true. 
            return
         end if
      end do

      deallocate(contouri)
      deallocate(wgfi)
      




      return
      end subroutine read_green





!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------


! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine ts_iohs(task, gamma, onlyS, no_u, no_s, Enspin,   &
      		 indxuo, maxnh, numh, listhptr, listh, H, S, qtot, &
		 temp, xij, fnlength, fname, na_u, lasto, isa, ef, &
		 ucell, ts_kscell_file, ts_kdispl_file,        	   &
		 ts_gamma_scf_file, xa, istep, ia1)

! *********************************************************************
! Saves the hamiltonian and overlap matrices, and other data required
! to obtain the bands and density of states
! Writen by J.Soler July 1997.
! Note because of the new more compact method of storing H and S
! this routine is NOT backwards compatible
! Modified by M.Paulsson 2009 to:
! 1: To include information of which FC step for phonon calculations
! 2: To only save the overlap matrix if onlyS flag is set
!    (Used for e-ph coupling calculations)
! 3: File format changed to unify Copenhagen/Barcelona Transiesta vers.
! 4: Smaller files by writing arrays directly instead of element wise
! *************************** INPUT **********************************
! character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
! logical       gamma         : Is only gamma point used?
! logical       onlyS         : Should only overlap matrix be saved?
! ******************** INPUT or OUTPUT (depending on task) ***********
! integer no_u                : Number of basis orbitals per unit cell
! integer no_s                : Number of basis orbitals per supercell
! integer Enspin              : Spin polarization (1 or 2)
! integer indxuo(no_s)        : Index of orbitals in supercell
! integer maxnh               : First dimension of listh, H, S and
!                               second of xij
! integer numh(nuo)           : Number of nonzero elements of each row
!                               of hamiltonian matrix
! integer listhptr(nuo)       : Pointer to the start of each row (-1)
!                               of hamiltonian matrix
! integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
!                               indexes for each matrix row
! real*8  H(maxnh,Enspin)     : Hamiltonian in sparse form
! real*8  S(maxnh)            : Overlap in sparse form
! real*8  qtot                : Total number of electrons
! real*8  temp                : Electronic temperature for Fermi smearing
! real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
!                               (not read/written if only gamma point)
! TSS Begin
! ********************* ADDED ARGUMENTS FOR TRANSIESTA ****************
! integer fnlength            : file name length
! character(fnlength)  fname  : file nema for input or output
! integer na_u                : Number of atoms per unit cell
! integer istep, ia1          : Force constant step and atom number
! TSS End
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! *********************************************************************

!
!  Modules
!
      use precision,    only : dp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb, &
                               GlobalToLocalOrb, GetNodeOrbs
      
      use files,        only : slabel, label_length
      use sys,          only : die
      use m_ts_kpoints, only: ts_gamma_scf, ts_kscell, ts_kdispl
! TSS End
#ifdef MPI
      use mpi_siesta
#endif
! TSS Begin
! SIESTA Modules
      use m_spin, only: nspin       
! TSS End

      implicit          none

      character(8) :: task
      character(len=label_length+3) :: paste
      logical           gamma, onlyS, onlySfile
      integer           maxnh, no_u, no_s, Enspin
! TSS
!      integer           indxuo(no_s), listh(maxnh), numh(*), listhptr(*)
      integer, dimension(:), pointer :: indxuo, listh, numh, listhptr 
! TSS
!      real(dp)          H(maxnh,Enspin), S(maxnh), &
!                        qtot, temp, xij(3,maxnh)
      real(dp), dimension(:,:), pointer :: H, xij
      real(dp), dimension(:), pointer :: S
      real(dp) :: qtot, temp 
      integer :: istep, ia1
      external          io_assign, io_close, paste
! TSS Begin
! Added arguments for TRANSIESTA
      integer :: fnlength, na_u
      character(fnlength) :: fname
      integer, dimension(:), pointer :: lasto, isa
      integer, dimension(3,3) :: ts_kscell_file
      real(dp) :: ef
      real(dp), dimension(3,3) :: ucell
      real(dp), dimension(:,:), pointer :: xa
      real(dp), dimension(3) :: ts_kdispl_file 
      logical :: ts_gamma_scf_file
! TSS End

! Internal variables and arrays
      integer    is, iu, k, ns
      integer    ih,hl, maxnhtot, maxhg
      integer, dimension(:), allocatable :: numhg
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
      real(dp), dimension(:,:), allocatable :: buffer2
#endif
      logical    found, gammaonfile
! TSS Begin
! Aux Variables
      integer :: i, j, i2, l
! TSS End



! #################### READ ######################################
! Choose between read or write
      if (task.eq.'read' .or. task.eq.'READ') then

        ! Note that only the master node is suppossed to
        ! call this routine in "read" mode...
    
! Check if input file exists
          inquire( file=fname, exist=found)

        if (found) then

! Open file
            call io_assign( iu )
            open( iu, file=fname, form='unformatted', status='old' )      

! Read dimensions
            read(iu) na_u, no_u, no_s, Enspin, maxnh

! Check dimensions

            if (Enspin  .ne. nspin) then
               if (node == 0) then
                 write(6,"(a,i10,a,i10)") "nspin mismatch: In "   &
                      // trim(fname) // " : ",                    &
                      Enspin, ". Current run: ", nspin
               endif
               call die('ts_iohs: TSHS file contains a wrong nspin')
            endif
           
! TSS Begin
! Allocate arrays that are going to be read now
          nullify(xa,isa)
          allocate(xa(3,na_u))
          allocate(isa(na_u)) 

! Read Geometry information
           read(iu) xa
           read(iu) isa   
       	   read(iu) ucell  

! Read k-point sampling information
           read(iu) gammaonfile   
	   read(iu) onlySfile
           read(iu) ts_gamma_scf_file       
           read(iu) ts_kscell_file
           read(iu) ts_kdispl_file  
	   read(iu) istep, ia1

! Check whether file is compatible from a gamma point of view
           if (onlySfile) then
             write(*,*) 'TSHS file does not contain Hamiltonian'
             call die('ts_iohs: onlyS flag, no H in file')
           endif
           if (.not.gamma.and.gammaonfile) then
             write(*,*) 'System Gamma:',gamma
             write(*,*) 'Electrode Gamma:',gammaonfile
             call die('ts_iohs: Non-TSgamma information not present')
           endif
           if (.not.ts_gamma_scf .and. ts_gamma_scf_file) then
             write(*,*) 'TS_Electrode Gamma:',ts_gamma_scf_file
             write(*,*) 'TS_System Gamma:',ts_gamma_scf
             call die('ts_iohs: Incompatible k samplings!') 
           end if
           do i = 1,3
             if ( ts_kdispl_file(i) /= ts_kdispl(i) ) then
               write(*,'(a,2F8.4)') 'TS_Electrode k displacements:',(ts_kdispl_file(j),j=1,2)
               write(*,'(a,2F8.4)') 'TS_System k displacements:', (ts_kdispl(j),j=1,2)
               call die('ts_iohs: Incompatibles k displacements') 
             end if
           enddo
           do i = 1,3
             do j = 1,3
               if ( ts_kscell(i,j)  /= ts_kscell_file(i,j) ) then
                 do i2=1,2
                   write(*,*) 'TS_Electrode k points:', (ts_kscell_file(i2,l),l=1,3)
                 enddo 
                 do i2=1,2
                   write(*,*) 'TS_System k points:', (ts_kscell(i2,l),l=1,3)
                 enddo
                 call die('ts_iohs: Incompatible k-point samplings')
               end if
             enddo
           enddo


! Read sparse listings
          nullify(lasto)
          allocate(lasto(0:na_u))
            read(iu) lasto


          if (.not.gamma) then
! Allocate arrays that are going to be read now
            nullify(indxuo)
            allocate(indxuo(1:no_s))
              read(iu) (indxuo(ih),ih=1,no_s)

          endif

! Allocate local array for global numh
          nullify(numh)
          allocate(numh(no_u))
          call memory('A','I',no_u,'iohs')

! Read numh and send to appropriate Node
	    read(iu) numh(1:no_u)



! Read Electronic Structure Information
            read(iu) qtot,temp
            read(iu) ef


! Create listhptr
          nullify(listhptr)
          allocate(listhptr(no_u))
          listhptr(1) = 0
! TSS nuo->no_u
! TSS hl->hi
          do ih = 2,no_u
            listhptr(ih) = listhptr(ih-1) + numh(ih-1)
          enddo

! Read listh
! Allocate lish
          nullify(listh)
          allocate(listh(maxnh))

           do ih = 1,no_u
             read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih))
           enddo

! Read Overlap matrix
! Allocate S
          nullify(S)
          allocate(S(maxnh))
           do ih = 1,no_u
             read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
           enddo

! Read Hamiltonian
! Allocate H
          nullify(H)
          allocate(H(maxnh,Enspin))
           do is = 1,Enspin
             do ih = 1,no_u
               read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
             enddo
           enddo

          if (.not.gamma) then
! Read interorbital vectors for K point phasing
! Allocate xij
            nullify(xij)
            allocate(xij(3,maxnh))
             do ih = 1,no_u
               read(iu) (xij(k,listhptr(ih)+1:listhptr(ih)+numh(ih)),k=1,3)
             enddo

          endif


! Close file
             call io_close( iu )

        else
            write(*,*) 'iohs: ERROR: file not found: ', fname
            call die('iohs: ERROR: file not found')
        endif
! FDN Temp Fim
! #################### WRITE ######################################
      elseif (task.eq.'write' .or. task.eq.'WRITE') then

! Find total numbers over all Nodes
#ifdef MPI
         call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum, &
         MPI_Comm_World,MPIerror)
#else
         maxnhtot = maxnh
#endif

        if (Node.eq.0) then
! Open file
          call io_assign( iu )
          open( iu, file=fname, form='unformatted', status='unknown' )      
! Write Dimensions Information
          write(iu) na_u, no_u, no_s, Enspin, maxnhtot

! Write Geometry information
          write(iu) xa(1:3,1:na_u)
          write(iu) isa(1:na_u)          
          write(iu) ucell

! Write k-point samplung information
          write(iu) gamma
          write(iu) onlyS
          write(iu) ts_gamma_scf_file
          write(iu) ts_kscell_file
          write(iu) ts_kdispl_file
	  write(iu) istep, ia1

! Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

! Write sparse listings
            write(iu) lasto(0:na_u)

          if (.not.gamma) then
            write(iu) (indxuo(ih),ih=1,no_s)
          endif

        endif

! Create globalised numh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            numhg(ih) = numh(hl)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(numh(hl),1,MPI_integer, &
              0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numhg(ih),1,MPI_integer, &
              BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

        if (Node.eq.0) then
! Write numh ( Sparse listings )
          maxhg = 0
          do ih = 1,no_u
            maxhg = max(maxhg,numhg(ih))
          enddo
	  write(iu) numhg(1:no_u)
#ifdef MPI
          allocate(buffer(maxhg))
          call memory('A','D',maxhg,'iohs')
          allocate(ibuffer(maxhg))
          call memory('A','I',maxhg,'iohs')
#endif



! Write Electronic Structure Information
          write(iu) qtot,temp
          write(iu) ef
        endif
! Write listh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
	    write(iu) listh(listhptr(hl)+1:listhptr(hl)+numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numhg(ih),MPI_integer,BNode,1, &
             MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(listh(listhptr(hl)+1),numh(hl),MPI_integer, &
             0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
                write(iu) ibuffer(1:numhg(ih))
            endif
          endif
#endif
        enddo
#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'iohs')
          deallocate(ibuffer)
        endif
#endif

! Write Overlap matrix
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
	    write(iu) S(listhptr(hl)+1:listhptr(hl)+numh(hl))
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
             BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(S(listhptr(hl)+1),numh(hl), &
             MPI_double_precision,0,1,MPI_Comm_World, &
             Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
	      write(iu) buffer(1:numhg(ih))
            endif
          endif
#endif
        enddo

	if (.not. onlyS) then
! Write Hamiltonian	 
        do is=1,Enspin	 
          do ih=1,no_u	 
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu) H(listhptr(hl)+1:listhptr(hl)+numh(hl),is)
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numhg(ih),MPI_double_precision, &
                BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(H(listhptr(hl)+1,is),numh(hl), &
               MPI_double_precision,0,1,MPI_Comm_World, &
               Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
	        write(iu) buffer(1:numhg(ih))
              endif
            endif
#endif
          enddo
        enddo
	endif  ! onlyS

#ifdef MPI
          if (Node .eq. 0) then
! Free buffer array
             call memory('D','D',size(buffer),'iohs')
             deallocate(buffer)
          endif
#endif


        if (.not.gamma) then
#ifdef MPI
! Allocate buffer array
          if (Node .eq. 0) then
             allocate(buffer2(3,maxhg))
             call memory('A','D',3*maxhg,'iohs')
          endif
#endif
          do ih = 1,no_u
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu) (xij(k,listhptr(hl)+1:listhptr(hl)+numh(hl)),k=1,3)
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer2(1,1),3*numhg(ih), &
               MPI_double_precision,BNode,1,MPI_Comm_World, &
               Request,MPIerror) 
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(xij(1,listhptr(hl)+1),3*numh(hl), &
               MPI_double_precision,0,1,MPI_Comm_World, &
               Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                write(iu) (buffer2(k,1:numhg(ih)),k=1,3)
              endif
            endif
#endif
          enddo
#ifdef MPI
          if (Node .eq. 0) then
! Free buffer array
             call memory('D','D',size(buffer2),'iohs')
             deallocate(buffer2)
          endif
#endif
        endif   ! not gamma

        if (Node.eq.0) then
! Deallocate local array for global numh
          call memory('D','I',size(numhg),'iohs')
          deallocate(numhg)
! Close file
          call io_close( iu )
        endif

      endif


      end subroutine ts_iohs


END MODULE m_ts_io
