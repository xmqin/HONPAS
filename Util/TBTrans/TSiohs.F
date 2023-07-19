C FDN kscell and kdispl added as dummy
      subroutine TSiohs(
     .task, finname,gamma,nua,nuotot,notot,nspin,
     .                 maxnh,numh, listhptr, listh, H, S, qtot, temp, 
     .                 xij, indxuo, efs, ucell, isa, lasto, xa,
     .                 kscell,kdispl  )
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C Note because of the new more compact method of storing H and S
C this routine is NOT backwards compatible
C *************************** INPUT **********************************
C character*(*) task          : 'read'/'READ' or 'write'/'WRITE'
C logical       gamma         : Is only gamma point used?
C ******************** INPUT or OUTPUT (depending on task) ***********
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
C integer nuotot              : Number of basis orbitals per supercell
C integer nua                 : No. atoms in unitcell
C integer isa(nua)            : Species index of each atom
C integer lasto(0:nua)        : Index of last orbital of each atom
C real*8 xa(3,nua)            : Atomic coordinates (Bohr)
C real*8 efs(nspin)           : Fermi energy for each spin
c real*8  cell(3,3)           : Unit cell vectors in real space cell(ixyz,ivec)
C integer indxuo(notot)       : Index of equivalent orbital in unit cell
C                             Unit cell orbitals must be the first in
C                             orbital lists, i.e. indxuo.le.nuo, with
C                             nuo the number of orbitals in unit cell


C *************************** UNITS ***********************************
C Units should be consistent between task='read' and 'write'
C *********************************************************************

      use parallel, only : IOnode, Node
      use sys, only : die
#ifdef MPI
      use mpi_siesta, only : MPI_Comm_World
      use m_mpi_utils, only : broadcast
#endif MPI

      implicit          none
C      character         task*(*),finname*33
      character task*4,finname*33    
      logical           gamma, ts_gamma_scf_file
      integer           maxnh, nuotot, nspin
c      integer           listh(maxnh),  listhptr(*)
c      double precision  H(maxnh,nspin), S(maxnh),
c     .                  xij(3,maxnh)
      double precision  qtot, temp
      external          io_assign, io_close

      integer, dimension (:), pointer:: listh,listhptr,
     .                        numh,indxuo,lasto,isa
      double precision, dimension (:,:), pointer:: H,xij,xa
      double precision, dimension (:), pointer:: S,efs

c Internal variables and arrays
      integer    im, is, iu, j, ju, k, mnh, ns, i,ic
      integer    ih,hl,nut
      logical    found
cccc
      integer nua                
      real*8  cell(3,3)            
      integer notot

C FDN
      integer kscell(3,3)
      real*8  ucell(3,3)
      real*8  kdispl(3)
      double precision   ef
! Added to TSHS file by MPN
      logical onlyS
      integer istep, ia1
C FDN

#ifdef MPI
      integer MPIerror
#endif MPI

cccc
C May be ....
      inquire(file=finname,opened=found,number=iu)
      if (found) call io_close(iu)

C Choose between read or write
      
      if (task.eq.'read' .or. task.eq.'READ') then

C Check if input file exists
        if (IOnode) inquire( file=finname, exist=found )
#ifdef MPI
        call broadcast(found)
#endif MPI

        if (found) then
C Open file
            if (IOnode) then
               call io_assign( iu )
              open( iu, file=finname, status='old' ,
     .       form = 'unformatted')
            end if


            if (IOnode) read(iu) nua, nuotot, notot, nspin, maxnh
#ifdef MPI
            call broadcast(nua)
            call broadcast(nuotot)
            call broadcast(notot)
            call broadcast(nspin)
            call broadcast(maxnh)
#endif MPI
       


            nullify(isa)
            allocate(xa(3,nua))
            allocate(isa(nua)) 

            if (IOnode) read(iu) xa
            if (IOnode) read(iu) isa  
            if (IOnode) read(iu) ucell
#ifdef MPI
            call broadcast(xa)
            call broadcast(isa)
            call broadcast(ucell)
#endif MPI
            

            if (IOnode) read(iu) gamma
! Variable added by MPN to be used by Phonons codes, not used by tbtrans
            if (IOnode) read(iu) onlyS
            if (IOnode) read(iu) ts_gamma_scf_file
#ifdef MPI
            call broadcast(gamma)
            call broadcast(ts_gamma_scf_file)
#endif MPI
       
            if (IOnode) read(iu) kscell
            if (IOnode) read(iu) kdispl 
#ifdef MPI
            call broadcast(kscell)
            call broadcast(kdispl)
#endif MPI
! Variable added by MPN to be used by Phonons codes, not used by tbtrans
            if (IOnode) read(iu) istep, ia1


            nullify(lasto)
            allocate(lasto(0:nua))
            if (IOnode) read(iu) lasto
#ifdef MPI
            call broadcast(lasto)
#endif MPI


            if (.not.gamma) then
               nullify(indxuo)
               allocate(indxuo(1:notot))
               if (IOnode) read(iu) (indxuo(ih),ih=1,notot)
#ifdef MPI
               call broadcast(indxuo)
#endif MPI
            end if


           nullify(numh)
           allocate(numh(nuotot))
           if (IOnode) then
              read(iu) numh(1:nuotot)
           end if ! IOnode
#ifdef MPI
           call broadcast(numh)
#endif MPI


           if (IOnode) read(iu) qtot,temp
#ifdef MPI
           call broadcast(qtot)
           call broadcast(temp)
#endif MPI

           allocate(efs(nspin))
           if (IOnode) read(iu) efs(1)
#ifdef MPI
           call broadcast(efs)
#endif MPI
           efs(nspin)=efs(1)

C Create listhptr
           nullify(listhptr)
           allocate(listhptr(nuotot))
           listhptr(1) = 0

           do ih = 2,nuotot
              listhptr(ih) = listhptr(ih-1) + numh(ih-1)
           enddo

C Read listh
C Allocate lish
           nullify(listh)
           allocate(listh(maxnh))

           if (IOnode) then
              do ih = 1,nuotot
                 read(iu) listh(listhptr(ih)+1:listhptr(ih)+numh(ih)) 
              enddo
           end if ! IOnode
#ifdef MPI
           call broadcast(listh)
#endif MPI


           nullify(S)
           allocate(S(maxnh))

           if (IOnode) then
              do ih = 1,nuotot
                 read(iu) S(listhptr(ih)+1:listhptr(ih)+numh(ih))
              enddo
           end if ! IOnode
#ifdef MPI
           call broadcast(S)
#endif MPI


           nullify(H)
           allocate(H(maxnh,nspin))

           if (IOnode) then
              do is = 1,nspin
                 do ih = 1,nuotot
                    read(iu) H(listhptr(ih)+1:listhptr(ih)+numh(ih),is)
                 enddo
              enddo
           end if ! IOnode
#ifdef MPI
           call broadcast(H)
#endif MPI


           if (.not.gamma) then

             nullify(xij)
             allocate(xij(3,maxnh))
             if (IOnode) then
                do ih = 1,nuotot
                   read(iu) 
     . (xij(k,listhptr(ih)+1:listhptr(ih)+numh(ih)),k=1,3)
                enddo
             end if ! IOnode
#ifdef MPI
             call broadcast(xij)
#endif MPI

           endif



            
c Close file
          
          if (IOnode) call io_close( iu )
        else
            if (IOnode) then
              write(6,*) 'iohs: ERROR: file not found: ', finname
              call die('iohs: ERROR: file not found')
            end if
#ifdef MPI
       call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI            
     
        endif
 
       endif

#ifdef MPI
       call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif MPI

      end
