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
      subroutine mulliken(iopt,nspin,natoms,nbasistot,maxnh,numh,
     .                    listhptr,listh,s,dm,isa,lasto,iaorb,
     .                    iphorb)
C ********************************************************************
C Subroutine to perform Mulliken population analysis.
C (Overlap and total populations, both for orbitals and for atoms)
C The density matrix (d.m.) and overlap matrix are passed in sparse form
C (both with the same sparse structure)
C There is no output. The populations are printed to the output.
C
C Written by P.Ordejon, October'96
C Non-collinear spin added by J.M.Soler, May 1998. 
C Symmetry label for each orbital included by DSP, Oct. 1998.
C Label with the principal quantum number introduced by DSP, July 1999.
C ************************** INPUT ************************************
C integer iopt                : Work option: 1 = atomic and orbital charges
C                                            2 = 1 + atomic overlap pop.
C                                            3 = 2 + orbital overlap pop.
C integer nspin               : Number of spin components
C integer natoms              : Number of atoms
C integer nbasistot           : Number of basis orbitals over all nodes
C integer maxnh               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements
C integer numh(nbasis)        : First Control vector of d.m. and overlap
C integer listhptr(nbasis)    : Second Control vector of d.m. and overlap
C integer listh(maxnh)        : Third Control vector of d.m. and overlap
C real*8  s(maxnh)            : Overlap matrix in sparse form
C real*8  dm(maxnh,nspin)     : Density matrix in sparse form 
C integer isa(natoms)         : Species index of each atom
C integer lasto(0:maxa)       : Index of last orbital of each atom
C                               (lasto(0) = 0) 
C integer iaorb(nbasis)       : Atomic index of each orbital 19684
C integer iphorb(nbasis)      : Orbital index of each orbital in its atom
C ************************* OUTPUT *************************************
C No output. The results are printed to standard output
C **********************************************************************
C
C  Modules
C
      use precision,    only : dp
      use parallel,     only : IOnode, Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer
     .  iopt,natoms,nbasistot,maxnh,nspin

      integer
     .  numh(*),lasto(0:natoms),listh(maxnh),listhptr(*),
     .  iphorb(*), isa(natoms), iaorb(*)

      real(dp)
     .  dm(maxnh,nspin),s(maxnh)

      external
     .  memory
C Internal parameters ..................................................
C Number of columns in printout.  Must be smaller than 20
      integer ncol, nbasis
      parameter (ncol = 8)

#ifdef MPI
      integer  :: MPIerror, mpistatus(MPI_STATUS_SIZE)
      real(dp) :: pb(ncol)
      real(dp) :: temp1(5), temp2(5) 
      real(dp), dimension(:), allocatable :: qb
      character(len=200) ::  mpibuff ! 200 is sufficiently long for all buffered writes
#endif

      integer i,ia,ib,ii,imax,in,io,ior,ip,is,j,ja,jja,jo,jor,
     .  ind, nao, nblock, ns, ispec, irow, nrow,
     .  config(ncol), iNode, itot

      real(dp)
     .  p(ncol),qa,qas(4),qts(4),qtot,qtot_temp,stot,svec(3)

      real(dp), dimension(:),   allocatable :: qo
      real(dp), dimension(:,:), allocatable :: qos

      character sym_label(ncol)*7, atm_label*20

      character(len=8), allocatable :: orb_label(:)
      character(len=13) :: atomtitle

C ......................

C Get Node numberReduce
#ifdef MPI
C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(nbasistot,Node,Nodes,nbasis)
#else
      nbasis = nbasistot
#endif

      if (iopt.eq.0) then
C iopt = 0 implies no analysis
        return
      elseif (iopt.lt.0 .or. iopt.gt.3) then
        if (ionode) then
          write(6,"(a)") 'mulliken: ERROR: Wrong iopt'
        endif
        return
      endif 

C Allocate local memory
      allocate(qos(nspin,nbasis))
      call memory('A','D',nspin*nbasis,'mulliken')

      ns=0
      do i = 1,natoms
        ns=max(ns,isa(i))
      enddo 

C Compute and print Overlap Populations for Orbitals ....................
      if (iopt .eq. 3) then
        if (ionode) then
          write(6,*) 
          write(6,"(a)")
     .         'mulliken: Mulliken Overlap Populations between Orbitals'
        endif
        nblock = nbasistot / ncol
        ip=1
        if (nblock*ncol .eq. nbasistot) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = nbasistot - nblock * ncol  
          do ii=1,imax 
             sym_label(ii)=symfio(isa(iaorb((ib-1)*ncol+ii)),
     .                iphorb((ib-1)*ncol+ii))
             config(ii)=cnfigfio(isa(iaorb((ib-1)*ncol+ii)),
     .                iphorb((ib-1)*ncol+ii))
          enddo 
          if (ionode) then
            write(6,*) 
            write(6,'(14x,20(2x,i4,3x))')((ib-1)*ncol+ii,ii=1,imax) 
            write(6,'(17x,20(1x,i1,a7))') 
     .              (config(ii),sym_label(ii),ii=1,imax)
          endif

! Following algorithm is designed to ensure that terms relating
! to orbitals are printed out in the correct sequence. It's even
! less efficient than the algorithm that used to be here - on the
! other hand, it does actually work properly.
! (Briefly: iterate over orbitals - 
!   if it's on the current node, and we are node 0 then print immediately; 
!   else if it's not, print the output to a buffer & send it to node 0.
!   if it's not on the current node, then if we're node 0, wait for the
!   message from whichever node it is on, and print the buffer. 
!   Otherwise do nothing.
          do itot = 1,nbasistot
            call GlobalToLocalOrb(itot,Node,Nodes,i)
            if (i.gt.0) then
              sym_label(1)=symfio(isa(iaorb(itot)),iphorb(itot))
              config(1)=cnfigfio(isa(iaorb(itot)),iphorb(itot))
              do ii = 1,imax
                p(ii) = 0._dp
              enddo
              do in = 1,numh(i)
                ind = listhptr(i)+in
                j = listh(ind)
                ii = j - (ib - 1) * ncol
                if (ii .ge. 1 .and. ii .le. imax) then
                  p(ii) = 0._dp
                  do is = 1,min(nspin,2)
                    p(ii) = p(ii) + dm(ind,is) * s(ind)
                  enddo
                endif
              enddo
              if (ionode) then
                write(6,15) itot,config(1), sym_label(1),
     .               (p(ii),ii=1,imax)
#ifdef MPI
              else
                write(mpibuff,15) itot,config(1), sym_label(1),
     .                (p(ii),ii=1,imax)
                call mpi_send(mpibuff, 200, MPI_Character, 0, itot, 
     .             MPI_Comm_world, mpierror)
              endif
            else
              if (ionode) then
                call mpi_recv(mpibuff, 200, MPI_Character, 
     .                MPI_ANY_SOURCE, itot,
     .                MPI_Comm_world, mpistatus, mpierror)
                write(6,'(a)') trim(mpibuff)
#endif
              endif
            endif
          enddo
        enddo
      endif
C ...................

C Compute and print Overlap Populations for Atoms ....................
      if (iopt .ge. 2) then
        if (ionode) then
          write(6,*) 
          write(6,"(a)")
     .          'mulliken: Mulliken Overlap Populations between Atoms'
        endif
        nblock = natoms / ncol
        ip=1
        if (nblock*ncol .eq. natoms) ip=0
        do ib = 1,nblock+ip
          imax = ncol
          if (ib .eq. nblock+1) imax = natoms - nblock * ncol
          if (ionode) then
            write(6,*) 
            write(6,10) ((ib-1)*ncol+ii,ii=1,imax)
          endif
          do i = 1,natoms
            do ii = 1,imax
              p(ii) = 0.0
            enddo
            do ior = lasto(i-1)+1,lasto(i)
              call GlobalToLocalOrb(ior,Node,Nodes,itot)
              if (itot.gt.0) then
                do in = 1,numh(itot)
                  ind = listhptr(itot)+in
                  jor = listh(ind)
                  do jja = 1,natoms
                    if (lasto(jja) .ge. jor) then
                      ja = jja
                      goto 100
                    endif
                  enddo
                  goto 110
 100              ii = ja - (ib - 1) * ncol
                  if (ii .ge. 1 .and. ii .le. imax) then
                    do is = 1,min(nspin,2)
                      p(ii) = p(ii) + dm(ind,is) * s(ind)
                    enddo
                  endif
 110              continue
                enddo
              endif
            enddo

#ifdef MPI
C Global sum of values stored in p
            pb(1:imax)=p(1:imax)
            call MPI_Reduce(pb,p,imax,MPI_double_precision,MPI_sum,
     .        0,MPI_Comm_World,MPIerror)
#endif

            if (ionode) then
              write(6,11) i,(p(ii),ii=1,imax)
            endif
          enddo
        enddo
      endif
C ....................

C Compute and print Mulliken Orbital and Atomic Populations ..........
      if (iopt .ge. 1) then
        if (ionode) then
          write(6,*) 
          write(6,"(a)")
     .         'mulliken: Mulliken Atomic and Orbital Populations:'
        endif
        if (nspin .le. 2) then
C We only write mulliken to cml file if we are unpolarized or collinear 
C and (currently) don't bother in the non-collinear case.
          if (cml_p) then
             call cmlStartPropertyList(xf=mainXML, 
     .            dictRef='siesta:mulliken',
     .            title='Mulliken Population Analysis')
          endif
          do is = 1,nspin
            if (nspin .eq. 2) then
              if(is .eq. 1) then
                if (ionode) write(6,'(/,a)') 'mulliken: Spin UP '
                if (cml_p) call cmlStartPropertyList(xf=mainXML, 
     .               title='SpinUp')
              elseif(is .eq. 2) then
                if (ionode) write(6,'(/,a)') 'mulliken: Spin DOWN '
                if (cml_p) call cmlStartPropertyList(xf=mainXML, 
     .               title='SpinDown')
              endif
            else
              if (cml_p) call cmlStartPropertyList(xf=mainXML, 
     .               title='Unpolarized')
            endif
            qtot = 0._dp
            do ispec = 1,ns  

             atm_label = labelfis(ispec)
             if (ionode) then
               write(6,'(/2a)')'Species: ', atm_label 
               write(6,'(a4,a7,a6)') 'Atom', 'Qatom', 'Qorb'
             endif
C DSP, Writing symmetries for each orbital. 
C DSP, Orbitals with a 'P' belong to the polarization shell
             nao = nofis(ispec)
             allocate(orb_label(nao))
             do io=1,nao
               write(orb_label(io),'(i1,a7)')
     .              cnfigfio(ispec,io), symfio(ispec,io)
             enddo
             nrow=nao/ncol
             io=1
             do irow=1,nrow 
               if (ionode) then
                 write(6,'(15x,8(a8))') orb_label(io:io+ncol-1)
                 io = io+ncol
               endif
             enddo 
             if (ionode) then
               write(6,'(15x,8(a8))') orb_label(io:nao)
             endif
             if (cml_p) then
               call cmlStartPropertyList(xf=mainXML,
     .              title=trim(atm_label))
               call cmlAddProperty(xf=mainXML,
     .              title='Orbital list', value=orb_label)
             endif
             deallocate(orb_label)
             
             do ia = 1,natoms
               if (isa(ia).eq.ispec) then
C             Compute charge in each orbital of atom ia
                 nao = lasto(ia) - lasto(ia-1)
                 allocate(qo(nao))
                 qa = 0.0_dp
                 do io = lasto(ia-1)+1,lasto(ia)
                   nao = io - lasto(ia-1)
                   qo(nao) = 0.0_dp
                   call GlobalToLocalOrb(io,Node,Nodes,itot)
                   if (itot.gt.0) then
                     do in = 1,numh(itot)
                       ind = listhptr(itot) + in
                       qo(nao) = qo(nao) + dm(ind,is) * s(ind)
                     enddo
                     qa = qa + qo(nao)
                   endif
                 enddo
                 qtot = qtot + qa
#ifdef MPI
C Global sum of values stored in qo
                 allocate(qb(nao))
                 qb = qo
                 call MPI_Reduce(qb,qo,nao,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
                 qb(1)=qa
                 call MPI_Reduce(qb(1),qa,1,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
                 deallocate(qb)
#endif
                 if (ionode) then
                   write(6,'(i4,f7.3,8f8.3,(/11x,8f8.3))')
     .               ia, qa, qo 
                 endif
                 if (cml_p) then
                   write(atomtitle,'(a,i8)') 'atom ',ia
                   call cmlStartPropertyList(xf=mainXML,
     .                  title=trim(atomtitle))
                   call cmlAddProperty(xf=mainXML,
     .                  title="Atomic Charge", value=qa, 
     .                  units="siestaUnits:e")
                   call cmlAddProperty(xf=mainXML,
     .                  title="Orbital Charges", value=qo,
     .                  units="siestaUnits:e")
                   call cmlEndPropertyList(xf=mainXML) !atomtitle
                 endif 
                 deallocate(qo)
               endif
             enddo
             if (cml_p) call cmlEndPropertyList(xf=mainXML) !atm_label
            enddo

#ifdef MPI
C Global sum of total charge
            qtot_temp = qtot
            call MPI_Reduce(qtot_temp,qtot,1,MPI_double_precision,
     $           MPI_sum, 0,MPI_Comm_World,MPIerror)
#endif
            if (ionode) then
              write(6,"(/a,f12.3)") 'mulliken: Qtot = ', qtot
            endif
            if (cml_p) then
              call cmlAddProperty(xf=mainXML,
     .             title="Total Charge", value=qtot, 
     .             units="siestaUnits:e")
              call cmlEndPropertyList(xf=mainXML) ! spins
            endif
          enddo
          if (cml_p) call cmlEndPropertyList(xf=mainXML) ! mulliken

        elseif (nspin .eq. 4) then
          do is = 1,nspin
            qts(is) = 0.0_dp
            do io = 1,nbasis
              qos(is,io) = 0.0_dp
                enddo
          enddo
          do is = 1,nspin
            do io = 1,nbasis
              do in = 1,numh(io)
                ind = listhptr(io)+in
                jo = listh(ind)
                qos(is,io) = qos(is,io) + dm(ind,is)*s(ind)/2
                if (jo.le.nbasis)
     .            qos(is,jo) = qos(is,jo) + dm(ind,is)*s(ind)/2
              enddo
            enddo
          enddo 

          do ispec=1,ns
            atm_label=labelfis(ispec)
            if (ionode) then
              write(6,'(/2a)')'Species: ', atm_label
              write(6,'(/,a4,a9,4x,2a10,3x,a8,/,a)')
     .             'Atom', 'Orb', 'Charge', 'Spin', 'Svec',
     .             repeat('-', 64)
            endif

            do ia = 1,natoms 
              if (isa(ia).eq.ispec) then 
                qas(1:nspin) = 0._dp
                do iNode = 0,Nodes-1
                  if (iNode.eq.Node) then
                    do io = lasto(ia-1)+1,lasto(ia)  
                      call GlobalToLocalOrb(io,Node,Nodes,itot)
                      if (itot.gt.0) then
C DSP, Writing symmetries for each orbital.
C DSP, Orbitals with a 'P' belong to the polarization shell
                        sym_label(1)=symfio(ispec,iphorb(io))
                        config(1)=cnfigfio(ispec,iphorb(io))
                        do is = 1,nspin
                          qas(is) = qas(is) + qos(is,itot)
                          qts(is) = qts(is) + qos(is,itot)
                        enddo
                        call spnvec( nspin, qos(1,itot), qtot,stot,svec)
                        if (ionode) then
                          write(6,'(i4,i5,i1,a7,2f10.5,3x,3f8.3)')
     .                         ia, io,config(1),sym_label(1), 
     .                         qtot, stot, svec
#ifdef MPI
                        else
                          write(mpibuff,'(i4,i5,i1,a7,2f10.5,3x,3f8.3)')
     .                         ia, io,config(1),sym_label(1), 
     .                         qtot, stot, svec
                          call mpi_send(mpibuff, 200, MPI_Character, 
     .                         0, io, MPI_Comm_world, mpierror)
                        endif
                      else
                        if (ionode) then
                          call mpi_recv(mpibuff, 200, MPI_Character, 
     .                          MPI_ANY_SOURCE, io,
     .                          MPI_Comm_world, mpistatus, mpierror)
                          write(6,'(a)')  trim(mpibuff)
#endif
                        endif
                      endif
                    enddo
                  endif
                enddo
                call spnvec( nspin, qas, qtot, stot, svec )
#ifdef MPI
C Global reduction of terms
                temp1(1:3) = svec
                temp1(4)   = qtot
                temp1(5)   = stot
                call MPI_Reduce(temp1,temp2,5,MPI_double_precision,
     $               MPI_sum,0, MPI_Comm_World,MPIerror)
                svec = temp2(1:3)
                qtot = temp2(4)
                stot = temp2(5)
#endif
                if (ionode) then
                  write(6,'(i4,4x,a6,3x,2f10.5,3x,3f8.3,/)')
     .                 ia, 'Total', qtot, stot, svec 
                endif
              endif
            enddo
          enddo 
          call spnvec( nspin, qts, qtot, stot, svec )
#ifdef MPI
C Global reduction of terms
          temp1(1:3) = svec
          temp1(4)   = qtot
          temp1(5)   = stot
          call MPI_Reduce(temp1,temp2,5,MPI_double_precision,
     $         MPI_sum,0,MPI_Comm_World,MPIerror)
          svec = temp2(1:3)
          qtot = temp2(4)
          stot = temp2(5)
#endif
          if (ionode) then
            write(6,'(a,/,2a8,1x,2f10.5,3x,3f8.3,/)')
     .        repeat('-',64), 'Total', 'Total', qtot, stot, svec
          endif
        endif
      endif

C Deallocate local memory
      call memory('D','D',size(qos),'mulliken')
      deallocate(qos)
C ...................

10    format(12x,20(2x,i4,2x))
11    format(i12,20(1x,f7.3))
15    format(i4,1x,i1,a7,20(1x,f8.3))
      return
      end



      subroutine spnvec( ns, qs, qt, st, sv )
c ********************************************************************
c Finds the spin vector components from the spin density matrix
c Written by J.M.Soler, May 1998.
c ******* Input ******************************************************
c integer ns     : Number of components in spin density matrix
c real*8  qs(ns) : Spin density matrix elements with the convention
c                  is=1 => Q11; is=2 => Q22; is=3 => Real(Q12);
c                  is=4 => Imag(Q12)
c ******* Output *****************************************************
c real*8  qt    : Total charge
c real*8  st    : Total spin
c real*8  sv(3) : Spin vector
c ********************************************************************

      use precision

      implicit          none
      integer           ns
      real(dp)          qs(ns), qt, st, sv(3)
      real(dp)          cosph, costh, sinph, sinth, tiny
      parameter ( tiny = 1.e-12_dp )

      if (ns .eq. 1) then
        qt = qs(1)
        st = 0.0_dp
        sv(1) = 0.0_dp
        sv(2) = 0.0_dp
        sv(3) = 0.0_dp
      elseif (ns .eq. 2) then
        qt = qs(1) + qs(2)
        st = qs(1) - qs(2)
        sv(1) = 0.0_dp
        sv(2) = 0.0_dp
        sv(3) = st
      elseif (ns .eq. 4) then
        qt = qs(1) + qs(2)
        st = sqrt( (qs(1)-qs(2))**2 + 4._dp*(qs(3)**2+qs(4)**2) )
        costh = ( qs(1) - qs(2) ) / ( st + tiny )
        sinth = sqrt( 1.0_dp - costh**2 )
        cosph =  qs(3) / ( sqrt( qs(3)**2 + qs(4)**2 ) + tiny )
        sinph = -qs(4) / ( sqrt( qs(3)**2 + qs(4)**2 ) + tiny )
        sv(1) = st * sinth * cosph
        sv(2) = st * sinth * sinph
        sv(3) = st * costh
      else
        write(6,*) 'spnvec: ERROR: invalid argument ns =', ns
        return
      endif
      end
