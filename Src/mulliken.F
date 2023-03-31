! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine mulliken(iopt,natoms,nbasistot,maxnh,numh,
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
C integer natoms              : Number of atoms
C integer nbasistot           : Number of basis orbitals over all nodes
C integer maxnh               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements
C integer numh(nbasis)        : First Control vector of d.m. and overlap
C integer listhptr(nbasis)    : Second Control vector of d.m. and overlap
C integer listh(maxnh)        : Third Control vector of d.m. and overlap
C real*8  s(maxnh)            : Overlap matrix in sparse form
C real*8  dm(maxnh,spin%DM)   : Density matrix in sparse form 
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
      use parallelsubs, only : LocalToGlobalOrb
      use atmfuncs,     only : symfio, cnfigfio, labelfis, nofis
      use alloc,        only : re_alloc, de_alloc
      use siesta_cml
      use m_spin,       only : spin
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in) :: 
     .  iopt,natoms,nbasistot,maxnh

      integer
     .  numh(*),lasto(0:natoms),listh(maxnh),listhptr(*),
     .  iphorb(*), isa(natoms), iaorb(*)

      real(dp)
     .  dm(maxnh,spin%DM),s(maxnh)

C Internal parameters ..................................................
C Number of columns in printout.  Must be smaller than 20
      integer ncol, nbasis
      parameter (ncol = 8)

#ifdef MPI
      integer  :: MPIerror, mpistatus(MPI_STATUS_SIZE)
      real(dp) :: pb(ncol)
      real(dp),  pointer :: qb(:)=>null()
      character(len=128) ::  mpibuff ! 128 is sufficiently long
#endif

      integer i,ia,ib,ii,imax,in,io,ior,ip,is,j,ja,jja,jo,jor,
     .  ind, nao, nblock, ns, ispec, irow, nrow,
     .  config(ncol), itot

      real(dp)
     .     p(ncol),qa,qas(4),qts(4),qtot,qtot_temp,stot,svec(3)
      real(dp) :: dm_aux(4)

      real(dp), pointer :: qo(:)=>null(), qos(:,:)=>null()

      character sym_label(ncol)*7, atm_label*20

      character(len=8), pointer :: orb_label(:)=>null()
      character(len=13) :: atomtitle

C ......................
#ifdef DEBUG
      call write_debug( '    PRE mulliken' )
#endif

C Get Node numberReduce
#ifdef MPI
C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(nbasistot,Node,Nodes,nbasis)
#else
      nbasis = nbasistot
#endif

      if (iopt.eq.0) then
C iopt = 0 implies no analysis
#ifdef DEBUG
      call write_debug( '    POS mulliken iopt==0' )
#endif
        return
      elseif (iopt.lt.0 .or. iopt.gt.3) then
        if (ionode) then
          write(6,"(a)") 'mulliken: ERROR: Wrong iopt'
        endif
#ifdef DEBUG
      call write_debug( '    POS mulliken iopt<0 or iopt>3' )
#endif
        return
      endif 


      ns=0
      do i = 1,natoms
        ns=max(ns,isa(i))
      enddo 

C Compute and print Overlap Populations for Orbitals ....................
      if (iopt .eq. 3) then
        if (ionode) then
          write(6,*) 
          write(6,"(a)")'mulliken: Overlap Populations between Orbitals'
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
                  do is = 1,spin%Spinor
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
                call mpi_send(mpibuff, 128, MPI_Character, 0, itot, 
     .             MPI_Comm_world, mpierror)
              endif
            else
              if (ionode) then
                call mpi_recv(mpibuff, 128, MPI_Character, 
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
          write(6,"(a)")'mulliken: Overlap Populations between Atoms'
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
                    do is = 1,spin%Spinor
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
          write(6,"(a)")'mulliken: Atomic and Orbital Populations:'
        endif
        if (spin%DM .le. 2) then
C We only write mulliken to cml file if we are unpolarized or collinear 
C and (currently) don't bother in the non-collinear case.
          if (cml_p) then
             call cmlStartPropertyList(xf=mainXML, 
     .            dictRef='siesta:mulliken',
     .            title='Mulliken Population Analysis')
          endif
          do is = 1,spin%DM
            if (spin%DM .eq. 2) then
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
             call re_alloc( orb_label, 1, nao, 'orb_label', 'mulliken' )
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
     $              dictref='siesta:mull_orb_list',
     .              title='Orbital list', value=orb_label)
             endif
             call de_alloc( orb_label, 'orb_label', 'mulliken' )

             do ia = 1,natoms
               if (isa(ia).eq.ispec) then
C             Compute charge in each orbital of atom ia
                 nao = lasto(ia) - lasto(ia-1)
                 call re_alloc( qo, 1, nao, 'qo', 'mulliken' )
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
                 call re_alloc( qb, 1, nao, 'qb', 'mulliken' )
                 qb = qo
                 call MPI_Reduce(qb,qo,nao,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
                 qb(1)=qa
                 call MPI_Reduce(qb(1),qa,1,MPI_double_precision,
     .             MPI_sum,0,MPI_Comm_World,MPIerror)
                 call de_alloc( qb, 'qb', 'mulliken' )
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
     $                  dictref='siesta:mull_atm_charge',
     .                  title="Atomic Charge", value=qa, 
     .                  units="siestaUnits:e")
                   call cmlAddProperty(xf=mainXML,
     $                  dictref='siesta:mull_orb_charge',
     .                  title="Orbital Charges", value=qo,
     .                  units="siestaUnits:e")
                   call cmlEndPropertyList(xf=mainXML) !atomtitle
                 endif 
                 call de_alloc( qo, 'qo', 'mulliken' )
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
     $             dictref='siesta:mull_total_charge',
     .             title="Total Charge", value=qtot, 
     .             units="siestaUnits:e")
              call cmlEndPropertyList(xf=mainXML) ! spins
            endif
          enddo
          if (cml_p) call cmlEndPropertyList(xf=mainXML) ! mulliken

        elseif (spin%DM .ge. 4) then
          call re_alloc( qos, 1, spin%Grid, 1, nbasistot,
     &          'qos', 'mulliken' )
          do is = 1,spin%Grid
            qts(is) = 0.0_dp
            do io = 1,nbasistot
              qos(is,io) = 0.0_dp
            enddo
          enddo
          do io = 1,nbasis
             call LocalToGlobalOrb(io,Node, Nodes, itot)
             do in = 1,numh(io)
                ind = listhptr(io)+in
                jo = listh(ind)
                dm_aux(1:4) = dm(ind,1:4)
                if (spin%SO) then
                   dm_aux(3) = 0.5_dp*(dm_aux(3)+dm(ind,7))
                   dm_aux(4) = 0.5_dp*(dm_aux(4)+dm(ind,8))
                endif
                do is = 1,spin%Grid
                   qos(is,itot) = qos(is,itot) + dm_aux(is)*s(ind)
                enddo
             enddo
          enddo 
#ifdef MPI
          call re_alloc( qb, 1, spin%Grid*nbasistot,
     &         'qb', 'mulliken' )
          call MPI_Reduce(qos(1,1), qb(1), spin%Grid*nbasistot,
     &         Mpi_Double_Precision, MPI_Sum, 
     &         0, MPI_Comm_world, MPIerror)
          itot = 0
          do io = 1 , nbasistot
             do is = 1 , spin%Grid
                itot = itot + 1
                qos(is, io) = qb(itot)
             end do
          end do 
          call de_alloc(qb, 'qb', 'mulliken')
#endif

          ! The above reduction means that the
          ! following mulliken calculation actually is
          ! serial. 
          ! We could have performed other analysis, but
          ! that would prove additional "small" communication.
          ! This should aptly provide speed and flexibility.
          ! Also, this entire loop may be encapsulated
          ! In (if ionode)... However, this will also
          ! streamline the processors.
          do ispec=1,ns
            atm_label=labelfis(ispec)
            if (ionode) then
              write(6,'(/2a)')'Species: ', atm_label
              write(6,'(/,a4,1x,a8,4x,2a10,3x,a8,/,a)')
     .             'Atom', 'Orb', 'Charge', 'Spin', 'Svec',
     .             repeat('-', 64)
            endif

            do ia = 1,natoms
               
              ! Check for same species
              if ( isa(ia) /= ispec ) cycle
                  
              ! Initialize atomic charges (per spin)
              qas(1:spin%Grid) = 0._dp
              do io = lasto(ia-1)+1, lasto(ia)
                
C DSP, Writing symmetries for each orbital.
C DSP, Orbitals with a 'P' belong to the polarization shell
                 sym_label(1) = symfio(ispec,iphorb(io))
                 config(1) =  cnfigfio(ispec,iphorb(io))

                 ! Reduce total charge per-atom
                 do is = 1 , spin%Grid
                    qas(is) = qas(is) + qos(is,io)
                 enddo
                 
                 call spnvec( spin%Grid, qos(1:spin%Grid,io), 
     $                qtot, stot, svec)
                 if ( IONode ) then
                    write(6,'(i5,i3,tr1,i1,a7,2f10.5,3x,3f8.3)')
     .                   ia,io-lasto(ia-1),config(1),sym_label(1), 
     .                   qtot, stot, svec
                 end if

              end do ! atomic orbitals

              call spnvec( spin%Grid, qas, qtot, stot, svec )
              if ( IONode ) then
                 write(6,'(i5,5x,a5,2x,2f10.5,3x,3f8.3,/)')
     .                ia, 'Total', qtot, stot, svec 
              end if

              ! Reduce total charge in system
              do is = 1 , spin%Grid
                 qts(is) = qts(is) + qas(is)
              end do

            end do ! loop atoms

            call spnvec( spin%Grid, qts, qtot, stot, svec )
            if ( IONode ) then
               write(6,'(a,/,a5,12x,2f10.5,3x,3f8.3,/)')
     .              repeat('-',64), 'Total', qtot, stot, svec
            end if

          end do ! species

          call de_alloc( qos, 'qos', 'mulliken' )
          
        end if ! spin%DM .ge. 4
      end if

C ...................

#ifdef DEBUG
      call write_debug( '    POS mulliken' )
#endif

10    format(12x,20(2x,i4,2x))
11    format(i12,20(1x,f7.3))
15    format(i4,1x,i1,a7,20(1x,f8.3))
      return
      end subroutine mulliken

      subroutine spnvec( ns, qs, qt, st, sv )
c ********************************************************************
c Finds the spin vector components from the spin density matrix
c Written by J.M.Soler, May 1998.
c ******* Input ******************************************************
c integer ns     : Number of components in spin density matrix
c real*8  qs(ns) : Spin density matrix elements with the convention
c                  is=1 => Q11; is=2 => Q22; is=3 => Real(Q12);
c                  is=4 => -Imag(Q12)  (*Note 
c ******* Output *****************************************************
c real*8  qt    : Total charge
c real*8  st    : Total spin moment
c real*8  sv(3) : Spin moment vector
c ********************************************************************
!
! Spin magnetic **moments** in bohr magnetons:
! m_x = 2 Re{Q12} = 2*qs(3)
! m_y = -2 Im{Q12} = 2*qs(4)
! m_z = Q11-Q22 = qs(1)-qs(2)
!
! The spin moment for a spin angular momentum of 1/2 (hbar) is (very close to) 1 bohr magneton.
!
      use precision

      implicit          none
      integer, intent(in)   ::          ns
      real(dp), intent(in)  ::          qs(ns)
      real(dp), intent(out) ::          qt, st, sv(3)

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
         sv(1) = 2.0_dp * qs(3)
         sv(2) = 2.0_dp * qs(4)
         sv(3) = qs(1) - qs(2)
         st = sqrt(sv(1)**2+sv(2)**2+sv(3)**2)
      else
        write(6,*) 'spnvec: ERROR: invalid argument ns =', ns
        return
      endif
      end
