! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_ksv
      public :: KSV_pol
      CONTAINS
      subroutine KSV_pol( nua, na, xa, rmaxo, scell, ucell, nuotot,
     .                    nuo, no, nspin, qspin, maxnh, 
     .                    maxkpol, numh, listhptr, listh, H, S, 
     .                    xijo, indxuo, isa, iphorb, iaorb, 
     .                    lasto, shape, nkpol,kpol,
     .                    wgthpol, polR, polxyz )
C *********************************************************************
C Finds polarization using the method of King-Smith and Vanderbilt
C ( Geometric Berry phase).
C Written by DSP, March 1999.
C **************************** INPUT **********************************
C integer nua                 : Number of atoms in the unit cell
C integer na                  : Number of atoms 
C real*8  xa(3,na)            : Atomic positions in cartesian coordinates
C real*8  rmaxo               : Maximum cutoff for atomic orbitals
C real*8  scell(3,3)          : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  ucell(3,3)          : Unit cell vectors
C integer nuotot              : No. of basis orbitals in the unit cell (global)
C integer nuo                 : No. of basis orbitals in the unit cell (local)
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C real*8  qspin(2)            : Total population of spin up and down
C integer maxnh               : Maximum number of orbitals interacting  
C                               with any orbital
C integer maxk                : Last dimension of kpoint 
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row of
C                               hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C real*8  H(maxnh,spin)       : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xijo(3,maxnh)       : Vectors between orbital centers (sparse)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C integer isa(na)             : Species index for each atom
C integer iphorb(no)          : Orbital index within the atom for each
C                               orbital
C integer lasto(0:na)      : Last orbital index of each atom 
C *************************** INPUT/OUTPUT ****************************
C integer  nkpol              :  Maximum number of grid points for the
C                                bidimensional integrals
C *************************** OUTPUT **********************************
C real*8   kpol(3,maxkpol)   : Auxiliar array to store the kpoints
C                              for the bidimensional integrations.
C real*8   wgthpol(maxkpol)  : Auxiliar array to store the weigths of
C                              the kpoints.            
C real*8   polR(3,nspin)     : Macroscopic polarization per unit cell
C                              along the vectors which define the unit 
C                              cell in real space.
C real*8   polxyz(3,nspin)   : Macroscopic polarization per unit cell
C                              along each of the cartesian coordinates.
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************
C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : IOnode
      use sys,           only : die
      use atmfuncs,      only : zvalfis
      use densematrix,   only : allocDenseMatrix, resetDenseMatrix
      use densematrix,   only : Haux, Saux, psi
      use alloc,         only : re_alloc, de_alloc
      USE m_ksvinit,     only : repol

      implicit          none

      integer           maxkpol, maxnh, nuo, nuotot, no, nspin, na
      integer           indxuo(no), listh(maxnh), numh(nuo), nkpol
      integer           listhptr(nuo), isa(na), iphorb(no), iaorb(no) 
      integer           nua, lasto(0:na)
      real(dp)          ddot, 
     .                  H(maxnh,nspin), kpol(3,maxkpol), 
     .                  S(maxnh), xijo(3,maxnh),
     .                  polR(3,nspin), polxyz(3,nspin), ucell(3,3),
     .                  wgthpol(maxkpol), 
     .                  scell(3,3), rmaxo,  xa(3,na)
C *********************************************************************
C Internal variables 
      real(dp)  ntote_real   !!** Needed to deal with synthetics
      integer
     .  ik, il, io, ispin, iuo, ix, i,
     .  iy, npl, jo, is, ia, ntote, nocc(2),
     .  nk2D(3), kscell(3,3), igrd, 
     .  nk, nmeshk(3,3), Nptot, 
     .  notcal(3), nhs, npsi

      real(dp)
     .  difA, pi, rcell(3,3), uR(3,3),  
     .  displ(3), dsp(3), cutoff, dk(3), detr, deti, 
     .  prodtr, prodti, kint(3), volcel, dkmod,
     .  dkxij, ckxij, skxij, dmod, polion(3),
     .  tiny, phase, ph(3,2), Debye,
     .  vaux(3,2), area, J, qspin(2), dq, phaseold(2)

      real(dp), dimension(:), pointer ::  ek => null()

      parameter (Debye  = 0.393430d0)  

      character         shape*10

      external          ddot, volcel, reclat, memory

      integer, dimension(:), pointer ::  muo => null()
      real(dp), dimension(:), pointer :: psi1 => null()
      real(dp), dimension(:), pointer :: psiprev => null()
      real(dp), dimension(:), pointer :: aux => null()
      
      parameter (  tiny= 1.0d-8  )

C Start time counter 
      call timer( 'KSV_pol', 1 )

!! jjunquer
!      write(6,*)' Node, Nodes = ', Node, Nodes
!! end jjunquer

C Reading unit cell and calculate the reciprocal cell
      call reclat( ucell, rcell, 1 )

C Find the integration grids in reciprocal space   
      call repol(nmeshk,dsp) 
      nk=nmeshk(1,1)*nmeshk(2,1)*nmeshk(3,1)+
     .   nmeshk(1,2)*nmeshk(2,2)*nmeshk(3,2)+
     .   nmeshk(1,3)*nmeshk(2,3)*nmeshk(3,3)

C Exit if no calculation is going to be performed
      if (nk.eq.0) goto 999

      if (nspin.eq.1) then 
C Total number of valence electrons in the unit cell
        ntote_real=0.0_dp
        do ia=1,nua
          is=isa(ia) 
          ntote_real=ntote_real+zvalfis(is)
        enddo
        if ((nint(ntote_real) - ntote_real) .gt. 1.0e-6_dp) then
           if (IOnode)
     .       write(6,'(/,a,/,a,/a)')
     .      'KSV_pol: Non-integer number of electrons',
     .      'KSV_pol: This is hardly an insulator',
     .      'KSV_pol: No polarization calculation performed' 
           goto 999
        endif
        ntote = nint(ntote_real)
        if (mod(ntote,2).ne.0) then 
           if (IOnode) then
              write(6,'(/,a,/,a,/a)')
     .      'KSV_pol: Odd total number of electrons',
     .      'KSV_pol: This is hardly an insulator',
     .      'KSV_pol: No polarization calculation performed' 
           endif
          goto 999
        else 
C Number of occupied bands for each spin component
          nocc(1)=ntote/2
        endif
c        write(6,'(/a,i5)') 
c     .    'KSV_pol: Total number of electrons ',ntote
c         write(6,'(a,i5)')
c     .    'KSV_pol: Number of occupied bands ', nocc(1)
      else
        dq=dabs(qspin(1)-qspin(2))
        if (dabs(dq-nint(dq)).gt.1.0d-2) then 
           if (IOnode) then
              write(6,'(/,a,/,a,/a)')
     .    'KSV_pol: Spin polarization should have an integer value',
     .    'KSV_pol: This is not an insulator for both spin components',
     .    'KSV_pol: No polarization calculation performed'
           endif
          goto 999
        endif 
        nocc(1)=nint(qspin(1))
        nocc(2)=nint(qspin(2))
      endif

C Check parameter maxkpol 
      if (nkpol .gt. maxkpol) then
         if (IONode) then
            write(6,'(/,a,/,a)')
     .    'KSV_pol: WARNING: parameter maxkpol too small',
     .    'KSV_pol: No polarization calculation performed'
         endif
        goto 999
      endif

C Allocate local memory
      nhs = 2*nuotot*nuo
      npsi = 2*nuotot*nuo
      call allocDenseMatrix(nhs, nhs, npsi)

      nullify( muo, ek, psi1, psiprev, aux )
      call re_alloc( muo,     1, nuotot, 'muo',     'KSV_pol' )
      call re_alloc( ek,      1, nuotot, 'ek',      'KSV_pol' )
      call re_alloc( psi1,    1, npsi,   'psi1',    'KSV_pol' )
      call re_alloc( psiprev, 1, npsi,   'psiprev', 'KSV_pol' )
      call re_alloc( aux, 1 , maxnh, 'aux', 'KSV_pol' )

C Initialise psi
      do io = 1,npsi
        psi(io) = 0.0d0
        psi1(io) = 0.0d0
        psiprev(io) = 0.0d0
      enddo

C Check indxuo 
      do iuo = 1,nuotot
        muo(iuo) = 0
      enddo
      do io = 1,no
        iuo = indxuo(io)
        if (indxuo(io).le.0 .or. indxuo(io).gt.no) then
           if (IOnode) then
              write(6,*) 'KSV_pol: invalid index: io, indxuo =',
     .             io, indxuo(io)
           endif
          call die('KSV_pol: invalid indxuo')
        endif
        muo(iuo) = muo(iuo) + 1
      enddo
      do iuo = 1,nuo
        if (muo(iuo) .ne. muo(1)) then
           if (IOnode)
     .        write(6,'(/,2a,3i6)') 
     .          'KSV_pol: ERROR: inconsistent indxuo.',
     .          ' iuo, muo(iuo), muo(1) =', iuo, muo(iuo), muo(1)
          call die('KSV_pol: ERROR: inconsistent indxuo.')
        endif
      enddo
      
C******************PI*************************************
      pi=dacos(-1.0d0)

C  Calculation  of the three components of the macroscopic
C  polarization

C Nuclear component

      do igrd=1,3
        dmod=ddot(3,ucell(1,igrd),1,ucell(1,igrd),1)
        dmod=dsqrt(dmod)
        do ix=1,3
          uR(ix,igrd)=ucell(ix,igrd)/dmod
        enddo
      enddo
   
      do ix = 1,3 
        do ispin = 1,nspin
          polxyz(ix,ispin) = 0.0d0
          polR(ix,ispin) = 0.0d0
          ph(ix,ispin) = 0.0d0
        enddo 
        polion(ix) = 0.0d0
        do ia = 1,nua
          is = isa(ia)
          dmod = ddot(3,rcell(1,ix),1,xa(1,ia),1)*zvalfis(is)/nspin
          polion(ix) = polion(ix) + 
     .     dmod*dsqrt(ddot(3,ucell(1,ix),1,ucell(1,ix),1))/(2.0d0*pi) 
        enddo 
      enddo 
    
      do igrd = 1,3 
        notcal(igrd) = 0
        Nptot = nmeshk(1,igrd)*nmeshk(2,igrd)*nmeshk(3,igrd)

        if ( Nptot.gt. 0 ) then

C Obtain the points and weigths for the bidimensional integration
          cutoff=0.0d0
          do ix=1,3
            do iy= 1,3
              kscell(ix,iy)=0
            enddo
          enddo
          nk2D(igrd)=1
          do ix=1,3
            if (ix.ne.igrd) then
              kscell(ix,ix)= nmeshk(ix,igrd)
              displ(ix)=dsp(igrd)
              nk2D(igrd)=nk2D(igrd)*nmeshk(ix,igrd)
            else
              kscell(ix,ix)=1
              displ(ix)=0.0d0
            endif
          enddo
          nk=nkpol
          call kgridinit( ucell, kscell, displ, cutoff, nk )
          call kgrid( ucell, kscell, displ,
     .                nk, kpol, wgthpol )


C Direction of the polarization (path for the line integral)  
          npl=nmeshk(igrd,igrd) 
          dkmod=0.0d0
          do ix=1,3
            dk(ix)=rcell(ix,igrd)/npl 
            dkmod=dkmod+dk(ix)**2 
          enddo
          dkmod=dsqrt(dkmod)

C Calculation of the Jacobian
          i=0
          do iy=1,3
            if (iy.ne.igrd) then 
              i=i+1
              do ix=1,3
                vaux(ix,i)=rcell(ix,iy)
              enddo
            endif 
          enddo
          area=dsqrt(
     .      (vaux(2,1)*vaux(3,2)-vaux(3,1)*vaux(2,2))**2 +
     .      (vaux(3,1)*vaux(1,2)-vaux(1,1)*vaux(3,2))**2 +
     .      (vaux(1,1)*vaux(2,2)-vaux(2,1)*vaux(1,2))**2 )
          J=dabs(area/volcel(rcell))
          difA=(3.0d0-dble(nspin))*J

C Construction of the matrix elements of the scalar product dk*r 
          call phirphi(nua, na, nuo, no, scell, xa, rmaxo,
     .              maxnh, lasto, iphorb, isa,
     .              numh, listhptr, listh, dk, aux) 

C Begin the bidimensional integration over the path integrals
          do ik = 1, nk
            do ispin = 1,nspin 
              prodtr = 1.0d0
              prodti = 0.0d0
              phase = 0.0d0
              do il = 0, npl
            
                do ix = 1,3
                  kint(ix) = kpol(ix,ik) + il*dk(ix)
                enddo            

C Find Wavefunctions 
                if (il.ne.npl) then 
                  if (nspin.le.2) then 
                    call diagpol( ispin, nspin, nuo, no,
     .                  nuotot, maxnh, numh, listhptr, 
     .                  listh, H, S, xijo, indxuo, kint,
     .                  ek, psi, 2, Haux, Saux )
                  elseif (nspin.eq.4) then  
               call die('KSV_pol: ERROR: nspin=4 not yet implemented')
                  else
                  call die('KSV_pol: ERROR: incorrect value of nspin')
                  endif
                endif  

C In the first point we just store the wavefunctions      
                if (il.eq.0) then 
                  iuo = 0
                  do io = 1,nuo
                    do jo = 1,nuotot
                      ia = iaorb(jo)
                      dkxij = ddot(3,dk,1,xa(1,ia),1)
                      ckxij = dcos(npl*dkxij)
                      skxij = dsin(npl*dkxij)
                      psi1(iuo+1) = psi(iuo+1)*ckxij + psi(iuo+2)*skxij 
                      psi1(iuo+2) = psi(iuo+2)*ckxij - psi(iuo+1)*skxij
                      iuo = iuo + 2
                    enddo 
                  enddo 
                  detr = 1.0d0
                  deti = 0.0d0 

C Store wavefunction for the next point
                  call savepsi(psiprev,psi,nuo,nuotot,nocc(ispin))

                elseif (il.ne.npl) then 
C Calculate the determinant of the overlap matrix between the 
C periodic Bloch functions in this k point and in the previous one.   
                  call detover(psiprev, psi, S, aux,
     .              numh, listhptr, listh, indxuo, no, nuo, xijo, 
     .              maxnh, nuotot, nocc(ispin), kint, dk, 
     .              detr, deti )
 
C Store wavefunction for the next point
                  call savepsi(psiprev,psi,nuo,nuotot,nocc(ispin))

                else 
C Calculate the determinant of the overlap matrix between the
C periodic Bloch functions in the last k point and the first one
                  call detover(psiprev, psi1, S, aux,
     .              numh, listhptr, listh, indxuo, no, nuo, xijo, 
     .              maxnh, nuotot, nocc(ispin), kint, dk, 
     .              detr, deti )
    
                endif

C Product of all the determinants along the integration path
                prodtr = prodtr*detr - prodti*deti
                prodti = prodti*detr + prodtr*deti

C Phase { i.e. Im[ln(det)] } of the individual determinants 
                if (detr.ne.0.0d0) then
                  if (detr.gt.0d0) then
                    phase = phase + datan(deti/detr) 
                  elseif (deti/detr.gt.0.0d0) then 
                    phase = phase + datan(deti/detr)-pi  
                  else 
                    phase = phase + pi + datan(deti/detr) 
                  endif 
                else
                  if (dabs(deti).gt.tiny) then
                    if (deti.lt.0.0d0) then
                      phase=phase-pi/2.0d0
                    else 
                      phase=phase+pi/2.0d0
                    endif
                  else
                    if (IOnode)
     $               write(6,*) 'KSV_pol: ERROR!!!!, determinant zero'
                  endif
                endif
                phase=phase-nint(phase/(2.0d0*pi))*2.0d0*pi        
              enddo 

C Now phase is defined in the interval [-pi,pi]. It is
C convinient to bring it back to the interval [-2pi:0] to
C avoid abrupt jumps near the value pi or -pi.
       
              if ((shape.ne.'molecule').and.(phase.lt.0.0d0))
     .            phase = 2.0d0*pi + phase

C Continuity of the phase in the plane******************
C
       if(ik.gt.1.and.
     .       dabs(phase-phaseold(ispin)).gt.pi) then
              if((phase-phaseold(ispin)).lt.0.0d0) then
                   phase=phase+2.0d0*pi
              else
                   phase=phase-2.0d0*pi
              endif
       endif
       phaseold(ispin)=phase
C
C*********************************************************

C Calculating the Berry phase from the global phase of the 
C product of all the determinants 
              if (prodtr.ne.0.0d0) then 
                if (prodtr.gt.0.0d0) then 
                  polR(igrd,ispin)=polR(igrd,ispin)+
     .              difA*datan(prodti/prodtr)*wgthpol(ik) 
                elseif(prodti/prodtr.gt.0.0d0) then 
                  polR(igrd,ispin)=polR(igrd,ispin)+
     .              difA*(datan(prodti/prodtr)-pi)*wgthpol(ik)  
                else
                  polR(igrd,ispin)=polR(igrd,ispin)+
     .              difA*(datan(prodti/prodtr)+pi)*wgthpol(ik)
                endif 
              else
                if (dabs(prodti).gt.tiny) then
                  if (prodti.lt.0.0d0) then
                    polR(igrd,ispin)=-difA*pi*wgthpol(ik)/2.0d0
     .                +polR(igrd,ispin)
                  else
                    polR(igrd,ispin)=difA*pi*wgthpol(ik)/2.0d0
     .                +polR(igrd,ispin)
                  endif
                else
                   if (IOnode)
     $              write(6,*) 'KSV_pol: ERROR!!!!, determinant zero'
                endif
              endif

C We can also calculate the Berry phase as the sum of 
C all the individual phases of the determinant. In principle,
C this should be equal to the previous calculation. 
C In practice I think it is more accurate.
              ph(igrd,ispin)=ph(igrd,ispin)+phase*wgthpol(ik)*difA 
       
            enddo
          enddo  
        
        else
          if (igrd.eq.1) then  
            if (IOnode) write(6,'(/,2a)')
     .        'KSV_pol: Polarization not calculated along the',
     .        ' first lattice vector'
          else 
            if (IOnode) write(6,'(/,a,i3,a)')
     .        'KSV_pol: Polarization not calculated along the',
     .        igrd,'th lattice vector'
          endif 
          notcal(igrd)=1
        endif
  
      enddo 
      
      do ispin=1,nspin
 
C Usually adding the individual phases is more accurate than 
C the phase of the product of the determinants, so we drop the
C results obtained from the the product.
        do ix=1,3
          polR(ix,ispin)=ph(ix,ispin)
        enddo 

C Adding the nuclear component of the polarization
        do ix=1,3 
          polR(ix,ispin)=polR(ix,ispin)+polion(ix) 
        enddo 

        do ix=1,3
          if (notcal(ix).eq.1) then 
            polR(ix,ispin)=0.0d0
          endif
        enddo 
       
C Transform to cartesian coordinates
        do ix=1,3
          polxyz(ix,ispin)=uR(ix,1)*polR(1,ispin)+
     .      uR(ix,2)*polR(2,ispin)+uR(ix,3)*polR(3,ispin)
        enddo 
        
c      if(nspin.gt.1) then 
c        if(ispin.eq.1) write(6,'(\a)')
c    . 'KSV_pol: Polarization for spin Up:'
c        if(ispin.eq.2) write(6,'(\a)')
c    . 'KSV_pol: Polarization for spin Down:'
c      endif

c      write(6,'(\a)') 
c    . 'KSV_pol: Macroscopic polarization per unit cell (a.u.):' 
c       write(6,'(\a,\3f12.6)')
c    .   'KSV_pol: Along the lattice vectors  ',
c    .       (polR(ix,ispin),ix=1,3) 
c       write(6,'(\a,\3f12.6)') 
c    .   'KSV_pol: Along cartesian directions ',
c    .    (polxyz(ix,ispin),ix=1,3)
c       write(6,'(\a)')
c    . 'KSV_pol: Macroscopic polarization per unit cell (Debye):'
c       write(6,'(\a,\3f12.6)')
c    .   'KSV_pol: Along the lattice vectors  ',
c    .       (polR(ix,ispin)/Debye,ix=1,3)
c       write(6,'(\a,\3f12.6)')
c    .   'KSV_pol: Along cartesian directions ',
c    .    (polxyz(ix,ispin)/Debye,ix=1,3)

      enddo 

c      if(nspin.gt.1) then 
c        write(6,'(\a,/,a,3f12.6)')
c    .   'KSV_pol: Sum along cartesian directions (a.u.): ',
c    .   'KSV_pol: ',(polxyz(ix,1)+polxyz(ix,2),ix=1,3)
c        write(6,'(\a,/,a,3f12.6)')
c    .   'KSV_pol: Sum along cartesian directions (Debye): ',
c    .   'KSV_pol: ',((polxyz(ix,1)+polxyz(ix,2))/Debye,ix=1,3)
c      endif

C This is the only exit point 
  999 continue

C Deallocate local memory
      call de_alloc( muo,     'muo',     'KSV_pol' )
      call de_alloc( ek,      'ek',      'KSV_pol' )
      call de_alloc( psi1,    'psi1',    'KSV_pol' )
      call de_alloc( psiprev, 'psiprev', 'KSV_pol' )
      call de_alloc( aux,     'aux',     'KSV_pol' )
      call resetDenseMatrix()
      
      if (nkpol.gt.0.and.IOnode) then
        do ispin = 1,nspin
          if (nspin.gt.1) then
            if (ispin.eq.1) write(6,'(/,a)')
     .       'siesta: Macroscopic polarization for spin Up:'
            if (ispin.eq.2) write(6,'(/,a)')
     .       'siesta: Macroscopic polarization for spin Down:'
          endif
          write(6,'(/,a)')
     .     'siesta: Macroscopic polarization per unit cell (a.u.):'
          write(6,'(a,3f12.6)')
     .     'siesta: Along the lattice vectors  ',
     .       (polR(ix,ispin),ix=1,3)
          write(6,'(a,3f12.6)')
     .     'siesta: Along cartesian directions ',
     .      (polxyz(ix,ispin),ix=1,3)
          write(6,'(/,a)')
     .     'siesta: Macroscopic polarization per unit cell (Debye):'
          write(6,'(a,3f12.6)')
     .     'siesta: Along the lattice vectors  ',
     .       (polR(ix,ispin)/Debye,ix=1,3)
          write(6,'(a,3f12.6)')
     .     'siesta: Along cartesian directions ',
     .      (polxyz(ix,ispin)/Debye,ix=1,3) 
        enddo 
        if (nspin.gt.1) then 
C Modified so that compiler is happy when nspin = 1 and bounds checking
C is turned on. JDG
          write(6,'(/,a,/a,3f12.6)')
     .      'siesta: Sum along cartesian directions (a.u.): ',
     .      'siesta: ',(polxyz(ix,1)+polxyz(ix,min(nspin,2)),ix=1,3)
          write(6,'(/,a,/a,3f12.6)')
     .      'siesta: Sum along cartesian directions (Debye): ',
     .      'siesta: ',((polxyz(ix,1)+polxyz(ix,min(nspin,2)))/Debye,
     .      ix=1,3)
        endif
      endif

      call timer( 'KSV_pol', 2 )

      end subroutine KSV_pol
      end module m_ksv
