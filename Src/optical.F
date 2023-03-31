! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine optical( nua, na, xa, scell, ucell, nuotot, nuo, no, 
     &                    nspin, qspin, maxnh, numh, 
     &                    listhptr, listh, H, S, 
     &                    xijo, indxuo, ek, efermi, temp,
     &                    isa, iphorb, iphkb, lasto, lastkb, shape )
C *********************************************************************
C Calculates the imaginary part of the dielectric function using
C first order perturbation theory.
C Written by DSP, August 1999.
C Restyled for f90 version by JDG. June 2004
C Small modification by DSP, May 2010
C **************************** INPUT **********************************
C integer nua                 : Number of atoms in the unit cell
C integer na                  : Number of atoms 
C real*8  xa(3,na)            : Atomic positions in cartesian coordinates
C real*8  scell(3,3)          : Supercell vectors scell(ixyz,ivect)
C real*8  ucell(3,3)          : Unit cell vectors
C integer nuotot              : No. of basis orbitals in the unit cell (global)
C integer nuo                 : No. of basis orbitals in the unit cell (local)
C integer no                  : Number of basis orbitals
C integer nspin               : Number of spin components
C real*8  qspin(2)            : Total population of spin up and down
C integer maxnh               : Maximum number of orbitals interacting  
C                               with any orbital
C integer numh(nuo)           : Number of nonzero elements of each row 
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element  
C                               column indexes for each matrix row
C integer listhptr(nuo)       : Pointer to start of row in listh
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  xijo(3,*)           : Vectors between orbital centers (sparse)
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C                               Unit cell orbitals must be the first in
C                               orbital lists, i.e. indxuo.le.nuo, with
C                               nuo the number of orbitals in unit cell
C real*8  ek(nuo)             : Array to store the eigenvalues
C real*8  efermi              : Fermi energy
c real*8  temp                : Electronic temperature 
C integer isa(na)             : Species index for each atom
C integer iphorb(no)          : Orbital index within the atom for each
C                               orbital
C integer iphKB(nokb)         : KB projector index within the atom for 
C                               each KB projector
C integer lasto(0:na)         : Last orbital index of each atom 
C integer lastkb(0:na)        : Last KB projector of each atom
C *************************** INTERNAL ********************************
C integer  nkopt              : Maximum number of grid points in the
C                               BZ integration
C real*8   kopt(3,nkopt)      : Auxiliar array to store the kpoints
C                               for the BZ integration
C real*8   wgthopt(nkopt)     : Auxiliar array to store the weigths of
C                               the kpoints.            
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C Energies in Rydbergs.
C *********************************************************************

      use precision,    only : dp
      use atomlist,     only : qtot
      use densematrix,  only : allocDenseMatrix, resetDenseMatrix
      use densematrix,  only : Haux, Saux, psi
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb
      use alloc
      use sys,          only : die
      use m_fermid,     only : stepf
      use files,        only : slabel, label_length
      use units,        only : pi
#ifdef MPI
      use mpi_siesta
#endif
      implicit          none
C     Passed variables
      integer           maxnh, nuotot, nuo, no, nspin
      integer           indxuo(no), listh(maxnh), numh(nuo)
      integer           listhptr(nuo)
      integer           na, isa(na), iphorb(no)
      integer           nua, lasto(0:na),lastkb(0:na)
      integer           iphKB(*)
      real(dp)          ek(nuotot), qspin(2), H(maxnh,nspin), 
     &                  S(maxnh), xijo(3,*),
     &                  ucell(3,3), efermi, temp,
     &                  scell(3,3), xa(3,na)

C     Internal variables 
      real(dp) ntote_real      !!** Note: needed to deal with synthetics
      integer
     &  iie, ik, io, ispin, iuo, ix, iy,  is, ia, ntote,
     &  nbands, ixmin, ije, j, maxhs, maxpsi,
     &  kscell(3,3), nk, nmeshk(3), npol, ipol, 
     &  ie, je, iw, nwi, nwf, nw, iu, ngamma, ig

#ifdef MPI
      integer                     :: MPIerror
#endif
      integer                     :: nkopt
         
      logical                     :: DoOptical
      logical                     :: gamma

      real(dp)                    :: efermiloc
      real(dp),        parameter  :: tiny = 1.0d-9
      real(dp),        parameter  :: ediffmin = 1.0d-3
      real(dp), pointer           :: epsimg(:,:) => null()
      real(dp), pointer           :: ekloc(:) => null()
      real(dp), pointer           :: kopt(:,:) => null()
      real(dp), pointer           :: wgthopt(:) => null()
#ifdef MPI
      real(dp), pointer           :: epsbuff(:,:) => null()
#endif
      real(dp), pointer           :: intraband(:,:) => null()
      real(dp), pointer           :: aux(:,:) => null()
      real(dp)
     &  cutoff, displ(3), dmod, dsp, dw, ediff, efield(3),
     &  fsum(2), f1, f2, inplane(3,2), kint(3),
     &  smooth, T, volcel, volume, xv(3), w,
     &  wmin, wmax, wi, wf, escissor, omega_p2(2),aux_p2(2)

      character(len=label_length+7) :: filename
      character(len=10)             :: shape
      character(len=12)             :: option
      character(len=1)              :: matrix

      external          io_assign, io_close, volcel

#ifdef DEBUG
      call write_debug( '  PRE optical' )
#endif
C     Start time counter 
      call timer( 'optical', 1 )

C     Initialise nbands to maximum number possible 
      nbands = nuotot

C     Read optical adsorption parameters
      call reoptical(DoOptical,nmeshk,dsp,wmin,wmax,smooth,
     &               efield,nbands,option,escissor)

      if (.not.DoOptical) goto 999

C     Unit cell volume
      volume = volcel(ucell)

C     Set initial value of nkopt and gamma
      nkopt = 0
      ngamma = 2
      gamma = .false.

C     Find the integration grid in reciprocal space   
      nk = nmeshk(1)*nmeshk(2)*nmeshk(3)

C     Use only gamma for atoms and molecules
C     unless otherwise requested
       if (((shape.eq.'atom').or.(shape.eq.'molecule'))
     &     .and.(nk.eq.0)) then 
        do ix=1,3
          nmeshk(ix)=1
        enddo 
        nk=1
	dsp=0.0d0
       endif

C     If nk.ne.0
      if (nk.ne.0) then
        do ix = 1,3
          do iy = 1,3
            kscell(ix,iy) = 0
          enddo
        enddo
        do ix = 1,3
          kscell(ix,ix) = nmeshk(ix)
          displ(ix) = dsp
        enddo
  
        cutoff = 0.0d0
C    DSP, May 2010 call to kgridinit has been eliminated because
C        the results are not coherent with those of kgrid
c       call kgridinit( ucell, kscell,  displ, cutoff, nkopt )
        nkopt=nk
C       Resize arrays to fit nkopt
        call re_alloc( kopt, 1, 3, 1, nkopt, 'kopt','optical' )
        call re_alloc( wgthopt, 1, nkopt, 'wgthopt', 'optical')

C       Obtain the points and weigths for the kspace integration
        do ix = 1,3
          do iy = 1,3
            kscell(ix,iy) = 0
          enddo
        enddo
        do ix = 1,3
          kscell(ix,ix) = nmeshk(ix)
          displ(ix) = dsp
        enddo
        call kgrid( ucell, kscell, displ, nkopt, kopt, wgthopt )
        ngamma=2
        if (nkopt.eq.1 .and. abs(kopt(1,1)).lt.tiny .and.
     .                  abs(kopt(2,1)).lt.tiny .and.
     .                  abs(kopt(3,1)).lt.tiny) then
          gamma = .true.
          ngamma = 1
        endif
        
C If nk.eq.0
      else
        nkopt = 0
      endif

C     Check whether nkopt is greater than 0
      if (nkopt.eq.0) then
         if (Node .eq. 0) then
          write(6,'(/a)') 
     &     'optical: no grid specified for optical calculation' 
          write(6,'(/a)')
     &     'optical: no calculation performed'
         endif
         goto 999
      endif 


      if (nspin.eq.1) then 
C       Total number of valence electrons in the unit cell
! (Taking into account net charge...)
        ntote_real = qtot
        if ((nint(ntote_real) - ntote_real) .gt. 1.0e-6_dp) then
           call die('Optical: Non-integer number of electrons')
        endif
        ntote = nint(ntote_real)
        if (nbands.le.nint(ntote/2.0d0)) then 
          if (Node.eq.0) then
            write(6,'(/a,i5)') 
     &        'Optical: You need include at least', nint(ntote/2.0d0)
            write(6,'(a)') 
     &        'Optical: bands for the occupied states,' 
            write(6,'(a)')
     &        'Optical: and some excited states'       
            write(6,'(a)')
     &        'Optical: No calculation performed'
          endif
          goto 999
        endif 
      else
        if (nbands.le.nint(qspin(1))) then
          if (Node.eq.0) then
            write(6,'(/a,i5)')
     &        'Optical: You need include at least', nint(qspin(1))
            write(6,'(a)')
     &        'Optical: bands for the spin Up occupied states,' 
            write(6,'(a)')
     &        'Optical: and some excited states'
            write(6,'(a)')
     &        'Optical: No calculation performed'
          endif
          goto 999
        endif 
        if (nbands.le.nint(qspin(2))) then
          if (Node.eq.0) then
            write(6,'(/a,i5)')
     &       'Optical: You need include at least', nint(qspin(2))
            write(6,'(a)')
     &       'Optical: bands for the spin Up occupied states,'
            write(6,'(a)')
     &       'Optical: and some excited states'
            write(6,'(a)')
     &       'Optical: No calculation performed'
          endif
          goto 999
        endif 
      endif

C     Find the size of the array to hold the imaginary dielectric constant
      dw = smooth/10
      nw = nint((wmax-wmin)/dw) + 1
      call re_alloc( epsimg, 1, nw, 1, 2, 'epsimg', 'optical')

C     Check size of remaining arrays
      call re_alloc( ekloc, 1, nuotot, 'ekloc', 'optical')
      
C     The calculation  of the imaginary part of the dielectric 
C     function starts here


C     Depending on the selected option, the meaning of the vector
C     efield changes 
      if (Node.eq.0) then
        write(6,'(/a,a)')   'Optical: Performing optical calculation: '
        write(6,'(/a,a)')     'Optical: Polarization type       = ',
     &    option
        write(6,'(a,f8.4,a)') 'Optical: Minimum of energy range = ',
     &    wmin,' Ry '
        write(6,'(a,f8.4,a)') 'Optical: Maximum of energy range = ',
     &    wmax,' Ry '
        write(6,'(a,f8.4,a)') 'Optical: Gaussian broadening     = ',
     &    smooth,' Ry '
        write(6,'(a,f8.4,a)') 'Optical: Scissor operator        = ',
     &    escissor,' Ry '
        write(6,'(a,i8)')     'Optical: Number of bands         = ',
     &    nbands
        write(6,'(a,i8)')     'Optical: Number of electrons     = ',
     .    ntote
        write(6,'(a,3i4)')    'Optical: BZ mesh dimensions      = ',
     &    (nmeshk(j),j=1,3)
      endif

      if ( option.eq.'polarized') then
C       For a polarized field, the vector efield gives the direction
C       of the electric field 
        dmod = 0.0d0
        do ix = 1,3
          dmod = dmod + efield(ix)**2
        enddo 
        if (dmod.gt.tiny) then 
          do ix = 1,3 
            efield(ix) = efield(ix)/dsqrt(dmod)
          enddo
        else
          call die('Optical: ERROR: zero electric field')
        endif

        if (Node.eq.0) then
          write(6,'(/a)') 
     &       'Optical: electric field direction'
          write(6,'(a,3f12.5)') 'Optical: ',
     &              (efield(ix),ix=1,3)
        endif

      elseif (option.eq.'unpolarized') then 
C       For the unpolarized option we have to calculate two vectors
C       perpendicular to the direcction of the vector efield ( which
C       now stands for the propagation direction.
        dmod = 0.0d0
        do ix = 1,3
          dmod = dmod + efield(ix)**2
        enddo 
        if (dmod.gt.tiny) then 
          do ix = 1,3 
            efield(ix) = efield(ix)/dsqrt(dmod)
          enddo
        else
          call die('Optical: ERROR: no propagation direction given')
        endif

        do ix = 1,3
          inplane(ix,1) = 0.0d0
          inplane(ix,2) = 0.0d0
        enddo 
        do ix = 1,3
          if (dabs(efield(ix)).lt.tiny) then 
            inplane(ix,1) = 1.0d0
            goto 10
          endif
        enddo
        dmod = 1.0d0
        do ix = 1,3
          if (dabs(efield(ix)).lt.dmod) then 
            dmod = dabs(efield(ix))
            ixmin = ix  
          endif
        enddo 
        inplane(ixmin,1) = 1.0d0
        do ix = 1,3
          inplane(ix,1) = inplane(ix,1) - dmod*efield(ix)
        enddo 
        dmod = 0.0d0
        do ix = 1,3
          dmod = dmod + inplane(ix,1)**2
        enddo
        do ix = 1,3
          inplane(ix,1) = inplane(ix,1)/dsqrt(dmod)
        enddo
10      continue    
        inplane(1,2) = inplane(2,1)*efield(3) - inplane(3,1)*efield(2)
        inplane(2,2) = inplane(3,1)*efield(1) - inplane(1,1)*efield(3)
        inplane(3,2) = inplane(1,1)*efield(2) - inplane(2,1)*efield(1)
        dmod = 0.0d0
        do ix = 1,3
          dmod = dmod + inplane(ix,2)**2
        enddo
        do ix = 1,3
          inplane(ix,2) = inplane(ix,2)/dsqrt(dmod)
        enddo
        if (Node.eq.0) then
          write(6,'(/a)')
     &'Optical: Vectors defining the plane that contains the Efield:'
             write(6,'(a,3f12.5)')
     &           'Optical: ', (inplane(ix,1),ix=1,3)
             write(6,'(a,3f12.5)')
     &           'Optical: ', (inplane(ix,2),ix=1,3)
        endif 
      endif 

C     Allocate Haux/Saux/psi arrays to required size
      if (nspin.le.2 .and. gamma) then
        maxhs = nuotot*nuo
        maxpsi = nuotot*nuo
      else
        maxhs = 2*nuotot*nuo
        maxpsi = 2*nuotot*nuo
      endif

      call allocDenseMatrix(maxhs, maxhs, maxpsi)

      ! Allocate auxiliary array for calculating phirphi
      call re_alloc( aux,  1, maxnh, 1, 3, 'aux', 'optical')

C     What kind of matrix elements will be used.
C     Real space is possible for molecules and atoms
C     For infinte media we will need caculate everything
C     using the momentum
      matrix = 'P' 
      if (((shape.eq.'atom').or.(shape.eq.'molecule'))
     &    .and.gamma) matrix = 'R' 
C     parameter npol
      if (option.eq.'polarized') then
        npol = 1
      elseif (option.eq.'unpolarized') then
        npol = 2
      else
        npol = 3
      endif
C     Allocate intraband array to calculate the Drude contribution
C     for metals
      call re_alloc( intraband, 1, ngamma, 1, nuo,
     &               'intraband', 'optical' )

C     Construction of the matrix elements of the scalar product n*r 
        
      if ( option.eq.'polarized') then 
        call phirphi_opt( nua, na, nuo, no, scell, xa, maxnh, 
     &                    lasto, lastkb, iphorb, iphKB, isa, numh,
     &                    listhptr, listh, efield, matrix, aux(1,1) )
        do io = 1,maxnh
          aux(io,2) = 0.0d0
          aux(io,3) = 0.0d0
        enddo
      elseif( option.eq.'unpolarized') then
        call phirphi_opt( nua, na, nuo, no, scell, xa, maxnh, lasto,
     &                    lastkb, iphorb, iphKB, isa, numh, listhptr,
     &                    listh, inplane(1,1), matrix, aux(1,1) )
        call phirphi_opt( nua, na, nuo, no, scell, xa, maxnh, lasto,
     &                    lastkb, iphorb, iphKB, isa, numh, listhptr,
     &                    listh, inplane(1,2), matrix, aux(1,2) )
        do io = 1,maxnh
          aux(io,3) = 0.0d0
        enddo
      else
         do ix = 1 , 3
            xv = 0._dp
            xv(ix) = 1._dp
            call phirphi_opt( nua, na, nuo, no, scell, xa, maxnh, 
     &           lasto, lastkb, iphorb, iphKB, isa, numh, 
     &           listhptr, listh, xv, matrix, aux(1,ix) )
         end do
      endif

C     Initialize array for the dielectric function
      do ix = 1,nw
        epsimg(ix,1) = 0.0_dp
        epsimg(ix,2) = 0.0_dp
      enddo
C...  For Drude formula.....
      omega_p2(1) = 0.0_dp  
      omega_p2(2) = 0.0_dp

C     Begin the BZ integration
      do ik = 1,nkopt
        do ispin = 1,nspin 
          do ix = 1,3
            kint(ix) = kopt(ix,ik)
          enddo            
             
C         Find Wavefunctions 
        
          if (nspin.le.2) then 
            call diagpol( ispin, nspin, nuo, no, nuotot,
     &                    maxnh, numh, listhptr, listh, H, S,
     &                    xijo, indxuo, kint, ek, psi, 
     &                    ngamma, Haux, Saux )
          elseif (nspin.eq.4) then  
            call die('Optical: ERROR: nspin=4 not yet implemented')
          else
            call die('Optical: ERROR: incorrect value of nspin')
          endif

C         Copy eigenvalues to an auxiliary array
          efermiloc = efermi
          do iuo = 1,nuotot
            ekloc(iuo)=ek(iuo)
          enddo

C         Apply scissor operator if required
          if (escissor.gt.0.0d0) then
            efermiloc = efermiloc + 0.5d0*escissor
            do iuo = 1,nuotot
              if (ekloc(iuo).gt.efermi) then
                ekloc(iuo) = ekloc(iuo) + escissor
              endif
            enddo
          endif

          do ipol = 1,npol
C Matrix elements for the dipolar transition are stored in Haux
C DSP May 2010, modified call to subroutine transition_rate 
C information wmin and wmax
C passed to subroutine transition_rate so only relevant
C transition probabilities are explicitly calculated.
             call transition_rate( ngamma, psi, ek, efermi,
     .            ekloc, efermiloc, temp,
     .            smooth, wmin, wmax,
     .            aux(1,ipol), Haux, Saux, numh, 
     .            listhptr, listh, indxuo, no, 
     .            nuo, nuotot, xijo, maxnh, nbands, 
     .            kint, matrix,intraband)

            do ie = 1,nbands
              call GlobalToLocalOrb(ie,Node,Nodes,iie)
              if (iie.gt.0) then
C               Drude term for metals
                T=0.0_dp
                do ig=1,ngamma
                  T = T+ intraband(ig,iie)**2
                enddo
                T = T/dble(npol)
                T = 4.0d0*pi*T/volume
                T = 2.0d0*T*wgthopt(ik)/dble(nspin)
                ediff=ekloc(ie)-efermiloc
C               factor of two because smooth is in Ry rather than in Ha
                omega_p2(ispin)= omega_p2(ispin) +
     &          2.0d0*dexp(-(ediff/smooth)**2)*T/(dsqrt(pi)*smooth)

                f1 = stepf((ekloc(ie)-efermiloc)/temp)
                do je = 1,nbands 
                  if (dabs(ekloc(ie)-ekloc(je)).gt.ediffmin) then 
                    f2 = stepf((ekloc(je)-efermiloc)/temp)
                    if (f1*(1.0d0-f2).gt.tiny) then 
                      ediff = ekloc(je) - ekloc(ie)
                      if ((ediff.ge.max(wmin-2.0d0*smooth,0.0d0)).and.
     &                        (ediff.le.wmax+2.0d0*smooth)) then
                        ije = ngamma*(je-1)+ngamma*nuotot*(iie-1)
                        if(ngamma.eq.2) then 
                          T = Haux(1+ije)**2+Haux(2+ije)**2 
                        else
                          T = Haux(1+ije)**2
                        endif
                        T = T*f1*(1.0d0-f2)/dble(npol)
                        T = 4.0d0*(pi**2)*T/volume
                        T = 2.0d0*T*wgthopt(ik)/dble(nspin)
C Modified April 2010
                        wi = max(wmin,ediff-20.0d0*smooth) 
                        wf = min(wmax,ediff+20.0d0*smooth) 
                        nwi = int((wi-wmin)/dw)+1
                        nwf = int((wf-wmin)/dw)+1 
                        do iw = nwi,nwf
                          w = wmin + (iw-1)*dw
C                         The factor of two comes from the fact that
C                         smooth is in Ry instead of Ha.
                          epsimg(iw,ispin) = epsimg(iw,ispin) +
     &                        2.0d0*dexp(-((ediff-w)/smooth)**2)*T/
     &                        (dsqrt(pi)*smooth)
                        enddo
                      endif
                    endif
                  endif
                enddo
              endif
            enddo
C           Closing loop on polarization options
          enddo 
C         Closing loop on spin orientations
        enddo 
C       Closing loop on k_points
      enddo 

#ifdef MPI
      call re_alloc( epsbuff, 1, nw, 1, 2, 'epsbuff', 'optical' )
      call MPI_Reduce(epsimg(1,1),epsbuff(1,1),2*nw,
     &  MPI_double_precision,MPI_sum,0,MPI_Comm_World,MPIerror)
      epsimg(1:nw,1:2) = epsbuff(1:nw,1:2)
      call de_alloc( epsbuff, 'epsbuff', 'optical' )

      call MPI_Reduce(omega_p2(1),aux_p2(1),nspin,
     &  MPI_double_precision,MPI_sum,0,MPI_Comm_World,MPIerror)
      omega_p2(1)=aux_p2(1)
      if(nspin.gt.1) omega_p2(2)=aux_p2(2)
#endif

C Output results
      if (Node.eq.0) then
        call io_assign(iu)
        filename = trim(slabel) // '.EPSIMG'
        open(iu,file=filename, status='unknown')
        nwf = int((wmax-wmin)/dw) + 1
        do ispin = 1,nspin
          dmod = 0.0d0
          do iw = 1,nwf 
            w = wmin + (iw-1)*dw
            dmod = dmod + epsimg(iw,ispin)*w*dw*volume
          enddo 
C A factor of 4 due to the integration in Ry rather than in a.u.
          fsum(ispin) = dmod/4.0d0
        enddo 
      
        write(iu,'(a)') '## Minimum and maximum energy in eV'
        write(iu,'(a,2f10.4)') 
     &    '## ', wmin*13.6058d0, (wmin+(nwf-1)*dw)*13.6058d0
        write(iu,'(a)') '## Number of spin components'
        write(iu,'(a,i3)') '## ',nspin
        if (nspin.gt.1) then 
          write(iu,'(a)') '## f-sum rule spin 1' 
          write(iu,'(a,f10.4)') '## ',fsum(1)/2.0d0/pi**2/qspin(1) 
          write(iu,'(a)') '## f-sum rule spin 2'
          write(iu,'(a,f10.4)') '## ',fsum(2)/2.0d0/pi**2/qspin(2) 
          write(iu,'(a)') '## wp^2 in Ha^2 for Drude, spin 1'
          write(iu,'(a,f10.4)') '## ',omega_p2(1)
          write(iu,'(a)') '## wp^2 in Ha^2 for Drude, spin 2'
          write(iu,'(a,f10.4)') '## ',omega_p2(2)
        else
          write(iu,'(a)') '## f-sum rule'
          write(iu,'(a,f10.4)') '## ',fsum(1)/2.0d0/pi**2/dble(ntote)
          write(iu,'(a)') '## wp^2 in Ha^2 for Drude'
          write(iu,'(a,f10.4)') '## ',omega_p2(1)
        endif   
        do iw = 1,nwf
          w = wmin + (iw-1)*dw
          write(iu,*) w*13.6058d0, (epsimg(iw,ispin),ispin=1,nspin)
        enddo 
        call io_close(iu)
      endif
  
      if (nspin.eq.2) ntote = nint(qspin(1)+qspin(2))
      if (nspin.eq.2) fsum(1) = fsum(1) + fsum(2)
      if (Node.eq.0) then
        write(6,'((/a,f10.6),3(/,2a))') 
     &   'Optical: Checking f-sum rule',
     &       fsum(1)/2.0d0/pi**2/dble(ntote),
     &   'Optical: ',
     &   'For insulators and closed shell systems', 
     &   'Optical: this number should be close',
     &   ' to 1 if an enough',
     &   'Optical: number of unoccupied bands have been',
     &   ' included' 
      endif
       
C     This is the only exit point 
  999 continue

C     Deallocate internal arrays
      call de_alloc( epsimg,    'epsimg',    'optical' )
      call de_alloc( intraband, 'intraband', 'optical' )
      call de_alloc( ekloc,     'ekloc',     'optical' )
      call de_alloc( kopt,      'kopt',      'optical' )
      call de_alloc( wgthopt,   'wgthopt',   'optical' )
      call de_alloc( aux,       'aux',       'optical' )
      call resetDenseMatrix()
      
      call timer( 'optical', 2 )

#ifdef DEBUG
      call write_debug( '  POS optical' )
#endif
      end
