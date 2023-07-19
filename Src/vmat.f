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
      subroutine vmat( no, np, dvol, nspin, V, nvmax, 
     .                 numVs, listVsptr, listVs, Vs, 
     .                 nuo, nuotot, iaorb, iphorb, isa )

C ********************************************************************
C Finds the matrix elements of the potential.
C First version written by P.Ordejon.
C Name and interface modified by J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of vmat that use a direct algorithm to save memory.
C Modified by J.D.Gale, November'99
C *********************** INPUT **************************************
C integer no              : Number of basis orbitals
C integer np              : Number of columns in C (local)
C real*8  dvol            : Volume per mesh point
C integer nspin           : Number of spin components
C real*4  V(np,nspin)     : Value of the potential at the mesh points
C integer nvmax           : First dimension of listV and Vs, and maxim.
C                           number of nonzero elements in any row of Vs
C integer numVs(nuo)      : Number of non-zero elements in a row of Vs
C integer listVsptr(nuo)  : Pointer to the start of rows in listVs
C integer listVs(nvmax)   : List of non-zero elements of Vs
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(nvmax,nspin) : Value of nonzero elements in each row 
C                           of Vs to which the potential matrix 
C                           elements are summed up
C *********************************************************************

C  Modules
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo
      use listsc_module, only: listsc
      use mesh, only: dxa, nsp, xdop, xdsp
      use meshdscf
      use meshphi
      use parallel,      only: Nodes
      use alloc,         only: re_alloc, de_alloc

      implicit none

C Argument types and dimensions
      integer
     .   no, np, nvmax, nuo, nuotot, iaorb(*), nspin,
     .   iphorb(*), isa(*), numVs(nuo), listVsptr(nuo), listVs(nvmax)
      real(grid_p), intent(in)  ::  V(nsp,np,nspin)
      real(dp)
     .   dvol, Vs(nvmax,nspin)

C Internal variables and arrays

      integer, parameter ::
     .  minloc = 100,  ! Min buffer size for local copy of Dscf
     .  maxoa  = 100   ! Max # of orbitals per atom
      integer
     .  i, ia, ic, ii, ijl, il, imp, ind, iop, ip, iphi, io,
     .  is, isp, ispin, iu, iul, ix, j, jc, jl, 
     .  last, lasta, lastop, maxloc, maxloc2, nc, nlocal, 
     .  nphiloc, nvmaxl

      integer,  pointer :: ilc(:), ilocal(:), iorb(:)

      logical  ParallelLocal
      real(dp)
     .  Vij, dxsp(3), phia(maxoa,nsp), r2cut(nsmax), r2sp

      real(dp), pointer :: VClocal(:)
      real(dp), pointer :: Clocal(:,:), Vlocal(:,:)

C  Start time counter
      call timer('vmat',1)

C  Set algorithm logical
      ParallelLocal = (Nodes.gt.1)

C  Find value of maxloc

      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )

      !  Allocate local memory
      nullify ( ilocal )
      call re_alloc( ilocal, 1, no, name='ilocal', routine='vmat' )

      nullify ( ilc )
      call re_alloc( ilc, 1, maxloc2, name='ilc', routine='vmat' )

      nullify (iorb )
      call re_alloc( iorb, 0, maxloc, name='iorb', routine='vmat' )

      ijl = (maxloc+1)*(maxloc+2)/2

      nullify( Vlocal )
      call re_alloc( Vlocal, 1, ijl, 1, nspin, name='Vlocal',
     &               routine='vmat' )

      nullify( Clocal )
      call re_alloc( Clocal, 1, nsp, 1, maxloc2, name='Clocal',
     &               routine='vmat' )

      nullify( VClocal )
      call re_alloc( VClocal, 1, nsp, name='VClocal', routine='vmat' )

      if (ParallelLocal) then
        if (nrowsDscfL.gt.0) then
          nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        else
          nvmaxl = 1
        endif
        call re_alloc( DscfL, 1, nvmaxl, 1, nspin,
     &                 name='DscfL',  routine='vmat' )
        DscfL(1:nvmaxl,1:nspin) = 0.0_dp
      endif

C  Full initializations done only once
      ilocal(1:no) = 0
      iorb(0:maxloc) = 0
      Vlocal(:,:) = 0.0_dp
      last = 0

C  Find atomic cutoff radii
      r2cut(:) = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Loop over grid points
      do ip = 1,np

C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)

C  Find new required size of Vlocal
        nlocal = last
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
        enddo

C  If overflooded, add Vlocal to Vs and reinitialize it
        if (nlocal .gt. maxloc) then
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            if (ParallelLocal) then
              iul = NeedDscfL(iu)
              if (i .eq. iu) then
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listdl(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin) 
                  enddo
                enddo
              else
                do ii = 1, numdl(iul)
                  ind = listdlptr(iul)+ii
                  j = listsc( i, iu, listdl(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin) 
                  enddo
                enddo
              endif
            else
              if (i .eq. iu) then
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listVs(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin) 
                  enddo
                enddo
              else
                do ii = 1, numVs(iu)
                  ind = listVsptr(iu)+ii
                  j = listsc( i, iu, listVs(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nspin
                    Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .                Vlocal(ijl,ispin) 
                  enddo
                enddo
              endif
            endif
          enddo
          ilocal(iorb(1:last)) = 0
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vlocal(1:ijl,1:nspin) = 0.0_dp
          last = 0
        endif

C  Look for required orbitals not yet in Vlocal
        if (nlocal .gt. last) then
          do ic = 1,nc
            imp = endpht(ip-1) + ic
            i = lstpht(imp)
            if (ilocal(i) .eq. 0) then
              last = last + 1
              ilocal(i) = last
              iorb(last) = i
            endif
          enddo
        endif

C  Loop on first orbital of mesh point
        lasta = 0
        lastop = 0
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)
          ilc(ic) = il

C  Generate or retrieve phi values
          if (DirectPhi) then
            if (ia.ne.lasta .or. iop.ne.lastop) then
              lasta = ia
              lastop = iop
              do isp = 1,nsp
                do ix = 1,3
                  dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
                enddo
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                else
                  phia(:,isp) = 0.0_dp
                endif
              enddo
            endif
            iphi = iphorb(i)
            do isp = 1,nsp
              Clocal(isp,ic) = phia(iphi,isp)
            enddo
          else
            do isp = 1,nsp
              Clocal(isp,ic) = phi(isp,imp)
            enddo
          endif

C  Pre-multiply V and Clocal(,ic)
          do ispin = 1,nspin
            do isp = 1,nsp
              VClocal(isp) = V(isp,ip,ispin) * Clocal(isp,ic)
            enddo

C  Loop on second orbital of mesh point (only for jc.le.ic)
            do jc = 1,ic
              jl = ilc(jc)

C  Loop over sub-points
              Vij = 0.0_dp
              do isp = 1,nsp
                Vij = Vij + VClocal(isp) * Clocal(isp,jc)
              enddo

              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              if (ic.ne.jc.and.il.eq.jl) then
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + 2.0_dp*Vij
              else
                Vlocal(ijl,ispin) = Vlocal(ijl,ispin) + Vij
              endif

            enddo
          enddo
        enddo
      enddo

C  Add final Vlocal to Vs
      do il = 1,last
        i = iorb(il)
        iu = indxuo(i)
        if (ParallelLocal) then
          iul = NeedDscfL(iu)
          if (i .eq. iu) then
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listdl(ind)
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
              enddo
            enddo
          else
            do ii = 1, numdl(iul)
              ind = listdlptr(iul)+ii
              j = listsc( i, iu, listdl(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
              enddo
            enddo
          endif
        else
          if (i .eq. iu) then
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listVs(ind)
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
              enddo
            enddo
          else
            do ii = 1, numVs(iu)
              ind = listVsptr(iu)+ii
              j = listsc( i, iu, listVs(ind) )
              jl = ilocal(j)
              if (il.gt.jl) then
                ijl = il*(il+1)/2 + jl + 1
              else
                ijl = jl*(jl+1)/2 + il + 1
              endif
              do ispin = 1,nspin
                Vs(ind,ispin) = Vs(ind,ispin) + dVol * 
     .            Vlocal(ijl,ispin) 
              enddo
            enddo
          endif
        endif
      enddo

  999 continue

      if (ParallelLocal) then
C Redistribute Hamiltonian from mesh to orbital based distribution
        call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, 
     .      nuotot, nspin, DscfL, Vs )
C Free memory 
        call de_alloc( DscfL,  name='DscfL' )
      endif

      !  Free local memory

      call de_alloc( VClocal,  name='VClocal' )
      call de_alloc( Clocal,  name='Clocal' )
      call de_alloc( Vlocal,  name='Vlocal' )
      call de_alloc( iorb,  name='iorb' )
      call de_alloc( ilc,  name='ilc' )
      call de_alloc( ilocal,  name='ilocal' )

      call timer('vmat',2)
      return
      end
