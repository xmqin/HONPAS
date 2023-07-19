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
      subroutine rhoofdsp( no, np, maxnd, numd, listdptr, 
     .                     listd, nspin, Dscf, rhoscf, nuo, nuotot, 
     .                     iaorb, iphorb, isa, q )
C ********************************************************************
C Finds the SCF density at the mesh points from the density matrix.
C Written by P.Ordejon and J.M.Soler. May'95.
C Re-ordered so that mesh is the outer loop and the orbitals are
C handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
C Version of rhoofd that optionally uses a direct algorithm to save 
C memory. Modified by J.D.Gale, November'99
C Spiral version written by V. M. Garcia-Suarez. June 2002.
C *********************** InpUT **************************************
C integer no              : Number of basis orbitals
C integer np              : Number of mesh points
C integer maxnd           : First dimension of listD and Dscf, and
C                           maximum number of nonzero elements in
C                           any row of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of rows in listd
C integer listd(maxnd)    : List of nonzero elements in each row of Dscf
C integer nspin           : Number of spin components
C real*8  Dscf(maxnd)     : Rows of Dscf that are non-zero 
C integer nuo             : Number of orbitals in unit cell locally
C integer nuotot          : Number of orbitals in unit cell in total
C integer iaorb(*)        : Pointer to atom to which orbital belongs
C integer iphorb(*)       : Orbital index within each atom
C integer isa(*)          : Species index of all atoms
C real*8 q(3)             : Wave vector for spiral configuration.
C *********************** OUTPUT **************************************
C real    rhoscf(nsp,np)  : SCF density at mesh points
C *********************************************************************

C  Modules
      use precision, only: dp, grid_p
      use atmfuncs,  only: rcut, all_phi
      use atm_types, only: nsmax=>nspecies
      use atomlist,  only: indxuo
      use listsc_module, only: listsc
      use mesh, only: nsp, dxa, xdop, xdsp, cmesh, nmeshg, nsm
      use meshdscf
      use meshphi
      use sys, only : die
      use alloc,     only: re_alloc, de_alloc
      implicit none

C Argument types and dimensions
      integer
     .   no, np, nspin, maxnd, nuo, nuotot, iaorb(*),
     .   iphorb(*), isa(*), numd(nuo), listdptr(nuo), listd(maxnd)

      real(grid_p), intent(out) ::  rhoscf(nsp,np,nspin)

      real(dp)
     .   Dscf(maxnd,nspin), q(3)

      external
     .   memory, timer, ipack

C Internal variables and arrays
      integer, parameter ::
     .  minloc = 100,  ! Min buffer size for local copy of Dscf
     .  maxoa  = 100   ! Max # of orbitals per atom

      integer
     .  i, ia, ic, ii, ijl, il, imp, ind, ispin, io, ix,
     .  iop, ip, iphi, is, isp, iu, iul, j, jc, jl, nsd, 
     .  last, lasta, lastop, maxloc, maxloc2, maxndl, nc, nphiloc,
     .  iii(3)

      integer, dimension(:), pointer, save :: 
     .  ilc, ilocal, iorb

      logical
     .  Parallel_Flag

      real(dp)
     .  Cij(nsp), Dij, Dij1, Dij2, Dij3, Dij4,
     .  dxsp(3), phia(maxoa,nsp), r2cut(nsmax), r2sp,
     .  xr(3), Rdi(3), qRdi, cqRdi, sqRdi

      real(dp), dimension(:,:), pointer ::
     .  Clocal, Dlocal

      real(dp), dimension(:,:,:), pointer, save ::
     .  DlocalSp

C  Start time counter
      call timer('rhoofdsp',1)

C  Set algorithm logical
      Parallel_Flag = (nuo .ne. nuotot)

C  Find size of buffers to store partial copies of Dscf and C
      maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
      maxloc = maxloc2 + minloc
      maxloc = min( maxloc, no )

C  If spiral, the diagonal elements of Vlocal do not change
      nsd = 2

C  Allocate local memory
      nullify ( ilocal )
      call re_alloc( ilocal, 1, no, name='ilocal', routine='rhoofdsp' )

      nullify ( ilc )
      call re_alloc( ilc, 1, maxloc2, name='ilc', routine='rhoofdsp' )

      nullify (iorb )
      call re_alloc( iorb, 0, maxloc, name='iorb', routine='rhoofdsp' )

      ijl = (maxloc+1)*(maxloc+2)/2

      nullify( Dlocal )
      call re_alloc( Dlocal, 1, ijl, 1, nsd, name='Dlocal',
     &               routine='rhoofdsp' )

      nullify( DlocalSp )
      call re_alloc( DlocalSp, 0, maxloc, 0, maxloc, 1, nsd, 
     &               name    = 'DlocalSp',
     &               routine = 'rhoofdsp' )

      nullify( Clocal )
      call re_alloc( Clocal, 1, nsp, 1, maxloc2, name='Clocal',
     &               routine='rhoofdsp' )

      if (Parallel_Flag) then
        if (nrowsDscfL.gt.0) then
          maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        else
          maxndl = 1
        endif
        call re_alloc( DscfL, 1, maxndl, 1, nsd+2, name='DscfL',
     &               routine='rhoofdsp' )
C Redistribute Dscf to DscfL form
        call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, nuotot,
     .    nsd, Dscf, DscfL )
      endif

C  Find atomic cutoff radiae
      r2cut(:) = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Initializations
      rhoscf(:,:,:) = 0.0_grid_p
      Dlocal(:,:) = 0.0_dp
      DlocalSp(:,:,:) = 0.0_dp
      ilocal(:) = 0
      iorb(:) = 0
      last = 0

C  Loop over grid points
      do ip = 1,np

C  Find point coordinates
        call ipack(-1,3,nmeshg/nsm,iii,ip)
        do ix = 1,3
          xr(ix) = iii(1) * cmesh(ix,1) + iii(2)*cmesh(ix,2) +
     .             iii(3) * cmesh(ix,3)
        enddo
 
C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)

C  iorb(il)>0 means that row il of Dlocal must not be overwritten
C  iorb(il)=0 means that row il of Dlocal is empty
C  iorb(il)<0 means that row il of Dlocal contains a valid row of 
C             Dscf, but which is not required at this point
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          if (il.gt.0) iorb(il) = i
        enddo

C  Look for required rows of Dscf not yet stored in Dlocal
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) then

C           Look for an available row in Dlocal
            do il = 1,maxloc
C             last runs circularly over rows of Dlocal
              last = last + 1
              if (last .gt. maxloc) last = 1
              if (iorb(last) .le. 0) goto 10
            enddo
            call die('rhoofdsp: no slot available in Dlocal')
   10       continue

C  Copy row i of Dscf into row last of Dlocal
            j = abs(iorb(last))
            if (j.ne.0) ilocal(j) = 0
            ilocal(i) = last
            iorb(last) = i
            il = last
            iu = indxuo(i)
            if (Parallel_Flag) then
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
                  do ispin = 1,nsd
                    Dij = DscfL(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                    Dij = DscfL(ind,ispin+2)
                    DlocalSp(il,jl,ispin) = Dij
                    DlocalSp(jl,il,ispin) = Dij
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
                  do ispin = 1,nsd
                    Dij = DscfL(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                    Dij = DscfL(ind,ispin+2)
                    DlocalSp(il,jl,ispin) = Dij
                    DlocalSp(jl,il,ispin) = Dij
                  enddo
                enddo
              endif
            else
              if (i .eq. iu) then
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listd(ind)
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    Dij = Dscf(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                    Dij = Dscf(ind,ispin+2)
                    DlocalSp(il,jl,ispin) = Dij
                    DlocalSp(jl,il,ispin) = Dij
                  enddo
                enddo
              else
                do ii = 1, numd(iu)
                  ind = listdptr(iu)+ii
                  j = listsc( i, iu, listd(ind) )
                  jl = ilocal(j)
                  if (il.gt.jl) then
                    ijl = il*(il+1)/2 + jl + 1
                  else
                    ijl = jl*(jl+1)/2 + il + 1
                  endif
                  do ispin = 1,nsd
                    Dij = Dscf(ind,ispin)
                    Dlocal(ijl,ispin) = Dij
                    Dij = Dscf(ind,ispin+2)
                    DlocalSp(il,jl,ispin) = Dij
                    DlocalSp(jl,il,ispin) = Dij
                  enddo
                enddo
              endif
            endif
          endif
        enddo

        lasta=0
        lastop=0

C  Generate or retrieve phi values for all orbitals up to nc
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)
          ilc(ic) = il

          if (DirectPhi) then
            if (ia.ne.lasta .or. iop.ne.lastop) then
              lasta = ia
              lastop = iop
              do isp = 1,nsp
                dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                r2sp = sum(dxsp**2)
                if (r2sp.lt.r2cut(is)) then
                  call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                else
                  phia(:,isp) = 0.0_dp
                endif
              enddo
            endif
            iphi = iphorb(i)
            Clocal(:,ic) = phia(iphi,:)
          else
            Clocal(:,ic) = phi(:,imp)
          endif

        enddo

C  Loop on first orbital of mesh point
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          iu = indxuo(i)
          ia = iaorb(i)
          is = isa(ia)
          iop = listp2(imp)

C  Calculate spiral phase
          Rdi(1:3) = xr(1:3) - xdop(1:3,iop) + dxa(1:3,ia)
          qRdi = q(1) * Rdi(1) + q(2) * Rdi(2) + q(3) * Rdi(3)
          cqRdi = cos(qRdi)
          sqRdi = sin(qRdi)

C  Loop on second orbital of mesh point (NOT only for jc.le.ic)
          do jc = 1,nc
            jl = ilc(jc)

            do isp = 1, nsp
              Cij(isp) = Clocal(isp,ic) * Clocal(isp,jc)
            enddo

            if (il.gt.jl) then
              ijl = il*(il+1)/2 + jl + 1
            else
              ijl = jl*(jl+1)/2 + il + 1
            endif

            if (jc.le.ic) then
              if (ic .eq. jc) then
                Dij1 = Dlocal(ijl,1)
                Dij2 = Dlocal(ijl,2)
              else
                Dij1 = 2*Dlocal(ijl,1)
                Dij2 = 2*Dlocal(ijl,2)
              endif
            endif
            Dij3 = DlocalSp(il,jl,1)
            Dij4 = DlocalSp(il,jl,2)

C  Loop over sub-points
            do isp = 1, nsp
              if (jc.le.ic) then
                rhoscf(isp,ip,1) = rhoscf(isp,ip,1) +
     .            Dij1*Cij(isp)
                rhoscf(isp,ip,2) = rhoscf(isp,ip,2) +
     .            Dij2*Cij(isp)
              endif
              rhoscf(isp,ip,3) = rhoscf(isp,ip,3) +
     .          (Dij3*cqRdi - Dij4*sqRdi)*Cij(isp)
              rhoscf(isp,ip,4) = rhoscf(isp,ip,4) +
     .          (Dij4*cqRdi + Dij3*sqRdi)*Cij(isp)
            enddo

          enddo

        enddo

C  Restore iorb for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          il = ilocal(i)
          iorb(il) = -i
        enddo

      enddo

C  Free local memory
      call de_alloc( ilocal,  name='ilocal' )
      call de_alloc( ilc,  name='ilc' )
      call de_alloc( iorb,  name='iorb' )
      call de_alloc( Clocal,  name='Clocal' )
      call de_alloc( Dlocal,  name='Dlocal' )
      call de_alloc( DlocalSp,  name='DlocalSp' )
     
      if (Parallel_Flag) then
         call de_alloc( DscfL,  name='DscfL' )
      endif

      call timer('rhoofdsp',2)
      end
