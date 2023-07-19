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
      subroutine dfscf( ifa, istr, na, no, nuo, nuotot, np, nspin,
     .                  indxua, isa, iaorb, iphorb,
     .                  maxnd, numd, listdptr, listd, Dscf, Datm,
     .                  Vscf, Vatm, dvol, VolCel, Fal, Stressl )

C ********************************************************************
C Adds the SCF contribution to atomic forces and stress.
C Written by P.Ordejon, J.M.Soler, and J.Gale.
C Last modification by J.M.Soler, October 2000.
C *********************** INPUT **************************************
C integer ifa             : Are forces required? (1=yes,0=no)
C integer istr            : Is stress required? (1=yes,0=no)
C integer na              : Number of atoms
C integer no              : Number of basis orbitals
C integer nuo             : Number of orbitals in unit cell (local)
C integer nuotot          : Number of orbitals in unit cell (global)
C integer np              : Number of mesh points (total is nsp*np)
C integer nspin           : Number of spin components
C integer indxua(na)      : Index of equivalent atom in unit cell
C integer isa(na)         : Species index of each atom
C integer iaorb(no)       : Atom to which orbitals belong
C integer iphorb(no)      : Index of orbital within its atom
C integer maxnd           : First dimension of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of row in listd
C integer listd(maxnd)    : List of nonzero elements of Dscf
C real*8  Dscf(maxnd,nspin): Value of nonzero elemens of density matrix
C real*8  Datm(nuotot)    : Occupations of basis orbitals in free atom
C real*4  Vscf(nsp,np,nspin): Value of SCF potential at the mesh points
C real*4  Vatm(nsp,np)    : Value of Harris potential (Hartree potential
C                           of sum of atomic desities) at mesh points
C                           Notice single precision of Vscf and Vatm
C real*8  dvol            : Volume per mesh point
C real*8  VolCel          : Unit cell volume
C *********************** INPUT and OUTPUT ****************************
C real*8  Fal(3,*)       : Atomic forces (contribution added on output)
C real*8  Stressl(3,3)    : Stress tensor
C *********************************************************************
C    6  10        20        30        40        50        60        7072

C  Modules
      use precision, only: dp, grid_p
      use atmfuncs, only: rcut, all_phi
      use atm_types, only: nsmax=>nspecies
      use atomlist, only: indxuo
      use listsc_module, only: listsc
      use mesh, only: dxa, nsp, xdop, xdsp
      use meshphi, only: endpht, lstpht, listp2
      use meshdscf, only: DscfL, nrowsDscfL, needDscfL
      use meshdscf, only: listDl, listDlPtr, numdL
      use alloc,    only: re_alloc, de_alloc, alloc_default,
     $                    allocDefaults
      use parallel, only: Nodes
      use sys,      only: die

      implicit none

C  Passed arguments
      integer, intent(in) ::
     .   ifa, istr, na, no, nuo, nuotot, np, nspin,  
     .   indxua(na), isa(na), iaorb(no), iphorb(no), 
     .   maxnd, numd(nuo), listdptr(nuo), listd(maxnd)

      real(grid_p), intent(in) ::
     .   Vscf(nsp,np,nspin), Vatm(nsp,np)

      real(dp), intent(in) ::
     .   Datm(nuotot), Dscf(maxnd,nspin), dvol, VolCel

      real(dp), intent(inout) :: Fal(3,*), Stressl(9)

C Internal variables
      integer, parameter ::
     .   minb  = 100,  ! Min buffer size for local copy of Dscf
     .   maxoa = 100   ! Max # of orbitals per atom
      integer
     .   i, ia, ib, ibuff(no), ic, ii, imp, ind, iop, ip, iphi, io,
     .   is, isp, ispin, iu, iua, iul, ix, ix1, ix2, iy,
     .   j, jb, jc, last, lasta, lastop, maxb, maxc, maxndl,
     .   nc, nphiloc
      real(dp)
     .   CD(nsp), CDV(nsp), DF(12), Dji, dxsp(3,nsp),  
     .   gCi(12,nsp), grada(3,maxoa,nsp),
     .   phia(maxoa,nsp), rvol, r2sp, r2cut(nsmax), V(nsp,nspin)
!
      integer, pointer, save ::  ibc(:), iob(:)
      real(dp), pointer, save :: C(:,:), D(:,:,:),
     .                           gC(:,:,:), xgC(:,:,:)
      logical ::           Parallel_Run, nullified=.false.
      type(allocDefaults) oldDefaults

C  Start time counter
      call timer('dfscf',1)

C  Nullify pointers
      if (.not.nullified) then
        nullify( C, D, gC, ibc, iob, xgC )
        nullified = .true.
      end if

C  Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults,
     .                    copy=.false., shrink=.false.,
     .                    imin=1, routine='dfscf' )
      
C  Allocate buffers to store partial copies of Dscf and C
      maxc = maxval(endpht(1:np)-endpht(0:np-1))
      maxb = maxc + minb
      maxb = min( maxb, no )
      call re_alloc( C,   1,nsp,  1,maxc,          name='C'   )
      call re_alloc( D,   0,maxb, 0,maxb, 1,nspin, name='D'   )
      call re_alloc( gC,  1,3,    1,nsp,  1,maxc,  name='gC'  )
      call re_alloc( ibc, 1,maxc,                  name='ibc' )
      call re_alloc( iob, 0,maxb,                  name='iob' )
      call re_alloc( xgC, 1,9,    1,nsp,  1,maxc,  name='xgC' )

C  Set logical that determines whether we need to use parallel or serial mode
      Parallel_Run = (Nodes.gt.1)

C  If parallel, allocate temporary storage for Local Dscf
      if (Parallel_Run) then
        if (nrowsDscfL.gt.0) then
          maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
        else
          maxndl = 1
        endif
        call re_alloc( DscfL, 1, maxndl, 1, nspin,
     &                 name='DscfL',  routine='dfscf' )
C Redistribute Dscf to DscfL form
        call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, nuotot,
     .                   nspin, Dscf, DscfL )

      endif

C  Find range of a single array to hold force and stress derivatives
C  Range 1-3 for forces
      if (ifa.eq.1) then
        ix1 = 1
      else
        ix1 = 4
      end if
C  Range 4-12 for stress
      if (istr.eq.1) then
        ix2 = 12
      else
        ix2 = 3
      end if

C  Initialise variables
      D(:,:,:) = 0.0_dp
      ibuff(:) = 0
      iob(:) = 0
      last = 0

C  Find atomic cutoff radii
      r2cut(:) = 0.0_dp
      do i = 1,nuotot
        ia = iaorb(i)
        is = isa(ia)
        io = iphorb(i)
        r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
      enddo

C  Evaluate constants
      rvol = 1.0_dp / VolCel

C  Loop over grid points
      do ip = 1,np

C  Find number of nonzero orbitals at this point
        nc = endpht(ip) - endpht(ip-1)

C  iob(ib)>0 means that row ib of D must not be overwritten
C  iob(ib)=0 means that row ib of D is empty
C  iob(ib)<0 means that row ib of D contains a valid row of 
C             Dscf, but which is not required at this point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          ib = ibuff(i)
          if (ib.gt.0) iob(ib) = i
        enddo

C  Look for required rows of Dscf not yet stored in D
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ibuff(i) .eq. 0) then

C  Look for an available row in D
            do ib = 1,maxb
C  last runs circularly over rows of D
              last = last + 1
              if (last .gt. maxb) last = 1
              if (iob(last) .le. 0) goto 10
            enddo
            call die('rhoofd: no slot available in D')
   10       continue

C  Copy row i of Dscf into row last of D
            j = abs(iob(last))
            if (j.ne.0) ibuff(j) = 0
            ibuff(i) = last
            iob(last) = i
            ib = last
            iu = indxuo(i)
            if (Parallel_Run) then
              iul = NeedDscfL(iu)
              do ii = 1, numdl(iul)
                ind = listdlptr(iul)+ii
                j = listdl(ind)
                if (i.ne.iu) j = listsc( i, iu, j )
                jb = ibuff(j)
                D(ib,jb,1:nspin) = DscfL(ind,1:nspin)
                D(jb,ib,1:nspin) = DscfL(ind,1:nspin)
              enddo
            else
              do ii = 1, numd(iu)
                ind = listdptr(iu)+ii
                j = listd(ind)
                if (i.ne.iu) j = listsc( i, iu, j )
                jb = ibuff(j)
                D(ib,jb,1:nspin) = Dscf(ind,1:nspin)
                D(jb,ib,1:nspin) = Dscf(ind,1:nspin)
              enddo
            endif
          endif
          ibc(ic) = ibuff(i)
        enddo

C  Restore iob for next point
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          ib = ibuff(i)
          iob(ib) = -i
        enddo

C  Calculate all phi values and derivatives at all subpoints
        lasta = 0
        lastop = 0
        do ic = 1,nc

          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iphi = iphorb(i)
          is = isa(ia)
          iua = indxua(ia)
          iop = listp2(imp)
          if (ia.ne.lasta .or. iop.ne.lastop) then
            lasta = ia
            lastop = iop
            do isp = 1,nsp
              dxsp(1:3,isp) = xdop(1:3,iop)+xdsp(1:3,isp)-dxa(1:3,ia)
              r2sp = dxsp(1,isp)**2 + dxsp(2,isp)**2 + dxsp(3,isp)**2
              if (r2sp.lt.r2cut(is)) then
                call all_phi( is,+1, dxsp(:,isp), nphiloc,
     .                        phia(:,isp), grada(:,:,isp))
              else
                phia(:,isp) = 0.0_dp
                grada(1:3,:,isp) = 0.0_dp
              endif
            enddo
          endif
          C(1:nsp,ic) = phia(iphi,1:nsp)
          gC(1:3,1:nsp,ic) = grada(1:3,iphi,1:nsp)

C  If stress required. Generate stress derivatives
          if (istr.eq.1) then
            do isp = 1,nsp
              ii = 0
              do ix = 1,3
                do iy = 1,3
                  ii = ii + 1
                  xgC(ii,isp,ic) = dxsp(iy,isp) * gC(ix,isp,ic) * rvol
                enddo
              enddo
            enddo
          endif
        enddo

C  Copy potential to a double precision array
        V(1:nsp,1:nspin) = Vscf(1:nsp,ip,1:nspin)

C  Factor two for nondiagonal elements for non-collinear spin
        V(1:nsp,3:nspin) = 2.0_dp * V(1:nsp,3:nspin)

C  Loop on first orbital of mesh point
        do ic = 1,nc

          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          iu = indxuo(i)
          ia = iaorb(i)
          iua = indxua(ia)
          ib = ibc(ic)

C  Copy force and stress gradients to a single array
          gCi(1:3, 1:nsp) =  gC(1:3,1:nsp,ic)
          gCi(4:12,1:nsp) = xgC(1:9,1:nsp,ic)

C  Find Sum_{j,spin}(Vscf*Dij*Cj) for every subpoint
C  Some loops are not done using f90 form as this
C  leads to much slower execution on machines with stupid f90
C  compilers at the moment
          CDV(1:nsp) = 0.0_dp
          do ispin = 1,nspin

C  Loop on second orbital of mesh point
            CD(1:nsp) = 0.0_dp
            do jc = 1,nc
              jb = ibc(jc)
              Dji = D(jb,ib,ispin)
              do isp = 1,nsp
                CD(isp) = CD(isp) + C(isp,jc) * Dji
              end do
            end do

            do isp = 1,nsp
              CDV(isp) = CDV(isp) + CD(isp) * V(isp,ispin)
            end do
          enddo
          do isp = 1,nsp
            CDV(isp) = CDV(isp) - 
     .                 C(isp,ic) * Datm(iu) * Vatm(isp,ip)
            CDV(isp) = 2.0_dp * dVol * CDV(isp)
          end do

C  Add 2*Dscf_ij*<Cj|Vscf|gCi> to forces
          do ix = ix1,ix2
            DF(ix) = 0.0_dp
            do isp = 1, nsp
              DF(ix) = DF(ix) + gCi(ix,isp) * CDV(isp)
            enddo
          enddo

C  Add force and stress to output arrays
          if (ifa.eq.1) then
            Fal(1:3,iua) = Fal(1:3,iua) + DF(1:3)
          endif
          if (istr.eq.1) then
            Stressl(1:9) = Stressl(1:9) + DF(4:12)
          endif

C  End of first orbital loop
        enddo

C  End of mesh point loop
      enddo
  
C  Deallocate local memory
      call de_alloc( xgC, name='xgC'  )
      call de_alloc( iob, name='iob' )
      call de_alloc( ibc, name='ibc'  )
      call de_alloc( gC,  name='gC'   )
      call de_alloc( D,   name='D'    )
      call de_alloc( C,   name='C'    )
      if (Parallel_Run) then
        call de_alloc( DscfL,  name='DscfL' )
      endif

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

      call timer('dfscf',2)
      end
