! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_dfscf
   
      public :: dfscf
    
      contains

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
C integer nspin           : spin%Grid (i.e., min(nspin,4))
C                           nspin=1 => Unpolarized, nspin=2 => polarized
C                           nspin=4 => Noncollinear spin/SOC
C integer indxua(na)      : Index of equivalent atom in unit cell
C integer isa(na)         : Species index of each atom
C integer iaorb(no)       : Atom to which orbitals belong
C integer iphorb(no)      : Index of orbital within its atom
C integer maxnd           : First dimension of Dscf
C integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
C integer listdptr(nuo)   : Pointer to start of row in listd
C integer listd(maxnd)    : List of nonzero elements of Dscf
C real*8  Dscf(maxnd,*)   : Value of nonzero elemens of density matrix
C                         : It is used in "hermitified" here
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
      use precision,     only: dp, grid_p
      use atmfuncs,      only: rcut, all_phi
      use atm_types,     only: nsmax=>nspecies
      use atomlist,      only: indxuo
      use listsc_module, only: listsc
      use mesh,          only: dxa, nsp, xdop, xdsp
      use meshphi,       only: endpht, lstpht, listp2
      use meshdscf,      only: matrixOtoM
      use meshdscf,      only: DscfL, nrowsDscfL, needDscfL
      use meshdscf,      only: listDl, listDlPtr, numdL
      use alloc,         only: re_alloc, de_alloc, alloc_default,
     $                         allocDefaults
      use parallel,      only: Nodes, Node
      use sys,           only: die
      use parallelsubs,  only: GlobalToLocalOrb

      implicit none

C  Passed arguments  

      integer, intent(in) ::
     .   ifa, istr, na, no, nuo, nuotot, np, 
     .   nspin,
     .   indxua(na), isa(na), iaorb(no), iphorb(no), 
     .   maxnd, numd(nuo), listdptr(nuo), listd(maxnd)

      real(grid_p), intent(in) ::
     .   Vscf(nsp,np,nspin), Vatm(nsp,np)

      real(dp), intent(in) ::  Datm(nuotot), dvol, VolCel
      real(dp), intent(in), target :: Dscf(:,:)

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
     .     phia(maxoa,nsp), rvol, r2sp, r2cut(nsmax)

      ! Allocate
      real(dp), pointer :: V(:,:) => null()
      real(dp), pointer :: DM_spbherm(:,:) => null()
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

      if ( size(Dscf, 2) == 8 ) then
         if (nspin /= 4) call die("Spin size inconsistency in dfscf")
         !     Prepare "spin-box hermitian" form of DM for work below
         !     We could re-use similar work in rhoofd
         !     This is the relevant part of the DM in view of the structure
         !     of the potential Vscf.
         call re_alloc(DM_spbherm,1,size(Dscf,1),1,4,"DM_spbherm")   
         DM_spbherm(:,1) = Dscf(:,1)
         DM_spbherm(:,2) = Dscf(:,2)
         DM_spbherm(:,3) = 0.5_dp * (Dscf(:,3)+Dscf(:,7))
         DM_spbherm(:,4) = 0.5_dp * (Dscf(:,4)+Dscf(:,8))
      else
         DM_spbherm => Dscf    ! Just use passed Dscf
      end if

C  Allocate buffers to store partial copies of Dscf and C
      maxc = maxval(endpht(1:np)-endpht(0:np-1))
      maxb = maxc + minb
      maxb = min( maxb, no )
      call re_alloc( C, 1, nsp, 1, maxc, 'C', 'dfscf' )
      call re_alloc( D, 0, maxb, 0, maxb, 1, nspin, 'D', 'dfscf' )
      call re_alloc( gC, 1, 3, 1, nsp, 1, maxc, 'gC', 'dfscf' )
      call re_alloc( ibc, 1, maxc, 'ibc', 'dfscf' )
      call re_alloc( iob, 0, maxb, 'iob', 'dfscf' )
      call re_alloc( xgC, 1, 9, 1, nsp, 1, maxc, 'xgC', 'dfscf' )
      call re_alloc( V, 1, nsp, 1, nspin, 'V', 'dfscf' )
      V = 0._dp

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
     .                 'DscfL', 'dfscf' )
C Redistribute Dscf to DscfL form
        call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo,
     .                   nspin, DM_spbherm, DscfL )

      endif

!     Note: Since this routine is only called to get forces AND stresses,
!     ifa = 1 and istr = 1 always. We could get rid of ix1 and ix2 and
!     unroll the relevant loops below for more efficiency.
!     Also, it might be worth unrolling some of the nsp loops below
!     if nsp is always 8.
      
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
            call die('dfscf: no slot available in D')
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
                ! i-j symmetry due to the spin-box hermiticity
                D(ib,jb,:) = DscfL(ind,:)
                D(jb,ib,:) = DscfL(ind,:)
              enddo
            else
              call GlobalToLocalOrb( iu, Node, Nodes, iul )
              do ii = 1, numd(iul)
                ind = listdptr(iul)+ii
                j = listd(ind)
                if (i.ne.iu) j = listsc( i, iu, j )
                jb = ibuff(j)
                D(ib,jb,:) = DM_spbherm(ind,:)
                D(jb,ib,:) = DM_spbherm(ind,:)
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
                call all_phi( is,+1, dxsp(:,isp), maxoa, nphiloc,
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

C     Copy potential to a double precision array

        V(1:nsp,1:nspin) = Vscf(1:nsp,ip,1:nspin)

C     Factor two for nondiagonal elements for non-collinear spin (and SO)
        if ( nspin == 4 ) then
           V(1:nsp,3:4) = 2.0_dp * V(1:nsp,3:4)
        end if

C     Loop on first orbital of mesh point
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
      call de_alloc( xgC, 'xgC', 'dfscf' )
      call de_alloc( iob, 'iob', 'dfscf' )
      call de_alloc( ibc, 'ibc', 'dfscf' )
      call de_alloc( gC, 'gC', 'dfscf' )
      call de_alloc( D, 'D', 'dfscf' )
      call de_alloc( C, 'C', 'dfscf' )
      call de_alloc( V, 'V', 'dfscf' )
      if (Parallel_Run) then
        call de_alloc( DscfL, 'DscfL', 'dfscf' )
      endif
      if (size(Dscf,2) == 8) then
         call de_alloc( DM_spbherm, 'DM_spbherm', 'dfscf' )
      endif

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

      call timer('dfscf',2)
 
      end subroutine dfscf

      end module m_dfscf
