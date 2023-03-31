! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_rhoofdsp

  implicit none

  private
  public :: rhoofdsp

contains
  subroutine rhoofdsp( no, nml, nmpl, maxnd, numd, listdptr, &
      listd, nspin, Dscf, rhoscf, nuo, nuotot, &
      iaorb, iphorb, isa, q )
! ********************************************************************
! Finds the SCF density at the mesh points from the density matrix.
! Written by P.Ordejon and J.M.Soler. May'95.
! Re-ordered so that mesh is the outer loop and the orbitals are
! handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
! Version of rhoofd that optionally uses a direct algorithm to save
! memory. Modified by J.D.Gale, November'99
! Spiral version written by V. M. Garcia-Suarez. June 2002.
! Fixed for parallel runs by Nick Papior, June 2020.
! *********************** InpUT **************************************
! integer no              : Number of basis orbitals
! integer nml             : Local numer of mesh divisions in each cell vector
! integer nmpl            : Local numer of mesh divisions
! integer maxnd           : First dimension of listD and Dscf, and
!                           maximum number of nonzero elements in
!                           any row of Dscf
! integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
! integer listdptr(nuo)   : Pointer to start of rows in listd
! integer listd(maxnd)    : List of nonzero elements in each row of Dscf
! integer nspin           : Number of spin components
! real*8  Dscf(maxnd)     : Rows of Dscf that are non-zero
! integer nuo             : Number of orbitals in unit cell locally
! integer nuotot          : Number of orbitals in unit cell in total
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! real*8 q(3)             : Wave vector for spiral configuration.
! *********************** OUTPUT **************************************
! real    rhoscf(nsp,nmpl,nspin): SCF density at mesh points
! *********************************************************************

!  Modules
    use precision, only: dp, grid_p
    use parallel, only: Nodes
    use atmfuncs,  only: rcut, all_phi
    use atm_types, only: nsmax=>nspecies
    use atomlist,  only: indxuo
    use listsc_module, only: listsc
    use mesh, only: nsp, dxa, xdop, xdsp, cmesh, nsm, meshLim
    use meshdscf
    use meshphi
    use sys, only : die
    use alloc,     only: re_alloc, de_alloc

! Argument types and dimensions
    integer :: &
        no, nml(3), nmpl, nspin, maxnd, nuo, nuotot, iaorb(*), &
        iphorb(*), isa(*), numd(nuo), listdptr(nuo), listd(maxnd)

    real(grid_p), intent(out) :: rhoscf(nsp,nmpl,nspin)

    real(dp) :: Dscf(maxnd,nspin), q(3)

    external :: memory, timer, ipack

    ! Internal variables and arrays
    integer, parameter :: minloc = 1000 ! Min buffer size for local copy of Dscf
    integer, parameter :: maxoa  = 100  ! Max # of orbitals per atom

    integer :: triang, &
        i, ia, ic, ii, ijl, il, imp, ind, ispin, io, ix, &
        iop, ip, iphi, is, isp, iu, iul, j, jc, jl, &
        last, lasta, lastop, maxloc, maxloc2, maxndl, nc, nphiloc, &
        iii(3)

    integer, pointer :: ilc(:) => null()
    integer, pointer :: ilocal(:) => null()
    integer, pointer :: iorb(:) => null()
    real(dp), pointer :: phia(:,:) => null()

    real(dp), pointer :: Clocal(:,:) => null()
    real(dp), pointer :: Dlocal(:,:) => null()

    logical :: ParallelLocal

    real(dp) :: Cij, Dij(4), &
        dxsp(3), r2cut(nsmax), r2sp, &
        xr(3), Rdi, qRdi, cqRdi(nsp), sqRdi(nsp)

    !  Start time counter
    call timer('rhoofdsp',1)

    ! Find atomic cutoff radiae
    r2cut(:) = 0.0_dp
    do i = 1,nuotot
      ia = iaorb(i)
      is = isa(ia)
      io = iphorb(i)
      r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
    end do

    ! Set algorithm logical
    ParallelLocal = (Nodes > 1)

    ! Find value of maxloc
    maxloc2 = maxval(endpht(1:nmpl)-endpht(0:nmpl-1))
    maxloc = maxloc2 + minloc
    maxloc = min( maxloc, no )
    triang  = (maxloc+1)*(maxloc+2)/2
    if ( ParallelLocal ) then
      if ( nrowsDscfL > 0 ) then
        maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
      else
        maxndl = 1
      end if
    end if

    ! Allocate local memory
    nullify(ilocal, ilc, iorb, phia, Clocal, Dlocal, DscfL)

    call re_alloc( ilocal, 1, no, 'ilocal', 'rhoofdsp' )
    call re_alloc( ilc, 1, maxloc2, 'ilc', 'rhoofdsp' )
    call re_alloc( iorb, 0, maxloc, 'iorb', 'rhoofdsp' )

    call re_alloc( Dlocal, 1, nspin, 1, triang, 'Dlocal', 'rhoofdsp' )
    call re_alloc( Clocal, 1, nsp, 1, maxloc2, 'Clocal', 'rhoofdsp' )
    if ( DirectPhi ) then
      call re_alloc(phia, 1, maxoa, 1, nsp, "phia", "rhoofdsp")
    end if

    if ( ParallelLocal ) then
      call re_alloc( DscfL, 1, maxndl, 1, nspin, 'DscfL', 'rhoofdsp' )

      ! Redistribute Dscf to DscfL form
      call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, &
          nspin, Dscf, DscfL )
    endif

    !  Initializations
    rhoscf(:,:,:) = 0.0_grid_p
    Dlocal(:,:) = 0.0_dp
    ilocal(:) = 0
    iorb(:) = 0
    last = 0

    ! Loop over grid points
    do ip = 1,nmpl

      ! Find point coordinates
      call ipack(-1,3,nml,iii,ip)
      ! Correct for processor offset
      iii(:) = iii(:) + meshLim(1,:) - 1
      do ix = 1,3
        xr(ix) = iii(1)*cmesh(ix,1) + iii(2)*cmesh(ix,2) + iii(3)*cmesh(ix,3)
      end do
 
      ! Find number of nonzero orbitals at this point
      nc = endpht(ip) - endpht(ip-1)

      !  iorb(il)>0 means that row il of Dlocal must not be overwritten
      !  iorb(il)=0 means that row il of Dlocal is empty
      !  iorb(il)<0 means that row il of Dlocal contains a valid row of
      !             Dscf, but which is not required at this point
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        il = ilocal(i)
        if ( il > 0 ) iorb(il) = i
      end do

      ! Look for required rows of Dscf not yet stored in Dlocal
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        if ( ilocal(i) == 0 ) then

          ! Look for an available row in Dlocal
          do il = 1,maxloc
            ! last runs circularly over rows of Dlocal
            last = last + 1
            if (last > maxloc) last = 1
            if (iorb(last) <= 0) goto 10
          enddo
          call die('rhoofdsp: no slot available in Dlocal')
10        continue

          !  Copy row i of Dscf into row last of Dlocal
          j = abs(iorb(last))
          if ( j /= 0 ) ilocal(j) = 0
          ilocal(i) = last
          iorb(last) = i
          il = last
          iu = indxuo(i)
          if (ParallelLocal) then
            iul = NeedDscfL(iu)
            if (i .eq. iu) then
              do ii = 1, numdl(iul)
                ind = listdlptr(iul)+ii
                j = listdl(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                do ispin = 1,nspin
                  Dlocal(ispin,ijl) = DscfL(ind,ispin)
                enddo
              enddo
            else
              do ii = 1, numdl(iul)
                ind = listdlptr(iul)+ii
                j = listsc( i, iu, listdl(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                do ispin = 1,nspin
                  Dlocal(ispin,ijl) = DscfL(ind,ispin)
                enddo
              enddo
            endif
          else
            if (i .eq. iu) then
              do ii = 1, numd(iu)
                ind = listdptr(iu)+ii
                j = listd(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                do ispin = 1,nspin
                  Dlocal(ispin,ijl) = Dscf(ind,ispin)
                enddo
              enddo
            else
              do ii = 1, numd(iu)
                ind = listdptr(iu)+ii
                j = listsc( i, iu, listd(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                do ispin = 1,nspin
                  Dlocal(ispin,ijl) = Dscf(ind,ispin)
                enddo
              enddo
            endif
          endif
        endif
      enddo

      lasta = 0
      lastop = 0

      ! Generate or retrieve phi values for all orbitals up to nc
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
                call all_phi( is, +1, dxsp, maxoa, nphiloc, phia(:,isp) )
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

      ! Loop on first orbital of mesh point
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        il = ilocal(i)
        iu = indxuo(i)
        ia = iaorb(i)
        is = isa(ia)
        iop = listp2(imp)

        ! Calculate spiral phase
        do isp = 1, nsp
          qRdi = 0._dp
          do ix = 1, 3
            Rdi = xr(ix) + xdsp(ix,isp) - xdop(ix,iop) + dxa(ix,ia)
            qRdi = qRdi + q(ix) * Rdi
          end do
          cqRdi(isp) = cos(qRdi)
          sqRdi(isp) = sin(qRdi)
        end do

        ! Loop on second orbital of mesh point, here only jc < ic
        do jc = 1, ic - 1
          jl = ilc(jc)

          ijl = idx_ijl(il,jl)

          Dij(1) = 2*Dlocal(1,ijl)
          Dij(2) = 2*Dlocal(2,ijl)
          Dij(3) = Dlocal(3,ijl)
          Dij(4) = Dlocal(4,ijl)

          ! Loop over sub-points
          do isp = 1, nsp
            Cij = Clocal(isp,ic) * Clocal(isp,jc)
            rhoscf(isp,ip,1) = rhoscf(isp,ip,1) + Dij(1)*Cij
            rhoscf(isp,ip,2) = rhoscf(isp,ip,2) + Dij(2)*Cij
            rhoscf(isp,ip,3) = rhoscf(isp,ip,3) + &
                (Dij(3)*cqRdi(isp) + Dij(4)*sqRdi(isp))*Cij
            rhoscf(isp,ip,4) = rhoscf(isp,ip,4) + &
                (Dij(4)*cqRdi(isp) - Dij(3)*sqRdi(isp))*Cij
          end do

        enddo

        ! onsite orbital
        ijl = idx_ijl(il,il)
        Dij(1) = Dlocal(1,ijl)
        Dij(2) = Dlocal(2,ijl)
        Dij(3) = Dlocal(3,ijl)
        Dij(4) = Dlocal(4,ijl)
        do isp = 1, nsp
          Cij = Clocal(isp,ic) * Clocal(isp,ic)
          rhoscf(isp,ip,1) = rhoscf(isp,ip,1) + Dij(1)*Cij
          rhoscf(isp,ip,2) = rhoscf(isp,ip,2) + Dij(2)*Cij
          rhoscf(isp,ip,3) = rhoscf(isp,ip,3) + &
              (Dij(3)*cqRdi(isp) + Dij(4)*sqRdi(isp))*Cij
          rhoscf(isp,ip,4) = rhoscf(isp,ip,4) + &
              (Dij(4)*cqRdi(isp) - Dij(3)*sqRdi(isp))*Cij
        end do

        ! Loop on second orbital of mesh point (here only jc > ic)
        do jc = ic + 1, nc
          jl = ilc(jc)
          ijl = idx_ijl(il,jl)
          Dij(3) = Dlocal(3,ijl)
          Dij(4) = Dlocal(4,ijl)
          do isp = 1, nsp
            Cij = Clocal(isp,ic) * Clocal(isp,jc)
            rhoscf(isp,ip,3) = rhoscf(isp,ip,3) + &
                (Dij(3)*cqRdi(isp) + Dij(4)*sqRdi(isp))*Cij
            rhoscf(isp,ip,4) = rhoscf(isp,ip,4) + &
                (Dij(4)*cqRdi(isp) - Dij(3)*sqRdi(isp))*Cij
          end do
        end do

      enddo

      ! Restore iorb for next point
      do imp = 1+endpht(ip-1), endpht(ip)
        i = lstpht(imp)
        il = ilocal(i)
        iorb(il) = -i
      enddo

    enddo

    ! Free local memory
    call de_alloc( ilocal, 'ilocal', 'rhoofdsp' )
    call de_alloc( ilc,'ilc', 'rhoofdsp' )
    call de_alloc( iorb,'iorb', 'rhoofdsp' )
    call de_alloc( Clocal, 'Clocal', 'rhoofdsp' )
    call de_alloc( Dlocal, 'Dlocal', 'rhoofdsp' )
    if ( DirectPhi ) then
      call de_alloc( phia, 'phia', 'vmatsp' )
    end if

    if (ParallelLocal) then
      call de_alloc( DscfL, 'DscfL', 'rhoofdsp' )
    endif

    call timer('rhoofdsp',2)

  contains

    ! In any case will the compiler most likely inline this
    ! small routine. So it should not pose any problem.
    pure function idx_ijl(i,j) result(ij)
      integer, intent(in) :: i,j
      integer :: ij
      if ( i > j ) then
        ij = i * (i + 1)/2 + j + 1
      else
        ij = j * (j + 1)/2 + i + 1
      end if
    end function idx_ijl

  end subroutine rhoofdsp

end module m_rhoofdsp
