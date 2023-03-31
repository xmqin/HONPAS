! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_vmatsp

  implicit none

  private
  public :: vmatsp

contains

  subroutine vmatsp( no, nml, nmpl, dvol, nspin, V, nvmax, &
      numVs, listVsptr, listVs, Vs, &
      nuo, nuotot, iaorb, iphorb, isa, q )

! ********************************************************************
! Finds the matrix elements of the potential.
! First version written by P.Ordejon.
! Name and interface modified by J.M.Soler. May'95.
! Re-ordered so that mesh is the outer loop and the orbitals are
! handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
! Version of vmat that use a direct algorithm to save memory.
! Modified by J.D.Gale, November'99
! Spiral version written by V. M. Garcia-Suarez. June 2002.
! Fixed for parallel runs by Nick Papior, June 2020.
! *********************** INPUT **************************************
! integer no              : Number of basis orbitals
! integer nml             : Local numer of mesh divisions in each cell vector
! integer nmpl            : Local numer of mesh divisions
! real*8  dvol            : Volume per mesh point
! integer nspin           : Number of spin components
! real*4  V(nsp,nmpl,nspin): Value of the potential at the mesh points
! integer nvmax           : First dimension of listV and Vs, and maxim.
!                           number of nonzero elements in any row of Vs
! integer numVs(nuo)      : Number of non-zero elements in a row of Vs
! integer listVsptr(nuo)  : Pointer to the start of rows in listVs
! integer listVs(nvmax)   : List of non-zero elements of Vs
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! real*8 q(3)             : Wave vector for spiral configuration.
! ******************** INPUT AND OUTPUT *******************************
! real*8  Vs(nvmax,nspin) : Value of nonzero elements in each row
!                           of Vs to which the potential matrix
!                           elements are summed up
! *********************************************************************

!  Modules
    use precision, only: dp, grid_p
    use parallel, only: Nodes
    use atmfuncs,  only: rcut, all_phi
    use atm_types, only: nsmax=>nspecies
    use atomlist,  only: indxuo
    use listsc_module, only: listsc
    use mesh, only: dxa, nsp, xdop, xdsp, cmesh, idop, nsm, meshLim
    use meshdscf
    use meshphi
    use alloc, only: re_alloc, de_alloc

! Argument types and dimensions
    integer :: no, nml(3), nmpl, nvmax, nuo, nuotot, iaorb(*), nspin, &
        iphorb(*), isa(*), numVs(nuo), listVsptr(nuo), listVs(nvmax)
    real(grid_p), intent(in) :: V(nsp,nmpl,nspin)
    real(dp) :: dvol, Vs(nvmax,nspin), q(3)

! Internal variables and arrays
    integer, parameter :: minloc = 1000 ! Min buffer size for local copy of Dscf
    integer, parameter :: maxoa  = 100  ! Max # of orbitals per atom
    integer, parameter :: nsd = 2  ! spin-box size
    integer :: &
        i, ia, ic, ii, ijl, il, imp, ind, iop, ip, iphi, io, &
        is, isp, ispin, iu, iul, ix, j, jc, jl, triang, &
        last, lasta, lastop, maxloc, maxloc2, nc, nlocal,  &
        nphiloc, nvmaxl, iii(3)

    logical :: ParallelLocal
    real(dp) :: Vij(4), dxsp(3), &
        r2cut(nsmax), r2sp, xr(3), Rdi, qRdi, cqRdi, sqRdi

    integer, pointer :: ilc(:) => null()
    integer, pointer :: ilocal(:) => null()
    integer, pointer :: iorb(:) => null()
    real(dp), pointer :: phia(:,:) => null()
    real(dp), pointer :: VClocal(:,:) => null()
    real(dp), pointer :: Clocal(:,:) => null()
    real(dp), pointer :: Vlocal(:,:) => null()
    real(dp), pointer :: VlocalSp(:,:,:) => null()

    ! Start time counter
    call timer('vmatsp',1)

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
        nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
      else
        nvmaxl = 1
      end if
    end if

    ! Allocate local memory
    nullify(ilocal, ilc, iorb, phia, Vlocal, VlocalSp, Clocal, VClocal, DscfL)

    call re_alloc( ilocal, 1, no, 'ilocal', 'vmatsp' )
    call re_alloc( ilc, 1, maxloc2, 'ilc', 'vmatsp' )
    call re_alloc( iorb, 0, maxloc, 'iorb', 'vmatsp' )

    call re_alloc( Vlocal, 1, nsd, 1, triang, 'Vlocal', 'vmatsp' )
    call re_alloc( VlocalSp, 1, nsd, 0, maxloc, 0, maxloc, 'VlocalSp', 'vmatsp' )
    call re_alloc( Clocal, 1, nsp, 1, maxloc2, 'Clocal', 'vmatsp' )
    call re_alloc( VClocal, 1, nspin, 1, nsp, 'VClocal', 'vmatsp' )
    if ( DirectPhi ) then
      call re_alloc(phia, 1, maxoa, 1, nsp, "phia", "vmatsp")
    end if
    if ( ParallelLocal ) then
      call re_alloc( DscfL, 1, nvmaxl, 1, nspin, 'DscfL', 'vmatsp' )
      DscfL(1:nvmaxl,1:nspin) = 0.0_dp
    end if

    ! Full initializations done only once
    ilocal(1:no) = 0
    iorb(0:maxloc) = 0
    Vlocal(:,:) = 0.0_dp
    VlocalSp(:,:,:) = 0.0_dp
    last = 0

    ! Loop over (big) grid points
    do ip = 1 , nmpl

      ! Find point coordinates
      call ipack(-1,3,nml,iii,ip)
      ! Correct for processor offset
      iii(:) = iii(:) + meshLim(1,:) - 1
      do ix = 1,3
        xr(ix) = iii(1)*cmesh(ix,1) + iii(2)*cmesh(ix,2) + &
            iii(3)*cmesh(ix,3)
      end do

      ! Find number of nonzero orbitals at this point
      nc = endpht(ip) - endpht(ip-1)

      ! Find new required size of Vlocal
      nlocal = last
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        if ( ilocal(i) == 0 ) nlocal = nlocal + 1
      end do

      ! If overflooded, add Vlocal to Vs and reinitialize it
      if ( nlocal > maxloc .and. last > 0 ) then
        if ( ParallelLocal ) then
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            iul = NeedDscfL(iu)
            if ( i == iu ) then
              do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j = listdl(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                do ispin = 1, nsd
                  DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * Vlocal(ispin,ijl)
                  DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
                end do
              end do
            else
              do ii = 1, numdl(iul)
                ind = listdlptr(iul) + ii
                j = LISTSC( i, iu, listdl(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il, jl)
                do ispin = 1,nsd
                  DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * Vlocal(ispin,ijl)
                  DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
                end do
              end do
            end if
          end do
        else ! ParallelLocal
          do il = 1,last
            i = iorb(il)
            iu = indxuo(i)
            !iul = iu ! Nodes == 1
            if ( i == iu ) then
              do ii = 1, numVs(iu)
                ind = listVsptr(iu) + ii
                j = listVs(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il, jl)
                do ispin = 1, nsd
                  Vs(ind,ispin) = Vs(ind,ispin) + dVol * Vlocal(ispin,ijl)
                  Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
                end do
              end do
            else
              do ii = 1, numVs(iu)
                ind = listVsptr(iu) + ii
                j = LISTSC( i, iu, listVs(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il, jl)
                do ispin = 1,nsd
                  Vs(ind,ispin) = Vs(ind,ispin) + dVol * Vlocal(ispin,ijl)
                  Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
                end do
              end do
            end if
          end do
        end if
        ! Reset local arrays
        do i = 1, last
          ilocal(iorb(i)) = 0
          iorb(i) = 0
        end do
        ijl = (last+1)*(last+2)/2
        Vlocal(1:nsd,1:ijl) = 0.0_dp
        VlocalSp(1:nsd,0:last,0:last) = 0.0_dp
        last = 0
      endif

      !  Look for required orbitals not yet in Vlocal
      if (nlocal > last) then
        do ic = 1, nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if ( ilocal(i) == 0 ) then
            last = last + 1
            ilocal(i) = last
            iorb(last) = i
          end if
        end do
      end if

      if ( DirectPhi ) then
        lasta = 0
        lastop = 0

        !  Generate or retrieve phi values for all orbitals up to nc
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          il = ilocal(i)
          ia = iaorb(i)
          iop = listp2(imp)
          ilc(ic) = il

          if ( ia /= lasta .or. iop /= lastop ) then
            lasta = ia
            lastop = iop
            is = isa(ia)
            do isp = 1,nsp
              do ix = 1,3
                dxsp(ix) = xdsp(ix,isp) + xdop(ix,iop) - dxa(ix,ia)
              enddo
              r2sp = sum(dxsp**2)
              if ( r2sp < r2cut(is) ) then
                call all_phi( is, +1, dxsp, maxoa, nphiloc, phia(:,isp) )
              else
                phia(:,isp) = 0.0_dp
              end if
            end do
          end if
          iphi = iphorb(i)
          do isp = 1,nsp
            Clocal(isp,ic) = phia(iphi,isp)
          end do

        end do

      else ! DirectPhi
        !  Generate or retrieve phi values for all orbitals up to nc
        do ic = 1,nc
          imp = endpht(ip-1) + ic
          il = ilocal(lstpht(imp))
          ilc(ic) = il

          Clocal(:,ic) = phi(:,imp)
        end do

      end if

      ! Loop on first orbital of mesh point
      do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        il = ilocal(i)
        iu = indxuo(i)
        ia = iaorb(i)
        is = isa(ia)
        iop = listp2(imp)

        ! Pre-multiply V and Clocal(,ic)
        do isp = 1,nsp
          ! Calculate spiral phase
          qRdi = 0._dp
          do ix = 1, 3
            Rdi = xr(ix) + xdsp(ix,isp) - xdop(ix,iop) + dxa(ix,ia)
            qRdi = qRdi + q(ix) * Rdi
          end do
          cqRdi = cos(qRdi)
          sqRdi = sin(qRdi)

          VClocal(1,isp) = V(isp,ip,1) * Clocal(isp,ic)
          VClocal(2,isp) = V(isp,ip,2) * Clocal(isp,ic)
          VClocal(3,isp) = (V(isp,ip,3)*cqRdi - V(isp,ip,4)*sqRdi) * Clocal(isp,ic)
          VClocal(4,isp) = (V(isp,ip,4)*cqRdi + V(isp,ip,3)*sqRdi) * Clocal(isp,ic)
        enddo

        ! Loop on second orbital of mesh point (NOT only for jc.le.ic)
        do jc = 1,nc
          jl = ilc(jc)

          ! Loop over sub-points
          Vij(:) = 0.0_dp
          do isp = 1,nsp
            Vij(:) = Vij(:) + VClocal(:,isp) * Clocal(isp,jc)
          end do

          ijl = idx_ijl(il, jl)

          ! Here we need the spiral contribution for jc > ic
          ! So an explicit check is required
          if ( jc == ic ) then
            Vlocal(1,ijl) = Vlocal(1,ijl) + Vij(1)
            Vlocal(2,ijl) = Vlocal(2,ijl) + Vij(2)
          else if ( jc < ic ) then
            if ( il == jl ) then
              Vlocal(1,ijl) = Vlocal(1,ijl) + 2.0_dp*Vij(1)
              Vlocal(2,ijl) = Vlocal(2,ijl) + 2.0_dp*Vij(2)
            else
              Vlocal(1,ijl) = Vlocal(1,ijl) + Vij(1)
              Vlocal(2,ijl) = Vlocal(2,ijl) + Vij(2)
            end if
          end if

          VlocalSp(1,jl,il) = VlocalSp(1,jl,il) + Vij(3)
          VlocalSp(2,jl,il) = VlocalSp(2,jl,il) + Vij(4)

        end do

      end do
    end do

    ! Add final Vlocal to Vs
    if ( ParallelLocal .and. last > 0 ) then

      do il = 1,last
        i = iorb(il)
        iu = indxuo(i)
        iul = NeedDscfL(iu)
        if ( i == iu ) then
          do ii = 1, numdl(iul)
            ind = listdlptr(iul) + ii
            j = listdl(ind)
            jl = ilocal(j)
            ijl = idx_ijl(il, jl)
            do ispin = 1, nsd
              DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * Vlocal(ispin,ijl)
              DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
            end do
          end do
        else
          do ii = 1, numdl(iul)
            ind = listdlptr(iul)+ii
            j = LISTSC( i, iu, listdl(ind) )
            jl = ilocal(j)
            ijl = idx_ijl(il, jl)
            do ispin = 1, nsd
              DscfL(ind,ispin) = DscfL(ind,ispin) + dVol * Vlocal(ispin,ijl)
              DscfL(ind,ispin+2) = DscfL(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
            end do
          end do
        end if
      end do

    else if ( last > 0 ) then ! ParallelLocal

      do il = 1 , last
        i = iorb(il)
        iu = indxuo(i)
        if ( i == iu ) then
          do ii = 1, numVs(iu)
            ind = listVsptr(iu) + ii
            j = listVs(ind)
            jl = ilocal(j)
            ijl = idx_ijl(il, jl)
            do ispin = 1 , nsd
              Vs(ind,ispin) = Vs(ind,ispin) + dVol * Vlocal(ispin,ijl)
              Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
            end do
          end do
        else
          do ii = 1, numVs(iu)
            ind = listVsptr(iu) + ii
            j = LISTSC( i, iu, listVs(ind) )
            jl = ilocal(j)
            ijl = idx_ijl(il, jl)
            do ispin = 1,nsd
              Vs(ind,ispin) = Vs(ind,ispin) + dVol * Vlocal(ispin,ijl)
              Vs(ind,ispin+2) = Vs(ind,ispin+2) + dVol * VlocalSp(ispin,jl,il)
            end do
          end do
        end if
      end do
    end if

    !  Free local memory
    call de_alloc( ilocal, 'ilocal', 'vmatsp' )
    call de_alloc( ilc, 'ilc', 'vmatsp' )
    call de_alloc( iorb, 'iorb', 'vmatsp' )
    call de_alloc( Vlocal, 'Vlocal', 'vmatsp' )
    call de_alloc( VlocalSp, 'VlocalSp', 'vmatsp' )
    call de_alloc( Clocal, 'Clocal', 'vmatsp' )
    call de_alloc( VClocal, 'VClocal', 'vmatsp' )
    if ( DirectPhi ) then
      call de_alloc( phia, 'phia', 'vmatsp' )
    end if

    if ( ParallelLocal ) then
      ! Redistribute Hamiltonian from mesh to orbital based distribution
      call matrixMtoO( nvmaxl, nvmax, numVs, listVsptr, nuo, &
          nspin, DscfL, Vs)
      ! Free memory
      call de_alloc( DscfL, 'DscfL', 'vmatsp' )
    endif

    call timer('vmatsp',2)

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

  end subroutine vmatsp

end module m_vmatsp
