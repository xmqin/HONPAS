! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_delk

  implicit none

  private

  public :: delk

contains

  subroutine delk( iexpikr, no, np, dvol, nvmax, &
       numVs, listVsptr, listVs, &
       nuo, nuotot, iaorb, iphorb, isa )

! ********************************************************************
! Finds the matrix elements of a plane wave
! < \phi_{\mu} | exp^( iexpikr * i * \vec{k} \cdot \vec{r} ) | \phi_{\nu} >
! First version written by J. Junquera in Feb. 2008
! Adapted from an existing version of vmat after the parallelization 
! designed by BSC in July 2011.
! *********************** INPUT **************************************
! integer iexpikr         : Prefactor of the dot product between the
!                           the k-vector and the position-vector in exponent.
!                           it might be +1 or -1
! integer no              : Number of basis orbitals
! integer np              : Number of columns in C (local)
! real*8  dvol            : Volume per mesh point
! integer nvmax           : First dimension of listV and Vs, and maxim.
!                           number of nonzero elements in any row of delkmat
! integer numVs(nuo)      : Number of non-zero elements in a row of delkmat
! integer listVsptr(nuo)  : Pointer to the start of rows in listVs
! integer listVs(nvmax)   : List of non-zero elements of delkmat
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! *****************************  OUTPUT *******************************
! complex*16  delkmat(nvmax) : Value of nonzero elements in each row
!                              of the matrix elements of exp(i*\vec{k}\vec{r})
!                              this variable is defined in the module
!                              m_dimensionsefield (file m_dimefield.f90)
! *********************************************************************

!  Modules
    use precision,     only: dp, grid_p
    use atmfuncs,      only: rcut, all_phi
    use atm_types,     only: nsmax=>nspecies
    use atomlist,      only: indxuo, indxua
    use siesta_geom,   only: xa
    use listsc_module, only: listsc
    use mesh,          only: dxa, nsp, xdop, xdsp, ne, nem
    use mesh,          only: cmesh, ipa, idop, nmsc, iatfold
    use mesh,          only: meshLim
    use meshdscf,      only: matrixMtoOC
    use meshdscf,      only: needdscfl, listdl, numdl, nrowsdscfl, listdlptr
    use meshphi,       only: directphi, endpht, lstpht, listp2, phi
    use parallel,      only: Nodes, node
    use alloc,         only: re_alloc, de_alloc
    use parallelsubs,  only: GlobalToLocalOrb
    use m_planewavematrixvar, only: delkmat, wavevector
#ifdef MPI
    use mpi_siesta
#endif
#ifdef _OPENMP
    use omp_lib
#endif

! Argument types and dimensions
    integer :: iexpikr, no, np, nvmax, nuo, nuotot
    integer :: iaorb(*), iphorb(*), isa(*), numVs(nuo)
    integer :: listVsptr(nuo), listVs(nvmax)
    real(dp) :: dvol
! Internal variables and arrays
    integer, parameter :: minloc = 1000 ! Min buffer size
    integer, parameter :: maxoa  = 100  ! Max # of orb/atom
    integer :: i, ia, ic, ii, ijl, il, imp, ind, iop
    integer :: ip, iphi, io, is, isp, iu, iul
    integer :: ix, j, jc, jl, last, lasta, lastop
    integer :: maxloc, maxloc2, nc, nlocal, nphiloc
    integer :: nvmaxl, triang, lenx, leny, lenz,lenxy
    logical :: ParallelLocal
    real(dp) :: Vij(2), r2sp, dxsp(3), VClocal(2,nsp)

    integer, pointer :: ilc(:), ilocal(:), iorb(:)
    real(dp), pointer :: DscfL(:,:),  t_DscfL(:,:,:), Clocal(:,:)
    real(dp), pointer :: Vlocal(:,:), phia(:,:), r2cut(:)
    integer :: NTH, TID

! Variables to compute the matrix element of the plane wave
! (not included in the original vmat subroutine)
    integer :: irel, iua, irealim, inmp(3)
    real(dp) :: kxij, aux(2), dist(3), kpoint(3)
    real(dp) :: dxp(3), displaat(3)
    real(dp) :: dxpgrid(3,nsp)
    complex(dp), pointer :: delkmats(:), t_delkmats(:,:)

#ifdef _TRACE_
    integer :: MPIerror
#endif

#ifdef DEBUG
    call write_debug( '    PRE delk' )
#endif

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 4 )
#endif

!   Start time counter
    call timer('delk',1)

!   Initialize the matrix elements of exp(i*\vec{k} \vec{r})
    delkmat(:) = 0.0_dp
    kpoint(:)  = dble(iexpikr) * wavevector(:)

!! For debugging
!      kpoint(:) = 0.0_dp
!! End debugging

!   Find atomic cutoff radii
    nullify(r2cut)
    call re_alloc( r2cut, 1, nsmax, 'r2cut', 'delk' )
    r2cut = 0.0_dp
    do i = 1,nuotot
       ia = iaorb(i)
       is = isa(ia)
       io = iphorb(i)
       r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
    end do

!   Set algorithm logical
    ParallelLocal = (Nodes > 1)
    lenx  = meshLim(2,1) - meshLim(1,1) + 1
    leny  = meshLim(2,2) - meshLim(1,2) + 1
    lenz  = meshLim(2,3) - meshLim(1,3) + 1
    lenxy = lenx*leny

!   Find value of maxloc
    maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
    maxloc  = maxloc2 + minloc
    maxloc  = min( maxloc, no )
    triang  = (maxloc+1)*(maxloc+2)/2
    if ( ParallelLocal ) then
       if ( nrowsDscfL > 0 ) then
          nvmaxl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
       else
          nvmaxl = 1
       end if
    end if

!   Allocate local memory
!$OMP parallel default(shared), &
!$OMP&shared(NTH,t_DscfL,t_delkmats), &
!$OMP&private(TID,last,delkmats,irealim), &
!$OMP&private(ip,nc,nlocal,ic,imp,i,il,iu,iul,ii,ind,j,ijl,jl), &
!$OMP&private(lasta,lastop,ia,is,iop,isp,ix,dxsp,r2sp,nphiloc,iphi,jc), &
!$OMP&private(Vij,VClocal,DscfL,ilocal,ilc,iorb,Vlocal,Clocal,phia), &
!$OMP&private(irel,inmp,dxp,dxpgrid,dist,kxij,iua,displaat,aux)

!$OMP single
#ifdef _OPENMP
    NTH = omp_get_num_threads( )
#else
    NTH = 1
#endif
!$OMP end single ! implicit barrier, IMPORTANT

#ifdef _OPENMP
    TID = omp_get_thread_num( ) + 1
#else
    TID = 1
#endif

    nullify(Clocal,phia,ilocal,ilc,iorb,Vlocal)
!$OMP critical
! Perhaps the critical section is not needed,
! however it "tells" the OS to allocate per
! thread, possibly waiting for each thread to
! place the memory in the best position.
    allocate( Clocal(nsp,maxloc2) )
    allocate( ilocal(no) , ilc(maxloc2) , iorb(maxloc) )
    allocate( Vlocal(triang,2) )
    if ( DirectPhi ) allocate( phia(maxoa,nsp) )
!$OMP end critical

!$OMP single
    if ( ParallelLocal ) then
       nullify( t_DscfL )
       call re_alloc( t_DscfL, 1, nvmaxl, 1, 2, 1, NTH, &
            'DscfL',  'delk' )
    else
       if ( NTH > 1 ) then
          nullify( t_delkmats )
          call re_alloc( t_delkmats, 1, nvmax, 1, NTH, &
               'delkmats',  'delk' )
       end if
    end if
!$OMP end single

    if ( ParallelLocal ) then
       DscfL => t_DscfL(1:nvmaxl,1:2,TID)
       DscfL(1:nvmaxl,1:2) = 0.0_dp
    else
       if ( NTH > 1 ) then
          delkmats => t_delkmats(1:nvmax,TID)
       else
          delkmats => delkmat
       end if
    end if

!   Full initializations done only once
    ilocal(1:no)         = 0
    iorb(1:maxloc)       = 0
    Vlocal(1:triang,1:2) = 0.0_dp
    last                 = 0

!   Loop over grid points
!$OMP do
    do ip = 1,np
!      Find number of nonzero orbitals at this point
       nc = endpht(ip) - endpht(ip-1)
!      Find new required size of Vlocal
       nlocal = last
       do ic = 1,nc
          imp = endpht(ip-1) + ic
          i = lstpht(imp)
          if (ilocal(i) .eq. 0) nlocal = nlocal + 1
       end do

!      If overflooded, add Vlocal to delkmat and reinitialize it
       if (nlocal > maxloc .and. last > 0) then
          if ( ParallelLocal ) then
             do il = 1,last
                i = iorb(il)
                iu = indxuo(i)
                iul = NeedDscfL(iu)
                if ( i == iu ) then
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul)+ii
                      j = listdl(ind)
                      jl = ilocal(j)
                      ijl = idx_ijl(il,jl)
! The variables we want to compute in this subroutine are complex numbers
! Here, when irealim =1 we refer to the real part, and
! when irealim = 2 we refer to the imaginary part
                      DscfL(ind,:) = DscfL(ind,:) + Vlocal(ijl,:) * dVol
                   end do
                else
                   ia  = iaorb(i)
                   iua = indxua(ia)
                   do ix = 1, 3
                      displaat(ix) = &
                           (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                           (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                           (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                   end do
                   dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                   kxij    = kpoint(1) *  dist(1)  + &
                        kpoint(2) *  dist(2)  + &
                        kpoint(3) *  dist(3)
                   do ii = 1, numdl(iul)
                      ind = listdlptr(iul)+ii
                      j   = listsc( i, iu, listdl(ind) )
                      jl  = ilocal(j)
                      ijl = idx_ijl(il,jl)
                      aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                           Vlocal(ijl,2) * dsin(kxij)
                      aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                           Vlocal(ijl,1) * dsin(kxij)
                      DscfL(ind,:) = DscfL(ind,:) + aux(:) * dVol
                   end do
                end if
             end do
          else 
             
             do il = 1,last
                i = iorb(il)
                iu = indxuo(i)
                call GlobalToLocalOrb( iu, Node, Nodes, iul )
                if (i == iu) then
                   do ii = 1, numVs(iul)
                      ind = listVsptr(iul)+ii
                      j = listVs(ind)
                      jl = ilocal(j)
                      ijl = idx_ijl(il,jl)
                      delkmats(ind) = delkmats(ind) + &
                           cmplx(Vlocal(ijl,1), Vlocal(ijl,2), kind=dp) * dVol

                   end do
                else
                   ia  = iaorb(i)
                   iua = indxua(ia)
                   do ix = 1, 3
                      displaat(ix) = &
                           (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                           (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                           (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
                   end do
                   dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
                   kxij    = kpoint(1) *  dist(1)  + &
                        kpoint(2) *  dist(2)  + &
                        kpoint(3) *  dist(3)
                   do ii = 1, numVs(iul)
                      ind = listVsptr(iul)+ii
                      j = listsc( i, iu, listVs(ind) )
                      jl = ilocal(j)
                      ijl = idx_ijl(il,jl)
                      aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                           Vlocal(ijl,2) * dsin(kxij)
                      aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                           Vlocal(ijl,1) * dsin(kxij)

                      delkmats(ind) = delkmats(ind) + &
                           cmplx(aux(1),aux(2),kind=dp) * dVol

                   end do
                end if
             end do

          end if

!         Reset local arrays
          do ii = 1, last
             ilocal(iorb(ii)) = 0
          end do
          iorb(1:last) = 0
          ijl = (last+1)*(last+2)/2
          Vlocal(1:ijl,1:2) = 0.0_dp
          last = 0
       end if

!      Look for required orbitals not yet in Vlocal
       if ( nlocal > last ) then
          do ic = 1,nc
             imp = endpht(ip-1) + ic
             i = lstpht(imp)
             if (ilocal(i) == 0) then
                last = last + 1
                ilocal(i) = last
                iorb(last) = i
             end if
          end do
       end if

       if ( DirectPhi ) then
          lasta = 0
          lastop = 0
          do ic = 1 , nc
             imp     = endpht(ip-1) + ic
             i       = lstpht(imp)
             il      = ilocal(i)
             ia      = iaorb(i)
             iop     = listp2(imp)
             ilc(ic) = il
             
!            Localize the position of the mesh point
             irel = idop(iop) + ipa(ia)
             call ipack( -1, 3, nem, inmp, irel )
             inmp(:) = inmp(:) + (meshLim(1,:) - 1)
             inmp(:) = inmp(:) - 2 * ne(:)
             
             dxp(:) = cmesh(:,1) * inmp(1) + &
                  cmesh(:,2) * inmp(2) + &
                  cmesh(:,3) * inmp(3)
             
             do isp = 1, nsp
                dxpgrid(:,isp) = dxp(:) + xdsp(:,isp)
             end do
             
!            Generate phi values
             if ( ia /= lasta .or. iop /= lastop ) then
                lasta = ia
                lastop = iop
                is = isa(ia)
                do isp = 1,nsp
                   dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                   r2sp = sum(dxsp**2)
                   if ( r2sp < r2cut(is) ) then
!$OMP critical
                      call all_phi( is, +1, dxsp, maxoa, nphiloc, phia(:,isp) )
!$OMP end critical
                   else
                      phia(:,isp) = 0.0_dp
                   end if
                end do
             end if
             iphi = iphorb(i)
             Clocal(:,ic) = phia(iphi,:)

!            Pre-multiply V and Clocal(,ic)
             Vij(:) = 0._dp
             do isp = 1 , nsp
                kxij = kpoint(1) * dxpgrid(1,isp) + &
                     kpoint(2) * dxpgrid(2,isp) + &
                     kpoint(3) * dxpgrid(3,isp)
                VClocal(1,isp) = dcos(kxij) * Clocal(isp,ic)
                VClocal(2,isp) = dsin(kxij) * Clocal(isp,ic)
                Vij(:) = Vij(:) + VClocal(:,isp) * Clocal(isp,ic)
             end do
             
!            ic == jc, hence we simply do
             ijl = idx_ijl(il,ilc(ic))
             Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:)
             
!            Loop on second orbital of mesh point
             do jc = 1 , ic - 1
                jl = ilc(jc)
                
!               Loop over sub-points
                Vij(:) = 0.0_dp
                do isp = 1 , nsp
                   Vij(:) = Vij(:) + VClocal(:,isp) * Clocal(isp,jc)
                end do
                
                ijl = idx_ijl(il,jl)
                if ( il == jl ) then
                   Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:) * 2._dp
                else
                   Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:)
                end if
                
             end do
             
          end do
          
       else

!         Loop on first orbital of mesh point
          do ic = 1 , nc
             imp     = endpht(ip-1) + ic
             i       = lstpht(imp)
             il      = ilocal(i)
             ia      = iaorb(i)
             iop     = listp2(imp)
             ilc(ic) = il

!            Localize the position of the mesh point
             irel = idop(iop) + ipa(ia)
             call ipack( -1, 3, nem, inmp, irel )
             inmp(:) = inmp(:) + (meshLim(1,:) - 1)
             inmp(:) = inmp(:) - 2 * ne(:)

             dxp(:) = cmesh(:,1) * inmp(1) + &
                  cmesh(:,2) * inmp(2) + &
                  cmesh(:,3) * inmp(3)

             do isp = 1, nsp
                dxpgrid(:,isp) = dxp(:) + xdsp(:,isp)
             end do

!            Retrieve phi values
             Clocal(:,ic) = phi(:,imp)

!            Pre-multiply V and Clocal(,ic)
             Vij(:) = 0._dp
             do isp = 1 , nsp
                kxij = kpoint(1) * dxpgrid(1,isp) + &
                     kpoint(2) * dxpgrid(2,isp) + &
                     kpoint(3) * dxpgrid(3,isp)
                VClocal(1,isp) = dcos(kxij) * Clocal(isp,ic)
                VClocal(2,isp) = dsin(kxij) * Clocal(isp,ic)
                Vij(:) = Vij(:) + VClocal(:,isp) * Clocal(isp,ic)
             end do

!            ic == jc, hence we simply do
             ijl = idx_ijl(il,ilc(ic))
             Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:)

!            Loop on second orbital of mesh point
             do jc = 1 , ic - 1
                jl = ilc(jc)

!               Loop over sub-points
                Vij(:) = 0.0_dp
                do isp = 1 , nsp
                   Vij(:) = Vij(:) + VClocal(:,isp) * Clocal(isp,jc)
                end do

                ijl = idx_ijl(il,jl)
                if ( il == jl ) then
                   Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:) * 2._dp
                else
                   Vlocal(ijl,:) = Vlocal(ijl,:) + Vij(:)
                end if

             end do

          end do
       end if

    end do
!$OMP end do nowait

!$OMP barrier
    
!   Add final Vlocal to delkmat
    if ( ParallelLocal .and. last > 0 ) then

       do il = 1 , last
          i = iorb(il)
          iu = indxuo(i)
          iul = NeedDscfL(iu)
          if ( i == iu ) then
             do ii = 1, numdl(iul)
                ind = listdlptr(iul)+ii
                j = listdl(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                DscfL(ind,:) = DscfL(ind,:) + Vlocal(ijl,:) * dVol
             end do
          else
             ia  = iaorb(i)
             iua = indxua(ia)
             do ix = 1, 3
                displaat(ix) = &
                     (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                     (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                     (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
             end do
             dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
             kxij    = kpoint(1) * dist(1) + &
                  kpoint(2) * dist(2) + &
                  kpoint(3) * dist(3)
             do ii = 1, numdl(iul)
                ind = listdlptr(iul)+ii
                j = listsc( i, iu, listdl(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                     Vlocal(ijl,2) * dsin(kxij)
                aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                     Vlocal(ijl,1) * dsin(kxij)
                DscfL(ind,:) = DscfL(ind,:) + aux(:) * dVol
             end do
          end if
       end do

    else if ( last > 0 ) then

       do il = 1 , last
          i = iorb(il)
          iu = indxuo(i)
          if ( i == iu ) then
             do ii = 1, numVs(iu)
                ind = listVsptr(iu)+ii
                j = listVs(ind)
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                delkmats(ind) = delkmats(ind) + &
                     cmplx(Vlocal(ijl,1), Vlocal(ijl,2),kind=dp) * dVol
             end do
          else
             ia  = iaorb(i)
             iua = indxua(ia)
             do ix = 1, 3
                displaat(ix) = &
                     (iatfold(1,ia)*nmsc(1))*cmesh(ix,1)+ &
                     (iatfold(2,ia)*nmsc(2))*cmesh(ix,2)+ &
                     (iatfold(3,ia)*nmsc(3))*cmesh(ix,3)
             end do
             dist(:) = xa(:,iua) - xa(:,ia) - displaat(:)
             kxij    = kpoint(1) * dist(1) + &
                  kpoint(2) * dist(2) + &
                  kpoint(3) * dist(3)
             do ii = 1, numVs(iu)
                ind = listVsptr(iu)+ii
                j = listsc( i, iu, listVs(ind) )
                jl = ilocal(j)
                ijl = idx_ijl(il,jl)
                aux(1) = Vlocal(ijl,1) * dcos(kxij) - &
                     Vlocal(ijl,2) * dsin(kxij)
                aux(2) = Vlocal(ijl,2) * dcos(kxij) + &
                     Vlocal(ijl,1) * dsin(kxij)
                delkmats(ind) = delkmats(ind) + &
                     cmplx(aux(1),aux(2),kind=dp) * dVol
             end do
          end if
       end do

    end if

!$OMP barrier

    if ( ParallelLocal .and. NTH > 1 ) then
       do irealim = 1, 2
!$OMP do
          do ind = 1, nvmaxl
             do ii = 2, NTH
                t_DscfL(ind,irealim,1) = t_DscfL(ind,irealim,1) + &
                     t_DscfL(ind,irealim,ii)
             end do
          end do
!$OMP end do nowait
       end do
!$OMP barrier
    else if ( NTH > 1 ) then
!$OMP do
       do ind = 1, nvmax
          do ii = 1, NTH
             delkmat(ind) = delkmat(ind) + t_delkmats(ind,ii)
          end do
       end do
!$OMP end do
    end if

!   Free local memory
    deallocate(Clocal,ilocal,ilc,iorb,Vlocal)
    if ( DirectPhi ) deallocate( phia )

!$OMP master
    if ( ParallelLocal ) then
!      Redistribute Hamiltonian from mesh to orbital based distribution
       DscfL => t_DscfL(1:nvmaxl,1:2,1)
       call matrixMtoOC( nvmaxl, nvmax, numVs, listVsptr, nuo, DscfL, delkmat )
       call de_alloc( t_DscfL, 'DscfL', 'delk' )
    else if ( NTH > 1 ) then
       call de_alloc( t_delkmats, 'delkmats', 'delk' )
    end if
!$OMP end master

!$OMP end parallel

    call de_alloc( r2cut, 'r2cut', 'delk' )

#ifdef _TRACE_
    call MPI_Barrier( MPI_Comm_World, MPIerror )
    call MPItrace_event( 1000, 0 )
#endif

    call timer('delk',2)

#ifdef DEBUG
    call write_debug( '    POS delk' )
#endif

  contains

!   In any case will the compiler most likely inline this
!   small routine. So it should not pose any problem.
    pure function idx_ijl(i,j) result(ij)
      integer, intent(in) :: i,j
      integer :: ij
      if ( i > j ) then
         ij = i * (i + 1)/2 + j + 1
      else
         ij = j * (j + 1)/2 + i + 1
      end if
    end function idx_ijl

  end subroutine delk

end module m_delk
