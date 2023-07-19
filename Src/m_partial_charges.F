      module m_partial_charges
      public :: compute_partial_charges
      logical, public :: want_partial_charges = .false.
      private
      CONTAINS
!
      subroutine compute_partial_charges(DRho,rhoatm,
     .                  nspin, iaorb, iphorb, 
     .                  isa, nmpl,dvol)


C Calculates the Hirshfeld and Voronoi charges 
C Coded by P. Ordejon, August 2004
C Integrated into dhscf workflow by A. Garcia, October 2012
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nspin         : Number of different spin polarisations
C                         nspin=1 => Unpolarized, nspin=2 => polarized
C                         nspin=4 => Noncollinear spin -NOT IMPLEMENTED-
C integer iaorb(no_s)   : Atom to which each orbital belongs
C integer iphorb(no_s)  : Orbital index (within atom) of each orbital
C integer indxuo        : Index of equivalent orbital in unit cell
C integer isa(na_s)     : Species index of all atoms in supercell
C real(grid_p) Drho(:,:): Electron charge density (in cluster form)
C real(grid_p) rhoatm(:): Superposition of atomic Electron charge densities
C ----------------------------------------------------------------------
C Output : None
C ----------------------------------------------------------------------

C  Modules
C
      use precision,    only:  dp, grid_p
      use parallel
      use atomlist, only: indxua, indxuo, no_s, no_u, Datm
      use siesta_geom, only: na_u, na_s
      use atmfuncs, only: rcut, rcore, symfio, cnfigfio, labelfis, 
     .                    nofis, zvalfis
      use sys
!      use mesh
      use siesta_options, only: hirshpop, voropop
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in)     :: nspin,
     .  iaorb(no_s), iphorb(no_s), isa(na_s), nmpl

      real(dp), intent(In)         :: dvol
      real(grid_p), dimension(:), intent(in)   :: rhoatm
      real(grid_p), dimension(:,:), intent(in) :: DRho

C
C Internal variables
C  real*8 q(na) : atomic charges 
C
      integer ::   ia, ispin, ns, ispec, is

#ifdef MPI
      integer
     .  MPIerror
#endif

      real(dp), dimension(:), allocatable       :: q, qtmp

      character(len=20) ::  atm_label

      external  memory

C ----------------------------------------------------------------------
C General initialisation
C ----------------------------------------------------------------------

      if (nspin .eq. 4) then
        call message("Hirshfeld and Voronoi populations" //
     $               " not implemented for non-collinear spin")
        return
      endif

C ----------------------------------------------------------------------
      allocate(q(na_u))
      call memory('A','D',na_u,'hirsh')
      allocate(qtmp(na_u))
      call memory('A','D',na_u,'hirsh')

C ----------------------------------------------------------------------
C Find Hirshfeld charges
C ----------------------------------------------------------------------

      if (hirshpop) then
        call hirshfeld( no_s, na_s, na_u, nspin, indxuo, indxua, 
     .                  nmpl, datm, rhoatm, DRho, iaorb, iphorb, isa, 
     .                  q )

#ifdef MPI
        call MPI_AllReduce(q,qtmp,na_u,
     .    MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
        q(1:na_u) = qtmp(1:na_u)
#endif

        if (Node.eq.0) then
           write(6,"(/,a)")  'Hirshfeld Net Atomic Populations:'
           write(6,'(a6,2x,a7,2x,a)') 'Atom #', 'Qatom', 'Species'
           do ia = 1,na_u
              is = isa(ia)
              atm_label = labelfis(is)
              write(6,'(i6,2x,f7.3,2x,a20)'), ia,
     $             zvalfis(is)-q(ia)*dvol, atm_label
           enddo
        endif
      endif

C ----------------------------------------------------------------------
C Find Voronoi charges
C ----------------------------------------------------------------------

      if (voropop) then
        call voronoi( no_s, na_s, na_u, nspin, indxua, 
     .                nmpl, rhoatm, DRho, iaorb,
     .                q )
#ifdef MPI
        call MPI_AllReduce(q,qtmp,na_u,
     .    MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
        q(1:na_u) = qtmp(1:na_u)
#endif

        if (Node.eq.0) then
           write(6,"(/,a)")  'Voronoi Net Atomic Populations:'
           write(6,'(a6,2x,a7,2x,a)') 'Atom #', 'Qatom', 'Species'
           do ia = 1,na_u
              is = isa(ia)
              atm_label = labelfis(is)
              write(6,'(i6,2x,f7.3,2x,a20)'), ia,
     $             -q(ia)*dvol, atm_label
           enddo
        endif
      endif

      call memory('D','D',size(qtmp),'hirsh')
      deallocate(qtmp)
      call memory('D','D',size(q),'hirsh')
      deallocate(q)

      end subroutine compute_partial_charges


      subroutine hirshfeld( no_s, na_s, na_u, nspin, indxuo, indxua, np,
     .                   Datm, rhoatm, rhoscf, iaorb, iphorb, isa, qh )
C ********************************************************************
C Finds the Hirshfeld atomic charges  
C Hirshfeld, Theo Chem Acta 44, 129 (1977)
C See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
C 
C Written by P.Ordejon. August'04.
C *********************** InpUT **************************************
C integer no_s                : Number of basis orbitals
C integer na_s                : Number of atoms
C integer na_u               : Number of atoms in unit cell
C integer nspin             : Number of spins
C integer indxuo(no_s)        : Index of equivalent orbital in unit cell
C integer indxua(na_s)        : Index of equivalent atom in unit cell
C integer np                : Number of mesh points
C real*8 Datm(no_s)           : Occupations of basis orbitals in free atom
C integer iaorb(*)          : Pointer to atom to which orbital belongs
C integer iphorb(*)         : Orbital index within each atom
C integer isa(*)            : Species index of all atoms
C real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
C real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
C *********************** OUTPUT **************************************
C real*8 qh(na_u)            : Hirshfeld charge on each orbital
C *********************************************************************
C
C  Modules
C
      use precision, only: dp, grid_p
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshphi
C
      implicit none

      integer          no_s, np, nspin, na_s, na_u
      integer          indxuo(no_s), indxua(na_s),
     .                 iaorb(*), iphorb(*), isa(*)

      real(grid_p) ::  rhoatm(nsp,np),rhoscf(nsp,np,nspin)

      real(dp)         Datm(no_s), phip, qh(na_u)

      integer          i, ip, isp, iu, kn, io, iop, is, iphi, ia, ix, 
     .                 ispin, iua
      real(dp)         Ci, gradCi(3), r2o, r2sp, dxsp(3), Qi

      ! A very small number to avoid division by zero
      real(grid_p)  :: rhoatm_tolerance = 1.0e-12_grid_p


C  Initialise Hirshfeld charges

      do iua = 1, na_u
        qh(iua) = 0.0_dp
      enddo

C  Loop on mesh points
      do ip = 1,np

C  Loop on orbitals of mesh point
        do kn = 1+endpht(ip-1), endpht(ip)
          i = lstpht(kn)
          iu = indxuo(i)

C  Generate phi value and loop on subpoints
          iphi = iphorb(i)
!!          if (iphi .gt. no_s) call die("hirshfeld: iphi error")
          ia = iaorb(i)
          iua = indxua(ia)
          is = isa(ia)
          r2o = rcut(is,iphi)**2
          iop = listp2(kn)
          do isp = 1,nsp
             if (abs(rhoatm(isp,ip)) <= rhoatm_tolerance) then
                ! Atomic charge is really too small here... so
                ! we do not count this point
                CYCLE
c$$$                 write(6,"(a,g14.6,/,10x,a,2g14.6)")
c$$$     $               "rhoatm: ", rhoatm(isp,ip), "rho: ",
c$$$     $                rhoscf(isp,ip,:)
              endif

            do ix = 1,3
              dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
            enddo
            r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
            if (r2sp.lt.r2o) then
              call phiatm(is,iphi,dxsp,phip,gradCi)
              Ci = phip
              Qi = Datm(iu) * Ci * Ci
              do ispin = 1, nspin
                qh(iua) = qh(iua) 
     .                + Qi * rhoscf(isp,ip,ispin)/rhoatm(isp,ip)
              enddo
            endif
          enddo

        enddo

      enddo

      end subroutine hirshfeld


      subroutine voronoi( no_s, na_s, na_u, nspin, indxua, np,
     .                   rhoatm, rhoscf, iaorb, qv )
C ********************************************************************
C Finds the Voronoi atomic charges
C Bickelhaupt et al, Organometallics 15, 2923 (1996)
C See Fonseca et al, J. Comp. Chem. 25, 189 (2003)
C 
C Written by P.Ordejon. August'04.
C *********************** InpUT **************************************
C integer no_s                : Number of basis orbitals
C integer na_s                : Number of atoms
C integer na_u               : Number of atoms in unit cell
C integer nspin             : Number of spins
C integer indxua(na_s)        : Index of equivalent atom in unit cell
C integer np                : Number of mesh points
C integer iaorb(*)          : Pointer to atom to which orbital belongs
C real rhoatm(nsp,np)       : Harris (sum of atoms) density at mesh points
C real rhoscf(nsp,np,nspin) : Selfconsistent charge density at mesh points
C *********************** OUTPUT **************************************
C real*8 qv(na_u)              : Voronoi charge on each orbital
C *********************************************************************
C
C  Modules
C
      use precision, only: dp, grid_p
      use atmfuncs, only: rcut, phiatm
      use mesh,     only: nsp, dxa, xdop, xdsp
      use meshphi
C
      implicit none

      integer          no_s, np, nspin, na_s, na_u
      integer          indxua(na_s),
     .                 iaorb(no_s)

      real(grid_p) ::  rhoatm(nsp,np),rhoscf(nsp,np,nspin)

      real(dp)     ::  phip, qv(na_u)

      integer          i, ip, isp, iu, kn, io, iop, is, ia, ix, 
     .                 ispin, iua, ntimes(na_u), maxeq, neq, ieq

      real(dp)         Ci, gradCi(3), r2o, r2sp, dxsp(3), Qi, huge, 
     .                 rmin, qtot, tol

      logical          eq

      parameter (huge = 1.0e30_dp)
      parameter (maxeq = 20)
      parameter (tol = 1.0e-2_dp)

      integer          iatom(maxeq)


C  Initialise Hirshfeld charges

      do iua = 1, na_u
        ntimes(iua) = 0
        qv(iua) = 0.0_dp
      enddo


C  Loop on mesh points
      do ip = 1,np

C  Loop on mesh subpoints
        do isp = 1,nsp
          rmin = huge
          eq = .false.
          neq = 1
          do ieq = 1, maxeq
            iatom(ieq) = 0
          enddo

C  Loop on orbitals of mesh point, to check which atom is closest to subpoint
          do kn = 1+endpht(ip-1), endpht(ip)
            i = lstpht(kn)
            ia = iaorb(i)

C  iua is the index of the atom in unit cell
            iua = indxua(ia)

            iop = listp2(kn)
            do ix = 1,3
              dxsp(ix) = xdop(ix,iop) + xdsp(ix,isp) - dxa(ix,ia)
            enddo
            r2sp = dxsp(1)**2 + dxsp(2)**2 + dxsp(3)**2
            r2sp = dsqrt(r2sp)

C  If distance is equal to the previous minimum, determine if it is
C  another atom, or the same. If it is a different one, then
C  handle multiplicity of nearest atoms to grip point
C  Consider distances equal if within tolerance

            if (dabs(rmin-r2sp) .lt. tol) then
              eq = .false.
              do ieq = 1,neq
                if (iua .eq. iatom(ieq)) eq = .true.
              enddo
              if (.not. eq) then
                neq = neq+1
                if (neq .gt. maxeq) stop 'voronoi: increase maxeq'
                iatom(neq) = iua
              endif
              goto 100
            endif

            if (r2sp .lt. rmin) then
              neq = 1
              iatom(neq) = iua
              rmin = r2sp
            endif

100         continue

          enddo

C  Assign charge to atom iatom; if no atom was found, then the charge
C  is zero, so move to next subpoint

          if (iatom(1) .ne. 0) then
            do ieq = 1,neq
              ntimes(iatom(ieq)) = ntimes(iatom(ieq))+1
            enddo
            do ispin = 1, nspin
              do ieq = 1,neq
                qv(iatom(ieq)) = qv(iatom(ieq)) 
     .            + (rhoscf(isp,ip,ispin) - (rhoatm(isp,ip)/nspin))/neq
              enddo
            enddo
          else
C  Check that charge is actually zero
            qtot=0.0_dp
            do ispin=1,nspin
              qtot = qtot + rhoscf(isp,ip,ispin) + rhoatm(isp,ip)
            enddo
            if (qtot .gt. 0.0d0) stop 'voronoi: Error in grid charge'
          endif

        enddo

      enddo

      end subroutine voronoi
      end module m_partial_charges
