!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine moments(iopt, natoms, no_u, maxnh, numh,
     .                    listhptr,listh,s,dm,isa,lasto,iaorb,
     .                    iphorb, indxuo)
C ********************************************************************
C Subroutine to calculate orbital moments on atoms.
C The density matrix (d.m.) and overlap matrix are passed in sparse form
C (both with the same sparse structure)
C There is no output. The populations are printed to the output.
C
C A modification of mulliken.F
C
C ************************** INPUT ************************************
C integer iopt                : Work option: 1 = atomic and orbital charges
C integer natoms              : Number of atoms in unit cell
C integer no_u           : Number of basis orbitals over all nodes
C integer maxnh               : First dimension of d.m. and overlap, and its
C                               maximum number of non-zero elements
C integer numh(no_l)        : First Control vector of d.m. and overlap
C integer listhptr(no_l)    : Second Control vector of d.m. and overlap
C integer listh(maxnh)        : Third Control vector of d.m. and overlap
C real*8  s(maxnh)            : Overlap matrix in sparse form
C real*8  dm(maxnh,8): Density matrix in sparse form 
C integer isa(natoms)         : Species index of each atom
C integer lasto(0:natoms)     : Index of last orbital of each atom
C                               (lasto(0) = 0) 
C integer iaorb(no_u)       : Atomic index of each orbital
C integer iphorb(no_u)      : Orbital index of each orbital in its atom
C integer indxuo(*)        : Index of equivalent unit-cell orbital
C ************************* OUTPUT *************************************
C No output. The results are printed to standard output
! Written by Jaime Ferrer (February 2009).
! Retain only the orbital angular moments (Alberto Garcia, April 2018)
C **********************************************************************
C
!
! JFR notes
!
! The density matrix is written as
!      | n(1)         n(3) - i n(4) |       | n   +   m_z   m_x - i m_y |
!  n = |                            | = 1/2 |                           |
!      | n(3) + i n(4)     n(2)     |       | m_x + i m_y   n   -   m_z |
!
! Therefore vec{m}(r) = (mx, my, mz) = Tr [ n vec{tau} ] 
!                     = ( 2 n(3), 2 n(4), n(1) - n(2) )
! where tau are the three Pauli spin matrices.
!
! The expectation value of the spin is 
! <vec{S}> = 1/2 int dr vec{m}(r) 
!          = sum S_{mn} ( rho_{mn}(3), rho_{mn}(4), 1/2 (rho_{mn}(1)-rho_{mn}(2)) )
! where S is the overlap matrix and rho is SIESTA's density matrix.
!
! The Zeeman interaction term in the hamiltonian is
! E_{zeeman} = mu_B (g <vec{S}> + <vec{L}> ) vec{B}
!
! As a consequence, the magnetic moments are
! vec{M} = - partial E_{zeeman} / partial vec{B} = 
!        = vec{M}_S + vec{M}_L = - mu_B ( g <vec{S}> + <vec{L}>
!
! This subroutine writes the magnetic moments. The expectation value
! of the spin is obtained by dividing by 2.
!
! The expectation values of the orbital angular moments of all the electrons in the system,
!  <vec{L}_a>, with respect to an origin centered at atom "a" are 
! <L_a> = Sum_n f_n <psi_n | L_a | psi_n >
!       = sum_{mu,nu} <phi_mu | L_a | phi_nu> (Dscf_{nu,mu}^{up,up}+Dscf_{nu,mu}^{dn,dn})
!  
! Since <phi_mu | L_a |phi_nu > are purely imaginary, only the imaginary components
! of the density matrix are selected.
!
! Instead of <L_a>, this subroutine computes the "atomic" orbital moments, where
! <L> is computed for each single atom "b", using an an on-site approximation, 
! where only the orbitals at atom "b" are considered. Further, in this version,
! only orbitals with same z-number are taken into account.
! In a way, the approximation consists of assuming that each electron contributes
! to only a few atoms, and that all <L_b> can in the end be summed up even though
! they use different origins.
!
! Notice that within DFT <S^2> = 3 hbar/ 4 N, where N is the number of electrons.
! This rather poor estimate of <S^2> comes about because of the use of Slater determinats.
! As an alternative, this subroutine computes (<S_x>_a,<S_y>_a,<S_z>_a) for each atom "a", 
! by summing all Mulliken populations for atom "a", and then uses the estimate 
! <S_a^2> ~ <S_x>_a^2 + <S_y>_a^2 + <S_z>_a^2 
! Also, the total spin of the system is computed by summing the individual spins:
! <vec S> = (<S_x>,<S_y>,<S_z>) = sum_a (<S_x>_a,<S_y>_a,<S_z>_a) 
! and <S^2> ~ <S_x>^2 + <S_y>^2 + <S_z>^2
!
! Notice that <S^2> =/ sum_a <S_a^2> because of crossing terms.
!
C  Modules
C
      use precision, only: dp
      use parallel,  only: IOnode, Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, GetNodeOrbs,
     $                         LocalToGlobalOrb
      use atmfuncs, only: symfio, cnfigfio, labelfis, nofis,
     $                    mofio, lofio, zetafio
      use spinorbit, only: int_so_ang
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in) ::  iopt, natoms, no_u,maxnh

      integer, intent(in) ::
     .  numh(*),lasto(0:natoms),listh(maxnh),listhptr(*),
     .  iphorb(*), isa(natoms), iaorb(*), indxuo(*)

      real(dp), intent(in) :: DM(maxnh,8), S(maxnh)

      external memory

      integer :: no_l, io_u, io_l, li, mi, zi

#ifdef MPI
      integer MPIerror, mpistatus(MPI_STATUS_SIZE)
      real(dp), dimension(3) :: qb, pb
      character(len=200)     :: mpibuff
#endif
      integer i, ioa, ia, in, io, is, ind, ns, ispec, config, 
     .        jo, jo_u, joa, ja, lj, mj, zj, iNode

      real(dp) :: L(3)
      real(dp) l_atm(3), l_sqr, lvec(3), l_tot(3)
      real(dp), dimension(:,:), allocatable ::  l_orb


      character sym_label*7, atm_label*20

C ......................

#ifdef MPI
C Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)
#else
      no_l = no_u
#endif

      if (iopt.eq.0) then
C iopt = 0 implies no analysis
        return
      elseif (iopt.lt.0 .or. iopt.gt.1) then
        if (ionode) then
          write(6,"(a)") 'moments: ERROR: Wrong flag'
        endif
        return
      endif

C Allocate local memory
      allocate(l_orb(3,no_l))
      call memory('A','D',no_l*3,'moments')

      ns=0
      do i = 1,natoms
        ns=max(ns,isa(i))
      enddo 

C Compute Orbital and Spin moments by orbitals.........
      if (ionode) then
        write(6,*) 
        write(6,"(a)")
     .   'moments: Magnetic moments from orbital angular momenta:'
      endif
      l_orb(:,:) = 0.d0

C --- accummulate orbital and spin moments in each orbital:

      do io = 1,no_l
        call LocalToGlobalOrb(io,Node,Nodes,io_u)
        ia = iaorb(io_u)
        ioa = iphorb(io_u)
        is = isa(ia)
        li = lofio (is,ioa)
        if (li == 0) CYCLE  ! No contribution from l=0 orbs
        mi = mofio (is,ioa)
        zi = zetafio(is,ioa)
        do in = 1,numh(io)
          ind = listhptr(io)+in
          jo = listh(ind)
          jo_u = indxuo(jo)
          if (jo .ne. jo_u) CYCLE ! Not in the unit cell
          ja = iaorb(jo)
          if (ja .ne. ia) CYCLE ! Not in the same atom
          joa = iphorb(jo)
          lj=lofio(is,joa)
          if (li /= lj) CYCLE ! Different l
          zj = zetafio(is,joa)
          if (zi /= zj) CYCLE ! Different z
          mj=mofio(is,joa)
          if (mi == mj) CYCLE ! Same m
          call int_so_ang(li, mj, mi, L(:))
          l_orb(1,io)=l_orb(1,io)+L(3)*(dm(ind,5)+dm(ind,6))
          l_orb(2,io)=l_orb(2,io)-L(2)*(dm(ind,5)+dm(ind,6))
          l_orb(3,io)=l_orb(3,io)-L(1)*(dm(ind,5)+dm(ind,6))
        enddo
      enddo

C .... printout in a loop over species, and atoms in species

      l_tot(:) = 0.0_dp

      do ispec=1,ns

        atm_label=labelfis(ispec)

        if (ionode) then
         write(6,'(/2a)')'Species: ', atm_label
         write(6,'(/,a4,a7,8x,a12,a18,/,100("-"))')
     .            'Atom', 'Orb', 'sqrt(<L>^2)', '<(Lx,Ly,Lz)>'
        endif

        do ia = 1,natoms 
         if (isa(ia) /= ispec) CYCLE  
         do iNode = 0, Nodes-1
          if (iNode .ne. Node) CYCLE
          l_atm(:) = 0.0_dp  
           do io = lasto(ia-1)+1,lasto(ia)  
            call GlobalToLocalOrb(io,Node,Nodes,io_l)
            if (io_l .gt. 0) then 
              l_atm(1:3) = l_atm(1:3) + l_orb(1:3,io_l)
              l_tot(1:3) = l_tot(1:3) + l_orb(1:3,io_l)
              sym_label=symfio(ispec,iphorb(io))
              config=cnfigfio(ispec,iphorb(io))
              lvec(1) = l_orb(1,io_l)
              lvec(2) = l_orb(2,io_l)
              lvec(3) = l_orb(3,io_l)
              l_sqr = sqrt(lvec(1)**2+lvec(2)**2+lvec(3)**2)
              if (ionode) then
                write(6,
     .            '(i4,i5,2x,i1,a6,1x,f8.3,2x,3f7.3)')
     .             ia, io, config, sym_label, l_sqr, lvec(1:3)
#ifdef MPI
              else
                write(mpibuff,
     .            '(i4,i5,2x,i1,a6,1x,f8.3,2x,3f7.3)')
     .             ia, io, config, sym_label, l_sqr, lvec(1:3)
                call mpi_send(mpibuff, 200, MPI_Character,
     .                        0, io, MPI_Comm_world, mpierror)
              endif
            else
              if (ionode) then
                call mpi_recv(mpibuff, 200, MPI_Character,
     .                          MPI_ANY_SOURCE, io,
     .                          MPI_Comm_world, mpistatus, mpierror)
                          write(6,'(a)')  trim(mpibuff)
#endif
              endif
            endif
           enddo ! io
           lvec(1:3) = l_atm(1:3)
         enddo ! iNodes

#ifdef MPI
C Global reduction of terms
         qb(1:3)=lvec(1:3)
         call MPI_Reduce(qb,pb,3,MPI_double_precision,MPI_sum,0,
     .                   MPI_Comm_World,MPIerror)
         lvec(1:3)=pb(1:3)
#endif
         if (ionode) then
           l_sqr = sqrt(lvec(1)**2+lvec(2)**2+lvec(3)**2)
           write(6,'(100("-"))')
           write(6,'(i4,5x,a6,4x,f8.3,2x,3f7.3/)')
     .           ia, 'Total', l_sqr, lvec
         endif
        enddo    !   do ia = 1,natoms 
      enddo    !   do ispec=1,ns

      lvec(1:3) = l_tot(1:3)

#ifdef MPI
C Global reduction of terms
      qb(1:3)=lvec(1:3)
      call MPI_Reduce(qb,pb,3,MPI_double_precision,MPI_sum,0,
     .                MPI_Comm_World,MPIerror)
      lvec(1:3)=pb(1:3)
#endif

      if (ionode) then
           l_sqr = sqrt(lvec(1)**2+lvec(2)**2+lvec(3)**2)
           write(6,'(100("-"))')
           write(6,'(a4,15x,f8.3,2x,3f7.3/)')
     .           'Sum', l_sqr,lvec
      endif
  
C Deallocate local memory
      call memory('D','D',size(l_orb),'moments')
      deallocate(l_orb)

      end subroutine moments

