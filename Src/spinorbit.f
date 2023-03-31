! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

C *********************************************************************
C
C SPIN-ORBIT INTERACTION ---  ON-SITE APPROXIMATION 
C
C The spin-orbit hamiltonian has the form:
C 
C                        | Lz             Lx - i Ly |
C  <i|HSO|j> = <i| V(r)  |                          | |j>
C                        | Lx + i Ly      - Lz      |
C
C                        | i*L1(li,mi,mj)         L2(li,mi,mj)-i*L3(li,mi,mj) |
C  <i|HSO|j> = 0.5 M(li) |                                                    |
C                        | -L2(li,mi,mj)-i*L3(li,mi,mj)       -i*L1(li,mi,mj) |
C
C where M(li) is the radial part, and Li(li,mi,mj) are the radial bits
C 
C Hence,  <i|Lz|j> =  i L1
C         <i|Lx|j> = -i L3
C         <i|Ly|j> =  i L2
C
C *********************************************************************
C
C   Written by Jaime Ferrer, Lucas Fernandez Seivane, Universidad de Oviedo
C              Miguel Angel Oliveria, Trinity Colleg, Dublin
C
C   Modified by Alberto Garcia
C
      module spinorbit

      use precision, only: dp
      use atmfuncs
      use atm_types
      use fdf, only: fdf_get
      use parallel, only: Node, Nodes
      use parallelsubs, only: LocalToGlobalOrb

      implicit none
!
!     Interim data structures to hold Vso information
!     which was known only to the master node
!
!     Wasteful, as we hold values for r > rcut
!     It would be better to make vso into some
!     kind of extra projector (with the forthcoming
!     indexing technology)

      logical, save  :: vso_setup = .false.
      real(dp), save :: so_strength = 1.0_dp

      integer, pointer, save   :: nr(:)
      real(dp), pointer, save  :: vso(:,:,:)
      real(dp), pointer, save  :: r(:,:)
      real(dp), pointer, save  :: drdi(:,:)

      public :: spinorb, int_so_ang
      private

      CONTAINS

      subroutine init_vso()

      use basis_types,     only: basis_parameters
      use atm_types,       only: nspecies
      use pseudopotential, only: pseudopotential_t
      use atmparams,       only: lmaxd
      use parallel,        only: Node
      use m_mpi_utils,     only: broadcast
      use sys,             only: die

      integer  :: is, mx_nrval, li, ir, iup
      real(dp) :: a, b, rpb, ea
      logical  :: there_are_so_potentials
      type(pseudopotential_t), pointer :: p

      allocate(nr(nspecies))
      mx_nrval = 0
      if (Node .eq. 0) then
         do is = 1, nspecies
            p=> basis_parameters(is)%pseudopotential
            nr(is) = p%nrval
            ! Could get maximum rcut for is's orbitals
            ! and compute maximum effective nrval
            mx_nrval = max(mx_nrval,nr(is))
         enddo
      endif

      call broadcast(mx_nrval)
      call broadcast(nr)
      allocate(vso(mx_nrval,1:lmaxd,nspecies))
      allocate(r(mx_nrval,nspecies))
      allocate(drdi(mx_nrval,nspecies))

      vso = 0.0_dp
      r   = 0.0_dp
      drdi = 0.0_dp

      if (Node .eq. 0) then
         write(6,"(/,a)")
     $        "Initializing spin-orbit part of the Hamiltonian"
         there_are_so_potentials = .false.
         do is = 1, nspecies
            p=> basis_parameters(is)%pseudopotential
            if ((p%irel == "rel") .and. (p%npotu > 0)) then
               write(6,"(a)") "  Adding spin-orbit effects for "
     $                         // p%name
               there_are_so_potentials = .true.
               do iup = 1, p%npotu
                  li = p%lup(iup)
                  vso(1:nr(is),li,is) = p%vup(li,:)
               enddo
            else
               ! No spin-orbit components for this species
               vso(1:nr(is),:,is) = 0.0_dp
            endif
            r(1:nr(is),is) = p%r(:)
            a = p%a      
            b = p%b
            rpb=b
            ea=exp(a)
            do ir=1,nr(is)
               drdi(ir,is)=a*rpb
               rpb=rpb*ea
            enddo
         enddo
         write(6,"(/)")
         if (.not. there_are_so_potentials) then
            call die("No spin-orbit components for any species!")
         endif
      endif

      call broadcast(vso)
      call broadcast(r)
      call broadcast(drdi)
      
      so_strength = fdf_get('SpinOrbitStrength',1.0_dp)

      end subroutine init_vso

!------------------------------------------------

      subroutine spinorb(no_u,no_l,iaorb,iphorb,isa,indxuo,
     $                   maxnh,numh,listhptr,listh,Dscf,H,Eso)
     $                   
C *********************************************************************
C Spin-orbit contributions to matrix elements.
C Energies in Ry. Lengths in Bohr.
C **************************** INPUT **********************************
C integer no_u             : Number of orbitals in unit cell
C integer no_l             : Number of orbitals in unit cell local to node
C integer maxnh            : First dimension of H and listh (and D)
C integer iphorb(no_u)     : Orbital index of each orbital in its atom
C integer iaorb(no_u)      : Atom to which each orbital belongs
C integer indxuo(*)        : Index of equivalent unit-cell orbital
C integer isa(*)           : Species index of each atom
C integer numh(no_l)       : Number of nonzero elements of each row
C                            of the hamiltonian matrix
C integer listhptr(no_l)   : Pointer to the start of rows (-1) of
C                            the hamiltonian matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements
C                            of each row of the hamiltonian matrix
C integer Dscf(maxnh,8)    : Density matrix
C **************************** INPUT / OUTPUT *********************************
C real*8 Eso:		   : Spin-orbit energy
C **************************** INPUT / OUTPUT *********************************
C real*8 H(maxnh,3:8)      : Spin-Orbit H matrix elements
C *********************************************************************
C
      use m_mpi_utils, only: globalize_sum
      implicit none

C Arguments

      integer, intent(in) :: no_u, no_l, maxnh,
     .     iphorb(no_u), iaorb(no_u),
     .     listh(maxnh), numh(no_l), listhptr(no_l),
     .     indxuo(:), isa(:)

      real(dp), intent(in)    :: Dscf(maxnh,8)
      real(dp), intent(inout) :: H(maxnh,3:8), Eso


C Internal variables

      integer  ::  ia, ja, ioa, is, joa, j, ind, 
     .             li, lj, mi, mj
      integer  ::  io_l, io_u, jo_s, jo_u, ih

      real(dp) :: int_rad, int_ang(1:3), buffer
      real(dp) :: Hso_ji(3)

!------------------------------------------- BEGIN

      if (.not. vso_setup) then
         call init_vso()
         vso_setup = .true.
      endif

      call timer( 'spinorb', 1 )

!  AG
!     On-site approximation: Only matrix elements between orbitals
!     on the same atom are considered.
!     Of these, only those with the same l
!     Of these, only those with different m's
!          (Aside: What happens to double-z orbitals??)
!
!     So Hso will be VERY sparse.

      Eso = 0.0_dp
      do io_l = 1, no_l
         call LocalToGlobalOrb(io_l,Node,Nodes,io_u)
         ia = iaorb(io_u)
         ioa = iphorb(io_u)
         is = isa(ia)
         li=lofio(is,ioa)
         if (li == 0) CYCLE  ! No contribution from l=0 orbs
         mi=mofio(is,ioa)
         do j = 1,numh(io_l)
            ind = listhptr(io_l) + j
            jo_s = listh(ind)
            jo_u = indxuo(jo_s)
            if (jo_s .ne. jo_u) CYCLE ! Not in the unit cell
            ja = iaorb(jo_u)
            if (ja .ne. ia) CYCLE ! Not in the same atom
            joa = iphorb(jo_u)
            lj=lofio(is,joa)
            if (li /= lj) CYCLE ! Different l
            mj=mofio(is,joa)
            if (mi == mj) CYCLE ! Same m
            
            call int_so_rad(is, li, joa, ioa, int_rad)
            call int_so_ang(li, mj, mi, int_ang(:))
            Hso_ji(:)=so_strength*int_rad*int_ang(:)

            H(ind,3) = H(ind,3) + Hso_ji(2)
            H(ind,4) = H(ind,4) + Hso_ji(3)
            H(ind,5) = H(ind,5) + Hso_ji(1)
            H(ind,6) = H(ind,6) - Hso_ji(1)
            H(ind,7) = H(ind,7) - Hso_ji(2)
            H(ind,8) = H(ind,8) - Hso_ji(3)

            Eso = Eso + Hso_ji(2) * (-Dscf(ind,3)+Dscf(ind,7)) +
     .                  Hso_ji(3) * (-Dscf(ind,4)+Dscf(ind,8)) +
     .                  Hso_ji(1) * (-Dscf(ind,5)+Dscf(ind,6))

         enddo
      enddo

      ! Globalzie Eso
      call globalize_sum(Eso, buffer)
      Eso = buffer
      
      call timer( 'spinorb', 2 )
      end subroutine spinorb


C *********************************************************************
C
C Subroutine to calculate the spin-orbit angular integral
C Calculates L1(li,mi,mj), L2(li,mi,mj) and L3(li,mi,,mj) 
C
C *********************************************************************

      subroutine int_so_ang(li, mi, mj, L)

      implicit none

      integer, intent(in)  :: li, mi, mj
      real(dp),intent(out) :: L(3)

      real(dp) ::  La, Lb, Lc
      real(dp), parameter :: one = 1._dp, two = 2._dp
      real(dp), parameter :: six = 6._dp



      L(1:3)= 0.0_dp

      La = sqrt(li*(li+1._dp)/2._dp)
      Lb = sqrt(li*(li+1._dp)-2._dp)/2._dp
      Lc = 0.0_dp
      if (li .ge. 3) Lc = sqrt(li*(li+1._dp)-6._dp)/2._dp

      if((mi+mj).EQ.0) L(1) = real(mj, dp)

      select case ( mi )
      case ( 0 )
         select case ( mj )
         case ( -1 )
            L(3) = La
         case ( 1 )
            L(2) = La
         end select
      case ( 1 )
         select case ( mj )
         case ( -2 )
            L(3) =  Lb
         case ( 0 )
            L(2) = -La
         case ( 2 )
            L(2) =  Lb
         end select
      case ( -1 )
         select case ( mj )
         case ( -2 )
            L(2) =  Lb
         case ( 0 )
            L(3) = -La
         case ( 2 )
            L(3) = -Lb
         end select
      case ( 2 )
         select case ( mj )
         case ( -3 )
            L(3) =  Lc
         case ( -1 )
            L(3) =  Lb
         case ( 1 )
            L(2) = -Lb
         case ( 3 )
            L(2) =  Lc
         end select
      case ( -2 )
         select case ( mj )
         case ( -3 )
            L(2) =  Lc
         case ( -1 )
            L(2) = -Lb
         case ( 1 )
            L(3) = -Lb
         case ( 3 )
            L(3) = -Lc
         end select
      case ( 3 )
         select case ( mj )
         case ( -2 )
            L(3) =  Lc
         case ( 2 )
            L(2) = -Lc
         end select
      case ( -3 )
         select case ( mj )
         case ( -2 )
            L(2) = -Lc
         case ( 2 )
            L(3) = -Lc
         end select
      end select

      end subroutine int_so_ang

C *********************************************************************
C
C Subroutine to calculate the radial spin-orbit integral
C Calculates 0.5*M(li)=0.5*<i|V_{li}^{SO}|i>
C
C *********************************************************************

      subroutine int_so_rad(is, li, ioa, joa, result)
      use atmfuncs, only: rcut

      implicit none

      integer,intent(in) :: is, li, ioa, joa
      real(dp),intent(out) :: result

      integer nrval, ir
      real(dp) fi, fj, grad, rr, rmax

      nrval = nr(is)
      rmax = rcut(is,ioa)
      rmax = max(rmax,rcut(is,joa))

      result = 0.0_dp
      do ir = 2, nrval
        rr = r(ir,is)
        if (rr > rmax) exit
        call rphiatm(is,ioa,rr,fi,grad)
        call rphiatm(is,joa,rr,fj,grad)
        result = result+fi*vso(ir,li,is)*fj*rr*drdi(ir,is)
      enddo

      result = 0.5_dp*result

      end subroutine int_so_rad

      end module spinorbit
