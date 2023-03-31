! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine phirphi(nua, na, nuo, no, scell, xa, rmaxo,
     .                   maxnh, lasto, iphorb, isa, 
     .                   numh, listhptr, listh, dk, S)
C *********************************************************************
C Finds the matrix elements of the dk*r between basis
C orbitals, where dk is a given vector and r is the position operator.
C
C     S_a,b(R_a,R_b) = 0.5*( <phi_a(r-R_a)|dk*(r-R_b)|phi_b(r-R_b)>
C                      +      <phi_b(r-R_b)|dk*(r-R_a)|phi_a(r-R_a)> )
C
C Energies in Ry. Lengths in Bohr.
C Writen by DSP March 1999 ( Based in routine overfsm )
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer nuo              : Number of basis orbitals in unit cell
C integer no               : Number of basis orbitals 
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C real*8  rmaxo            : Maximum cutoff for atomic orbitals
C integer maxnh            : First dimension of listh
C integer lasto(0:na)      : Last orbital index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer isa(na)          : Species index of each atom
C integer numh(nuo)        : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listhptr(nuo)    : Pointer to the start of each row 
C                            of the overlap matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C real*8  dk(3)            : Vector in k-space
C **************************** OUTPUT *********************************
C real*8  S(maxnh)         : Sparse overlap matrix
C *********************************************************************

      use precision
      use atmfuncs,     only : rcut, orb_gindex
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb
      use alloc,        only : re_alloc, de_alloc
      use neighbour,    only : jna=>jan, xij, r2ij
      use neighbour,    only : mneighb, reset_neighbour_arrays
      use m_new_matel,  only : new_matel

      implicit none

C Passed variables 
      integer  maxnh, na, no, nua, nuo

      integer
     .  iphorb(no), isa(na), lasto(0:na), 
     .  listh(maxnh), listhptr(nuo), numh(nuo)
      
      real(dp)
     .  scell(3,3), rmaxo, dk(3), S(maxnh), xa(3,na)

C Internal variables 
      integer
     .  ia, ind, iio, io, ioa, is, ix, ig, jg,
     .  j, ja, jn, jo, joa, js, nnia

      real(dp)
     .  grSij(3), rij, Sij, xinv(3)

      real(dp), dimension(:), pointer ::  Si

      real(dp), parameter  :: tiny = 1.0d-9

C ......................
      call timer( 'phirphi', 1 )

C Initialize neighb subroutine 
      call mneighb( scell, 2.d0*rmaxo, na, xa, 0, 0, nnia )

C Allocate local memory
      nullify( Si )
      call re_alloc( Si, 1, no, name='Si',
     &               routine='subroutine phirphi' )

      do jo = 1,no
        Si(jo) = 0.0d0
      enddo
      do ia = 1,nua 
        is = isa(ia)
        call mneighb( scell, 2.0d0*rmaxo, na, xa, ia, 0, nnia )
           
        do iio = lasto(ia-1)+1,lasto(ia)  
          call GlobalToLocalOrb(iio,Node,Nodes,io)
          if (io .gt. 0) then
            ioa = iphorb(iio)
            ig = orb_gindex(is,ioa)
            do jn = 1,nnia 
              do ix = 1,3
                xinv(ix) = - xij(ix,jn)
              enddo 
              ja = jna(jn)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja) 
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then  
                   jg = orb_gindex(js,joa)

                  if (abs(dk(1)).gt.tiny) then
                    call new_MATEL('X', ig, jg, xij(1:3,jn),
     .                          Sij, grSij ) 
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(1)  
 
                    call new_MATEL('X', jg, ig, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(1)  
                  endif
                     
                  if (abs(dk(2)).gt.tiny) then
                    call new_MATEL('Y', ig, jg, xij(1:3,jn),
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(2) 
                
                    call new_MATEL('Y', jg, ig, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(2)  
                  endif
 
                  if (abs(dk(3)).gt.tiny) then
                    call new_MATEL('Z', ig, jg, xij(1:3,jn),
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(3) 
 
                    call new_MATEL('Z', jg, ig, xinv,
     .                          Sij, grSij )
                    Si(jo) = Si(jo) + 0.5d0*Sij*dk(3) 
                  endif
                endif
              enddo
            enddo
            do j = 1,numh(io)
              ind = listhptr(io)+j
              jo = listh(ind)
              S(ind) = Si(jo) 
              Si(jo) = 0.0d0 
            enddo
          endif
        enddo
      enddo

C Deallocate local memory
!      call new_MATEL('Z', 0, 0, 0, 0, xinv, Sij, grSij )
      call reset_neighbour_arrays( )
      call de_alloc( Si, name='Si' )

      call timer( 'phirphi', 2 )
      end
