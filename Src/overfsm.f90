! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_overfsm
  
  use precision,     only : dp
  use parallel,      only : Node, Nodes
  use parallelsubs,  only : GlobalToLocalOrb
  use atmfuncs,      only : rcut, orb_gindex
  use neighbour,     only : jna=>jan, r2ij, xij, mneighb, reset_neighbour_arrays
  use alloc,         only : re_alloc, de_alloc
  use m_new_matel,   only : new_matel
  use t_spin, only: tSpin
  
  implicit none

  private
  public :: overfsm
  
contains
  
  subroutine overfsm(na_u, na_s, no_l, no_s, scell, xa, indxua, rmaxo, &
      lasto, iphorb, isa, maxnd, numd, listdptr, listd, &
      spin, Escf, fa, stress)
    
    ! *********************************************************************
    ! Overlap matrix and orthonormalization contribution to forces and stress.
    ! Energies in Ry. Lengths in Bohr.
    ! Writen by J.Soler and P.Ordejon, August-October'96, June'98
    ! **************************** INPUT **********************************
    ! integer na_u             : Number of atoms in unit cell
    ! integer na_s             : Number of atoms in supercell
    ! integer no_l             : Number of orbitals in unit-cell (local)
    ! integer no_s             : Number of orbitals in supercell
    ! real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
    ! real*8  xa(3,na_s)       : Atomic positions in cartesian coordinates
    ! integer indxua(na_s)     : Index of equivalent atom in unit cell
    ! real*8  rmaxo            : Maximum cutoff for atomic orbitals
    ! integer maxnd            : First dimension of Escf and listd
    ! integer lasto(0:na_s)    : Last orbital index of each atom
    ! integer iphorb(no_s)     : Orbital index of each orbital in its atom
    ! integer isa(na_s)        : Species index of each atom
    ! integer numd(no_l)       : Number of nonzero elements of each row
    !                            of Escf
    ! integer listdptr(no_l)   : Pointer to start of rows (-1) of Escf
    ! integer listd(maxnd)     : Column indexes of the nonzero elements  
    !                            of each row of Escf
    ! type(tSpin) spin         : Spin configuration
    ! integer Escf(maxnd,spin%EDM): Energy-Density matrix
    ! ********************** INPUT and OUTPUT *****************************
    ! real*8  fa(3,na_u)       : Atomic forces (Orthog. part added to input)
    ! real*8  stress(3,3)      : Stress tensor (Orthog. part added to input)
    ! *********************************************************************
    
    
    integer, intent(in) :: na_u, na_s, no_l, no_s
    integer, intent(in) :: indxua(na_s), iphorb(no_s), isa(na_s), lasto(0:na_s)
    integer, intent(in) :: maxnd
    integer, intent(in) :: listd(maxnd), numd(no_l), listdptr(no_l)
    
    real(dp), intent(inout) :: fa(3,na_u), stress(3,3)
    type(tSpin), intent(in) :: spin
    real(dp), intent(in) :: scell(3,3), Escf(maxnd,spin%EDM), rmaxo, xa(3,na_s)

    ! Internal variables ......................................................
  
    integer :: ia, ind, io, ioa, is, ispin, ix, iio, ig, jg
    integer :: j, ja, jn, jo, joa, js, jua, jx, nnia
    
    real(dp) :: fij(3), grSij(3), rij, Sij, volcel, volume

    real(dp), dimension(:), pointer ::  Di => null()

    external :: timer

    ! Start timer
    call timer( 'overfsm', 1 )
    
    volume = na_u * volcel(scell) / na_s
    
    ! Initialize neighb subroutine 
    call mneighb( scell, 2.d0*rmaxo, na_s, xa, 0, 0, nnia)
    
    ! Allocate local memory
    call re_alloc( Di, 1, no_s, 'Di', 'overfsm' )
    
    Di(1:no_s) = 0.0d0
    
    do ia = 1, na_u
      is = isa(ia)
      call mneighb( scell, 2.d0*rmaxo, na_s, xa, ia, 0, nnia )
      do io = lasto(ia-1)+1,lasto(ia)
        
        ! Is this orbital on this Node?
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        if ( iio > 0 ) then
          
          ! Valid orbital 
          ioa = iphorb(io)
          ig = orb_gindex(is,ioa)
          do j = 1,numd(iio)
            ind = listdptr(iio)+j
            jo = listd(ind)
            do ispin = 1, spin%spinor  
              Di(jo) = Di(jo) + Escf(ind,ispin)
            end do
          end do
          
          do jn = 1,nnia
            ja = jna(jn)
            jua = indxua(ja)
            rij = sqrt( r2ij(jn) )
            do jo = lasto(ja-1) + 1 , lasto(ja)
              joa = iphorb(jo)
              js = isa(ja)
              
              if ( rcut(is,ioa) + rcut(js,joa) > rij ) then
                jg = orb_gindex(js,joa)
                call new_MATEL( 'S', ig, jg, xij(1:3,jn), Sij, grSij )
                do ix = 1,3
                  fij(ix) = - Di(jo) * grSij(ix)
                  fa(ix,ia)  = fa(ix,ia)  + fij(ix)
                  fa(ix,jua) = fa(ix,jua) - fij(ix)
                  do jx = 1,3
                    stress(jx,ix) = stress(jx,ix) + xij(jx,jn) * fij(ix) / volume
                  end do
                end do
              end if
              
            end do
          end do

          ! Reset arrays for next orbital interaction
          do j = 1,numd(iio)
            jo = listd(listdptr(iio)+j)
            Di(jo) = 0._dp
          end do
          
        end if
        
      end do
    end do
    
    ! Deallocate local memory
    call reset_neighbour_arrays( )
    call de_alloc( Di, 'Di', 'overfsm' )
    
    ! Finish timer
    call timer( 'overfsm', 2 )
    
  end subroutine overfsm

end module m_overfsm




