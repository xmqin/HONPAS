! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine kgrid( cell, kscell, displ, nk, points, weight )
c **********************************************************************
c Finds Monkhost-Pack k-point coordinates and weights.
c This version assumes no symmetry except time reversal, i.e.
c inversion in reciprocal space.
c Refs: H.J.Monkhorst and J.D.Pack, Phys Rev B 13, 5188 (1976)
c       J.Moreno and J.M.Soler, Phys Rev B 45, 13891 (1992)
c Written by J.M.Soler. July 1997.
c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c integer kscell(3,3): Supercell reciprocal of k-grid unit cell
c                      scell(ix,i) = sum_j cell(ix,j)*kscell(j,i)
c real*8  displ(3)   : Grid origin in k-grid-vector coordinates:
c                      origin(ix) = sum_j gridk(ix,j)*displ(j)
c integer nk         : Dimension of arrays points and weight
c ***************** OUTPUT *********************************************
c real*8  points(3,nk) : K-point cartesian coordinates
c                        Only if input_nk .ge. output_nk
c real*8  weight(nk)   : K-point weights 
c                        Only if input_nk .ge. output_nk
c ***************** UNITS **********************************************
c points returned in units reciprocal of cell
c **********************************************************************
C
C  Modules
C
      use units,        only : Ang
      use precision
      use parallel,     only : Node
      use fdf

      implicit          none

      integer, intent(inout)  ::          kscell(3,3), nk
      real(dp), intent(in)    ::          cell(3,3)
      real(dp), intent(inout) ::          displ(3)
      real(dp), intent(out)   ::          points(3,*), weight(*)
c ----------------------------------------------------------------------

c Internal variables

      character(len=24), save        :: accum

      integer           i, i1, i2, i3, igmax(3), igmin(3),
     .                  ir, ik, ix, j,
     .                  kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3),
     .                  ng(3), ni, nkr(3), nktot, proj(3,3)
      real(dp)          d(3), dkg(3), dkx(3), dscell(3,3),
     .                  gridk(3,3), gscell(3,3), pi,
     .                  scell(3,3), tiny, w1, wtot

      type(block_fdf) :: bfdf

      external          idiag, reclat
      logical           spiral
      parameter (tiny   = 1.d-12)

C Find out if there is spiral arrangement of spins
      spiral = fdf_isblock('SpinSpiral')

C Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )
C Find k-grid supercell
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) +
     .                  cell(ix,2) * kscell(2,i) +
     .                  cell(ix,3) * kscell(3,i)
        enddo
      enddo

C Find equivalent diagonal supercell
      call idiag( 3, kscell, kdsc, ml, mr, maux )
      proj(:,:) = 0  ! Possible sign changes
      do i = 1, 3
         proj(i,i) = 1
         if (kdsc(i,i) < 0) proj(i,i) = -1
      enddo
      kdsc = matmul(kdsc,proj)
      mr = matmul(mr,proj)
      dscell = matmul(scell,mr)

C Find k-grid unit vectors
      call reclat( dscell, gridk, 1 )
 
C Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) +
     .            gscell(ix,2) * displ(2) +
     .            gscell(ix,3) * displ(3)
      enddo

C Find grid origin in gridk coordinates
      pi = 4.d0 * atan(1.d0)
      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) +
     .             dkx(2) * dscell(2,i) +
     .             dkx(3) * dscell(3,i) ) / (2*pi)
      enddo

C Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
        igmin(j) = -( (ng(j)-1) / 2)
        igmax(j) = ng(j) / 2
      enddo

C Find number of points with time-reversal (inversion) symmetry,
C (if not spiral) after reflection on each alternative plane
      if (spiral) then
C       If spin-spiral, use all k points
        nk = nktot
      else
        do j = 1,3
          ni = ng(j)
          if (abs(dkg(j)) .lt. tiny) then
            ni = ng(j)/2 + 1
          elseif (abs(dkg(j)-0.5d0) .lt. tiny) then
            ni = (ng(j)-1)/2 + 1
          endif
          ! To work around an Intel_12 compiler bug
          write(accum,"(3i8)") ni,nktot,kdsc(j,j)
          nkr(j) = ni * nktot / kdsc(j,j)
        enddo

C Select reflection plane
        ir = 3
        if (nkr(2) .lt. nkr(ir)) ir = 2
        if (nkr(1) .lt. nkr(ir)) ir = 1
        igmin(ir) = 0
        if (abs(dkg(ir)-0.5d0) .lt. tiny)
     .    igmax(ir) = (ng(ir)-1)/2
        nk = nkr(ir)
      endif

C Find k points and weights
      w1 = 1.0d0 / nktot
      nk = 0
      do i3 = igmin(3),igmax(3)
      do i2 = igmin(2),igmax(2)
      do i1 = igmin(1),igmax(1)
        nk = nk + 1
        d(1) = i1 + dkg(1)
        d(2) = i2 + dkg(2)
        d(3) = i3 + dkg(3)
        if (d(1) .gt. 0.5d0*ng(1)+tiny) d(1) = d(1) - ng(1)
        if (d(2) .gt. 0.5d0*ng(2)+tiny) d(2) = d(2) - ng(2)
        if (d(3) .gt. 0.5d0*ng(3)+tiny) d(3) = d(3) - ng(3)
        do ix = 1,3
          points(ix,nk) = gridk(ix,1)*d(1) + 
     .                    gridk(ix,2)*d(2) +
     .                    gridk(ix,3)*d(3)
        enddo
        if (spiral) then
          weight(nk) = w1
        else
          if ( abs(d(ir))              .lt. tiny .or.
     .         abs(d(ir)-0.5d0*ng(ir)) .lt. tiny) then
            weight(nk) = w1
          else
            weight(nk) = 2.0d0 * w1
          endif
        endif
      enddo
      enddo
      enddo

C Check that weight is normalised
      wtot = 0.0d0
      do ik = 1,nk
        wtot = wtot + weight(ik)
      enddo
      if (abs(wtot-1.0d0) .gt. nk*tiny) then
        w1 = dble(nk)/wtot
        do ik = 1,nk
          weight(ik) = w1*weight(ik)
        enddo
      endif

      return
      end
