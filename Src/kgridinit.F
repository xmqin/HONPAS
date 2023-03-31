! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine kgridinit( cell, kscell, displ, cutoff, nk )

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
c real*8  cutoff     : Minimum k-grid cutoff required
c                      Not used unless det(kscell)=0
c integer nk         : Dimension of arrays points and weight
c ***************** OUTPUT *********************************************
c integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
c                        Only if det(kscell)=0 on input
c real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
c                        Only if det(kscell)=0 on input
c real*8  cutoff       : Actual k-grid cutoff
c integer nk           : Actual number of irreducible k-points
c ***************** UNITS **********************************************
c cutoff must be in the same units as cell
c ***************** BEHAVIOUR ******************************************
c - If det(kscell).ne.0, input cutoff is not used
c - If det(kscell).eq.0, kscell and displ are generated according to
c   input cutoff
c - If det(kscell).eq.0 .AND. cutoff.le.0.d0, they are readed
c   from the input fdf data file.
c   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
c   If both are present, kgrid_Monkhorst_Pack has priority. If none is
c   present, the cutoff defect is cero, producing only the gamma point.
c   Examples of fdf data specifications:
c     kgrid_cutoff  50. Bohr
c     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
c     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
c     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
c     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
c     %endblock kgrid_Monkhorst_Pack
c **********************************************************************
C
C  Modules
C
      use precision,  only : dp
      use units,      only : Ang
      use parallel,   only : Node
      use m_minvec,   only : minvec
      use sys,        only : die
      use fdf

      implicit          none

      character(len=24), save         :: accum

C Passed variables
      integer, intent(inout)   :: kscell(3,3)
      integer, intent(out)     :: nk
      real(dp), intent(in)     :: cell(3,3)
      real(dp), intent(inout)  :: cutoff, displ(3)

      external                 ::  idiag, reclat
      real(dp), external       ::  volcel

C Internal variables
      integer           i, ir, ix, j,
     .                  kdsc(3,3), maux(3,3,2), ml(3,3), mr(3,3),
     .                  ng(3), ni, nkr(3), nktot, factor(3,3),
     $                  proj(3,3)
      real(dp)          defcut, dkg(3), dkx(3), dscell(3,3),
     .                  gridk(3,3), gscell(3,3), huge, pi,  
     .                  scell(3,3), scmin(3,3), tiny, vmod
      real(dp)          ctransf(3,3)
      logical           genlogic, spiral

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

      parameter (defcut = 0.d0)
      parameter (huge   = 1.d30)
      parameter (tiny   = 1.d-12)

      genlogic = .false.
C Find out if there is spiral arrangement of spins 
      spiral = fdf_isblock('SpinSpiral')

C Find total number of points (determinant of kscell)
      nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .             kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .             kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .             kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .             kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .             kscell(3,1) * kscell(2,2) * kscell(1,3) )

C Look for kscell or cutoff in input fdf file
      if ( nktot.eq.0 .and. cutoff.lt.tiny ) then
        genlogic = fdf_block('kgrid_Monkhorst_Pack',bfdf)
        if ( genlogic ) then
          do i= 1, 3
            if (.not. fdf_bline(bfdf,pline))
     .        call die('kgridinit: ERROR in kgrid_Monkhorst_Pack block')
            kscell(1,i) = fdf_bintegers(pline,1)
            kscell(2,i) = fdf_bintegers(pline,2)
            kscell(3,i) = fdf_bintegers(pline,3)
            if ( fdf_bnvalues(pline) > 3 ) then
               displ(i) = mod(fdf_bvalues(pline,4), 1._dp)
            else
               displ(i) = 0._dp
            end if
          enddo
          call fdf_bclose(bfdf)
          nktot = abs( kscell(1,1) * kscell(2,2) * kscell(3,3) +
     .                 kscell(2,1) * kscell(3,2) * kscell(1,3) +
     .                 kscell(3,1) * kscell(1,2) * kscell(2,3) -
     .                 kscell(1,1) * kscell(3,2) * kscell(2,3) -
     .                 kscell(2,1) * kscell(1,2) * kscell(3,3) -
     .                 kscell(3,1) * kscell(2,2) * kscell(1,3) )
        else
C The second argument is the default value
          cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
        endif
      endif


C Find kscell from required cutoff (i.e., if nktot is still zero)
C If kgrid_cutoff is zero, we will end up with nktot = 1
C and kscell as the unit matrix.
C
      if (nktot .eq. 0) then
        call minvec( cell, scmin, ctransf )
!        print "(a,/,3(3f15.8,/))", "Rounded unit cell vectors: ",
!     $         ((scmin(i,j)/Ang,i=1,3),j=1,3)
        nktot = 1
        do j = 1,3
          factor(j,1:3) = 0
          vmod = sqrt( scmin(1,j)**2 + scmin(2,j)**2 + scmin(3,j)**2 )
          factor(j,j) = int(cutoff/(vmod/2.d0)) + 1
CC          if (mod(factor(j),2) .eq. 0) displ(j) = 0.5d0
          nktot = nktot * factor(j,j)
        enddo
C
C       Generate actual supercell skeleton
C
!!        print "(a,3i4)", "Factors from rounded cell: ", 
!!     $         (factor(j,j),j=1,3)
        kscell = matmul(ctransf, factor)
!!        print "(a,/,3(3i4,/))", "Kscell from required cutoff: ",
!!     $         ((kscell(i,j),j=1,3),i=1,3)

      endif
C 
C Now we have kscell, either from MP input, or from kgrid_cutoff input
C Find k-grid supercell
C
!!      scell = matmul(cell,kscell)
      do i = 1,3
        do ix = 1,3
          scell(ix,i) = cell(ix,1) * kscell(1,i) +
     .                  cell(ix,2) * kscell(2,i) +
     .                  cell(ix,3) * kscell(3,i)
        enddo
        vmod = sqrt( scell(1,i)**2 + scell(2,i)**2 + scell(3,i)**2 )
!!        print *, i, vmod/Ang
      enddo
!!      print "(a,/,3(3f15.8,/))",  "K grid supercell vectors: ",
!!     $         ((scell(i,j)/Ang,i=1,3),j=1,3)

C Find actual cutoff
      call minvec( scell, scmin, ctransf )
      cutoff = huge
      do i = 1,3
        vmod = sqrt( scmin(1,i)**2 + scmin(2,i)**2 + scmin(3,i)**2 )
        cutoff = min( cutoff, vmod/2.d0 )
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
!!      print *, "Eq. diagonal supercell: ", (kdsc(i,i),i=1,3)
!!      print *, "Mright: ",
!!     $         ((mr(i,j),j=1,3),i=1,3)
!!      print *, "Mleft: ",
!!     $         ((ml(i,j),j=1,3),i=1,3)
!!      print *, "check:", matmul(ml,matmul(kscell,mr))
C
C     Set the displacements
C
      do j = 1, 3
         if (mod(kdsc(j,j),2) .eq. 0) then
            if (genlogic .and. displ(j) /= 0.5d0) then
               if (Node .eq. 0)
     $           write(6,"(a,i4,a,2f8.2)")
     $            "k-point displ. along", j, " input, could be: ",
     $           displ(j), 0.5d0
            else
               displ(j) = 0.5d0
            endif
         else
            if (genlogic .and. displ(j) /= 0.0d0) then
               if (Node .eq. 0)
     $           write(6,"(a,i4,a,2f8.2)")
     $            "k-point displ. along", j, " input, could be: ",
     $           displ(j), 0.0d0
            else
               displ(j) = 0.0d0
            endif
         endif
      enddo
C
C     Equivalent supercell DA with the property that there exists
C     a primitive cell (pa') such that DA_i = N_i*pa'_i
C     (See Moreno and Soler)
C
      dscell = matmul(scell,mr)
!!      print *, "Dscell cell vectors: ",
!!     $         ((dscell(i,j)/Ang,i=1,3),j=1,3)

C Find k-grid unit vectors
      call reclat( dscell, gridk, 1 )

C Find grid origin in cartesian coordinates
      call reclat( scell, gscell, 1 )
      do ix = 1,3
        dkx(ix) = gscell(ix,1) * displ(1) +
     .            gscell(ix,2) * displ(2) +
     .            gscell(ix,3) * displ(3)
      enddo
!!      print *, "Grid origin in cartesian coords ", (dkx(i), i=1,3)

C Find grid origin in gridk coordinates
      pi = 4.d0 * atan(1.d0)
      do i = 1,3
        dkg(i) = ( dkx(1) * dscell(1,i) +
     .             dkx(2) * dscell(2,i) +
     .             dkx(3) * dscell(3,i) ) / (2*pi)
      enddo
!!      print *, "Grid origin in grid coords ", (dkg(i), i=1,3)

C Some printout for debugging
!!      write(6,'(/,a,/,(3f12.6,i6,f12.6))') 'kgrid: gridK,ng,dg =',
!!     .  ((gridk(ix,i),ix=1,3),kdsc(i,i),dkg(i),i=1,3)

C Find total range of grid indexes
      do j = 1,3
        ng(j) = kdsc(j,j)
      enddo

C Find number of points with time-reversal (inversion) symmetry,
C (if not spiral) after reflection on each alternative plane
      if (spiral) then
C If spiral, use all k points
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

C Find k points and weights
        nk = nkr(ir)
      endif
!!      print *, "nk: ", nk

      return
      end
