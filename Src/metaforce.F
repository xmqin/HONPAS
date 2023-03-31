! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module metaforce

      use precision, only : dp

      implicit none

      integer,           save :: nGauss = 0
      logical,           save :: lMetaForce = .false.
      integer,  pointer, save :: nGaussPtr(:,:) => null()
      real(dp), pointer, save :: GaussK(:) => null()
      real(dp), pointer, save :: GaussR0(:) => null()
      real(dp), pointer, save :: GaussZeta(:) => null()

      CONTAINS

      subroutine initmeta()
C
C  Reads data for metadynamics Gaussians
C
C  Julian Gale, NRI, Curtin University, Feb. 2006
C
      use parallel,     only : Node
      use alloc,        only : re_alloc
      use sys,          only : die
      use fdf

      implicit none

      integer            :: k
      integer            :: ni
      integer            :: nr

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C     Find if there are any Gaussians specified
      nGauss = fdf_block_linecount('MetaForce', 'iirrr')
      lMetaForce = nGauss > 0
      if ( lMetaForce ) then

C Allocate arrays for Gaussian data
        call re_alloc( nGaussPtr, 1, 2, 1, nGauss, 'nGaussPtr',
     &                 'initmeta' )
        call re_alloc( GaussK, 1, nGauss, 'GaussK', 'initmeta' )
        call re_alloc( GaussR0, 1, nGauss, 'GaussR0', 'initmeta' )
        call re_alloc( GaussZeta, 1, nGauss, 'GaussZeta', 'initmeta' )

C Read Gaussians from block
        lMetaForce = fdf_block('MetaForce', bfdf)

        do while ( fdf_bline(bfdf, pline) )
C Read and parse data line
          ni = fdf_bnintegers(pline)
          nr = fdf_bnreals(pline)

C Check that correct info is given
          if ((ni .ge. 2) .and. (nr .ge. 3)) then
            nGauss = nGauss + 1
            ! Atom 1
            nGaussPtr(1,nGauss) = abs(fdf_bintegers(pline,1))
            ! Atom 2
            nGaussPtr(2,nGauss) = abs(fdf_bintegers(pline,2))
            ! Gaussian pre-factor for the force
            ! Unit: Ry
            GaussK(nGauss) = fdf_breals(pline,1)
            ! Distance in r - r_0, where r is the distance between
            ! atom i and j
            ! Unit: Ang
            GaussR0(nGauss) = fdf_breals(pline,2)/0.529177_dp
            ! Gaussian prefactor
            ! Unit: 1/Bohr**2
            GaussZeta(nGauss) = fdf_breals(pline,3)
          endif
        enddo

      endif

      end subroutine initmeta

      subroutine meta(xa,na,ucell,Emeta,fa,stress,forces,stresses)

      implicit none

C
C  Compute energy and forces due to Gaussians
C
C  Julian Gale, NRI, Curtin University, Feb. 2006
C

C Passed variables
      integer,    intent(in)    :: na 
      logical,    intent(in)    :: forces
      logical,    intent(in)    :: stresses
      real(dp),   intent(in)    :: ucell(3,3)
      real(dp),   intent(in)    :: xa(3,na)
      real(dp),   intent(out)   :: Emeta
      real(dp),   intent(inout) :: fa(3,na)
      real(dp),   intent(inout) :: stress(3,3)

C Local variables
      integer                   :: ii
      integer                   :: i
      integer                   :: j
      integer                   :: k
      integer                   :: ixcell(27)
      integer                   :: iycell(27)
      integer                   :: izcell(27)
      real(dp)                  :: dtrm
      real(dp)                  :: exptrm
      real(dp)                  :: r
      real(dp)                  :: r2
      real(dp)                  :: r2min
      real(dp)                  :: rvol
      real(dp)                  :: vol
      real(dp)                  :: volcel
      real(dp)                  :: x
      real(dp)                  :: y
      real(dp)                  :: z
      real(dp)                  :: xc
      real(dp)                  :: yc
      real(dp)                  :: zc
      real(dp)                  :: xmin
      real(dp)                  :: ymin
      real(dp)                  :: zmin
      real(dp)                  :: xcell(27)
      real(dp)                  :: ycell(27)
      real(dp)                  :: zcell(27)

      if (stresses) then       
        vol = volcel(ucell)          
        rvol = 1.0_dp/vol           
      endif
C Initialise energy
      Emeta = 0.0_dp

C Find cell images
      call cellimagelist(ucell,xcell,ycell,zcell,ixcell,iycell,izcell)

C Loop over Gaussians 
      do k = 1,nGauss

C Set pointers to atoms
        i = nGaussPtr(1,k)
        j = nGaussPtr(2,k)

C Set initial coordinates
        x = xa(1,j) - xa(1,i)
        y = xa(2,j) - xa(2,i)
        z = xa(3,j) - xa(3,i)

C Get coordinates of nearest image
        r2min = 100000.0_dp
        do ii = 1,27
          xc = x + xcell(ii)
          yc = y + ycell(ii)
          zc = z + zcell(ii)

C Compute distance
          r2 = xc*xc + yc*yc + zc*zc
          if (r2.lt.r2min) then
            r2min = r2
            xmin = xc
            ymin = yc
            zmin = zc
          endif

        enddo

        r = sqrt(r2min)

C Compute energy and forces
        exptrm = exp(-GaussZeta(k)*(r - GaussR0(k))**2)
        Emeta = Emeta + 0.5_dp*GaussK(k)*exptrm
        dtrm = GaussK(k)*exptrm*GaussZeta(k)*(r - GaussR0(k))/r
        if (forces) then
          fa(1,i) = fa(1,i) - dtrm*xmin
          fa(2,i) = fa(2,i) - dtrm*ymin
          fa(3,i) = fa(3,i) - dtrm*zmin
          fa(1,j) = fa(1,j) + dtrm*xmin
          fa(2,j) = fa(2,j) + dtrm*ymin
          fa(3,j) = fa(3,j) + dtrm*zmin
        endif
        if (stresses) then
          dtrm = dtrm*rvol
          stress(1,1) = stress(1,1) + dtrm*xmin*xmin
          stress(2,1) = stress(2,1) + dtrm*ymin*xmin
          stress(3,1) = stress(3,1) + dtrm*zmin*xmin
          stress(1,2) = stress(1,2) + dtrm*xmin*ymin
          stress(2,2) = stress(2,2) + dtrm*ymin*ymin
          stress(3,2) = stress(3,2) + dtrm*zmin*ymin
          stress(1,3) = stress(1,3) + dtrm*xmin*zmin
          stress(2,3) = stress(2,3) + dtrm*ymin*zmin
          stress(3,3) = stress(3,3) + dtrm*zmin*zmin
        endif

      enddo

      end subroutine meta

      end module metaforce
