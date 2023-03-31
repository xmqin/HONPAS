! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_ksvinit
      USE precision, only : dp

      private

      integer, public, save           :: nkpol
      real(dp), pointer, public, save :: wgthpol(:) 
      real(dp), pointer, public, save :: kpol(:,:)
      logical, public, save           :: gamma_polarization

      public :: estimate_pol_kpoints, repol

      CONTAINS

      subroutine estimate_pol_kpoints(ucell)
      USE alloc, only: re_alloc

      implicit none

      real(dp), intent(in)   :: ucell(3,3)

      nullify(kpol,wgthpol)

      ! Find the grid for the calculation of the polarization
      nkpol = 1
      call re_alloc(kpol,1,3,1,nkpol,name='kpol',routine='siesta',
     .              copy=.false.)
      call re_alloc(wgthpol,1,nkpol,name='wgthpol',routine='siesta',
     .              copy=.false.)

      call KSV_init(ucell, maxkpol=0)

      call re_alloc(kpol,1,3,1,nkpol,name='kpol',routine='siesta',
     .              shrink=.false.,copy=.false.)
      call re_alloc(wgthpol,1,nkpol,name='wgthpol',routine='siesta',
     .              shrink=.false.,copy=.false.)

      gamma_polarization = (nkpol == 0) 
!!      print *, "nkpol, gamma_pol: ", nkpol, gamma_polarization
      end subroutine estimate_pol_kpoints
!
!------------------------------------------------------------------------
      subroutine KSV_init( ucell, maxkpol )
C *********************************************************************
C Finds polarization using the method of King-Smith and Vanderbilt
C ( Geometric Berry phase). Initialisation routine.
C Written by DSP, March 1999.
C **************************** INPUT **********************************
C real*8   ucell(3,3)         : Unit cell vectors
C integer  maxkpol            : Last dimension of kpoint 
C *************************** UNITS ***********************************
C Lengths in atomic units (Bohr).
C k vectors in reciprocal atomic units.
C *********************************************************************
      implicit          none
      integer           maxkpol
      double precision  ucell(3,3)
C *********************************************************************

C Internal variables 
      integer  ix, iy, kscell(3,3), igrd, nk, nmeshk(3,3)
         
      double precision  displ(3), dsp(3), cutoff

      external  timer

C Start time counter 
      call timer( 'KSV_init', 1 )

C Find the integration grids in reciprocal space   
C In the first call this is dome to set the rigth 
C dimension for the arrays which depend on parameter
C nkpol
      call repol(nmeshk,dsp) 
      
      nk=nmeshk(1,1)*nmeshk(2,1)*nmeshk(3,1)+
     .  nmeshk(1,2)*nmeshk(2,2)*nmeshk(3,2)+
     .  nmeshk(1,3)*nmeshk(2,3)*nmeshk(3,3) 

      if (nk.ne.0) then  
        nkpol=1
        do igrd=1,3
          cutoff=0.0d0
          do ix=1,3
            do iy= 1,3
              kscell(ix,iy)=0
            enddo
          enddo
          do ix=1,3
            if (ix.ne.igrd) then
              kscell(ix,ix)= nmeshk(ix,igrd)
              displ(ix)=dsp(igrd)
            else 
              kscell(ix,ix)=1
              displ(ix)=0.0d0
            endif
          enddo 
          nk=kscell(1,1)*kscell(2,2)*kscell(3,3) 
          if (nk.gt.0) then 
            call kgridinit( ucell, kscell, displ, cutoff, nk )
            if (maxkpol.gt.nk) then
              call kgrid( ucell, kscell,  displ,
     .                    nk, kpol, wgthpol )
            endif
          endif 
          nkpol=max(nk,nkpol)       
        enddo  
      else
        nkpol=0
      endif  

      call timer( 'KSV_init', 2 )

      end subroutine ksv_init
      
!----------------------------------------------------------------
      subroutine repol(nmeshk,displ)

c *******************************************************************
c Reading k_space grids for the calculation of the macroscopic 
c polarization
c
c Reads fdf block. All three grids must be specified, therefore 
c a 3x3 matrix of integer numbers must be given:
C    For example: 
C        %block PolarizationGrids
C              100 11   6  yes 
C              10  120  7  no 
C              12  15   85 yes
C        %endblock PolarizationGrids  
C
c   *  First grid will be used to calculate the polarization
c      along the direction of the first lattice vector.
c   *  Second grid  will be used for the calculation along the 
C      the direction of the second lattice vector.
C   *  Third grid  will be used for the calculation along the 
C      the direction of the third lattice vector. 
C   *  The last column specifies if the bidimensional grids
C      are going to be diplaced from the origin or not as in the
C      Monkhorst-Pack algorithm (PHYS.REV. B13, 5188 (1976)).
C      This last column is optional
C
C   The number in the diagonal of the matrix specified the number
C   of points to be used in the one dimensional line integrals along 
C   the different directions. The other number specified the mesh 
c   used in the surface integrals. For more details see, 
c   R.D. King-Smith, and D. Vanderbilt, PRB 47, 1651 (1993).
c   If the number of point in one of the grids is zero, the calculation
C   will not be done for this particular direction. 
c      
c Written by D. Sanchez-Portal. March 1999
c ********* OUTPUT **************************************************
c integer   nmeshk(3,3)     : Grids for the integration in reciprocal 
c                             space
c real*8    displ(3)        : Displacement of the bidimensional 
c                             integration grid from the origin
c *******************************************************************

C
C  Modules
C
      use precision
      use parallel,  only : Node
      use fdf
      use sys

      implicit          none

C Passed variables and arrays
      integer           nmeshk(3,3)
      real(dp)          displ(3)

C Internal variables and arrays
 
      character         yiorno
      integer           ni, nn, igrd, ix

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Set outputs to zero
      do igrd = 1,3
        do ix = 1,3
          nmeshk(ix,igrd) = 0
        enddo 
        displ(igrd) = 0.0d0
      enddo

C Check for block and read

      if ( fdf_block('PolarizationGrids',bfdf) ) then

        do igrd= 1, 3
          if (.not. fdf_bline(bfdf,pline)) then
            call die('repol: ERROR not enough lines in ' //
     &               'PolarizationGrids block')
          endif
          nn = fdf_bnnames(pline)
          ni = fdf_bnintegers(pline)
          if (ni .eq. 3) then
            nmeshk(1,igrd) = fdf_bintegers(pline,1)
            nmeshk(2,igrd) = fdf_bintegers(pline,2)
            nmeshk(3,igrd) = fdf_bintegers(pline,3)
               
            if (nn .ge. 1) then 
              yiorno = fdf_bnames(pline,1)
              if ( yiorno.eq.'yes'  .or.
     .             yiorno.eq.'YES'  .or.
     .             yiorno.eq.'y'    .or.
     .             yiorno.eq.'Y'    .or.
     .             yiorno.eq.'true' .or.
     .             yiorno.eq.'TRUE' .or.
     .             yiorno.eq.'t'    .or.
     .             yiorno.eq.'T' ) then 
                      
                displ(igrd) = 0.5d0
                           
              endif 
            endif         
          else
            call die('repol: Check syntax in PolarizationGrids block')
          endif

        enddo
        call fdf_bclose(bfdf)
      endif

      end subroutine repol

      end module m_ksvinit
