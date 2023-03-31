! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine cgvc( na, xa, fa, cell, stress, dxmax, 
     .                 tp, ftol, strtol, varcel, relaxd, usesavecg )
c ***************************************************************************
c Variable-cell conjugate-gradient geometry optimization
c
c   Energy minimization including atomic coordinates and unit cell vectors.
c   It allows an external target stress:
c              %block MD.TargetStress
c                  3.5  0.0  0.0  0.0  0.0  0.0
c              %endblock MD.TargetStress
c   corresponding to xx, yy, zz, xy, xz, yz.
c   In units of (-MD.TargetPressure)
c   Default: hydrostatic pressure: -1, -1, -1, 0, 0, 0
c
c   Based on E({xa},stress), with {xa} in fractional coor
c   The gradient of the energy given by {cfa} forces (fractional!)
c   The fractional coordinates are multiplied by the initial cell vectors
c   to get them in Bohr for dxmax and preconditioning.
c      
c Written by E. Artacho. November 1998. 
c ******************************** INPUT ************************************
c integer na            : Number of atoms in the simulation cell
c real*8 fa(3,na)       : Atomic forces
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
c logical varcel        : true if variable cell optimization
c logical usesavecg     : true if we're using saved CG files.
c *************************** INPUT AND OUTPUT ******************************
c real*8 xa(3,na)       : Atomic coordinates
c                         input: current step; output: next step
c real*8 cell(3,3)      : Matrix of the vectors defining the unit cell 
c                         input: current step; output: next step
c                         cell(ix,ja) is the ix-th component of ja-th vector
c real*8 stress(3,3)    : Stress tensor components
c real*8 strtol         : Maximum stress tolerance
c ******************************** OUTPUT ***********************************
c logical relaxd        : True when converged
c ***************************************************************************

C
C  Modules
C
      use precision, only : dp
      use parallel,    only : Node
      use sys,         only : die
      use fdf

      use m_mpi_utils, only:  broadcast
      use units, only: kBar, Ang
!
!     Use old-style conjugate-gradient routine
!     Use the new one only for the Zmatrix case
!
      use m_conjgr_old, only: conjgr
      use alloc,        only: re_alloc, de_alloc

      implicit none

! Subroutine arguments:

      integer, intent(in) :: na
      real(dp), intent(in) :: fa(3,na), dxmax,
     .                        tp, ftol, strtol
      logical, intent(in) :: varcel, usesavecg
      real(dp), intent(inout) :: xa(3,na), stress(3,3), cell(3,3)
      logical, intent(out) :: relaxd

c Internal variables and arrays

      real(dp)            :: new_volume, trace, ftol_tmp, volume

      logical           found
      integer           ia, i, j, n, indi

      real(dp) ::  celli(3,3), sxx, syy, szz, sxy, sxz, syz
      real(dp) ::  stress_dif(3,3)

      real(dp), pointer :: gxa(:)=>null(), gfa(:)=>null()
      real(dp), pointer, save :: cgaux(:)=>null()

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline=>null()

! Saved internal variables:

      logical, save :: frstme = .true.
      logical, save :: tarstr = .false.
      logical, save :: constant_volume
      real(dp), save :: initial_volume

      real(dp), save :: cgcntr(0:20) = 0.0_dp

      integer, save :: ndeg,
     .                 linmin

      real(dp), save :: tstres(3,3),
     .                  modcel(3),
     .                  precon,
     .                  strain(3,3),
     .                  cellin(3,3)
      
      real(dp) :: volcel
      external :: volcel
c ---------------------------------------------------------------------------

      volume = volcel(cell)

C Allocate local memory
      call re_alloc( gfa, 1, na*3+6, 'gfa', 'cgvc' )
      call re_alloc( gxa, 1, na*3+6, 'gxa', 'cgvc' )
      if (.not.associated(cgaux))
     &  call re_alloc( cgaux, 1, na*6+12, 'cgaux', 'cgvc' )

C If first call to cgvc, check dim and get target stress --------------------

      if ( frstme ) then
  
        if ( varcel ) then

C Check if we want a constant-volume simulation
          constant_volume = fdf_get("MD.ConstantVolume", .false.)

C Look for target stress and read it if found, otherwise generate it --------

          tarstr = fdf_block('MD.TargetStress',bfdf)

          if (tarstr) then
            if (Node.eq.0) then
              write(6,'(/a,a)') 'cgvc: Reading %block MD.TargetStress',
     .                          ' (units of MD.TargetPressure).'
            endif
            if (.not. fdf_bline(bfdf,pline))
     .        call die('cgvc: ERROR in MD.TargetStress block')
            sxx = fdf_bvalues(pline,1)
            syy = fdf_bvalues(pline,2)
            szz = fdf_bvalues(pline,3)
            sxy = fdf_bvalues(pline,4)
            sxz = fdf_bvalues(pline,5)
            syz = fdf_bvalues(pline,6)
            call fdf_bclose(bfdf)
            
            tstres(1,1) = - sxx * tp
            tstres(2,2) = - syy * tp
            tstres(3,3) = - szz * tp
            tstres(1,2) = - sxy * tp
            tstres(2,1) = - sxy * tp
            tstres(1,3) = - sxz * tp
            tstres(3,1) = - sxz * tp
            tstres(2,3) = - syz * tp
            tstres(3,2) = - syz * tp
          else
            if (Node.eq.0) then
              write(6,'(/a,a)') 'cgvc: No target stress found, ',
     .             'assuming hydrostatic MD.TargetPressure.'
            endif
            do i= 1, 3
              do j= 1, 3
                tstres(i,j) = 0._dp
              enddo
              tstres(i,i) = - tp
            enddo
          endif

C Write target stress down --------------------------------------------------

          if (constant_volume) then
            tstres(:,:) = 0.0_dp
            if (Node.eq.0) then
              write(6,"(a)") "***Target stress set to zero " //
     $            "for constant-volume calculation"
            endif
          endif
          if (Node.eq.0) then
            write(6,"(/a)") 'cgvc: Target stress (kBar)'
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(1,1)/kBar, tstres(1,2)/kBar,
     .       tstres(1,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(2,1)/kBar, tstres(2,2)/kBar,
     .       tstres(2,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(3,1)/kBar, tstres(3,2)/kBar,
     .       tstres(3,3)/kBar
          endif


C Moduli of original cell vectors for fractional coor scaling back to au ---

          do n = 1, 3
            modcel(n) = 0.0_dp
            do j = 1, 3
              modcel(n) = modcel(n) + cell(j,n)*cell(j,n)
            enddo
            modcel(n) = dsqrt( modcel(n) )
          enddo

C Scale factor for strain variables to share magnitude with coordinates 
C ---- (a length in Bohrs typical of bond lengths ..) ------------------

          precon = fdf_get('MD.PreconditionVariableCell',
     .                     9.4486344d0,'Bohr')

C Dimension of space where E is minimized ------------------------------

          ndeg = na*3 + 6

C Initialize absolute strain and save initial cell vectors -------------
C Initialization to 1. for numerical reasons, later substracted --------

          strain(1:3,1:3) = 1.0_dp
          cellin(1:3,1:3) = cell(1:3,1:3)
          initial_volume = volcel(cellin)

        else

          ndeg = na*3

        endif

C Initialize and read cgaux and cgcntr if present and wanted ---------------

        if (usesavecg) then
          if (Node.eq.0) then
            call iocg( 'read', ndeg*2, cgaux, cgcntr, relaxd, found )
          endif
          call broadcast(cgaux)
          call broadcast(cgcntr)
          call broadcast(relaxd)
          call broadcast(found)

          if ( found ) then
            linmin = cgcntr(1)
! Simple bugfix for fixed cell restarts. In the case of a tightend ftol
! after a sucessful CG run, relaxed is set to .true. by the read from
! iocg. This results in the test for convergence being skipped within
! conjgr. Setting relaxed to .false. here avoids this. In the case of 
! a variable cell relaxed is set below anyway and the converence test
! is always performed. A better solution (moving the test out into a
! seperate routine) will follow in a later version.  -- AMW 8/7/2008
            relaxd = .false.
          else
            if (Node.eq.0) then
              write(6,'(/,a)') 'cgvc: WARNING: CG file not found'
            endif
            relaxd = .false.
            cgcntr(0) = 0
            linmin = 1
          endif
        else
          relaxd = .false.
          cgcntr(0) = 0
          linmin = 1
        endif

        frstme = .false.
      endif

C Variable cell -------------------------------------------------------------

      if ( varcel ) then

C Inverse of matrix of cell vectors  (transpose of) ------------------------

        call reclat( cell, celli, 0 )

C Transform coordinates and forces to fractional ---------------------------- 
C but scale them again to Bohr by using the (fix) moduli of the original ----
C lattice vectors (allows using maximum displacement as before) -------------
C convergence is checked here for input forces as compared with ftol --------

        relaxd = .true.
        do ia = 1, na
          do n = 1, 3
            indi = 3*(ia - 1) + n
            gxa(indi) = 0.0_dp
            gfa(indi) = 0.0_dp
            relaxd = relaxd .and. ( abs(fa(n,ia)) .lt. ftol )
            do i = 1, 3
              gxa(indi) = gxa(indi) + celli(i,n) * xa(i,ia) * modcel(n)
              gfa(indi) = gfa(indi) + cell(i,n) * fa(i,ia) / modcel(n)
            enddo
          enddo
        enddo

C Symmetrizing the stress tensor --------------------------------------------

        do i = 1, 3
           do j = i+1, 3
              stress(i,j) = 0.5_dp*( stress(i,j) + stress(j,i) )
              stress(j,i) = stress(i,j)
           enddo
        enddo

C Subtract target stress

        stress_dif = stress - tstres
!
!       Take 1/3 of the trace out here if constant-volume needed
!
        if (constant_volume) then
           trace = stress_dif(1,1) + stress_dif(2,2) + stress_dif(3,3)
           do i=1,3
              stress_dif(i,i) = stress_dif(i,i) - trace/3.0_dp
           enddo
        endif

C Append excess stress and strain to gxa and gfa ------ 
C preconditioning: scale stress and strain so as to behave similarly to x,f -

        indi = 3*na
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              gfa(indi) = -stress_dif(i,j)*volume/precon
              gxa(indi) = strain(i,j) * precon
           enddo
        enddo

C Check stress convergence --------------------------------------------------

        do i = 1, 3
           do j = 1, 3
              relaxd = relaxd .and. 
     .          ( abs(stress_dif(i,j)) .lt. abs(strtol) )
           enddo
        enddo

C Call conjugate gradient minimization -------------------------------------- 

        if ( .not. relaxd ) then
           ftol_tmp = 0.0_dp
!!           do i=1,ndeg 
!!              print *, "gxa, gfa ", i, gxa(i), gfa(i)
!!           enddo
           call conjgr(ndeg,gxa,gfa,dxmax, ftol_tmp ,cgcntr,cgaux )
        endif


C Fixed cell ----------------------------------------------------------------

      else

         if (.not. relaxd) then
!!            do i=1,ndeg 
!!               print *, "gxa, gfa ", i, gxa(i), gfa(i)
!!            enddo
            call conjgr( 3*na, xa, fa, dxmax, ftol, cgcntr, cgaux )
            relaxd = (int(cgcntr(0)) .eq. 0)           !! ???
         endif

      endif

C Checking line minimizations and convergence -------------------------------

      if (nint(cgcntr(1)) .ne. linmin) then
        if (Node.eq.0) then
          write(6,'(/a,i4,a,f10.4)')
     .      'cgvc: Finished line minimization ', linmin,
     .      '.  Mean atomic displacement =', cgcntr(18)/sqrt(dble(na))
        endif
        linmin = nint(cgcntr(1))
      endif

C Transform back if variable cell ------------------------------------------- 

      if ( varcel .and. (.not. relaxd) ) then

C New cell ------------------------------------------------------------------

        indi = 3*na
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              strain(i,j) = gxa(indi) / precon
              strain(j,i) = strain(i,j)
           enddo
        enddo

        cell = cellin + matmul(strain-1.0_dp,cellin)
        if (constant_volume) then
           new_volume = volcel(cell)
           if (Node.eq.0) write(6,"(a,f12.4)")
     $          "Volume before coercion: ",  new_volume/Ang**3
           cell =  cell * (initial_volume/new_volume)**(1.0_dp/3.0_dp)
        endif

C Output fractional coordinates to cartesian Bohr, and copy to xa ----------- 

        do ia = 1, na
          do i = 1, 3
            xa(i,ia) = 0.0_dp
            do n = 1, 3
              indi = 3*(ia - 1) + n
              xa(i,ia) = xa(i,ia) + cell(i,n) * gxa(indi) / modcel(n)
            enddo
          enddo
        enddo

      endif

C Save cgaux ----------------------------------------------------------------

      if (Node.eq.0) 
     .  call iocg( 'write', ndeg*2, cgaux, cgcntr, relaxd, found )

C Deallocate local memory
      call de_alloc( gxa, 'gxa', 'cgvc' )
      call de_alloc( gfa, 'gfa', 'cgvc' )

      return
      end

