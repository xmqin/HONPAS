! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine cgvc_zmatrix( na, xa, fa, cell, stress, dxmax, 
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
c real*8 stress(3,3)    : Stress tensor components
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
c real*8 strtol         : Maximum stress tolerance
c logical varcel        : true if variable cell optimization
c logical usesavecg     : true if we're using saved CG files
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
      use precision,   only : dp
      use parallel,    only : Node
      use sys,         only : die
      use fdf
      use m_mpi_utils, only:  broadcast
      use units, only: kBar, Ang

      use zmatrix,     only : VaryZmat, lUseZmatrix, Zmat
      use zmatrix,     only : CartesianForce_to_ZmatForce
      use zmatrix,     only : ZmatForce,ZmatForceVar
      use zmatrix,     only : iZmattoVars,ZmatType
      use zmatrix,     only : Zmat_to_Cartesian
      use zmatrix,     only : coeffA, coeffB, iNextDept
      use Zmatrix,     only : ZmatForceTolLen, ZmatForceTolAng
      use Zmatrix,     only : ZmatMaxDisplLen, ZmatMaxDisplAng
      use m_conjgr,    only : conjgr
      use alloc,       only : re_alloc, de_alloc

      implicit none

C Subroutine arguments:
      integer,  intent(in)   :: na
      logical,  intent(in)   :: varcel, usesavecg
      logical,  intent(out)  :: relaxd
      real(dp), intent(in)   :: fa(3,na), dxmax, tp, ftol
      real(dp), intent(inout):: xa(3,na), stress(3,3), strtol, cell(3,3)

C Internal variables and arrays
      integer             :: ia, i, j, n, indi,indi1,vi,k,ndegi,ndi
      logical             :: found
      real(dp)            :: volume, new_volume, trace
      real(dp)            :: sxx, syy, szz, sxy, sxz, syz
      real(dp)            :: celli(3,3)
      real(dp)            :: stress_dif(3,3)
      real(dp), pointer   :: gxa(:) =>null(), gfa(:)=>null()
      real(dp)            :: force, force1

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

C Saved internal variables:
      integer,       save :: ndeg, linmin
      logical,       save :: frstme = .true.
      logical,       save :: tarstr = .false.
      logical,       save :: constant_volume
      real(dp),      save :: initial_volume

      real(dp),      save :: cgcntr(0:20) = 0.0_dp
      real(dp),      save :: tstres(3,3) 
      real(dp),      save :: modcel(3) 
      real(dp),      save :: precon 
      real(dp),      save :: strain(3,3)
      real(dp),      save :: cellin(3,3)

      real(dp), pointer, save :: cgaux(:)
      real(dp), pointer       :: ftoln(:)=>null(), dxmaxn(:)=>null(),
     &                           ftoln_tmp(:)=>null()

      real(dp), external :: volcel

c ---------------------------------------------------------------------------

      volume = volcel(cell)

C If first call to cgvc_zmatrix, check dim and get target stress --------------------

      if ( frstme ) then
  
C Look for target stress and read it if found, otherwise generate it --------

        if ( varcel ) then
     
C Check if we want a constant-volume simulation
          constant_volume = fdf_get("MD.ConstantVolume", .false.)

          tarstr = fdf_block('MD.TargetStress',bfdf)

          if (tarstr) then
            if (Node.eq.0) then
              write(6,'(/a,a)')
     $            'cgvc_zmatrix: Reading %block MD.TargetStress',
     .                     ' (units of MD.TargetPressure).'
            endif
            if (.not. fdf_bline(bfdf,pline))
     $        call die('cgvc_zmatrix: ERROR in MD.TargetStress block')
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
              write(6,'(/a,a)')
     $            'cgvc_zmatrix: No target stress found, ',
     .            'assuming hydrostatic MD.TargetPressure.'
            endif
            do i= 1, 3
              do j= 1, 3
                tstres(i,j) = 0.0_dp
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
            write(6,"(/a)") 'cgvc_zmatrix: Target stress (kBar)'
            do i = 1, 3
              write(6,"(a,2x,3f12.3)") 
     .             'cgvc_zmatrix:', (tstres(1,j)/kBar,j=1,3)
            enddo
          endif

C Moduli of original cell vectors for fractional coor scaling back to au ---
          do n= 1, 3
            modcel(n) = 0.0_dp
            do j= 1, 3
              modcel(n) = modcel(n) + cell(j,n)*cell(j,n)
            enddo
            modcel(n) = dsqrt( modcel(n) )
          enddo

C Scale factor for strain variables to share magnitude with coordinates -----
C ---- (a length in Bohrs typical of bond lengths ..) -----------------------

          precon = fdf_get('MD.PreconditionVariableCell',
     .                     9.4486344d0,'Bohr')

C Initialize absolute strain and save initial cell vectors -----------------
C Initialization to 1. for numerical reasons, later substracted ------------
       
          strain(1:3,1:3) = 1.0_dp
          cellin(1:3,1:3) = cell(1:3,1:3)
          initial_volume = volcel(cellin)

        endif             ! varcel

        ! Find number of degrees of freedom
        if (lUseZmatrix) then
           ndeg = 0
           do ia = 1,na
              do n = 1,3
                 indi = 3*(ia-1) + n
                 if (VaryZmat(indi)) then
                    ndeg = ndeg + 1
                 endif
              enddo
           enddo
        else
           ndeg = 3*na
        endif
        if (varcel) then
           ndeg = ndeg + 6
        endif

C Initialize and read cgaux and cgcntr if present and wanted ---------------

        if (.not.associated(cgaux))
     &     call re_alloc( cgaux, 1, 2*ndeg, 'cgaux', 'cgvc_zmatrix' )

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
          else
            if (Node.eq.0) then
              write(6,'(/,a)')
     $              'cgvc_zmatrix: WARNING: CG file not found'
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

      endif      ! First time

C Allocate local memory
      call re_alloc( gfa, 1, ndeg, 'gfa', 'cgvc_zmatrix' )
      call re_alloc( gxa, 1, ndeg, 'gxa', 'cgvc_zmatrix' )
      call re_alloc( ftoln, 1, ndeg, 'ftoln', 'cgvc_zmatrix' )
      call re_alloc( ftoln_tmp, 1, ndeg, 'ftoln_tmp', 'cgvc_zmatrix' )
      call re_alloc( dxmaxn, 1, ndeg, 'dxmaxn', 'cgvc_zmatrix' )

C Variable cell -------------------------------------------------------------

      if ( varcel ) then

C Inverse of matrix of cell vectors  (transverse of) ------------------------

        call reclat( cell, celli, 0 )

C Transform coordinates and forces to fractional ---------------------------- 
C but scale them again to Bohr by using the (fix) moduli of the original ----
C lattice vectors (allows using maximum displacement as before) -------------
C convergence is checked here for input forces as compared with ftol --------

        relaxd = .true.
        ndegi = 0
        if (lUseZmatrix) then
C Loop over degrees of freedom, scaling cartesian coordinates to fractional
          do ia = 1,na
            do n = 1,3
              indi = 3*(ia-1) + n
              if (VaryZmat(indi)) then
                ndegi = ndegi + 1
                if (ZmatType(indi).gt.1) then
                  ftoln(ndegi) = ZmatForceTolLen
                  dxmaxn(ndegi) = ZmatMaxDisplLen
                else
                  ftoln(ndegi) = ZmatForceTolAng
                  dxmaxn(ndegi) = ZmatMaxDisplAng
                endif
                vi = iZmattoVars(indi)
                if (vi.eq.0) then
                  force = ZmatForce(indi)
                else
                  force = ZmatForceVar(vi)
                endif 
                relaxd=relaxd.and.(dabs(force).lt.ftoln(ndegi))
                if (ZmatType(indi).gt.2) then
C Cartesian coordinate                
                  gxa(ndegi) = 0.0_dp
                  gfa(ndegi) = 0.0_dp
                  do i = 1,3
                    indi1 = 3*(ia-1)+i
                    gxa(ndegi) = gxa(ndegi)+
     .                          celli(i,n)*Zmat(indi1)*modcel(n)
                    if (i.eq.n) then
                      force1 = force
                    else
                      force1 = ZmatForce(indi1)
                    endif
                    gfa(ndegi) = gfa(ndegi)+ 
     .                          cell(i,n)*force1/modcel(n)
                  enddo
                else
                  gxa(ndegi) = Zmat(indi)
                  gfa(ndegi) = force
                endif
              endif
            enddo
          enddo
        else
C No Z-Matrix
          do ia = 1,na
            do n = 1,3
              ndegi = ndegi + 1
              ftoln(ndegi) = ftol
              dxmaxn(ndegi) = dxmax
              gxa(ndegi) = 0.0_dp
              gfa(ndegi) = 0.0_dp
              relaxd = relaxd .and. ( dabs(fa(n,ia)) .lt. ftoln(ndegi) )
              do i = 1, 3
                gxa(ndegi) = gxa(ndegi)+celli(i,n)*xa(i,ia)*modcel(n)
                gfa(ndegi) = gfa(ndegi)+ cell(i,n)*fa(i,ia)/modcel(n)
              enddo
            enddo
          enddo
        endif


C Symmetrizing the stress tensor --------------------------------------------
        do i = 1,3
          do j = i+1,3
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

C Append stress (substracting target stress) and strain to gxa and gfa ------ 
C preconditioning: scale stress and strain so as to behave similarly to x,f -
        do i = 1,3
          do j = i,3
            ndegi = ndegi + 1
            gfa(ndegi) = -stress_dif(i,j)*volume/precon
            gxa(ndegi) = strain(i,j) * precon
            if (lUseZmatrix) then
              dxmaxn(ndegi) = ZmatMaxDisplLen
            else
              dxmaxn(ndegi) = dxmax
            endif
          enddo
        enddo

C  testing
!!          print *, "ndeg, ndegi:", ndeg, ndegi
        
C Check stress convergence --------------------------------------------------
        strtol = dabs(strtol)
        do i = 1,3
          do j = 1,3
            relaxd = relaxd .and. 
     .        ( dabs(stress_dif(i,j)) .lt. abs(strtol) )
          enddo
        enddo

C Call conjugate gradient minimization -------------------------------------- 
        if ( .not. relaxd ) then
          do i=1,ndeg 
            ftoln_tmp(i) = 0.d0
!!            print *, "gxa, gfa ", i, gxa(i), gfa(i)
          enddo
          call conjgr(ndeg,gxa,gfa,dxmaxn,ftoln_tmp,cgcntr,cgaux)
        endif

      else     ! Fixed cell
C Set number of degrees of freedom & variables
        relaxd = .true.
        ndegi = 0
        if (lUseZmatrix) then
          do i = 1,na
            do k = 1,3
              indi = 3*(i-1)+k
              if (VaryZmat(indi)) then
                ndegi = ndegi + 1
                gxa(ndegi) = Zmat(indi)
                vi = iZmattoVars(indi)
                if (vi.eq.0) then
                  force = ZmatForce(indi)
                else
                  force = ZmatForceVar(vi)
                endif
                gfa(ndegi) = force
                if (ZmatType(indi).eq.1) then
                  ftoln(ndegi) = ZmatForceTolAng
                  dxmaxn(ndegi) = ZmatMaxDisplAng
                else
                  ftoln(ndegi) = ZmatForceTolLen
                  dxmaxn(ndegi) = ZmatMaxDisplLen
                endif
                relaxd=relaxd.and.(abs(force).lt.ftoln(ndegi))
              endif
            enddo
          enddo
        else
          do i = 1,na
            do k = 1,3
              ndegi = ndegi + 1
              gxa(ndegi) = xa(k,i)
              gfa(ndegi) = fa(k,i)
              dxmaxn(ndegi) = dxmax
              ftoln(ndegi) = ftol
              relaxd = relaxd .and. ( dabs(fa(k,i)) .lt. ftoln(ndegi) )
            enddo
          enddo
        endif

        if (.not. relaxd) then
!!           do i=1,ndeg 
!!              print *, "gxa, gfa ", i, gxa(i), gfa(i)
!!           enddo
           call conjgr( ndeg, gxa, gfa, dxmaxn, ftoln, cgcntr, cgaux)
           relaxd = (int(cgcntr(0)) .eq. 0)
        endif

      endif
C End of fixed cell

C Checking line minimizations and convergence -------------------------------

      if (nint(cgcntr(1)) .ne. linmin) then
        if (Node.eq.0) then
          write(6,'(/a,i4,a,f10.4)')
     .      'cgvc_zmatrix: Finished line minimization ', linmin,
     .      '.  Mean atomic displacement =', cgcntr(18)/sqrt(dble(na))
        endif
        linmin = nint(cgcntr(1))
      endif

C Transform back if variable cell ------------------------------------------- 
      if ( varcel .and. (.not. relaxd)) then

C New cell ------------------------------------------------------------------
        indi = ndeg-6
        do i = 1,3
          do j = i,3
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
        ndegi = 0
        if (lUseZmatrix) then
          do ia=1,na
            do k=1,3
              indi = 3*(ia-1)+k
              if (VaryZmat(indi)) then
                ndegi = ndegi + 1
                j = indi
                do while (j.ne.0) 
                  if (ZmatType(j).gt.2) then
                    Zmat(j) = 0.0_dp
                    do n=1,3
                      indi1 = 3*(ia-1)+n
                      ! Assume all three coords of this atom
                      ! are variables
                      ndi = ndegi + n - k   ! correct index
                      Zmat(j)=Zmat(j)+cell(k,n)*gxa(ndi)/modcel(n)
                    enddo
                  else
                    Zmat(j) = gxa(ndegi)
                  endif
                  Zmat(j) = Zmat(j)*coeffA(j) + coeffB(j)
                  j = iNextDept(j)
                enddo
              endif
            enddo
          enddo  
          call Zmat_to_Cartesian(xa)
        else
          do ia = 1,na
            do k = 1,3
              ndegi = ndegi+1
              xa(k,ia) = 0.0_dp
              do n = 1,3
                indi1 = 3*(ia-1) + n
                ndi = ndegi + n - k ! correct index == indi1
                xa(k,ia) = xa(k,ia) + cell(k,n) * gxa(ndi) / modcel(n)
              enddo
            enddo
          enddo
        endif

      else    ! fixed cell

        ! Return coordinates to correct arrays 
        ndegi = 0
        if (lUseZmatrix) then
          do ia = 1,na
            do k = 1,3
              indi = 3*(ia-1)+k
              if (VaryZmat(indi)) then
                ndegi = ndegi + 1
                j = indi
                do while (j.ne.0)
                  Zmat(j) = gxa(ndegi)*coeffA(j)+coeffB(j)
                  j = iNextDept(j)
                enddo
              endif
            enddo
          enddo
          call Zmat_to_Cartesian(xa)
        else
          do ia = 1,na
            do k = 1,3
              ndegi = ndegi + 1
              xa(k,ia) = gxa(ndegi)
            enddo
          enddo
        endif

      endif

C Save cgaux ----------------------------------------------------------------

      if (Node.eq.0) 
     .  call iocg( 'write', ndeg*2, cgaux, cgcntr, relaxd, found )

C Deallocate local memory
      call de_alloc( dxmaxn, 'dxmaxn', 'cgvc_zmatrix' )
      call de_alloc( ftoln_tmp, 'ftoln_tmp', 'cgvc_zmatrix' )
      call de_alloc( ftoln, 'ftoln', 'cgvc_zmatrix' )
      call de_alloc( gxa, 'gxa', 'cgvc_zmatrix' )
      call de_alloc( gfa, 'gfa', 'cgvc_zmatrix' )

      return
      end
