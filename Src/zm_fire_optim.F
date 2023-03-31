! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_zm_fire_optim
      public :: zm_fire_optimizer
      private

      CONTAINS
!
      subroutine zm_fire_optimizer( na, xa, fa, cell, stress, dxmax,
     .                              tp, ftol, strtol, varcel, relaxd )

c ******************************** INPUT ************************************
c integer na            : Number of atoms in the simulation cell
c real*8 fa(3,na)       : Atomic forces
c real*8 stress(3,3)    : Stress tensor components
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
c real*8 strtol         : Maximum stress tolerance
c logical varcel        : true if variable cell optimization
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
      use parallel,    only : Node, IONode
      use fdf,         only : fdf_get
      use units, only: Ang

      use m_fire
      use siesta_options, only: dt

      use zmatrix,     only : VaryZmat, Zmat
      use zmatrix,     only : CartesianForce_to_ZmatForce
      use zmatrix,     only : ZmatForce,ZmatForceVar
      use zmatrix,     only : iZmattoVars,ZmatType
      use zmatrix,     only : Zmat_to_Cartesian
      use zmatrix,     only : coeffA, coeffB, iNextDept
      use Zmatrix,     only : ZmatForceTolLen, ZmatForceTolAng
      use Zmatrix,     only : ZmatMaxDisplLen, ZmatMaxDisplAng

      implicit none

C Subroutine arguments:
      integer,  intent(in)   :: na
      logical,  intent(in)   :: varcel
      logical,  intent(out)  :: relaxd
      real(dp), intent(in)   :: fa(3,na), dxmax, tp, ftol
      real(dp), intent(inout):: xa(3,na), stress(3,3), strtol, cell(3,3)

C Internal variables and arrays
      integer             :: ia, i, j, n, indi,indi1,vi,k
      real(dp)            :: volume, new_volume, trace
      real(dp)            :: celli(3,3)
      real(dp)            ::  stress_dif(3,3)
      real(dp), dimension(:), allocatable :: gxa, gfa
      real(dp)            :: force, force1

C Saved internal variables:
      integer,       save :: ndeg
      logical,       save :: frstme = .true.
      logical,       save :: constant_volume
      real(dp),      save :: initial_volume

      real(dp),      save :: tstres(3,3) 
      real(dp),      save :: cellin(3,3) 
      real(dp),      save :: modcel(3) 
      real(dp),      save :: precon 
      real(dp),      save :: strain(3,3)  ! Special treatment !!

      real(dp), dimension(:), allocatable    :: ftoln, dxmaxn

      type(fire_t), save  :: b

      integer        :: ndegi, ndi

      logical, save  :: fire_debug
      real(dp), save :: fire_mass
      real(dp)       :: fire_dt, fire_dt_inc,
     $                  fire_dt_dec, fire_alphamax,
     $                  fire_alpha_factor, fire_dtmax
      integer        :: fire_nmin

      external          memory
      real(dp), external :: volcel

c ---------------------------------------------------------------------------

      volume = volcel(cell)

      if ( frstme ) then
         
        fire_mass = fdf_get("MD.FIRE.Mass", 1.0_dp)
        fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
        fire_dt_inc = fdf_get("MD.FIRE.TimeInc", FIRE_DEF_dt_inc)
        fire_dt_dec = fdf_get("MD.FIRE.TimeDec", FIRE_DEF_dt_dec)
        fire_nmin = fdf_get("MD.FIRE.Nmin", FIRE_DEF_nmin)
        fire_alphamax = fdf_get("MD.FIRE.AlphaMax",
     $       FIRE_DEF_alphamax)
        fire_alpha_factor = fdf_get("MD.FIRE.AlphaFactor",
     &       FIRE_DEF_alpha_factor)
        fire_dtmax = fdf_get("MD.FIRE.MaxTimeStep",
     $       FIRE_DEF_dtmax)
        fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
        fire_debug = fdf_get("MD.FIRE.Debug", .false.)

        ! Find number of variables
        ndeg = 0
        do ia = 1,na
           do n = 1,3
              indi = 3*(ia-1) + n
              if (VaryZmat(indi)) then
                 ndeg = ndeg + 1
              endif
           enddo
        enddo
        if ( varcel ) then
           ndeg = ndeg + 6      ! Including stress
        endif

        if (Ionode) then
           write(6,'(a,i6)') "Fire_optim: No of elements: ", ndeg
        endif
         call fire_setup(b, n=ndeg, dt=fire_dt,
     $                   debug=fire_debug,
     $                   dt_inc=fire_dt_inc, dt_dec=fire_dt_dec,
     $                   alphamax=fire_alphamax,
     $                   alpha_factor=fire_alpha_factor,
     $                   nmin=fire_nmin)


        if ( varcel ) then

          constant_volume = fdf_get("MD.ConstantVolume", .false.)

          call get_target_stress(tp,tstres)

            if (constant_volume) then
               tstres(:,:) = 0.0_dp
               write(6,"(a)") "***Target stress set to zero " //
     $              "for constant-volume calculation"
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
C ---- (a length in Bohrs typical of bond lengths ..) 

          ! AG: This could better be volume^(1/3) by default
          precon = fdf_get('MD.PreconditionVariableCell',
     .                           9.4486344d0,'Bohr')

C Initialize absolute strain and save initial cell vectors
C Initialization to 1 for numerical reasons, later substracted
       
          strain(1:3,1:3) = 1.0_dp
          cellin(1:3,1:3) = cell(1:3,1:3)
          initial_volume = volcel(cellin)

        endif     ! varcel


        frstme = .false.
      endif                 ! First time

C Allocate local memory

      allocate(gfa(ndeg))
      call memory('A','D',ndeg,'cgvc_zmatrix')
      allocate(gxa(ndeg))
      call memory('A','D',ndeg,'cgvc_zmatrix')
      allocate(ftoln(ndeg))
      call memory('A','D',ndeg,'cgvc_zmatrix')
      allocate(dxmaxn(ndeg))
      call memory('A','D',ndeg,'cgvc_zmatrix')

      if ( varcel ) then

        ! Inverse of matrix of cell vectors  (transpose of)
        call reclat( cell, celli, 0 )

C Transform coordinates and forces to fractional ---------------------------- 
C but scale them again to Bohr by using the (fix) moduli of the original ----
C lattice vectors (allows using maximum displacement as before) -------------
C convergence is checked here for input forces as compared with ftol --------

        relaxd = .true.
        ndegi = 0
        ! Loop over degrees of freedom, scaling 
        ! only cartesian coordinates to fractional
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
            dxmaxn(ndegi) = ZmatMaxDisplLen
          enddo
        enddo

C Check stress convergence
        strtol = dabs(strtol)
        do i = 1,3
          do j = 1,3
            relaxd = relaxd .and. 
     .        ( dabs(stress_dif(i,j)) .lt. abs(strtol) )
          enddo
        enddo

      else   ! FIXED CELL

C Set number of degrees of freedom & variables
         relaxd = .true.
        ndegi = 0
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
                relaxd=relaxd.and.(dabs(force).lt.ftoln(ndegi))
              endif
            enddo
          enddo

      endif

      if (relaxd) RETURN   ! Will leave work arrays allocated

      call fire_step(b,gfa,gxa,dxmaxn)


      ! Transform back
      if ( varcel) then

        ! New cell
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
                      ndi = ndegi + n - k
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
C Fixed cell
C Return coordinates to correct arrays 
        ndegi = 0
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
      endif


C Deallocate local memory

      call memory('D','D',size(dxmaxn),'cgvc_zmatrix')
      deallocate(dxmaxn)
      call memory('D','D',size(ftoln),'cgvc_zmatrix')
      deallocate(ftoln)
      call memory('D','D',size(gxa),'cgvc_zmatrix')
      deallocate(gxa)
      call memory('D','D',size(gfa),'cgvc_zmatrix')
      deallocate(gfa)

      end subroutine zm_fire_optimizer

      end module m_zm_fire_optim
