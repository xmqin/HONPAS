! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_fire_optim
!
!     FIRE geometry optimization
!      
      use precision, only : dp
      use fdf,          only: fdf_block
      use m_mpi_utils, only:  broadcast
      use units, only: kBar, eV, Ang
      use alloc, only: re_alloc, de_alloc
      use fdf, only: fdf_get
      
      use m_fire
      use siesta_options, only: dt
      use parallel, only : ionode

      implicit none

      public :: fire_optimizer
      private

      CONTAINS
      subroutine fire_optimizer( na, xa, fa, cell, stress,
     $     dxmax, tp, ftol, strtol, varcel, relaxd)
c ***************************************************************************
c FIRE geometry optimization
c
c ******************************** INPUT ************************************
c integer na            : Number of atoms in the simulation cell
c real*8 fa(3,na)       : Atomic forces
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
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

! Subroutine arguments:

      integer, intent(in) :: na
      real(dp), intent(in) :: fa(3,na), dxmax,
     .                        tp, ftol, strtol
      logical, intent(in) :: varcel
      real(dp), intent(inout) :: xa(3,na), stress(3,3), cell(3,3)
      logical, intent(out) :: relaxd

c Internal variables and arrays

      real(dp)            :: new_volume, trace, volume

      integer           ia, i, j, n, indi

      real(dp) ::  celli(3,3)
      real(dp) ::  stress_dif(3,3)

      real(dp), dimension(:), pointer       :: gxa, gfa, deltamax
      real(dp), dimension(:), pointer       :: rold, rdiff

! Saved internal variables:

      logical, save :: frstme = .true.
      logical, save :: constant_volume
      real(dp), save :: initial_volume


      real(dp), save :: tstres(3,3),
     .                  modcel(3),
     .                  precon,
     .                  strain(3,3),
     .                  cellin(3,3)

      type(fire_t), save  :: b

      integer, save  :: numel

      logical, save  :: fire_debug
      real(dp), save :: fire_mass
      real(dp)       :: fire_dt, fire_dt_inc,
     $                  fire_dt_dec, fire_alphamax,
     $                  fire_alpha_factor, fire_dtmax
      integer        :: fire_nmin

      real(dp), parameter ::  fovermp = 0.009579038 * Ang**2 / eV
      real(dp), parameter ::  dstrain_max = 0.1_dp

      real(dp) :: volcel
      external :: volcel
c ---------------------------------------------------------------------------

      volume = volcel(cell)

C If first call to optim, check dim and get target stress --------------------

      if ( frstme ) then
         
         fire_mass = fdf_get("MD.FIRE.Mass", 1.0_dp)
         fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
         fire_dt_inc = fdf_get("MD.FIRE.TimeInc", FIRE_DEF_dt_inc)
         fire_dt_dec = fdf_get("MD.FIRE.TimeDec", FIRE_DEF_dt_dec)
         fire_nmin = fdf_get("MD.FIRE.Nmin", FIRE_DEF_nmin)
         fire_alphamax = fdf_get("MD.FIRE.AlphaMax",
     $        FIRE_DEF_alphamax)
         fire_alpha_factor = fdf_get("MD.FIRE.AlphaFactor",
     &        FIRE_DEF_alpha_factor)
         fire_dtmax = fdf_get("MD.FIRE.MaxTimeStep",
     $        FIRE_DEF_dtmax)
         fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
         fire_debug = fdf_get("MD.FIRE.Debug", .false.)

         if (varcel ) then
            numel = 3*na + 6
         else
            numel = 3*na
         endif

         if (Ionode) then
           write(6,'(a,i6)') "FIRE: No of elements: ", numel
         endif
         call fire_setup(b, n=numel, dt=fire_dt,
     $                   debug=fire_debug,
     $                   dt_inc=fire_dt_inc, dt_dec=fire_dt_dec,
     $                   alphamax=fire_alphamax,
     $                   alpha_factor=fire_alpha_factor,
     $                   nmin=fire_nmin)

        if ( varcel ) then

          constant_volume = fdf_get("MD.ConstantVolume", .false.)

          call get_target_stress(tp,tstres)

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
     .                           9.4486344d0,'Bohr')

C Initialize absolute strain and save initial cell vectors -------------
C Initialization to 1. for numerical reasons, later substracted --------

          strain(1:3,1:3) = 1.0_dp
          cellin(1:3,1:3) = cell(1:3,1:3)
          initial_volume = volcel(cellin)

        else

        ! Do nothing if not variable cell

        endif

        frstme = .false.
      endif


C Variable cell -------------------------------------------------------------

      if ( varcel ) then

        nullify( gfa )
        call re_alloc( gfa, 1, numel, name='gfa',
     &                 routine='fire_optimizer' )
        nullify( gxa )
        call re_alloc( gxa, 1, numel, name='gxa',
     &                 routine='fire_optimizer' )
        nullify( deltamax )
        call re_alloc( deltamax, 1, numel, name='deltamax',
     &                 routine='fire_optimizer' )

C Inverse of matrix of cell vectors  (transpose of) ------------------------

        call reclat( cell, celli, 0 )

C Transform coordinates and forces to fractional 
C but scale them again to Bohr by using the (fixed) moduli of the original
C lattice vectors (allows using maximum displacement as before) 
C convergence is checked here for input forces as compared with ftol

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
              deltamax(indi) = dxmax
            enddo
          enddo
        enddo

C Symmetrizing the stress tensor 

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
              deltamax(indi) = dstrain_max
           enddo
        enddo

C Check stress convergence --------------------------------------------------

        do i = 1, 3
           do j = 1, 3
              relaxd = relaxd .and. 
     .          ( abs(stress_dif(i,j)) .lt. abs(strtol) )
           enddo
        enddo

        if (relaxd) RETURN

        ! Pre-condition forces
!!        gfa = fovermp * gfa / fire_mass
        call fire_step(b,gfa,gxa,deltamax)

C Fixed cell ----------------------------------------------------------------

      else

        relaxd = .true.
        do ia = 1, na
          do n = 1, 3
            relaxd = relaxd .and. ( abs(fa(n,ia)) .lt. ftol )
          enddo
        enddo
        if (relaxd) RETURN

          nullify( rold )
          call re_alloc( rold, 1, numel, name='rold',
     &                   routine='fire_optimizer' )
          nullify( rdiff )
          call re_alloc( rdiff, 1, numel, name='rdiff',
     &                   routine='fire_optimizer' )

           indi = 0
           do i = 1, na
              rold(indi+1:indi+3) = xa(1:3,i)
              rdiff(indi+1:indi+3) = fa(1:3,i)
              indi = indi + 3
           enddo

!!!!!!!           rold = reshape(xa, (/1, 3*na /))
!!!!!!!           rdiff = reshape(fa, (/1, 3*na /))

              ! Pre-condition forces
!!              rdiff = fovermp * rdiff / fire_mass
              call fire_step(b,rdiff,rold,(/ (dxmax, i=1,numel) /) )
              xa(1:3,1:na) = reshape(rold, (/ 3, na /))

      endif

C Transform back if variable cell

      if ( varcel ) then

      ! New cell 

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

      ! Deallocate local memory

        call de_alloc( gxa, name='gxa' )
        call de_alloc( gfa, name='gfa' )
        call de_alloc( deltamax, name='deltamax' )

      else

        call de_alloc( rold, name='rold' )
        call de_alloc( rdiff, name='rdiff' )

      endif ! variable cell

      end subroutine fire_optimizer
      end module m_fire_optim

