! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_cell_fire_optim
!
!     Fire geometry optimization
!      
      use precision, only : dp
      use sys, only: die
      use alloc
      use fdf, only: fdf_get

      use m_fire
      use siesta_options, only: dt

      use m_memory, only: memory, mem_stat
      use parallel, only : ionode

      implicit none

      public :: cell_fire_optimizer
      private

      CONTAINS
      subroutine cell_fire_optimizer( na, xa, cell, stress,
     $           tp, strtol, varcel, relaxd)
c ***************************************************************************
c Fire geometry optimization (Cell Only)
c
c Written by Alberto Garcia, April 2007

c ******************************** INPUT ************************************
c na                    : number of atoms (will need to keep fractional coords)
c real*8 tp             : Target pressure
c logical varcel        : true if variable cell optimization
c *************************** INPUT AND OUTPUT ******************************
c real*8 cell(3,3)      : Matrix of the vectors defining the unit cell 
c                         input: current step; output: next step
c                         cell(ix,ja) is the ix-th component of ja-th vector
c real*8 xa(3,na)       : atomic coordinates
c real*8 stress(3,3)    : Stress tensor components
c real*8 strtol         : Maximum stress tolerance
c ******************************** OUTPUT ***********************************
c logical relaxd        : True when converged
c ***************************************************************************

C
C  Modules
C

      implicit none
! Subroutine arguments:

      real(dp), intent(in) :: tp,  strtol
      logical, intent(in) :: varcel
      real(dp), intent(inout) :: stress(3,3), cell(3,3)
      integer, intent(in)  :: na
      real(dp), intent(inout) :: xa(3,na)
      logical, intent(out) :: relaxd

c Internal variables and arrays

      real(dp)       :: volume, new_volume, trace

      integer           i, j, indi

      real(dp) ::  celli(3,3)
      real(dp) ::  stress_dif(3,3)

      real(dp), dimension(:), allocatable       :: gxa, gfa
      real(dp), dimension(:), pointer       :: deltamax


! Saved internal variables:

      logical, save :: frstme = .true.
      logical, save :: constant_volume
      real(dp), save :: initial_volume


      real(dp), save :: target_stress(3,3),
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
      real(dp), parameter ::  dstrain_max = 0.1_dp

      real(dp) :: volcel
      external :: volcel
c ---------------------------------------------------------------------------

      volume = volcel(cell)

      if ( frstme ) then
         fire_mass = fdf_get("MD.FIRE.Mass", 1.0_dp)
         fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
         fire_dt_inc = fdf_get("MD.FIRE.TimeInc", FIRE_DEF_dt_inc)
         fire_dt_dec = fdf_get("MD.FIRE.TimeDec", FIRE_DEF_dt_dec)
         fire_nmin = fdf_get("MD.FIRE.Nmin", FIRE_DEF_nmin)
         fire_alphamax = fdf_get("MD.FIRE.AlphaMax", FIRE_DEF_alphamax)
         fire_alpha_factor = fdf_get("MD.FIRE.AlphaFactor",
     &        FIRE_DEF_alpha_factor)
         fire_dtmax = fdf_get("MD.FIRE.MaxTimeStep", FIRE_DEF_dtmax)
         fire_dt = fdf_get("MD.FIRE.TimeStep", dt)
         fire_debug = fdf_get("MD.FIRE.Debug", .false.)

         
         if (varcel ) then
            numel =  6
         else
            call die("no varcel?")
         endif
         if (Ionode) then
           write(6,'(a,i6)') "Cell_Fire_optim: No of elements: ",
     $                numel
         endif

         call fire_setup(b, n=numel, dt=fire_dt,
     $                   debug=fire_debug,
     $                   dt_inc=fire_dt_inc, dt_dec=fire_dt_dec,
     $                   alphamax=fire_alphamax,
     $                   alpha_factor=fire_alpha_factor,
     $                   nmin=fire_nmin)

        if ( varcel ) then

C Check if we want a constant-volume simulation
          constant_volume = fdf_get("MD.ConstantVolume", .false.)

          call get_target_stress(tp,target_stress)
          if (constant_volume) target_stress = 0.0_dp


C Scale factor for strain variables to share magnitude with coordinates 
C ---- (a length in Bohrs typical of bond lengths ..) ------------------

          precon = fdf_get('MD.PreconditionVariableCell',
     .                           9.4486344d0,'Bohr')

C Initialize absolute strain and save initial cell vectors -------------
C Initialization to 1. for numerical reasons, later substracted --------

          strain(1:3,1:3) = 1.0_dp
          cellin(1:3,1:3) = cell(1:3,1:3)
          initial_volume = volcel(cellin)

        endif

        frstme = .false.
      endif                     ! First call

C Variable cell -------------------------------------------------------------

      if ( varcel ) then

         allocate(gfa(numel), stat=mem_stat)
         call memory('A','D',numel,'Fire_optim',
     $        stat=mem_stat,id="gfa")
         allocate(gxa(numel), stat=mem_stat)
         call memory('A','D',numel,'Fire_optim',
     $        stat=mem_stat,id="gxa")
        nullify( deltamax )
        call re_alloc( deltamax, 1, numel, name='deltamax',
     &                 routine='fire_optimizer' )


        relaxd = .true.

C Symmetrizing the stress tensor 

        do i = 1, 3
           do j = i+1, 3
              stress(i,j) = 0.5_dp*( stress(i,j) + stress(j,i) )
              stress(j,i) = stress(i,j)
           enddo
        enddo

C Subtract target stress

        stress_dif = stress - target_stress
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

        indi = 0
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

C Call Fire step

        call fire_step(b,gfa,gxa,deltamax)

      endif

C Transform back if variable cell

      if ( varcel ) then

      ! New cell 

        indi = 0
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              strain(i,j) = gxa(indi) / precon
              strain(j,i) = strain(i,j)
           enddo
        enddo

        ! Inverse of matrix of cell vectors  (transpose of)
        call reclat( cell, celli, 0 )

!       Update cell

        cell = cellin + matmul(strain-1.0_dp,cellin)
        if (constant_volume) then
           new_volume = volcel(cell)
           cell =  cell * (initial_volume/new_volume)**(1.0_dp/3.0_dp)
        endif

        call rescale_coordinates(na,xa, celli_oldcell=celli,
     $                                  new_cell=cell)


      ! Deallocate local memory

        deallocate (gxa, stat=mem_stat)
        call memory('D','D',numel,'Fire_optim',
     $       stat=mem_stat,id="gxa")
        deallocate (gfa, stat=mem_stat)
        call memory('D','D',numel,'Fire_optim',
     $       stat=mem_stat,id="gfa")
        call de_alloc( deltamax, name='deltamax' )

      endif ! variable cell

      end subroutine cell_fire_optimizer
!
!--------------------------------------------------------------

      subroutine rescale_coordinates(na,xa,
     $                            celli_oldcell,new_cell)
      use precision, only : dp
      use zmatrix

      integer, intent(in)     :: na
      real(dp), intent(inout) :: xa(3,na)
      real(dp), intent(in)    :: celli_oldcell(3,3)
      real(dp), intent(in)    :: new_cell(3,3)


      real(dp), dimension(3) :: r, frac
      integer  :: ifirst, imol, icart, i, j

      !  NOTE: We have to be careful here if using a Zmatrix
      if (lUseZmatrix) then

        !     re-scale only the positions of the leading atoms
        !     in each molecule,
        !     plus any cartesian block,
        !     and recompute the cartesian coordinate array
        !     
           do imol = 1, nZmol
              ifirst = nZmolStartAtom(imol)
              r(1:3) = Zmat(3*(ifirst-1)+1:3*(ifirst-1)+3)
              frac(1:3) = matmul(transpose(celli_oldcell),r(1:3))
              r(1:3) = matmul(new_cell,frac(1:3))
              Zmat(3*(ifirst-1)+1:3*(ifirst-1)+3) = r(1:3)
           enddo
           do icart = 1, nZcart
             ! Process cartesian block
              ifirst = nZcartStartAtom(icart)
              do j = ifirst, ifirst + nZcartAtoms(icart) - 1
                 r(1:3) = Zmat(3*(j-1)+1:3*(j-1)+3)
                 frac(1:3) = matmul(transpose(celli_oldcell),r(1:3))
                 Zmat(3*(j-1)+1:3*(j-1)+3) = matmul(new_cell,frac(1:3))
              enddo
           enddo
           call Zmat_to_Cartesian(xa)
        else  
           ! No Zmatrix
           ! Rescale coordinates for all atoms
           do i = 1, na
              xa(:,i) = matmul(transpose(celli_oldcell),xa(:,i))
              xa(:,i) = matmul(new_cell,xa(:,i))
           enddo
           
        endif

      end subroutine rescale_coordinates


      end module m_cell_fire_optim

