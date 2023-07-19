! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_broyden_optim
!
!     Broyden geometry optimization
!      
      use precision, only : dp
      use fdf,          only: fdf_block
      use m_fdf_global, only: fdf_global_get
      use m_mpi_utils, only:  broadcast
      use units, only: kBar
!
      use m_broyddj_nocomm, only: broyden_t
      use m_broyddj_nocomm, only: broyden_reset, broyden_step,
     $                     broyden_destroy, broyden_init,
     $                     broyden_is_setup

      use m_memory, only: memory, mem_stat
      use parallel, only : ionode

      implicit none

      public :: broyden_optimizer
      private

      CONTAINS
      subroutine broyden_optimizer( na, xa, fa, cell, stress,
     $     dxmax, tp, ftol, strtol, varcel, relaxd)
c ***************************************************************************
c Broyden geometry optimization
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
c Written by Alberto Garcia, December 2005

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

C
C  Modules
C

! Subroutine arguments:

      integer, intent(in) :: na
      real(dp), intent(in) :: fa(3,na), dxmax,
     .                        tp, ftol, strtol
      logical, intent(in) :: varcel
      real(dp), intent(inout) :: xa(3,na), stress(3,3), cell(3,3)
      logical, intent(out) :: relaxd

c Internal variables and arrays

      real(dp)            :: volume, new_volume, trace
      real(dp)            :: Max_step_strain, Max_step_coordinates

      logical           found
      integer           iu, ia, i, j, n, indi

      real(dp) ::  celli(3,3), sxx, syy, szz, sxy, sxz, syz
      real(dp) ::  stress_dif(3,3)

      real(dp), dimension(:), allocatable       :: gxa, gfa, max_step
      real(dp), dimension(:), allocatable       :: rold, rnew, rdiff

! Saved internal variables:

      logical, save :: frstme = .true.
      logical, save :: tarstr = .false.
      logical, save :: constant_volume
      real(dp), save :: initial_volume


      real(dp), save :: tstres(3,3),
     .                  modcel(3),
     .                  precon,
     .                  strain(3,3),
     .                  cellin(3,3)

      type(broyden_t), save  :: br
      logical, save           :: initialization_done = .false.
      integer, save  :: numel

      real(dp), save :: jinv0
      integer, save  :: maxit
      logical, save  :: cycle_on_maxit, variable_weight
      logical, save  :: broyden_debug
      real(dp)       :: weight

      real(dp) :: volcel
      external :: volcel
c ---------------------------------------------------------------------------

      volume = volcel(cell)

      if (.not. initialization_done) then

          call fdf_global_get(maxit,"MD.Broyden.History.Steps",5)
          call fdf_global_get(cycle_on_maxit,
     $              "MD.Broyden.Cycle.On.Maxit",.true.)
          call fdf_global_get(variable_weight,
     $         "MD.Broyden.Variable.Weight",.false.)
          call fdf_global_get(broyden_debug,
     $         "MD.Broyden.Debug",.false.)
          call fdf_global_get(jinv0,
     $         "MD.Broyden.Initial.Inverse.Jacobian",1.0_dp)

         if (ionode) then
          write(6,*) 
          write(6,'(a,i3)')
     $         "Broyden_optim: max_history for broyden: ", maxit
          write(6,'(a,l1)')
     $         "Broyden_optim: cycle on maxit: ", cycle_on_maxit
          if (variable_weight) write(6,'(a)')
     $         "Broyden_optim: Variable weight not implemented yet"
          write(6,'(a,f8.4)')
     $         "Broyden_optim: initial inverse jacobian: ", jinv0
          write(6,*) 
        endif

         call broyden_init(br,debug=broyden_debug)
         initialization_done = .true.

      endif


C If first call to optim, check dim and get target stress --------------------

      if ( frstme ) then
  
         if (varcel ) then
            numel = 3*na + 6
         else
            numel = 3*na
         endif
         if (Ionode) then
           write(6,'(a,i6)') "Broyden_optim: No of elements: ", numel
         endif

        if ( varcel ) then

C Check if we want a constant-volume simulation
          call fdf_global_get(constant_volume,
     .                        "MD.ConstantVolume", .false.)

C Look for target stress and read it if found, otherwise generate it --------

          if (ionode) then
            tarstr = fdf_block('MD.TargetStress',iu)

            if (tarstr) then
               write(6,'(/a,a)')
     $              'Broyden_optim: Reading %block MD.TargetStress',
     .                        ' (units of MD.TargetPressure).'
               read(iu,*, end=50) sxx, syy, szz, sxy, sxz, syz
               tstres(1,1) = - sxx * tp
               tstres(2,2) = - syy * tp
               tstres(3,3) = - szz * tp
               tstres(1,2) = - sxy * tp
               tstres(2,1) = - sxy * tp
               tstres(1,3) = - sxz * tp
               tstres(3,1) = - sxz * tp
               tstres(2,3) = - syz * tp
               tstres(3,2) = - syz * tp
   50          continue
            else
              write(6,'(/a,a)')
     $              'Broyden_optim: No target stress found, ',
     .                'assuming hydrostatic MD.TargetPressure.'
              do i = 1, 3
                do j = 1, 3
                  tstres(i,j) = 0._dp
                enddo
                tstres(i,i) = - tp
              enddo
            endif

C Write target stress down --------------------------------------------------

            write(6,"(/a)") 'Broyden_optim: Target stress (kBar)'
            write(6,"(a,2x,3f12.3)") 
     .       'Broyden_optim:', tstres(1,1)/kBar, tstres(1,2)/kBar, 
     .       tstres(1,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'Broyden_optim:', tstres(2,1)/kBar, tstres(2,2)/kBar, 
     .       tstres(2,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'Broyden_optim:', tstres(3,1)/kBar, tstres(3,2)/kBar, 
     .       tstres(3,3)/kBar
          endif  ! node0

          call broadcast(tstres(1:3,1:3))


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

          call fdf_global_get(precon,'MD.PreconditionVariableCell',
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


!     Broyden section
!

      Max_step_strain = dxmax
      Max_step_coordinates = dxmax

C Variable cell -------------------------------------------------------------

      if ( varcel ) then

         allocate(gfa(numel), stat=mem_stat)
         call memory('A','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="gfa")
         allocate(gxa(numel), stat=mem_stat)
         call memory('A','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="gxa")
         allocate(max_step(numel), stat=mem_stat)
         call memory('A','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="max_step")
         allocate (rnew(numel), stat=mem_stat)
         call memory('A','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="rnew")

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
              max_step(indi) = Max_step_coordinates
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
              max_step(indi) = Max_step_strain
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

C Call Broyden step

        if (.not. broyden_is_setup(br)) 
     $    call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $                       jinv0,0.01_dp,max_step)

           weight = 1.0_dp
           call broyden_step(br,gxa,gfa,w=weight,newx=rnew)

C Fixed cell ----------------------------------------------------------------

      else

        relaxd = .true.
        do ia = 1, na
          do n = 1, 3
            relaxd = relaxd .and. ( abs(fa(n,ia)) .lt. ftol )
          enddo
        enddo
        if (relaxd) RETURN


           allocate (rold(numel), stat=mem_stat)
           call memory('A','D',numel,'Broyden_optim',
     $          stat=mem_stat,id="rold")
           allocate (rnew(numel), stat=mem_stat)
           call memory('A','D',numel,'Broyden_optim',
     $          stat=mem_stat,id="rnew")
           allocate (rdiff(numel),stat=mem_stat)
           call memory('A','D',numel,'Broyden_optim',
     $          stat=mem_stat,id="rdiff")
           allocate (max_step(numel),stat=mem_stat)
           call memory('A','D',numel,'Broyden_optim',
     $          stat=mem_stat,id="max_step")

           weight = 1.0_dp
           indi = 0
           do i = 1, na
              rold(indi+1:indi+3) = xa(1:3,i)
              rdiff(indi+1:indi+3) = fa(1:3,i)
              max_step(indi+1:indi+3) = Max_step_coordinates
              indi = indi + 3
           enddo

           if (.not. broyden_is_setup(br)) then

              call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $             jinv0,0.01_dp,max_step)

           endif

!!!!!!!           rold = reshape(xa, (/1, 3*na /))
!!!!!!!           rdiff = reshape(fa, (/1, 3*na /))
           
           call broyden_step(br,rold,rdiff,w=weight,newx=rnew)
           xa(1:3,1:na) = reshape(rnew, (/ 3, na /))

      endif

C Transform back if variable cell

      if ( varcel ) then

      ! New cell 

        indi = 3*na
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              strain(i,j) = rnew(indi) / precon
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
              xa(i,ia) = xa(i,ia) + cell(i,n) * rnew(indi) / modcel(n)
            enddo
          enddo
        enddo

      ! Deallocate local memory

        deallocate (rnew, stat=mem_stat)
        call memory('D','D',numel,'Broyden_optim',
     $       stat=mem_stat,id="rnew")
        deallocate (gxa, stat=mem_stat)
        call memory('D','D',numel,'Broyden_optim',
     $       stat=mem_stat,id="gxa")
        deallocate (gfa, stat=mem_stat)
        call memory('D','D',numel,'Broyden_optim',
     $       stat=mem_stat,id="gfa")
        deallocate (max_step, stat=mem_stat)
        call memory('D','D',numel,'Broyden_optim',
     $       stat=mem_stat,id="max_step")

      else

         deallocate (rold, stat=mem_stat)
         call memory('D','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="rold")
         deallocate (rnew, stat=mem_stat)
         call memory('D','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="rnew")
         deallocate (rdiff, stat=mem_stat)
         call memory('D','D',numel,'Broyden_optim',
     $        stat=mem_stat,id="rdiff")

      endif ! variable cell

      end subroutine broyden_optimizer
      end module m_broyden_optim

