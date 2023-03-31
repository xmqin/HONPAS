! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_cell_broyden_optim
!
!     Broyden geometry optimization
!      
      use precision, only : dp
      use fdf
      use units, only: kBar
      use sys, only: die
!
      use m_broyddj_nocomm, only: broyden_t
      use m_broyddj_nocomm, only: broyden_reset, broyden_step,
     $                     broyden_destroy, broyden_init,
     $                     broyden_is_setup

      use parallel, only : ionode

      implicit none

      public :: cell_broyden_optimizer
      private

      CONTAINS
      subroutine cell_broyden_optimizer( na, xa, cell, stress,
     $                tp, strtol, varcel, relaxd)
c ***************************************************************************
c Broyden geometry optimization (Cell Only)
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

! Subroutine arguments:

      real(dp), intent(in) :: tp,  strtol
      logical, intent(in) :: varcel
      real(dp), intent(inout) :: stress(3,3), cell(3,3)
      integer, intent(in)  :: na
      real(dp), intent(inout) :: xa(3,na)
      logical, intent(out) :: relaxd

c Internal variables and arrays
      real(dp) :: volume, new_volume, trace
      integer  :: i, j, indi
      real(dp) :: celli(3,3), sxx, syy, szz, sxy, sxz, syz
      real(dp) :: stress_dif(3,3)

      real(dp), allocatable ::  rnew(:)
      real(dp), dimension(:), allocatable  :: gxa, gfa, max_step

      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

! Saved internal variables:

      logical, save :: frstme = .true.
      logical, save :: constant_volume
      real(dp), save :: initial_volume


      real(dp), save :: target_stress(3,3),
     .                  precon,
     .                  strain(3,3),
     .                  cellin(3,3)

      type(broyden_t), save  :: br
      logical, save           :: initialization_done = .false.
      integer, save  :: numel

      real(dp), save :: jinv0, max_step_strain
      integer, save  :: maxit
      logical, save  :: cycle_on_maxit, variable_weight
      logical, save  :: broyden_debug
      real(dp)       :: weight
      logical        :: tarstr

      real(dp) :: volcel
      external :: volcel
c ---------------------------------------------------------------------------

      volume = volcel(cell)

      if (.not. initialization_done) then

        maxit           = fdf_get("MD.Broyden.History.Steps",5)
        cycle_on_maxit  = fdf_get("MD.Broyden.Cycle.On.Maxit",.true.)
        variable_weight = fdf_get("MD.Broyden.Variable.Weight",.false.)
        broyden_debug   = fdf_get("MD.Broyden.Debug",.false.)
        jinv0           = fdf_get("MD.Broyden.Initial.Inverse.Jacobian",
     $                            1.0_dp)
        Max_step_strain = fdf_get(
     $         "MD.MaxCGDispl",0.02_dp,"Bohr")

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
          numel =  6
        else
          call die("no varcel?")
        endif
        if (Ionode) then
          write(6,'(a,i6)') "Cell_Broyden_optim: No of elements: ",
     $                      numel
        endif

        if ( varcel ) then

C Check if we want a constant-volume simulation
          constant_volume = fdf_get("MD.ConstantVolume", .false.)

C Look for target stress and read it if found, otherwise generate it --------

          tarstr = fdf_block('MD.TargetStress',bfdf)

          if (tarstr) then
            if (ionode) then
              write(6,'(/a,a)')
     $            'Broyden_optim: Reading %block MD.TargetStress',
     .                      ' (units of MD.TargetPressure).'
            endif
            if (.not. fdf_bline(bfdf,pline))
     $        call die('cell_broyden_optimizer: ERROR in ' //
     .                 'MD.TargetStress block')
            sxx = fdf_bvalues(pline,1)
            syy = fdf_bvalues(pline,2)
            szz = fdf_bvalues(pline,3)
            sxy = fdf_bvalues(pline,4)
            sxz = fdf_bvalues(pline,5)
            syz = fdf_bvalues(pline,6)
            call fdf_bclose(bfdf)
            
            target_stress(1,1) = - sxx * tp
            target_stress(2,2) = - syy * tp
            target_stress(3,3) = - szz * tp
            target_stress(1,2) = - sxy * tp
            target_stress(2,1) = - sxy * tp
            target_stress(1,3) = - sxz * tp
            target_stress(3,1) = - sxz * tp
            target_stress(2,3) = - syz * tp
            target_stress(3,2) = - syz * tp
          else
            if (ionode) then
              write(6,'(/a,a)')
     $            'Broyden_optim: No target stress found, ',
     .            'assuming hydrostatic MD.TargetPressure.'
            endif
            target_stress(:,:) = 0.0_dp
            do i= 1, 3
              target_stress(i,i) = - tp
            enddo
          endif

C Write target stress down --------------------------------------------------

          if (ionode) then
            write(6,"(/a)") 'Broyden_optim: Target stress (kBar)'
            do i=1, 3
              write(6,"(a,2x,3f12.3)") 
     .            'Broyden_optim:', (target_stress(i,j)/kBar, j=1,3)
            enddo
          endif

C Scale factor for strain variables to share magnitude with coordinates 
C ---- (a length in Bohrs typical of bond lengths ..) ------------------

          precon = fdf_get('MD.PreconditionVariableCell',
     .                     9.4486344d0,'Bohr')

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

         allocate(gfa(1:numel),gxa(1:numel))
         allocate(rnew(1:numel),max_step(1:numel))

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

        if (.not. broyden_is_setup(br)) then
           call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $                       jinv0,0.01_dp,max_step)
        endif

        weight = 1.0_dp
        call broyden_step(br,gxa,gfa,w=weight,newx=rnew)

      endif

C Transform back if variable cell

      if ( varcel ) then

      ! New cell 

        indi = 0
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              strain(i,j) = rnew(indi) / precon
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
        deallocate(rnew,gxa,gfa,max_step)

      endif ! variable cell

      end subroutine cell_broyden_optimizer
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


      end module m_cell_broyden_optim

