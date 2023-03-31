! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_rhog

  ! This module implements most of the functionality necessary to
  ! to mix the Fourier components of the charge density
  !
  ! Theory: as presented in Kresse and Furthmuller
  !  The most basic form is the Kerker mixing, with a parameter related
  !  to the Thomas Fermi screening wavevector.
  !  There is also a Pulay (DIIS) mixer.

  ! This module is not completely self-contained. It interacts with:
  !
  ! - dhscf_init: the allocation of the rhog arrays is done there once
  ! the mesh-related parameters are known.
  !
  ! - dhscf: It can optionally use rhog_in instead of Dscf as starting
  ! point. It also computes and stores rhog by default.
  !
  ! - compute_energies: Apart from correcting EKS when mixing the charge,
  ! it generates rhog from DM_out.
  !
  ! - siesta_forces: In transiesta runs the history is reset upon starting
  !   transiesta
  !
  ! - setup_hamiltonian: it saves rhog_in.
  !

  use precision
  use class_dData1D
  use class_Pair_dData1D
  use class_Fstack_Pair_dData1D

  use m_spin, only: nspin

  implicit none

  ! The (complex) fourier components of rho(G) as stored in a real
  ! array, with dimensions:
  !   complex (2)
  !   np = product(ntml) (local mesh divisions)
  !   nspin
  real(grid_p), pointer, public  :: rhog(:,:,:) => null()
  real(grid_p), pointer, public  :: rhog_in(:,:,:) => null()

  ! these are auxiliary arrays
  real(dp), pointer              :: g2(:) => null()
  logical, pointer               :: g2mask(:) => null()
  integer, pointer               :: gindex(:) => null()
  integer, pointer               :: star_index(:) => null()

  ! Small cutoff for mixing with DIIS/Pulay method
  real(dp)               :: rhog_cutoff 
  integer                :: ng_diis   !  Number of G's in Pulay scheme
                                    !  (in this node)
  real(dp), pointer      :: g2_diis(:) => null()

  logical                :: using_diis_for_rhog
  ! rg_in and rg_diff are one-dimensional arrays holding the
  ! fourier components for selected G-vectors (those with a
  ! norm smaller than rhog_cutoff
  real(dp), pointer              :: rg_in(:) => null()
  real(dp), pointer              :: rg_diff(:) => null()


  real(dp)               :: q0sq    !  Thomas-Fermi K2 for damping
  real(dp)               :: q1sq    !  For scalar product

  type(Fstack_Pair_dData1D), save :: rhog_stack

  integer :: jg0   ! Index of G=0 vector

  ! Callable routines:

  ! Sets up the auxiliary arrays
  public :: order_rhog
  ! Performs the mixing
  public :: mix_rhog
  ! Computes the mismatch between in and out charge densities 
  ! (in Fourier space)
  public :: compute_charge_diff
  ! Deallocates rhog storage
  public :: resetRhog

  private

CONTAINS


  subroutine mix_rhog(iscf)
    use siesta_options, only: wmix
    use precision,      only: dp
    use fdf,            only: fdf_get
    use m_diis,         only: diis

    integer, intent(in) :: iscf

    integer  :: j, i
    real(dp) :: alpha
    real(dp), allocatable :: coeff(:)
    logical  :: mix_first_iter

    if (using_diis_for_rhog) then

       ! Do not use the first diff in the cycle for DIIS
       if (iscf>1) call add_rhog_in_and_diff_to_stack()
       !call print_type(rhog_stack)

       if (n_items(rhog_stack) > 1) then
          allocate(coeff(n_items(rhog_stack)))

          ! Get the DIIS coefficients
          call diis(rhog_stack,scalar_product,coeff)

          ! This will replace the small-G set coefficients
          ! by the DIIS-optimal ones
          call get_optimal_rhog_in()

          deallocate(coeff)
       endif
    endif
           
    ! Do Kerker mixing on the whole fourier series

    ! Check whether we want to mix in the first step
    mix_first_iter = fdf_get("SCF.MixCharge.SCF1",.false.)

    if ((iscf == 1) .and. (.not. mix_first_iter)) then
       ! Do not mix. Take the output density
       rhog_in(:,:,:) = rhog(:,:,:)
    else
       do j = 1, size(rhog_in,dim=2)
          !!! maybe           if (g2(j) <= certain cutoff) then
          if (.true.) then
             alpha = wmix * g2(j) / (g2(j) + q0sq)
             if (alpha == 0) alpha = wmix  ! for G=0
             rhog_in(:,j,:) = alpha * rhog(:,j,:) +  &
                  (1.0_dp - alpha) * rhog_in(:,j,:)
          else
             rhog_in(:,j,:) = rhog(:,j,:)
          endif
       enddo
    endif

  CONTAINS

    subroutine add_rhog_in_and_diff_to_stack()

    ! Store rho_in(G) and rho_diff(G) as single vectors
    ! in a circular stack of the appropriate size

    type(dData1D)      :: vin, vdiff
    type(Pair_dData1D) :: pair
    integer :: ip, i, j, ispin
    character(len=20) :: msg

    ip = 0
    do ispin = 1, nspin
       do j = 1, size(rhog_in,dim=2)
          if (g2mask(j)) then   ! if this G is treated within DIIS
             do i = 1, 2
                ip = ip + 1
                rg_in(ip) = rhog_in(i,j,ispin)
                rg_diff(ip) = rhog(i,j,ispin) - rg_in(ip)
             enddo
          endif
       enddo
    enddo

    write(msg,"(a,i3)") "scf step: ",iscf
    call newdData1D(vin,rg_in,name="(rhog_in -- " // trim(msg) //")")
    call newdData1D(vdiff,rg_diff,name="(rhog_diff -- " // trim(msg) //")")
    call new(pair,vin,vdiff,"(pair in-diff -- " // trim(msg) //")")

    call push(rhog_stack,pair)    ! Store in stack

    call delete(vin)
    call delete(vdiff)
    call delete(pair)

  end subroutine add_rhog_in_and_diff_to_stack

  subroutine get_optimal_rhog_in()

    ! Synthesize the DIIS-optimal rho_in(G) and rho_out(G)
    ! from the DIIS coefficients.

    real(dp), dimension(:), pointer :: vin, vdiff
    type(Pair_dData1D), pointer     :: pairp
    type(dData1D),  pointer         :: vp
    integer :: ip, i, j, ispin, k

    ! zero-out the components treated with DIIS
    do ispin = 1, nspin
       do j = 1, size(rhog_in,dim=2)
          if (g2mask(j)) then
             rhog_in(:,j,ispin) = 0.0_dp
             rhog(:,j,ispin)    = 0.0_dp
          endif
       enddo
    enddo

    do k = 1, n_items(rhog_stack)
       pairp => get_pointer(rhog_stack,k)
       call firstp(pairp,vp)
       vin => val(vp)
       call secondp(pairp,vp)
       vdiff => val(vp)

       ip = 0
       do ispin = 1, nspin
          do j = 1, size(rhog_in,dim=2)
             if (g2mask(j)) then
                do i = 1, 2
                   ip = ip + 1
                   rhog_in(i,j,ispin) = rhog_in(i,j,ispin) + &
                                               coeff(k)* vin(ip)
                   rhog(i,j,ispin)  =   rhog(i,j,ispin) +    &
                                             coeff(k)* (vdiff(ip) + vin(ip))
                enddo
             endif
          enddo
       enddo
    enddo

  end subroutine get_optimal_rhog_in

  end subroutine mix_rhog

!----------------------------------------------------------------
  subroutine compute_charge_diff(drhog)
    use precision
    use m_spin, only: nspin
    use fdf,    only: fdf_get
    use parallel, only: Node

#ifdef MPI
    use m_mpi_utils,           only: globalize_max
#endif

    real(dp), intent(out) :: drhog
#ifdef MPI
    real(dp) :: buffer1
#endif
    integer  :: j, js, i, ispin
    logical  :: debug_stars
    real(dp) :: ss

    ! Note that we now use the complex norm instead of the
    ! abs-norm...
    drhog = -huge(1.0_dp)
    do ispin = 1, nspin
       do j = 1, size(rhog_in,dim=2)
          ss = 0.0_dp
          do i = 1, 2
             ss = ss + (rhog(i,j,ispin)-rhog_in(i,j,ispin))**2
             !drhog = max(drhog,abs(rhog(i,j,ispin)-rhog_in(i,j,ispin)))
          enddo
          drhog = max(drhog,ss)
       enddo
    enddo
    drhog = sqrt(drhog)
    if (Node == 0) print "(a,f12.6)", " Max |\Delta rho(G)|: ", drhog

    ! Print info about the first 10 stars
    debug_stars = fdf_get("SCF.DebugRhogMixing",.false.)
    if (debug_stars .and. Node == 0) then
       print "(a8,2x,a20,4x,a20,2x,a8)", "G2", &
            "rho_in(G) (R, C)", "Diff_rho(G) (R, C)", "damping"
       do ispin = 1, nspin
          do js = 1, 10
             j = gindex(star_index(js))
             print "(f8.4,2x,2f10.5,4x,2f10.5,2x,f8.4)", g2(j), &
                  rhog_in(:,j,ispin) , (rhog(:,j,ispin)-rhog_in(:,j,ispin)), &
                  g2(j) / (g2(j) + q0sq)

          enddo
       enddo
    endif

#ifdef MPI
!     Ensure that drhog is the same on all nodes 
      call globalize_max(drhog,buffer1)
      drhog = buffer1
#endif
    end subroutine compute_charge_diff

  subroutine order_rhog(cell, n1, n2, n3, mesh, nsm)

      ! Sets up indexes for the handling of rho(G)

      use parallel,    only : Node, Nodes, ProcessorY
      use alloc, only: re_alloc
      use fdf,   only: fdf_get, fdf_defined
      use sorting, only: ordix

      use m_mpi_utils,           only: globalize_max
      use m_mpi_utils,           only: globalize_min
      use m_mpi_utils,           only: globalize_sum

      use atomlist,              only: qtot

      implicit none

      real(dp), intent(in):: cell(3,3)
      integer, intent(in) :: n1, n2, n3
      integer, intent(in) :: mesh(3)
      integer, intent(in) :: nsm

      ! Local variables

      real(dp)              :: B(3,3), g(3), celvol
      integer               :: I, I1, I2, I3, IX, J, J1, J2, J3, JX,  &
                               NP, NG, NG2, NG3,                      &
                               ProcessorZ, Py, Pz, J2min, J2max,      &
                               J3min, J3max, J2L, J3L, NRemY, NRemZ,  &
                               BlockSizeY, BlockSizeZ
      external :: reclat
      real(dp), external :: volcel

      real(dp), parameter :: tiny = 1.e-10_dp   ! for regularization
      integer :: n_rhog_depth
      integer :: ng_diis_min, ng_diis_max, ng_diis_sum

      real(dp) :: pi, qtf2, length, length_max, q0_size
      
!
!     Find the genuine thomas fermi k0^2. This does not
!     seem to be too relevant for the preconditioning, as
!     it tends to be too big.
!
      pi = 4*atan(1.0_dp)
      celvol = volcel(cell)
      qtf2 = 4*(3*qtot/(pi*celvol))**(1.0_dp/3.0_dp)

      ! Find the maximum length of a lattice vector
      ! This could set a better scale for the Kerker preconditioning

      length_max = 0.0_dp
      do i = 1, 3
         length = sqrt(dot_product(cell(:,i),cell(:,i)))
         length_max = max(length,length_max)
      enddo
      q0_size = (2*pi/length_max)

      if (Node == 0 ) then
         print "(a,f12.6)", "Thomas-Fermi K2 (Ry):", qtf2
         print "(a,f12.6)", "L max (bohr):", length_max
         print "(a,f12.6)", "q0_size = 2pi/L (Bohr^-1):", q0_size
         print "(a,f12.6)", "q0_size^2 (Ry) :", q0_size**2
      endif

      if (fdf_defined("Thomas.Fermi.K2")) then
         if (Node == 0 ) then
            print "(a,f12.6)", "Please use 'SCF.Kerker.q0sq' " // &
                 "instead of 'Thomas.Fermi.K2'"
         endif
      endif

      q0sq = fdf_get("SCF.Kerker.q0sq", 0.0_dp, "Ry")   ! in Ry
      if (Node == 0 ) then
         print "(a,f12.6)", "Kerker preconditioner q0^2 (Ry):", q0sq
      endif

      q0sq = q0sq  + tiny 

      call re_alloc(g2, 1, n1*n2*n3, "g2", "order_rhog")
      call re_alloc(g2mask, 1, n1*n2*n3, "g2mask", "order_rhog")
      call re_alloc(gindex, 1, n1*n2*n3, "gindex", "order_rhog")

      rhog_cutoff = fdf_get("SCF.RhoG-DIIS-Cutoff", 9.0_dp, "Ry")   ! in Ry
      ng_diis = 0

!     Find reciprocal lattice vectors
      call reclat(CELL, B, 1 )

!     Work out processor grid dimensions

      jg0 = -1

      ProcessorZ = Nodes/ProcessorY
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1

      NG2 = Mesh(2)
      NG3 = Mesh(3)
      BlockSizeY = ((NG2/NSM)/ProcessorY)*NSM
      BlockSizeZ = ((NG3/NSM)/ProcessorZ)*NSM
      NRemY = (NG2 - BlockSizeY*ProcessorY)/NSM
      NRemZ = (NG3 - BlockSizeZ*ProcessorZ)/NSM
      J2min = (Py-1)*BlockSizeY + NSM*min(Py-1,NRemY)
      J2max = J2min + BlockSizeY - 1
      if (Py-1.lt.NRemY) J2max = J2max + NSM
      J2max = min(J2max,NG2-1)
      J3min = (Pz-1)*BlockSizeZ + NSM*min(Pz-1,NRemZ)
      J3max = J3min + BlockSizeZ - 1
      if (Pz-1.lt.NRemZ) J3max = J3max + NSM
      J3max = min(J3max,NG3-1)

      do J3 = J3min,J3max
        if (J3.gt.NG3/2) then
          I3 = J3 - NG3
        else
          I3 = J3
        endif
        do J2 = J2min,J2max
          if (J2.gt.NG2/2) then
            I2 = J2 - NG2
          else
            I2 = J2
          endif
          do J1 = 0,N1-1
            if (J1.gt.N1/2) then
              I1 = J1 - N1
            else
              I1 = J1
            endif
            G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
            G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
            G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3

            J2L = J2 - J2min 
            J3L = J3 - J3min 
            J = 1 + J1 + N1 * J2L + N1 * N2 * J3L

            G2(J) = G(1)**2 + G(2)**2 + G(3)**2

            if (max(abs(i1),abs(i2),abs(i3)) == 0) then
               jg0 = j
            endif

            ! Do not include G=0 in the DIIS subset
            g2mask(j) = ((g2(j) <= rhog_cutoff) .and. (g2(j) > 0))

            !   if (g2mask(j)) then
            !      print "(i5,f10.4,4x,3i6)", j, g2(j), i1, i2, i3
            !   endif

         enddo
      enddo
   enddo
   ! This will work only in serial form for now
   ! Sort by module of G
   call ordix(g2,1,n1*n2*n3,gindex)

   ! Get index of star representatives
   call get_star_reps(g2,gindex,star_index)

   ng_diis = count(g2mask)

   call globalize_min(ng_diis,ng_diis_min)
   call globalize_max(ng_diis,ng_diis_max)
   call globalize_sum(ng_diis,ng_diis_sum)

   if (Node == 0 ) then
      print "(a,i10)", "Number of G's in DIIS: ",  ng_diis_sum
#ifdef MPI
      print "(a,3i10)", "Distrib of G's in DIIS (min, ave, max): ", &
           ng_diis_min, nint(dble(ng_diis_sum)/Nodes), ng_diis_max
#endif
   endif
      
   n_rhog_depth = fdf_get("SCF.RhoG-DIIS-Depth",0)
   
   ! Note that some nodes might not have any Gs in the DIIS procedure
   ! But we still go ahead
   using_diis_for_rhog = ((ng_diis_max > 1) .and. (n_rhog_depth > 1))

   if (using_diis_for_rhog)  call set_up_diis()

   CONTAINS
     
   subroutine get_star_reps(a,aindex,star_index)
     ! Determine star representatives by checking
     ! when the (sorted) modulus of G changes...

     real(dp), intent(in)  :: a(:)
     integer, intent(in)   :: aindex(:)
     integer, pointer      :: star_index(:)

     integer :: j, ng, jj, jp1, ns
     
     ng = size(a)
     
     ns = 0
     do jj = 1, ng-1
        j = aindex(jj)
        jp1 = aindex(jj+1)
        if (abs(a(j) - a(jp1)) > 1.e-7_dp)   ns = ns + 1
     enddo
     !
     call re_alloc(star_index, 1, ns, "star_index", "m_rhog")
     !
     ns = 0
     do jj = 1, ng-1
        j = aindex(jj)
        jp1 = aindex(jj+1)
        if (abs(a(j) - a(jp1)) > 1.e-7_dp) then
           ns = ns + 1
           star_index(ns) = jj
        endif
     enddo
     if (Node == 0) print *, "Number of stars: ", ns

   end subroutine get_star_reps
     
   subroutine set_up_diis()
#ifdef MPI
      use m_mpi_utils, only: globalize_max, globalize_min
#endif

      integer :: ip, ispin, n
      real(dp):: max_g2, min_g2, q1sq_def

#ifdef MPI
      real(dp):: global_max, global_min
#endif

    if (Node == 0) print "(/,a)", "Setting up DIIS for rho(G) -------"
    if (Node == 0) print "(a,i6)",  &
                    "Number of G's treated with DIIS in Node 0: ", ng_diis

    call re_alloc(g2_diis, 1, 2*ng_diis*nspin, "g2_diis", "order_rhog")
    call re_alloc(rg_in, 1, 2*ng_diis*nspin, "rg_in", "order_rhog")
    call re_alloc(rg_diff, 1, 2*ng_diis*nspin, "rg_diff", "order_rhog")

   ! Create the g2_diis array holding |G|^2 for each G in the DIIS set
   ip = 0
   do ispin = 1, nspin
      do j = 1, size(rhog_in,dim=2)
         if (g2mask(j)) then
             do i = 1, 2         ! A bit superfluous
                ip = ip + 1
                g2_diis(ip) = g2(j)
             enddo
          endif
       enddo
    enddo

    min_g2 = minval(g2_diis)
    max_g2 = maxval(g2_diis)
#ifdef MPI
    call globalize_min(min_g2,global_min)
    min_g2 = global_min
    call globalize_max(max_g2,global_max)
    max_g2 = global_max
#endif

    if (Node == 0) then
       print "(a,2f10.3)",  &
             "Minimum and maximum g2 in DIIS: ", min_g2, max_g2
    endif
    ! KF parameter for scalar-product weight
    ! Weight smallest g 20 times more than largest g
    ! (when possible -- if not, set factor n<20)
    n = floor(max_g2/min_g2) - 1
    n = min(n,20)
    q1sq_def = (n-1)*max_g2 / ((max_g2/min_g2) - n)
    if (Node == 0) then
       print "(a,f10.3,a,i3)", &
            "Metric preconditioner cutoff default (Ry): ", q1sq_def, &
            " n:",n
    endif
    q1sq = fdf_get("SCF.RhoG.Metric.Preconditioner.Cutoff",q1sq_def, "Ry")
    if (Node == 0) then
       print "(a,f10.3)", &
                      "Metric preconditioner cutoff (Ry): ", q1sq
       print "(a,2f10.4)", "Max and min weights: ", &
                            (min_g2 + q1sq)/min_g2, &
                            (max_g2 + q1sq)/max_g2
    endif


   call new(rhog_stack,n_rhog_depth,"(rhog DIIS stack)")

  end subroutine set_up_diis

 end subroutine order_rhog

!---
 subroutine resetRhog(continuation)
   use alloc, only: de_alloc
   !< If .true. we don't de-allocate anything, we only reset the history
   logical, intent(in), optional :: continuation
   logical :: lcontinuation
   
   lcontinuation = .false.
   if ( present(continuation) ) lcontinuation = continuation
   
   if ( lcontinuation ) then
     call reset(rhog_stack)
   else
     call de_alloc(rhog_in)
     call de_alloc(rhog)
     call de_alloc(g2)
     call de_alloc(g2mask)
     call de_alloc(gindex)
     call de_alloc(star_index)
     call de_alloc(g2_diis)
     call de_alloc(rg_in)
     call de_alloc(rg_diff)
     call delete(rhog_stack)
   end if
   
 end subroutine resetRhog
 
 ! Auxiliary function to perform the scalar product of residual
 ! vectors in the DIIS method
 function scalar_product(a,b) result (sp)
   real(dp), intent(in) :: a(:), b(:)
   real(dp) :: sp

   integer :: ip, ispin, j
   real(dp):: weight

   sp = 0
   ip = 1
   ! Use standard definition of the scalar product as conjg(A)*B,
   ! but take the real part, as per DIIS equations
   ! If q1sq is not zero, g2_diis(ip) determines the weight, as in KF
   do ispin = 1, nspin
     do j = 1, ng_diis
       weight = (g2_diis(ip) + q1sq) / g2_diis(ip)
       sp = sp + (a(ip)*b(ip) + a(ip+1)*b(ip+1)) * weight
       ip = ip + 2
     enddo
   enddo
   ! Note: This is a "local scalar product"
   ! For efficiency reasons, the All_reduce call is done
   ! on the whole matrix in m_diis
 end function scalar_product

end module m_rhog
