! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_broyden_mixing

      use precision, only: dp

      use m_broyddj, only: broyden_t, broyden_init, broyden_is_setup
      use m_broyddj, only: broyden_reset, broyden_step
      use fdf
      use alloc,     only: re_alloc, de_alloc

      use parallel, only: ionode
      use m_mpi_utils, only: globalize_sum, globalize_max

      public :: broyden_mixing
      private

      CONTAINS

      subroutine broyden_mixing(iscf,mix_scf1,nbasis,maxnd,numd,
     .                   listdptr,nspin,alpha,nkick,alpha_kick,
     $                   dmnew,dmold,dmax)

C ************************** INPUT **************************************

C integer iscf               : Current SCF iteration
C integer nbasis             : Number of atomic orbitals stored locally
C integer maxnd              : First dimension of D.M., and 
C                              maximum number of nonzero elements of D.M.
C integer numd(:)            : Control vector of D.M.
C                              (number of nonzero elements of each row)
C integer nspin              : Spin polarization (1=unpolarized, 2=polarized,
C                               4=Non-collinear)
C real*8 alpha               : Mixing parameter (for linear mixing)
C integer nkick              : A kick is given every nkick iterations
C real*8 alpha_kick          : Mixing parameter for kicks
C logical mix_scf1           : Mix on first iteration?
C ********************* INPUT AND OUTPUT*********************************
C real*8 dmnew(maxnd)        : Density Matrix
C                              Input: d.m. output in current SCF step
C                              Output: d.m. input for next SCF iteration
C real*8 dmold(maxnd)        : Density matrix
C                              Input: d.m. input in current SCF step
C                              Output: d.m. input for next SCF iteration
C ************************** OUTPUT *************************************
C real*8 dmax                : Maximum change of a DM element between 
C                              input and output

      implicit none

      integer, intent(in) :: iscf,maxnd,nbasis,nspin,nkick
      logical, intent(in) :: mix_scf1
      integer, intent(in) ::  numd(:), listdptr(:)
      real(dp), intent(in) :: alpha, alpha_kick
      real(dp), intent(inout) :: dmnew(maxnd,nspin),
     $                           dmold(maxnd,nspin)
      real(dp), intent(out) ::  dmax


      logical, save           :: initialization_done = .false.


      integer ::   i0,i,is,j, numel,k, global_numel

      real(dp), dimension(:), pointer  :: rold, rnew, rdiff 

      type(broyden_t), save ::  br

      real(dp), save :: jinv0
      integer, save  :: maxit
      logical, save  :: cycle_on_maxit, variable_weight
      logical, save  :: do_not_mix, broyden_debug

      real(dp) :: global_dnorm, global_dmax,  dnorm, diff, weight

      numel = nspin * sum(numd(1:nbasis))
      call Globalize_sum(numel,global_numel)

      if (.not. initialization_done) then

        if (ionode) then
          print *, "Broyden: No of relevant DM elements: ", global_numel
        endif
        maxit = fdf_integer("DM.Number.Broyden",5)
        cycle_on_maxit =
     $          fdf_boolean("DM.Broyden.Cycle.On.Maxit",.true.)
        variable_weight =
     $          fdf_boolean("DM.Broyden.Variable.Weight",.true.)
        broyden_debug = 
     $          fdf_boolean("DM.Broyden.Debug",.false.)

        jinv0 = fdf_double("DM.Broyden.Initial.Mixing",alpha)
        if (ionode) then
          print *, "maxit for broyden: ", maxit
          print *, "cycle on maxit: ", cycle_on_maxit
          print *, "variable weight: ", variable_weight
          print *, "initial alpha: ", jinv0
        endif

        call broyden_init(br,broyden_debug)

        initialization_done = .true.

      endif

      do_not_mix = (iscf == 1 .and. .not. mix_scf1)
      if (kick_is_due(iscf,nkick) .or. do_not_mix) then

          ! Kick without saving any history for later
         if (broyden_debug .and. ionode) then
            if (do_not_mix) then
               print *, "No mix in first iteration"
            else
               print *, "Kick"
            endif
         endif

        ! Linear mixing with alpha_kick (or no mixing)

         dmax = 0.0_dp
         dnorm = 0.0_dp
         do is = 1,nspin
            do i = 1,nbasis
               do j = 1,numd(i)
                  k = listdptr(i) + j
                  diff = dmnew(k,is) - dmold(k,is)
                  dmax = max(dmax, abs(diff))
                  dnorm = dnorm + diff**2
                  if (.not. do_not_mix) then
                     dmnew(k,is) = dmold(k,is) + alpha_kick*diff
                  endif
                  dmold(k,is) = dmnew(k,is)
               enddo
            enddo
         enddo
         call Globalize_sum(dnorm,global_dnorm)
         call Globalize_max(dmax,global_dmax)
!
         global_dnorm = sqrt(global_dnorm)
         dmax = global_dmax
         if (broyden_debug .and. ionode)
     $                   print *, "Dnorm = ", global_dnorm

         if (broyden_is_setup(br)) then
            if (broyden_debug .and. ionode)
     $           print *, "Resetting history after kick or 1st."
            call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $           jinv0,0.01_dp)
         endif

         RETURN

      endif       ! Linear mixing
!
!     Broyden section
!
      nullify( rold )
      call re_alloc( rold, 1, numel, name='rold',
     &               routine='broyden_mixing' )
      nullify( rnew )
      call re_alloc( rnew, 1, numel, name='rnew',
     &               routine='broyden_mixing' )
      nullify( rdiff )
      call re_alloc( rdiff, 1, numel, name='rdiff',
     &               routine='broyden_mixing' )
!
!          Copy input to auxiliary arrays
!          (memory will be saved by inlining the whole thing
!           later on)
!
           i0 = 0
           dmax = 0.0_dp
           dnorm = 0.0_dp
           do is = 1,nspin
             do i = 1,nbasis
               do j = 1,numd(i)
                 i0 = i0 + 1
                 k = listdptr(i) + j
                 rold(i0) = dmold(k,is)
                 rdiff(i0) = dmnew(k,is) -  dmold(k,is)
                 dmax = max(dmax, abs(rdiff(i0)))
                 dnorm = dnorm + rdiff(i0)**2
               enddo
             enddo
           enddo
           call Globalize_sum(dnorm,global_dnorm)
           call Globalize_max(dmax,global_dmax)
!
           global_dnorm = sqrt(global_dnorm)
           dmax = global_dmax
           if (broyden_debug .and. ionode)
     $                  print *, "Dnorm = ", global_dnorm
!

           if (iscf == 1) then
              if (broyden_debug .and. ionode) then
              ! Forget about history upon new scf cycle
                 print *, "Resetting history for new SCF cycle."
                 print *, "Broyden: No of relevant DM elements: ",
     $                global_numel
              endif
              call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $             jinv0,0.01_dp)

           endif

           if (.not. broyden_is_setup(br)) then
              if (broyden_debug .and. ionode) then
                 print *, "Broyden: No of relevant DM elements: ",
     $                global_numel
              endif
              call broyden_reset(br,numel,maxit,cycle_on_maxit,
     $             jinv0,0.01_dp)

           endif

!          Broyden step
!
           if (variable_weight) then
!
!            Heuristic weight 1 < w < 300
!
              weight = exp(1.0_dp/(global_dmax+0.20))
              if (broyden_debug .and. ionode)
     $                     print *, "weight: ", weight
           else
              weight = 1.0_dp
           endif

           call broyden_step(br,rold,rdiff,w=weight,newx=rnew)

!
!         Copy back the results
!
           i0 = 0
           do is = 1,nspin
             do i = 1,nbasis
               do j = 1,numd(i)
                 i0 = i0 + 1
                 k = listdptr(i)+j
                 dmnew(k,is) = rnew(i0)
                 dmold(k,is) = dmnew(k,is)
               enddo
             enddo
           enddo

           call de_alloc( rold, name='rold' )
           call de_alloc( rnew, name='rnew' )
           call de_alloc( rdiff, name='rdiff' )

      end subroutine broyden_mixing

!---------------------------------------------------------
      function kick_is_due(n,step) result(due)
      logical :: due
      integer, intent(in) :: n, step
      
      due = .false.
      if (step .le. 0) RETURN
      due = (modulo(n,step) == 0)

      end function kick_is_due
!---------------------------------------------------------
      
      end module m_broyden_mixing

