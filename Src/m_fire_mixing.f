! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_fire_mixing

      use precision, only: dp

      use m_fire_para
      use fdf, only: fdf_get
      use alloc,     only: re_alloc, de_alloc

      use parallel, only: ionode
      use m_mpi_utils, only: broadcast, globalize_sum, globalize_max

      public :: fire_mixing
      private

      CONTAINS

      subroutine fire_mixing(iscf,mix_scf1,nbasis,maxnd,numd,
     .                   listdptr,nspin,alpha,nkick,alpha_kick,
     $                   dmnew,dmold,dmax)

C ************************** INPUT **************************************

C integer iscf               : Current SCF iteration
C integer nbasis             : Number of atomic orbitals stored locally
C integer maxnd              : First dimension of D.M., and 
C                              maximum number of nonzero elements of D.M.
C integer numd(:)            : Control vector of D.M.
C                              (number of nonzero elements of each row)
C integer nspin              : Spin polarization (1=unpolarized, 2=polarized, 4=non-collinear)
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

      real(dp), dimension(:), pointer  :: rold, rdiff

      type(fire_t), save ::  b

      logical  :: do_not_mix, fire_debug

      real(dp), save :: fire_mass
      real(dp)       :: fire_dt, fire_dt_inc,
     $                  fire_dt_dec, fire_alphamax,
     $                  fire_alpha_factor, fire_dtmax
      integer        :: fire_nmin

      real(dp) :: global_dnorm, global_dmax,  dnorm, diff

      numel = nspin * sum(numd(1:nbasis))
      call Globalize_sum(numel,global_numel)

      if (.not. initialization_done) then

         fire_dt = alpha        ! Fix this, unless the user disagrees

         fire_mass = fdf_get("DM.FIRE.Mass", 1.0_dp)
         fire_dt_inc = fdf_get("DM.FIRE.TimeInc", FIRE_DEF_dt_inc)
         fire_dt_dec = fdf_get("DM.FIRE.TimeDec", FIRE_DEF_dt_dec)
         fire_nmin = fdf_get("DM.FIRE.Nmin", FIRE_DEF_nmin)
         fire_alphamax = fdf_get("DM.FIRE.AlphaMax", FIRE_DEF_alphamax)
         fire_alpha_factor = fdf_get("DM.FIRE.AlphaFactor",
     &        FIRE_DEF_alpha_factor)
         fire_dtmax = fdf_get("DM.FIRE.MaxTimeStep", FIRE_DEF_dtmax)
         fire_debug = fdf_get("DM.FIRE.Debug", .true.)

         if (ionode) then
            print *, "Fire: No of relevant DM elements: ",
     $           numel, global_numel
         endif
         
         call fire_setup(b, n=numel, dt=fire_dt,
     $        debug=fire_debug,
     $        dt_inc=fire_dt_inc, dt_dec=fire_dt_dec,
     $        alphamax=fire_alphamax,
     $        alpha_factor=fire_alpha_factor,
     $        nmin=fire_nmin)

         initialization_done = .true.
          
      endif

      do_not_mix = (iscf == 1 .and. .not. mix_scf1)
      if (do_not_mix) then

         dmax = 0.0_dp
         dnorm = 0.0_dp
         do is = 1,nspin
            do i = 1,nbasis
               do j = 1,numd(i)
                  k = listdptr(i) + j
                  diff = dmnew(k,is) - dmold(k,is)
                  dmax = max(dmax, abs(diff))
                  dnorm = dnorm + diff**2
                  dmold(k,is) = dmnew(k,is)
               enddo
            enddo
         enddo

         call Globalize_sum(dnorm,global_dnorm)
         call Globalize_max(dmax,global_dmax)
!
         global_dnorm = sqrt(global_dnorm)
         dmax = global_dmax
         if (fire_debug .and. ionode)
     $                   print *, "Dnorm = ", global_dnorm

         RETURN

         endif



!     Fire section
!
      nullify( rold )
      call re_alloc( rold, 1, numel, name='rold',
     &               routine='fire_mixing' )
      nullify( rdiff )
      call re_alloc( rdiff, 1, numel, name='rdiff',
     &               routine='fire_mixing' )
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
           if (fire_debug .and. ionode)
     $                  print *, "Dnorm = ", global_dnorm
!

!          Fire step
!
           call fire_step(b,rdiff,rold, (/ (dmax, i=1,numel) /) )
!
!         Copy back the results
!
           i0 = 0
           do is = 1,nspin
             do i = 1,nbasis
               do j = 1,numd(i)
                 i0 = i0 + 1
                 k = listdptr(i)+j
                 dmnew(k,is) = rold(i0)
                 dmold(k,is) = dmnew(k,is)
               enddo
             enddo
           enddo

           call de_alloc( rold, name='rold' )
           call de_alloc( rdiff, name='rdiff' )

      end subroutine fire_mixing

      end module m_fire_mixing

