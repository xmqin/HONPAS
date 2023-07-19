c
      subroutine compat_params(str)
c
c     Some internal parameters tend to change over time... This routine
c     assigns values to them based on a set of directives, specified
c     at the top of the input file as:
c     
c     %define COMPAT_UCB
c     %define NEW_CC
c  
c

c       
c     COMPAT_UCB :  Revert to the standard circa 1990 UCB values. Note
c                that these correspond to the first released version
c                of Jose Luis Martins code, not to the old Froyen
c                version (that would be 'froyen' --to be implemented)
c     (The default is: to use a denser grid up to larger radii.
c                Use a larger value for the ps's ecuts.
c                Use the Soler-Balbas XC package)
c
c     NEW_CC     : New core-correction scheme
c     OLD_CC     : Old core-correction scheme  (see wrapup.f)
c
c     The default is to use the new CC scheme only for GGA calculations.
c
c
c     For compatibility with an interim scheme using strings in the
c     input file, this routine accepts an argument "str". If not
c     empty, the user is warned that those strings now carry no
c     weight.
c

      character*(*) str

      include 'compat.h'

      logical leqi, defined
      external leqi, defined
c
      if (str .ne. " ") then
         write(6,'(a,a20)')
     $        '** WARNING: Compatibility string obsolete: ',
     $        str
         stop 'COMPAT'
      endif

      if (defined('COMPAT_UCB')) then

         write(6,'(a)') '*** UCB compatibility mode ***'
         aa_def = 6.d0
         bb_def = 40.d0
         rmax_def = 80.d0
         ecuts = 1.2d-4
         use_excorr = .true.
c
      else

         aa_def = 6.d0
         bb_def = 80.d0
         rmax_def = 120.d0 
         ecuts = 1.d-3
         use_excorr = .false.

      endif
c
c     Flag to use the old excorr subroutine.
c
      if (defined('USE_OLD_EXCORR')) use_excorr = .true.
c
c     Avoid cutting off ionic unscreened pseudopotentials
c
      if (defined('NO_PS_CUTOFFS')) ecuts = 0.d0
c
      use_old_cc = defined('OLD_CC')
      use_new_cc = defined('NEW_CC')

      end

