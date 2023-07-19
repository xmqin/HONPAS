c
c compat.h
c
c     Variables and flags to modify the behavior of the program without
c     explicit reference in the input file.
c
      real*8 aa_def, bb_def, rmax_def, ecuts
      logical use_excorr, use_old_cc, use_new_cc
c
      common /rcompat/ aa_def, bb_def, rmax_def, ecuts
      common /lcompat/ use_excorr, use_old_cc, use_new_cc
      save /rcompat/, /lcompat/
c---------

