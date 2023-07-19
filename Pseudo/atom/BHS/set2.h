c---
c    Set2.h - Implementation of a set of pairs
c
      integer nmax
      parameter (nmax=50)
c
      integer nels
      character*30 el_name(nmax)
      double precision el_val(nmax)
c
      common /set_num/ el_val, nels
      common /set_chr/ el_name
c
      save /set_num/, /set_chr/
c---
