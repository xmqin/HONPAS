c---
c     Variables for communication between input and header
c     (Old input was split: the new input does only
c      reading, and header does the writing).
c
      double precision aa, bb, zion
      character ray(5)*10, title(5)*10
      integer nval
      common /input_s/ aa, bb, zion, nval
      common /input_chr/ ray, title
c
