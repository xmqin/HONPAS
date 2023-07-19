c
      subroutine Header
c
c     Writes out the header information
c     read by subroutine Input.
c
c     Most of the information is in
c     the standard program common blocks.
c     A few things in the input.h blocks.
c
c     Alberto Garcia, July 8, 2002
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'compat.h'
      include 'input.h'
      include 'version.h'
c
      character*3 name
      integer i
      double precision xji, zero
      
      parameter(zero = 0.d0)
c
c   find jobname and date and printout.
c
      ray(1) = atom_id
      call cal_date(ray(2))
c
c   printout
c
      write(6,9090) ray(1), ray(2), title, '&v&d'
 9090 format(1x,a10,a10,1x,5a10,1x,a4/60('-'),/)
      if (job .eq. 0) then
         write(6,9100) nameat
      else if (job .lt. 4) then
         write(6,9110) nameat
      else if (job .eq. 4) then
         write(6,9120) nameat
      else if (job .eq. 5) then
         write(6,9130) nameat
      end if
 9100 format(1x,a2,' all electron calculation ',/1x,27('-'),/)
 9110 format(1x,a2,' pseudopotential generation',/1x,29('-'),/)
 9120 format(1x,a2,' pseudopotential test',/1x,23('-'),/)
 9130 format(1x,a2,' pseudo test + charge mod ',/1x,27('-'),/)
      if (ispp .eq. 'r') then
         write(6,9140)
 9140    format(' r e l a t i v i s t i c ! !',/)
         name = '   '
      else if (ispp .eq. ' ') then
         name = 'non'
      else
         name = '   '
      end if
      write(6,9150) icorr, name
 9150 format(' correlation = ',a2,3x,a3,'spin-polarized',/)
      write(6,9160) znuc, ncore, nval, zel, zion
 9160 format(' nuclear charge             =',f10.6,
     &      /' number of core orbitals    =',i3,
     &      /' number of valence orbitals =',i3,
     &      /' electronic charge          =',f10.6,
     &      /' ionic charge               =',f10.6,//)
      if (zsh .gt. 0.00001D0) write(6,9170) zsh, rsh
 9170 format(' shell charge =',f6.2,' at radius =',f6.2,//)
      write(6,9180)
 9180 format(' input data for orbitals',
     &      //'  i    n    l    s     j     occ',/)
      xji = zero
      do 140 i = 1, norb
         if (ispp .eq. 'r') xji = lo(i) + so(i)
         write(6,9190) i, no(i), lo(i), so(i), xji, zo(i)
 9190    format(1x,i2,2i5,2f6.1,f10.4)
  140 continue
      if (job .lt. 4) write(6,9200) r(2), nr, r(nr), aa, bb
 9200 format(//' radial grid parameters',//' r(1) = .0 , r(2) =',e9.3,
     &      ' , ... , r(',i4,') =',f8.3,/' a =',f7.3,'  b =',f8.3,/)
c
      return
c
      end



