! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      real function ray( g )
c *******************************************************************
c Returns the light intensity (brightness) arriving from a 
c surface pixel with a given orientation.
c Written by J.M.Soler, April 1998
c ********** Input **************************************************
c real g(3) : Vector perpendicular to the surface, pointing 
c             inwards, but not necessarily normalized
c *******************************************************************
C    6  10        20        30        40        50        60        7072

      implicit none

c maxl : Maximun number of light sources
      integer maxl
      real    eps
      parameter ( maxl = 10 )
      parameter ( eps  = 1.e-12 )

c Internal variables
      logical  frstme, yes
      integer  il, nl
      real     ambient, cl(maxl), cosnl, dx2, diff,
     .         dl(maxl), dlight(maxl), g(3), gmod,
     .         pi, plight(maxl), refl, reflect,
     .         tlight(maxl), wlight(maxl), wl2(maxl),
     .         xi(3), xl(3,maxl), xn(3)
      save     ambient, dl, dlight, frstme, nl, plight, reflect, 
     .         tlight, wlight, wl2, xl
      data     frstme /.true./

c Default illumination data
c reflect : Surface reflectivity (between 0 and 1)
c ambient : Ambient light intensity
c nl      : Number of directional light sources
c dlight  : Intensity of directional light source
c tlight  : Polar angle (theta) between the light source and viewpoint
c plight  : Azimutal angle (phi) between the light source and viewpoint
c wlight  : Angular light-width (all angles in degrres)
      data 
     . reflect / 0.3 /,
     . ambient / 0.2 /,
     . nl      / 1   /,
     . (dlight(il),wlight(il),tlight(il),plight(il),il=1,1) /
     .   0.8, 10., 30., 90. /

c Read illumination data
      if (frstme) then
        inquire( file='light.dat', exist=yes )
        if (yes) then
          open(1,file='light.dat',status='old')
          read(1,*) reflect
          read(1,*) ambient
          read(1,*) nl
          if (nl .gt. maxl) stop 'ray: dimension maxl too small'
          do il = 1,nl
            read(1,*) dlight(il), wlight(il), tlight(il), plight(il)
          enddo
          close(1)
        endif
      endif

c Find illumination vector
      if (frstme) then
        pi = 4. * atan(1.)
        do il = 1,nl
          xl(1,il) = - sin(tlight(il)*pi/180.) * cos(plight(il)*pi/180.)
          xl(2,il) = - sin(tlight(il)*pi/180.) * sin(plight(il)*pi/180.)
          xl(3,il) = - cos(tlight(il)*pi/180.)
          wl2(il) = (wlight(il)*pi/180.)**2 + eps
          dl(il) = dlight(il) / pi / wl2(il)
        enddo
      endif

c Add light contributions
      ray = ambient
      do il = 1,nl

c       Find normalized surface normal
        gmod = sqrt( g(1)**2 + g(2)**2 + g(3)**2 ) + eps
        xn(1) = g(1) / gmod
        xn(2) = g(2) / gmod
        xn(3) = g(3) / gmod

c       Find difuse (scattered) light intensity
        cosnl = xn(1)*xl(1,il) + xn(2)*xl(2,il) + xn(3)*xl(3,il)
        cosnl = max( cosnl, 0. )
        diff  = (1.-reflect) * dlight(il) * cosnl

c       Find reflected light intensity
c       Find incident direction, which reflects towards (0,0,1)
        xi(1) = -2.*xn(3)*xn(1)
        xi(2) = -2.*xn(3)*xn(2)
        xi(3) = -2.*xn(3)*xn(3) + 1.
        dx2 = (xi(1)-xl(1,il))**2 +
     .        (xi(2)-xl(2,il))**2 +
     .        (xi(3)-xl(3,il))**2
        refl = reflect * dl(il) * exp(-dx2/wl2(il))
      
        ray = ray + diff + refl
      enddo

c Saturate ray intensity to one
      ray = min( ray, 0.999999 )

      frstme = .false.
      end



