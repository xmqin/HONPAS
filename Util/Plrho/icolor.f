! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      integer function icolor( bright, ncolor, color )

c **********************************************************************
c Finds a color index previously set by routine redblue.
c Written by J.M.Soler. Nov'97.
c ********** Input *****************************************************
c real    bright : Color intensity in the range (0.,1.)
c integer ncolor : Number of color values 
c real    color  : Value of function(s) which determines color
c **********************************************************************

      implicit none
      integer  ncolor
      real     b, bright, color

      integer  icb, icmax, icmin, icr, iseed, nc
      real     cb, colblu, colred, colwht, cr, dcb, dcr, ran1
      common /comcol/ colred, colwht, colblu, icmin, nc
      data iseed /-6176601/

      if (ncolor .eq. 0) then
        b = max( bright, 0. )
        b = min( b, 1. )
        icolor = icmin + b * nc
      elseif (ncolor .eq. 1) then
        if (color .lt. colwht) then
          cr = 1.
          cb = 1. - min( 1., abs(color-colwht)/abs(colred-colwht) )
        else
          cb = 1.
          cr = 1. - min( 1., (color-colwht)/(colblu-colwht) )
        endif
        dcr = nc * cr * bright * 0.999999
        dcb = nc * cb * bright * 0.999999
        icr = int(dcr)
        icb = int(dcb)
        dcr = dcr - icr
        dcb = dcb - icb
        if (ran1(iseed) .lt. dcr) icr = icr + 1
        if (ran1(iseed) .lt. dcb) icb = icb + 1
        icolor = icmin + icr + (nc+1) * icb
      else
        stop 'icolor: ERROR: ncolor value not allowed'
      endif
      end


      integer function icback( ncolor )

c **********************************************************************
c Returnrs background color
c Written by J.M.Soler. Nov'97.
c **********************************************************************

      implicit none
      integer  ncolor

      integer  icmin, nc
      real     colblu, colred, colwht
      common /comcol/ colred, colwht, colblu, icmin, nc
 
      if (ncolor .eq. 0) then
        icback = icmin - 1 + nc
      else
        icback = icmin - 1 + (nc+1)**2
      endif

*     icback = 0
      end

      subroutine grays()

c **********************************************************************
c Sets up a grayscale color table using pgplot library.
c Written by J.M.Soler. Nov'97.
c **********************************************************************

      implicit none
      integer  ic, icmax, icmin, nc
      real     b, colblu, colred, colwht
      common /comcol/ colred, colwht, colblu, icmin, nc

      call pgqcir( icmin, icmax )
      nc = icmax - icmin + 1
      write(6,*) 'grays: color index range =', icmin, icmax
      do ic = 1,nc
        b = real(ic) / nc
        call pgscr( icmin+ic, b, b, b )
      enddo
      end



      subroutine redblu( colmin, col0, colmax )

c **********************************************************************
c Sets up a red->white->blue color table using pgplot library.
c Written by J.M.Soler. Nov'97.
c ********** Input *****************************************************
c real colmin, colmax : Range of the function to represent by color.
c                       colmin will correspond to pure red
c                       colmax will correspond to pure blue
c                       zero will correspont to pure white
c **********************************************************************

      implicit none
      real     colmin, colmax, col0
      external pgscr

c Internal variables
      integer  ic, icb, icmax, icmin, icr, icrang, nc
      real     blue, colblu, colred, colwht, green, red
      common /comcol/ colred, colwht, colblu, icmin, nc

c Copy color range to common variables
      colred = colmin
      colblu = colmax
      colwht = col0

c Find the range of color indexes (number of simultaneous colors
c                                  accepted by the screen used)
      call pgqcir( icmin, icmax )
      icrang = icmax - icmin + 1
*     write(6,*) 'redblu: color index range =', icmin, icmax

c Find the number of colors for red and/or blue
      nc = sqrt( real(icrang) ) - 1

c Set color index table
      ic = icmin - 1
      do icb = 0,nc
        do icr = 0,nc
          ic = ic + 1
          red  = real(icr) / real(nc)
          blue = real(icb) / real(nc)
          green = min( red, blue )
          call pgscr( ic, red, green, blue )
        enddo
      enddo
      end




      FUNCTION RAN1(IDUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END


