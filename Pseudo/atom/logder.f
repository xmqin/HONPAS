c
      subroutine logder(nfirst,nlast,flag)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'elecpot.h'
      include 'energy.h'
c
      include 'ode_blk.h'
      include 'ode_path.h'
c
c     Computes the logarithmic derivative (LD) as a function of energy.
c     The all-electron and pseudopotential results can be compared
c     to assess the transferability of the pseudopotential.
c
c     Note:  Non-relativistic wave functions are used.
c
c     Alberto Garcia, April 27, 1991
c
c-----
c     October 1995: Removed output to unit 6.
c
c     Number of energy values at which the LD is computed and
c     half energy range around the eigenvalue (~2 rydberg)
c
      integer npt_energ
      double precision half_range
      parameter (npt_energ=51,half_range=2.d0)
      double precision eps
      parameter (eps=1.d-7)
c
c     Orbitals to be processed: nfirst <= i <= nlast
c
      integer nfirst, nlast
c
c     Flag: 'AE' for all-electron
c           'PS' for pseudo
c
      character flag*2
c
      double precision y(2)
c
      integer i, j, jr0, k, lp, llp, nok, nbad
      double precision emin, emax, eigv, step, r0, ld, h1, hmin
      character filename*9
c
      external ode, rkqc
c
c    Ode integration operational parameters
c
      h1 = 0.1d0
      hmin = 0.d0
      kmax = 200
      dxsav = 0.1d0
c
ccccc      write(6,'(/,1x,a,1x,a2,/)') 'Logarithmic derivatives', flag
c
      do 100 i = nfirst, nlast
c
        eigv = ev(i)
        emin = eigv - half_range
        emax = eigv + half_range
        step = 2 * half_range / (npt_energ - 1)
c
         lp = lo(i) + 1
         llp = lo(i)*lp
c
         write(filename,9900) flag, 'LOGD', lo(i)
 9900    format(a2,a4,i1)
         open(unit=3,file=filename,form='formatted',status='unknown')
c
c        Radius at which the LD will be calculated. Fix it at
c        a grid point. (Numerical Recipes, p.90)
c
         r0 = logder_radius
         call locate(r,nr,r0,jr0)
         r0 = r(jr0)
c
ccccccc         write(6,9000) i, no(i), lo(i), so(i), lo(i)+so(i), r0
 9000    format(/,1x,i2,2i5,2f6.1,2x,'r0:',f5.3,/)
c
         if (down(i)) then
            do 50 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j) + llp/r(j)**2
   50       continue
         else
            do 60 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j) + llp/r(j)**2
   60       continue
         end if
c
        do 30 k = 0, npt_energ - 1
c
          energ = emin + step * k
c
c         Call the integration routine. Adaptive stepsize.
c
c         y(1): wavefunction   ( = r R )
c         y(2): first derivative
c
c         Initial values (at r(2), near r=0 and so R ~ r^l )
c
          y(1) = r(2) ** lp
          y(2) = lp * r(2)**(lp-1)
c
          call odeint(y,2,r(2),r0,eps,h1,hmin,nok,nbad,ode,rkqc)
c
c         Write solution to file 20
c
c$$$          do 200 j = 1, kmax
c$$$             write(20,8000)  xp(j),(yp(m,j),m=1,2)
c$$$  200     continue
c$$$ 8000     format(1x,f10.6,4x,2(2x,f14.6))
c$$$c 
c$$$          write(20,'(/,1x,a,/)') '---------------------'
c
          ld = y(2) / y(1)
ccccccccc          write(6,9100) energ, y(2), y(1), ld, nok, nbad
          write(3,8000) energ, ld, eigv
 8000     format(1x,f12.6,3x,f12.6,3x,f12.6)
 9100     format(1x,f12.6,3x,2(2x,f10.6),3x,f10.6,2(2x,i4))
c
   30   continue
c
          close(3)
c
  100 continue  
c
      return
c
      end

