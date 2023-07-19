       subroutine blypxc(nspin,dens,gdens,EX,EC,
     .                   dEXdd,dECdd,dEXdgd,dECdgd) 
c ***************************************************************
c Implements Becke gradient exchange functional (A.D. 
c Becke, Phys. Rev. A 38, 3098 (1988)) and Lee, Yang, Parr
c correlation functional (C. Lee, W. Yang, R.G. Parr, Phys. Rev. B
c 37, 785 (1988)), as modificated by Miehlich,Savin,Stoll and Preuss,
c Chem. Phys. Lett. 157,200 (1989). See also Johnson, Gill and Pople,
c J. Chem. Phys. 98, 5612 (1993). Some errors were detected in this
c last paper, so not all of the expressions correspond exactly to those
c implemented here.
c Written by Maider Machado. July 1998.
c **************** INPUT ******************************************** 
c integer nspin          : Number of spin polarizations (1 or 2)
c real*8  dens(nspin)    : Total electron density (if nspin=1) or
c                           spin electron density (if nspin=2)
c real*8  gdens(3,nspin) : Total or spin density gradient
c ******** OUTPUT *****************************************************
c real*8  ex             : Exchange energy density
c real*8  ec             : Correlation energy density
c real*8  dexdd(nspin)   : Partial derivative
c                           d(DensTot*Ex)/dDens(ispin),
c                           where DensTot = Sum_ispin( DENS(ispin) )
c                          For a constant density, this is the
c                          exchange potential
c real*8  decdd(nspin)   : Partial derivative
c                           d(DensTot*Ec)/dDens(ispin),
c                           where DensTot = Sum_ispin( DENS(ispin) )
c                          For a constant density, this is the
c                          correlation potential
c real*8  dexdgd(3,nspin): Partial derivative
c                           d(DensTot*Ex)/d(GradDens(i,ispin))
c real*8  decdgd(3,nspin): Partial derivative
c                           d(DensTot*Ec)/d(GradDens(i,ispin))
c ********* UNITS ****************************************************
c Lengths in Bohr
c Densities in electrons per Bohr**3
c Energies in Hartrees
c Gradient vectors in cartesian coordinates
c ********************************************************************
 
      implicit none
      integer nspin
      double precision  dens(nspin), gdens(3,nspin), EX, EC,
     .                  dEXdd(nspin), dECdd(nspin), dEXdgd(3,nspin),
     .                  dECdgd(3,nspin)

c Internal variables
      integer is,ix,ois
      double precision pi, beta, thd, tthd, thrhlf, half, fothd,
     .                 d(2),gd(3,2),dmin, ash,gdm(2),denmin,dt, 
     .                 g(2),x(2),a,b,c,dd,onzthd,gdmin, 	     
     .                 ga, gb, gc,becke,dbecgd(3,2),
     .                 dgdx(2), dgdxa, dgdxb, dgdxc,dgdxd,dbecdd(2),
     .                 den,omega, domega, delta, ddelta,cf,
     .                 gam11, gam12, gam22, LYPa, LYPb1,
     .                 LYPb2,dLYP11,dLYP12,dLYP22,LYP,
     .                 dd1g11,dd1g12,dd1g22,dd2g12,dd2g11,dd2g22,
     .                 dLYPdd(2),dg11dd(3,2),dg22dd(3,2),
     .                 dg12dd(3,2),dLYPgd(3,2)
  
c Lower bounds of density and its gradient to avoid divisions by zero
      parameter ( denmin=1.d-8 )
      parameter (gdmin=1.d-8)
      parameter (dmin=1.d-5)

c Fix some numerical parameters 
      parameter ( thd = 1.d0/3.d0, tthd=2.d0/3.d0 )
      parameter ( thrhlf=1.5d0, half=0.5d0,
     .            fothd=4.d0/3.d0, onzthd=11.d0/3.d0)

c Empirical parameter for Becke exchange functional (a.u.)
      parameter(beta= 0.0042d0) 

c Constants for LYP functional (a.u.) 
      parameter(a=0.04918d0, b=0.132d0, c=0.2533d0, dd=0.349d0)

       pi= 4*atan(1.d0)
       

c Translate density and its gradient to new variables
      if (nspin .eq. 1) then
        d(1) = half * dens(1)
        d(1) = max(denmin,d(1))
        d(2) = d(1)
        dt = max( denmin, dens(1) )
        do ix = 1,3
          gd(ix,1) = half * gdens(ix,1)    
          gd(ix,2) = gd(ix,1)
        enddo 
      else
        d(1) = dens(1)
        d(2) = dens(2)
        do is=1,2
         d(is) = max (denmin,d(is))
        enddo
        dt = max( denmin, dens(1)+dens(2) )  
        do ix = 1,3
          gd(ix,1) = gdens(ix,1)
          gd(ix,2) = gdens(ix,2)
        enddo
      endif

      gdm(1) = sqrt( gd(1,1)**2 + gd(2,1)**2 + gd(3,1)**2 )
      gdm(2) = sqrt( gd(1,2)**2 + gd(2,2)**2 + gd(3,2)**2 )
 
      do is=1,2
      gdm(is)= max(gdm(is),gdmin)
      enddo

c Find Becke exchange energy
       ga = -thrhlf*(3.d0/4.d0/pi)**thd
      do is=1,2
       if(d(is).lt.dmin) then
        g(is)=ga
       else
        x(is) = gdm(is)/d(is)**fothd
        gb = beta*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1)) 
        gc = 1+6*beta*x(is)*ash        
        g(is) = ga-gb/gc
       endif
      enddo

c   Density of energy 
      becke=(g(1)*d(1)**fothd+g(2)*d(2)**fothd)/dt

      
c Exchange energy derivatives
       do is=1,2
        if(d(is).lt.dmin)then
         dbecdd(is)=0.
         do ix=1,3
          dbecgd(ix,is)=0.
         enddo
        else
        dgdxa=6*beta**2*x(is)**2
        ash=log(x(is)+sqrt(x(is)**2+1))
        dgdxb=x(is)/sqrt(x(is)**2+1)-ash
        dgdxc=-2*beta*x(is)
        dgdxd=(1+6*beta*x(is)*ash)**2
        dgdx(is)=(dgdxa*dgdxb+dgdxc)/dgdxd
        dbecdd(is)=fothd*d(is)**thd*(g(is)-x(is)*dgdx(is))
        do ix=1,3
         dbecgd(ix,is)=d(is)**(-fothd)*dgdx(is)*gd(ix,is)/x(is)
        enddo 
        endif
       enddo

c  Lee-Yang-Parr correlation energy
      den=1+dd*dt**(-thd)
      omega=dt**(-onzthd)*exp(-c*dt**(-thd))/den
      delta=c*dt**(-thd)+dd*dt**(-thd)/den
      cf=3.*(3*pi**2)**tthd/10.
      gam11=gdm(1)**2
      gam12=gd(1,1)*gd(1,2)+gd(2,1)*gd(2,2)+gd(3,1)*gd(3,2)
      gam22=gdm(2)**2
      LYPa=-4*a*d(1)*d(2)/(den*dt)
      LYPb1=2**onzthd*cf*a*b*omega*d(1)*d(2)
      LYPb2=d(1)**(8./3.)+d(2)**(8./3.)
      dLYP11=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)
     .*d(1)/dt)-d(2)**2)
      dLYP12=-a*b*omega*(d(1)*d(2)/9.*(47.-7.*delta)
     .-fothd*dt**2)
      dLYP22=-a*b*omega*(d(1)*d(2)/9.*(1.-3.*delta-(delta-11.)*
     .d(2)/dt)-d(1)**2)

c    Density of energy
      LYP=(LYPa-LYPb1*LYPb2+dLYP11*gam11+dLYP12*gam12
     .+dLYP22*gam22)/dt

c   Correlation energy derivatives
       domega=-thd*dt**(-fothd)*omega*(11.*dt**thd-c-dd/den)
       ddelta=thd*(dd**2*dt**(-5./3.)/den**2-delta/dt)

c   Second derivatives with respect to the density
       dd1g11=domega/omega*dLYP11-a*b*omega*(d(2)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2))

       dd1g12=domega/omega*dLYP12-a*b*omega*(d(2)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)

      dd1g22=domega/omega*dLYP22-a*b*omega*(1./9.*d(2)
     . *(1.-3.*delta-(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3.+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2)-2*d(1))

       
      dd2g22=domega/omega*dLYP22-a*b*omega*(d(1)/9.*
     . (1.-3.*delta-2*(delta-11.)*d(2)/dt)-d(1)*d(2)/9.*
     . ((3+d(2)/dt)*ddelta-(delta-11.)*d(2)/dt**2))
      
 
      dd2g12=domega/omega*dLYP12-a*b*omega*(d(1)/9.*
     . (47.-7.*delta)-7./9.*d(1)*d(2)*ddelta-8./3.*dt)
      
      dd2g11=domega/omega*dLYP11-a*b*omega*(1./9.*d(1)
     . *(1.-3.*delta-(delta-11.)*d(1)/dt)-d(1)*d(2)/9.*
     . ((3.+d(1)/dt)*ddelta-(delta-11.)*d(1)/dt**2)-2*d(2))


        dLYPdd(1)=-4*a/den*d(1)*d(2)/dt*
     . (thd*dd*dt**(-fothd)/den
     . +1./d(1)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(2)*(onzthd*
     . d(1)**(8./3.)+d(2)**(8./3.)))+dd1g11*gam11+
     . dd1g12*gam12+dd1g22*gam22


       dLYPdd(2)=-4*a/den*d(1)*d(2)/dt*(thd*dd*dt**(-fothd)/den
     . +1./d(2)-1./dt)-2**onzthd*cf*a*b*(domega*d(1)*d(2)*
     . (d(1)**(8./3.)+d(2)**(8./3.))+omega*d(1)*(onzthd*
     . d(2)**(8./3.)+d(1)**(8./3.)))+dd2g22*gam22+
     . dd2g12*gam12+dd2g11*gam11


c   Second derivatives with respect to the density gradient

        do is=1,2
          do ix=1,3
           dg11dd(ix,is)=2*gd(ix,is)
           dg22dd(ix,is)=2*gd(ix,is)
          enddo
        enddo
        do ix=1,3
          dLYPgd(ix,1)=dLYP11*dg11dd(ix,1)+dLYP12*gd(ix,2)
          dLYPgd(ix,2)=dLYP22*dg22dd(ix,2)+dLYP12*gd(ix,1)
        enddo




c    Set output arguments
       EX=becke
       EC=LYP
       do is=1,nspin
        dEXdd(is)=dbecdd(is)
        dECdd(is)=dLYPdd(is)
        do ix=1,3
         dEXdgd(ix,is)=dbecgd(ix,is)
         dECdgd(ix,is)=dLYPgd(ix,is)
        enddo
       enddo
       end 


