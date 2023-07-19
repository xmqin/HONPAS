c
c $Id: excorr.f,v 1.3 1999/02/26 14:26:43 wdpgaara Exp $
c
      subroutine excorr(id,cdd,cdu,cdc,vod,vou,vxc,vc,exc,ec)
c
      implicit none
c
      include 'radial.h'
      include 'param.h'
c
c     Compute the LDA exchange and correlation energy and potential.
c
c     Revised by Alberto Garcia
c
c    The only major modification is that the constants for the
c    ceperly-alder 'ca' method are placed in parameter
c    statements, this was done so non-opt compiliers
c    would minimize the number of calculations.
c
C     .. Parameters ..
c
      double precision tiny_charge
      parameter (tiny_charge=1.d-12)
c
      double precision zero, one, pfive, opf, pnn
      parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0,pnn=.99D0)
      double precision pthree, psevf, c0504
      parameter (pthree=0.3D0,psevf=0.75D0,c0504=0.0504D0)
      double precision c0254, c014, c0406
      parameter (c0254=0.0254D0,c014=0.014D0,c0406=0.0406D0)
      double precision c15p9, c0666, c11p4
      parameter (c15p9=15.9D0,c0666=0.0666D0,c11p4=11.4D0)
      double precision c045, c7p8, c88, c20p592
      parameter (c045=0.045D0,c7p8=7.8D0,c88=0.88D0,c20p592=20.592D0)
      double precision c3p52, c0311, c0014
      parameter (c3p52=3.52D0,c0311=0.0311D0,c0014=0.0014D0)
      double precision c0538, c0096, c096
      parameter (c0538=0.0538D0,c0096=0.0096D0,c096=0.096D0)
      double precision c0622, c004, c0232
      parameter (c0622=0.0622D0,c004=0.004D0,c0232=0.0232D0)
      double precision c1686, c1p3981, c2611
      parameter (c1686=0.1686D0,c1p3981=1.3981D0,c2611=0.2611D0)
      double precision c2846, c1p0529, c3334
      parameter (c2846=0.2846D0,c1p0529=1.0529D0,c3334=0.3334D0)
      double precision con1, con2, con3
      parameter (con1=1.D0/6,con2=0.008D0/3,con3=0.3502D0/3)
      double precision con4, con5, con6
      parameter (con4=0.0504D0/3,con5=0.0028D0/3,con6=0.1925D0/3)
      double precision con7, con8, con9
      parameter (con7=0.0206D0/3,con8=9.7867D0/6,con9=1.0444D0/3)
      double precision con10, con11
      parameter (con10=7.3703D0/6,con11=1.3336D0/3)
C     ..
C     .. Scalar Arguments ..
      double precision vxc, vc, exc, ec
      character id*1
C     ..
C     .. Array Arguments ..
      double precision cdd(nrmax), cdu(nrmax), cdc(nrmax), 
     &                 vod(nrmax), vou(nrmax)
c
C     .. Local Scalars ..
      double precision a0, alb, aln, alp,  be, beta,
     &                 cdsum, ecf, ecp, ect, excf, excp,
     &                 exct, exf, exp_var, ftrd, fz, fzp, pi, rs,
     &                 rslog, sb, sqrs, te, tftm, trd, vcd, vcf,
     &                 vcp, vcu, vxcd, vxcf, vxcp, vxcu, vxf, 
     &                 lda_xpot, x, z
      integer i, ll
C     ..
C     .. External Subroutines ..
      external ext
C     ..
C     .. Intrinsic Functions ..
      intrinsic atan, log, sqrt
C     ..
      logical leqi
      external leqi
C     ..
c
      pi = 4*atan(one)
c
c
      trd = one/3
      ftrd = 4*trd
      tftm = 2**ftrd - 2
      a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
      alp = one
      if (.not. leqi(icorr,'xa')) alp = 2*trd
c
c     Initialize
c
      vxc = zero
      vc = zero
      exc = zero
      ec = zero
c
c      start loop (at the second point...see below)
c
      ll = 4
      do 70 i = 2, nr
c
         cdsum = cdd(i) + cdu(i)
         if (ifcore .ge. 1) cdsum = cdsum + cdc(i)
c
         vxcd = 0.d0
         vxcu = 0.d0
         vcd = 0.d0
         vcu = 0.d0
         exct = 0.d0
         ect = 0.d0
c
cag****!!!!!!         if (cdsum .le. tiny_charge) go to 100
c
         rs = (3*r(i)**2/cdsum)**trd
c
c        Spin variables
c
         z = zero
         fz = zero
         fzp = zero
         if (leqi(id,'s')) then
            z = (cdd(i)-cdu(i))/cdsum
            fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
            fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
         end if
c
c      exchange (only use (xa))
c
         lda_xpot = -3*alp/(pi*a0*rs)
         exp_var = 3*lda_xpot/4
c
c        Relativistic correction to exchange
c
         if (leqi(id,'r')) then
            beta = c014/rs
            sb = sqrt(1+beta*beta)
            alb = log(beta+sb)
            lda_xpot = lda_xpot*(-pfive+opf*alb/(beta*sb))
            exp_var = exp_var*(one-opf*((beta*sb-alb)/beta**2)**2)
         end if
c
   60    continue
c
         vxf = 2**trd*lda_xpot
         exf = 2**trd*exp_var
         vcp = zero
         ecp = zero
         vcf = zero
         ecf = zero
c
         if (leqi(icorr,'ca')) then
c          ceperly-alder (ca)
c          The Perdew-Zunger parameterization is used.
c          See Phys. Rev. B 23 5075 (1981).
            if (rs .gt. one) then
               sqrs = sqrt(rs)
               te = one + con10*sqrs + con11*rs
               be = one + c1p0529*sqrs + c3334*rs
               ecp = -c2846/be
               vcp = ecp*te/be
               te = one + con8*sqrs + con9*rs
               be = one + c1p3981*sqrs + c2611*rs
               ecf = -c1686/be
               vcf = ecf*te/be
            else
               rslog = log(rs)
               ecp = (c0622+c004*rs)*rslog - c096 - c0232*rs
               vcp = (c0622+con2*rs)*rslog - con3 - con4*rs
               ecf = (c0311+c0014*rs)*rslog - c0538 - c0096*rs
               vcf = (c0311+con5*rs)*rslog - con6 - con7*rs
            end if
c
         else if (leqi(icorr,'xa')) then
c
c          correlation
c
         else if (leqi(icorr,'wi')) then
c
c          wigner (wi)
            vcp = -(c3p52*rs+c20p592)/(3*(rs+c7p8)**2)
            ecp = -c88/(rs+c7p8)
c
         else if (leqi(icorr,'hl')) then
c          hedin-lundqvist (hl)
            x = rs/21
            aln = log(1+1/x)
            vcp = -c045*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c045*ecp
c
         else if (leqi(icorr,'gl')) then
c          gunnarson-lundqvist-wilkins (gl)
            x = rs/c11p4
            aln = log(1+1/x)
            vcp = -c0666*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0666*ecp
            x = rs/c15p9
            aln = log(1+1/x)
            vcf = -c0406*aln
            ecf = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecf = ((con1/x-pthree)/x+psevf)/x
            ecf = -c0406*ecf
c
         else if (leqi(icorr,'bh')) then
c          von barth - hedin (bh)
            x = rs/30
            aln = log(1+1/x)
            vcp = -c0504*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0504*ecp
            x = rs/75
            aln = log(1+1/x)
            vcf = -c0254*aln
            ecf = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecf = ((con1/x-pthree)/x+psevf)/x
            ecf = -c0254*ecf
c
         else if (leqi(icorr,'gr')) then
c          von barth - hedin + gradient corrections (gr)
            x = rs/30
            aln = log(1+1/x)
            vcp = -c0504*aln
            ecp = aln + (x**3*aln-x*x) + x/2 - trd
            if (x .gt. 500*one) ecp = ((con1/x-pthree)/x+psevf)/x
            ecp = -c0504*ecp
c
         else
c
            write(6,9050) icorr
 9050       format('error in velect - icorr =',a2,' not implemented')
            call ext(400)
c
         end if
c
         vxcp = lda_xpot + vcp
         vxcf = vxf + vcf
         vxcd = vxcp
         vxcu = vxcp
         excp = exp_var + ecp
         excf = exf + ecf
         vcd = vcp
         vcu = vcp
         exct = excp
         ect = ecp
c
         if (z .ne. zero) then
            vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
            vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
            vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
            vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
            exct = exct + fz*(excf-excp)
            ect = ect + fz*(ecf-ecp)
         end if
c
 100     continue
c
         vod(i) = vod(i) + vxcd
         vou(i) = vou(i) + vxcu
c
c        Add to the integrated value
c
         vxc = vxc + ll*(cdd(i)*vxcd+cdu(i)*vxcu)*rab(i)
         vc = vc + ll*(cdd(i)*vcd+cdu(i)*vcu)*rab(i)
         exc = exc + ll*cdsum*exct*rab(i)
         ec = ec + ll*cdsum*ect*rab(i)
         ll = 6 - ll
c
   70 continue
c
      vxc = vxc/3
      vc = vc/3
      ec = ec/3
      exc = exc/3
c
c     Extrapolate backwards for the first point...
c
      vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
      vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))
c
      return
c
      end
