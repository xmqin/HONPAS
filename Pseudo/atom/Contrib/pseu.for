c
      subroutine pseudo(itype,icorr,ispp,nr,a,b,r,rab,nameat,norb,ncore,
     &                  no,lo,so,zo,znuc,zel,cdd,cdu,cdc,viod,viou,vid,
     &                  viu,vod,vou,etot,ev,ek,ep)
c
c *************************************************************
c *                                                           *
c *    pseudo generates the pseudo potential using            *
c *  the scheme of Hamann, Schluter and Chiang -              *
c *  Phys. Rev. Lett. 43, 1494 (1979).                        *
c *                                                           *
c *************************************************************
c
c  njtj  *** modifications  ***
c    The only major modifications are in the spin-polarized
c    treatment of the el-el unscreening of the pseudopotential
c    A spin-polarized pseudopotential is unscreened
c    with a spin-polarized valence charge.  This was not done
c    in pseudo or pseudok in earlier versions of this
c    program.
c  njtj  *** modifications  ***
c
C     .. Parameters ..
      double precision zero, ecuts, tpfive, one
      parameter (zero=0.D0,ecuts=2.0D-4,tpfive=2.5D0,one=1.D0)
      double precision small, small2, small3, pzfive
      parameter (small=1.D-13,small2=1.D-10,small3=1.D-18,pzfive=.05D0)
      double precision pfive, small4, ai
      parameter (pfive=0.5D0,small4=1.D-6,ai=2*137.0360411D0)
      integer lmax, norbmx, nrmax
      parameter (lmax=4,norbmx=40,nrmax=1000)
C     ..
C     .. Scalar Arguments ..
      double precision a, b, zel, znuc
      integer itype, ncore, norb, nr
      character ispp*1, icorr*2, nameat*2
C     ..
C     .. Array Arguments ..
      double precision cdc(*), cdd(*), cdu(*), ek(*), ep(*), etot(10),
     &                 ev(*), r(*), rab(*), so(*), vid(*), viod(lmax,*),
     &                 viou(lmax,*), viu(*), vod(*), vou(*), zo(*)
      integer lo(*), no(*)
C     ..
C     .. Local Scalars ..
      double precision aa, ac, ag, arp, arpm, bc, cdcp, cfac, cl, dcl,
     &                 delta, dev, devold, ecut, eviae, fcut, fjm1,
     &                 gamma, gg, gpp, pi, rbnew, rbold, rcfac, rextr,
     &                 rmind, rminu, rra, rrc, rrp, rzero, tanb, viodj,
     &                 viouj, vp2z, vps, vpsdm, vpsum, zeld, zelt, zelu,
     &                 zion, zot, zratio, zval
      integer i, icore, ifcore, iflag, ifull, ist, j, j3rc, jcut, k, ka,
     &        ll, llp, lp, ncp, nextr, noi, npotd, npotu
      character blank*1, irel*3, nicore*4
C     ..
C     .. Local Arrays ..
      double precision rc(5), rcut(10)
      integer indd(5), indu(5)
      character il(5)*1, iray(6)*10, ititle(7)*10
C     ..
C     .. External Subroutines ..
      external difnrl, difrel, dsolv2, etotal, ext, potran, potrv,
     &         potrw, velect, wtrans, zedate
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, atan, exp, sin, sqrt
C     ..
      double precision ar(nrmax), arps(nrmax), br(nrmax), f(nrmax),
     &                 g(nrmax), v(nrmax), wk1(nrmax), wk2(nrmax),
     &                 wk3(nrmax), wk4(nrmax), wk5(nrmax), wk6(nrmax),
     &                 wk7(nrmax), wkb(6*nrmax)
      integer nops(norbmx)
C     ..
C     ..
C     .. Data statements ..
c
      data il/'s', 'p', 'd', 'f', 'g'/
C     ..
      do 10 i = 1, 5
         indd(i) = 0
         indu(i) = 0
   10 continue
      if (ncore .eq. norb) return
      if (itype .ne. 1 .and. itype .ne. 2 .and. itype .ne. 3) return
      ifcore = itype - 1
      pi = 4*atan(one)
c
c  Spin-polarized potentails should be unscreened with
c  a spin-polarized valence charge.  This was not
c  done in pseudo and pseudk in earlier versions
c  of this program.
c
      if (ispp .eq. 's') then
         blank = 's'
      else
         blank = ' '
      end if
c
c  read rc(s),rc(p),rc(d),rc(f),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,9000) (rc(i),i=1,4), cfac, rcfac
 9000 format(6f10.5)
      if (cfac .eq. zero) cfac = one
c
c   Reset vod and vou to zero.  They are here used to store
c   the pseudo valence charge density.
c
      do 20 i = 1, nr
         vod(i) = zero
         vou(i) = zero
   20 continue
c
c  Print the heading.
c
      write(6,9010) nameat
 9010 format(//a2,' Pseudopotential HSC generation',/1x,35('-'),
     &      //' nl    s    eigenvalue',6x,'rc',4x,6x,'cl',9x,'gamma',7x,
     &      'delta',/)
c
c      start loop over valence orbitals
c
      ncp = ncore + 1
      do 220 i = ncp, norb
         lp = lo(i) + 1
         llp = lo(i)*lp
         if (so(i) .lt. 0.1D0) then
            if (indd(lp) .ne. 0) then
               write(6,9020) lp - 1
               call ext(800+lp)
            else
               indd(lp) = i
            end if
         else
            if (indu(lp) .ne. 0) then
               write(6,9030) lp - 1
               call ext(810+lp)
            else
               indu(lp) = i
            end if
         end if
 9020    format(//
     &         'error in pseudo - two down spin orbitals of the same ',
     &         /'angular momentum (',i1,') exist')
 9030    format(//'error in pseudo - two up spin orbitals of the same ',
     &         /'angular momentum (',i1,') exist')
c
c      find all electron wave function
c
         do 30 j = 1, nr
            ar(j) = zero
   30    continue
         if (so(i) .lt. 0.1D0) then
            do 40 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
   40       continue
         else
            do 50 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
   50       continue
         end if
         if (ispp .ne. 'r') then
            do 60 j = 2, nr
               v(j) = v(j) + llp/r(j)**2
   60       continue
         end if
         if (ispp .ne. 'r') then
            call difnrl(0,i,v,ar,br,nr,a,b,r,rab,norb,no,lo,so,znuc,
     &                  viod,viou,vid,viu,ev,iflag)
         else
            call difrel(0,i,v,ar,br,nr,a,b,r,rab,norb,no,lo,so,znuc,
     &                  viod,viou,vid,viu,ev)
         end if
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an useful number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
         ist = 1
         if (ar(nr-85) .lt. zero) ist = -1
         call potrw(ar,r,nr-85,lo(i),1,ist,rc(lp))
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum.
c
         ka = lo(i) + 1
         if (so(i) .lt. 0.1D0 .and. lo(i) .ne. 0) ka = -lo(i)
         nextr = no(i) - lo(i)
         rzero = zero
         arp = br(2)
c
         if (ispp .eq. 'r') then
            if (so(i) .lt. 0.1D0) then
               arp = ka*ar(2)/r(2) + (ev(i)-viod(lp,2)/r(2)-vid(2)+
     &               ai*ai)*br(2)/ai
            else
               arp = ka*ar(2)/r(2) + (ev(i)-viou(lp,2)/r(2)-viu(2)+
     &               ai*ai)*br(2)/ai
            end if
         end if
c
         do 70 j = 3, nr
            if (nextr .eq. 0) go to 80
            if (ar(j-1)*ar(j) .le. zero) rzero = (ar(j)*r(j-1)-
     &          ar(j-1)*r(j))/(ar(j)-ar(j-1))
            arpm = arp
            arp = br(j)
c
            if (ispp .eq. 'r') then
               if (so(i) .lt. 0.1D0) then
                  arp = ka*ar(j)/r(j) + (ev(i)-viod(lp,j)/r(j)-vid(j)+
     &                  ai*ai)*br(j)/ai
               else
                  arp = ka*ar(j)/r(j) + (ev(i)-viou(lp,j)/r(j)-viu(j)+
     &                  ai*ai)*br(j)/ai
               end if
            end if
c
            if (arp*arpm .gt. zero) go to 70
            rextr = (arp*r(j-1)-arpm*r(j))/(arp-arpm)
            nextr = nextr - 1
   70    continue
c
c  Check rc, if outside bounds reset.
c
   80    continue
         if (rzero .lt. r(2)) rzero = r(2)
         if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) go to 90
         if (rc(lp) .ge. zero) then
            write(6,9040) rc(lp), rextr
         end if
 9040    format(/'Warning, the Core radius =',f5.2,
     &         /' is larger then wave function',' extrema position =',
     &         f5.2,/)
         if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c  Reset the n quantum numbers.
c
   90    continue
         do 100 j = 1, norb
            nops(j) = 0
  100    continue
         nops(i) = lp
c
c  njtj  ***  modification start  ***
c    Sset up the functions f(r/rc) and g(r/rc) and
c  modify the ionic potential.
c
         aa = 4*one
         dcl = -6*one*lp
         cl = dcl
c
         do 110 j = 1, nr
            rrc = r(j)/rc(lp)
            rra = rrc**aa
            f(j) = zero
            if (rra .lt. 88*one) f(j) = exp(-rra)
            g(j) = rrc**lp*f(j)
            fjm1 = one - f(j)
            if (fjm1 .lt. small4) fjm1 = (one-pfive*rra)*rra
            if (so(i) .lt. 0.1D0) then
               viod(lp,j) = fjm1*viod(lp,j) - f(j)*r(j)*vid(j) +
     &                      dcl*r(j)*f(j)
            else
c
c bug fix Alberto Garcia 5/11/90
c
               viou(lp,j) = fjm1*viou(lp,j) - f(j)*r(j)*viu(j) +
     &                      dcl*r(j)*f(j)
            end if
            if (rrc .lt. 3*one) j3rc = j
  110    continue
         dcl = dcl/2
c
c   Start the iteration loop to find cl.
c
         eviae = ev(i)
         devold = zero
         do 150 j = 1, 100
            call dsolv2(j,2,blank,ifcore,nr,a,b,r,rab,norb,ncore,nops,
     &                  lo,so,zo,znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev,
     &                  ek,ep)
            dev = eviae - ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c   division by zero errors in the calculation of
c   dcl = -dev*dcl / (dev-devold).
c
            if ((abs(dev).lt.small2.or.abs(dev-devold).lt.small3) .and.
     &          j .ne. 1) then
c
               go to 160
c
            else
               if (j .gt. 20 .or. abs(dev) .lt. 0.001D0) then
c
c   Use newton raphson iteration to change cl.
c
                  dcl = -dev*dcl/(dev-devold)
               else
                  if (dev*dcl .lt. zero) then
                     dcl = -dcl/3
                  end if
               end if
            end if
c
c  njtj  ***  modification end  ***
c
c  Find the new potential.
c
  120       continue
            if (so(i) .lt. 0.1D0) then
               do 130 k = 2, nr
                  viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
  130          continue
            else
               do 140 k = 2, nr
                  viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
  140          continue
            end if
            cl = cl + dcl
            devold = dev
  150    continue
c
c  End the iteration loop for cl.
c
         call ext(820+lp)
c
c   Find the pseudo-wavefunction.
c
  160    continue
         if (so(i) .lt. 0.1D0) then
            do 170 j = 2, nr
               v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
  170       continue
         else
            do 180 j = 2, nr
               v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
  180       continue
         end if
         call difnrl(0,i,v,arps,br,nr,a,b,r,rab,norb,nops,lo,so,znuc,
     &               viod,viou,vid,viu,ev,iflag)
c
c  Compute delta and gamma.
c
         gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
         ag = zero
         gg = zero
         ll = 4
         do 190 j = 2, nr
            ag = ag + ll*arps(j)*g(j)*rab(j)
            gg = gg + ll*g(j)*g(j)*rab(j)
            ll = 6 - ll
  190    continue
         ag = ag/3
         gg = gg/3
         delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c     Modify the pseudo-wavefunction and pseudo-potential and
c   add to charge density.
c
         if (so(i) .lt. 0.1D0) then
            do 200 j = 2, nr
               arps(j) = gamma*(arps(j)+delta*g(j))
               vod(j) = vod(j) + zo(i)*arps(j)*arps(j)
               if (arps(j) .lt. small .and.
     &             r(j) .gt. one) arps(j) = small
               rrp = r(j)/rc(lp)
               gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &               r(j)**2
               viod(lp,j) = viod(lp,j) + gamma*delta*
     &                      ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  200       continue
         else
            do 210 j = 2, nr
               arps(j) = gamma*(arps(j)+delta*g(j))
               vou(j) = vou(j) + zo(i)*arps(j)*arps(j)
               if (arps(j) .lt. small .and.
     &             r(j) .gt. one) arps(j) = small
               rrp = r(j)/rc(lp)
               gpp = (llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)*g(j)/
     &               r(j)**2
               viou(lp,j) = viou(lp,j) + gamma*delta*
     &                      ((ev(i)-v(j))*g(j)+gpp)*r(j)/arps(j)
  210       continue
         end if
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a useful number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
         ist = 1
         if (arps(nr-85) .lt. zero) ist = -1
         call potrw(arps,r,nr-85,lo(i),0,ist,rc(lp))
         call wtrans(arps,r,nr,lo(i),ist,wk1,wk2,wk3)
c
c  njtj  ***  user should adjust for their needs  ***
c
         write(6,9050) nops(i), il(lp), so(i), ev(i), rc(lp), cl,
     &     gamma, delta
 9050    format(1x,i1,a1,f6.1,5f12.6)
  220 continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 230 i = ncp, norb
         nops(i) = lo(i) + 1
         zval = zval + zo(i)
  230 continue
      zion = zval + znuc - zel
      if (zval .ne. zero) zratio = zion/zval
      do 240 i = 1, nr
         cdd(i) = vod(i)
  240 continue
      do 250 i = 1, nr
         cdu(i) = vou(i)
  250 continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
         ac = zero
         bc = zero
         icore = 1
         if (cfac .le. zero .or. zratio .eq. zero) then
            write(6,9060) r(icore), ac, bc
         else
            if (rcfac .le. zero) then
               do 260 i = nr, 2, -1
                  if (cdc(i) .gt. cfac*zratio*
     &                (cdd(i)+cdu(i))) go to 280
  260          continue
            else
               do 270 i = nr, 2, -1
                  if (r(i) .le. rcfac) go to 280
  270          continue
            end if
  280       continue
            icore = i
            cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))
            tanb = cdc(icore)/(r(icore)*cdcp-cdc(icore))
            rbold = tpfive
            do 300 i = 1, 50
               rbnew = pi + atan(tanb*rbold)
               if (abs(rbnew-rbold) .lt. .00001D0) then
                  bc = rbnew/r(icore)
                  ac = cdc(icore)/(r(icore)*sin(rbnew))
                  do 290 j = 1, icore
                     cdc(j) = ac*r(j)*sin(bc*r(j))
  290             continue
                  write(6,9060) r(icore), ac, bc
c
                  go to 310
c
               else
                  rbold = rbnew
               end if
  300       continue
            write(6,9070)
            call ext(830)
         end if
      end if
 9060 format(//' core correction used',/' pseudo core inside r =',f6.3,
     &      /' ac =',f6.3,' bc =',f6.3,/)
 9070 format(//' error in pseudo - noncovergence in finding ',
     &      /'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
  310 continue
      if (ispp .eq. 's') then
         blank = 's'
      else
         blank = ' '
      end if
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,cdd,cdu,cdc,vod,
     &            vou,etot)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      write(6,9080)
 9080 format(/)
      ecut = ecuts
      do 380 i = ncp, norb
         lp = lo(i) + 1
         if (so(i) .lt. 0.1D0) then
            do 320 j = 2, nr
               viod(lp,j) = viod(lp,j) + (vid(j)-vod(j))*r(j)
               vp2z = viod(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
  320       continue
            rcut(i-ncore) = r(jcut)
            do 330 j = jcut, nr
               fcut = exp(-5*(r(j)-r(jcut)))
               viod(lp,j) = -2*zion + fcut*(viod(lp,j)+2*zion)
  330       continue
            do 340 j = 2, nr
               v(j) = viod(lp,j)/r(j)
  340       continue
c
c  njtj  ***  plotting routines ***
c
            call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
            call potrv(v,r,nr-120,lo(i),zion)
c
c  njtj  ***  user should adjust for their needs  ***
c
         else
            do 350 j = 2, nr
               viou(lp,j) = viou(lp,j) + (viu(j)-vou(j))*r(j)
               vp2z = viou(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
  350       continue
            rcut(i-ncore) = r(jcut)
            do 360 j = jcut, nr
               fcut = exp(-5*(r(j)-r(jcut)))
               viou(lp,j) = -2*zion + fcut*(viou(lp,j)+2*zion)
  360       continue
            do 370 j = 2, nr
               v(j) = viou(lp,j)/r(j)
  370       continue
c
c  njtj  ***  plotting routines ***
c
            call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
            call potrv(v,r,nr-120,lo(i),zion)
c
c  njtj  ***  user should adjust for their needs  ***
c
         end if
  380 continue
      write(6,9080)
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
cag removed
c      write(3,9090)
c      write(3,9100) zion
c 9090 format(1x,'marker zio')
c 9100 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's') then
         do 410 i = ncp, norb, 2
            lp = lo(i) + 1
            zot = zo(i) + zo(i+1)
            if (zot .ne. zero) then
               do 390 j = 2, nr
                  viod(lp,j) = (viod(lp,j)*zo(i)+viou(lp,j)*zo(i+1))/zot
                  viou(lp,j) = viod(lp,j)
  390          continue
            else
               do 400 j = 2, nr
                  viod(lp,j) = viod(lp,j)/2 + viou(lp,j)/2
                  viou(lp,j) = viod(lp,j)
  400          continue
            end if
  410    continue
      end if
c
      do 420 i = 2, nr
         vid(i) = vod(i)
         viu(i) = vou(i)
  420 continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,nr,a,b,r,rab,norb-ncore,0,nops(ncp),
     &            lo(ncp),so(ncp),zo(ncp),znuc,cdd,cdu,cdc,viod,viou,
     &            vid,viu,ev(ncp),ek(ncp),ep(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,9110) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,9120) (ev(i),i=ncp,norb)
 9110 format(//' test of eigenvalues',//' rcut =',8(2x,a1,f7.2))
 9120 format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,9130)
 9130 format(///' l    vps(0)    vpsmin      at r',/)
      do 450 i = 1, lmax
         if (indd(i)+indu(i) .eq. 0) go to 450
         if (indd(i) .ne. 0) then
            vpsdm = zero
            do 430 j = 2, nr
               if (r(j) .lt. .00001D0) go to 430
               vps = viod(i,j)/r(j)
               if (vps .lt. vpsdm) then
                  vpsdm = vps
                  rmind = r(j)
               end if
  430       continue
            write(6,9140) il(i), viod(i,2)/r(2), vpsdm, rmind
         end if
         if (indu(i) .ne. 0) then
            vpsum = zero
            do 440 j = 2, nr
               if (r(j) .lt. .00001D0) go to 440
               vps = viou(i,j)/r(j)
               if (vps .lt. vpsum) then
                  vpsum = vps
                  rminu = r(j)
               end if
  440       continue
            write(6,9140) il(i), viou(i,2)/r(2), vpsum, rminu
         end if
 9140    format(1x,a1,3f10.3)
  450 continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,nops(ncp),lo(ncp),so(ncp),
     &            zo(ncp),etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = '   Hamann,'
      iray(4) = ' Schluter '
      iray(5) = 'and Chiang'
      iray(6) = ' potential'
c
c  Encode the title array.
c
      do 460 i = 1, 7
         ititle(i) = '          '
  460 continue
      do 470 i = 1, lmax
         if (indd(i) .eq. 0 .and. indu(i) .eq. 0) go to 470
         zelu = zero
         zeld = zero
         if (indd(i) .ne. 0) then
            noi = no(indd(i))
            zeld = zo(indd(i))
         end if
         if (indu(i) .ne. 0) then
            noi = no(indu(i))
            zelu = zo(indu(i))
         end if
         zelt = zeld + zelu
         if (ispp .ne. 's') then
            write(ititle(2*i-1),9150) noi, il(i), zelt
            write(ititle(2*i),9160) ispp, rc(i)
 9150       format(i1,a1,'(',f6.2,')')
 9160       format(a1,' rc=',f5.2)
         else
            write(ititle(2*i-1),9170) noi, il(i), zeld
            write(ititle(2*i),9180) zelu, ispp, rc(i)
 9170       format(i1,a1,'  (',f4.2,',')
 9180       format(f4.2,')',a1,f4.2)
         end if
  470 continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
         if (indu(1) .eq. 0) go to 490
         indd(1) = indu(1)
         indu(1) = 0
         do 480 j = 2, nr
            viod(1,j) = viou(1,j)
            viou(1,j) = zero
  480    continue
  490    continue
         do 510 i = 2, lmax
            if (indd(i) .eq. 0 .or. indu(i) .eq. 0) go to 510
            do 500 j = 2, nr
               viodj = viod(i,j)
               viouj = viou(i,j)
               viod(i,j) = ((i-1)*viodj+i*viouj)/(2*i-1)
               viou(i,j) = 2*(viouj-viodj)/(2*i-1)
  500       continue
  510    continue
      end if
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 520 i = 1, lmax
         if (indd(i) .ne. 0) npotd = npotd + 1
         if (indu(i) .ne. 0) npotu = npotu + 1
  520 continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
         if (ifull .eq. 0) then
            nicore = 'pcec'
         else
            nicore = 'fcec'
         end if
      else if (ifcore .eq. 2) then
         if (ifull .eq. 0) then
            nicore = 'pche'
         else
            nicore = 'fche'
         end if
      else
         nicore = 'nc  '
      end if
      if (ispp .eq. 's') then
         irel = 'isp'
      else if (ispp .eq. 'r') then
         irel = 'rel'
      else
         irel = 'nrl'
      end if
      rewind 1
      write(1) nameat, icorr, irel, nicore, (iray(i),i=1,6),
     &  (ititle(i),i=1,7), npotd, npotu, nr - 1, a, b, zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 530 i = 1, lmax
         if (indd(i) .eq. 0) go to 530
         write(1) i - 1, (viod(i,j),j=2,nr)
  530 continue
      do 540 i = 1, lmax
         if (indu(i) .eq. 0) go to 540
         write(1) i - 1, (viou(i,j),j=2,nr)
  540 continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .ne. 1) then
         write(1) (zero,i=2,nr)
      else
         write(1) (cdc(i),i=2,nr)
      end if
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
c
      end
