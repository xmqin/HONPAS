c
      subroutine wrapup(pot_id)
c
      implicit none
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'ion.h'
      include 'charge.h'
      include 'elecpot.h'
      include 'energy.h'
      include 'compat.h'
      include 'version.h'
      include 'pseudowave.h'
c
c     ecuts now set in compat_params...
c
      double precision zero, tpfive, one
      parameter (zero=0.d0,tpfive=2.5d0, one=1.d0)
c
      integer i, icore, j, jcut, lp, noi, npotd, npotu, ifull
      integer nops(norbmx), position, iunit
      character ray(6)*10, title*70, pot_id*40, id*1
c
      double precision zval, zratio, zion, ac, bc, cdcp, tanb, rbold,
     &                 rbnew, pi, ecut, vp2z, fcut, zot, vpsdm,
     &                 vps, rmind, vpsum, rminu, zelu, zeld, zelt,
     &                 viodj, viouj, cc
      double precision rcut(10), v(nrmax)
c
      double precision absval, minabs, maxabs, norm1, norm2
      double precision fourier_area(5), qc(5), dummy_real, dummy_qc
      double precision fourier_eps
      parameter (fourier_eps = 1.0d-2)

      integer n_channels, lun

      double precision cutoff_function, force_underflow
      external cutoff_function, force_underflow

      logical new_scheme
c
      external logder
c
      pi = 4*atan(one)
c
c     Do not use relativity for what follows
c
      if (polarized) then
         id = 's'
      else
         id = ' '
      endif
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
c
      do 5 i = 1, norb
         nops(i) = 0
    5 continue
c
      do 10 i = ncp, norb
         nops(i) = lo(i) + 1
         zval = zval + zo(i)
   10 continue
      zion = zval + znuc - zel
      if (zval .ne. zero) zratio = zion/zval
c
      do 20 i = 1, nr
         cdd(i) = vod(i)
         cdu(i) = vou(i)
   20 continue
c
c=====================================================================
c  If a core correction is indicated construct pseudo core charge
c  if cfac < 0 or the valence charge is zero the full core is used
c
      ifcore = job - 1
      if (ifcore .ne. 0) then
         ac = zero
         bc = zero
         cc = zero
         icore = 1
         if (cfac .le. zero .or. zratio .eq. zero) then
            write(6,9000) r(icore), ac, bc, cc
            write(6,'(a)') '(Full core used)'
            call coreq
         else
            if (rcfac .le. zero) then
               do 30 i = nr, 2, -1
                  if (cdc(i) .gt. cfac*zratio*
     &                (cdd(i)+cdu(i))) go to 50
   30          continue
            else
               do 40 i = nr, 2, -1
                  if (r(i) .le. rcfac) go to 50
   40          continue
            end if
   50       continue
            icore = i
C
C--------------------------------------------------------------------
C Choice of methods to compute the core correction:
C
C 1. Traditional 'Froyen-Louie-Cohen' with cdc(r) = Arsin(Br)
C    and value and first-derivative matching. It is the default
c    for LDA calculations. See compat_params.f
C    and input.f for info on how to force its use from the input file.
C
C 2. New by Jose Luis Martins' group, using a Kerker-like exp( ) function
C    and matching also the second derivative. This is the default with
C    the 'mons' compatibility mode for GGA calculations.
C

            new_scheme = .false.
            if (is_gga) then
               write(6,'(a)') 'Note: GGA calculation ==> New CC scheme'
               new_scheme = .true.
            endif

            if ( (use_old_cc) .and. (is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using old-style core corrections',
     $            ' despite this being a GGA calculation.',
     $              'I hope you know what you are doing...'
               new_scheme = .false.
            endif

            if ( (use_new_cc) .and. (.not. is_gga)) then

               write(6,'(/,2a,/a)')
     $              'WARNING: Using new core corrections',
     $            ' despite this being an LDA calculation.',
     $         'Results will not be compatible with older versions.'
               new_scheme = .true.
            endif

            if (.not. new_scheme) then

c           Fit to  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c
c           Find derivative at core radius. Use a five-point formula
c           instead of the old two-point method:
c              cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))
c           Since the r abscissae are not equally spaced, we perform
c           the derivative using the chain rule:
c
c           r(i) = a [ exp(b*(i-1)) - 1 ]
c           r(x) = a [ exp(b*x) - 1 ]
c
c           f'(r) = f'(x) / r'(x)
c
c           r'(x) = a*b*exp(b*x) = (r(x)+a)*b
c           r'(i) = a*b*exp[b*(i-1)] = rab(i)
c
c           To compute f'(x), use eq. 25.3.6 
c           (p. 883) in Abramowitz-Stegun, with h=1
c
c           f'(0) = 1/12 f(-2) - 2/3 f(-1) + 2/3 f(1) - 1/12 f(2)
c
cold            cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))

            cdcp =   cdc(icore-2) - 8*cdc(icore-1) +
     $             8*cdc(icore+1) -   cdc(icore+2)
            cdcp = cdcp / 12.d0 / rab(icore)
c
c           Now fit ac and bc using the function and derivative
c           information. 
c
c           g(r) = Arsin(Br) ==>  g / (rg'-g) = tanBr/(Br)
c
c           Use an iterative method to find Br and from that ac and bc.
c           Start near Br=2.5 so that we are in an invertible region 
c           of the tanx/x function.
c 
c
            tanb = cdc(icore)/(r(icore)*cdcp-cdc(icore))
            rbold = tpfive
            do 70 i = 1, 50
               rbnew = pi + atan(tanb*rbold)
               if (abs(rbnew-rbold) .lt. .00001D0) then
                  bc = rbnew/r(icore)
                  ac = cdc(icore)/(r(icore)*sin(rbnew))
                  do 60 j = 1, icore
                     cdc(j) = ac*r(j)*sin(bc*r(j))
   60             continue
                  write(6,9000) r(icore), ac, bc, cc
c
                  call coreq
                  go to 80
c
               else
                  rbold = rbnew
               end if
   70       continue
            write(6,9010)
            call ext(830)

            else

c                 Use subroutine provided by JLM to fit
c                 cdc(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
c
                  CALL PCC_EXP(NR,ICORE,AC,BC,CC,R,CDC)
                  write(6,9000) r(icore), ac, bc, cc

c
           endif

           call coreq

         end if
      end if
C---------------------------------------------------------------------
 9000 format(//' Core correction used',/' Pseudo core inside r =',f6.3,
     &      /' ac =',f6.3,' bc =',f6.3,' cc =',f6.3,/)
 9010 format(//' Error in pseudo - nonconvergence in finding ',
     &      /'pseudo-core values')
c
c  End the pseudo core charge.
c======================================================================
c
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentials should be unscreened with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
   80 continue
c
      call Velect(0,1,id,zval)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius.
c
c  Note that the cutting off of the pseudopotentials (making
c  them approach -2*Zion/r faster) is not strictly necessary.
c  It might even be argued that it should be left to "client"
c  programs to decide what to do.

cag
c
c     On the issue of plotting:
c
c     For non-relativistic, non-spin-polarized calculations, all
c     the orbitals are considered as "down".
c     For relativistic calculations, the "s" orbitals are considered
c     as "up". The actual things plotted on the files (which only
c     record l, not the up/down character) depend on the order of
c     enumeration of the orbitals. Since these are always "down/up",
c     the potentials plotted are always the "up" ones (except of
c     course for scalar calculations and for "s" states in relativistic
c     calculations, for which the distinction is irrelevant.

      write(6,9020)
 9020 format(/)
      ecut = ecuts
      do 150 i = ncp, norb
         lp = lo(i) + 1
         if (down(i)) then
            do 90 j = 2, nr
               v(j) = viod(lp,j)/r(j) + vid(j)
               viod(lp,j) = viod(lp,j) + (vid(j)-vod(j))*r(j)
               vp2z = viod(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
   90       continue
cag
c           Plot screened ionic potential
c
            call potrvs(v,r,nr-120,lo(i))
cag
c           Default cutoff function: f(r)=exp(-5*(r-r_cut)). It damps
c           down the residual of rV+2*Zion.
c           Should be made smoother... Vps ends up with a kink at rcut.
c           Maybe use one of the Vanderbilt generalized gaussians.
cag
            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 100 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viod(lp,j) = -2*zion + fcut*(viod(lp,j)+2*zion)
  100       continue
            do 110 j = 2, nr
               v(j) = viod(lp,j)/r(j)
  110       continue
c
c           Dwon potentials are  always generated
c
            call potran(lo(i)+1,v,r,nr,zion,fourier_area(lo(i)+1),
     $                  fourier_eps,qc(lo(i)+1))
            call potrv(v,r,nr-120,lo(i),zion)
c
         else
            do 120 j = 2, nr
               v(j) = viou(lp,j)/r(j) + viu(j)
               viou(lp,j) = viou(lp,j) + (viu(j)-vou(j))*r(j)
               vp2z = viou(lp,j) + 2*zion
               if (abs(vp2z) .gt. ecut) jcut = j
  120       continue
cag
c           Plot screened ionic potential
c
            call potrvs(v,r,nr-120,lo(i))
cag
            rcut(i-ncore) = r(jcut)
            if (rcut(i-ncore) .lt. rc(lp)) then
               write(6,'(a,2f8.4)') 'Vps rcut point moved out to rc: ',
     $              rcut(i-ncore), rc(lp)
               rcut(i-ncore) = rc(lp)
            endif
            do 130 j = jcut, nr
cag               fcut = exp(-5*(r(j)-r(jcut)))
               fcut = cutoff_function(r(j)-r(jcut))
               viou(lp,j) = -2*zion + fcut*(viou(lp,j)+2*zion)
  130       continue
            do 140 j = 2, nr
               v(j) = viou(lp,j)/r(j)
  140       continue
c
c           Up potentials are not always generated
c
            call potran(lo(i)+1,v,r,nr,zion,dummy_real,
     $           fourier_eps,dummy_qc)
            call potrv(v,r,nr-120,lo(i),zion)
c
         end if
c
  150 continue
c
c     Write out the Fourier area for each pseudo channel
c
      call get_unit(lun)
      open(unit=lun,file="FOURIER_AREA",form="formatted",
     $     status="unknown")
      rewind(lun)
c
c     Compute also the minimum, maximum, mean, and root-mean-square.
c
      n_channels = 0
      maxabs = -1.0d0
      minabs = 1.0d10
      norm1 = 0.0d0
      norm2 = 0.0d0
      do j=lo(ncp) + 1, lo(norb) + 1
         absval = fourier_area(j)
         if (absval .gt. maxabs) maxabs = absval
         if (absval .lt. minabs) minabs = absval
         norm1 = norm1 + absval
         norm2 = norm2 + absval*absval
         n_channels =  n_channels + 1
      enddo
      norm1 = norm1 / n_channels
      norm2 = sqrt( norm2 / n_channels)
      write(lun,"(i4)") n_channels
      write(lun,"(5f10.5)") (fourier_area(j),j=lo(ncp)+1,lo(norb)+1)
      write(lun,"(4f10.5)") minabs, maxabs, norm1, norm2
      close(lun)
c
c     Write out the Fourier threshold for each channel
c
      call get_unit(lun)
      open(unit=lun,file="FOURIER_QMAX",form="formatted",
     $     status="unknown")
      rewind(lun)
c
c     Compute also the minimum, maximum, mean, and root-mean-square.
c
      n_channels = 0
      maxabs = -1.0d0
      minabs = 1.0d10
      norm1 = 0.0d0
      norm2 = 0.0d0
      do j=lo(ncp) + 1, lo(norb) + 1
         absval = qc(j)
         if (absval .gt. maxabs) maxabs = absval
         if (absval .lt. minabs) minabs = absval
         norm1 = norm1 + absval
         norm2 = norm2 + absval*absval
         n_channels =  n_channels + 1
      enddo
      norm1 = norm1 / n_channels
      norm2 = sqrt( norm2 / n_channels)
      write(lun,"(i4)") n_channels
      write(lun,"(5f14.5)") (qc(j),j=lo(ncp)+1,lo(norb)+1)
      write(lun,"(4f14.5)") minabs, maxabs, norm1, norm2
      close(lun)
c
      write(6,9020)
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (polarized) then
         do 180 i = ncp, norb, 2
            lp = lo(i) + 1
            zot = zo(i) + zo(i+1)
            if (zot .ne. zero) then
               do 160 j = 2, nr
                  viod(lp,j) = (viod(lp,j)*zo(i)+viou(lp,j)*zo(i+1))/zot
                  viou(lp,j) = viod(lp,j)
  160          continue
            else
               do 170 j = 2, nr
                  viod(lp,j) = viod(lp,j)/2 + viou(lp,j)/2
                  viou(lp,j) = viod(lp,j)
  170          continue
            end if
  180    continue
      end if
c
      do 190 i = 2, nr
         vid(i) = vod(i)
         viu(i) = vou(i)
  190 continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,id,ncp,norb,0,nops)
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,9030) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,9040) (ev(i),i=ncp,norb)
 9030 format(//' test of eigenvalues',//' rcut =',8(2x,a1,f7.2))
 9040 format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,9050)
 9050 format(///' l    vps(0)    vpsmin      at r',/)
c
      do 220 i = 1, lmax
         if (indd(i)+indu(i) .eq. 0) go to 220
         if (indd(i) .ne. 0) then
            vpsdm = zero
            do 200 j = 2, nr
               if (r(j) .lt. .00001D0) go to 200
               vps = viod(i,j)/r(j)
               if (vps .lt. vpsdm) then
                  vpsdm = vps
                  rmind = r(j)
               end if
  200       continue
            write(6,9060) il(i), viod(i,2)/r(2), vpsdm, rmind
         end if
         if (indu(i) .ne. 0) then
            vpsum = zero
            do 210 j = 2, nr
               if (r(j) .lt. .00001D0) go to 210
               vps = viou(i,j)/r(j)
               if (vps .lt. vpsum) then
                  vpsum = vps
                  rminu = r(j)
               end if
  210       continue
            write(6,9060) il(i), viou(i,2)/r(2), vpsum, rminu
         end if
 9060    format(1x,a1,3f10.3)
  220 continue
c
c   Print out the energies from etotal. (Valence only...)
c
      call etotal(ncp,norb)
c
c   Compute the logarithmic derivative as a function of energy 
c
      if (logder_radius .gt. 0.d0) call logder(ncp,norb,'PS')
c
c  Find the jobname and date.
c
      ray(1) = atom_id
      call cal_date(ray(2))
c  
      read(pot_id,'(4a10)') (ray(i),i=3,6)
c
c  Encode the title array.
c
      title = ' '
      position = 1
      do 240 i = 1, lmax
         if (indd(i) .eq. 0 .and. indu(i) .eq. 0) go to 240
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
         if ( .not. polarized) then
            write(title(position:),9070) noi, il(i), zelt, ispp, rc(i)
 9070       format(i1,a1,f5.2,a1,' r=',f5.2,'/')
            position = position + 17
         else
            write(title(position:),9090)
     $                             noi, il(i), zeld, zelu, ispp, rc(i)
 9090       format(i1,a1,f4.2,',',f4.2,a1,f4.2,'/')
            position = position + 17
         end if
  240 continue
c
c  Construct relativistic sum and difference potentials.
c
      if (relativistic) then
c
c   ***  The s potential is from now on considered as "down", even
c        though s=0.5 in the relativistic case.
c
         if (indu(1) .eq. 0) go to 260
         indd(1) = indu(1)
         indu(1) = 0
         do 250 j = 2, nr
            viod(1,j) = viou(1,j)
            viou(1,j) = zero
  250    continue
  260    continue
         do 280 i = 2, lmax
            if (indd(i) .eq. 0 .or. indu(i) .eq. 0) go to 280
            do 270 j = 2, nr
               viodj = viod(i,j)
               viouj = viou(i,j)
               viod(i,j) = ((i-1)*viodj+i*viouj)/(2*i-1)
               viou(i,j) = 2*(viouj-viodj)/(2*i-1)
  270       continue
  280    continue
      end if
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 290 i = 1, lmax
         if (indd(i) .ne. 0) npotd = npotd + 1
         if (indu(i) .ne. 0) npotu = npotu + 1
  290 continue
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
      if (polarized) then
         irel = 'isp'
      else if (relativistic) then
         irel = 'rel'
      else
         irel = 'nrl'
      end if

      open(unit=1,file='VPSOUT',status='unknown',form='unformatted')
      rewind 1
      open(unit=2,file='VPSFMT',status='unknown',form='formatted')
      rewind 2
      open(unit=3,file='PSWFFMT',status='unknown',form='formatted')
      rewind 3

      write(1) nameat, icorr, irel, nicore, (ray(i),i=1,6), title,
     &         npotd, npotu, nr - 1, a, b, zion
      write(1) (r(i),i=2,nr)
c
      do 700 iunit=2,3
         write(iunit,8005) nameat, icorr, irel, nicore
         write(iunit,8010) (ray(j),j=1,6), title
         write(iunit,8015) npotd, npotu, nr-1, a, b, zion
         write(iunit,8040) 'Radial grid follows' 
         write(iunit,8030) (r(j),j=2,nr)
 700  continue
c
 8000 format(1x,i2)
 8005 format(1x,a2,1x,a2,1x,a3,1x,a4)
 8010 format(1x,6a10,/,1x,a70)
 8015 format(1x,2i3,i5,3g20.12)
 8030 format(4(g20.12))
 8040 format(1x,a)
c
c  Write the potentials to file (unit=1, and 2).
c
      do 300 i = 1, lmax
         if (indd(i) .eq. 0) go to 300
         write(1) i - 1, (viod(i,j),j=2,nr)
         write(2,8040) 'Down Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (force_underflow(viod(i,j)),j=2,nr)
  300 continue
      do 310 i = 1, lmax
         if (indu(i) .eq. 0) go to 310
         write(1) i - 1, (viou(i,j),j=2,nr)
         write(2,8040) 'Up Pseudopotential follows (l on next line)'
         write(2,8000) i-1
         write(2,8030) (force_underflow(viou(i,j)),j=2,nr)
 310  continue
c
c  Write the charge densities      
c  Note that this charge density is the "pseudo" one.
c
      write(2,8040) 'Core charge follows'

      if (ifcore .ne. 1) then
         write(1) (zero,i=2,nr)
         write(2,8030) (zero,i=2,nr)
      else
         write(1) (cdc(i),i=2,nr)
         write(2,8030) (force_underflow(cdc(i)),i=2,nr)
      end if
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
      write(2,8040) 'Valence charge follows'
      write(2,8030) (force_underflow(zratio*(cdd(i)+cdu(i))),i=2,nr)

      close(1)
      close(2)
c
c  Write the pseudo-wavefunctions (only in formatted form, in
c  auxiliary file) 'u_n,l (r) = 1/r R_n,l (r)'
c
 8045 format(1x,i2,i2)
      do 600 i = 1, lmax
         if (indd(i) .eq. 0) go to 600
         write(3,8040)
     $      'Down Pswavefunction (R/r) follows (l,n on next line)'
         write(3,8045) i-1, no(indd(i))
         write(3,8030) (force_underflow(pswfnrd(i,j)),j=2,nr)
  600 continue
      do 620 i = 1, lmax
         if (indu(i) .eq. 0) go to 620
         write(3,8040)
     $      'Up Pswavefunction (R/r) follows (l,n on next line)'
         write(3,8045) i-1, no(indu(i))
         write(3,8030) (force_underflow(pswfnru(i,j)),j=2,nr)
  620 continue
c
      close(3)
c
      open(unit=3,file='PSCHARGE',form='formatted',status='unknown')
c
c     NOTE: We no longer put "zratio" here!!!
c     (We still do in the ps file for compatibility with PW and
c      SIESTA) 
c     (Only affects plots for ionic configurations)
c     BEAR THIS IN MIND IF YOU ARE USING THE HEURISTIC CORE CORRECTION
c     CRITERION: If you specify a given "pc_weight" in the input file,
c     do not be surprised if the plot does not show rcore
c     in the place you expect it to be.
c
      if (ifcore .ne. 1) then
         do 400 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), zero
 400     continue
      else
         do 410 j = 2, nr
            write(3,9900) r(j), cdu(j), cdd(j), cdc(j)
 410     continue
      endif
c
 9900 format(1x,f15.10,3x,3f15.8)
c
      close(unit=3)
c
c     Write the pseudopotential in XML format
c
      call pseudoXML( ray, npotd, npotu, zion, zratio )

c
      return
c
      end

      double precision function cutoff_function(r)
      implicit none
c
c     Generalized cutoff function
c
      double precision r
c
c     Standard cutoff
c
      cutoff_function = exp(-5.d0*r)

      end










