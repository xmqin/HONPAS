c
      subroutine Input
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
c  subroutine to read input parameters
c
C     .. Parameters ..
      double precision one, zero, pfive
      parameter (one=1.D0,zero=0.D0,pfive=0.5D0)
c
C     .. Local Scalars ..
      double precision rmax, sc, si, zcore, zd, zu,  zval
      integer i, j, li, ni
      character type*2, flavor*3, compat_str*20
C     ..
C     .. Local Arrays ..
      integer lc(19), nc(19), nomin(0:4)
C     ..
C     .. External Functions ..
      double precision nucl_z
      logical leqi
      external nucl_z, leqi
C     ..
C     .. External Subroutines ..
      external ext, cal_date
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs, exp
C     ..
c
c  data for orbitals:
c 
c              1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p, 
c                     7s, 5f, 6d, 7p
c
      data nc /1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6,
     $         7, 5, 6, 7/
      data lc /0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1,
     $         0, 3, 2, 1/
c
c     This statement now goes in a block data subprogram  
cccc      data il /'s', 'p', 'd', 'f', 'g'/
C     ..
c
      do 10 i = 0, 4
         nomin(i) = 10
   10 continue
      do 20 i = 1, norbmx
         no(i) = 0
         lo(i) = 0
         so(i) = zero
         zo(i) = zero
   20 continue
c
c
c  read the type of calculation and title card
c   job =
c   ae = 0 all electron calculation
c   pg = 1 pseudopotential generation w/o core correction
c   pe = 2 pseudopotential generation w/  core correction exchange
c   ph = 3 pseudopotential generation w/  core correction hartree
c   pt = 4 pseudopotential test
c   pm = 5 pseudopotential test + valence charge modify
c
      job = -1
      read(5,9000,err=999,end=999) type, title
 9000 format(3x,a2,5a10)
c
c  if type = ' ' , no more data, program ends
c
      if (type .eq. 'ae') then
         job = 0
      else if (type .eq. 'pg') then
         job = 1
      else if (type .eq. 'pe') then
         job = 2
      else if (type .eq. 'ph') then
         job = 3
      else if (type .eq. 'pt') then
         job = 4
      else if (type .eq. 'pm') then
         job = 5
      else
         job = -1
c
         return
c
      end if
c
      ifcore = job - 1
c
c  njtj  ***  major modification  start  ***
c  There are seven ways to generate the pseudopotential :
c    flavor = van Vanderbilt
c    flavor = tam Troullier and Martins
c    flavor = ker (yes) Kerker
c    flavor = hsc (no)  Hamann Schluter and Chiang
c    flavor = min (oth) datafile made for minimization
c    flavor = bhs Bachelet, Hamann and Schluter
c    flavor = tm2 Improved Troullier and Martins
c                      
c
c     Flavor should not be necessary for the tests...
c
      if ((job .gt. 0) .and. (job .lt. 4)) then
         read(5,9010) flavor, logder_radius
 9010    format(8x,a3,f9.3)
         if (leqi(flavor,'tm2')) then
            scheme = 6
         else if (leqi(flavor,'bhs')) then
            scheme = 5
         else if (leqi(flavor,'oth') .or. leqi(flavor,'min')) then
            scheme = 4
         else if (leqi(flavor,'tbk') .or. leqi(flavor,'tam')) then
            scheme = 2
         else if (leqi(flavor,'yes') .or. leqi(flavor,'ker')) then
            scheme = 1
         else if (leqi(flavor,'no ') .or. leqi(flavor,' no')
     &            .or. leqi(flavor,'hsc') ) then
            scheme = 0
         else
            write(6,9020) flavor
            call ext(150)
         end if
      end if
 9020 format(//'error in input - flavor =',a3,' unknown')
c  njtj  ***  major modification end  ***
c
c   read element name, correlation type, polarization flag...
c   ispp = ' ' - nonspin calculation
c   ispp = s  - spin polarized calculation
c   ispp = r  - relativistic calculation
c
c   ... and a compatibility string (obsolete -- to be removed)
c
      read(5,9030) nameat, icorr, ispp, compat_str
 9030 format(3x,a2,3x,a2,a1,1x,a20)
c
      call compat_params(compat_str)
c
      if (.not.( leqi(ispp,'s') .or. leqi(ispp,'r') ) ) ispp = ' '
c
c     Not all the correlation schemes can be used for a
c     spin-polarized calculation.
c
      if ( leqi(ispp,'s') .and. 
     &     ( leqi(icorr,'xa') .or. leqi(icorr,'wi') .or.
     &       leqi(icorr,'hl') .or. leqi(icorr,'gr'))
     &   )   ispp = ' '
c
      normal = leqi(ispp,' ')
      relativistic = leqi(ispp,'r')
      polarized = leqi(ispp,'s')
c
c  njtj   ***  major modification start  ***
c   Floating point comparison error modification.
c   Read the atomic number (nuclear charge),
c   shell charge and radius (added to the nuclear potential),
c   and radial grid parameters.
c
      read(5,9040) znuc, zsh, rsh, rmax, aa, bb
 9040 format(6f10.3)
      if (abs(znuc) .le. 0.00001D0) znuc = nucl_z(nameat)
      if (job .lt. 4) then
c
c   set up grid
c
         if (abs(rmax) .lt. 0.00001D0) rmax = rmax_def
         if (abs(aa) .lt. 0.00001D0) aa = aa_def
         if (abs(bb) .lt. 0.00001D0) bb = bb_def
         a = exp(-aa)/znuc
         b = 1/bb
         do 30 i = 1, nrmax
            r(i) = a*(exp(b*(i-1))-1)
            rab(i) = (r(i)+a)*b
            if (r(i) .gt. rmax) go to 40
   30    continue
c
         write(6,9050)
 9050    format(/' error in input - arraylimits',
     &      ' for radial array exceeded',/)
c
         call ext(100)
c
   40    continue
c
         nr = i - 1
c
      end if
c  njtj  ***  major modification end  ***
c
c   read the number of core and valence orbitals
c
c
      read(5,9060) ncore, nval
 9060 format(2i5,2f10.3)
      if (ncore .gt. 19) then
         write(6,9070)
         call ext(101)
      end if
 9070 format(//'error in input - max number of core orbitals','is 19')
c
c   compute occupation numbers and orbital energies for the core
c
      zcore = zero
      if (ncore .eq. 0) go to 70
      sc = zero
      if (ispp .ne. ' ') sc = -pfive
      norb = 0
      do 60 i = 1, ncore
         do 50 j = 1, 2
            if (ispp .eq. ' ' .and. j .eq. 2) go to 50
            norb = norb + 1
            no(norb) = nc(i)
            lo(norb) = lc(i)
            so(norb) = sc
            zo(norb) = 2*lo(norb) + 1
            if (ispp .eq. ' ') zo(norb) = 2*zo(norb)
            if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc) + 1
            zcore = zcore + zo(norb)
c
c           If lo==0 and sc=-0.5, zo=0. This means that
c           s orbitals are "up" in a relativistic calculation.
c
            if (abs(zo(norb)) .lt. 0.1D0) norb = norb - 1
            if (ispp .ne. ' ') sc = -sc
   50    continue
   60 continue
      ncore = norb
cc
cc     Example of manual fiddling
cc     Ionize an electron from 1s...
cc
cc      zo(1) = 1.d0
cc      zcore = zcore - 1.d0
cc
c   for the valence orbitals
c
   70 continue
      if (job .ge. 4) ncore = 0
      norb = ncore
      zval = zero
      if (nval .eq. 0) go to 130
      do 90 i = 1, nval
         read(5,9060) ni, li, zd, zu
         si = zero
         if (ispp .ne. ' ') si = pfive
         do 80 j = 1, 2
            if (ispp .eq. ' ' .and. j .eq. 2) go to 80
            norb = norb + 1
            if (ispp .ne. ' ') si = -si
            no(norb) = ni
            lo(norb) = li
            so(norb) = si
            zo(norb) = zd + zu
            if (ispp .eq. 's') then
               if (si .lt. 0.1D0) then
                  zo(norb) = zd
               else
                  zo(norb) = zu
               end if
            else if (ispp .eq. 'r') then
               zo(norb) = zo(norb)*(2*(li+si)+1)/(4*li+2)
            end if
            zval = zval + zo(norb)
            if (ispp .eq. 'r' .and. li+si .lt. zero) norb = norb - 1
            if (norb .eq. 0) go to 80
            if (nomin(lo(norb)) .gt. no(norb)) 
     &                               nomin(lo(norb)) = no(norb)
   80    continue
   90 continue
c
c   Compute 'down' flag and
c   abort if two orbitals are equal
c
      nval = norb - ncore
      do 110 i = 1, norb
         down(i) = (so(i) .lt. 0.1D0)
         do 100 j = 1, norb
            if (i .le. j) go to 100
            if (no(i) .ne. no(j)) go to 100
            if (lo(i) .ne. lo(j)) go to 100
            if (abs(so(i)-so(j)) .gt. 0.001D0) go to 100
            write(6,9080) i
            call ext(110+i)
  100    continue
  110 continue
 9080 format(//'error in input - orbital ',i2,'is already occupied')
c
c   reduce n quantum number if pseudoatom
c
      if (job .ge. 4) then
         do 120 i = 1, nval
            no(i) = no(i) - nomin(lo(i)) + lo(i) + 1
  120    continue
      end if
  130 continue
      zion = znuc - zcore - zval
      zel = zval
      if (job .lt. 4) then
         zel = zel + zcore
      else
         znuc = znuc - zcore
      end if
c
c   find jobname and date and printout.
c
      ray(1) = atom_id
      call cal_date(ray(2))
cag
      ncp = ncore + 1
cag
 999  continue
      return
c
      end
c
      block data orb_init
      implicit none
      include 'orbital.h'
      data il /'s', 'p', 'd', 'f', 'g'/
      end

