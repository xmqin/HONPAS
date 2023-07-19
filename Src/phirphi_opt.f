! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine phirphi_opt( nua, na, nuo, no, scell, xa,
     .                        maxnh, lasto, 
     .                        lastkb, iphorb, iphKB, isa, 
     .                        numh, listhptr, listh, indxuo, 
     .                        indxua, iaorb, dk, matrix, Sp)
C *********************************************************************
C If matrix='R'
C Finds the matrix elements of the dk*r between basis
C orbitals, where dk is a given vector and r is the position operator,
C or the momentum operator plus a correction due to the non-local 
C potential in the case of an infinite solid.
C
C  S_a,b(R_a,R_b) =  <phi_a(r-R_a-t_a)|dk*r|phi_b(r-R_b-t_b)>
C
C  
C  Being R_a and R_b lattice vectors, and t_a and t_b the coordiantes
C  of atoms a and b in the unit cell
C
C If matrix='P'
C Finds the matrix elements of the -i*dk*p (i.e. minus dk dot gradient )
C between basis orbitals, phi_a(r-R_a) and phi_b(r-R_b),
C where dk is a given vector and p is the momentum operator,
C a correction must be included due to the use of non-local pseudoptentials
C VNL=|f(r-R_c)><f(r'-R_c)|
C
C  S_a,b(R_a,R_b) = -i*2.0d0*<phi_a(r-R_a-t_a)|dk*p|phi_b(r-R_b-t_b)>+
C    <phi_a(r-R_a-t_a)|f(r-R_c)><f(r'-R_c)|dk*(r'-R_c)|phi_b(r'-R_b-t_b)>
C   -<phi_a(r-R_a-t_a)|dk*(r-R_c)|f(r-R_c)><f(r'-R_c)|phi_b(r'-R_b-t_b)>
C
C  The factor of two in front of dk*p is related to the use of Ry and
C  Bohr instead of Ha and Ry. The factor will be canceled out
C  afterwards when we divide by (E_i-E-j), the energy denominator.
C  Being R_a and R_b lattice vectors, and t_a and t_b the coordiantes
C  of atoms a and b in the unit cell
C

C Energies in Ry. Lengths in Bohr.
C Written by DSP August 1999 ( Based in routine overfsm and nlefsm)
C Restyled for f90 version by JDG. June 2004
C **************************** INPUT **********************************
C integer nua              : Number of atoms in unit cell
C integer na               : Number of atoms in supercell
C integer nuo              : Number of basis orbitals in unit cell
C integer no               : Number of basis orbitals in super cell
C real*8  scell(3,3)       : Supercell vectors SCELL(IXYZ,IVECT)
C real*8  xa(3,na)         : Atomic positions in cartesian coordinates
C integer maxnh            : First dimension of listh
C integer lasto(0:na)      : Last orbital index of each atom
C integer lastkb(0:na)     : Last KB porjector index of each atom
C integer iphorb(no)       : Orbital index of each orbital in its atom
C integer iphKB(nokb)      : KB proj. index of each KB proj. in its atom
C integer isa(na)          : Species index of each atom
C integer numh(nuo)        : Number of nonzero elements of each row
C                            of the overlap matrix
C integer listh(maxnh)     : Column indexes of the nonzero elements  
C                            of each row of the overlap matrix
C integer indxuo(no)       : Index of equivalent orbital in unit cell
C integer indxua(na)   : Index of equivalent atom in unit cell
C integer iaorb(no)        : Atom to which orbitals belong
C real*8  dk(3)            : Vector in k-space
C character*1 matrix       : 'R' for molecules or atoms, the matrix
C                             elements of the position operator are 
C                             calculated
C                            'P' for solids, the matrix elements of 
C                             the momentum operator are calculated.
C **************************** OUTPUT *********************************
C real*8  Sp(maxnh)        : Sparse overlap matrix
C *********************************************************************

      use precision,    only : dp
      use parallel,     only : Node, Nodes
      use atmparams,    only : lmx2, nzetmx, nsemx
      use atmfuncs,     only : epskb, lofio, mofio, rcut, rphiatm
      use atm_types,    only : nspecies
      use parallelsubs, only : GlobalToLocalOrb, LocalToGlobalOrb
      use alloc,        only : re_alloc
      use sys,          only : die
      use neighbour,    only : jna=>jan, xij, r2ij
      use neighbour,    only : mneighb

      implicit none

C Passed variables
      integer  maxnh, na, no, nua, nuo

      integer
     .  iphorb(no), isa(na), lasto(0:na), lastkb(0:na),
     .  listh(maxnh), listhptr(nuo), numh(nuo), iaorb(no), indxuo(no),
     .  iphKB(*), indxua(na)
      
      real(dp)
     .  scell(3,3), dk(3), Sp(maxnh), xa(3,na)

      character 
     .  matrix*1

C Internal variables
C maxkba = maximum number of KB projectors of one atom
C maxno  = maximum number of basis orbitals overlapping a KB projector
        
      integer
     .  ia, iio, ind, io, ioa, is, ix, 
     .  j, ja, jn, jo, joa, js, nnia, norb, ka

      real(dp)
     .  grSij(3), rij, Sij, xinv(3), sum

      integer
     .  ikb, ina, ino, jno, ko, koa, ks,
     .  nkb, nna, nno, ilm1, ilm2, npoints, 
     .  ir

      real(dp)
     .  epsk, grSki(3), 
     .  rki, rmax, rmaxkb, rmaxo,
     .  Sik, Sjk, Sikr, Sjkr, 
     .  dintg2(3), dint1, dint2, dintgmod2, 
     .  dintg1(3), dintgmod1,
     .  phi1, phi2, dphi1dr, dphi2dr, Sir0, r
      
      integer,  pointer, save :: iano(:)
      integer,  pointer, save :: iono(:)
      integer,           save :: maxkba = 25
      integer,           save :: maxno = 1000
      logical,  pointer, save :: calculated(:,:,:)
      logical,           save :: frstme = .true.
      logical,  pointer, save :: listed(:)
      logical,  pointer, save :: needed(:)
      logical                 :: within
      real(dp),          save :: dx = 0.01d0
      real(dp), pointer, save :: Pij(:,:,:)
      real(dp), pointer, save :: Si(:)
      real(dp), pointer, save :: Ski(:,:,:)
      real(dp),          save :: tiny = 1.0d-9
      real(dp), pointer, save :: Vi(:)

      external  matel
      
C Start timer
      call timer('phirphiopt',1)

C Check input matrix
      if(matrix.ne.'P'.and.matrix.ne.'R')
     $   call die('phirphi_opt: matrix only can take values R or P')
 
C Nullify pointers
      if (frstme) then
        nullify(Si,Vi,listed,needed)
        nullify(Ski,iano,iono)
        nullify(Pij,calculated)
        frstme = .false.
      endif

C Allocate arrays
      call re_alloc(listed,1,no,name='listed',routine='phirphi_opt')
      call re_alloc(needed,1,no,name='needed',routine='phirphi_opt')
      call re_alloc(Si,1,no,name='Si',routine='phirphi_opt')
      call re_alloc(Vi,1,no,name='Vi',routine='phirphi_opt')
      call re_alloc(iano,1,maxno,name='iano',routine='phirphi_opt')
      call re_alloc(iono,1,maxno,name='iono',routine='phirphi_opt')
      call re_alloc(Ski,1,2,1,maxkba,1,maxno,name='Ski',
     .              routine='phirphi_opt')
      norb = lmx2*nzetmx*nsemx
      call re_alloc(Pij,1,norb,1,norb,1,nspecies,name='Pij',
     .              routine='phirphi_opt')
      call re_alloc(calculated,1,norb,1,norb,1,nspecies,
     .              name='calculated',routine='phirphi_opt')

C Initialise quantities
      do io = 1,maxnh
        Sp(io) = 0.0d0
      enddo 
      do ka = 1,nspecies
        do io = 1,norb
          do jo = 1,norb
            calculated(jo,io,ka) = .false.
          enddo 
        enddo 
      enddo

C Find maximum range
      rmaxo = 0.0d0
      rmaxkb = 0.0d0
      do ia = 1,na
        is = isa(ia)
        do ikb = lastkb(ia-1)+1,lastkb(ia)
          ioa = iphKB(ikb)
          rmaxkb = max( rmaxkb, rcut(is,ioa) )
        enddo
        do io = lasto(ia-1)+1,lasto(ia)
          ioa = iphorb(io)
          rmaxo = max( rmaxo, rcut(is,ioa) )
        enddo
      enddo
      rmax = rmaxo + rmaxkb

C Correction due to the non-locality of the potential have to be 
C calculated for the momentum operator matrix elements

      if (matrix.eq.'P') then 

C Initialize arrayd Vi only once
        no = lasto(na)
        do jo = 1,no
          Vi(jo) = 0.0d0
          listed(jo) = .false.
          needed(jo) = .false.
        enddo

C Find out which orbitals are needed on this node
        do iio = 1,nuo
          call LocalToGlobalOrb(iio,Node,Nodes,io)
          needed(io) = .true.
          do j = 1,numh(iio)
            ind = listhptr(iio) + j
            jo = listh(ind)
            needed(jo) = .true.
          enddo
        enddo

C Loop on atoms with KB projectors
        do ka = 1,na
          ks = isa(ka)
          nkb = lastkb(ka) - lastkb(ka-1)
          if (nkb.gt.maxkba) then
            maxkba = nkb + 10
            call re_alloc(Ski,1,2,1,maxkba,1,maxno,copy=.true.,
     .        name='Ski',routine='phirphi_opt')
          endif

C Find neighbour atoms
          call mneighb( scell, rmax, na, xa, ka, 0, nna )

C Find neighbour orbitals
          nno = 0
          do ina = 1,nna
            ia = jna(ina)

            is = isa(ia)
            rki = sqrt(r2ij(ina))
            do io = lasto(ia-1)+1,lasto(ia)
              if (needed(io)) then
                call GlobalToLocalOrb(io,Node,Nodes,iio)
                ioa = iphorb(io)

C Find if orbital is within range
                within = .false.
                do ko = lastkb(ka-1)+1,lastkb(ka)
                  koa = iphKB(ko)
                  if ( rki .lt. rcut(is,ioa)+rcut(ks,koa) ) 
     .              within = .true.
                enddo

C Find overlap between neighbour orbitals and KB projectors
                if (within) then
                  nno = nno + 1
                  if (nno.gt.maxno) then
                    maxno = nno + 500
                    call re_alloc(iano,1,maxno,copy=.true.,
     .                name='iano',routine='phirphi_opt')
                    call re_alloc(iono,1,maxno,copy=.true.,
     .                name='iono',routine='phirphi_opt')
                    call re_alloc(Ski,1,2,1,maxkba,1,maxno,
     .                copy=.true.,name='Ski',routine='phirphi_opt')
                  endif
                  iono(nno) = io
                  iano(nno) = ia
                  ikb = 0
                  do ko = lastkb(ka-1)+1,lastkb(ka)
                    ikb = ikb + 1
                    koa = iphKB(ko)
                    do ix = 1,3
                     xinv(ix) = - xij(ix,ina)
                    enddo 
                    call matel('S', is, ks, ioa, koa, xinv,
     .                         Ski(1,ikb,nno), grSki)
              
                    sum = 0.0d0
                    if (abs(dk(1)).gt.tiny) then
                      call matel('X', is, ks, ioa, koa, xinv,
     .                           Sik, grSki)
                      sum = sum + Sik*dk(1) 
                    endif
                    if (abs(dk(2)).gt.tiny) then
                      call matel('Y', is, ks, ioa, koa, xinv,
     .                           Sik, grSki)
                      sum = sum + Sik*dk(2)
                    endif
                    if (abs(dk(3)).gt.tiny) then
                      call matel('Z', is, ks, ioa, koa, xinv,
     .                           Sik, grSki)
                      sum = sum + Sik*dk(3)
                    endif
                    Ski(2,ikb,nno) = sum
                  enddo
                endif
              endif
            enddo

          enddo

C Loop on neighbour orbitals
          do ino = 1,nno
            io = iono(ino)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            if (iio .gt. 0) then
              ia = iano(ino)
              if (ia .le. nua) then

C Scatter filter of desired matrix elements
                do j = 1,numh(iio)
                  ind = listhptr(iio) + j
                  jo = listh(ind)
                  listed(jo) = .true.
                enddo

C Find matrix elements with other neighbour orbitals
                do jno = 1,nno
                  jo = iono(jno)
                  if (listed(jo)) then

C Loop on KB projectors
                    ikb = 0
                    do ko = lastkb(ka-1)+1,lastkb(ka)
                      ikb = ikb + 1
                      koa = iphKB(ko)
                      epsk = epskb(ks,koa)
                      Sik = Ski(1,ikb,ino)
                      Sikr= Ski(2,ikb,ino)
                      Sjk = Ski(1,ikb,jno)
                      Sjkr= Ski(2,ikb,jno)
                      Vi(jo) = Vi(jo) + epsk * (Sik * Sjkr - Sikr * Sjk)
                    enddo
  
                  endif
                enddo

C Pick up contributions to H and restore Di and Vi
                do j = 1,numh(iio)
                  ind = listhptr(iio) + j
                  jo = listh(ind) 
                  Sp(ind) = Sp(ind) + Vi(jo)
                  Vi(jo) = 0.0d0
                  listed(jo) = .false.
                enddo
  
              endif
            endif
          enddo
        enddo

      endif  

C Initialize neighb subroutine 
      call mneighb( scell, 2.0d0*rmaxo, na, xa, 0, 0, nnia )

      do jo = 1,no
        Si(jo) = 0.0d0
      enddo

      do ia = 1,nua 
        
        is = isa(ia)
        call mneighb( scell, 2.0d0*rmaxo, na, xa, ia, 0, nnia )
           
        do io = lasto(ia-1)+1,lasto(ia)  
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          if (iio .gt. 0) then
            ioa = iphorb(io)
            do jn = 1,nnia 
              do ix = 1,3
                xinv(ix) = - xij(ix,jn)
              enddo 
              ja = jna(jn)
              rij = sqrt( r2ij(jn) )
              do jo = lasto(ja-1)+1,lasto(ja)
                joa = iphorb(jo)
                js = isa(ja)  
 
                if (rcut(is,ioa)+rcut(js,joa) .gt. rij) then  

                  if (matrix.eq.'R') then 
 
                    call matel('X', js, is, joa, ioa, xinv,
     .                         Sij, grSij )
                    Si(jo) = Sij*dk(1)  
                     
                    call matel('Y', js, is, joa, ioa, xinv,
     .                         Sij, grSij )
                    Si(jo) = Si(jo) + Sij*dk(2)  
 
                    call matel('Z', js, is, joa, ioa, xinv,
     .                         Sij, grSij )
                    Si(jo) = Si(jo) + Sij*dk(3) 
           
                    call matel('S', is, js, ioa, joa, xij(1,jn),
     .                         Sij, grSij )
                    Si(jo) = Si(jo) + Sij*(
     .                   xa(1,ia)*dk(1)
     .                 + xa(2,ia)*dk(2)
     .                 + xa(3,ia)*dk(3))  
   
                  else
                            
                    if (rij.lt.tiny) then 
C Perform the direct computation of the matrix element of the momentum 
C within the same atom
                      if (.not.calculated(joa,ioa,is)) then 
                        ilm1 = lofio(is,ioa)**2 + lofio(is,ioa) + 
     .                    mofio(is,ioa) + 1
                        ilm2 = lofio(is,joa)**2 + lofio(is,joa) + 
     .                    mofio(is,joa) + 1
                        call intgry(ilm1,ilm2,dintg2)
                        call intyyr(ilm1,ilm2,dintg1)
                        dintgmod1 = dintg1(1)**2 + dintg1(2)**2 + 
     .                    dintg1(3)**2
                        dintgmod2 = dintg2(1)**2 + dintg2(2)**2 + 
     .                    dintg2(3)**2
                        Sir0 = 0.0d0
                        if ((dintgmod2.gt.tiny).or.(dintgmod1.gt.tiny)) 
     .                      then 
                          dint1 = 0.0d0
                          dint2 = 0.0d0
                          npoints = int(max(rcut(is,ioa),rcut(is,joa))
     .                      /dx) + 2
                          do ir = 1,npoints
                            r = dx*(ir-1)
                            call rphiatm(is,ioa,r,phi1,dphi1dr)
                            call rphiatm(is,joa,r,phi2,dphi2dr)
                            dint1 = dint1 + dx*phi1*dphi2dr*r**2
                            dint2 = dint2 + dx*phi1*phi2*r
                          enddo 
C The factor of two because we use Ry for the Hamiltonian
                          Sir0 =
     .                 -2.0d0*(dk(1)*(dint1*dintg1(1)+dint2*dintg2(1))+
     .                         dk(2)*(dint1*dintg1(2)+dint2*dintg2(2))+
     .                         dk(3)*(dint1*dintg1(3)+dint2*dintg2(3))) 
                        endif 
                        Pij(ioa,joa,is) = - Sir0
                        Pij(joa,ioa,is) =   Sir0
                        calculated(ioa,joa,is) = .true.
                        calculated(joa,ioa,is) = .true.
                      endif 
                      Si(jo) = Pij(joa,ioa,is)

                    else
C Matrix elements between different atoms are taken from the 
C gradient of the overlap 
                      call matel('S', is, js, ioa, joa, xij(1,jn),
     .                           Sij, grSij )
C The factor of two because we use Ry for the Hamiltonian
                      Si(jo) =
     .                  2.0d0*(grSij(1)*dk(1)
     .             +           grSij(2)*dk(2)
     .             +           grSij(3)*dk(3))

                    endif 
  
                  endif 
                endif
              enddo
            enddo
            do j = 1,numh(iio)
              ind = listhptr(iio) + j
              jo = listh(ind)
              Sp(ind) = Sp(ind) + Si(jo)
              Si(jo) = 0.0d0 
            enddo
          endif
        enddo
      enddo

C Reduce size of arrays
      call re_alloc(calculated,1,1,1,1,1,1,name='calculated',
     .              routine='phirphi_opt')
      call re_alloc(Pij,1,1,1,1,1,1,name='Pij',
     .              routine='phirphi_opt')
      call re_alloc(Ski,1,2,1,1,1,1,name='Ski',routine='phirphi_opt')
      call re_alloc(iono,1,1,name='iono',routine='phirphi_opt')
      call re_alloc(iano,1,1,name='iano',routine='phirphi_opt')
      call re_alloc(Vi,1,1,name='Vi',routine='phirphi_opt')
      call re_alloc(Si,1,1,name='Si',routine='phirphi_opt')
      call re_alloc(needed,1,1,name='needed',routine='phirphi_opt')
      call re_alloc(listed,1,1,name='listed',routine='phirphi_opt')

C Stop timer
      call timer('phirphiopt',2)

      return
      end


      subroutine Intgry( ilm1, ilm2, dintg)
C *******************************************************************
C Returns a vector with the value of the integral of a spherical 
C harmonic with the gradient of another spherical harmonic.
C written by DSP. August 1999 (from subroutine ylmexp by J. Soler)
C ************************* INPUT ***********************************
C integer  ilm1                     : first spherical harmonic index: 
C                                     ilm1=L1*L1+L1+M1+1.
C integer  ilm2                     : second spherical harmonic index
C ************************* OUTPUT **********************************
C real*8   dintg(3)                 : Vector with the value of the 
C                                    (angular) integral of the product
C                                     Yl1m1*grad(Yl2m2)
C *******************************************************************
C external rlylm(lmax,rvec,y,grady) : Returns spherical harmonics,
C                                     times r**l
C *******************************************************************

      use precision, only: dp
      use alloc
      use spher_harm,   only : gauleg, rlylm

      implicit none

C Passed variables
      integer           ilm1, ilm2
      real(dp)          dintg(3)

C Internal variables
      integer
     .  lmax, l2, ix, isp, iz, nsp, im,
     .  maxlm, maxsp
      real(dp)
     .  phi, pi, theta, dl, r(3)

      integer,           save :: maxl = 8
      logical,           save :: frstme = .true.
      real(dp), pointer, save :: gry(:,:)
      real(dp), pointer, save :: w(:)
      real(dp), pointer, save :: wsp(:)
      real(dp), pointer, save :: x(:,:)
      real(dp), pointer, save :: y(:)
      real(dp), pointer, save :: z(:)

      external
     .  dot

C Nullify pointers and initialise arrays on first call
      if (frstme) then
        nullify(gry,w,wsp,x,y,z)
        maxlm = (maxl+1)*(maxl+1)
        maxsp = (maxl+1)*(2*maxl+1)
        call re_alloc(gry,1,3,0,maxlm,name='gry',routine='intgry')
        call re_alloc(w,1,maxl+1,name='w',routine='intgry')
        call re_alloc(wsp,1,maxsp,name='wsp',routine='intgry')
        call re_alloc(x,1,3,1,maxsp,name='x',routine='intgry')
        call re_alloc(y,0,maxlm,name='y',routine='intgry')
        call re_alloc(z,1,maxl+1,name='z',routine='intgry')
        frstme = .false.
      endif

C Find special points and weights for gaussian quadrature 
      if (max(ilm1,ilm2).eq.1) then 
        lmax = 1
      else
        dl = dble(max(ilm1,ilm2)) - 1.0d0 + 1.0d-2
        lmax = int(dsqrt(dl)) + 1
      endif

C Check array dimensions
      if (lmax.gt.maxl) then
        maxl = lmax
        maxlm = (maxl+1)*(maxl+1)
        maxsp = (maxl+1)*(2*maxl+1)
        call re_alloc(gry,1,3,0,maxlm,name='gry',routine='intgry')
        call re_alloc(w,1,maxl+1,name='w',routine='intgry')
        call re_alloc(wsp,1,maxsp,name='wsp',routine='intgry')
        call re_alloc(x,1,3,1,maxsp,name='x',routine='intgry')
        call re_alloc(y,0,maxlm,name='y',routine='intgry')
        call re_alloc(z,1,maxl+1,name='z',routine='intgry')
      endif
C
      call gauleg( -1.0d0, 1.0d0, z, w, lmax+1 )
      pi = 4.0d0*atan(1.0d0)
      nsp = 0
      do iz = 1,lmax+1
        theta = acos( z(iz) )
        do im = 0,2*lmax
          nsp = nsp + 1
          phi = im*2.0d0*pi/(2*lmax+1)
          x(1,nsp) = sin(theta)*cos(phi)
          x(2,nsp) = sin(theta)*sin(phi)
          x(3,nsp) = cos(theta)
          wsp(nsp) = w(iz)*(2.0d0*pi) / (2*lmax+1)
        enddo
      enddo

C Find the integrals
      if (ilm2.eq.1) then
        l2 = 0
      else
        dl = dble(ilm2) - 1.0d0 + 1.0d-2
        l2 = int(dsqrt(dl))
      endif
      do ix = 1,3
        dintg(ix) = 0.0d0
      enddo 
      if (l2.eq.0) return
      do isp = 1,nsp
        r(1) = x(1,isp)
        r(2) = x(2,isp)
        r(3) = x(3,isp)
        call rlylm( lmax, r, y, gry )
        do ix = 1,3
          dintg(ix) = dintg(ix) + wsp(isp)*y(ilm1-1)*
     .      (gry(ix,ilm2-1) - l2*y(ilm2-1)*x(ix,isp))
        enddo
      enddo
      
      end


      subroutine Intyyr( ilm1, ilm2, dintg)
C *******************************************************************
C Returns a vector with the value of the integral of the product of spherical 
C two spherical harmonics and the radial unit vector. 
C ************************* INPUT ***********************************
C integer  ilm1                     : first spherical harmonic index: 
C                                     ilm1=L1*L1+L1+M1+1.
C integer  ilm2                     : second spherical harmonic index
C ************************* OUTPUT **********************************
C real*8   dintg(3)                 : Vector with the value of the 
C                                    (angular) integral of the product
C                                     Yl1m1*Yl2m2*(x/r,y/r,z/r)
C *******************************************************************
C external rlylm(lmax,rvec,y,grady) : Returns spherical harmonics,
C                                     times r**l
C *******************************************************************

      use precision, only: dp
      use alloc
      use spher_harm,   only : gauleg, rlylm

      implicit none

C Passed variables
      integer           ilm1, ilm2
      real(dp)          dintg(3)

C Internal variables 
      integer
     .  lmax,ix, isp, iz, nsp, im, maxlm, maxsp
      real(dp)
     .  phi, pi, theta, dl, r(3)

      integer,           save :: maxl = 8
      logical,           save :: frstme = .true.
      real(dp), pointer, save :: gry(:,:)
      real(dp), pointer, save :: w(:)
      real(dp), pointer, save :: wsp(:)
      real(dp), pointer, save :: x(:,:)
      real(dp), pointer, save :: y(:)
      real(dp), pointer, save :: z(:)

      external
     .  dot

C Nullify pointers and initialise arrays on first call
      if (frstme) then
        nullify(gry,w,wsp,x,y,z)
        maxlm = (maxl+1)*(maxl+1)
        maxsp = (maxl+1)*(2*maxl+1)
        call re_alloc(gry,1,3,0,maxlm,name='gry',routine='intyyr')
        call re_alloc(w,1,maxl+1,name='w',routine='intyyr')
        call re_alloc(wsp,1,maxsp,name='wsp',routine='intyyr')
        call re_alloc(x,1,3,1,maxsp,name='x',routine='intyyr')
        call re_alloc(y,0,maxlm,name='y',routine='intyyr')
        call re_alloc(z,1,maxl+1,name='z',routine='intyyr')
        frstme = .false.
      endif

C Find special points and weights for Gaussian quadrature
      if (max(ilm1,ilm2).eq.1) then 
        lmax = 1
      else
        dl = dble(max(ilm1,ilm2)) - 1.0d0 + 1.0d-2
        lmax = int(dsqrt(dl)) + 1
      endif

C Check array dimensions
      if (lmax.gt.maxl) then
        maxl = lmax
        maxlm = (maxl+1)*(maxl+1)
        maxsp = (maxl+1)*(2*maxl+1)
        call re_alloc(gry,1,3,0,maxlm,name='gry',routine='intgry')
        call re_alloc(w,1,maxl+1,name='w',routine='intgry')
        call re_alloc(wsp,1,maxsp,name='wsp',routine='intgry')
        call re_alloc(x,1,3,1,maxsp,name='x',routine='intgry')
        call re_alloc(y,0,maxlm,name='y',routine='intgry')
        call re_alloc(z,1,maxl+1,name='z',routine='intgry')
      endif
C
      call gauleg( -1.0d0, 1.0d0, z, w, lmax+1 )
      pi = 4.0d0*atan(1.0d0)
      nsp = 0
      do iz = 1,lmax+1
        theta = acos(z(iz))
        do im = 0,2*lmax
          nsp = nsp + 1
          phi = im * 2.0d0*pi/(2*lmax+1)
          x(1,nsp) = sin(theta)*cos(phi)
          x(2,nsp) = sin(theta)*sin(phi)
          x(3,nsp) = cos(theta)
          wsp(nsp) = w(iz)*(2.0d0*pi)/(2*lmax+1)
        enddo
      enddo

C Find the integrals
      do ix = 1,3
        dintg(ix) = 0.0d0
      enddo 
      do isp = 1,nsp
        r(1) = x(1,isp)
        r(2) = x(2,isp)
        r(3) = x(3,isp)
        call rlylm( lmax, r, y, gry )
        do ix = 1,3
          dintg(ix) = dintg(ix) +
     .      wsp(isp)*y(ilm1-1)*y(ilm2-1)*x(ix,isp) 
        enddo
      enddo
 
      end
