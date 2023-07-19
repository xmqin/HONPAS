C
      subroutine etotal(nfirst,nlast)
c
      implicit none
c
      include 'orbital.h'
      include 'param.h'
      include 'energy.h'
c
c  etotal computes the total energy from the
c  electron charge density.
c
c  Only the orbitals between nfirst and nlast are considered.
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
C     .. Parameters ..
      double precision zero
      parameter (zero=0.D0)
C     ..
      integer nfirst, nlast
c
C     .. Local Scalars ..
      double precision vsum
      integer i
      character*2 id
C     ..
C     .. Intrinsic Functions ..
      intrinsic abs
C     ..
      etot(1) = zero
      etot(2) = zero
      etot(3) = zero
      do 10 i = nfirst, nlast
         etot(1) = etot(1) + zo(i)*ev(i)
         etot(2) = etot(2) + zo(i)*ek(i)
         etot(3) = etot(3) + zo(i)*ep(i)
   10 continue
c
c   kinetic energy
c
      etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c   potential energy
c
      etot(9) = etot(3) + etot(4) + etot(7)
c
c      total energy
c
      etot(10) = etot(1) - etot(4) - etot(5) + etot(7)
c
c   printout
c
      write(6,9000) nameat
 9000 format(//1x,a2,' output data for orbitals',/1x,28('-'),
     &      //' nl    s      occ',9x,'eigenvalue',4x,'kinetic energy',
     &      6x,'pot energy',/)
c
      id = " "
      do 20 i = nfirst, nlast
         if (i .ge. ncp) id="&v"
         write(6,9010) no(i), il(lo(i)+1), so(i), zo(i), ev(i),
     &     ek(i), ep(i), id
 9010    format(1x,i1,a1,f6.1,f10.4,3f17.8,2x,a2)
   20 continue
c
      write(6,'(a)') '---------------------------- &v'
      write(6,9020) (etot(i),i=1,10)
 9020 format(//' total energies',/1x,14('-'),
     &      //' sum of eigenvalues        =',f18.8,
     &      /' kinetic energy from ek    =',f18.8,
     &      /' el-ion interaction energy =',f18.8,
     &      /' el-el  interaction energy =',f18.8,
     &      /' vxc    correction         =',f18.8,
     &      /' virial correction         =',f18.8,
     &      /' exchange + corr energy    =',f18.8,
     &      /' kinetic energy from ev    =',f18.8,
     &      /' potential energy          =',f18.8,/1x,45('-'),
     &      /' total energy              =',f18.8)
c
      if (job .ge. 4 .or. abs(zsh) .gt. 0.00001D0) return
c
c   virial theorem
c
      vsum = 2*etot(8) + etot(9) + etot(6)
      write(6,9030) 2*etot(8), etot(9), etot(6), vsum
 9030 format(//' virial theorem(nonrelativistic)',/1x,14('-'),
     &      //' kinetic energy  *  2      =',f18.8,
     &      /' potential energy          =',f18.8,
     &      /' virial correction         =',f18.8,/1x,45('-'),
     &      /' virial sum                =',f18.8)
c
      return
c
      end
