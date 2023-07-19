c
c $Id: pseudo.f,v 1.4 1999/02/26 14:26:46 wdpgaara Exp $
c
      subroutine Pseudo(pot_id,headline,ps_generator)
c
c     Generates the pseudopotential. The particular flavor is
c     determined by ps_generator.
c
      implicit none
c
      character pot_id*40, headline*79
      external ps_generator
c
      include 'radial.h'
      include 'orbital.h'
      include 'param.h'
      include 'elecpot.h'
c
C     .. Parameters ..
      double precision zero, one
      parameter (zero=0.D0,one=1.D0)
C     ..
C     .. Local Scalars ..
      integer i, llp, lp
C     ..
C     .. Local Arrays ..
      double precision ar(nrmax), br(nrmax)
C     ..
C     .. External Subroutines ..
      external ext, wf, wrapup, defined
      logical  defined
C     ..
c
      do 10 i = 1, 5
         indd(i) = 0
         indu(i) = 0
   10 continue
c
      if (ncore .eq. norb) return
      if (job .ne. 1 .and. job .ne. 2 .and. job .ne. 3) return
c
c
c  read rc(s),rc(p),rc(d),rc(f),cfac,rcfac
c
c  cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c
c  rcfac is used for the pseudocore cut off radius.  If set
c  to less than or equal to zero cfac is used.  cfac must be
c  set to greater than zero.
c
c     Optionally use free-format input here
c     The user should make sure that all six numbers
c     are present.
c
      if (defined('FREE_FORMAT_RC_INPUT')) then
         read(5,*) (rc_input(i),i=1,4), cfac, rcfac
      else
         read(5,9000) (rc_input(i),i=1,4), cfac, rcfac
      endif
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
      write(6,9010) nameat, pot_id, headline
 9010 format(//1x,a2,' pseudopotential generation: ',a,//,a,/)
c
c      start loop over valence orbitals
c
      ncp = ncore + 1
c
      do 220 i = ncp, norb
c
         lp = lo(i) + 1
         llp = lo(i)*lp
         if (down(i)) then
            if (indd(lp) .ne. 0) then
               write(6,9020) 'down', lp - 1
               call ext(800+lp)
            else
               indd(lp) = i
            end if
         else
            if (indu(lp) .ne. 0) then
               write(6,9020) 'up',lp - 1
               call ext(810+lp)
            else
               indu(lp) = i
            end if
         end if
 9020    format(//' error in pseudo - two ',a4,
     &            ' spin orbitals of the same ',
     &         /'angular momentum (',i1,') exist')
c
c
c      Find all electron wave function and its nodes and
c      extrema.
c
         call Wf(i,ar,br)
c
c      Find the pseudopotential
c
         call ps_generator(i,ar,br)
c
  220 continue
c
c  End loop over valence orbitals.
c
      call wrapup(pot_id)
c
      return
c
      end
