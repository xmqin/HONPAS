      subroutine pswf_store(arps,l,spin_down)
c
c     Stores the valence pseudowavefunction at the appropriate
c     place in the pseudowave data structures for further manipulation.
c
c     Alberto Garcia, after Javier Junquera, 2004-2005
c
      include "radial.h"
      include "ion.h"
      include "pseudowave.h"

      double precision     arps(*)            ! Pseudo-wavefunction
      integer              l                  ! angular momentum
      logical              spin_down          ! is it a down-one?

      integer              ir, lp

      lp = l + 1

      if (spin_down) then
         do ir = 1, nr
            pswfnrd(lp,ir) = arps(ir)
         enddo
      else
         do ir = 1, nr
            pswfnru(lp,ir) = arps(ir)
         enddo
      endif

      end


      
      
