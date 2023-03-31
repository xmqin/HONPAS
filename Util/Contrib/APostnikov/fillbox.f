C...............................................................
C
      subroutine fillbox(is1,obox,rbox,rinv,cc_ang,nat,coor_ang,nbox)
C
C     selects the atoms of the primitive cell, which -
C     including translations - fall into the selected box. 
C
C     Input     
C             is1 : read/write unit for temporary data
C            obox : origin of output box,
C            rbox : its three spanning vectors,
C            rinv : inversion of the matrix of spanning vectors
C          cc_ang : lattice vectors (in Ang)
C            nat  : number of atoms in the primitive cell 
C        coor_ang : atom coordinates in the cell (Cartesian, in Ang)
C     Output
C            nbox :   number of atoms in the output box 
C     The types and coordinates of selected atoms are written
C     into unit  is1. 
C
      implicit none

      integer nat, is1,limit
      double precision obox(3),rbox(3,3),rinv(3,3),coort(3),
     .                 cc_ang(3,3),coor_ang(3,nat)
      parameter (limit=5)  !  tried translations along each lattice vector
      integer iat,nbox,ishift,jshift,kshift,jj
      logical hit
      external hit

C --- write crystal structure data.
C     write selected atoms first into a scratch file (ios), for the case
C     there are zero. Then the label 'ATOMS' with no  atoms following
C     will crash XCrySDen.
      nbox = 0
      rewind is1
      do iat=1,nat
C       try displacements along lattice vectors to check whether
C       atoms fall into the given grid box. Hopefully [-limit,limit]
C       would be enough... 
        do ishift=-limit,limit   
        do jshift=-limit,limit   
        do kshift=-limit,limit   
          do jj=1,3
            coort(jj) = coor_ang(jj,iat) + 
     +             ishift*cc_ang(jj,1) +  
     +             jshift*cc_ang(jj,2) +  
     +             kshift*cc_ang(jj,3)   
          enddo
          if (hit(coort,obox,rinv)) then
            nbox = nbox + 1
            write (is1,'(i4,3f20.8)') iat, (coort(jj),jj=1,3)
          endif
        enddo
        enddo
        enddo
      enddo

      return
      end
