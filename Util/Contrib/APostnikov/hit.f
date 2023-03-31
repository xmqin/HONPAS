C ...........................................................
C
      logical function hit(coort,obox,rinv)
C
C     checks whether translated vector coort()
C     falls inside a box with the origin obox(),
C     spanned by three vectors, the inverse of which are rinv()
C
C     Note that rbox of the main program contains COLUMNWISE the vectors 
C     which span the box,
C     hence rinv contains LINEWISE the vectors of the reciprocal box.
C     The orthogonality relation:
C     Sum_k [rinv(i,k)*rbox(k,j)] = delta(i,j)
C
      implicit none
      integer ii,jj
      double precision coort(3),obox(3),rinv(3,3),rela
C
C     The test is done by transforming the (cartesian) displaced
C     coordinate coort-obox into relative coordinates within the box.
C     Then check whether all coordinates fall within [0,1]
C
      hit = .false.
      do ii=1,3
        rela = 0.0
        do jj=1,3
          rela = rela + rinv(ii,jj)*(coort(jj)-obox(jj))
        enddo
        if (rela.lt.0..or.rela.gt.1.) return
      enddo
      hit=.true.
      return
      end
