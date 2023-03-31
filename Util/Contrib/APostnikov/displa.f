C...............................................................
C
      subroutine displa(nat,ivmin,ivmax,mass,dstmin,evr,disp)
C     Recovers displacements for each atom in each mode
C
      implicit none
      integer iat,nat,ii,iev,ivmin,ivmax
      double precision mass(nat),freq(ivmin:ivmax),
     .      evr(3,nat,ivmin:ivmax),disp(3,nat,ivmin:ivmax),
     .      dlen,dstmin,dmax,dscale,dfact

  101 write (6,701,advance="no")
      read (5,*,err=101) dscale

      do iev=ivmin,ivmax
      write (6,*)'  iev=',iev
        dmax = 0.d0
        do iat=1,nat
          dlen = 0.d0
          do ii=1,3
            disp(ii,iat,iev)=evr(ii,iat,iev)/sqrt(mass(iat))
	    dlen = dlen + disp(ii,iat,iev)**2
          enddo
          dlen = sqrt(dlen) 
          if (dlen.gt.dmax) dmax=dlen
C         write (6,301) iat,(evr(ii,iat,iev),ii=1,3),
C    .                      (disp(ii,iat,iev),ii=1,3)
C 301   format(' iat=',i4,' evr=',3d10.4,' disp=',3e10.4)
        enddo   !  do iat=1,nat
C   scale displacement to a convenient length;
C   take some fraction of smallest interatomic distance dstmin:
        dfact = dstmin * dscale / dmax
        do iat=1,nat
          do ii=1,3
            disp(ii,iat,iev)=disp(ii,iat,iev)*dfact
          enddo
        enddo
      enddo

      return
  701 format (' Specify magnitude of max. displacement',
     .        ' in terms of minimal interatomic distance',/
     .        ' (something like 0.2 - 0.5; or negative : ')
      end
