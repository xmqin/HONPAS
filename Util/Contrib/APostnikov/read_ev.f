C...............................................................
C
      subroutine read_ev(ii2,nat,ivmin,ivmax,evr,freq)
C
C     reads selected modes only from ii2 ('vector' file),
C     check for consistency, and keeps selected modes only
C     (ivmin to ivmax)
C
      implicit none
      integer ii2,ii,nat,iat,ivmin,ivmax,iev,idum
      double precision qq(3),dummy(3),small,
     .                 freq(ivmin:ivmax),evr(3,nat,ivmin:ivmax)
      data small/10.e-4/
C
      rewind (ii2)
      read (ii2,100,err=301)   !  leading line
      read (ii2,101,err=301) (qq(ii),ii=1,3)  
      if (sqrt(qq(1)**2+qq(2)**2+qq(3)**2).gt.small) then
        write (6,*) ' Cannot vizualize for q.ne.0 !'
        stop
      endif
      do iev=1,nat*3
        read (ii2,102,err=302,end=401) idum
        if (idum.ne.iev) then
          write (6,*) ' Stop in read_ev: expected eigenvector ',iev,
     .                ' but found ',idum
          stop
        endif
        if (iev.lt.ivmin.or.iev.gt.ivmax) then  ! -- skip these modes
          read (ii2,103,err=303) dummy(1)
          read (ii2,100)   !  Eigenvector, real part follows
          do iat=1,nat
            read (ii2,104,err=304) (dummy(ii),ii=1,3)
          enddo
          read (ii2,100)   !  Eigenvector, imag part follows
          do iat=1,nat
            read (ii2,104,err=305) (dummy(ii),ii=1,3)
          enddo
	    if ((dummy(1)**2+dummy(2)**2+dummy(3)**2).gt.small) then
              print *,' For iev =',iev,' iat =',iat,
     .                ' Imag. part of eigenvector is not zero!',
     .                ' - Unsuited for visualization'
              stop
            endif
       else  ! --- selected modes... 
          read (ii2,103,err=303) freq(iev)
          read (ii2,100)   !  Eigenvector, real part follows
          do iat=1,nat
            read (ii2,104,err=304) (evr(ii,iat,iev),ii=1,3)
          enddo
          read (ii2,100)   !  Eigenvector, imag part follows
          do iat=1,nat
            read (ii2,104,err=305) (dummy(ii),ii=1,3)
          enddo
	    if ((dummy(1)**2+dummy(2)**2+dummy(3)**2).gt.small) then
              print *,' For iev =',iev,' iat =',iat,
     .                ' Imag. part of eigenvector is not zero!',
     .                ' - Unsuited for visualization'
              stop
            endif
       endif  
      enddo  
      return
  301 print *,' Error in 1st or 2d line of .vector file'
      stop
  302 print *,' Error reading Eigenvector line, last iev=',iev
      stop
  303 print *,' Error reading Frequency line, iev=',iev
      stop
  304 print *,' Error reading real part of eigenvector ',iev,
     .        ' iat=',iat
      stop
  305 print *,' Error reading imag part of eigenvector ',iev,
     .        ' iat=',iat
      stop
  401 print *,' End of file while expecting eigenvector ',iev
      stop
  100 format()
  101 format(15x,3f12.6)
  102 format(14x,i6)
  103 format(14x,f13.6)
  104 format(3e12.4)
      end
