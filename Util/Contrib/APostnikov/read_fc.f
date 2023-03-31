C...............................................................
C
      subroutine read_fc(ii1,nat,fc)
C
C     reads force constants file .FC from ii1 
C
      implicit none
      integer ii1,nat,ix,ipm,jat,iline
      double precision fc(:,:,:,:,:)
C
      rewind (ii1)
      iline = 0
      do iat=1,nat   !   atom which is displaced
      do ix =1,3     !   Cartesian displacements
      do ipm =1,2    !   displace + and -
        do jat=1,nat   !   atom at which force is induced
          iline = iline + 1
          read (ii1,101,err=301,end=302) (fc(jx,jat,ix,iat,ipm),jx=1,3)
        enddo
      enddo
      enddo
      enddo
      return

  301 continue
      print *,' Error reading FC file, line ',iline
      stop
  302 continue
      print *,' Unexpected end of FC file, iat,ix,ipm,jat=',
     .         iat,ix,ipm,jat     
      stop
  101 format (3f15.7)
      end

