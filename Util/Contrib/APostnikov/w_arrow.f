C...............................................................
C
      subroutine write_arrow(io1,is1,nbox,ivmin,ivmax,iev,
     .                       cc_ang,nat,iz,freq,disp)
C
C     write down atom coordinates and displacement vectors
C     in the format for Xcrysden
      implicit none
      integer io1,is1,nat,nbox,ibox,iev,ivmin,ivmax,iat,ii,jj,iz(nat)
      double precision cc_ang(3,3),coort(3),
     .       disp(3,nat,ivmin:ivmax),freq(ivmin:ivmax),
     .       btoang,dscal,twopi
      data btoang/0.529177/  !  Bohr to Angstroem
      data dscal /0.10/      !  A convenient scale factor for arrows' length

      write (io1,211) iev,freq(iev)
C --- header as for periodic structure
C     write (io1,'(A)') 'CRYSTAL' 
C     write (io1,'(A)') 'PRIMVEC' 
C     do ii=1,3
C       write (io1,'(3f16.9)') (cc_ang(ii,jj),jj=1,3)
C     enddo
C     write (io1,'(A)') 'PRIMCOORD' 
C     write (io1,'(2i5)') nat,1
C --- read coordinates of atoms in the box from  is1:
      rewind is1
C --- header as for molecule (= selected atoms in the box):
      write (io1,'(A)') 'ATOMS' 
      do ibox=1,nbox
        read  (is1,'(i4,3f20.8)') iat, (coort(jj),jj=1,3)
        write (io1,'(i4,3f12.7,2x,3f12.7)') iz(iat), 
     .              (coort(jj),jj=1,3),
     .              (disp(jj,iat,iev)*dscal,jj=1,3)
      enddo
      return

  211 format ('# ---- XSF block for ---- iev =',i6,
     .        '  freq = ',f14.6,' cm-1')
      end
