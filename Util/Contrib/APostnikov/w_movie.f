C...............................................................
C
      subroutine write_movie(io2,is1,nbox,ivmin,ivmax,iev,
     .                       cc_ang,nat,iz,freq,disp,nsteps)
C
C     write down atom coordinates and animation of the vibration mode
C     in the format for Xcrysden
      implicit none
      integer io2,is1,nat,iat,nbox,ibox,iev,ivmin,ivmax,ii,jj,iz(nat),
     .        nsteps,istep
      double precision cc_ang(3,3),coort(3),
     .       disp(3,nat,ivmin:ivmax),freq(ivmin:ivmax),
     .       btoang,dscal,twopi
      data btoang/0.529177/  !  Bohr to Angstroem

      twopi=8.d0*atan(1.d0)
      write (io2,212) iev,freq(iev)
      write (io2,'(A,i5)') 'ANIMSTEPS',nsteps 
      if (nsteps.le.0) return
C       write (io2,'(A,i5)') 'CRYSTAL'
C       write (io2,'(A)') 'PRIMVEC' 
C       do ii=1,3
C         write (io2,'(3f16.9)') (cc_ang(ii,jj),jj=1,3)
C       enddo
C     assume fixed-cell animation, because we deal with phonons
      do istep=1,nsteps
        rewind is1
        write (io2,'(A,i5)') 'ATOMS',istep
        do ibox=1,nbox
          read  (is1,'(i4,3f20.8)') iat, (coort(jj),jj=1,3)
          write (io2,'(i4,3f12.7,2x,3f12.7)') iz(iat), 
     .                (coort(jj) +
     +                 sin(twopi*istep/nsteps) *
     *                 disp(jj,iat,iev),jj=1,3)
        enddo
      enddo
      return

  212 format ('# --- AXSF block for ---- iev =',i6,
     .        '  freq = ',f14.6,' cm-1')
      end

