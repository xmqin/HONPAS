C...............................................................
C
      subroutine write_axsf2(ii2,io1,nstep,mdfirst,mdlast,mdstep,
     .                       nat,nz,mdmod,varcel,cc_ang,coord,veloc)
C
C     writes output info into the AXSF file for the case 
C     when no output box is selected.
C     The AXSF is as for periodic structure, with SIESTA native box,
C     either from MD (if it is written with variable cell), or from XV.
C     Read again nstep MD steps from ii2,
C     write those from mdfirst till mdlast, step mdstep.
C
      implicit none
      integer ii2,io1,mdfirst,mdlast,mdstep,nat,nz(nat),mdmod,
     .        ii,jj,istep,nstep,iat,idum
      double precision cc_ang(3,3),cc_bohr(3,3),cc_velo(3,3),cc_md(3,3),
     .                 coord(3,nat),veloc(3,nat),b2ang,coorf(3),coorx(3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      logical varcel
      character lett*1,symbol*2

      rewind (ii2) 
      rewind (io1) 
      write (io1,"('ANIMSTEPS',i6)") (mdlast-mdfirst)/mdstep+1
      if (mdmod.eq.1) then  ! (read data from MD, all in Bohr)
        write (io1,"('CRYSTAL')")
        if (.not.varcel) then
C         write cell vectors from XV (only once)
          write (io1,"('PRIMVEC')")
          do ii=1,3
            write (io1,'(3x,3f18.9)') (cc_ang(jj,ii),jj=1,3)
          enddo
        endif
        do istep = 1,nstep
          read (ii2) idum,coord,veloc            ! -- read atom coordinates...
          if (varcel) read (ii2) cc_bohr,cc_velo ! -- and cell only if (varcel)
          if ((istep.ge.mdfirst).and.
     .        (istep.le.mdlast).and.
     .        (mod(istep-mdfirst,mdstep).eq.0)) then
            if (varcel) then
              cc_ang = cc_bohr*b2ang
              write (io1,"('PRIMVEC',i8)") (istep-mdfirst)/mdstep+1
              do ii=1,3
                write (io1,'(3x,3f18.9)') (cc_ang(jj,ii),jj=1,3)
              enddo
            endif  ! if (varcel) 
            write (io1,"('PRIMCOORD',i6)") (istep-mdfirst)/mdstep+1
            write (io1,"(i6,'  1')") nat
            do iat=1,nat
              write (io1,'(i4,3f13.6)') nz(iat),
     .                                  (coord(ii,iat)*b2ang,ii=1,3)
C             coordinates in MD are in Bohr, and converted to Ang
            enddo
          endif  !  if write this MD step
        enddo  !  do istep

      elseif (mdmod.eq.2) then  ! (read data from MD_CAR, coord. in Ang)
        write (io1,"('CRYSTAL')")
        if (.not.varcel) then
C         write cell vectors from XV (only once)
          write (io1,"('PRIMVEC')")
          do ii=1,3
            write (io1,'(3x,3f18.9)') (cc_ang(jj,ii),jj=1,3)
          enddo
        endif
        do istep = 1,nstep
          read (ii2,'(a1)') lett  !  "---" // trim(slabel) //"---"
          read (ii2,'(a1)') lett  !  "       1.0"
          read (ii2,'(3f16.9)') cc_md  ! -- lattice vectors in Ang (each step)
          read (ii2,'(a1)') lett  !  List of atoms: modify!!!
          read (ii2,'(a1)') lett  !  "Direct"
          if ((istep.ge.mdfirst).and.
     .        (istep.le.mdlast).and.
     .        (mod(istep-mdfirst,mdstep).eq.0)) then
            if (varcel) then
              write (io1,"('PRIMVEC',i8)") (istep-mdfirst)/mdstep+1
              write (io1,'(3x,3f18.9)') cc_md
            endif  ! if (varcel) 
            write (io1,"('PRIMCOORD',i6)") (istep-mdfirst)/mdstep+1
            write (io1,"(i6,'  1')") nat
            do iat=1,nat
              read (ii2,'(3f16.9)') coorf ! -- fractional coordinates 
! ---         reconstruct cartesian coordinates in Ang:
              coorx(:) = coorf(1)*cc_md(:,1) +
     +                   coorf(2)*cc_md(:,2) +
     +                   coorf(3)*cc_md(:,3)
              write (io1,'(i4,3f13.6)') nz(iat),coorx
            enddo
          else   !  skip this MD step 
            do iat=1,nat
              read (ii2,'(3f16.9)') coorf ! -- dummy read, not used 
            enddo
          endif  !  if write this MD step
        enddo  !  do istep

      elseif (mdmod.eq.3) then  ! (read data from ANI, coord. in Ang)
        write (io1,"('CRYSTAL')")
C       write cell vectors from XV (only once)
        write (io1,"('PRIMVEC')")
        do ii=1,3
          write (io1,'(3x,3f18.9)') (cc_ang(jj,ii),jj=1,3)
        enddo
        do istep = 1,nstep
          read (ii2,301) nat
          do iat=1,nat
            read (ii2,'(a2,2x,3f12.6)') symbol,(coord(ii,iat),ii=1,3)
C           coordinates in ANI are in Angstroem, and so written to AXSF
          enddo
          if ((istep.ge.mdfirst).and.
     .        (istep.le.mdlast).and.
     .        (mod(istep-mdfirst,mdstep).eq.0)) then
            write (io1,"('PRIMCOORD',i6)") (istep-mdfirst)/mdstep+1
            write (io1,"(i6,'  1')") nat
            do iat=1,nat
              write (io1,'(i4,3f13.6)')  nz(iat),(coord(ii,iat),ii=1,3)
            enddo
          endif  !  if write this MD step
        enddo

      endif

      return
  301 format(i5,/)

      end
