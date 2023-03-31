! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

      program plsts

! Simple explorer of STS data files (wip)
!     This is too simple for now. Ideas are welcome.

      use m_gridfunc, only: gridfunc_t, read_gridfunc
      use m_gridfunc, only: grid_p, get_values_at_point

       implicit none

       type(gridfunc_t) :: gf

       integer, parameter :: dp = selected_real_kind(10,100)

! Internal variables
       character
     .   name*75, fname*80, fname_aux*80, oname*80
       integer
     .      i, ip, is, j, mesh(3), np, idum, nt, iz, iy, i1, i2, ne
       integer :: k1, k2, ic, ix
       real ::    x, y, e
       integer :: iostat, n3
       real(dp) :: cell(3,3), origin(3)
       real(dp) :: dxdm(2,2)  ! 2D plane grid vectors

       real(grid_p), allocatable :: vals(:)
       real(dp) :: zmax
       real(dp) :: point(3)
       
!     Read plot data
       print *, "Root of filename: "
       read(5,*) name

       fname = trim(name)//'.STS'
       fname_aux = trim(name)//'.STS_AUX'
       write(6,*) 'Reading STS data from file ',trim(fname)

       call read_gridfunc(fname,gf)
       
       np = product(gf%n(1:3))
       ne = gf%nspin
       cell = gf%cell
       mesh = gf%n
       origin = gf%origin
       
       write(6,*)
       write(6,*) 'Cell vectors (bohr)'
       write(6,*)
       write(6,*) cell(:,1)
       write(6,*) cell(:,2)
       write(6,*)
       write(6,*) 'Grid mesh: ',mesh(1),'x',mesh(2)
       write(6,*)
       write(6,*) 'ne = ', ne
       write(6,*)
       write(6,"(a,3f10.5)") 'Box origin (bohr): ', origin(1:3)

       ! Find the values at an arbitrary point
       ! Use z=0 to make sure that we get some info
       ! In the general case, one should be careful with possible
       ! assumptions about periodicity. Hence check below, which
       ! could be in the routine itself.

       print *, "Enter coordinates of sample point (fractional):"
       read(*,*) point(1:3)
       !!!                     point = [ 0.3_dp, 0.5_dp, 0.0_dp ]
       if (.not. gf%is_periodic(3)) then
          zmax = (gf%n(3) - 1) * (1.0_dp / gf%n(3))
          if (point(3) > zmax) then
             stop 'z coord beyond last meaningful plane'
          endif
       endif
       allocate ( vals(ne) )
       call get_values_at_point(gf,point,vals)
       
       open(1,file=fname_aux,form="formatted", status="old",
     $      position="rewind",action="read")
       open(2,file=trim(name)//'.sts.dat',form="formatted",
     $      status="unknown", position="rewind",action="write")

       write(6,*) 'Writing STS I(E) data at point ' //
     $            ' in file: ' // trim(name)//'.sts.dat'

       write(2,'(a)') '#  Energy (eV)  I(E)  '
       do i = 1, ne
          read(1,*) idum, e
          write(2,*) e, vals(idum)
       enddo
       close(1)
       close(2)
       deallocate(vals)

       oname = trim(name)//'.STM.Sample.CH'

       ! Write 2D file info
       
       OPEN( unit=2, file=oname )
       write(6,*) 'Writing STM(E) contour image in file ', trim(oname)

          DO IC = 1,2
             DO IX = 1,2
                DXDM(IX,IC) = CELL(IX,IC) / MESH(IC)
             ENDDO
          ENDDO

       ! Write an example of contour, at an arbitrary energy
       do k1 = 0, mesh(1) - 1
          do k2 = 0, mesh(2) - 1
             x = dxdm(1,1) * k1 + dxdm(1,2) * k2
             y = dxdm(2,1) * k1 + dxdm(2,2) * k2
             write(2,*) x, y, gf%val(k1+1,k2+1,1,ne/2)
          enddo
          write(2,*)
       enddo
       close(2)
          
      end

