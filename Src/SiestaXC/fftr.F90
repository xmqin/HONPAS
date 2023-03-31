!!@LICENSE
!
!******************************************************************************
! MODULE fftr
! Handles 3-D fast Fourier transform of real functions
! Written by J.M.Soler. Feb.2008
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fftr2k : 3-D FFT transform from real to reciprocal space of a real function
! fftk2r : 3-D FFT transform from reciprocal to real space
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use m_fft3d, only: fft3d         ! 3D FFT of a complex function
! use mesh3D,  only: myMeshBox     ! Returns the mesh box of my processor
! use mesh3D,  only: setMeshDistr  ! Sets a new mesh distribution
! use mesh3D,  only: fftMeshDistr  ! Sets a mesh distr. for FFTs
! use mesh3D,  only: freeMeshDistr ! Frees a mesh distribution ID
! use mesh3D,  only: copyMeshData  ! Copies data in a box of mesh points
! use alloc,   only: de_alloc      ! De-allocation routine
! use alloc,   only: re_alloc      ! Re-allocation routine
! use sys,     only: die           ! Termination routine
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
! use precision,  only: gp=>grid_p      ! Real precision type of mesh arrays
!
!   EXTERNAL procedures used:
! none
!
!******************************************************************************
!
! SUBROUTINE fftr2k( nMesh, rDistr, f )

! 3-D fast Fourier transform from real to reciprocal space of a real function.
! On input, f must contain a real function in a box of grid points.
! On output, it contains the real part of its discrete Fourier transform 
! (if grid_index(k)<grid_index(-k)) or its imaginary part (if 
! grid_index(k)>grid_index(-k)) in a box of k mesh points. The real space box
! is determined by rDistr. The k-box is determined by kDistr, as returned by
! routine fftMeshDistr. See module mesh3D for details. In serial execution, 
! f contains the function at all mesh points in both r and k space
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)       : Total num. of mesh points in each axis
! integer rDistr         : Distribution ID of real-space mesh
!----------------------- INPUT and OUTPUT -------------------------------------
! real(gp) f(0:,0:,0:,:) : Functions to be Fourier transformed
!                          Array shape must be enough to contain
!                          both real- and k-space mesh boxes
!----------------------------- UNITS ------------------------------------------
! Units of f are arbitrary
!----------------------------- USAGE ------------------------------------------
!   Sample usage to find the electrostatic Hartree potential:
! use precision, only: gp=>grid_p
! use fftr,      only: fftk2r, fftr2k
! use mesh3D,    only: fftMeshDistr, myMeshBox, setMeshDistr
! integer:: kBox(2,3), kDistr, kMesh(3), maxMesh(3), nMesh(3), &
!           rBox(2,3), rDistr, rMesh(3)
! real(dp):: rCell(3,3), kCell(3,3), kVec(3), pi
! real(gp),allocatable:: f(:,:,:,:)
! ...Find rCell, kCell, nMesh, pi
! call setMeshDistr( rDistr, nMesh )
! call fftMeshDistr( nMesh, kDistr )
! call myMeshBox( nMesh, rDistr, rBox )
! call myMeshBox( nMesh, kDistr, kBox )
! rMesh(:) = rBox(2,:)-rBox(1,:)+1
! kMesh(:) = kBox(2,:)-kBox(1,:)+1
! maxMesh(:) = max( rMesh(:), kMesh(:) )
! allocate( f(0:maxMesh(1)-1,0:maxMesh(2)-1,0:maxMesh(3)-1,1) )
! f(0:rMesh(1)-1,0:rMesh(2)-1,0:rMesh(3)-1,1) = my_density
! call fftr2k( nMesh, rDistr, f )
! do k3 = kBox(1,3),kBox(2,3)
! do k2 = kBox(1,2),kBox(2,2)
! do k1 = kBox(1,1),kBox(2,1)
!   i1 = k1 - kBox(1,1)
!   i2 = k2 - kBox(1,2)
!   i3 = k3 - kBox(1,3)
!   kVec(:) = kCell(:,1)*k1 + kCell(:,2)*k2 + kCell(:,3)*k3
!   f(i1,i2,i3,1) = f(i1,i2,i3,1) * 4*pi/sum(kVec**2)
! end do
! end do
! end do
! call fftk2r( nMesh, rDistr, f )
! my_Vhartree = f(0:rMesh(1)-1,0:rMesh(2)-1,0:rMesh(3)-1,1)
!---------------------------- BEHAVIOUR ---------------------------------------
! If shape of f is not enough to contain my processor's mesh box in both real
!   and reciprocal space, it stops with an error warning.
!--------------------------- ALGORITHMS ---------------------------------------
! The real function array is copied into a complex array.
! The complex array data is redistributed from rDistr to kDistr using 
!   copyMeshData
! The complex array is Fourier transformed using fft3d
! The real imaginary parts of f, in reciprocal space, are stored in pairs of
!   (k,-k) values. The real part goes into whatever k or -k has a lower
!   mesh point index. The imaginary part to the higher-index one.
!
!******************************************************************************
!
! Subroutine fftk2r( nMesh, rDistr, f )
!
! 3-D fast Fourier transform from reciprocal to real space of a real function.
! On input, f must contain, in a box of reciprocal grid points, either the
! real part (if grid_index(k)<grid_index(-k)) or the imaginary part 
! (if grid_index(k)>grid_index(-k)) of a function that is real in r space. 
! On output, it contains its inverse Fourier transform, i.e. the function in 
! a box of real space mesh points.
!------------------------------ INPUT -----------------------------------------
! integer nMesh(3)    : Total num. of mesh points in each axis
! integer rDistr      : Distribution ID of real-space mesh
!----------------------- INPUT and OUTPUT -------------------------------------
! real(gp) f(0:,0:,0:,:) : Functions to be Fourier transformed
!                          Array shape must be enough to contain
!                          both real- and k-space mesh boxes
!----------------------------- UNITS ------------------------------------------
! Units of f are arbitrary
!----------------------------- USAGE ------------------------------------------
! See usage section in fftr2k
!---------------------------- BEHAVIOUR ---------------------------------------
! If shape of f is not enough to contain my processor's mesh box in both real
!   and reciprocal space, it stops with an error warning.
!--------------------------- ALGORITHMS ---------------------------------------
! The reciprocal-space values of f, at -k values, are obtained using 
!   copyMeshData, for the negative image of my processor's kBox
! The real and imaginary parts of f, at all my kBox points, is obtained from 
!   the values of f at k and -k.
! The complex function is back Fourier transformed to real space.
! The real part of the complex function, in real space, is redistributed and
!   copied to the output array within my processor's rBox, using copyMeshData.
!
!******************************************************************************

MODULE fftr

  use precision, only: dp            ! Double precision real type
  use precision, only: gp=>grid_p    ! Real precision type of mesh arrays
  use alloc,     only: de_alloc      ! De-allocates arrays
  use alloc,     only: re_alloc      ! Re-allocates arrays
  use m_fft3d,   only: fft3d         ! 3D FFT of a complex function
  use mesh3D,    only: associateMeshTask ! Associates commun. tasks and distr.
  use mesh3D,    only: myMeshBox     ! Returns the mesh box of my processor
  use mesh3D,    only: setMeshDistr  ! Sets a new mesh distribution
  use mesh3D,    only: fftMeshDistr  ! Sets a mesh distr. for FFTs
  use mesh3D,    only: freeMeshDistr ! Frees a mesh distribution ID
  use mesh3D,    only: copyMeshData  ! Copies data in a box of mesh points
  use sys,       only: die           ! Terminates execution
#ifdef DEBUG_XC
!  use debugXC,   only: udebug        ! File unit for debug output
#endif /* DEBUG_XC */

  implicit none

PUBLIC::    &
  fftk2r,   &! k->r (inverse) Fourier transform
  fftr2k     ! r->k (direct) Fourier transform

PRIVATE ! Nothing is declared public beyond this point

  integer,save:: k2r=-1         ! ID of task to change from k to r distribs.
  integer,save:: kDistr=-1      ! k-mesh distribution ID
  integer,save:: mk2k=-1        ! ID of task to change from -k to k data
  integer,save:: r2k=-1         ! ID of task to change from r to k distribs.
  integer,save:: theMesh(3)=0   ! Number of mesh points in each axis

CONTAINS

!------------------------------------------------------------------------------

subroutine fftr2k( nMesh, rDistr, f )

! 3-D fast Fourier transform from real to reciprocal space of a real function.

  implicit none

  integer, intent(in)   :: nMesh(3)    ! Total num. of mesh points in each axis
  integer, intent(in)   :: rDistr      ! Distribution ID of real-space mesh
  real(gp),intent(inout):: f(0:,0:,0:,:) ! Functions to be Fourier transformed
                                       ! Array shape must be enough to contain
                                       ! both real- and k-space mesh boxes

  character(len=*),parameter:: myName='fftr2k '
  integer:: ib1, ib2, ib3, ic, ic1, ic2, ic3, ip, jc1, jc2, jc3, jf, jp, &
            kBox(2,3), mkBox(2,3), myKmesh(3), myRmesh(3), rBox(2,3)
  real(gp),pointer:: fc(:,:,:,:)=>null()

! Find box of real-space mesh points corresponding to this node
  call myMeshBox( nMesh, rDistr, rBox )
  myRmesh(:) = rBox(2,:) - rBox(1,:) + 1

! Find box of k-space mesh points corresponding to this node
  call fftMeshDistr( nMesh, kDistr )
  call myMeshBox( nMesh, kDistr, kBox )
  myKmesh(:) = kBox(2,:) - kBox(1,:) + 1

! Associate communication tasks to mesh distributions
  call associateMeshTask( r2k, rDistr, kDistr )

! Check that input array is large enough
  do ic = 1,3
    if (size(f,ic) < max(myRmesh(ic),myKmesh(ic))) &
      call die(myName//'ERROR: size of input array f too small')
  end do

! Allocate complex array for fft
  call re_alloc( fc, kBox(1,1),kBox(2,1), &
                     kBox(1,2),kBox(2,2), &
                     kBox(1,3),kBox(2,3), 1,2, myName//'fc' )

! Loop on functions to be Fourier transformed
  do jf = 1,size(f,4)

    ! Copy real function to complex array
    call copyMeshData( nMesh, rDistr, f(:,:,:,jf:jf), kBox, fc(:,:,:,1:1), r2k )
    fc(:,:,:,2) = 0

    ! Perform Fourier transform
    call fft3d( fc, kDistr, nMesh, +1 )

    ! Store real and imaginary parts in pairs of k and -k vectors
    ! Real part goes to that of k,-k with lowest point index in cell
    ! Imaginary part goes to that with highest point index
    f(:,:,:,jf) = 0
    do ic3 = kBox(1,3),kBox(2,3)   ! Mesh indexes relative to cell origin
    do ic2 = kBox(1,2),kBox(2,2)
    do ic1 = kBox(1,1),kBox(2,1)
      ib1 = ic1 - kBox(1,1)        ! Mesh indexes relative to box origin
      ib2 = ic2 - kBox(1,2)
      ib3 = ic3 - kBox(1,3)
      ip = 1 + ic1 + nMesh(1)*ic2 + nMesh(1)*nMesh(2)*ic3  ! k point index
      ! Find indexes of -k vector in unit cell
      jc1 = modulo(-ic1,nMesh(1))
      jc2 = modulo(-ic2,nMesh(2))
      jc3 = modulo(-ic3,nMesh(3))
      jp = 1 + jc1 + nMesh(1)*jc2 + nMesh(1)*nMesh(2)*jc3  ! -k point index
      ! Fill f(k) with real or imaginary part of (k,-k) pair
      if (ip<=jp) then ! Store f(k) = Re(fc(k)) = Re(fc(-k))
        f(ib1,ib2,ib3,jf) = fc(ic1,ic2,ic3,1)
      else  ! (ip>jp) => Store f(k) = Im(fc(-k)) = -Im(fc(k))
        f(ib1,ib2,ib3,jf) = -fc(ic1,ic2,ic3,2)
      end if ! (ip<=jp)
    end do ! ic1
    end do ! ic2
    end do ! ic3

  end do ! jf

! Deallocate temporary array
  call de_alloc( fc, myname//'fc' )

end subroutine fftr2k

!------------------------------------------------------------------------------

subroutine fftk2r( nMesh, rDistr, f )

! 3-D fast Fourier transform from reciprocal to real space of a real function.

  implicit none

  integer, intent(in)   :: nMesh(3)    ! Total num. of mesh points in each axis
  integer, intent(in)   :: rDistr      ! Distribution ID of real-space mesh
  real(gp),intent(inout):: f(0:,0:,0:,:) ! Functions to be Fourier transformed
                                       ! Array shape must be enough to contain
                                       ! both real- and k-space mesh boxes

  character(len=*),parameter:: myName='fftk2r '
  integer:: ib1, ib2, ib3, ic, ic1, ic2, ic3, ip, jc1, jc2, jc3, jf, jp, &
            kBox(2,3), mkBox(2,3), myKmesh(3), myRmesh(3), rBox(2,3)
  real(gp),pointer:: fmk(:,:,:,:)=>null(), fc(:,:,:,:)=>null()

! Find box of real-space mesh points corresponding to this node
  call myMeshBox( nMesh, rDistr, rBox )
  myRmesh(:) = rBox(2,:) - rBox(1,:) + 1

! Find box of k-space mesh points corresponding to this node
  call fftMeshDistr( nMesh, kDistr )
  call myMeshBox( nMesh, kDistr, kBox )
  myKmesh(:) = kBox(2,:) - kBox(1,:) + 1

! Check that input array is large enough
  do ic = 1,3
    if (size(f,ic) < max(myRmesh(ic),myKmesh(ic))) &
      call die(myName//'ERROR: size of input array f too small')
  end do

! Allocate a complex array for FFT and another to hold f(-k)
  call re_alloc( fc, kBox(1,1),kBox(2,1), &
                     kBox(1,2),kBox(2,2), &
                     kBox(1,3),kBox(2,3), 1,2, myName//'fc' )
  if (kDistr/=0) then  ! f distributed over parallel processors
    mkBox(1,:) = -kBox(2,:)
    mkBox(2,:) = -kBox(1,:)
    call re_alloc( fmk, mkBox(1,1),mkBox(2,1), &
                        mkBox(1,2),mkBox(2,2), &
                        mkBox(1,3),mkBox(2,3), 1,1, myName//'fmk' )
  end if ! (kDistr/=0)

! Associate communication tasks to mesh distributions
  call associateMeshTask( k2r, rDistr, kDistr )
  call associateMeshTask( mk2k, kDistr )

! Loop on functions to be transformed
  do jf = 1,size(f,4)

    ! Get f vales at -k mesh points
    if (kDistr/=0) &
      call copyMeshData( nMesh, kDistr, f(:,:,:,jf:jf), mkBox, fmk, mk2k )

    ! Retrieve real and imaginary parts from pairs of k and -k vectors
    do ic3 = kBox(1,3),kBox(2,3)   ! Mesh indexes relative to cell origin
    do ic2 = kBox(1,2),kBox(2,2)
    do ic1 = kBox(1,1),kBox(2,1)
      ib1 = ic1 - kBox(1,1)        ! Mesh indexes relative to box origin
      ib2 = ic2 - kBox(1,2)
      ib3 = ic3 - kBox(1,3)
      ip = 1 + ic1 + nMesh(1)*ic2 + nMesh(1)*nMesh(2)*ic3
      ! Find indexes of -k vector in unit cell
      jc1 = modulo(-ic1,nMesh(1))
      jc2 = modulo(-ic2,nMesh(2))
      jc3 = modulo(-ic3,nMesh(3))
      jp = 1 + jc1 + nMesh(1)*jc2 + nMesh(1)*nMesh(2)*jc3
      ! Retrieve real and imaginary parts from (k,-k) pair
      if (ip==jp) then      ! -k=k+G  =>  fc(k) is real
        fc(ic1,ic2,ic3,1) = f(ib1,ib2,ib3,jf)
        fc(ic1,ic2,ic3,2) = 0
      else if (ip<jp) then  ! fc(k) = cmplx(f(k),f(-k))
        fc(ic1,ic2,ic3,1) = f(ib1,ib2,ib3,jf)
        if (kDistr==0) then ! f array not distributed
          fc(ic1,ic2,ic3,2) = f(jc1,jc2,jc3,jf)
        else                ! f array distributed over parallel processors
          fc(ic1,ic2,ic3,2) = fmk(-ic1,-ic2,-ic3,1)
        end if
      else ! (ip>jp)   =>     fc(k) = conjg(fc(-k)) = conjg(f(-k),f(k))
        if (kDistr==0) then
          fc(ic1,ic2,ic3,1) = f(jc1,jc2,jc3,jf)
        else
          fc(ic1,ic2,ic3,1) = fmk(-ic1,-ic2,-ic3,1)
        end if
        fc(ic1,ic2,ic3,2) = -f(ib1,ib2,ib3,jf)
      end if ! (ip==jp)
    end do ! ic1
    end do ! ic2
    end do ! ic3

    ! Perform inverse Fourier transform
    call fft3d( fc, kDistr, nMesh, -1 )

    ! Copy real function from complex to output array
    call copyMeshData( nMesh, kDistr, fc(:,:,:,1:1), rBox, f(:,:,:,jf:jf), k2r )

  end do ! jf

! Deallocate temporay arrays
  if (kDistr/=0) call de_alloc( fmk, myName//'fmk' )
  call de_alloc( fc, myname//'fc' )

end subroutine fftk2r

END MODULE fftr

