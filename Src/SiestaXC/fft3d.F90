!!@LICENSE
!
!******************************************************************************
! MODULE m_fft3d
! 3-D fast complex Fourier transform
! Written by J.M.Soler (Feb.2008).
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fft3d : 3-D complex FFT
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,    only: die           ! Termination routine
! use alloc,  only: de_alloc      ! De-allocation routine
! use alloc,  only: re_alloc      ! Re-allocation routine
! use fft1d,  only: gpfa          ! 1-D fast Fourier transform
! use fft1d,  only: setgpfa       ! Initializes gpfa routine
! use mesh3D, only: copyMeshData  ! Copies data in a box of mesh points
! use mesh3D, only: fftMeshDistr  ! Sets a mesh distr. for FFTs
! use mesh3D, only: myMeshBox     ! Returns the mesh box of my processor
! use mesh3D, only: redistributeMeshData ! Changes data distribution
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
! use precision,  only: gp=>grid_p      ! Real precision type of mesh arrays
!
!   EXTERNAL procedures used:
! gpfa    : 1D complex FFT
! setgpfa : Initializes gpfa
! timer   : CPU time counter
!
!******************************************************************************
!
! SUBROUTINE fft3d( dat, meshDistr, nMesh, r2k )
!
! 3D complex FFT
!------------------------------ INPUT -----------------------------------------
! integer meshDistr : Mesh distribution ID of dat array
! integer nMesh(3)  : Mesh divisions of each cell axis
! integer r2k       : Direction of Fourier transform:
!                       r2k=+1 => r->k; r2k=-1 => k->r
!----------------------- INPUT and OUTPUT -------------------------------------
! real(gp) dat(:,:,:,:) : Complex data to be Fourier transformed:
!                           dat(:,:,:,1) = Real part
!                           dat(:,:,:,2) = Imaginary part
!                         The dat array must be declared target.
!----------------------------- UNITS ------------------------------------------
! Units of dat are arbitrary
!----------------------------- USAGE ------------------------------------------
! use precision, only: gp=>grid_p
! use mesh3D,    only: fftMeshDistr, myMeshBox, setMeshDistr
! use m_fft3d,   only: fft3d
! integer:: myBox(2,3), myDistr, nMesh(3)
! real(gp),allocatable:: dat(:,:,:,:)
! ...Find nMesh
! call setMeshDistr( myDistr, nMesh )
!    Or, alternatively
! call fftMeshDistr( nMesh, myDistr )
! call myMeshBox( nMesh, myDistr, myBox )
! allocate( dat(myBox(1,1):myBox(2,1),   &
!               myBox(1,2):myBox(2,2),   &
!               myBox(1,3):myBox(2,3),2) )
! ...Find dat in myBox
! call fft3d( dat, myDistr, nMesh, +1 )
!---------------------------- BEHAVIOUR ---------------------------------------
! It stops with an error message in the following cases:
! - The dimension of dat array in indexes 1-3 is not enough to contain myBox 
!   in meshDistr
! - The dimension of dat in index 4 is different from 2
!--------------------------- ALGORITHMS ---------------------------------------
! 1) dat array is redistributed into kDistr, as returned by fftMeshDistr. 
! 2) dat array is redistributed again into axisDistr(1), as returned from 
!    fftMeshDistr. In this mesh distribution, each processor's box spans the
!    whole x axis (or, more properly, the whole first cell vector).
! 3) dat is Fourier transformed in the x axis
! 4) dat is redistributed again into kDistr
! 5) Steps 2-4 are repeated for the y and z axes.
! 6) dat is redistributed bak into the input meshDistr.
! In serial execution, none of these data redistributions really occur. 
!   Rather, pointers are used to point the 'redistributed' arrays to the
!   original ones.
!
!******************************************************************************

MODULE m_fft3d

  ! Used module procedures
  use sys,    only: die          ! Terminates execution
  use alloc,  only: de_alloc     ! Deallocates arrays
  use alloc,  only: re_alloc     ! (Re)allocates arrays
  use fft1d,  only: gpfa         ! 1-D fast Fourier transform
  use fft1d,  only: setgpfa      ! Initializes gpfa routine
  use mesh3D, only: associateMeshTask ! Associates commun. tasks to distr.
  use mesh3D, only: copyMeshData ! Copies a box of mesh data
  use mesh3D, only: fftMeshDistr ! Sets mesh distributions for FFTs
  use mesh3D, only: myMeshBox    ! Returns the mesh box of my node
  use mesh3D, only: redistributeMeshData ! Changes data distribution
  use m_timer,only: timer_start  ! Start counting CPU time
  use m_timer,only: timer_stop   ! Stop counting CPU time

  ! Used module parameters
  use precision,    only: dp           ! Real double precision type
  use precision,    only: gp=>grid_p   ! Real type of mesh array data
#ifdef DEBUG_XC
  use debugXC,      only: udebug    ! File unit for debug output
#endif /* DEBUG_XC */

  implicit none

PUBLIC:: &
  fft3d   ! 3D complex FFT

PRIVATE ! Nothing is declared public beyond this point

CONTAINS

!******************************************************************************

subroutine fft3d( dat, meshDistr, nMesh, r2k )

  implicit none

  ! Passed arguments
  real(gp),target,intent(inout):: dat(0:,0:,0:,:) ! Data to be Fourier 
                                 ! transformed (one complex function):
                                 !   dat(:,:,:,1) = Real part
                                 !   dat(:,:,:,2) = Imaginary part
  integer,intent(in):: meshDistr ! Mesh distribution ID of dat array
  integer,intent(in):: nMesh(3)  ! Mesh divisions of each cell axis
  integer,intent(in):: r2k       ! Direction of Fourier transform:
                                 ! r2k=+1 => r->k; r2k=-1 => k->r

  ! Internal parameters, variables, and arrays
  character(len=*),parameter:: myName = 'fft3d '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer,save:: a2my(3)=-1, aDistr(3)=-1, io2my=-1, my2a(3)=-1, &
                 myDistr=-1, my2io=-1, oldMesh(3)=0
  logical:: finished, newMesh
  integer:: aBox(2,3), aMesh(3), datShape(4), i1, i2, i3, ia, iter, &
            maxTrigs, meshBox(2,3), myBox(2,3), nTrigs(3)
  real(dp),pointer,save:: trigs(:,:)=>null()
  real(gp),pointer:: myDat(:,:,:,:)=>null(), aDat(:,:,:,:)=>null()

  ! Start time counter
  call timer_start( myName )

  ! Get my node's mesh box in input mesh distribution
  call myMeshBox( nMesh, meshDistr, meshBox )

  ! Check dat array shape
  datShape = shape(dat)
  if (any(datShape(1:3) < meshBox(2,:)-meshBox(1,:)+1)) &
    call die(errHead//'dat shape inconsistent with meshDistr')
  if (datShape(4)/=2) &
    call die(errHead//'incorrect shape of dat array')

  ! Check if this is a new mesh and store it for later calls
  newMesh = any(nMesh/=oldMesh)
  oldMesh = nMesh

  ! Set trigs array, used by 1D-FFT routine gpfa
  if (newMesh) then
    maxTrigs = 256
    do iter = 1,2
      call re_alloc( trigs, 1,maxTrigs, 1,3, myName//'trigs' )
      do ia=1,3  ! Loop on cell axes
        call setgpfa( trigs(:,ia), maxTrigs, nTrigs(ia), nMesh(ia) )
      end do ! ia
      if (all(nTrigs<=maxTrigs)) then
        finished = .true.
        exit ! iter loop
      else
        finished = .false.
        maxTrigs = maxval(nTrigs)
      end if
    end do ! iter
    if (.not.finished) &
      call die(errHead//'trigs iteration not converged')
  end if ! (newMesh)

  ! (Re)set FFT mesh distributions myDistr and aDistr.
  ! myDistr is a homogeneous 3D mesh distribution.
  ! aDistr(1:3) are three homogen. 2D distribs. In each one,
  ! all the node's mesh boxes span one of the whole cell axes
  if (newMesh) call fftMeshDistr( nMesh, myDistr, aDistr )

#ifdef DEBUG_XC
  ! Check that aDistr mesh boxes span whole cell axes
  if (newMesh) then
    do ia = 1,3
      call myMeshBox( nMesh, aDistr(ia), aBox )
      aMesh(:) = aBox(2,:) - aBox(1,:) + 1
      if (aMesh(ia)/=nMesh(ia)) call die(errHead//'incorrect aBox')
    end do ! ia
  end if ! (newMesh)
#endif /* DEBUG_XC */

  ! Associate communication tasks to distributions
  call associateMeshTask( io2my, meshDistr, myDistr )
  call associateMeshTask( my2io, meshDistr, myDistr )
  do ia = 1,3
    call associateMeshTask( my2a(ia), myDistr, aDistr(ia) )
    call associateMeshTask( a2my(ia), myDistr, aDistr(ia) )
  end do ! ia

  ! Redistribute data over myDistr
  ! Notice: if (meshDistr==myDistr) redistributeMeshData may
  ! just direct pointer myDat=>dat, without reallocating it
  call redistributeMeshData( meshDistr, dat, myDistr, &
                             dstData=myDat, task=io2my)

  ! FFT in first axis. Arguments of gpfa are:
  ! gpfa( realData(*), imagData(*), trigs(*),  
  !       dataSpan, vectorSpan, vectorSize, numVectors, fftDir )
  call myMeshBox( nMesh, aDistr(1), aBox )
  aMesh(:) = aBox(2,:) - aBox(1,:) + 1
  call redistributeMeshData( myDistr, myDat, aDistr(1), &
                             dstData=aDat, task=my2a(1) )
  call gpfa( aDat(:,:,:,1), aDat(:,:,:,2), trigs(:,1), &
             1, aMesh(1), aMesh(1), aMesh(2)*aMesh(3), -r2k )
  call redistributeMeshData( aDistr(1), aDat, myDistr, &
                             dstData=myDat, task=a2my(1) )
      
  ! FFT in second axis
  call myMeshBox( nMesh, aDistr(2), aBox )
  aMesh(:) = aBox(2,:) - aBox(1,:) + 1
  call redistributeMeshData( myDistr, myDat, aDistr(2), &
                             dstData=aDat, task=my2a(2) )
  do i3 = 0,aMesh(3)-1
    call gpfa( aDat(:,:,i3,1), aDat(:,:,i3,2), trigs(:,2), &
               aMesh(1), 1, aMesh(2), aMesh(1), -r2k )
  end do ! i3
  call redistributeMeshData( aDistr(2), aDat, myDistr, &
                             dstData=myDat, task=a2my(2) )
      
  ! FFT in third axis
  call myMeshBox( nMesh, aDistr(3), aBox )
  aMesh(:) = aBox(2,:) - aBox(1,:) + 1
  call redistributeMeshData( myDistr, myDat, aDistr(3), &
                             dstData=aDat, task=my2a(3) )
  call gpfa( aDat(:,:,:,1), aDat(:,:,:,2), trigs(:,3), &
             aMesh(1)*aMesh(2), 1, aMesh(3), aMesh(1)*aMesh(2), -r2k )
  call redistributeMeshData( aDistr(3), aDat, myDistr, &
                             dstData=myDat, task=a2my(3) )
      
  ! Copy data to output distr. (unless myDat and dat are the same array)
  if (.not.associated(myDat,dat)) &
    call copyMeshData( nMesh, myDistr, myDat, meshBox, dat, my2io )

  ! Divide by number of points
  if (r2k>0) dat = dat / product(nMesh)

  ! Deallocate internal arrays (possibly allocated within redistributeMeshData)
  if (associated(aDat,myDat)) then
    nullify(aDat)
  else
    call de_alloc( aDat,  name='redistributeMeshData dstData' )
  end if
  if (associated(myDat,dat)) then
    nullify(myDat)
  else
    call de_alloc( myDat, name='redistributeMeshData dstData' )
  end if
        
  ! Stop time counter
  call timer_stop( myName )

end subroutine fft3d

END MODULE m_fft3d

