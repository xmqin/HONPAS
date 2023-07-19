! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine cdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C real symmetric matrices.
C Written by G.Fabricius and J.Soler, March 1998
C Rewritten by Julian Gale, August 2004
C ************************** INPUT ******************************************
C complex*16 H(nml,nm)             : Symmetric H matrix
C complex*16 S(nml,nm)             : Symmetric S matrix
C integer n                        : Order of the generalized  system
C integer nm                       : Right hand dimension of H and S matrices
C integer nml                      : Left hand dimension of H and S matrices
C                                    which is greater than or equal to nm
C integer neigvec                  : No. of eigenvectors to calculate
C integer iscf                     : SCF cycle
C ************************** OUTPUT *****************************************
C real*8 w(nml)                    : Eigenvalues
C complex*16 Z(nml,nm)             : Eigenvectors
C integer ierror                   : Flag indicating success code for routine
C                                  :  0 = success
C                                  : -1 = repeat call as memory is increased
C                                  :  1 = fatal error
C ************************* PARALLEL ****************************************
C When running in parallel this routine now uses Scalapack to perform a
C parallel matrix diagonalisation. This requires Scalapack and Blacs to
C be installed first. Although globally a 1-D block cyclic data distribution
C is employed, locally 1 or 2-D distributions are allowed for. 
C The routine allows access to all the phases of diagonalisation for fuller
C control, and allows for parallel divide and conquer with reduced memory.
C The presence of eigenvalue clusters is checked for and the memory adjusted
C accordingly to try to guarantee convergence.
C Note that the blacs grid is only initialised on the first call to avoid
C exceeding MPI limits for split/group creation.
C ***************************************************************************
C
C  Modules
C
      use precision
      use fdf
      use parallel,   only : Node, Nodes, ProcessorY
      use parallel,   only : ParallelOverK, ResetFirstCall
!
      use parallel,   only : BlockSize_orbs => BlockSize
      use m_spin,     only : Mcoll              ! Mixed spins?
!
      use diagmemory, only : MemoryFactor
#ifdef MPI
      use mpi_siesta, only : mpi_bcast, mpi_logical, mpi_comm_world
#endif
      use alloc
      use sys,        only : die

      implicit          none

C Passed variables
      integer                 :: ierror
      integer                 :: iscf
      integer                 :: n
      integer                 :: neigvec
      integer                 :: nm
      integer                 :: nml
      real(dp)                :: w(nml)
      real(dp)                :: H(2,nml,nm)
      real(dp)                :: S(2,nml,nm)
      real(dp)                :: Z(2,nml,nm)

C Local variables
      type(allocDefaults) oldDefaults
#ifdef MPI
      integer                 :: BlockSize
      integer                 :: anb
      integer                 :: clustersize
      integer                 :: desch(9)
      integer                 :: MPIerror
      integer                 :: iacol
      integer                 :: iarow
      integer                 :: iceil
      integer,           save :: ictxt
      integer                 :: indxg2p
      integer                 :: lwork2
      integer,           save :: maxclustersize = 12
      integer                 :: mq0
      integer                 :: mycol
      integer                 :: myrow
      integer                 :: nc
      integer                 :: nn
      integer                 :: np
      integer                 :: np0
      integer                 :: npcol
      integer                 :: nprow
      integer                 :: nps
      integer                 :: nq
      integer                 :: nhegst_lwopt
      integer                 :: nhetrd_lwopt
      integer                 :: numroc
      integer                 :: nz
      integer                 :: oldmaxclustersize
      integer                 :: pjlaenv
      integer                 :: sqnpc
      integer,  pointer, save :: iclustr(:)
      real(dp), pointer, save :: gap(:)
      real(dp)                :: dscale
      logical                 :: BlacsOK
      logical                 :: lBtmp(6)
C Additional variables for a 2D proc grid
      integer                 :: np2d_row, np2d_col
      integer                 :: my2d_row, my2d_col
      integer,           save :: i2d_ctxt
      integer, dimension(9)   :: desc_h2d
      
C Matrices for 2D
      real(dp), pointer, save, dimension(:,:,:) ::  h2d, s2d, z2d
      integer ::  mat_2d_row, mat_2d_col
#endif
      character               :: jobz
      character               :: range
      integer                 :: ibtype
      integer                 :: ilaenv
      integer                 :: info
      integer                 :: liwork
      integer                 :: lrwork
      integer                 :: lwork
      integer                 :: nb
      integer                 :: neigok
      integer,  pointer, save :: ifail(:)
      integer,  pointer, save :: iwork(:)
      logical,           save :: AllInOne = .false.
      logical,           save :: DivideConquer = .true.
      logical,           save :: DoRead = .true.
      logical,           save :: FirstCall = .true.
      logical,           save :: NoExpert = .false.
      logical,           save :: PreRotate = .false.
      logical,           save :: SaveEigenvectors = .false.
      logical,           save :: Use2D = .true.
      logical,           save :: Serial = .true.
      logical                 :: ChangeAlgorithm
      real(dp)                :: abstol
      real(dp)                :: orfac
      real(dp)                :: vl
      real(dp)                :: vu
      real(dp), pointer, save :: rwork(:)
      real(dp), pointer, save :: work(:)
      real(dp), pointer, save :: Zsave(:,:,:)

      complex(dp), parameter  :: 
     .             zero_z = (0.0_dp, 0.0_dp)
     .           , one_z  = (1.0_dp, 0.0_dp)
C Start time count
      call timer('cdiag',1)

#ifdef MPI
!
! Adjust blocksize if we are dealing with mixed spins
!
      BlockSize = BlockSize_orbs * Mcoll
#endif

C Check whether first call needs to be reset due to a change of algorithm
      if (ResetFirstCall) then
        FirstCall = .true.
        ResetFirstCall = .false.
      endif

C Nullify pointers
      if (FirstCall) then
        nullify( ifail, iwork, rwork, work, Zsave )
#ifdef MPI
        nullify( gap, iclustr )
        nullify( h2d, s2d, z2d )
#endif
      endif

C Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults, copy=.false., shrink=.true.,
     .                    imin=1, routine='cdiag' )

C*******************************************************************************
C Setup                                                                        *
C*******************************************************************************

C Initialise error flag
      ierror = 0

! vl and vu are not currently used, but they must be initialized
      vl = 0
      vu = n

C Set algorithm logicals
      if (DoRead) then
        DoRead = .false.
        if (Node.eq.0) then
          AllInOne      = fdf_boolean('Diag.AllInOne',.false.)
          DivideConquer = fdf_boolean('Diag.DivideAndConquer',.true.)
          NoExpert      = fdf_boolean('Diag.NoExpert',.false.)
          PreRotate     = fdf_boolean('Diag.PreRotate',.false.)
          Use2D         = fdf_boolean('Diag.Use2D',.true.)
        endif
#ifdef MPI
        lBtmp(1) = AllInOne
        lBtmp(2) = DivideConquer
        lBtmp(3) = NoExpert
        lBtmp(4) = PreRotate
        lBtmp(5) = Use2D
        call MPI_Bcast(lBtmp,5,MPI_logical,0,MPI_Comm_World,MPIerror)
        AllInOne = lBtmp(1)
        DivideConquer = lBtmp(2) 
        NoExpert = lBtmp(3) 
        PreRotate = lBtmp(4) 
        Use2D = lBtmp(5)
#endif
      endif
      Serial = (Nodes.eq.1.or.ParallelOverK)
      if (Serial) Use2D = .false.
      if (DivideConquer.and..not.Serial) AllInOne = .false.
      if (AllInOne.and.Serial) DivideConquer = .false.
      if (PreRotate) SaveEigenvectors = .true.
      if (SaveEigenvectors) then
        if (Use2D) then
          call die('Zsave not dimensioned in cdiag')
        else
          call re_alloc( Zsave, 1,2, 1,n, 1,nm, name='Zsave', 
     .                   routine='cdiag' )
        endif
      endif

C Is it time to switch from standard diagonaliser?
      ChangeAlgorithm = (iscf.gt.2)

#ifdef MPI
      if (.not.Serial) then
C Get Blacs context and initialise Blacs grid for main grid
        nprow = 1
        npcol = Nodes
        if (FirstCall) then
          call blacs_get( -1, 0, ictxt )
          call blacs_gridinit( ictxt, 'C', nprow, npcol )
        endif

        if (Use2D) then
C If 2D approach is being used then setup secondary grid
          np2d_row = processorY
          np2d_col = Nodes/processorY
          if (FirstCall) then
            call blacs_get(ictxt, 10, i2d_ctxt)
            call blacs_gridinit(i2d_ctxt, 'R', np2d_row, np2d_col)
          endif
          call blacs_gridinfo(i2d_ctxt, np2d_row, np2d_col,
     .       my2d_row, my2d_col)
        endif

C Set up blacs descriptors for parallel case
        BlacsOK = .true.
        call descinit( desch, n, n, BlockSize, BlockSize, 0, 0, 
     .                 ictxt, n, info )
        if (info.ne.0) BlacsOK = .false.
        if (.not.BlacsOK) then
          call die('ERROR : Blacs setup has failed in cdiag!')
        endif

        if (Use2D) then
C Enquire size of local part of 2D matices
          mat_2d_row = numroc(n, BlockSize, my2d_row, 0, np2d_row)
          mat_2d_col = numroc(n, BlockSize, my2d_col, 0, np2d_col)

C Set up blacs descriptors for 2D case
          call descinit(desc_h2d, n, n, Blocksize, BlockSize, 0, 0,
     .                  i2d_ctxt,  mat_2d_row, info)
          if (info.ne.0) BlacsOK = .false.
          if (.not.BlacsOK) then
            call die('ERROR : Blacs setup has failed in cdiag!')
          endif
        endif

      endif
#endif

C Set general Lapack/Scalapack parameters
C NB: Divide and conquer doesnot have eigenvalue only case in Scalapack
      if (neigvec.gt.0.or.(DivideConquer.and..not.Serial)) then
        jobz   = 'V'
        if (neigvec.eq.n) then
          range  = 'A'
        else
          range  = 'I'
        endif
      else
        jobz   = 'N'
        range  = 'A'
      endif
      ibtype = 1
      abstol = 1.0d-8
      orfac  = 1.0d-3

C Calculate memory requirements
      if (Serial) then
        if (DivideConquer) then
          if (neigvec.gt.0) then
            lwork = 1 + 6*n 
            lrwork = 1 + 5*n + 2*n*n
            liwork = 3 + 5*n
          else
            lwork = 1 + 2*n
            lrwork = n
            liwork = 1
          endif
        else
          nb = max(ilaenv(1,'ZHETRD','U',n,-1,-1,-1),
     .             ilaenv(1,'ZUNMTR','U',n,-1,-1,-1))
          if (NoExpert) then
            lwork = (nb + 1)*n
            lrwork = max(1,3*n-2)
            liwork = 1
          else
            lwork = (nb + 1)*n
            lrwork = 7*n
            liwork = 5*n
          endif
        endif
#ifdef MPI
      else
        if (Use2D) then
          call blacs_gridinfo(desc_h2d(2),nprow,npcol,myrow,mycol)
          iarow = indxg2p(1,desc_h2d(6),myrow,desc_h2d(7),nprow)
          iacol = indxg2p(1,desc_h2d(5),mycol,desc_h2d(8),npcol)
        else
          call blacs_gridinfo(desch(2),nprow,npcol,myrow,mycol)
          iarow = indxg2p(1,desch(6),myrow,desch(7),nprow)
          iacol = indxg2p(1,desch(5),mycol,desch(8),npcol)
        endif
        nn = max(n,BlockSize)
        nn = max(nn,2)
        np0 = numroc(nn,BlockSize,0,0,nprow)
        mq0 = numroc(nn,BlockSize,0,0,npcol)
        np = numroc(n,BlockSize,myrow,iarow,nprow)
        nq = numroc(n,BlockSize,mycol,iacol,npcol)
        if (DivideConquer) then
          lwork = n + BlockSize*(np0 + mq0 + BlockSize)
          nhegst_lwopt = (2*np + nq)*BlockSize + 
     .                   BlockSize*BlockSize
          lwork = max(lwork,nhegst_lwopt)
C Note: Scalapack documentation fibs about the size of lrwork!
          lrwork = 1 + 10*n + 3*np*nq
          liwork = 7*n + 8*Nodes + 2
        else
          if (neigvec.gt.0) then
            lwork = n + (np0 + mq0 + BlockSize)*BlockSize

            anb = pjlaenv(ictxt,3,'PDSYTTRD','L',0,0,0,0)
            sqnpc = int(sqrt(dble(Nodes)))
            nps = max(numroc(n,1,0,0,sqnpc),2*anb)
            nhetrd_lwopt = 2*(anb + 1)*(4*nps + 2) + 
     .                     (nps + 1)*nps
            nhegst_lwopt = (2*np + nq)*BlockSize + 
     .                     BlockSize*BlockSize
            lwork2 = max(n + nhetrd_lwopt,nhegst_lwopt)
            lwork = max(lwork,lwork2)
            lwork2 = n + max(nn*(np0+1),3)
            lwork = max(lwork,lwork2)
            lrwork = 4*n + max(5*nn,np0*mq0) + iceil(neigvec,Nodes)*nn +
     .               maxclustersize*n
          else
            lwork = n + max(nn*(np0+1),3)
            nhegst_lwopt = (2*np + nq)*BlockSize + 
     .                     BlockSize*BlockSize
            lwork = max(lwork,nhegst_lwopt)
            lrwork = 5*nn + 4*n
          endif
          liwork = 6*max(n,Nodes+1,4)
        endif
#endif
      endif

C Scale memory by memory factor
      lwork = nint(MemoryFactor*dble(lwork))

C Allocate workspace arrays
      call re_alloc( work,    1,2*lwork, name='work', 
     .	routine='cdiag' )
      call re_alloc( rwork,   1,lrwork,  name='rwork',
     .  routine='cdiag')
      call re_alloc( iwork,   1,liwork,  name='iwork',
     .  routine='cdiag')
      call re_alloc( ifail,   1,n,       name='ifail',
     .  routine='cdiag')
#ifdef MPI
      call re_alloc( gap,     1,Nodes,   name='gap',
     .  routine='cdiag')
      call re_alloc( iclustr, 1,2*Nodes, name='iclustr',
     .  routine='cdiag')

      if (Use2D) then
C Set up workspace arrays for 2D versions of 1D arrays
        call re_alloc(h2d, 1,2, 1,mat_2d_row, 1,mat_2d_col, name='h2d',
     .    routine='cdiag')
        call re_alloc(s2d, 1,2, 1,mat_2d_row, 1,mat_2d_col, name='s2d',
     .    routine='cdiag')
        call re_alloc(z2d, 1,2, 1,mat_2d_row, 1,mat_2d_col, name='z2d',
     .    routine='cdiag')

C For 2D case, copy arrays to new distribution
        call pzgemr2d(n, n, H, 1, 1, desch, h2d, 1, 1, desc_h2d, ictxt)
        call pzgemr2d(n, n, S, 1, 1, desch, s2d, 1, 1, desc_h2d, ictxt)
      endif
#endif

      if (AllInOne) then
C*******************************************************************************
C Perform diagonalisation using simple driver                                  *
C*******************************************************************************
        if (Serial) then
          if (NoExpert) then
            call zhegv(1,jobz,'U',n,H,n,S,n,w,work,lwork,rwork,info)
            neigok = n
            call zcopy(n*n,H,1,Z,1)
          else
            call zhegvx(1,jobz,range,'U',n,H,n,S,n,vl,vu,1,neigvec,
     .                  abstol,neigok,w,Z,n,work,lwork,rwork,iwork,
     .                  ifail,info)
          endif
#ifdef MPI
        else
          if (Use2D) then
            call pzhegvx(1,jobz,range,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .                   desc_h2d,vl,vu,1,neigvec,abstol,neigok,nz,w,
     .                   orfac,z2d,1,1,desc_h2d,work,lwork,rwork,lrwork,
     .                   iwork,liwork,ifail,iclustr,gap,info)
            call pzgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
          else
            call pzhegvx(1,jobz,range,'U',n,H,1,1,desch,S,1,1,desch,
     .                   vl,vu,1,neigvec,abstol,neigok,nz,w,orfac,Z,
     .                   1,1,desch,work,lwork,rwork,lrwork,iwork,liwork,
     .                   ifail,iclustr,gap,info)
          endif

C If there is insufficient memory due to eigenvalue clustering increase memory
          if (info.eq.2) then
            oldmaxclustersize = maxclustersize 
            do nc = 1,Nodes 
              clustersize = iclustr(2*nc) - iclustr(2*nc-1)
              maxclustersize = max(maxclustersize,clustersize)
            enddo
            if (maxclustersize.gt.oldmaxclustersize) then
C If the memory has increased then try again - otherwise fail
              ierror = - 1
              goto 999
            endif
          endif
#endif
        endif

C Check error flag
        if (info.ne.0) then
          ierror = 1
          if (info.lt.0) then
            call die('Illegal argument to general eigensolver')
          elseif (info.gt.0)   then
            if (mod(info/2,2).ne.0) then
              if (Node.eq.0) then
                write(6,'(/,''Clustered eigenvectors not converged - '',
     .            ''more memory required'',/)')
              endif
            endif
            call die('Failure to converge general eigenproblem')
          endif
        endif
        if (neigok.lt.neigvec) then
          call die('Insufficient eigenvalues converged in cdiag')
        endif

      else
C*******************************************************************************
C Factorise overlap matrix                                                     *
C*******************************************************************************
        call timer('cdiag1',1)
        if (Serial) then
          call zpotrf('U',n,S,n,info)
#ifdef MPI
        else
          if (Use2D) then
            call pzpotrf('U',n,s2d,1,1,desc_h2d,info)
          else
            call pzpotrf('U',n,S,1,1,desch,info)
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in Cholesky factorisation in cdiag')
        endif
        call timer('cdiag1',2)

C*******************************************************************************
C Transform problem to standard eigenvalue problem                             *
C*******************************************************************************
        call timer('cdiag2',1)
        if (Serial) then
          call zhegst(ibtype,'U',n,H,n,S,n,info)
#ifdef MPI
        else
          if (Use2D) then
            call pzhengst(ibtype,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .                    desc_h2d,dscale,work,lwork,info)
          else
            call pzhengst(ibtype,'U',n,H,1,1,desch,S,1,1,
     .                    desch,dscale,work,lwork,info)
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in forward transformation in cdiag')
        endif
        call timer('cdiag2',2)

C*******************************************************************************
C Pre-rotate using trial eigenvectors                                          *
C*******************************************************************************
        if (PreRotate.and.ChangeAlgorithm) then
          if (Serial) then
C Rotate according to eigenvectors from previous cycle
            call zhemm('L','U',n,n,one_z,H,n,Zsave,n,zero_z,Z,n)
            call zgemm('T','N',n,n,n,one_z,Zsave,n,Z,n,zero_z,H,n)
#ifdef MPI
          else
            if (Use2D) then
              call pzhemm('L','U',n,n,one_z,h2d,1,1,desc_h2d,Zsave,1,1,
     .                     desc_h2d,zero_z,z2d,1,1,desc_h2d)
              call pzgemm('T','N',n,n,n,one_z,Zsave,1,1,desc_h2d,z2d,1,
     .                     1,desch,zero_z,h2d,1,1,desc_h2d)
            else
              call pzhemm('L','U',n,n,one_z,H,1,1,desch,Zsave,1,1,desch,
     .                     zero_z,Z,1,1,desch)
              call pzgemm('T','N',n,n,n,one_z,Zsave,1,1,desch,Z,1,1,
     .                     desch,zero_z,H,1,1,desch)
            endif
#endif
          endif
        endif

C*******************************************************************************
C Solve standard eigenvalue problem                                            *
C*******************************************************************************
        call timer('cdiag3',1)
        if (Serial) then
          if (DivideConquer) then
            call zheevds(jobz,'U',n,H,n,Z,n,w,work,lwork,rwork,lrwork,
     .                   iwork,liwork,info)
            neigok = n
          else
            call zheevx(jobz,range,'U',n,H,n,vl,vu,1,neigvec,abstol,
     .                  neigok,w,Z,n,work,lwork,rwork,iwork,ifail,info)
          endif
#ifdef MPI
        else
          if (DivideConquer) then
            if (Use2D) then
              call pzheevd(jobz,'U',n,h2d,1,1,desc_h2d,w,z2d,1,1,
     .                      desc_h2d,work,lwork,rwork,lrwork,iwork,
     .                      liwork,info)
            else
              call pzheevd(jobz,'U',n,H,1,1,desch,w,Z,1,1,desch,work,
     .                      lwork,rwork,lrwork,iwork,liwork,info)
            endif
            neigok = n
          else
            if (Use2D) then
              call pzheevx(jobz,range,'U',n,h2d,1,1,desc_h2d,vl,vu,1,
     .                     neigvec,abstol,neigok,nz,w,orfac,z2d,1,1,
     .                     desc_h2d,work,lwork,rwork,lrwork,iwork,
     .                     liwork,ifail,iclustr,gap,info)
            else
              call pzheevx(jobz,range,'U',n,H,1,1,desch,vl,vu,1,neigvec,
     .                     abstol,neigok,nz,w,orfac,Z,1,1,desch,work,
     .                     lwork,rwork,lrwork,iwork,liwork,ifail,
     .                     iclustr,gap,info)
            endif

C If there is insufficient memory due to eigenvalue clustering increase memory
            if (info.eq.2) then
              oldmaxclustersize = maxclustersize
              do nc = 1,Nodes
                clustersize = iclustr(2*nc) - iclustr(2*nc-1)
                maxclustersize = max(maxclustersize,clustersize)
              enddo
              if (maxclustersize.gt.oldmaxclustersize) then
C If the memory has increased then try again - otherwise fail
                ierror = - 1
                goto 999
              endif
            endif

          endif
#endif
        endif

C Check error flag
        if (info.ne.0) then
          ierror = 1
          if (info.lt.0) then
            call die('Illegal argument to standard eigensolver')
          elseif (info.gt.0)   then
            if (mod(info/2,2).ne.0) then
              if (Node.eq.0) then
                write(6,'(/,''Clustered eigenvectors not converged - '',
     .            ''more memory required'',/)')
              endif
            endif
            call die('Failure to converge standard eigenproblem')
          endif
        endif
        if (neigok.lt.neigvec) then
          call die('Insufficient eigenvalues converged in cdiag')
        endif

        call timer('cdiag3',2)

C*******************************************************************************
C Post-rotate using trial eigenvectors                                         *
C*******************************************************************************
        if (PreRotate.and.ChangeAlgorithm) then
          if (Serial) then
            call zgemm('N','N',n,n,n,one_z,Zsave,n,Z,n,zero_z,H,n)
            call zcopy(n*nm,H(1,1,1),1,Z(1,1,1),1)
#ifdef MPI
          else
            if (Use2D) then
              call pzgemm('N','N',n,n,n,one_z,Zsave,1,1,desc_h2d,z2d,1,
     .                     1,desc_h2d,zero_z,h2d,1,1,desc_h2d)
              call zcopy(n*nm,h2d(1,1,1),1,z2d(1,1,1),1)
            else
              call pzgemm('N','N',n,n,n,one_z,Zsave,1,1,desch,Z,1,1,
     .                     desch,zero_z,H,1,1,desch)
              call zcopy(n*nm,H(1,1,1),1,Z(1,1,1),1)
            endif
#endif
          endif
        endif
  
C Save eigenvectors if required
        if (SaveEigenvectors) then
#ifdef MPI
          if (Use2D) then
            call zcopy(n*nm,z2d(1,1,1),1,Zsave(1,1,1),1)
          else
#endif
            call zcopy(n*nm,Z(1,1,1),1,Zsave(1,1,1),1)
#ifdef MPI
          endif
#endif
        endif

C*******************************************************************************
C Back transformation                                                          *
C*******************************************************************************
        call timer('cdiag4',1)
        if (Serial) then
          if (neigvec.gt.0) then
            call ztrsm('Left', 'U', 'N', 'Non-unit', 
     .                 n, neigvec, one_z, S, n, Z, n)
          endif
#ifdef MPI
        else
          if (neigvec.gt.0) then
            if (Use2D) then
              call pztrsm('Left','U','N','Non-unit',n,neigvec,
     .                     one_z,s2d,1,1,desc_h2d,z2d,1,1,desc_h2d)
              call pzgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
            else
              call pztrsm('Left','U','N','Non-unit',n,neigvec,
     .                     one_z,S,1,1,desch,Z,1,1,desch)
            endif
          endif
          if (dscale.ne.1.0_dp) then
            call dscal(n,dscale,w,1)
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in back transformation in cdiag')
        endif
        call timer('cdiag4',2)

      endif

C*******************************************************************************
C Clean up                                                                     *
C*******************************************************************************

C Common exit point 
  999 continue

C Set first time call logical to false
      FirstCall = .false.

C Deallocate workspace arrays
#ifdef MPI
      call de_alloc( iclustr, name='iclustr')
      call de_alloc( gap,     name='gap'    )
#endif
      call de_alloc( ifail,   name='ifail'  )
      call de_alloc( iwork,   name='iwork'  )
      call de_alloc( work,    name='work'   )

C  Restore old allocation defaults
      call alloc_default( restore=oldDefaults )

C Stop time count
      call timer('cdiag',2)

      end
