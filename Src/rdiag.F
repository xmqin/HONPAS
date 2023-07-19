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
      subroutine rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)
C ***************************************************************************
C Subroutine  to solve all eigenvalues and eigenvectors of the
C real general eigenvalue problem  H z = w S z,  with H and S
C real symmetric matrices.
C Written by G.Fabricius and J.Soler, March 1998
C Rewritten by Julian Gale, August 2004
C ************************** INPUT ******************************************
C real*8 H(nml,nm)                 : Symmetric H matrix
C real*8 S(nml,nm)                 : Symmetric S matrix
C integer n                        : Order of the generalized  system
C integer nm                       : Right hand dimension of H and S matrices
C integer nml                      : Left hand dimension of H and S matrices
C                                    which is greater than or equal to nm
C integer neigvec                  : No. of eigenvectors to calculate
C integer iscf                     : SCF cycle
C ************************** OUTPUT *****************************************
C real*8 w(nml)                    : Eigenvalues
C real*8 Z(nml,nm)                 : Eigenvectors
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
      use parallel,   only : Node, Nodes, BlockSize, ProcessorY
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
      real(dp)                :: H(nml,nm)
      real(dp)                :: S(nml,nm)
      real(dp)                :: w(nml)
      real(dp)                :: Z(nml,nm)

C Local variables
      type(allocDefaults) oldDefaults
#ifdef MPI
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
      integer                 :: nz
      integer                 :: nn
      integer                 :: np
      integer                 :: np0
      integer                 :: npcol
      integer                 :: nprow
      integer                 :: nps
      integer                 :: nq
      integer                 :: nsygst_lwopt
      integer                 :: nsytrd_lwopt
      integer                 :: nt
      integer                 :: numroc
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
      real(dp), pointer, save, dimension(:,:) ::  h2d, s2d, z2d
      integer ::  mat_2d_row, mat_2d_col
#endif
      character               :: jobz
      character               :: range
      integer                 :: ilaenv
      integer                 :: info
      integer                 :: liwork
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
      real(dp)                :: vl
      real(dp)                :: vu
      real(dp), pointer, save :: work(:)
      real(dp), pointer, save :: Zsave(:,:)

      real(dp), parameter     :: abstol = 1.0e-8_dp
      real(dp), parameter     :: orfac = 1.0e-3_dp
      real(dp), parameter     :: zero = 0.0_dp
      real(dp), parameter     :: one = 1.0_dp

C Start time count
      call timer('rdiag',1)

C Nullify pointers
      if (FirstCall) then
        nullify( ifail, iwork, work, Zsave )
#ifdef MPI
        nullify( gap, iclustr )
        nullify( h2d, s2d, z2d )
#endif
      endif

C Get old allocation defaults and set new ones
      call alloc_default( old=oldDefaults, copy=.false., shrink=.true.,
     .                    imin=1, routine='rdiag' )

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
        if (Node.eq.0) then
           lBtmp(1) = AllInOne
           lBtmp(2) = DivideConquer
           lBtmp(3) = NoExpert
           lBtmp(4) = PreRotate
           lBtmp(5) = Use2D
        endif
        call MPI_Bcast(lBtmp,5,MPI_logical,0,MPI_Comm_World,MPIerror)
        AllInOne = lBtmp(1)
        DivideConquer = lBtmp(2)
        NoExpert = lBtmp(3)
        PreRotate = lBtmp(4)
        Use2D = lBtmp(5)
#endif
        Serial = (Nodes.eq.1)
        if (Serial) Use2D = .false.
        if (AllInOne.and.Serial) DivideConquer = .false.
        if (neigvec.ne.n) PreRotate = .false.
        if (PreRotate) SaveEigenvectors = .true.
        if (Nodes.eq.1) Use2D = .false.
        if (SaveEigenvectors) then
          if (Use2D) then
            call die('Zsave not dimensioned in rdiag')
          else
            call re_alloc( Zsave, 1, nml, 1, nm, name='Zsave', 
     .                     routine='rdiag' )
          endif
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
          call die('ERROR : Blacs setup has failed in rdiag!')
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
            call die('ERROR : Blacs setup has failed in rdiag!')
          endif
        endif

      endif
#endif

C Set general Lapack parameters
      if (neigvec.gt.0) then
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

C Calculate memory requirements
      if (Serial) then
        if (DivideConquer) then
          if (neigvec.gt.0) then
            lwork = 1 + 6*n + n*n
            liwork = 3 + 5*n
          else
            lwork = 1 + 2*n
            liwork = 1
          endif
        else
          nb = ilaenv(1,'DSYTRD','U',n,-1,-1,-1)
          if (NoExpert) then
            lwork = (nb+2)*n
            liwork = 1
          else
            lwork = max(8*n,(nb+3)*n)
            liwork = 5*n
          endif
        endif
#ifdef MPI
      else
        if (Use2D) then
          call blacs_gridinfo(desc_h2d(2),nprow,npcol,myrow,mycol)
        else
          call blacs_gridinfo(desch(2),nprow,npcol,myrow,mycol)
        endif
        if (DivideConquer) then
          np0 = numroc(n,BlockSize,0,0,nprow)
          mq0 = numroc(n,BlockSize,0,0,npcol)
          lwork = max(BlockSize*(np0+1),3*BlockSize)
          lwork2 = BlockSize*(2*np0 + mq0 + BlockSize)
          lwork = max(lwork,lwork2)
          iarow = indxg2p(1,BlockSize,myrow,desch(7),nprow)
          iacol = indxg2p(1,BlockSize,mycol,desch(8),npcol)
          np = numroc(n,BlockSize,myrow,iarow,nprow)
          nq = numroc(n,BlockSize,mycol,iacol,npcol)
          nt = 3*n + max(BlockSize*(np+1),3*BlockSize)
          lwork2 = max(1+6*n+2*np*nq,nt) + 2*n
          lwork = max(lwork,lwork2)
          liwork = 7*n + 8*Nodes + 2
        else
          nn = max(n,BlockSize)
          nn = max(nn,2)
          np0 = numroc(nn,BlockSize,0,0,nprow)
          mq0 = numroc(nn,BlockSize,0,0,npcol)
          if (neigvec.gt.0) then
            lwork = 5*n + max(5*nn,np0*mq0+2*BlockSize*BlockSize) +
     .              iceil(neigvec,Nodes)*nn + maxclustersize*n
            anb = pjlaenv(ictxt,3,'PDSYTTRD','L',0,0,0,0)
            sqnpc = int(sqrt(dble(Nodes)))
            nps = max(numroc(n,1,0,0,sqnpc),2*anb)
            nsytrd_lwopt = n + 2*(anb + 1)*(4*nps + 2) + (nps + 3)*nps
            nsygst_lwopt = (2*np0 + mq0)*BlockSize + BlockSize*BlockSize
            lwork2 = max(5*n + nsytrd_lwopt,nsygst_lwopt)
            lwork = max(lwork,lwork2)
          else
            lwork = 5*n + max(5*nn,BlockSize*(np0+1))
          endif
          liwork = 6*max(n,Nodes+1,4)
        endif
#endif
      endif

C Scale memory by memory factor
      lwork = nint(MemoryFactor*dble(lwork))

C Allocate workspace arrays
      call re_alloc( work,    1,lwork,   name='work'  )
      call re_alloc( iwork,   1,liwork,  name='iwork' )
      call re_alloc( ifail,   1,n,       name='ifail'   )
#ifdef MPI
      call re_alloc( gap,     1,Nodes,   name='gap'     )
      call re_alloc( iclustr, 1,2*Nodes, name='iclustr' )

      if (Use2D) then
C Set up workspace arrays for 2D versions of 1D arrays
        call re_alloc( h2d, 1, mat_2d_row, 1, mat_2d_col, name='h2d')
        call re_alloc( s2d, 1, mat_2d_row, 1, mat_2d_col, name='s2d')
        call re_alloc( z2d, 1, mat_2d_row, 1, mat_2d_col, name='z2d')

C For 2D case, copy arrays to new distribution
        call pdgemr2d(n, n, H, 1, 1, desch, h2d, 1, 1, desc_h2d, ictxt)
        call pdgemr2d(n, n, S, 1, 1, desch, s2d, 1, 1, desc_h2d, ictxt)
      endif
#endif

      if (AllInOne) then
C*******************************************************************************
C Perform diagonalisation using simple driver                                  *
C*******************************************************************************
        if (Serial) then
          if (NoExpert) then
            call dsygv(1,jobz,'U',n,H,n,S,n,w,work,lwork,info)
            neigok = n
            call dcopy(n*n,H,1,Z,1)
          else
            call dsygvx(1,jobz,range,'U',n,H,n,S,n,vl,vu,1,neigvec,
     .                  abstol,neigok,w,Z,n,work,lwork,iwork,ifail,
     .                  info)
          endif
#ifdef MPI
        else 
   10     continue
          if (Use2D) then
            call pdsygvx(1,jobz,range,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .                   desc_h2d,vl,vu,1,neigvec,abstol,neigok,nz,w,
     .                   orfac,z2d,1,1,desc_h2d,work,lwork,iwork,liwork,
     .                   ifail,iclustr,gap,info)
            call pdgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
          else
            call pdsygvx(1,jobz,range,'U',n,H,1,1,desch,S,1,1,desch,
     .                   vl,vu,1,neigvec,abstol,neigok,nz,w,orfac,Z,
     .                   1,1,desch,work,lwork,iwork,liwork,ifail,
     .                   iclustr,gap,info)
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
          call die('Insufficient eigenvalues converged in rdiag')
        endif

      else
C*******************************************************************************
C Factorise overlap matrix                                                     *
C*******************************************************************************
        call timer('rdiag1',1)
        if (Serial) then
          call dpotrf('U',n,S,n,info)
#ifdef MPI
        else
          if (Use2D) then
            call pdpotrf('U',n,s2d,1,1,desc_h2d,info)
          else
            call pdpotrf('U',n,S,1,1,desch,info)
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in Cholesky factorisation in rdiag')
        endif
        call timer('rdiag1',2)

C*******************************************************************************
C Transform problem to standard eigenvalue problem                             *
C*******************************************************************************
        call timer('rdiag2',1)
        if (Serial) then
          call dsygst(1,'U',n,H,n,S,n,info)
#ifdef MPI
        else
          if (Use2D) then
            call pdsyngst(1,'U',n,h2d,1,1,desc_h2d,s2d,1,1,
     .                    desc_h2d,dscale,work,lwork,info)
          else
            call pdsyngst(1,'U',n,H,1,1,desch,S,1,1,
     .                    desch,dscale,work,lwork,info)
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in forward transformation in rdiag')
        endif
        call timer('rdiag2',2)

C*******************************************************************************
C Pre-rotate using trial eigenvectors                                          *
C*******************************************************************************
        if (PreRotate.and.ChangeAlgorithm) then
          if (Serial) then
C Rotate according to eigenvectors from previous cycle
            call dsymm('L','U',n,n,one,H,n,Zsave,n,zero,Z,n)
            call dgemm('T','N',n,n,n,one,Zsave,n,Z,n,zero,H,n)
#ifdef MPI
          else
            if (Use2D) then
              call pdsymm('L','U',n,n,one,h2d,1,1,desc_h2d,Zsave,1,1,
     .                     desc_h2d,zero,z2d,1,1,desc_h2d)
              call pdgemm('T','N',n,n,n,one,Zsave,1,1,desc_h2d,z2d,1,
     .                     1,desch,zero,h2d,1,1,desc_h2d)
            else
              call pdsymm('L','U',n,n,one,H,1,1,desch,Zsave,1,1,desch,
     .                     zero,Z,1,1,desch)
              call pdgemm('T','N',n,n,n,one,Zsave,1,1,desch,Z,1,1,
     .                     desch,zero,H,1,1,desch)
            endif
#endif
          endif
        endif

C*******************************************************************************
C Solve standard eigenvalue problem                                            *
C*******************************************************************************
        call timer('rdiag3',1)
        if (Serial) then
          if (DivideConquer) then
            call dsyevds(jobz,'U',n,H,n,w,Z,n,work,lwork,iwork,liwork,
     .                   info)
            neigok = n
          else
            call dsyevx(jobz,range,'U',n,H,n,vl,vu,1,neigvec,abstol,
     .                  neigok,w,Z,n,work,lwork,iwork,ifail,info)
          endif
#ifdef MPI
        else
          if (DivideConquer) then
            if (Use2D) then
              call pdsyevd(jobz,'U',n,h2d,1,1,desc_h2d,w,z2d,1,1,
     .                     desc_h2d,work,lwork,iwork,liwork,info)
            else
              call pdsyevd(jobz,'U',n,H,1,1,desch,w,Z,1,1,desch,work,
     .                     lwork,iwork,liwork,info)
            endif
            neigok = n
          else
            if (Use2D) then
              call pdsyevx(jobz,range,'U',n,h2d,1,1,desc_h2d,vl,vu,1,
     .                     neigvec,abstol,neigok,nz,w,orfac,z2d,1,1,
     .                     desc_h2d,work,lwork,iwork,liwork,ifail,
     .                     iclustr,gap,info)
            else
              call pdsyevx(jobz,range,'U',n,H,1,1,desch,vl,vu,1,neigvec,
     .                     abstol,neigok,nz,w,orfac,Z,1,1,desch,work,
     .                     lwork,iwork,liwork,ifail,iclustr,gap,info)
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
                call timer('rdiag3',2)
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
          call die('Insufficient eigenvalues converged in rdiag')
        endif

        call timer('rdiag3',2)

C*******************************************************************************
C Post-rotate using trial eigenvectors                                         *
C*******************************************************************************
        if (PreRotate.and.ChangeAlgorithm) then
          if (Serial) then
            call dgemm('N','N',n,n,n,one,Zsave,n,Z,n,zero,H,n)
            call dcopy(n*nm,H(1,1),1,Z(1,1),1)
#ifdef MPI
          else
            if (Use2D) then
              call pdgemm('N','N',n,n,n,one,Zsave,1,1,desc_h2d,z2d,1,
     .                     1,desc_h2d,zero,h2d,1,1,desc_h2d)
              call dcopy(n*nm,h2d(1,1),1,z2d(1,1),1)
            else
              call pdgemm('N','N',n,n,n,one,Zsave,1,1,desch,Z,1,1,
     .                     desch,zero,H,1,1,desch)
              call dcopy(n*nm,H(1,1),1,Z(1,1),1)
            endif
#endif
          endif
        endif

C Save eigenvectors if required
        if (SaveEigenvectors) then
#ifdef MPI
          if (Use2D) then
            call dcopy(n*nm,z2d(1,1),1,Zsave(1,1),1)
          else
#endif
            call dcopy(n*nm,Z(1,1),1,Zsave(1,1),1)
#ifdef MPI
          endif
#endif
        endif

C*******************************************************************************
C Back transformation                                                          *
C*******************************************************************************
        call timer('rdiag4',1)
        if (Serial) then
          if (neigvec.gt.0) then
            call dtrsm('Left','U','N','Non-unit',n,neigvec,one,S,n,Z,n)
          endif
#ifdef MPI
        else
          if (neigvec.gt.0) then
            if (Use2D) then
              call pdtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .                     s2d,1,1,desc_h2d,z2d,1,1,desc_h2d)
              call pdgemr2d(n,n,z2d,1,1,desc_h2d,Z,1,1,desch,ictxt)
            else
              call pdtrsm('Left','U','N','Non-unit',n,neigvec,one,
     .                     S,1,1,desch,Z,1,1,desch)
            endif
            if (dscale.ne.one) then
              call dscal(n,dscale,w,1)
            endif
          endif
#endif
        endif
        if (info.ne.0) then
          call die('Error in back transformation in rdiag')
        endif
        call timer('rdiag4',2)

      endif

C*******************************************************************************
C Clean up                                                                     *
C*******************************************************************************

C Common exit point 
  999 continue

C Set first call flag to false
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
      call timer('rdiag',2)

      end
