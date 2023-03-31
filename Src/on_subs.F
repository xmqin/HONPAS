! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_on_subs
      public :: ctrans1, ctrans2, axb_build1, axb_build2
      CONTAINS
      subroutine ctrans1(no_cl,ntmax)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector cttoc.
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nbands         : Number of columns of full C
C integer ntmax          : maximum number of nonzero elements of each
C                           column of C
C **************************** OUTPUT *********************************
C integer numt(nc)       : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer listt(ntmax,nc) : Control vector of C transpose matrix
C                           (list of nonzero elements of each column of C)
C integer cttoc(ntmax,nc) : Map from C transpose to C indexing (pointer)
C *********************************************************************

      use alloc
      use on_core, only : numct, listct, cttoc
      use on_main, only : numc, listc, ncP2T, nbG2L, nbandsloc

      implicit none

      integer, intent(in)    :: no_cl
      integer, intent(inout) :: ntmax

C Internal variables ..................................................
      integer
     .  i,il,imu,mu,n

C  Initialize numt list 
      do il = 1,nbandsloc
        numct(il) = 0
      enddo

C  Construct information for transpose of C 
      ! numct(b_l): number of coefficients
      ! handled here (for 1:no_cl) for the
      ! local band b_l
      ! By definition, a band is local if
      ! this number is not zero
      do mu = 1,no_cl
        do imu = 1,numc(mu)
          i = listc(imu,mu)
          il = nbG2L(i)
          if (il.gt.0) then
            numct(il) = numct(il) + 1
            n = numct(il)
            ntmax = max(ntmax,n)
          endif
        enddo
      enddo

C  Re-size arrays
      call re_alloc(listct,1,ntmax,1,nbandsloc,name='listct')
      call re_alloc(cttoc,1,ntmax,1,nbandsloc,name='cttoc')

C  Reinitialize numt list 
      do il = 1,nbandsloc
        numct(il) = 0
      enddo

C  Construct information for transpose of C 
      do mu = 1,no_cl
        do imu = 1,numc(mu)
          i = listc(imu,mu)
          il = nbG2L(i)
          if (il.gt.0) then
            numct(il) = numct(il) + 1
            n = numct(il)
            listct(n,il) = mu
            cttoc(n,il) = imu
          endif
        enddo
      enddo

      end subroutine ctrans1


      subroutine ctrans2(no_cl,maxnf,maxnft,numf,listf)
C ********************************************************************
C Finds the C transpose matrix control vectors numt and listt,
C and the index vector cttoc.
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nbands         : Number of rows of C 
C integer nbasis         : Number of columns of full C
C integer maxnf          : First dimension of list and C, and maximum
C                           number of nonzero elements of each row of C
C integer maxnft         : maximum number of nonzero elements of each
C                           column of C
C integer numf(nr)        : Control vector of C matrix
C                           (number of nonzero elements of each row of C)
C integer listf(maxnf,nr) : Control vector of C matrix (pointer)
C                          (list of nonzero elements of each row of C)
C **************************** OUTPUT *********************************
C integer numft(nbasis)       : Control vector of C matrix
C                         (number of nonzero elements of each column of C)
C integer listft(maxnft,nbasis) : Control vector of C transpose matrix
C                           (list of nonzero elements of each column of C)
C integer fttof(maxnft,nbasis) : Map from C transpose to C indexing (pointer)
C *********************************************************************

      use alloc
      use on_core, only : numft, listft, fttof
      use on_main, only : ncG2L, nbL2G, nbandsloc

      implicit none

      integer, intent(in) ::
     .  no_cl,maxnf,listf(maxnf,nbandsloc),numf(nbandsloc)
      integer, intent(out) :: maxnft

C Internal variables ..................................................
      integer
     .  i,ii,imu,mu,mul,n

C  Initialize numft list
      do i = 1,no_cl
        numft(i) = 0
      enddo

C  Find maximum dimension for arrays
      do mul = 1,nbandsloc
        mu = nbL2G(mul)
        do imu = 1,numf(mul)
          i = listf(imu,mul)
          ii = ncG2L(i)
          if (ii.gt.0) then
            numft(ii) = numft(ii) + 1
          endif
        enddo
      enddo

C  Find maximum value of numft and re-zero numft
      maxnft = 0
      do i = 1,no_cl
        maxnft = max(maxnft,numft(i))
        numft(i) = 0
      enddo

C Resize arrays
      call re_alloc(listft,1,maxnft,1,no_cl,name='listft')
      call re_alloc(fttof,1,maxnft,1,no_cl,name='fttof')

C  Construct information for transpose of F
      do mul = 1,nbandsloc
        mu = nbL2G(mul)
        do imu = 1,numf(mul)
          i = listf(imu,mul)
          ii = ncG2L(i)
          if (ii.gt.0) then
            numft(ii) = numft(ii) + 1
            n = numft(ii)
            listft(n,ii) = mu
            fttof(n,ii) = imu
          endif
        enddo
      enddo

      end subroutine ctrans2


      subroutine axb_build1(nbasis,maxnct,numct,listct,nhmax,numh,
     .                      listh,listhptr,indon,maxnf)
C ********************************************************************
C Constructs control indexes of a C matrix in sparse form,
C C being the product of A and B (also in sparse form)
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nbands x nbasis
C               B is rectangular, and has dimension:  nbasis x nbasis
C and, as a result:
C               C is rectangular, and has dimension:  nbands x nbasis
C
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nbands             : Number of rows of A 
C integer nbasis             : Number of columns of A
C integer maxnct             : First dimension of A matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of A)
C integer numct(nbands)      : Control vector of A matrix
C                            (number of nonzero elements of each row of A)
C integer listct(maxnct,nbands)  : Control vector of A matrix
C                           (list of nonzero elements of each row of A)
C integer nbasis             : Number of rows of B
C integer nbasis             : Number of columns of B
C integer nhmax              : First dimension of B matrix in sparse form,
C                              as declared in calling routine
C                              (max. number of <>0 elements of each row of B)
C integer numh(nbasis)         : Control vector of B matrix
C                            (number of nonzero elements of each row of B)
C integer listh(nhmax)       : Control vector of B matrix
C                            (pointer to start of row in listb)
C integer listhptr(nbasis)   : Control vector of B matrix
C                            (list of nonzero elements of each row of B)
C integer indon(nbasis)      : Auxiliary array to build C in sparse form
C integer nindv(maxnf)       : Auxiliary array to store indexes of nonzero
C                              matrix elements of each row of C
C ************************INPUT/OUTPUT ********************************
C integer maxnf              : First dimension of C matrix in sparse form,
C                              as declared in calling routine
C                              (max. number of <>0 elements of each row of C)
C **************************** OUTPUT *********************************
C integer numc(nbands)       : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer listc(maxnf,nbands): Control vector of C matrix
C                            (list of nonzero elements of each row of C)
C *********************************************************************

      use alloc
      use on_core, only : nindv, numf, listf
      use on_main, only : ncT2P, nbandsloc, nbL2G, nbG2L

      implicit none

      integer, intent(in) ::
     .  nbasis,maxnct,nhmax,
     .  listct(maxnct,nbandsloc),listh(nhmax),
     .  listhptr(nbasis),numct(nbandsloc),numh(nbasis)
      integer, intent(inout) :: maxnf
      integer, intent(out)   :: indon(nbasis)

C Internal variables..................................................
      integer
     .  i,il,in,indk,j,k,kn,nind,nloc
      logical
     .  lbandloc

C Ensure arrays are large enough at entrance
      call re_alloc(nindv,1,max(maxnf,nbasis),name='nindv')

C Initialize internal variables .......................................
      nind = 0
      do i = 1,nbasis
        indon(i) = 0
      enddo
      do i = 1,maxnf
        nindv(i) = 0
      enddo

C Initialise number of local bands
      nloc = 0

C Find out C control vectors  - get number of elements
      do il = 1,nbandsloc
        lbandloc = .false.
        do in = 1,numct(il)
          k = ncT2P(listct(in,il))
          if (k.gt.0) then
            lbandloc = .true.
            indk = listhptr(k)
            do kn = 1,numh(k)
              j = listh(indk+kn)
              if (indon(j).eq.0) then
                indon(j) = 1
                nind = nind + 1
                nindv(nind) = j
              endif
            enddo
          endif
        enddo
        if (lbandloc) then
          nloc = nloc + 1
          numf(nloc) = nind
        endif
        maxnf = max(maxnf,nind)
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          indon(j) = 0
        enddo
        nind = 0
      enddo

C Resize arrays if necessary
      call re_alloc(nindv,1,max(maxnf,nbasis),name='nindv')
      call re_alloc(listf,1,maxnf,1,nbandsloc,name='listf')

C Find out C control vectors - store elements
      do il = 1,nbandsloc
        do in = 1,numct(il)
          k = ncT2P(listct(in,il))
          if (k.gt.0) then
            indk = listhptr(k)
            do kn = 1,numh(k)
              j = listh(indk+kn)
              if (indon(j).eq.0) then
                indon(j) = 1
                nind = nind + 1
                nindv(nind) = j
              endif
            enddo
          endif
        enddo
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          indon(j) = 0
          listf(in,il) = j
        enddo
        nind = 0
      enddo

      end subroutine axb_build1

      subroutine axb_build2(nbands,maxnf,numf,listf,no_cl,
     .                      ncmax,numc,listc,indon,maxnhij)
C ********************************************************************
C Constructs control indexes of a C matrix in sparse form,
C C being the product of A and B (also in sparse form)
C
C              C = A x B
C
C In full form: A is rectangular, and has dimension:  nbands x nbasis
C               B is rectangular, and has dimension:  nbasis x nbands
C and, as a result:
C               C is rectangular, and has dimension:  nbands x nbands
C
C Written by P.Ordejon. October'96
C ***************************** INPUT *********************************
C integer nbands               : Number of rows of A 
C integer nbasis               : Number of columns of A
C integer maxnf             : First dimension of A matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of A)
C integer numf(nbands)         : Control vector of A matrix
C                            (number of nonzero elements of each row of A)
C integer listf(maxnf,nbands)  : Control vector of A matrix
C                           (list of nonzero elements of each row of A)
C integer nbasis               : Number of rows of B
C integer nbands               : Number of columns of B
C integer ncmax             : First dimension of B matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of B)
C integer numb(nbasis)          : Control vector of B matrix
C                            (number of nonzero elements of each row of B)
C integer listb(ncmax,nbasis)   : Control vector of B matrix
C                            (list of nonzero elements of each row of B)
C integer indon(nbands)     : Auxiliary array to build C in sparse form
C integer nindv(ncmax)      : Auxiliary array to store indexes of nonzero
C                             matrix elements of each row of C
C ************************INPUT/OUTPUT ********************************
C integer maxnhij             : First dimension of C matrix in sparse form,
C                             as declared in calling routine
C                             (max. number of <>0 elements of each row of C)
C **************************** OUTPUT *********************************
C integer numc(nbands)          : Control vector of C matrix
C                            (number of nonzero elements of each row of C)
C integer listc(ncmax,nbands)   : Control vector of C matrix
C                            (list of nonzero elements of each row of C)
C *********************************************************************

      use alloc
      use on_core, only : nindv, numhij, listhij
      use on_main, only : ncG2L, nbandsloc

      implicit none

      integer, intent(in) ::
     .  nbands,maxnf,ncmax,no_cl,
     .  listf(maxnf,nbandsloc),listc(ncmax,no_cl),
     .  numf(nbandsloc),numc(no_cl)
      integer, intent(inout) :: maxnhij
      integer, intent(out)   :: indon(nbands)

C Internal variables..................................................
      integer
     .  i,il,in,j,k,kn,nind

C Ensure arrays are large enough on entry
      call re_alloc(nindv,1,max(ncmax,nbands),name='nindv')

C Initialize internal variables
      nind = 0
      do i = 1,nbands
        indon(i) = 0
      enddo
      do i = 1,ncmax
        nindv(i) = 0
      enddo

C Find out C dimensions
      do il = 1,nbandsloc
        do in = 1,numf(il)
          k = ncG2L(listf(in,il))
          do kn = 1,numc(k)
            j = listc(kn,k)
            if (indon(j) .eq. 0) then
              indon(j) = 1
              nind = nind + 1
              nindv(nind) = j
            endif
          enddo
        enddo
        maxnhij = max(nind,maxnhij)
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          indon(j) = 0
        enddo
        nind = 0
      enddo

C  Re-size arrays
      call re_alloc(numhij,1,nbandsloc,name='numhij')
      call re_alloc(listhij,1,maxnhij,1,nbandsloc,name='listhij')

C Find out C control vectors 
      do il = 1,nbandsloc
        do in = 1,numf(il)
          k = ncG2L(listf(in,il))
          do kn = 1,numc(k)
            j = listc(kn,k)
            if (indon(j) .eq. 0) then
              indon(j) = 1
              nind = nind + 1
              nindv(nind) = j
            endif
          enddo
        enddo
        numhij(il) = nind
        do in = 1,nind
          j = nindv(in)
          nindv(in) = 0
          indon(j) = 0
          listhij(in,il) = j
        enddo
        nind = 0
      enddo

      end  subroutine axb_build2

      end module m_on_subs

