! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine extrapolon(istep,iord,nspin,nbasis,nbasisCloc,maxnc,
     .                      numc,listc,aux,numcold,listcold,cold,c,
     .                      nbasisloc)
C ******************************************************************************
C Subroutine to extrapolate a given matrix M (like the coefficients of the
C wave functions, or the density matrix) for the next MD step.
C The matrix M is given in sparse form.
C
C Writen by P.Ordejon, November'96.
C ******************************* INPUT ***************************************
C integer istep                : Time step of the simulation
C integer iord                 : Extrapolation order (0 or 1)
C                                0 = 0th order;  1 = 1st order
C integer nspin                : Number of spin polarizations (1 or 2)
C integer nbasis               : Number of rows of matrix M
C integer nbasisloc            : Number of rows of matrix M held locally
C integer nbasisCloc           : Maximum number of rows of matrix M (dimension)
C integer maxnc                : First dimension of M matrix, and maximum
C                                number of nonzero elements of each column of M
C integer numc(maxnc)          : Control vector 1 of M matrix at t
C integer listc(maxnc,nbasisCloc) : Control vector 2 of M matrix at t
C real*8 aux(2,nbasis)         : Auxiliary storage array
C ************************** INPUT AND OUTPUT *********************************
C integer numcold(maxnc)         : Input: Control vector 1 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 1 of M matrix at t
C integer listcold(maxnc,nbasisCloc) : Input: Control vector 2 of M matrix at t-dt
C                                       (if istep .ne. 1)
C                                Output: Control vector 2 of M matrix at t
C real*8 cold(maxnc,nbasisCloc,nspin): Input: matrix M at t-2dt
C                                Output: matrix M at t-dt
C real*8 c(maxnc,nbasisCloc,nspin): New matrix M (extrapolated)
C                                Input: matrix at t-dt
C                                Output: matrix at t
C                                If istep = 1, c returned uncahanged
C **************************** SCRATCH SPACE **********************************
C real*8 aux(2,nbasis)         : Auxiliary storage array
C **************************** BEHAVIOUR **************************************
C The routine allows for the sparse structure of the matrix M to change
C between MD time steps. On input, the matrices of former steps (c and cold) 
C have the structure of last step (t-dt): numold and listold; whereas the new
C (extrapolated) matrix has the structure of the current time step (which
C must be determined before calling this routine!!): num and list.
C On output, the routine updates the structure of c and cold, to that
C at the current (t) time steps respectively. Same with numold and listold
C 
C For the first MD time step (istep = 1), there is no extrapolation. 
C In that case, c is returned unchanged.
C Also, in that case numold and listold are only an output, and are set equal
C to num and list
C *****************************************************************************

      use precision, only : dp
      use parallel,  only : IONode
      use sys,       only : die

      implicit none

      integer, intent(in) :: iord
      integer, intent(in) :: istep
      integer, intent(in) :: nspin
      integer, intent(in) :: nbasis
      integer, intent(in) :: nbasisloc
      integer, intent(in) :: nbasisCloc
      integer, intent(in) :: maxnc
      integer, intent(in) :: numc(nbasisCloc)
      integer, intent(in) :: listc(maxnc, nbasisCloc)

      integer, intent(inout) :: numcold(nbasisloc)
      integer, intent(inout) :: listcold(maxnc, nbasisloc)
      real(dp), intent(inout) :: cold(maxnc, nbasisloc, nspin)
      real(dp), intent(inout) :: c(maxnc, nbasisCloc, nspin)

      real(dp), intent(out) :: aux(2,nbasis)
C  Internal variables .......................................................

      integer
     .  i,in,ispin,j

      real(dp) ::
     .  msave

      logical :: 
     .  changed
C ...........................................................................

      if (iord /= 0 .and. iord /= 1) then
        if (IONode) then
          call die ('extrapolon: Wrong iord: '//
     .              'only 0 and 1 order available')
        endif
      endif

C Just initialize numcold and listcold if istep = 1 ...........................
      if (istep .eq. 1) then
        do i = 1,nbasisloc
          numcold(i) = numc(i)
          do in = 1,numc(i)
            listcold(in,i) = listc(in,i)
            do ispin = 1,nspin
              cold(in,i,ispin) = 0.0_dp
            enddo
          enddo
        enddo

      else

C Check if sparse structure has changed .....................................
        changed = .false.

        do i = 1,nbasisloc
          if (numcold(i).ne.numc(i)) then
            changed = .true.
            exit
          endif
          do in = 1,numc(i)
            if (listcold(in,i).ne.listc(in,i)) then
              changed = .true.
              exit
            endif
          enddo
        enddo

C If sparse structure has changed, re-order c and cold 
C and change numcold and listcold to current ones .............................

        if (changed) then
          aux = 0.0_dp
          
          do i = 1,nbasisloc
            do ispin = 1,nspin
              do in = 1,numcold(i)
                j = listcold(in,i)
                aux(1,j) = c(in,i,ispin)
                aux(2,j) = cold(in,i,ispin)
              enddo
              do in = 1,numc(i)
                j = listc(in,i)
                c(in,i,ispin) = aux(1,j)
                cold(in,i,ispin) = aux(2,j)
              enddo
              do in = 1,numcold(i)
                j = listcold(in,i)
                aux(1,j) = 0.0_dp
                aux(2,j) = 0.0_dp
              enddo
            enddo
            numcold(i) = numc(i)
            do in = 1,numc(i)
              listcold(in,i) = listc(in,i)
            enddo
          enddo
        endif !changed

C Extrapolate matrix M ......................................................

        do ispin = 1,nspin
          do i = 1,nbasisloc
            do in = 1,numc(i)
              msave = c(in,i,ispin)
              if (iord == 1) then
                c(in,i,ispin) = 2.0_dp*c(in,i,ispin) - cold(in,i,ispin)
              endif
              cold(in,i,ispin) = msave
            enddo
          enddo
        enddo

      endif !istep==1

      end
