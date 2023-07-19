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
      SUBROUTINE TIMER (PROG,IOPT)
      use parallel, only: IOnode
      use sys, only: die
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*) prog
      integer iopt
c
c     Two-in-one: CPU plus elapsed times 
c
#ifdef MPI
      integer MPIerror
#endif
      call timer1(prog,iopt)
      if (IOnode) then
         call elap1(prog,iopt)
      endif
      end
!----------------------------------------------------------------

      SUBROUTINE TIMER1 (PROG,IOPT)


C  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
C   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
C   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.
C  ARGUMENTS:
C    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
C    IOPT: INPUT OPTION PARAMETER:
C      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
C      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
C      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
C      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
C  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
C  WRITTEN BY J.SOLER (JSOLER AT EMDUAM11) DEC/90

C
C  Modules
C
      use precision
      use parallel,   only : Node
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      real :: treal     ! Default real for call to cpu_time

      integer :: iopt, iprog

      integer, parameter :: NMAX=500
      real(dp), parameter :: ZERO=0.0_dp,HUNDRD=100.0_dp,
     $                       TIMMIN=1.0e-6_dp
      real(dp) :: TIME1(NMAX),TIMET(NMAX)

      integer, save  :: nprogs = 0
      real(dp), save :: time0  = 0.0_dp

      real(dp) :: time, timtot, timetl, avgtme, fractn
#ifdef MPI
      real(dp) ::  buffer1(3),buffer2(3)
      integer mpierror
#endif
      INTEGER NCALLS(NMAX)
      CHARACTER*10 PROGS(NMAX),PROG*(*)
      SAVE PROGS,TIME1,TIMET,NCALLS

!!      CALL CPUTIM (TIME)
      call cpu_time( treal)         ! Standard Fortran95
      TIME = treal

      IF (IOPT.EQ.0) THEN
         NPROGS=0
         TIME0=TIME
      ELSEIF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         DO 10 IPROG=1,NPROGS
            IF (PROGS(IPROG).EQ.PROG) GO TO 20
   10    CONTINUE
            NPROGS=NPROGS+1
            IF (NPROGS.GT.NMAX) THEN
               if (Node.eq.0) then
                 WRITE (6,*) 'timer: NMAX IS SATURATED. PROG = ',PROG
               endif
               RETURN
            ENDIF
            IPROG=NPROGS
            PROGS(IPROG)=PROG
            NCALLS(IPROG)=0
            TIMET(IPROG)=ZERO
   20    CONTINUE
         IF (IOPT.EQ.1) THEN
            NCALLS(IPROG)=NCALLS(IPROG)+1
            TIME1(IPROG)=TIME
         ELSE
            TIMET(IPROG)=TIMET(IPROG)+TIME-TIME1(IPROG)
         ENDIF
      ELSEIF (IOPT.EQ.3) THEN
         TIMTOT=TIME-TIME0

C Sum TIMTOT across all Nodes and ensure that all Nodes have same value here
#ifdef MPI
         buffer1(1)=timtot
         call MPI_AllReduce(buffer1(1),timtot,1,MPI_double_precision,
     .     MPI_sum,MPI_Comm_World,MPIerror)
#endif

         IF (TIMTOT.LT.TIMMIN) RETURN

         IF (PROG.EQ.'ALL' .OR. PROG.EQ.'all') THEN
           if (Node.eq.0) then
             WRITE (6,'(/,A)') 'timer: CPU execution times:'
             WRITE (6,'(A,2X,A10,A9,2A12,A9)') 'timer:',
     .         'Routine   ', 'Calls', 'Time/call', 'Tot.time', '%'
           endif
           DO 40 IPROG=1,NPROGS
             TIMETL=TIMET(IPROG)
             AVGTME=TIMET(IPROG)/NCALLS(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=AVGTME
             buffer1(3)=FRACTN
             call MPI_Reduce(buffer1,buffer2,3,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             AVGTME=buffer2(2)
             FRACTN=buffer2(3)
#endif

             if (Node.eq.0) then
               WRITE(6,'(A,2X,A10,I9,2F12.3,F9.2)') 'timer:',
     .           PROGS(IPROG),NCALLS(IPROG),AVGTME,TIMETL,FRACTN
             endif
   40      CONTINUE
           if (Node.eq.0) then
             WRITE(6,*) ' '
           endif
         ELSE
           DO 50 IPROG=1,NPROGS
             IF (PROGS(IPROG).NE.PROG) GOTO 50
             TIMETL=TIMET(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

C Sum values across Nodes
#ifdef MPI
             buffer1(1)=TIMETL
             buffer1(2)=FRACTN
             call MPI_Reduce(buffer1,buffer2,2,MPI_double_precision,
     .         MPI_sum,0,MPI_Comm_World,MPIerror)
             TIMETL=buffer2(1)
             FRACTN=buffer2(2)
#endif

             if (Node.eq.0) then
               WRITE(6,'(A,A10,I6,F12.3,F7.2)')
     .          'timer: Routine,Calls,Time,% = ',
     .           PROGS(IPROG),NCALLS(IPROG),TIMETL,FRACTN
             endif
   50      CONTINUE
         ENDIF
      ELSE
         if (Node.eq.0) then
           WRITE(6,*) 'timer: INVALID OPTION IOPT =',IOPT
         endif
      ENDIF
      END
C
c-----------------------------------------------------------------
C
      SUBROUTINE elap1 (PROG,IOPT)


C  FINDS AND PRINTS THE WALL CLOCK SPENT IN DIFFERENT ROUTINES AND/OR
C   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=0 AT THE BEGINNING
C   OF EACH ROUTINE AND WITH IOPT=1 AT THE END OF IT.

!  Only for the Master Node in an MPI run

C  ARGUMENTS:
C    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
C    IOPT: INPUT OPTION PARAMETER:
C      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
C      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
C      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
C      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
C  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
C  WRITTEN BY Alberto Garcia, Feb 2000, stealing code from J. Soler
!  Maximum wall-clock time support by Andrew Walker, 2008

      use precision
      use m_walltime
      use parallel, only: node
      use sys, only : die
      

      implicit none

      integer :: iopt, iprog
      integer, parameter :: NMAX=500
      real(dp), parameter :: ZERO=0.0_dp,HUNDRD=100.0_dp,
     $                       TIMMIN=1.0e-6_dp
      real(dp) :: TIME1(NMAX),TIMET(NMAX)
      integer ::  wt = 6  ! Use standard output

      integer, save  :: nprogs = 0
      real(dp), save :: time0  = 0.0_dp

      real(dp) :: time, timtot, timetl, avgtme, fractn

      INTEGER NCALLS(NMAX)
      CHARACTER*10 PROGS(NMAX),PROG*(*)
      SAVE PROGS,TIME1,TIMET,NCALLS,wt

      call wall_time(time)

      IF (IOPT.EQ.0) THEN
         NPROGS=0
         TIME0=TIME
      ELSEIF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
         DO 10 IPROG=1,NPROGS
            IF (PROGS(IPROG).EQ.PROG) GO TO 20
   10    CONTINUE
            NPROGS=NPROGS+1
            IF (NPROGS.GT.NMAX) THEN
               if (Node.eq.0) then
                 WRITE (wt,*) 'elap: NMAX IS SATURATED. PROG = ',PROG
               endif
               RETURN
            ENDIF
            IPROG=NPROGS
            PROGS(IPROG)=PROG
            NCALLS(IPROG)=0
            TIMET(IPROG)=ZERO
   20    CONTINUE
         IF (IOPT.EQ.1) THEN
            NCALLS(IPROG)=NCALLS(IPROG)+1
            TIME1(IPROG)=TIME
         ELSE
            TIMET(IPROG)=TIMET(IPROG)+TIME-TIME1(IPROG)
         ENDIF
      ELSEIF (IOPT.EQ.3) THEN
         TIMTOT=TIME-TIME0

         IF (TIMTOT.LT.TIMMIN) RETURN

         IF (PROG.EQ.'ALL' .OR. PROG.EQ.'all') THEN
             WRITE (wt,'(/,A)') 'elaps: ELAPSED times:'
             WRITE (wt,'(A,2X,A10,A9,2A12,A9)') 'elaps:',
     .         'Routine   ', 'Calls', 'Time/call', 'Tot.time', '%'
           DO 40 IPROG=1,NPROGS
             TIMETL=TIMET(IPROG)
             AVGTME=TIMET(IPROG)/NCALLS(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

               WRITE(wt,'(A,2X,A10,I9,2F12.3,F9.2)') 'elaps:',
     .           PROGS(IPROG),NCALLS(IPROG),AVGTME,TIMETL,FRACTN
   40      CONTINUE
             WRITE(wt,*) ' '
         ELSE
           DO 50 IPROG=1,NPROGS
             IF (PROGS(IPROG).NE.PROG) GOTO 50
             TIMETL=TIMET(IPROG)
             FRACTN=TIMET(IPROG)/TIMTOT*HUNDRD

               WRITE(wt,'(A,A10,I6,F12.3,F7.2)')
     .          'elaps: Routine,Calls,Wall,% = ',
     .           PROGS(IPROG),NCALLS(IPROG),TIMETL,FRACTN
   50      CONTINUE
         ENDIF
      ELSE
           WRITE(wt,*) 'elap: INVALID OPTION IOPT =',IOPT
      ENDIF

      END


