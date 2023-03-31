c**********************************************************************
c   pi.f - compute pi by integrating f(x) = 4/(1 + x**2)     
c     
c   Each node: 
c    1) receives the number of rectangles used in the approximation.
c    2) calculates the areas of it's rectangles.
c    3) Synchronizes for a global summation.
c   Node 0 prints the result.
c
c  Variables:
c
c    pi  the calculated result
c    n   number of points of integration.  
c    x           midpoint of each rectangle's interval
c    f           function to integrate
c    sum,pi      area of rectangles
c    tmp         temporary scratch space for global summation
c    i           do loop index
c****************************************************************************
      program main

      use mpi_siesta

      implicit none

      integer, parameter :: sp = selected_real_kind(6,30)
      integer, parameter :: dp = selected_real_kind(14,100)

      real(dp), parameter ::  PI25DT = 3.141592653589793238462643_dp

      real(sp) t1, t2, elapsed, tot_time
      real(dp)  mypi, pi, h, sum, x, f, a
      integer n, myid, numprocs, i, rc, ntotal, ierr
      integer count, rate, max, count2
c                                 function to integrate
      f(a) = 4.d0 / (1.d0 + a*a)

c
c     Some timing calls for 'proof of concept' work
c
      call system_clock(count,rate,max)
      write(6,*) 'Initial clock count: ',count
c
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c
      call system_clock(count2,rate,max)
      print *, 'Process ', myid, ' of ', numprocs,                        &
     &         ' Took ', count2-count, ' ticks to complete setup,'

 10   if ( myid .eq. 0 ) then
cc         write(6,98)
cc 98      format('Enter the number of intervals: (0 quits)')
         n= 500000000
      endif
      
      call system_clock(count,rate,max)
c
c     Send both n and h (redundant, but it tests the interface...)
c
      call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
      h = 1.0d0/n
      call MPI_BCAST(h,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
c
      call system_clock(count2,rate,max)
      print *, 'Process ', myid,                                         &
     &         ' Took ', count2-count, ' ticks to complete bcast.'


      call cpu_time(t1)
      sum  = 0.0d0
      do 20 i = myid+1, n, numprocs
         x = h * (dble(i) - 0.5d0)
         sum = sum + f(x)
 20   continue
      mypi = h * sum
      call cpu_time(t2)
      elapsed = t2 - t1

c                                 collect all the partial sums
      call system_clock(count,rate,max)

      call MPI_REDUCE(mypi,pi,1,mpi_double_precision,MPI_SUM,0,           
     &     MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(elapsed,tot_time,1,mpi_real,MPI_SUM,0,  
     &     MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(n,ntotal,1,MPI_INTEGER,MPI_SUM,0,         
     &     MPI_COMM_WORLD,ierr)

      call system_clock(count2,rate,max)
      print *, 'Process ', myid,                                
     &         ' Took ', count2-count, ' ticks to complete reduces.'


c                                 node 0 prints the answer.
      if (myid .eq. 0) then
         write(6,*) pi, abs(pi - PI25DT)
         write(6,*) 'Total time spent in computation:', tot_time
         write(6,*) 'Total number of points:', ntotal
      endif

 30   call MPI_FINALIZE(rc)
      stop
      end




