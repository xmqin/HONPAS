! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      SUBROUTINE BONDS( CELL, NA, ISA, XA, RMAX, filename )

C **********************************************************************
C Finds the bond lengths (actually the distances to neighbors)
C Based on shaper (by J.M.Soler)
C Alberto Garcia, Oct 2006
C ************ INPUT ***************************************************
C REAL*8  CELL(3,3) : Lattice vectors CELL(Ixyz,Ivector)
C INTEGER NA        : Number of atoms
C INTEGER ISA(NA)   : Species index of each atom
C REAL*8  XA(3,NA)  : Cartesian atomic coordinates
C REAL*8  RMAX      : Maximum distance considered
C CHARACTER(LEN=*)  filename : File name for output
C ************ UNITS ***************************************************
C CELL and XA must be in the same units
C **********************************************************************

      use precision, only : dp
      use atmfuncs,  only : labelfis
      use units,     only : Ang
      use sorting,   only : order, iorder, ordix
      use alloc,       only: re_alloc, de_alloc
      use neighbour,   only: jna=>jan, xij, r2ij, maxna=>maxnna
      use neighbour,   only: mneighb, reset_neighbour_arrays
      implicit          none

      INTEGER,  intent(in) ::     NA, ISA(NA)
      real(dp), intent(in) ::     CELL(3,3), XA(3,NA), RMAX
      character(len=*), intent(in) :: filename

      EXTERNAL             :: io_assign, io_close
      integer              ::  IA, IN, JA, JS, NNA, iu
      integer,     pointer ::  index(:)
      real(dp)             :: RIJ
      real(dp),  parameter :: tol = 1.0e-8_dp

      nullify(index)
      call re_alloc( index, 1, maxna, 'index', 'bonds' )

C Initialize neighbour-locater routine
      CALL MNEIGHB( CELL, RMAX, NA, XA, 0, 0, NNA )

C Main loop

      call io_assign( iu )
      open(iu,file=filename,form='formatted',
     $     status='replace', action="write")      

      do IA = 1,NA

C Find neighbours of atom IA
        CALL MNEIGHB( CELL, RMAX, NA, XA, IA, 0, NNA )
        ! This call will do nothing if the internal neighbor arrays
        ! did not grow. If tight size were important, one could
        ! resize to nna instead.
        call re_alloc( index, 1, maxna, 'index', 'bonds' )

           if (nna < 2 ) then
              write(iu , *) "Atom ", ia, 
     $             ": no bonds for rmax specified: ", rmax/Ang,
     $             " Ang."
              call pxfflush(iu)
              cycle   ! loop over ia
           endif
           ! Sort by distance
           call ordix(r2ij,1,nna,index)
           call iorder( jna, 1, nna, index )
           call order(  r2ij, 1, nna, index )
           call order(  xij, 3, nna, index )

           write(iu,fmt="(a,i5,1x,a,3f8.4)")
     $       "Neighbors of: ",
     $          ia, trim(labelfis(isa(ia))) // " at: ", xa(:,ia)

          do IN = 1, NNA
            JA = JNA(IN)
            JS = ISA(JA)
            RIJ = SQRT(R2IJ(IN))
            if (rij > 0.0001_dp) then
               write(iu,fmt="(i5,1x,a,f8.4,2x,a,3f8.4)")
     $           ja, labelfis(js), rij/Ang, "Ang. Really at: ",
     $           xa(:,ia)+xij(:,in)
            endif
            ! Note that we now print the real location of the
            ! neighbor atom in space, and not, as in earlier
            ! versions, the location of the equivalent representative
            ! in the unit cell.
          enddo

      enddo

      call reset_neighbour_arrays( )
      call de_alloc( index, 'index', 'bonds' )
      call io_close(iu)

      end

