! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
c----------------------------------------------------------------------
c State: stable          Date: 6/3/95
c----------------------------------------------------------------------

        subroutine pspltm1( nrow, ncol, mode, ja, ia, title, ptitle, 
     &                      size, munt, nlines, lines, iwa, compax, 
     &                      iunt )

c----------------------------------------------------------------------
c PSPLTM - PostScript PLoTer of a (sparse) Matrix. It draws the matrix
c          with a bounding box and centered.
c
c input arguments description :
c
c nrow   = number of rows in matrix
c
c ncol   = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored in 
c           CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2)
c           or MSC (mode=3)
c
c ja     = column indices of nonzero elements when matrix is
c          stored rowise. Row indices if stores column-wise.
c          If MSR/MSC format is used the argument must be
c          jlu(nrows+2).
c ia     = integer array of containing the pointers to the 
c          beginning of the rows/columns in ja. If MSR/MSC format
c          is used the argument must be jlu(1).
c
c title  = character*(*). a title of arbitrary length to be printed 
c          as a caption to the figure. Can be a blank character if no
c          caption is desired.
c
c ptitle = position of title; 0 under the drawing, else above
c
c size   = size of the drawing  
c
c munt   = units used for size : 'cm' or 'in'
c
c nlines = number of separation lines to draw for showing a partionning
c          of the matrix. Enter zero if no partition lines are wanted.
c
c lines  = integer array of length nlines containing the coordinates of 
c          the desired partition lines . The partitioning is symmetric: 
c          a horizontal line across the matrix will be drawn in 
c          between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
c          an a vertical line will be similarly drawn between columns
c          lines(i) and lines(i)+1 for i=1,2,...,nlines 
c
c iunt   = logical unit number where to write the matrix into.
c
c compax = Compression parameter. If compax = 0 => No compresion.
c          If compax > 0 => compression
c          Compax is the number of columns with zero value between two 
c          columns with non zero value drawed as entries on the 
c          matrix.
c
c iwa = integer work array. Its minimun lenght is nrow/ncol. 
c
c Note: Use of 'cm' assumes european format for paper size (21cm wide) 
c       and use of 'in' assumes american format (8.5in wide). The 
c       correct centering of the figure depends on the proper choice.
c----------------------------------------------------------------------

        implicit none
c------------------------------------------------------- Input variables
        integer nrow, ncol, mode, ptitle, iunt, nlines, compax

        integer ja(*), ia(*), lines(nlines), iwa(ncol+1) 

        real size

        character title*(*), munt*2 

C------------------------------------------------------- Local variables
        integer n, nr, nc, nnz, m, maxdim, istart, ilast, ii, k, ltit

        real lrmrgn, botmrgn, xtit, ytit, ytitof, 
     &       fnstit, siz

        real xl, xr, yb, yt, scfct, u2dot, frlw, delt, paperx, 
     &       xx, yy

        integer kk, isep, kol, mj, j
C       integer nii

        real haf, zero, conv
        data haf /0.5/, zero/0.0/, conv/2.54/

        integer LENSTR
        external LENSTR


C---------------------------------------------------------------- BEGIN
        siz = size
        nr = nrow
        nc = ncol

        if ( (mode .eq. 0) .OR. (mode .eq. 2) ) then
          n = nr
        else if ( (mode .eq. 1) .OR. (mode .eq. 3) ) then
          n = nc
        else
          write(*,*) 'Error in mode parameter'
          STOP 'ERROR'
        endif

        nnz = ia(n+1) - ia(1) 
        maxdim = MAX(nrow, ncol)
        m = 1 + maxdim
        nc = nc + 1
        nr = nr + 1

c. . .  units (cm) to dot conversion factor and paper size (cm)
        u2dot = 72.0/conv

        if ( (munt .eq. 'cm') .OR. (munt .eq. 'CM') ) then
            paperx = 21.0

        else if ( (munt .eq. 'in') .OR. (munt .eq. 'IN') ) then
            paperx = 8.5*conv
            siz = siz*conv
        else
            write(*,*) 'Error in munt parameter'
            STOP 'ERROR'
        endif

c. . .  left and right margins (drawing is centered)
        if (siz .gt. (paperx+2)) then
            write(*,*) 'Size of the drawing too big'
            STOP 'ERROR'
        endif

        lrmrgn = (paperx-siz)/2.0

c. . .  bottom margin : 2 cm
        botmrgn = 2.0

c. . .  scaling factor
        scfct = (siz*u2dot) / m

c. . .  matrix frame line witdh
        frlw = 0.25

c. . .  font size for title (cm)
        fnstit = 0.5
        ltit = LENSTR(title)

c. . .  position of title : centered horizontally
c                           at 1.0 cm vertically over the drawing
        ytitof = 1.0
        xtit = paperx / 2.0
        ytit = botmrgn + (siz*nr)/m + ytitof

c. . .  almost exact bounding box
        xl = lrmrgn*u2dot - (scfct*frlw)/2
        xr = (lrmrgn+siz)*u2dot + (scfct*frlw)/2
        yb = botmrgn*u2dot - (scfct*frlw)/2
        yt = (botmrgn + (siz*nr)/m) * u2dot + (scfct*frlw)/2
        if (ltit .gt. 0) then
            yt = yt + (ytitof+fnstit*0.70)*u2dot
        endif

c. . .  add some room to bounding box: 10 dpi = 3.5 mm
c       delt = 0.0
        delt = 10.0
        xl = xl - delt
        xr = xr + delt
        yb = yb - delt
        yt = yt + delt

c. . .  correction for title under the drawing
        if ( (ptitle .eq. 0) .AND. (ltit .gt. 0) ) then
            ytit = botmrgn + fnstit*0.3
            botmrgn = botmrgn + ytitof + fnstit*0.7
        endif

c.......BEGIN OF OUTPUT
        write(iunt,10) '%!'
        write(iunt,10) '%%Creator: PSPLTM routine'
        write(iunt,12) '%%BoundingBox:', xl, yb, xr, yt
        write(iunt,10) '%%EndComments'
        write(iunt,10) '/cm {72 mul 2.54 div} def'
        write(iunt,10) '/mc {72 div 2.54 mul} def'
        write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
        write(iunt,10) 'cvs print ( ) print} def'
        write(iunt,10)
     &          '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'

c. . .  We leave margins etc. in cm so it is easy to modify them if
c       needed by editing the output file
        write(iunt,10) 'gsave'

        if (ltit .gt. 0) then
            write(iunt,*) 
     &              '/Helvetica findfont',fnstit,' cm scalefont setfont'
            write(iunt,*) xtit,' cm',ytit,' cm moveto'
            write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
        endif

        write(iunt,*) lrmrgn, ' cm ', botmrgn, ' cm translate'
        write(iunt,*) siz, ' cm ', m, ' div dup scale '

c. . .  draw a frame around the matrix
        write(iunt,*) frlw,' setlinewidth'
        write(iunt,10) 'newpath'
        write(iunt,11) 0, 0, ' moveto'

c. . .  change .false. to .true. if you prefer a square frame around
c       a rectangular matrix

        if (.false.) then
            write(iunt,11) m, 0, ' lineto'
            write(iunt,11) m, m, ' lineto'
            write(iunt,11) 0, m, ' lineto'
        else
            write(iunt,11) nc,  0, ' lineto'
            write(iunt,11) nc, nr, ' lineto'
            write(iunt,11)  0, nr, ' lineto'
        endif

        write(iunt,10) 'closepath stroke'

c. . . drawing the separation lines 
 
       write(iunt,*)  ' 0.2 setlinewidth'
      do kol= 1, nlines 
        isep = lines(kol) 
c. . . . . horizontal lines 
        yy = REAL(nrow-isep) + haf 
        xx = REAL(ncol+1) 
        write(iunt,13) zero, yy, ' moveto '
        write(iunt,13) xx, yy, ' lineto stroke '
c. . . . . vertical lines 
        xx = REAL(isep) + haf 
        yy = REAL(nrow+1)  
        write(iunt,13) xx, zero,' moveto '
        write(iunt,13) xx, yy, ' lineto stroke '             
      enddo

c------------------------------------------ PLOTTING LOOP
      write(iunt,10) '1 1 translate'
      write(iunt,10) '0.8 setlinewidth'
      write(iunt,10) '/p {moveto 0 -.40 rmoveto '
      write(iunt,10) '           0  .80 rlineto stroke} def'
      write(iunt,10) '/bm {moveto -.40 0 rmoveto} def'
      write(iunt,10) '/sm {0 rlineto .40 0 rlineto stroke} def'

C     nii = nrow/10

      DO ii= 1, n

          do k= 1, nc
              iwa(k) = 0
          enddo

          istart = ia(ii)
          ilast  = ia(ii+1)-1 

          IF ( (mode .eq. 1) .OR. (mode .eq. 3) ) THEN
c. . . . . .  CSC or MSC
              do k= istart, ilast
                  write(iunt,11) ii-1, nrow-ja(k), ' p'
              enddo

c. . . . . .  add diagonal element if MSC mode.
              if (mode .eq. 3) then
                  write(iunt,11) ii-1, nrow-ii, ' p' 
              endif

          ELSE
c. . . . . .  CSR or MSR
              mj = 0
              do k= istart, ilast
                  kk = ja(k)
                  iwa(kk) = 1
                  mj = MAX( mj, kk )
              enddo
              mj = MIN( mj, ncol )

              k = 1
              dowhile( (k .le. mj) .AND. (compax .ge. 1) )
                  do while( (iwa(k) .eq. 0) .AND. (k .le. mj) )
                      k = k + 1
                  enddo

                  do while( (iwa(k) .eq. 1) .AND. (k .le. mj) )
                      k = k + 1
                  enddo
                  kk = k
                  do while( (iwa(kk) .eq. 0) .AND. (kk .le. mj) )
                      kk = kk + 1
                  enddo
                  kk = kk - 1

                  if ( ((kk-k+1).le.compax) .AND. (kk.le.mj) ) then
                      do j= k, kk
                          iwa(j) = 1
                      enddo
                  endif
                  k = kk + 1
              enddo 

              k = 1
              dowhile (k .le. mj)
                  do while( (iwa(k) .eq. 0) .AND. (k .le. mj) )
                      k = k + 1
                  enddo

                  kk = k
                  do while ( (iwa(kk) .eq. 1) .AND. (kk .le. mj) )
                      kk = kk + 1
                  enddo

                  if (k .le. mj) then 
                      if ((kk-1) .eq. k) then
                          write(iunt,11) k-1, nrow-ii, ' p'
                      else        
                          write(iunt,11) k-1, nrow-ii, ' bm'
                          write(iunt,14) kk-1-k,' sm'
                      endif
                  endif
                  k = kk
              enddo 

c. . . . . .  add diagonal element if MSR mode.
              if (mode .eq. 2) then
                  write(iunt,11) ii-1, nrow-ii, ' p' 
              endif
          ENDIF
      ENDDO
c------------------------------------ END PLOTING LOOP
      write(iunt,10) 'showpage'
      RETURN

10    FORMAT (A)
11    FORMAT (2I6,A)
12    FORMAT (A,4F9.2)
13    FORMAT (2F9.2,A)
14    FORMAT (I6,A)
c------------------------------------------------------------------END
      end
