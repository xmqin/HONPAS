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
c **********************************************************************
c Reads and imposes required constraints to atomic displacements by
c making zero the forces in those directions. Constraints are specified
c by the FDF data block GeometryConstraints (see example below).
c Only types position and routine implemented in this version.
c Written by J.M.Soler. Feb., 1998
c Modified by P. Ordejon to output the total number of constraints
C    imposed.  June, 2003
c *********** INPUT ****************************************************
c real*8  cell(3,3)    : Lattice vectors
c real*8  stress( 3,3) : Stress tensor
c integer na           : Number of atoms
c integer isa(na)      : Species indexes
c real*8  amass(na)    : Atomic masses
c real*8  xa(3,na)     : Atomic cartesian coordinates
c real*8  fa(3,na)     : Atomic forces
c *********** OUTPUT ***************************************************
c real*8  cstress( 3,3) : Constrained stress tensor
c real*8  cfa(3,na)     : Constrained atomic forces
c integer ntcon         : Total number of position constraints imposed
c *********** UNITS ****************************************************
c Units are arbitrary but cell and xa must be in the same units
c *********** BEHAVIOUR ************************************************
c cstress may be the same physical array as stress, and cfa the same 
c as fa, i.e. it is allowed:
c     call fixed( cell, stress, na, isa, amass, xa, fa, stress, fa, ntcon )
C NOTE: This is NOT good practice, and might be flagged by some compilers
c *********** USAGE ****************************************************
c Example: consider a diatomic molecule (atoms 1 and 2) above a surface, 
c represented by a slab of 5 atomic layers, with 10 atoms per layer.
c To fix the cell height, the slab's botom layer (last 10 atoms),
c the molecule's interatomic distance, its height above the surface
c (but not its azimutal orientation and lateral position), and the
c relative height of the two atoms:
c
c   %block GeometryConstraints
c   cellside   c 
c   cellangle  alpha  beta  gamma
c   position  from -1 to -10
c   rigid  1  2
c   center 1  2   0.0  0.0  1.0
c   routine constr
c   stress 1  2  3
c   %endblock GeometryConstraints
c
c where constr is the following user-written subroutine:
c
c      subroutine constr( cell, na, isa, amass, xa, stress, fa, ntcon )
cc real*8  cell(3,3)    : input lattice vectors (Bohr)
cc integer na           : input number of atoms
cc integer isa(na)      : input species indexes
cc real*8  amass(na)    : input atomic masses
cc real*8  xa(3,na)     : input atomic cartesian coordinates (Bohr)
cc real*8  stress( 3,3) : input/output stress tensor (Ry/Bohr**3)
cc real*8  fa(3,na)     : input/output atomic forces (Ry/Bohr)
cc integer ntcon        : total number of position constraints, accumulative
c      integer na, isa(na), ntcon
c      double precision amass(na), cell(3,3), fa(3,na),
c     .                 stress(3,3), xa(3,na), fz
c      fz =  fa(3,1) + fa(3,2) 
c      fa(3,1) = fz*amass(1)/(amass(1)+amass(2))
c      fa(3,2) = fz*amass(2)/(amass(1)+amass(2))
c      ntcon = ntcon+1
c      end
c **********************************************************************

      subroutine fixed( cell, stress, na, isa, amass, xa, fa,
     .                  cstress, cfa, ntcon )
      use precision
      use parallel,   only : Node
      use fdf
      use parsing
      use alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in)   :: na, isa(na)
      real(dp), intent(in)  :: amass(na), cell(3,3), 
     .                         fa(3,na), stress(3,3), xa(3,na)
      integer, intent(out)  :: ntcon
      real(dp), intent(out) :: cfa(3,na), cstress(3,3)

      ! Internal parameters
      ! maxl  : maximum number of input constraint lines
      ! maxw  : maximum number of items in an input constraint line
      integer, parameter :: maxl= 10000,  maxw=  1000


      ! Internal variables
      logical found
      character line*130, name1*10, name2*10, name3*10, 
     &          name4*10, names*130


      integer i, ia, ia1, ia2, ia3, ic, il, iu, ix,
     &        jx, lastch, lch(0:maxw), ni, nn, nr, nv, maxcnew

      character(len=10), pointer, save :: ctype(:)
      integer,           pointer, save :: iac(:)
      real(dp),          pointer, save :: xc(:,:)


      integer   integs(maxw)
      real(dp)  ddot, fxc, reals(maxw), values(maxw), xnorm


      external ddot, constr, memory

      logical,  save :: frstme
      data              frstme /.true./

      ! maxc  : maximum number of constraints
      integer, save :: maxc = 1000

      integer,  save :: nc
      real(dp), save :: xs(6)


#ifdef MPI
      integer MPIerror
#endif
      ! Read constraint data only the first time
      if (frstme) then
        ! Allocate memory for constraints - do only once and preserve contents
         nullify( ctype )
         call re_alloc( ctype, 1, maxc, name='ctype', routine='fixed' )

        nullify( iac )
        call re_alloc( iac, 1, maxc, name='iac', routine='fixed' )

        nullify( xc )
        call re_alloc( xc, 1, 3, 1, maxc, name='xc', routine='fixed' )

        ! Initialise stress constraints to unconstrained state
        xs(1:6) = 0.0d0

        nc = 0

        ! Look for constraints data block
        if (Node.eq.0) then
          found = fdf_block('GeometryConstraints',iu)
        endif
#ifdef MPI
        call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
        if (.not.found) goto 30

        ! Loop on data lines
        do 20 il = 1,maxl

          ! Read and parse data line
          if(Node.eq.0) then
            read(iu,'(a)',end=30) line
            lastch = index(line,'#') - 1
          endif
#ifdef MPI
       call MPI_Bcast(line,130,MPI_character,0,MPI_Comm_World,MPIerror)
       call MPI_Bcast(lastch,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          if (lastch .le. 0) lastch = len(line)
          call parse( line(1:lastch), nn, lch, names, nv, values,
     &                ni, integs, nr, reals )

          ! Check if constraints are finished
          name1 = names(lch(0)+1:lch(1))
          if (name1 .eq. '%end' .or.
     &        name1 .eq. '%endblock') then
            goto 30

          ! Select type of constraint
          elseif (name1 .eq. 'routine') then
            if (nn.eq.1 .or. names(lch(1)+1:lch(2)) .eq. 'constr') then
              nc = nc + 1
              if (nc.eq.maxc) then
                maxcnew = nc + nint(0.1*nc)

                call re_alloc( ctype, 1, maxcnew, name='ctype',            
     &                         routine='fixed', copy=.true. )
                call re_alloc( iac, 1, maxcnew, name='iac',            
     &                          routine='fixed', copy=.true. )
                call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                         routine='fixed', copy=.true. )
                maxc = maxcnew
              endif
              ctype(nc) = 'routine'
            else
              if (Node.eq.0) then
                write(6,*) 'fixed: ERROR: user-constraints routine',
     &                   ' must be called constr'
              endif
            endif

          elseif (name1 .eq. 'position') then

            ! Check syntax
            if (nr.ne.0 .and. nr.ne.3) then
              if (Node.eq.0) then
                write(6,'(a,/,a)')
     &            'fixed: syntax ERROR in %block GeometryConstraints:',
     &            line(1:lastch)
              endif
              goto 20
            endif

            ! Find constrained atoms
            if (nn .gt. 1) then
              ! Atoms specified by range. Make list of them.
              name2 = names(lch(1)+1:lch(2))
              name3 = names(lch(2)+1:lch(3))
              if (nn.eq.4) name4 = names(lch(3)+1:lch(4))
              if (name2.eq.'from' .and. name3.eq.'to') then
                ia1 = integs(1)
                ia2 = integs(2)
                if (ia1 .lt. 0) ia1 = na + ia1 + 1
                if (ia2 .lt. 0) ia2 = na + ia2 + 1
                if (nn.eq.4 .and. name4.eq.'step') then
                  ia3 = abs(integs(3))
                else
                  ia3 = 1
                endif
                ni = 0
                do ia = min(ia1,ia2),max(ia1,ia2),ia3
                  ni = ni + 1
                  integs(ni) = ia
                enddo
              else
                if (Node.eq.0) then
                  write(6,'(a,/,a)')
     &            'fixed: syntax ERROR in %block GeometryConstraints:',
     &            line(1:lastch)
                endif
                goto 20
              endif
            endif

            ! Store position constraints
            do i = 1,ni
              ia = integs(i)
              if (ia .lt. 0) ia = na + ia + 1
              if (nr .eq. 0) then
                ! Make one constraint for each cartesian coordinate
                do ix = 1,3
                  nc = nc + 1
                  if (nc.eq.maxc) then
                    maxcnew = nc + nint(0.1*nc)
                    call re_alloc( ctype, 1, maxcnew, name='ctype',            
     &                             routine='fixed', copy=.true. )
                    call re_alloc( iac, 1, maxcnew, name='iac',            
     &                             routine='fixed', copy=.true. )
                    call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                             routine='fixed', copy=.true. )
                    maxc = maxcnew
                  endif
                  ctype(nc) = 'position'
                  iac(nc) = ia
                  do jx = 1,3
                    xc(jx,nc) = 0.0d0
                  enddo
                  xc(ix,nc) = 1.0d0
                enddo
              elseif (nr .eq. 3) then
                ! Make only one constraint in the specified direction
                nc = nc + 1
                if (nc.eq.maxc) then
                  maxcnew = nc + nint(0.1*nc)
                  call re_alloc( ctype, 1, maxcnew, name='ctype',            
     &                           routine='fixed', copy=.true. )
                  call re_alloc( iac, 1, maxcnew, name='iac',            
     &                           routine='fixed', copy=.true. )
                  call re_alloc( xc, 1, 3, 1, maxcnew, name='xc',            
     &                         routine='fixed', copy=.true. )
                  maxc = maxcnew
                endif
                ctype(nc) = 'position'
                iac(nc) = ia
                xnorm = sqrt(reals(1)**2 + reals(2)**2 + reals(3)**2)
                do ix = 1,3
                  xc(ix,nc) = reals(ix) / xnorm
                enddo
              endif
            enddo

          elseif (name1 .eq. 'stress') then

            do ic = 1,nv
              nc = nc + 1
              ctype(nc) = 'stress'
              ix = nint(values(ic))
              if (ix.ge.1.and.ix.le.6) xs(ix) = 1.0d0
            enddo

          else
            if (Node.eq.0) then
              write(6,*) 'fixed: ERROR: sorry, constraint type ',
     &                  name1, ' not implemented yet'
            endif
          endif
   20   continue
   30   continue
        frstme = .false.
      endif

      ! Copy stress and forces to output arrays 
      do ix = 1,3
        do jx = 1,3
          cstress(jx,ix) = stress(jx,ix)
        enddo
      enddo
      do ia = 1,na
        do ix = 1,3
          cfa(ix,ia) = fa(ix,ia)
        enddo
      enddo

      ! Apply constraints
      ntcon = 0
      do ic = 1,nc
        if (ctype(ic) .eq. 'routine') then
          call constr( cell, na, isa, amass, xa, cstress, cfa, ntcon )
        elseif (ctype(ic) .eq. 'position') then
          ntcon = ntcon + 1
          ia = iac(ic)
          fxc = ddot(3,cfa(1,ia),1,xc(1,ic),1)
          do ix = 1,3
            cfa(ix,ia) = cfa(ix,ia) - fxc * xc(ix,ic)
          enddo
        endif
      enddo

      ! Apply stress constraints
      cstress(1,1) = cstress(1,1) - xs(1)*cstress(1,1)
      cstress(2,2) = cstress(2,2) - xs(2)*cstress(2,2)
      cstress(3,3) = cstress(3,3) - xs(3)*cstress(3,3)
      cstress(3,2) = cstress(3,2) - xs(4)*cstress(3,2)
      cstress(2,3) = cstress(2,3) - xs(4)*cstress(2,3)
      cstress(3,1) = cstress(3,1) - xs(5)*cstress(3,1)
      cstress(1,3) = cstress(1,3) - xs(5)*cstress(1,3)
      cstress(1,2) = cstress(1,2) - xs(6)*cstress(1,2)
      cstress(2,1) = cstress(2,1) - xs(6)*cstress(2,1)

      end subroutine fixed
