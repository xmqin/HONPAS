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
      subroutine eggbox(task,cell,na,isa,ntm,xa,fa,Etot,eggbox_block)

c **********************************************************************
c eggbox - subtracts the eggbox energy and forces in the atomic 
c          approximation, by means of species-dependent fourier series.
c If including the free-atom total energy of every species as zeroth 
c          component, the KS energy becomes cohesive energy.
c
c Written by T. Archer June 2003
c *********** INPUT ****************************************************
c real*8  cell(3,3)    : Lattice vectors
c integer na           : Number of atoms
c integer isa(na)      : Species indexes
c integer ntm(3)       : real space grid
c real*8  xa(3,na)     : Atomic cartesian coordinates
c real*8  fa(3,na)     : Atomic forces
c logical eggbox_block : if the eggboox block appears in the fdf
c *********** OUTPUT ***************************************************
c real*8  fa(3,na)     : Atomic forces
c real*8  energy 
c *********** UNITS ****************************************************
c Units of the fourier coefficents in eV internally and by default
c **********************************************************************

C  Modules

      use precision
      use parallel,    only : Node
      use fdf
      use parsing
#ifdef MPI
      use mpi_siesta
#endif

C External

      implicit          none

      integer           na, isa(na), ntm(3)
      real(dp)          cell(3,3), fa(3,na), xa(3,na), Etot
      character task*6 
      logical eggbox_block

c Internal parameters

c maxl     : maximum number block lines
c maxw     : maximum number of items in an input constraint line
c maxshell : maximum number of fouriour shell 
c maxispe  : maximum number of species
      integer maxl, maxw, maxshell, maxispe
      parameter ( maxl = 10000 )
      parameter ( maxw = 1000 )
      parameter ( maxshell = 10 )
      parameter ( maxispe = 15 )

C Internal variables  

      logical, save :: frstme = .true.

      integer  ni, nn, nr, nv, integers(maxw), lc(0:maxw),
     .     ia, ix, i, iu, a, b, c
      real(dp) reals(maxw), values(maxw),
     .     pi, xac(3), k(3), xa_frac(3,na), 
     .     celli(3,3), f(3,na), E, x, y, z, 
     .     fxyz(3,na), eV, escale,
     .     vreal(0:maxshell,0:maxshell,0:maxshell,maxispe)
      character         line*200, names*80, name1*10
#ifdef MPI
      integer  MPIerror, nvreal
#endif

      save vreal, escale, pi

c --------------------------------------------------------------

      if (frstme) then
         frstme = .false.
         eV = 13.60580_dp
         pi = 4.d0 * atan(1.d0)

C Initalise vreal
         do a=0,maxshell
           do b=0,maxshell
             do c=0,maxshell
               do ia=1,maxispe
                 vreal(a,b,c,ia)=0
               enddo
             enddo
           enddo
         enddo
         
C Read cooeficents from fdf
         if (Node.eq.0) then
           eggbox_block = fdf_block('EggBoxRemove',iu)
         endif
#ifdef MPI
         call MPI_Bcast(eggbox_block,1,MPI_logical,0,MPI_Comm_World,
     .                  MPIerror)
#endif
         if (eggbox_block) then
           if (Node.eq.0) then
             write(6,'(/a)') 'eggbox: Reading %block EggBoxRemove'
             do i=1,maxl
               read(iu,'(a)') line
               call parse(line, nn, lc, names, nv, values, ni, integers,
     .                    nr, reals)
C Check if block has finished
               name1 = names(lc(0)+1:lc(1))
               if (name1 .eq. '%end' .or.name1 .eq. '%endblock') goto 30
C vreal is set out as vreal(kx,ky,kz,species)
               vreal(integers(2),integers(3),integers(4),
     .              integers(1))=reals(1)
             enddo
 30          continue

C Get the scale (default 1 eV)
             escale = fdf_physical('EggBoxScale', 1.d0, 'eV')
             escale = escale/eV
           endif

C Pass data to other nodes
#ifdef MPI
           nvreal = maxispe*(maxshell+1)**3
           call MPI_Bcast(vreal(0,0,0,1),nvreal,MPI_double_precision,
     .                    0,MPI_Comm_World,MPIerror)
           call MPI_Bcast(escale,1,MPI_double_precision,0,
     .                    MPI_Comm_World,MPIerror)
#endif
         else
c          write(6,'a') "eggbox: The eggbox effect will not be removed"
           return
         endif
      endif
         

C Find reciprical lattice unit cell
      call reclat(cell, celli, 0)

C Convert to fractional coordinates
      do ia = 1,na
         do ix = 1,3
            xac(ix) = xa(ix,ia)
         enddo
         do ix = 1,3
            xa_frac(ix,ia) = celli(1,ix) * xac(1) +
     .              celli(2,ix) * xac(2) +
     .              celli(3,ix) * xac(3)
         enddo
      enddo
         
C Define the k vector describing the periodicity of the grid 
      k(1) = 2.0d0*pi*ntm(1)
      k(2) = 2.0d0*pi*ntm(2)
      k(3) = 2.0d0*pi*ntm(3)

C Loop over all atoms predict the egg box effect and subtract it off
      E=0
      do ia=1,na
C i is the atom type
         i=isa(ia)
         x=xa_frac(1,ia)
         y=xa_frac(2,ia)
         z=xa_frac(3,ia)

C The predicted energy
         do a=0,maxshell
           do b=0,maxshell
             do c=0,maxshell
               E=E+vreal(a,b,c,i)*cos(a*k(1)*x+b*k(2)*y+c*k(3)*z)
             enddo
           enddo
         enddo

C The predicted forces
         f(1,ia)=0.0_dp
         f(2,ia)=0.0_dp
         f(3,ia)=0.0_dp
         do a=0,maxshell
            do b=0,maxshell
               do c=0,maxshell

                  f(1,ia)=f(1,ia)-a*k(1)*(
     .               vreal(a,b,c,i)*sin(a*k(1)*x+b*k(2)*y+c*k(3)*z))

                  f(2,ia)=f(2,ia)-b*k(2)*(
     .               vreal(a,b,c,i)*sin(a*k(1)*x+b*k(2)*y+c*k(3)*z))


                  f(3,ia)=f(3,ia)-c*k(3)*(
     .               vreal(a,b,c,i)*sin(a*k(1)*x+b*k(2)*y+c*k(3)*z))

               enddo
            enddo
         enddo

C Convert forces to cartesian coordinates
         do ix = 1,3
            fxyz(ix,ia) = celli(ix,1) * f(1,ia) +
     .           celli(ix,2) * f(2,ia) +
     .           celli(ix,3) * f(3,ia)
         enddo
         fxyz(1,ia)=-fxyz(1,ia)*escale
         fxyz(2,ia)=-fxyz(2,ia)*escale
         fxyz(3,ia)=-fxyz(3,ia)*escale
            
C Finally remove the eggbox forces from fa
         if (task.eq.'forces') then
            fa(1,ia)=(fa(1,ia)-fxyz(1,ia))
            fa(2,ia)=(fa(2,ia)-fxyz(2,ia))
            fa(3,ia)=(fa(3,ia)-fxyz(3,ia))
         endif         
      enddo  

C Remove the eggbox energy
      if (task.eq.'energy') then
         Etot=Etot-E*escale
      endif
         
      end subroutine eggbox
