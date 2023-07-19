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
      subroutine setup_kscell( cell, kscell, displ, firm_displ )

c ***************** INPUT **********************************************
c real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
c ***************** OUTPUT *********************************************
c integer kscell(3,3)  : Supercell reciprocal of k-grid unit cell
c real*8  displ(3)     : Grid origin in k-grid-vector coordinates:
c logical firm_displ   : User-specified displacements

c   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
c   If both are present, kgrid_Monkhorst_Pack has priority. If none is
c   present, the cutoff default is zero, producing only the gamma point.
c   Examples of fdf data specifications:
c     kgrid_cutoff  50. Bohr
c     %block kgrid_Monkhorst_Pack  # Defines kscell and displ
c     4  0  0   0.50               # (kscell(i,1),i=1,3), displ(1)
c     0  4  0   0.50               # (kscell(i,2),i=1,3), displ(2)
c     0  0  4   0.50               # (kscell(i,3),i=1,3), displ(3)
c     %endblock kgrid_Monkhorst_Pack
c **********************************************************************
C
C  Modules
C
      use precision,  only : dp
      use units,      only : Ang
      use parallel,   only : Node
      use m_minvec,   only : minvec
      use fdf
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

C Passed variables
      real(dp), intent(in)   :: cell(3,3)
      integer, intent(out)   :: kscell(3,3)
      real(dp), intent(out)  :: displ(3)
      logical, intent(out)   :: firm_displ

C Internal variables
      integer           i, iu, j,  factor(3,3), expansion_factor
#ifdef MPI
      integer           MPIerror
#endif
      real(dp)          scmin(3,3),  vmod
      real(dp)          ctransf(3,3)
      logical           genlogic

      real(dp), parameter :: defcut = 0.0_dp
      integer, dimension(3,3), parameter :: unit_matrix = 
     .                    reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

      if (Node.eq.0) then

         genlogic = fdf_block('kgrid_Monkhorst_Pack',iu)
         if ( genlogic ) then
            do i = 1,3
               read(iu,*) (kscell(j,i),j=1,3), displ(i)
            enddo
            firm_displ = .true.

         else

            cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
            write(6,"(a,f10.5)") "Kgrid cutoff input: ", cutoff
            firm_displ = .false.
            
            ! Find equivalent rounded unit-cell
            call minvec( cell, scmin, ctransf )

            expansion_factor = 1
            do j = 1,3
               factor(j,1:3) = 0
               vmod = sqrt(dot_product(scmin(1:3,j),scmin(1:3,j)))
               factor(j,j) = int(2.0_dp*cutoff/vmod) + 1
               if (mod(factor(j,j),2) == 0) then
                  displ(j) = 0.5_dp
               else
                  displ(j) = 0.0_dp
               endif
               expansion_factor = expansion_factor * factor(j,j)
            enddo
            ! Generate actual supercell skeleton
            kscell = matmul(ctransf, factor)
            ! Avoid confusing permutations
            if (expansion_factor == 1) then
               kscell = unit_matrix
            endif
         endif
      endif
C
#ifdef MPI
      call MPI_Bcast(kscell(1,1),9,MPI_integer,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(displ,3,MPI_double_precision,0,MPI_Comm_World,
     .     MPIerror)
      call MPI_Bcast(firm_displ,1,MPI_logical,0,MPI_Comm_World,
     .     MPIerror)
#endif

      end subroutine setup_kscell

