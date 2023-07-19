MODULE m_ts_kpoints
!
! Routines that are related to TS kpoint sampling
!
!==============================================================================
! CONTAINS:
!          1) setup_ts_scf_kscell
!          2) ts_write_k_points


  USE precision, only : dp
 
  implicit none

  private


!===================== K-POINT RELATED VARIABLES ============================== 
!==============================================================================
! Contains data structures and routines to deal with the kpoint-grid
! for the self-consistent calculation with Green's Functions.
! Original verison modified so that the kpoint sampling is done for the
! direction perpendicular to the transport directins
! Version created to integrate TRANSIESTA in siesta2.3
!==============================================================================
! The kpoint mesh parameters that are public may be accessed in other 
! parts of the code. With respesct to the original names of the variables used 
! in the m_kpoint_grid module, a "ts_" was added. Also in the SIESTA module
! kscell and kdispl are not public, but are made public (ts_kscell and
! ts_kdispl) here.

!==============================================================================
! DETAILS: To obatin the kpoints for the GFs calculations it uses the same 
! scheme as for SIESTA but puts ts_kscell(3,3)=1 and ts_kdispl(3)=0.0
!
!==============================================================================
! OBS: So far if only a kgrid cutoff is specified the code stops. TO be 
!      changed in the near future
!==============================================================================
 
  logical, public, save   :: ts_scf_kgrid_first_time = .true.
  logical, public, save   :: ts_gamma_scf
  integer, public, save   :: ts_maxk              ! 
  integer, public, save   :: ts_nkpnt             ! Total number of k-points
  real(dp), public, save  :: ts_eff_kgrid_cutoff  ! Effective kgrid_cutoff

  real(dp), pointer, public, save :: ts_kweight(:) 
  real(dp), pointer, public, save :: ts_kpoint(:,:)

  integer,  public, dimension(3,3), save  :: ts_kscell = 0
  real(dp), public, dimension(3), save    :: ts_kdispl = 0.0_dp

  logical, public, save     :: ts_user_requested_mp = .false.
  logical, public, save     :: ts_user_requested_cutoff = .false.

  logical, public, save     :: ts_spiral = .false.
  logical, public, save     :: ts_firm_displ = .false.

  public :: setup_ts_scf_kscell, ts_write_k_points

  CONTAINS


!-----------------------------------------------------------------------

  subroutine setup_ts_scf_kscell( cell, firm_displ )

! ***************** INPUT **********************************************
! real*8  cell(3,3)  : Unit cell vectors in real space cell(ixyz,ivec)
! ***************** OUTPUT *********************************************
! logical firm_displ   : User-specified displacements (firm)?

!   The relevant fdf labels are kgrid_cutoff and kgrid_Monkhorst_Pack.
!   If both are present, kgrid_Monkhorst_Pack has priority. If none is
!   present, the cutoff default is zero, producing only the gamma point.
!   Examples of fdf data specifications:
!     kgrid_cutoff  50. Bohr
!     %block kgrid_Monkhorst_Pack  # Defines kscell and kdispl
!     4  0  0   0.50               # (kscell(i,1),i=1,3), kdispl(1)
!     0  4  0   0.50               # (kscell(i,2),i=1,3), kdispl(2)
!     0  0  4   0.50               # (kscell(i,3),i=1,3), kdispl(3)
!     %endblock kgrid_Monkhorst_Pack
! **********************************************************************

!  Modules

      use precision,  only : dp
      use parallel,   only : Node
      use m_minvec,   only : minvec
      use fdf
      use sys,        only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

! Passed variables
      real(dp), intent(in)   :: cell(3,3)
      logical, intent(out)   :: firm_displ

! Internal variables
      integer           i, iu, j,  factor(3,3), expansion_factor
#ifdef MPI
      integer           MPIerror
#endif
      real(dp)          scmin(3,3),  vmod, cutoff
      real(dp)          ctransf(3,3)
      logical           mp_input

      real(dp), parameter :: defcut = 0.0_dp
      integer, dimension(3,3), parameter :: unit_matrix =  &
                         reshape ((/1,0,0,0,1,0,0,0,1/), (/3,3/))

      if (Node.eq.0) then

         mp_input = fdf_block('kgrid_Monkhorst_Pack',iu)     
         if ( mp_input ) then
          
            ts_user_requested_mp = .true.
            do i = 1,3
               read(iu,*) (ts_kscell(j,i),j=1,3), ts_kdispl(i)
            enddo
            firm_displ = .true.

         else

            write(*,*) 'WARNING !!!'
            write(*,*) 'TS kgrid determined first with 3D cell !!!'
            write(*,*) 'Specifying only cutoff in Electrode AND Scattering calculations might lead to problems !!'

            cutoff = fdf_physical('kgrid_cutoff',defcut,'Bohr')
            if (cutoff /= defcut) then
               ts_user_requested_cutoff = .true.
            endif

            ts_kdispl(1:3) = 0.0_dp  ! Might be changed later
            
            ! Find equivalent rounded unit-cell
            call minvec( cell, scmin, ctransf )

            expansion_factor = 1
            do j = 1,3
               factor(j,1:3) = 0
               vmod = sqrt(dot_product(scmin(1:3,j),scmin(1:3,j)))
               factor(j,j) = int(2.0_dp*cutoff/vmod) + 1
               expansion_factor = expansion_factor * factor(j,j)
            enddo
            ! Generate actual supercell skeleton
            ts_kscell = matmul(ctransf, factor)
            ! Avoid confusing permutations
            ! Defer implementation to avoid diffs with reference version
            !!if (expansion_factor == 1) then
            !!   kscell = unit_matrix
            !!endif
         endif
      endif

! Modify the ts_kscell and ts_kdispl to obtain the 2D sampling
    ts_kscell(1:3,3)=0
    ts_kscell(3,1:3)=0
    ts_kscell(3,3)=1
    ts_kdispl(3)=0.0

#ifdef MPI
      call MPI_Bcast(ts_kscell(1,1),9,MPI_integer,0,MPI_Comm_World, MPIerror)
      call MPI_Bcast(ts_kdispl,3,MPI_double_precision,0,MPI_Comm_World, MPIerror)
      call MPI_Bcast(ts_firm_displ,1,MPI_logical,0,MPI_Comm_World, MPIerror)
      call MPI_Bcast(ts_user_requested_mp,1,MPI_logical,0,   &
                     MPI_Comm_World, MPIerror)
      call MPI_Bcast(ts_user_requested_cutoff,1,MPI_logical,0,   &
                     MPI_Comm_World, MPIerror)
#endif
    

    end subroutine setup_ts_scf_kscell

    subroutine ts_write_k_points()
      USE siesta_options, only: writek
      USE units, only: Ang
      USE siesta_cml

      implicit none

      integer  :: ik, ix, i
      external :: iokp

      if ( writek ) then
         write(6,'(/,a)') 'transiesta: ts_k-point coordinates (Bohr**-1) and weights:'
         write(6,'(a,i4,3f12.6,3x,f12.6)')                          &
              ('transiesta: ', ik, (ts_kpoint(ix,ik),ix=1,3), ts_kweight(ik), &
              ik=1,ts_nkpnt)
      else
         call iokp( ts_nkpnt, ts_kpoint, ts_kweight )
      endif
      write(6,'(/a,i6)')  'transiesta: ts_k-grid: Number of Transport k-points =', ts_nkpnt
      write(6,'(a)') 'transiesta: ts_k-grid: Supercell and displacements'
      write(6,'(a,3i4,3x,f8.3)') 'transiesta: ts_k-grid: ',        &
           (ts_kscell(i,1),i=1,3), ts_kdispl(1)
      write(6,'(a,3i4,3x,f8.3)') 'transiesta: ts_k-grid: ',        &
           (ts_kscell(i,2),i=1,3), ts_kdispl(2)
!      write(6,'(a,3i4,3x,f8.3)') 'transiesta: ts_k-grid: ',        &
!           (ts_kscell(i,3),i=1,3), ts_kdispl(3)
      write(6,*)
      if (cml_p) then
          call cmlStartPropertyList(xf=mainXML, title="Transiesta k-points", &
                 dictRef="siesta:ts_kpoints")
          call cmlAddProperty(xf=mainXML, value=ts_nkpnt, dictref='siesta:ts_nkpnt', &
                 units="cmlUnits:countable")
          do ik=1, ts_nkpnt
              call cmlAddKPoint(xf=mainXML, coords=ts_kpoint(:,ik), weight=ts_kweight(ik))
          enddo
          call cmlEndPropertyList(mainXML)
          call cmlAddProperty(xf=mainXML, value=ts_kscell, dictref='siesta:ts_kscell', &
                 units="siestaUnits:Ang")
          call cmlAddProperty(xf=mainXML, value=ts_kdispl, dictref='siesta:ts_kdispl', &
                 units="siestaUnits:Ang")
      endif
    end subroutine ts_write_k_points

END MODULE m_ts_kpoints
