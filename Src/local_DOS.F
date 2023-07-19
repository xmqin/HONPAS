
      MODULE m_local_DOS
      private
      public :: local_DOS

      CONTAINS

      subroutine local_DOS( )
      use m_energies

      use sparse_matrices
      USE siesta_options, only :occtol, g2cut, isolve, 
     &                          fixspin, temp, neigwanted 
      use siesta_geom
      use atomlist,       only: indxuo, indxua           
      use atomlist,       only: qtot, no_u, no_l
      use atomlist,       only: iphorb                   
      use atomlist,       only: datm, no_s, iaorb        
      use fdf,            only : fdf_block, fdf_convfac  
      use sys,            only: die                      
      use m_mpi_utils,    only : broadcast               
      use files,          only : slabel                  
      use Kpoint_grid
      use parallel,       only: IOnode                   
      use files,          only : label_length            
      use parse
      use m_ntm
      use m_forces,       only: fa
      use m_eo
      use m_spin,         only: nspin, qs, efs
      use m_gamma
      implicit none

      character(len=label_length+5), external :: paste

      type(parsed_line), pointer  :: p


      integer :: dummy_iscf = 1
      character(len=label_length+5) :: fildos

      real(dp):: e1  ! Lower bound of energy range
      real(dp):: e2  ! Upper bound of energy range

      real(dp)  :: dummy_str(3,3), dummy_strl(3,3)  ! for dhscf call
      real(dp)  :: dummy_dipol(3)

      integer   :: iu
      real(dp)  :: factor, g2max
      logical   :: genlogic   ! Generate local density of states?
      character(len=150) :: line    ! Data line read

! Find local density of states
      if (IOnode) then
        genlogic = fdf_block('LocalDensityOfStates',iu)
      endif
      call broadcast(genlogic)

      if ( genlogic ) then

! Find the desired energy range
        if (IOnode) then
          read(iu,'(a)') line
          p=>digest(line)
          if (.not. match(p,"vvn"))
     .       call die("Wrong format in LocalDensityofStates")
          factor = fdf_convfac( names(p,1), 'Ry' )
          e1 = values(p,1)*factor
          e2 = values(p,2)*factor
          call destroy(p)
        endif
        call broadcast(e1)
        call broadcast(e2)

        !Find the density matrix for states between e1 and e2
        if (isolve .eq. 0) then
          call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                numh, listhptr, listh, numh, listhptr, listh, 
     .                H, S, qtot, fixspin, qs, temp, e1, e2,
     .                gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                occtol, dummy_iscf, neigwanted)

          !Find the LDOS in the real space mesh
          fildos = paste( slabel, '.LDOS' )
          g2max = g2cut
          call dhscf( nspin, no_s, iaorb, iphorb, no_l,
     .              no_u, na_u, na_s, isa, xa_last, indxua, 
     .              ucell_last, mscell, g2max, ntm,
     .              0, 0, 0, fildos, ' ', ' ', ' ', ' ', ' ',
     .              maxnh, numh, listhptr, listh, Dscf, Datm,
     .              maxnh, numh, listhptr, listh, H,
     .              Enaatm, Enascf, Uatm, Uscf, DUscf, DUext, Exc, Dxc,
     .              dummy_dipol, fa, dummy_str, fa, dummy_strl )
                    ! next to last argument is dummy here,
                    ! as no forces are calculated
                    ! todo: make all these optional
        else
          if (IOnode)  write(6,*)
     .       'siesta: ERROR: LDOS implemented only with diagon'
        endif

      endif ! genlogic


      END subroutine local_DOS

      END module m_local_DOS
