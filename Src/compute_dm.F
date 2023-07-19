      MODULE m_compute_dm
      private
      public :: compute_dm

      CONTAINS

      subroutine compute_dm( iscf )
      use precision
      use units, only: eV
      USE siesta_options
      use sparse_matrices
      use siesta_geom
      use atomlist, only: qa, lasto, iphorb, iaorb, no_u, no_s, indxuo, 
     .                    qtot, no_l
      use sys, only: die
      use m_ordern,     only: ordern
      use Kpoint_grid
      use m_energies,   only: Ebs, DEHarr, Ecorrec, Ef, Entrop
      use m_steps, only: istp
      use m_rmaxh
      use m_eo
      use m_spin,  only: nspin, qs, efs
      use m_gamma
      use parallel, only: IONode
#ifdef MPI
      use m_mpi_utils, only: globalize_sum
#endif
#ifdef CDF
      use iodmhs_netcdf, only: write_dmh_netcdf
#endif
#ifdef TRANSIESTA
      use m_ts_options, only : mixH
      use m_ts_global_vars, only : TSiscf, TSinit, TSrun, Hold
      use m_ts_in_siesta, only : transiesta
#endif /* TRANSIESTA */


      implicit none

      integer, intent(in) :: iscf

      integer :: io, ispin
      real(dp):: const  ! Temporary factor used in calculation of Harris forces
      ! e1>e2 To signal that we do not want DOS weights  
      real(dp), parameter  :: e1 = 1.0_dp, e2 = -1.0_dp

#ifdef MPI
      real(dp):: buffer1   ! Temporary for the result of a globalize operation
#endif

!------------------------------------------------------------------------- BEGIN

! Solve eigenvalue problem 
          if (isolve .eq. 0) then
            call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                  numh, listhptr, listh, numh, listhptr, listh, 
     .                  H, S, qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                  eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                  occtol, iscf, neigwanted)
            Ecorrec = 0.0_dp
!
          elseif (isolve .eq. 1) then
            if (.not. gamma) call die("Cannot do O(N) with k-points.")
            call ordern(usesavelwf,ioptlwf,na_u,no_u,no_l,lasto,
     .                  isa,qa,rcoor,rmaxh,ucell,xa,iscf,istp,ncgmax,
     .                  etol,eta,qtot,maxnh,numh,listhptr,listh,H,S,
     .                  chebef,noeta,rcoorcp,beta,pmax,Dscf,Escf,
     .                  Ecorrec,nspin,qs)
            Entrop = 0.0_dp
#ifdef TRANSIESTA
! TSS Begin
          elseif (isolve .eq. 2 .and. TSinit) then
            call diagon(no_s, nspin, nspin, no_l, maxnh, maxnh, no_u,
     .                  numh, listhptr, listh, numh, listhptr, listh,
     .                  H, S, qtot, fixspin, qs, temp, e1, e2,
     .                  gamma, xijo, indxuo, nkpnt, kpoint, kweight,
     .                  eo, qo, Dscf, Escf, ef, efs, Entrop, no_u,
     .                  occtol, iscf, neigwanted)

            Ecorrec = 0._dp

         elseif (TSrun) then
            call timer('TS',1)
            TSiscf=TSiscf+1

            call transiesta(gamma,lasto,
     .        no_l, no_s, nspin, na_u, maxnh, numh, listhptr,
     .        listh, xijo, indxuo, H, S, Dscf, Escf, ef, no_u,
     .        TSiscf,qtot)
            call timer('TS',2)

            Ecorrec = 0._dp
! TSS End
#endif /* TRANSIESTA */
          else
            call die('siesta: ERROR: wrong solution method')
          endif

#ifdef TRANSIESTA
! Save present H matrix
          if(mixH) then
             do ispin = 1,nspin
                do io = 1,nh
                   Hold(io,ispin) = H(io,ispin)
                enddo
             enddo
          end if
#endif /* TRANSIESTA */

! Harris-functional energy 
          DEharr = 0.0_dp
          Ebs = 0.0_dp
          do ispin = 1,nspin
! const factor takes into account that there are two nondiagonal
! elements in non-collinear spin density matrix, stored as
! ispin=1 => D11; ispin=2 => D22, ispin=3 => Real(D12);
! ispin=4 => Imag(D12)
            const = 1._dp
            if (ispin .gt. 2) const = 2._dp
            do io = 1,nh
              DEharr = DEharr + H(io,ispin) * const * 
     .                     ( Dscf(io,ispin) - Dold(io,ispin) )
              Ebs = Ebs + H(io,ispin) * const * Dscf(io,ispin)
            enddo
          enddo
#ifdef MPI
! Global reduction of DEharr
          call globalize_sum(DEharr,buffer1)
          DEharr = buffer1
          call globalize_sum(Ebs,buffer1)
          Ebs = buffer1
#endif
! ..................
#ifdef CDF
          if ( writedmhs_cdf_history) then
             call write_dmh_netcdf( no_l, nh, nspin,
     $                              Dold, H, Dscf)
          else if ( writedmhs_cdf) then
             call write_dmh_netcdf( no_l, nh, nspin,
     $                              Dold, H, Dscf, overwrite=.true.)
          endif
#endif

! Print populations at each SCF step if requested before mixing ......

          if (muldeb ) then
             if (ionode) write (6,"(/a)")
     .             'siesta: Mulliken populations before mixing'
             call mulliken( mullipop, nspin, na_u, no_u, maxnh,
     .                      numh, listhptr, listh, S, Dscf, isa,
     .                      lasto, iaorb, iphorb )
          endif
!--------------------------------------------------------------------------- END
      END subroutine compute_dm
      END MODULE m_compute_dm
