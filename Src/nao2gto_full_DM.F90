      module get_DM

      use precision,     only : dp
      use atm_types,     only : nso, nco
      use atmfuncs,      only : lofio, mofio
      use atomlist,      only : indxuo
      use listsc_module, only : listsc
      use hfx_types 
      implicit none

      private
      public :: sparse_dense, get_pmax_shell, dense2sparse, get_um_cut, &
                build_indx

      contains

!-----------------------------------------------------------------------
      subroutine sparse_dense( nspin, nuo, nuotot, norb, maxnd, numd, &
                               listdptr, listd, M_sparse, M_dense )
!-----------------------------------------------------------------------
!  Parallel treatment:
!  nuo      : Local number of orbitals within the unitcell local in this node
!             nuo <= nuotot
!  norb     : Global number of orbitals within the supercell
!             norb = nuotot * ncells
!  maxnd    : Maximum size of Sparse matrix
!  numd(nuo)       : Number of nonzero elements of each row
!  listd(maxnd)    : Nonzero hamiltonian-matrix element
!  listdptr(nuo)   : Control vector of density matrix
!  M_sparse(mxnd)  : sparse maxtrix

! Input variables-------------------------------------------------------
      integer  nspin, nuo, nuotot, norb, maxnd, numd(nuo),  &
               listdptr(nuo), listd(*)

      real(dp)  M_sparse(maxnd, nspin)
! Output variable-------------------------------------------------------
      real(dp)  M_dense(norb,norb,nspin)

! Tempos, internal variables -------------------------------------------
      integer  ispin, io, iio, jo, i, j, iu, ju, ind
      real(dp), dimension(:,:), allocatable  ::  M_aux
!-----------------------------------------------------------------------

      allocate(M_aux(nuotot,norb))

      M_aux(1:nuotot,1:norb) = 0.0d0
      M_dense(1:norb,1:norb,1:nspin) = 0.0d0

      do ispin=1,nspin
         do io = 1,nuotot
            do j = 1,numd(io)
               ind = listdptr(io) + j
               jo = listd(ind)
               M_aux(io,jo) = M_sparse(ind,ispin)
               M_dense(io,jo,ispin)=M_aux(io,jo)
            enddo
         enddo

!         if(nuotot.eq.norb)  then
!         return
!         else
         do io=nuotot+1,norb
            iu=indxuo(io)
            do j=1,numd(iu)
               ind=listdptr(iu) +j
               ju=listd(ind)
               jo=listsc(io,iu,ju)
               M_dense(io,jo,ispin)=M_aux(iu,ju)
            enddo
         enddo
!        endif
      enddo

      deallocate(M_aux)

      end subroutine sparse_dense

!-----------------------------------------------------------------------
      subroutine get_pmax_shell( nspin, norb, iaorb, iphorb, nuo, na, &
                                 isa, maxnd, numd, listdptr, listd, &
                                 Dscf, Dscf_max )

! Input variables-------------------------------------------------------
      integer nspin, norb, nuo, iaorb(norb), iphorb(norb), na, &
              isa(na), maxnd, numd(nuo), listdptr(nuo), listd(*)

      real(dp)  Dscf(maxnd, nspin)
! Output variable-------------------------------------------------------
      real(dp)  Dscf_max(maxnd)

! Tempos, internal variables -------------------------------------------
      integer  ispin, ia, ja, io, ioa, is, jo, joa, js, i, j, iu, ju, ind
      integer il, im, jl, jm
      integer :: io_in_shell, jo_in_shell, ind_in_shell
      logical :: find_in_Dscf

!-----------------------------------------------------------------------
!        call memory('A','L',norb*norb,'hfx_potential')

      do ispin = 1, nspin
      do io  = 1, norb
         ia = iaorb(io)
         ioa = iphorb(io)
         is  = isa(ia)
         il  = lofio(is,ioa)
         im  = mofio(is,ioa)
         if(im .ne. -il) cycle

         iu = indxuo(io)
         do j=1, numd(iu)
            ind = listdptr(iu) +j
            ju  = listd(ind)
            jo  = listsc(io,iu,ju)
            ja = iaorb(jo)
            joa = iphorb(jo)
            js = isa(ja)
            jl  = lofio(js,joa)
            jm  = mofio(js,joa)
            if(jm .ne. -jl) cycle

                do io_in_shell = io, io-1+nso(il)
                 do jo_in_shell = jo, jo-1+nso(jl)
                    ind_in_shell = D2Sindx(io_in_shell,jo_in_shell)
                    if((ind_in_shell.ne.0) .and. (abs(Dscf(ind_in_shell,ispin)) .gt. Dscf_max(ind)) ) then
                       Dscf_max(ind)= abs(Dscf(ind_in_shell,ispin))
                    endif
                enddo
               enddo

         enddo
       enddo
      enddo

      end subroutine get_pmax_shell

      subroutine build_indx(nuotot, norb, maxnd, numd, listdptr,listd,D2Sindx)

       implicit none

       integer , intent(in)  :: nuotot,norb, maxnd
       integer , intent(in)  :: numd(nuotot), listdptr(nuotot),listd(maxnd)
       integer, intent(inout) :: D2Sindx(norb,norb)
      ! Local variables
       integer :: io, jo, j, iu, ju, ind
      
!      write(99,*) maxnd
       do io = 1, norb
          iu = indxuo(io)
!          write(99,*) "io for ",io, iu
        do j=1,numd(iu)
           ind=listdptr(iu) + j
           ju=listd(ind)
           jo=listsc(io,iu,ju)
           D2Sindx(io,jo) = ind
!           write(99,*) "j ",j, jo, ju, ind
!            write(100,*) io, jo, ind
         enddo
      enddo

!      do  io = 1, norb
!         do jo = 1, norb
!            write(31,*) io, indxuo(io),jo, D2Sindx(io,jo)
!         enddo
!      enddo

       end subroutine build_indx


      subroutine dense2sparse(io_input, jo_input, ind_output, find_in_Dscf, &
                             nuotot, norb, maxnd, numd, listdptr,listd)

       use precision, only: dp

       implicit none

       integer , intent(in)  :: io_input, jo_input
       integer , intent(out) :: ind_output
       logical , intent(out) :: find_in_Dscf
       integer , intent(in)  :: nuotot, norb, maxnd
       integer , intent(in)  :: numd(nuotot), listdptr(nuotot), listd(maxnd)

      ! Local variables
       integer :: io, jo, j, iu, ju, ind

        io = io_input
        iu = indxuo(io)
        do j=1,numd(iu)
           ind=listdptr(iu) + j
           ju=listd(ind)
           jo=listsc(io,iu,ju)
           if(jo.eq.jo_input)  then
             ind_output = ind
             find_in_Dscf = .true.
             return
           endif
        enddo
        ind_output = 0
        find_in_Dscf = .false.

       end subroutine dense2sparse


      subroutine get_um_cut(io_input, jo_input, um_cut,  &
                          nuotot, norb, maxnd, numd, listdptr,listd)

      use precision, only: dp

      implicit none

      integer , intent(in)  :: io_input, jo_input
      logical , intent(out) :: um_cut
      integer , intent(in)  :: nuotot, norb, maxnd
      integer , intent(in)  :: numd(nuotot), listdptr(nuotot), listd(maxnd)

      ! Local variables
      integer :: io, jo, j, iu, ju, ind

        um_cut = .true.
        io = io_input
        iu = indxuo(io)
        do j=1,numd(iu)
           ind=listdptr(iu) + j
           ju=listd(ind)
           jo=listsc(io,iu,ju)
           if(jo.eq.jo_input)  then
             um_cut = .false.
             return
           endif
        enddo

      end subroutine get_um_cut

      end module get_DM
