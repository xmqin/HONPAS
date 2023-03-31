 subroutine extrapolatedSpData2D(SpM1,SpM2,sp_out,SpMout)

   use class_dSpData2D
   use class_dData2D
   use class_Sparsity
   use class_OrbitalDistribution

   use m_restruct_SpData2D, only: restruct_dSpData2D

   type(dSpData2D), intent(in)    :: SpM1, SpM2
   type(Sparsity), intent(in)    :: sp_out
   ! Note!!  inout is essential to avoid memory leaks...
   type(dSpData2D), intent(inout) :: SpMout

   integer, parameter :: dp = selected_real_kind(10,100)

   type(dSpData2D)     :: SpM1_out, SpM2_out
   type(dData2D)     :: a2d_out

   integer                              :: dim2
   real(dp), dimension(:,:), pointer    :: a_1, a_2, a_out

   call restruct_dSpData2D(SpM1,sp_out,SpM1_out)
   call restruct_dSpData2D(SpM2,sp_out,SpM2_out)

   a_1  => val(SpM1_out)
   a_2  => val(SpM2_out)
   dim2 = size(a_1, dim=2)

   call newdData2D(a2d_out,nnzs(sp_out),dim2,"(new in extrapolate)")
   a_out => val(a2d_out)

   ! Simple linear extrapolation for now
   a_out(:,:) = 2.0_dp * a_2(:,:) - a_1(:,:)

   call newdSpData2D(sp_out,a2d_out,dist(SpM1), &
        SpMout,name="Extrapolated SpM")

   call delete(a2d_out)
   call delete(SpM1_out)
   call delete(SpM2_out)

 CONTAINS
        subroutine die(str)
          character(len=*), optional :: str
          if (present(str)) then
             print *, trim(str)
          endif
          stop
        end subroutine die

 end subroutine extrapolatedSpData2D
