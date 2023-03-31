#include "settings.inc"
#if DIM <= 1
subroutine ROUTINE2(put_gatt,VAR)(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  VAR_TYPE, intent(in) :: att DIMS
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_put_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '/ /trim(name)/ /' in file: '/ /this)
end subroutine ROUTINE2(put_gatt,VAR)

subroutine ROUTINE2(put_att,VAR)(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  VAR_TYPE, intent(in) :: att DIMS
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  if (this%define > -1 ) call ncdf_redef(this)
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '/ /trim(var)/ /' : '/ /trim(name)/ /' in file: '/ /this)
  call ncdf_err(nf90_put_att(this%id, id, name, att), &
       'Saving (VAR) '/ /trim(var)/ /' attribute: '/ /trim(name)/ /' in file: '/ /this)
end subroutine ROUTINE2(put_att,VAR)

subroutine ROUTINE2(get_att,VAR)(this, var, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: var, name
  VAR_TYPE, intent(out) :: att DIMS
  integer :: id
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_inq_varid(this%id, trim(var), id), &
              'Retrieving id from (VAR) '/ /trim(var)/ /' : '/ /trim(name)/ /' in file: '/ /this)
  call ncdf_err(nf90_get_att(this%id, id, name, att), &
       'Retrieving (VAR) '/ /trim(var)/ /' attribute: '/ /trim(name)/ /' in file: '/ /this)
end subroutine ROUTINE2(get_att,VAR)

subroutine ROUTINE2(get_gatt,VAR)(this, name, att)
  type(hNCDF), intent(inout) :: this
  character(len=*), intent(in) :: name
  VAR_TYPE, intent(out) :: att DIMS
  if ( .not. ncdf_participate(this) ) return
  call ncdf_err(nf90_get_att(this%id, NF90_GLOBAL, name, att), &
       'Saving global (VAR) attribute: '/ /trim(name)/ /' in file: '/ /this)
end subroutine ROUTINE2(get_gatt,VAR)

#endif