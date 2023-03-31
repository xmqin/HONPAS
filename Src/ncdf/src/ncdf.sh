#!/bin/bash

_vpath=.
if [ ! -z "$VPATH" ]; then
    _vpath=$VPATH
fi

if [ -z "$DIR_FDICT" ]; then
    var_dir=$_vpath/../fdict
else
    var_dir=$DIR_FDICT
fi
source $var_dir/default_settings.bash
[ -e $var_dir/settings.sh ] && source $var_dir/settings.bash

# Override any special settings in this file
[ -e $_vpath/settings.bash ] && source $_vpath/settings.bash

# The different variable types used in this (long does not exist)
vars=(h s d c z i)

function c_to_r {
    local var=$1 ; shift
    case $var in
	c) _ps s ;;
	z) _ps d ;;
    esac
}
	   
function has_att {
    local var=$1 ; shift
    case $var in
	h|c|z|l) _ps 0 ;;
	*) _ps 1 ;;
    esac
}

# Create the interface files
{
# Variable creation
for sub in put get ; do
_psnl "interface ncdf_${sub}_var"
for v in ${vars[@]} ; do
    for d in `seq 0 $(var_N $v)` ; do 
	_psnl "module procedure ${sub}_var_${v}${d}_name"
    done
done
_psnl "end interface ncdf_${sub}_var"
done
# Global/local attribute set
for ssub in put get ; do
for sub in ${ssub}_gatt ${ssub}_att ; do
_psnl "interface ncdf_${sub}"
_psnl "module procedure ${sub}"
for v in ${vars[@]} ; do
    [ $(has_att $v) -eq 0 ] && continue
    for d in `seq 0 $(var_N $v)` ; do 
	# Attributes does not allow dimensions larger than 1
	[ $d -gt 1 ] && continue
	_psnl "module procedure ${sub}_${v}${d}"
    done
done
_psnl "end interface ncdf_${sub}"
done
done
# variable and defining fill-values
for sub in inq_var def_fill ; do
_psnl "interface ncdf_$sub"
[ "$sub" == "inq_var" ] && _psnl "module procedure ncdf_${sub}_def"
for v in ${vars[@]} ; do
    d=0
    _psnl "module procedure ${sub}_${v}${d}"
done
_psnl "end interface ncdf_$sub"
done
} > netcdf_ncdf_interface_.inc

{
for v in ${vars[@]} ; do
    _psnl "#define VAR_TYPE $(var_name $v)"
    for d in `seq 0 $(var_N $v)` ; do
	if [ $d -eq 0 ]; then
	    _psnl "#define DIMS"
	else
	    _psnl "#define DIMS $(dim_to_size $d)"
	fi
	_psnl "#define VAR $v$d"
	_psnl "#define DIM $d"
	case $v in
	    c)
		_psnl "#define COMPLEX_DTYPE sp"
		_psnl "#define REAL_TYPE $(var_name $(c_to_r $v))"
		;;
	    z)
		_psnl "#define COMPLEX_DTYPE dp"
		_psnl "#define REAL_TYPE $(var_name $(c_to_r $v))"
		;;
	    *)
		_psnl "#define REAL_TYPE VAR_TYPE"
		;;
	esac
        # Attributes only allowed for dimensions larger than 1
	if [ $d -le 1 ] && [ $(has_att $v) -eq 1 ]; then
	    _psnl '#include "netcdf_ncdf_att_inc.inc"'
	fi
	_psnl '#include "netcdf_ncdf_var_inc.inc"'
	_psnl "#undef COMPLEX_DTYPE"
	_psnl "#undef REAL_TYPE"
	_psnl "#undef VAR"
	_psnl "#undef DIM"
	_psnl "#undef DIMS"
    done
    _psnl "#undef VAR_TYPE"
done
} > netcdf_ncdf_funcs_.inc
