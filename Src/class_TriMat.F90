! We define all variables in this one
! As the entries are so short this makes more sense

! Currently there is no reason to generate all the modules.
! They are not used...
#ifdef __CORRECT_FOR_TRIMAT_USAGE
module class_lTriMat
!========================
#define TYPE_NAME lTriMat
#define STR_TYPE_NAME "lTriMat"
#define TYPE_NAME_ lTriMat_
#define NEW_TYPE newlTriMat
#define VAR_TYPE logical
! DO NOT DEFINE PREC
#define VAR_INIT .false.
#include "class_TriMat.T90"
!========================
end module class_lTriMat

module class_iTriMat
!========================
#define TYPE_NAME iTriMat
#define STR_TYPE_NAME "iTriMat"
#define TYPE_NAME_ iTriMat_
#define NEW_TYPE newiTriMat
#define VAR_TYPE integer
! DO NOT DEFINE PREC
#define VAR_INIT 0
#include "class_TriMat.T90"
!========================
end module class_iTriMat


module class_sTriMat
!========================
#define TYPE_NAME sTriMat
#define STR_TYPE_NAME "sTriMat"
#define TYPE_NAME_ sTriMat_
#define NEW_TYPE newsTriMat
#define VAR_TYPE real
#define PREC sp
#define VAR_INIT 0._sp
#include "class_TriMat.T90"
!========================
end module class_sTriMat

module class_dTriMat
!========================
#define TYPE_NAME dTriMat
#define STR_TYPE_NAME "dTriMat"
#define TYPE_NAME_ dTriMat_
#define NEW_TYPE newdTriMat
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_TriMat.T90"
!========================
end module class_dTriMat

module class_cTriMat
!========================
#define TYPE_NAME cTriMat
#define STR_TYPE_NAME "cTriMat"
#define TYPE_NAME_ cTriMat_
#define NEW_TYPE newcTriMat
#define VAR_TYPE complex
#define PREC sp
#define VAR_INIT cmplx(0._sp,0._sp)
#include "class_TriMat.T90"
!========================
end module class_cTriMat

#endif

module class_zTriMat
!========================
#define TYPE_NAME zTriMat
#define STR_TYPE_NAME "zTriMat"
#define TYPE_NAME_ zTriMat_
#define NEW_TYPE newzTriMat
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT dcmplx(0._dp,0._dp)
#include "class_TriMat.T90"
!========================
end module class_zTriMat


