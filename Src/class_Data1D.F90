! We define all variables in this one
! As the entries are so short this makes more sense

module class_lData1D
!========================
#define TYPE_NAME lData1D
#define STR_TYPE_NAME "lData1D"
#define TYPE_NAME_ lData1D_
#define NEW_TYPE newlData1D
#define VAR_TYPE logical
#define VAR_INIT .false.
! DO NOT DEFINE PREC
#include "class_Data1D.T90"
!========================
end module class_lData1D

module class_iData1D
!========================
#define TYPE_NAME iData1D
#define STR_TYPE_NAME "iData1D"
#define TYPE_NAME_ iData1D_
#define NEW_TYPE newiData1D
#define VAR_TYPE integer
#define VAR_INIT 0
! DO NOT DEFINE PREC
#include "class_Data1D.T90"
!========================
end module class_iData1D


module class_sData1D
!========================
#define TYPE_NAME sData1D
#define STR_TYPE_NAME "sData1D"
#define TYPE_NAME_ sData1D_
#define NEW_TYPE newsData1D
#define VAR_TYPE real
#define PREC sp
#define VAR_INIT 0._sp
#include "class_Data1D.T90"
!========================
end module class_sData1D

module class_dData1D
!========================
#define TYPE_NAME dData1D
#define STR_TYPE_NAME "dData1D"
#define TYPE_NAME_ dData1D_
#define NEW_TYPE newdData1D
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_Data1D.T90"
!========================
end module class_dData1D

module class_gData1D
#ifdef GRID_SP
use class_sData1D, gData1D => sData1D, newgData1D => newsData1D
#else
use class_dData1D, gData1D => dData1D, newgData1D => newdData1D
#endif
end module class_gData1D

module class_cData1D
!========================
#define TYPE_NAME cData1D
#define STR_TYPE_NAME "cData1D"
#define TYPE_NAME_ cData1D_
#define NEW_TYPE newcData1D
#define VAR_TYPE complex
#define PREC sp
#define VAR_INIT cmplx(0._sp,0._sp)
#include "class_Data1D.T90"
!========================
end module class_cData1D

module class_zData1D
!========================
#define TYPE_NAME zData1D
#define STR_TYPE_NAME "zData1D"
#define TYPE_NAME_ zData1D_
#define NEW_TYPE newzData1D
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT dcmplx(0._dp,0._dp)
#include "class_Data1D.T90"
!========================
end module class_zData1D


