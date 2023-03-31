! We define all variables in this one
! As the entries are so short this makes more sense
module class_lData2D
!========================
#define TYPE_NAME lData2D
#define STR_TYPE_NAME "lData2D"
#define TYPE_NAME_ lData2D_
#define NEW_TYPE newlData2D
#define VAR_TYPE logical
#define VAR_INIT .false.
#include "class_Data2D.T90"
!========================
end module class_lData2D

module class_iData2D
!========================
#define TYPE_NAME iData2D
#define STR_TYPE_NAME "iData2D"
#define TYPE_NAME_ iData2D_
#define NEW_TYPE newiData2D
#define VAR_TYPE integer
#define VAR_INIT 0
#include "class_Data2D.T90"
!========================
end module class_iData2D


module class_sData2D
!========================
#define TYPE_NAME sData2D
#define STR_TYPE_NAME "sData2D"
#define TYPE_NAME_ sData2D_
#define NEW_TYPE newsData2D
#define VAR_TYPE real
#define PREC sp
#define VAR_INIT 0._sp
#include "class_Data2D.T90"
!========================
end module class_sData2D

module class_dData2D
!========================
#define TYPE_NAME dData2D
#define STR_TYPE_NAME "dData2D"
#define TYPE_NAME_ dData2D_
#define NEW_TYPE newdData2D
#define VAR_TYPE real
#define PREC dp
#define VAR_INIT 0._dp
#include "class_Data2D.T90"
!========================
end module class_dData2D

module class_cData2D
!========================
#define TYPE_NAME cData2D
#define STR_TYPE_NAME "cData2D"
#define TYPE_NAME_ cData2D_
#define NEW_TYPE newcData2D
#define VAR_TYPE complex
#define PREC sp
#define VAR_INIT cmplx(0._sp,0._sp)
#include "class_Data2D.T90"
!========================
end module class_cData2D

module class_zData2D
!========================
#define TYPE_NAME zData2D
#define STR_TYPE_NAME "zData2D"
#define TYPE_NAME_ zData2D_
#define NEW_TYPE newzData2D
#define VAR_TYPE complex
#define PREC dp
#define VAR_INIT dcmplx(0._dp,0._dp)
#include "class_Data2D.T90"
!========================
end module class_zData2D


