! This file holds all the different modules for sparse vectors
! This makes editing easier as we do not need to consider *all* files

! A logical sparse array
module class_lSpData1D
  use class_lData1D
!========================
#define TYPE_NAME  lSpData1D
#define STR_TYPE_NAME "lSpData1D"
#define TYPE_NAME_ lSpData1D_
#define NEW_TYPE newlSpData1D
#define VAR_TYPE lData1D
#define VAR_NEW_TYPE newlData1D
#define VAR_TYPE_TYPE logical
! Do not define precision
#include "class_SpData1D.T90"
!========================
end module class_lSpData1D

module class_iSpData1D
  use class_iData1D
!========================
#define TYPE_NAME  iSpData1D
#define STR_TYPE_NAME "iSpData1D"
#define TYPE_NAME_ iSpData1D_
#define NEW_TYPE newiSpData1D
#define VAR_TYPE iData1D
#define VAR_NEW_TYPE newiData1D
#define VAR_TYPE_TYPE integer
#define VAR_INIT 0
! Do not define precision
#include "class_SpData1D.T90"
!========================
end module class_iSpData1D

module class_sSpData1D
  use class_sData1D
!========================
#define TYPE_NAME  sSpData1D
#define STR_TYPE_NAME "sSpData1D"
#define TYPE_NAME_ sSpData1D_
#define NEW_TYPE newsSpData1D
#define VAR_TYPE sData1D
#define VAR_NEW_TYPE newsData1D
#define VAR_TYPE_TYPE real
#define PREC sp
#include "class_SpData1D.T90"
!========================
end module class_sSpData1D

module class_dSpData1D
  use class_dData1D
!========================
#define TYPE_NAME  dSpData1D
#define STR_TYPE_NAME "dSpData1D"
#define TYPE_NAME_ dSpData1D_
#define NEW_TYPE newdSpData1D
#define VAR_TYPE dData1D
#define VAR_NEW_TYPE newdData1D
#define VAR_TYPE_TYPE real
#define PREC dp
#include "class_SpData1D.T90"
!========================
end module class_dSpData1D

module class_gSpData1D
#ifdef GRID_SP
  use class_sSpData1D, newgSpData1D => newsSpData1D, gSpData1D => sSpData1D
#else
  use class_dSpData1D, newgSpData1D => newdSpData1D, gSpData1D => dSpData1D
#endif
end module class_gSpData1D

module class_cSpData1D
  use class_cData1D
!========================
#define TYPE_NAME  cSpData1D
#define STR_TYPE_NAME "cSpData1D"
#define TYPE_NAME_ cSpData1D_
#define NEW_TYPE newcSpData1D
#define VAR_TYPE cData1D
#define VAR_NEW_TYPE newcData1D
#define VAR_TYPE_TYPE complex
#define PREC sp
#include "class_SpData1D.T90"
!========================
end module class_cSpData1D

module class_zSpData1D
  use class_zData1D
!========================
#define TYPE_NAME  zSpData1D
#define STR_TYPE_NAME "zSpData1D"
#define TYPE_NAME_ zSpData1D_
#define NEW_TYPE newzSpData1D
#define VAR_TYPE zData1D
#define VAR_NEW_TYPE newzData1D
#define VAR_TYPE_TYPE complex
#define PREC dp
#include "class_SpData1D.T90"
!========================
end module class_zSpData1D

