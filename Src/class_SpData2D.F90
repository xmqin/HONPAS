! This file holds all the different modules for sparse vectors
! This makes editing easier as we do not need to consider *all* files

! A logical sparse array
module class_lSpData2D
  use class_lData2D
!========================
#define TYPE_NAME  lSpData2D
#define STR_TYPE_NAME "lSpData2D"
#define TYPE_NAME_ lSpData2D_
#define NEW_TYPE newlSpData2D
#define VAR_TYPE lData2D
#define VAR_NEW_TYPE newlData2D
#define VAR_TYPE_TYPE logical
! Do not define precision
#include "class_SpData2D.T90"
!========================
end module class_lSpData2D

module class_iSpData2D
  use class_iData2D
!========================
#define TYPE_NAME  iSpData2D
#define STR_TYPE_NAME "iSpData2D"
#define TYPE_NAME_ iSpData2D_
#define NEW_TYPE newiSpData2D
#define VAR_TYPE iData2D
#define VAR_NEW_TYPE newiData2D
#define VAR_TYPE_TYPE integer
! Do not define precision
#include "class_SpData2D.T90"
!========================
end module class_iSpData2D

module class_sSpData2D
  use class_sData2D
!========================
#define TYPE_NAME  sSpData2D
#define STR_TYPE_NAME  "sSpData2D"
#define TYPE_NAME_ sSpData2D_
#define NEW_TYPE newsSpData2D
#define VAR_TYPE sData2D
#define VAR_NEW_TYPE newsData2D
#define VAR_TYPE_TYPE real
#define PREC sp
#include "class_SpData2D.T90"
!========================
end module class_sSpData2D

module class_dSpData2D
  use class_dData2D
!========================
#define TYPE_NAME  dSpData2D
#define STR_TYPE_NAME  "dSpData2D"
#define TYPE_NAME_ dSpData2D_
#define NEW_TYPE newdSpData2D
#define VAR_TYPE dData2D
#define VAR_NEW_TYPE newdData2D
#define VAR_TYPE_TYPE real
#define PREC dp
#include "class_SpData2D.T90"
!========================
end module class_dSpData2D

module class_cSpData2D
  use class_cData2D
!========================
#define TYPE_NAME  cSpData2D
#define STR_TYPE_NAME  "cSpData2D"
#define TYPE_NAME_ cSpData2D_
#define NEW_TYPE newcSpData2D
#define VAR_TYPE cData2D
#define VAR_NEW_TYPE newcData2D
#define VAR_TYPE_TYPE complex
#define PREC sp
#include "class_SpData2D.T90"
!========================
end module class_cSpData2D

module class_zSpData2D
  use class_zData2D
!========================
#define TYPE_NAME  zSpData2D
#define STR_TYPE_NAME  "zSpData2D"
#define TYPE_NAME_ zSpData2D_
#define NEW_TYPE newzSpData2D
#define VAR_TYPE zData2D
#define VAR_NEW_TYPE newzData2D
#define VAR_TYPE_TYPE complex
#define PREC dp
#include "class_SpData2D.T90"
!========================
end module class_zSpData2D

