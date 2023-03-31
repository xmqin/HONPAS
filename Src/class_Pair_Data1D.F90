module class_Pair_sData1D
use class_sData1D
character(len=*), parameter :: mod_name="class_Pair_sData1D.F90"

#define PAIR_NAME Pair_sData1D
#define _T1_ sData1D
#define _T2_ sData1D
#include "Pair.T90"

end module class_Pair_sData1D

module class_Pair_dData1D
use class_dData1D
character(len=*), parameter :: mod_name="class_Pair_dData1D.F90"

#define PAIR_NAME Pair_dData1D
#define _T1_ dData1D
#define _T2_ dData1D
#include "Pair.T90"
end module class_Pair_dData1D

module class_Pair_gData1D
#ifdef GRID_SP
use class_Pair_sData1D, Pair_gData1D => Pair_sData1D
#else
use class_Pair_dData1D, Pair_gData1D => Pair_dData1D
#endif
end module class_Pair_gData1D


