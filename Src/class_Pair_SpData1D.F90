module class_Pair_sSpData1D

use class_sSpData1D

character(len=*), parameter :: mod_name="class_Pair_sSpData1D.F90"

#define PAIR_NAME Pair_sSpData1D
#define _T1_ sSpData1D
#define _T2_ sSpData1D

#include "Pair.T90"

end module class_Pair_sSpData1D

module class_Pair_dSpData1D

use class_dSpData1D

character(len=*), parameter :: mod_name="class_Pair_dSpData1D.F90"

#define PAIR_NAME Pair_dSpData1D
#define _T1_ dSpData1D
#define _T2_ dSpData1D

#include "Pair.T90"

end module class_Pair_dSpData1D

module class_Pair_gSpData1D
#ifdef GRID_SP
use class_Pair_sSpData1D, Pair_gSpData1D => Pair_sSpData1D
#else
use class_Pair_dSpData1D, Pair_gSpData1D => Pair_dSpData1D
#endif
end module class_Pair_gSpData1D


