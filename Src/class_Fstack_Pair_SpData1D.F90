module class_Fstack_Pair_sSpData1D

use class_Pair_sSpData1D

character(len=*), parameter :: mod_name="class_Fstack_sSpData1D.F90"

#define _T_ Pair_sSpData1D
#define FSTACK_NAME Fstack_Pair_sSpData1D

#include "Fstack.T90"

end module class_Fstack_Pair_sSpData1D

module class_Fstack_Pair_dSpData1D

use class_Pair_dSpData1D

character(len=*), parameter :: mod_name="class_Fstack_dSpData1D.F90"

#define _T_ Pair_dSpData1D
#define FSTACK_NAME Fstack_Pair_dSpData1D

#include "Fstack.T90"

end module class_Fstack_Pair_dSpData1D

module class_Fstack_Pair_gSpData1D
#ifdef GRID_SP
use class_Fstack_Pair_sSpData1D, &
  Fstack_Pair_gSpData1D => Fstack_Pair_sSpData1D
#else
use class_Fstack_Pair_dSpData1D, &
  Fstack_Pair_gSpData1D => Fstack_Pair_dSpData1D
#endif
end module class_Fstack_Pair_gSpData1D
