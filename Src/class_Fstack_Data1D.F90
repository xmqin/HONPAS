module class_Fstack_sData1D

use class_sData1D

character(len=*), parameter :: mod_name="class_Fstack_sData1D.F90"

#define _T_ sData1D
#define FSTACK_NAME Fstack_sData1D

#include "Fstack.T90"

end module class_Fstack_sData1D

module class_Fstack_dData1D

use class_dData1D

character(len=*), parameter :: mod_name="class_Fstack_dData1D.F90"

#define _T_ dData1D
#define FSTACK_NAME Fstack_dData1D

#include "Fstack.T90"

end module class_Fstack_dData1D

module class_Fstack_gData1D
#ifdef GRID_SP
use class_Fstack_sData1D, &
  Fstack_gData1D => Fstack_sData1D
#else
use class_Fstack_dData1D, &
  Fstack_gData1D => Fstack_dData1D
#endif
end module class_Fstack_gData1D
