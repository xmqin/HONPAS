#!/bin/sh
#
# Output comparison helper
# Please note that the Reference directory was populated
# with runs with the Intel compiler, which is notoriously flaky.
# 
REFERENCE=../../ref-version/Tests
#
tests="h2o h2o_reparam h2o_basis h2o_dos h2o_orderN \
       floating bessel mgco3 si2x1h force_2 born\
       var_cell constant_volume batio3 fe fe_broyden sih si64 \
       h2oZ sih_op_broyden h2o_op_broyden zmatrix md_pr md_npr \
       md_anneal md_verlet md_nose si_bandpoints sih_fire \
       graphite_c6 oxyn partial h2o_findp_bug h2o_radialgrid"

for i in $tests ;do
echo "---------$i start"
diff $i/$i.out ${REFERENCE}/$i/$i.out
echo "****---------$i done"
done | grep -v " timer: " | grep -v " elaps: " | grep -v " of run"




