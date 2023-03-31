#!/bin/sh
#
# Copies results of TS tests to Reference directory
#
if [ $# != 1 ] 
then
   echo "Usage: $0 Reference_Directory (no trailing /)"
   exit
fi
#
refdir=$1
#
rm -f _sources
cat > _sources <<EOF
ts_au/au_111_capacitor.TBT.TRANS_Left-Right                     
ts_au/au_111_capacitor.out
ts_au/au_111_capacitor_single.TBT.AVADOS_Left
ts_au/au_111_capacitor_single.out
ts_au/bulk_au_111.TBT.TRANS_Left-Right
ts_au/bulk_au_111.out
ts_au/elec_au_111_abc.out
ts_au/tbt_au_111_capacitor.out
ts_au/tbt_bulk_au_111.out
ts_au_100/au_100.TBT.TRANS_Left-Right
ts_au_100/au_100.out
ts_au_100/elec_au_100.out
ts_au_100/tbt_au_100.out
ts_au_100_0.25V/au_100.TBT.TRANS_Left-Right
ts_au_100_0.25V/au_100.out
ts_au_100_0.25V/elec_au_100.out
ts_au_100_0.25V/tbt_au_100.out
ts_au_100_repetition/au_100.TBT.TRANS_Left-Right
ts_au_100_repetition/au_100.out
ts_au_100_repetition/elec_au_100.out
ts_au_100_repetition/tbt_au_100.out
ts_au_100_repetition_0.25V/au_100.TBT.TRANS_Left-Right
ts_au_100_repetition_0.25V/au_100.out
ts_au_100_repetition_0.25V/elec_au_100.out
ts_au_100_repetition_0.25V/tbt_au_100.out
ts_au_repetition/au_111_capacitor.TBT.TRANS_Left-Right
ts_au_repetition/au_111_capacitor.out
ts_au_repetition/elec_au_111_abc.out
ts_au_repetition/tbt_au_111_capacitor.out
ts_graphene/elec.out
ts_graphene/graphene-UA-tbt.out
ts_graphene/graphene-UA.TBT.TRANS_Left-Right
ts_graphene/graphene-UA.out
ts_graphene/graphene-tbt.out
ts_graphene/graphene.TBT.TRANS_Left-Right
ts_graphene/graphene.out
ts_graphene/tbt_graphene-UA.out
ts_graphene/tbt_graphene.out
ts_term3/elec-x.out
ts_term3/elec-z.out
ts_term3/tbt_term3.out
ts_term3/term3-tbt.out
ts_term3/term3.TBT.AVADOS_el-3
ts_term3/term3.TBT.AVTRANS_el-1-el-3
ts_term3/term3.TBT.AVTRANS_el-1-el-4
ts_term3/term3.TBT.AVTRANS_el-3-el-4
ts_term3/term3.out
ts_term4/elec-x.out
ts_term4/elec-z.out
ts_term4/tbt_term4.out
ts_term4/term4-tbt.out
ts_term4/term4.TBT.AVADOS_el-1
ts_term4/term4.TBT.AVTRANS_el-1-el-2
ts_term4/term4.TBT.AVTRANS_el-1-el-3
ts_term4/term4.TBT.AVTRANS_el-1-el-4
ts_term4/term4.TBT.AVTRANS_el-2-el-3
ts_term4/term4.TBT.AVTRANS_el-2-el-4
ts_term4/term4.TBT.AVTRANS_el-3-el-4
ts_term4/term4.out
EOF
#
for i in `cat _sources`; do
 target=`echo $i | cut -d "/" -f 2`
 echo "-- Copying $i to $target"
 cp -p $i ${refdir}/${target}
done
