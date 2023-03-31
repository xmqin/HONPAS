#!/bin/sh
#
# Drives ioncat to plot orbitals, KB projs, and Vna
# Examples:
#
#  sh ionplot.sh -o 3 H.ion
# 
# for i in $(ioncat -i H.ion) ; do
#  ionplot.sh -Z -o $i H.ion
# done
#
# ----------
#
# Make sure you can find the ioncat executable
#
echo $0
IONCAT=$(basename $0)/ioncat
#
## set -x
string="$@"
$IONCAT $@  > .tmp_ioncat
cat > .plot.g << EOF
set title "$string"
plot ".tmp_ioncat" using 1:2 w l title "f"
replot ".tmp_ioncat" using 1:3 w l title "grad f"
EOF
#
gnuplot -persist .plot.g
#


