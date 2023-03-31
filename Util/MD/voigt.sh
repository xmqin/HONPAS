#!/bin/sh
#
# Extracts and plots stress data from an MD simulation.
# It assumes that gnuplot is installed.
# 
# Usage:
#        sh voigt.sh Outputfile.out
#
grep oigt OUT | awk '{$1=""; $2=""; print NR, $0}' > stress.dat
rm -f voigt.gnu
#echo "set data style lines" >> voigt.gnu
echo "plot \"stress.dat\" using 1:2" >> voigt.gnu
echo "replot \"stress.dat\" using 1:3" >> voigt.gnu
echo "replot \"stress.dat\" using 1:4" >> voigt.gnu
echo "replot \"stress.dat\" using 1:5" >> voigt.gnu
echo "replot \"stress.dat\" using 1:6" >> voigt.gnu
echo "replot \"stress.dat\" using 1:7" >> voigt.gnu
#
#
gnuplot -persist < voigt.gnu  &
