#!/bin/sh
#
# Extracts and plots cell shape data from a variable-cell dynamics simulation.
# It assumes that gnuplot is installed.
# 
# Usage:
#        sh shape.sh Outputfile.out
#
file=$1
#
grep "vector module" $file | awk '{print NR, $7, $8, $9 }' > vectors.dat
grep "angles" $file | awk '{print NR, $6, $7, $8 }' > angles.dat
#
rm -f vectors.gnu
echo "set data style lines" >> vectors.gnu
echo "plot \"vectors.dat\" using 1:2" >> vectors.gnu
echo "replot \"vectors.dat\" using 1:3" >> vectors.gnu
echo "replot \"vectors.dat\" using 1:4" >> vectors.gnu
#
rm -f angles.gnu
echo "set data style lines" >> angles.gnu
echo "plot \"angles.dat\" using 1:2" >> angles.gnu
echo "replot \"angles.dat\" using 1:3" >> angles.gnu
echo "replot \"angles.dat\" using 1:4" >> angles.gnu

#
gnuplot -persist  vectors.gnu &
gnuplot -persist  angles.gnu  &
