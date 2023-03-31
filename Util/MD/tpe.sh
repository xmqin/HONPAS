#!/bin/sh
#
# Extracts and plots T, p, and E data from an MD simulation.
# It assumes that gnuplot is installed.
# 
# Usage:
#        sh tpe.sh Outputfile.MDE
#
grep -v Step *.MDE | awk '{print NR, $2, $3, $4, $6}' > tpe.dat
rm -f tpe.gnu
#echo "set data style lines" >> tpe.gnu
echo "plot \"tpe.dat\" using 1:2" >> tpe.gnu
#echo "replot \"tpe.dat\" using 1:3" >> tpe.gnu
#echo "replot \"tpe.dat\" using 1:4" >> tpe.gnu
#echo "replot \"tpe.dat\" using 1:5" >> tpe.gnu
#
#
gnuplot -persist < tpe.gnu  &
