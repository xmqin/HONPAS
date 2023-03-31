#!/bin/sh

rm -f file.2d
#
grid2d < raw.in > file.2d
gnuplot -persist cont.gplot
