#!/bin/sh

for i in  2.5 3.0 3.5 4.0 4.5 5.0 
do
   sed "s/ALAT/$i/g" template >| he_$i.fdf
  ../../Src/siesta < he_$i.fdf | tee he_$i.out 
  done
