#!/bin/bash

if [ $# != 2 ] 
then
    echo "Usage: $0 siesta file"
    exit
fi

siesta=$1
input=$2

cp ../$input.fdf .
#Find the maximum/minum values of the displacement
grid=$(awk '/MeshCutoff/{print $2}' $input.fdf);
gridPeriod=$(echo  "scale=5; 3.14159/sqrt($grid)"|bc);
echo "grid:" $grid " Ry"
echo " corresponding dx:" $gridPeriod "Bohr"

displMax=$(echo "scale=5;2*0.529177*$gridPeriod"|bc);
#echo $displMax " Ang"
displMin=0.0;
dispdx=$(echo "scale=5;$gridPeriod/10.0"|bc);
steps=$(echo "($displMax-$displMin)/$dispdx" |bc);

let steps=steps+1;

#echo $dispdx,$steps;

rm -f *.dat;
rm -f *.ps;

let l=0;
POS=$displMin;

cp $input.fdf $input.in;

cat >> $input.in <<End-of-file
%block AtomicCoordinatesOrigin
POS 0.00 0.00 1
%endblock AtomicCoordinatesOrigin
End-of-file

Emin=0;
Emax=-1000000;
while [ $l -lt $steps ];
  do

      #echo "l=" $l
  

  echo "  POS="$POS
  sed  -e "s/POS/$POS/"  $input.in > $input.fdf;
  $siesta < $input.fdf >  $input-$POS.out;
  F=$(awk '/constr/{print $2}' $input-$POS.out);
  E=$(awk '/Total =/{print $4}' $input-$POS.out);
  Emin=$(echo "define min(x,y){scale=5;if (x <= y) {return (x)}; if (y < x) {return (y)}}; min($Emin,$E)"|bc);
  Emax=$(echo "define max(x,y){scale=5;if (x >= y) {return (x)}; if (y > x) {return (y)}}; max($Emax,$E)"|bc);
  #echo "Emin: "$Emin
  echo "  " $E $F;
  echo $POS $E  $F >>  "E-F-vs-disp.dat";
  #echo "Emin,Emax:" $Emin $Emax
  rm -rf INPUT*;
  
  let l=l+1
  POS=$(echo  "scale=5; $POS+$dispdx"|bc);

done

echo "----------------------------"

EggBoxSize=0;
Emax=$(echo "scale=5; -1.0*$Emax"|bc);
EggBoxSize=$(echo "scale=5; -1.0*($Emin+$Emax)"|bc);
echo "Eggbox size:" $EggBoxSize " eV"

