#!/bin/bash

if [ $# != 5 ] 
then
    echo "Usage: $0 siesta file Grid_Min Grid_Max Grid_Increment"
    exit
fi

siesta=$1
input=$2
GridMin=$3
GridMax=$4
GridIncrement=$5

cp ../$input.fdf $input.in

rm -f E-vs-mesh.dat

Grid=$GridMin;

#echo "Grid:" $Grid
#echo "Max:" $GridMax
#echo "Min:" $GridMin
#echo "Increment:" $GridIncrement;

while [ $Grid -le $GridMax ];
do
		echo "Grid:" $Grid " Ryd";
        sed s/rmesh/$Grid/ $input.in > $input.fdf;
        $siesta < $input.fdf >& $input-$Grid.out;
        E=$(echo $(awk '/siesta:         Total =/{print $4}' $input-$Grid.out));
        echo $Grid $E >> E-vs-mesh.dat;
		Grid=$(echo "scale=5; $Grid+$GridIncrement"|bc);
done
