#!/bin/sh
#
#  run_script
#
#  Basis optimization with bond-length as 
#  criterion (objective-function: abs(b-b_exp))
#
SIESTA=$HOME/bin/siesta-simplex
#
#
name=$1
#
if [ -d $name ] ; then
  echo "Directory $name exists..."
  exit
fi

mkdir $name
sed -f $name.sed TEMPLATE > $name/$name.fdf
cp *.psf $name
cd $name
$SIESTA < $name.fdf > $name.out
rm -f *.psf *.xml *.psdump      # Optional. To save space
#
# Process bond length info
#
head -2 *.BONDS_FINAL | tail -1  \
     | awk 'BEGIN {bond=0.9584} {print sqrt(($3-bond)**2)}' \
     > OPTIM_OUTPUT
cd ..
