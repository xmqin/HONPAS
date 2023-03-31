#!/bin/sh
#
# Script to run a basis generation calculation
# Do not move directory from Siesta tree, or modify
# the lines below to help the script find the gen-basis program
# (note also re-definition of prog further down)
#
BASIS_TUTORIAL_DIR=`dirname $0`
EXEC_DIR=${BASIS_TUTORIAL_DIR}/../../Util/Gen-basis
GEN_BASIS_PROGRAM=${EXEC_DIR}/gen-basis
prog=${GEN_BASIS_PROGRAM} 
#
#
if [ "$#" != 2 ] 
then
      echo "Usage: sh gen-basis.sh <ChemLabel.fdf> <Pseudo.psf> "
      echo "       where ChemLabel _must_ be the same as in the"
      echo "       information blocks in the file"
      echo "       The pseudopotential file can have any name"
      exit
fi
#
fdf_file=$1
pseudo_file=$2
label=`basename $fdf_file .fdf`
#
dirname=$label
#
if [ -d $dirname ] 
then
        echo "Directory $dirname exists. Please delete it first"
        exit
fi
#
mkdir $dirname ; cd $dirname

prog=../$prog
BASIS_TUTORIAL_DIR=../${BASIS_TUTORIAL_DIR}

cp ../$fdf_file .
#
#  Make the pseudo file have the same root as the run Label
#
ln -sf ../$pseudo_file $label.psf
#
$prog < $fdf_file > OUT
ln -fs OUT $label.OUT
#
#
echo " "
echo "==> Calculation completed for $label "
echo " "
echo "The files necessary to plot the Basis information for $label"
echo "are now in directory $dirname. You can use the Gnuplot"
echo "scripts to visualize the PAOs, KB projectors, and other "
echo "interesting functions."
echo " "
echo " *** (Be sure to use """gnuplot -persist""" if using X"
echo " *** For postscript output, use the .gps files"
#
#  Copy plotting scripts
#
cp -f ${BASIS_TUTORIAL_DIR}/*.gp* .
#
# Rename files
#
for i in ORB.S?.?.$label KB.L?.?.$label ; do
      name=`basename $i .$label`
      mv $i $name
done
#
mv CHLOCAL.$label CHLOCAL
mv RED_VLOCAL.$label RED_VLOCAL
mv VNA.$label VNA
#
if [ -f CHCORE.$label ]
then
      mv CHCORE.$label CHCORE
fi

#


