#!/bin/sh
#
#  run_script
#
#  Pseudopotential optimization prototype based
#  on a penalty function to enforce smoothness.
#
ATM=atm-3.2.5          # You need atm >= 3.2.5
#
# Sigmoid parameters
a=0.1
b=10.0
#
# Maximum qmax allowed
q0=14.0
#
name=$1
#
if [ -d $name ] ; then
  echo "Directory $name exists..."
  exit
fi

mkdir -p $name/gen
sed -f $name.sed PSGEN.TEMPLATE > $name/gen/INP
#
# In case there is trouble, use a large value
cd $name ; echo "1000.0" > OPTIM_OUTPUT
#
cd gen
$ATM
#
# Get maximum qmax
#
qmax=`tail -1 FOURIER_QMAX | awk '{print $2}'`
#
cp VPSOUT ../VPSIN
#
# Get ready for test
#
cd ..
ln -sf ../Test.inp ./INP
$ATM
#
# Process econf diffs info
#
max_ediff=`tail -1  ECONF_DIFFS  | awk '{print $1}' `

#
# Constraint on qmax as sigmoid-like penalty function
#
objective_value=`awk -v a=$a -v b=$b -v ediff=$max_ediff -v qmax=$qmax \
                 -v q0=$q0 \
              'BEGIN {print ediff+b/(1.+exp(-(qmax-q0)/a))}' `
#
echo $objective_value >| OPTIM_OUTPUT
#
# Record what went on in a summary file (one for each job)
#
echo -n "$name " > log
echo -n ' $rcs $rcp $rcd $rcf ' | sed -f ../$name.sed >> log
echo " Ediff:" $max_ediff "Qmax:" $qmax "F:" $objective_value >> log
#
cd ..

