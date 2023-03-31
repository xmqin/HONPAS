#!/bin/sh

#
# Select the appropriate run below (at the end).
#
# Take care to use the appropriate (e.g., parallel or serial) copy
# of siesta. To run parallel jobs, you have to figure out how to
# set the MPI environment. For a simple single-node interactive calculation
# (where possible, such as in a Rocks cluster), you can use the parallel.sh
# script.
#
# Make sure that you are using the right version of Siesta. 
# The SIESTA setting below is quite naive and might not work
# in all cases. You can call this script as:
#
#            SIESTA=/path/to/siesta test.sh
#
#

ROOT="../../../.."
PSEUDOS=${ROOT}/Tests/Pseudos
#
if [ -z "$SIESTA" ] ; then
      SIESTA=${ROOT}/Obj/siesta
fi
echo "Using Siesta executable: $SIESTA"
#

if [ -d work ] ; then
   echo "Work directory 'work' exists. Please delete it"
   exit
else
   mkdir work
fi

cp -p h2o.fast.fdf h2o.conv.fdf driver.dat work
#
cd work
cp ${PSEUDOS}/H.psf  .
cp ${PSEUDOS}/O.psf  .
#
ln -sf ${SIESTA} ./siesta

../Src/driver < driver.dat | tee driver.out
../Src/simple  | tee simple.out

#
# Make sure you edit Src/para.f90 to suit your system before
# you un-comment this
#
#../Src/para  | tee para.out

