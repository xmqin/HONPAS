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

ROOT=$(cd -P -- "../../.." && pwd -P)
PSEUDOS=${ROOT}/Tests/Pseudos
#
if [ -z "$SIESTA" ] ; then
    if [ ! -z "$OBJDIR" ] ; then
	SIESTA_DIR=${ROOT}/$OBJDIR
    fi
    if [ -z "$SIESTA_DIR" ] ; then
	SIESTA_DIR=${ROOT}/Obj
    fi
    echo "Using Siesta dir: ${SIESTA_DIR}"
    SIESTA=$(cd -P -- ${SIESTA_DIR} && pwd -P)/siesta
fi

if [ ! -e "$SIESTA" ] ; then
   echo "Cannot find $SIESTA"
   exit
fi
echo "Using Siesta executable: $SIESTA"
#

#rm -r work
if [ -d work ] ; then
   echo "Work directory 'work' exists. Please delete it"
   exit
fi

#
mkdir work
cd work
cp -p ../*.fdf .
cp ${PSEUDOS}/H.psf  .
cp ${PSEUDOS}/O.psf  .
ln -sf ${SIESTA} ./siesta
#

cat socket.fdf >> h2o.fdf

echo ""; echo "simple_sockets_serial"
../Src/simple_sockets_serial    | tee simple_sockets_serial.out
mv h2o.out siesta_sockets_serial.out


