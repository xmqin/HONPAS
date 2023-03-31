#!/bin/sh

##set -x
#
# Run "partial atom" example
#
# Argument: SIESTA execution string
#
# Optional environmental SYMBOL:  OBJDIR
# This is essential if the argument 'SIESTA' does not refer to a typical installation
#
# The script is passed the (probably relative) path to the siesta
# executable, and maybe with a "mpirun" prefix
#
SIESTA="$1"
# 
# Now we try to guess how to build the 'fractional' executable in Util

# Extract last component of the executable, in case of mpirun-style string
REL_PATH=$(echo ${SIESTA} | awk '{print $NF}')
EXEC_PREFIX=$(echo ${SIESTA} | awk '{$NF=""; print}')
REL_PATH=$(which ${REL_PATH})
NAME=$(basename ${REL_PATH})
EXEC_DIR=$(dirname ${REL_PATH})
#
# Find absolute path -------
pushd ${EXEC_DIR} ; ABS_EXEC_DIR=$(pwd) ; popd
#---------------------------
ABS=${ABS_EXEC_DIR}/${NAME}
COMPILATION_DIR=$(basename ${ABS_EXEC_DIR})
echo "Running script with SIESTA=$EXEC_PREFIX $ABS"
#
# Make sure we can use the program location info...
# Use the sentinel .siesta file in modern versions
#
if [  -f $ABS_EXEC_DIR/.siesta ] ; then
    OBJDIR=${OBJDIR:-$COMPILATION_DIR}
else
    if [ -z $OBJDIR ] ; then
      echo "${ABS_EXEC_DIR} does not look like a compilation dir"
      echo "Need to specify OBJDIR"
      exit
    fi
fi
echo "Using OBJDIR=$OBJDIR to build 'fractional'"
#

fractional=../../../../Util/VCA/fractional
#
cp ../../Pseudos/O.psf .

## if [ ! -x $fractional ] ; then
# Compile always
  echo "Compiling $fractional..."
  (cd ../../../../Util/VCA ; make OBJDIR="$OBJDIR" clean fractional) > /dev/null
##   echo -n "Please compile 'fractional' in "
##   echo "Util/VCA before running this test"
##   exit 1
##fi

echo "==> Running $fractional"
$fractional O  0.5

#
cat > Job.fdf << EOF
SystemName          Oxypartial
SystemLabel         oxypartial
NumberOfAtoms       2
NumberOfSpecies     1

MeshCutoff 200 Ry

%block ChemicalSpeciesLabel
 1  201  O-Fraction-0.50000      # Species index, atomic number, species label
%endblock ChemicalSpeciesLabel
%block SyntheticAtoms
 1  
 2 2 3 4
 1.0 2.0 0.0 0.0
%endblock SyntheticAtoms

%block PAO.basis
O-Fraction-0.50000  2                 # Species label, number of l-shells
 n=2   0   2                         # n, l, Nzeta 
 0.0 0.0   
   1.000      1.000   
 n=2   1   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   4.139      2.740   
   1.000      1.000   
%endblock PAO.Basis

DM.NumberPulay 4

Spin.Polarized T

AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 0.000  0.000  0.000  1
 0.000  0.000  1.200  1
%endblock AtomicCoordinatesAndAtomicSpecies

MD.TypeOfRun Broyden
MD.NumCGSteps 40
EOF
#
echo "==> Running $SIESTA"
$SIESTA < Job.fdf > Job.out
cp Job.out ../partial.out

