#!/bin/sh
#
# The script is passed the (probably relative) path to the transiesta
# executable
#
### set -x
TS_RAW="$1"
#
# Extract last component, in case of mpirun-style string
# 
TS_REL_PATH=$(echo ${TS_RAW} | awk '{print $NF}')
TS_EXEC_PREFIX=$(echo ${TS_RAW} | awk '{$NF=""; print}')
#
TS_NAME=$(basename ${TS_REL_PATH})
EXEC_DIR=$(dirname ${TS_REL_PATH})
#
# Find absolute path -------
#
pushd ${EXEC_DIR}
ABS_EXEC_DIR=$(pwd)
popd
#---------------------------
TS_ABS=${ABS_EXEC_DIR}/${TS_NAME}
TS="$TS_EXEC_PREFIX $TS_ABS"
OBJDIR=$(basename ${ABS_EXEC_DIR})
ROOT_DIR=$(dirname ${ABS_EXEC_DIR})
echo "Running script with TranSIESTA=$TS, compiled in OBJDIR=${OBJDIR}"

#
# Start with the electrode calculation
#
echo "Electrode Calculation"
mkdir Elec
cd Elec
ln ../../H.psf .
ln ../../elec.fast.fdf .
${TS} < elec.fast.fdf > elec.fast.out 
RETVAL=$?
if [ $RETVAL -ne 0 ]
then
   echo "The electrode calculation did not go well ..."
   exit
fi
cp elec.fast.out ../../ts_fast_elec.out
#
# Go back to base directory
#
cd ..

#
# Scattering region calculation
#
echo "Scattering Region Calculation"
mkdir Scat
cd Scat
ln ../../scat.fast.fdf .
ln ../../H.psf .
# Copy the electrode's .TSHS
ln ../Elec/elec.fast.TSHS .
$TS < scat.fast.fdf > scat.fast.out
RETVAL=$?
if [ $RETVAL -ne 0 ]
then
   echo "The scattering region calculation did not go well ..."
   exit
fi
cp scat.fast.out ../../ts_fast_scat.out
#
# Go back to base directory
#
cd ..

#
# TBTrans calculation
#
echo "TBTrans Calculation"
#
# TBT can be specified in the command line, and will override
# the default location in Util
#
if [ -z "$TBT" ] ; then
  TBT=${ROOT_DIR}/Util/TBTrans/tbtrans
  echo "Compiling $TBT..."
  # Clean in case the compiled version there is not compatible
  (cd "${ROOT_DIR}/Util/TBTrans" ; make OBJDIR="$OBJDIR" clean ;
       make OBJDIR="$OBJDIR" )
  TBT="${TS_EXEC_PREFIX} ${ROOT_DIR}/Util/TBTrans/tbtrans"
fi
#
echo "Running script with tbtrans=$TBT"
mkdir TBT
cd TBT
# Copy input files
ln ../Elec/elec.fast.TSHS .
ln ../Scat/scat.fast.TSHS .
ln ../../scat.fast.fdf .
$TBT < scat.fast.fdf > tbt.out
RETVAL=$?
if [ $RETVAL -ne 0 ]
then
   echo "The scattering region calculation did not go well ..."
   exit
fi
cp tbt.out ../../ts_fast_tbt.out
#
# Go back to base directory
#
cd ..

# If it gets here it's because it finished without error
touch ../completed
