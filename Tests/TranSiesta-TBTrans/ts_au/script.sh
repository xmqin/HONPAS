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
ELEC=elec_au_111_abc
echo "Electrode Calculation"
mkdir Elec
cd Elec
ln ../../Au.psf .
ln ../../${ELEC}.fdf .
${TS} < ${ELEC}.fdf > ${ELEC}.out
RETVAL=$?
if [ $RETVAL -ne 0 ]
then
   echo "The electrode calculation did not go well ..."
   exit
fi
cp ${ELEC}.out ../..
#
# Go back to base directory
#
cd ..

#
# Scattering region calculation
#
for SCAT in bulk_au_111 au_111_capacitor 
do
  echo "==> Scattering Region Calculation for $SCAT"
  mkdir Scat_$SCAT
  cd Scat_$SCAT
  ln ../../Au.psf .
  ln ../../${SCAT}.fdf .
  # Copy the electrode's .TSHS
  ln ../Elec/${ELEC}.TSHS .
  $TS < ${SCAT}.fdf > ${SCAT}.out
  RETVAL=$?
  if [ $RETVAL -ne 0 ]
      then
      echo "** The scattering region calculation for $SCAT did not go well ..."
      exit
  fi
  cp ${SCAT}.out ../..
#
# Go back to base directory
#
  cd ..

#
# TBTrans calculation
#
  echo "==> TBTrans Calculation for $SCAT"
#
# TBT can be specified in the command line, and will override
# the default location in Util
#
 if [ -z "$TBT" ] ; then
  TBT=${ROOT_DIR}/Util/TBTrans/tbtrans
  echo "==> (Compiling $TBT...)"
  # Clean in case the compiled version there is not compatible
  (cd "${ROOT_DIR}/Util/TBTrans" ; make OBJDIR="$OBJDIR" clean ;
       make OBJDIR="$OBJDIR" )
  TBT="${TS_EXEC_PREFIX} ${ROOT_DIR}/Util/TBTrans/tbtrans"
 fi
#
 echo "==> Running $SCAT with tbtrans=$TBT"
 mkdir TBT_$SCAT
 cd TBT_$SCAT
 # Copy input files
 ln ../Elec/${ELEC}.TSHS .
 ln ../Scat_$SCAT/${SCAT}.TSHS .
 ln ../../${SCAT}.fdf .
 $TBT < ${SCAT}.fdf  > tbt_${SCAT}.out
 RETVAL=$?
 if [ $RETVAL -ne 0 ]
 then
   echo "The scattering region calculation did not go well ..."
   exit
 fi
 cp tbt_${SCAT}.out ../..
 #
 # Go back to base directory
 #
 cd ..
done
# If it gets here it's because it finished without error
touch ../completed
