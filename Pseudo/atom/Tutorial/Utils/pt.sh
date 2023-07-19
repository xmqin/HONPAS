#!/bin/sh -f
#
# pt.sh -- Script to run pseudopotential test calculations
#
# Usage: pt.sh <ptname.inp> <psname.vps>
#
# Make sure that atm is in path
#
DEFAULT_DIR=../../../Utils
ATOM_UTILS_DIR=${ATOM_UTILS_DIR:-${DEFAULT_DIR}}
#
default="../../../../atm"
prog=${ATOM_PROGRAM:-$default}
#
if [ "$#" != 2 ] 
then
	echo "Usage: $0 <ptname.inp> <psname.vps>"
	exit
fi
#
file=$1
psfile=$2
ptname=`basename $file .inp`
psname=`basename $psfile .vps`
name="$ptname-$psname"
#
#
if [ -d $name ] 
then
	echo "Directory $name exists. Please delete it first"
	exit
fi
#
mkdir $name ; cd $name
cp ../$file ./INP
cp ../$psfile ./VPSIN
#
$prog
#
echo "==> Output data in directory $name"
#
#  Copy plotting scripts
#
for i in charge vcharge vspin pt  ; do
            cp -f ${ATOM_UTILS_DIR}/$i.gps .
            cp -f ${ATOM_UTILS_DIR}/$i.gplot .
done


