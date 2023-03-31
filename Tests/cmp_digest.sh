#!/bin/sh
#
# cmp_digest -- Generates and compares the digests of two Siesta output files
#               It uses a simple-minded scheme
#
# If only one file is given, the reference is sought in a standard directory which
# must exist (it is linked to by newer versions of Siesta). It can be re-defined
# by setting the environment variable REFERENCE_DIR.
#
# If any differences are found, a file 'OUT.diffs' is left in the
# working directory.
#
# A. Garcia, March 14, 2011
#
# ------------------------------------------ 
# Get absolute path of this script
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
#
if [ -z "${REFERENCE_DIR}" ] ; then
    REFERENCE_DIR=$srcdir/Reference
fi
if [ -z "${SCRIPTS_DIR}" ] ; then
    SCRIPTS_DIR=$srcdir
fi
#
DIGEST_CREATOR=${SCRIPTS_DIR}/out_digest.awk
#
if [ $# = 2 ]
then
    f1=$1
    f2=$2
elif [ $# = 1 ]
then
    f1=$1
    f2=${REFERENCE_DIR}/$1
    echo "Using reference directory: ${REFERENCE_DIR}"
else 
      echo "Usage: [ REFERENCE_DIR=Dir ] $0 file1.out [file2.out]"
      exit
fi
#
if [ ! -r $f1 ] ; then echo "No such file: $f1"; exit ; fi
if [ ! -r $f2 ] ; then echo "No reference file: $f2"; exit ; fi
#
rm -f .tmp_dig1 .tmp_dig2
#
# Extract info and compress whitespace
awk -f $DIGEST_CREATOR $f1 | tr -s ' ' > .tmp_dig1
awk -f $DIGEST_CREATOR $f2 | tr -s ' ' > .tmp_dig2
#
# Get the diffs side-by-side, ignoring whitespace
#
diff -y -w --suppress-common-lines .tmp_dig1 .tmp_dig2 > OUT.diffs
# Erase the diff file if it is empty
if [ ! -s OUT.diffs ]
then
    rm -f OUT.diffs .tmp_dig1 .tmp_dig2
else
    echo "$f1: ** DIFFERENCES FOUND **"
fi
