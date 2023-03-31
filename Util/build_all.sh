#!/bin/sh
#
# build all Utils
#
# Usage: [ OBJDIR=Objdir] sh build_all
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
topdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
#---------------------------------------------
#
if [ -z "$OBJDIR" ] ; then
    OBJDIR=Obj
fi
echo "Using: ${OBJDIR}/arch.make"
echo ""
echo ""
sleep 1

# Collect all failed directories to print them out in the end
failed=""

for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile ); do
    relpath=${i%/*}
    sub=${relpath##*/}
    # Skip these directories
    case $sub in
	fdict|ncdf)
	    continue
	    ;;
    esac
    
    echo ""
    echo "====> Processing $relpath ..."
    echo ""
    cd $relpath
    
    make OBJDIR=${OBJDIR} clean 
    make OBJDIR=${OBJDIR}
    if [ $? -ne 0 ]; then
	echo "*** COMPILATION FAILED in $relpath ***"
	failed="$failed $relpath"
    fi
    cd $topdir
done

if [ ! -z "$failed" ]; then
    echo "" #empty line
    echo "" #empty line
    echo " *** All failed directories:"
    echo " *** (Some programs have to be compiled after compiling Siesta)"
    for p in $failed ; do
	echo "   $p"
    done
fi
