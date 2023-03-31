#!/bin/bash
#
# re-creates all dependencies
#
# Usage: [ OBJDIR=Objdir] bash dep_all.sh
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
    f=${i##*/}
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

    # Check whether the "DO NOT DELETE THIS LINE" exists
    # in the makefile. If so, simply delete everything after.
    x=$(grep 'DO NOT DELETE THIS LINE' $f)
    x=${x// /}
    if [ ${#x} -gt 10 ]; then
	sed '/DO NOT DELETE THIS LINE/q' $f > $f.old
	mv $f.old $f
    fi
    make OBJDIR=${OBJDIR} dep 2>/dev/null
    if [ $? -ne 0 ]; then
	echo "*** DEPENDENCY FAILED in $relpath ***"
	failed="$failed $relpath"
    fi
    cd $topdir
done

if [ ! -z "$failed" ]; then
    echo "" #empty line
    echo "" #empty line
    echo " *** All failed directories (DEPENDENCIES):"
    for p in $failed ; do
	echo "   $p"
    done
fi
