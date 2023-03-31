#!/bin/bash

# Installation script for zlib, hdf5, netcdf-c and netcdf-fortran
# with complete CDF-4 support (in serial).
# This installation script has been written by:
#  Nick R. Papior, 2016-2017.
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL license.

# VERY BASIC installation script of required libraries
# for installing these packages:
#   flook
# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.
# Optionally you could use VENDOR=gnu|intel
# to use defaults for the specific compiler.

# If you have downloaded other versions edit these version strings
f_v=0.8.1

# Install path, change accordingly
# You can change this variable to control the installation path
# If you want the installation path to be a "packages" folder in
# your home directory, change to this:
# ID=$HOME/packages
if [ -z $PREFIX ]; then
    ID=$(pwd)/build
else
    ID=$PREFIX
fi

while [ $# -gt 0 ]; do
    opt=$1 ; shift
    case $opt in
	--prefix|-p)
	    ID=$1 ; shift
	    ;;
	--flook-version|-flook-v)
	    f_v=$1 ; shift
	    ;;
	--help|-h)
	    echo " $0 --help shows this message"
	    echo ""
	    echo "These options are available:"
	    echo ""
	    echo "  --prefix|-p <>: specify the installation directory of the library"
	    echo "  --flook-version|-flook-v <>: specify the flook version (default: $f_v)"
	    echo ""
	    echo "To customize compilers and flags please export these environment variables:"
	    echo "  CC"
	    echo "  FC"
	    echo "  CFLAGS"
	    echo "  FFLAGS"
	    echo ""
	    exit 0
	    ;;
    esac
done


echo "Installing libraries in folder: $ID"
mkdir -p $ID

# First we check that the user have downloaded the files
function file_exists {
    if [ ! -e $(pwd)/$1 ]; then
	echo "I could not find file $1..."
	echo "Please download the file and place it in this folder:"
	echo " $(pwd)"
	exit 1
    fi
}

# Download a file, if able and the file does not exist
which wget > /dev/null
if [ $? -eq 0 ]; then
    # success we can download using wget
    function _dwn_file {
	wget -O $1 $2
    }
else
    function _dwn_file {
	curl -o $1 $2
    }
fi

# Use download function
#  $1 is name of file
#  $2 is URL
function download_file {
    if [ ! -e $(pwd)/$1 ] ; then
	# Try and download
	_dwn_file $1 $2
    fi
}

# Check for function $?
function retval {
    local ret=$1
    local info="$2"
    shift 2
    if [ $ret -ne 0 ]; then
	echo "Error: $ret"
	echo "$info"
	exit 1
    fi
}

# Download files if they can
download_file flook-${f_v}.tar.gz https://github.com/ElectronicStructureLibrary/flook/archive/v$f_v.tar.gz

file_exists flook-${f_v}.tar.gz
unset file_exists

#################
# Install flook #
#################
[ -d $ID/flook/${f_v}/lib64 ] && flook_lib=lib64 || flook_lib=lib
if [ ! -d $ID/flook/${f_v}/$flook_lib ]; then
    rm -rf flook-${f_v}
    tar xfz flook-${f_v}.tar.gz
    cd flook-${f_v}
    mkdir -p obj ; cd obj
    {
	echo "# Setup script creation of setup.make"
	[ "x$FC" != "x" ] && \
	    echo "FC = $FC"
	[ "x$FCFLAGS" != "x" ] && \
	    echo "FFLAGS = $FCFLAGS"
	[ "x$FFLAGS" != "x" ] && \
	    echo "FFLAGS = $FFLAGS"
	[ "x$CC" != "x" ] && \
	    echo "CC = $CC"
	[ "x$CFLAGS" != "x" ] && \
	    echo "CFLAGS = $CFLAGS"
    } > setup.make
    {
	echo TOP_DIR=..
	echo include ../Makefile
    } > Makefile
    if [ "x$VENDOR" == "x" ]; then
	make liball
    else
	make VENDOR=$VENDOR liball
    fi
    retval $? "flook make liball"
    make install PREFIX=$ID/flook/${f_v}
    retval $? "flook make install"
    cd ../../
    rm -rf flook-${f_v}
    echo "Completed installing flook"
    [ -d $ID/flook/${f_v}/lib64 ] && flook_lib=lib64 || flook_lib=lib
else
    echo "flook directory already found."
fi

##########################
# Completed installation #
##########################

echo ""
echo "##########################"
echo "# Completed installation #"
echo "#    of flook package    #"
echo "#  and its dependencies  #"
echo "##########################"
echo ""
echo ""

echo "Please add the following to the BOTTOM of your arch.make file"
echo ""
echo "INCFLAGS += -I$ID/flook/${f_v}/include"
echo "LDFLAGS += -L$ID/flook/${f_v}/$flook_lib -Wl,-rpath,$ID/flook/${f_v}/$flook_lib"
echo "LIBS += -lflookall -ldl"
echo "COMP_LIBS += libfdict.a"
echo "FPPFLAGS += -DSIESTA__FLOOK"
echo ""
