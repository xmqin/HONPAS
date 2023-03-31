#!/bin/bash

# Script for easy cycling and checking the OUT.diffs
# produced my the make check in the tests folder.

# Script created by:
#  Nick R. Papior, 2016

# This script loops over all sub-folder files as: */OUT.diffs
# and prints them out (with the file-name)
# Subsequently one may delete or keep the file.

function parse_diff {
    local f=$1
    shift
    if [ ! -e $f ]; then
       return 0
    fi

    # 1. Clear screen
    clear

    # 2. print out file name:
    echo "File: $f"
    printf "%61s | %s\n" "NEW" "OLD"
    echo ""

    # 3. Print out diff
    cat $f

    # 4. a couple of new-lines
    echo ""
    echo ""

    # 5. Ask whether the user will delete or keep the file
    read -n 1 -p "Delete file $f [y/n]: " delete
    case $delete in
	y|Y|yes|YES)
	    rm $f
	    ;;
    esac
}

# If there are any inputs we only use those which are cycle
if [ $# -ne 0 ]; then
    
    while [ $# -gt 0 ]
    do
	f=$1
	shift
	parse_diff $f
    done
    
else
    
    for f in */OUT.diffs
    do
	parse_diff $f
    done
    
fi
