#!/bin/bash

# bash-script to create tar file for distribution

# create the fdict file
cd fdict
fdict=`./dist.sh`
cd ../

# Get tag branch and number of commits since
describe=`git describe HEAD`
# Remove hash
describe=${describe%-*}

# Save file name
file=ncdf-$describe.tar.gz
rm -f $file

# Create the archive (with prefix)
git archive --prefix ncdf-$describe/ \
    --format tar.gz \
    -o $file HEAD

# Ensure correct directory of the fdict library
# Here we need to extract the created archive,
# and put in the fdict library
mkdir .tmp
cd .tmp
tar xfz ../$file
cd ${file//.tar.gz/}/
# Clean fdict directory
rm -rf fdict
# Extract fdict
tar xfz ../../fdict/$fdict
mv ${fdict//.tar.gz/} fdict
cd ../
tar cfz $file ${file//.tar.gz/}
mv $file ../
cd ../
rm -rf .tmp


# echo file name created
echo $file
