#!/bin/sh
#
psfile=$1
#
if [ $# == 2 ]
then
      pagenum=$2
      echo "Extracting page $2 of file $psfile"
else
      pagenum=_1
      echo "Extracting last page of file $psfile"
fi
#
rm -rf tmp__aux.ps
#
# Pick last page of braindead GNUPLOT output
#
psselect -p${pagenum} $psfile tmp__aux.ps
#
name=`basename $psfile .ps`
#
# Convert to PDF
#
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$name.pdf tmp__aux.ps
#
rm -rf tmp__aux.ps
#
echo " ==> Generated $name.pdf"
