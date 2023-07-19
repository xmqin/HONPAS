#!/bin/sh

SIESTA="$1"
echo "Running script with SIESTA=$SIESTA"
#
fractional=../../../../Util/VCA/fractional
#
cp ../../Pseudos/O.psf .

if [ ! -x $fractional ] ; then
#  echo "Compiling $fractional..."
#  (cd ../../../../Util/VCA ; make fractional)
   echo -n "Please compile 'fractional' in "
   echo "Util/VCA before running this test"
   exit 1
fi

$fractional O  0.5

#
cat > Job.fdf << EOF
SystemName          Oxypartial
SystemLabel         oxypartial
NumberOfAtoms       2
NumberOfSpecies     1

MeshCutoff 200 Ry

%block ChemicalSpeciesLabel
 1  201  O-Fraction-0.50000      # Species index, atomic number, species label
%endblock ChemicalSpeciesLabel
%block SyntheticAtoms
 1  
 2 2 3 4
 1.0 2.0 0.0 0.0
%endblock SyntheticAtoms

%block PAO.basis
O-Fraction-0.50000  2                 # Species label, number of l-shells
 n=2   0   2                         # n, l, Nzeta 
 0.0 0.0   
   1.000      1.000   
 n=2   1   2 P   1                   # n, l, Nzeta, Polarization, NzetaPol
   4.139      2.740   
   1.000      1.000   
%endblock PAO.Basis

DM.NumberPulay 4

Spin.Polarized T

AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 0.000  0.000  0.000  1
 0.000  0.000  1.200  1
%endblock AtomicCoordinatesAndAtomicSpecies

MD.TypeOfRun Broyden
MD.NumCGSteps 40
EOF
#
$SIESTA < Job.fdf > Job.out
cp Job.out ../partial.out

