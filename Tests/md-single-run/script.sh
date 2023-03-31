#!/bin/sh
#
# A single run of Verlet molecular dynamics
# The job runs for 80 steps.
# Compare with the sequence of jobs in test md-cont.
#
SIESTA="$1"
echo "Running script with SIESTA=$SIESTA"
#
cp ../../Pseudos/O.psf .
cp ../../Pseudos/H.psf .

cat > Job.fdf << EOF
SystemName          Water molecule -- md verlet
SystemLabel         h2o
NumberOfAtoms       3
NumberOfSpecies     2

MeshCutoff  100 Ry

%block ChemicalSpeciesLabel
 1  8  O      # Species index, atomic number, species label
 2  1  H
%endblock ChemicalSpeciesLabel

LatticeConstant 8.0 Ang
%block LatticeVectors
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 0.8
%endblock LatticeVectors

AtomicCoordinatesFormat  Ang
%block AtomicCoordinatesAndAtomicSpecies
 0.000  0.000  0.000  1
 0.757  0.586  0.000  2
-0.757  0.586  0.000  2
%endblock AtomicCoordinatesAndAtomicSpecies

Solution.Method       diagon
MeshCutoff             100 Ry

WriteCoorStep      .true.
WriteForces        .true.
WriteMDHistory     .true.

MD.UseSaveXV       T
MD.TypeOfRun         Verlet
MD.InitialTemperature 600 K
MD.Initial.Time.Step      1
MD.Final.Time.Step        80
MD.Length.Time.Step       0.2 fs
EOF
#
$SIESTA < Job.fdf > md-single-run.out
cp md-single-run.out ..


