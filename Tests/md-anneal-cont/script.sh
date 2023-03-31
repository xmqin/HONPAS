#!/bin/sh
#
# An annealing run followed by a Nose run
#
SIESTA="$1"
echo "Running script with SIESTA=$SIESTA"
#
cp ../../Pseudos/O.psf .
cp ../../Pseudos/H.psf .

cat > Anneal.fdf << EOF
SystemName          Water molecule -- md anneal
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

MD.TypeOfRun         anneal
MD.AnnealOption  Temperature
MD.InitialTemperature 100 K
MD.TargetTemperature 600 K
MD.LengthTimeStep  0.2 fs
MD.TauRelax       25 fs
MD.FinalTimeStep  200
EOF
#
$SIESTA < Anneal.fdf > Anneal.out
cp Anneal.out ..
#
# -------------- Nose job
#
cat > Nose.fdf << EOF
SystemName          Water molecule -- md Nose
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

MD.UseSaveXV    T
MD.TypeOfRun         nose
MD.TargetTemperature 600 K
MD.LengthTimeStep  0.2 fs
MD.Initial.Time.Step      1
MD.Final.Time.Step        100
EOF
#
$SIESTA < Nose.fdf > Nose.out
cp Nose.out ..
#------------------------------------


