for a in 4.46 4.48 4.50 4.52 4.54 4.56 4.58 4.60 4.62
do
mkdir $a
cp  *psf *fdf $a
cd $a
cat > coord.fdf <<EOF
LatticeConstant      $a Ang

%block LatticeVectors
  0.00  0.50  0.50
  0.50  0.00  0.50
  0.50  0.50  0.00
%endblock LatticeVectors

AtomicCoordinatesFormat     Fractional
%block AtomicCoordinatesAndAtomicSpecies
  0.00  0.00  0.00        1
  0.25  0.25  0.25        2
%endblock AtomicCoordinatesAndAtomicSpecies

EOF
mpirun -np 12 ~/honpas_v1.0/Obj/honpas < BP.fdf |tee BP.out
wait

E_eV=`grep "Total" BP.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" BP.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`

echo "Tol" $a  $E_eV  $Force >> ../E_a.dat
cd ..

done

