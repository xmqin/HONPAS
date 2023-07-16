for a in 3.52 3.53 3.54 3.55 3.56 3.57 3.58 3.59 3.60
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
  0.25  0.25  0.25        1
%endblock AtomicCoordinatesAndAtomicSpecies

EOF
mpirun -np 24 /public/home/xmqin/honpas_2023/honpas_v1.2/Obj/honpas < Diamond.fdf |tee Diamond.out

E_eV=`grep "Total" Diamond.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" Diamond.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`

echo "Tol" $a  $E_eV  $Force >> ../E_a.dat
cd ..

done

