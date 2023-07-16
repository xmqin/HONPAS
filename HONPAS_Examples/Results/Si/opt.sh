for a in 5.43 # 5.40 5.41 5.42 5.43 5.44 5.442 5.444 5.446 5.448 5.452 5.454 5.46 5.47 5.48 
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
mpirun -np 24 /public/home/xmqin/Paper/honpas_v1.2/Obj/honpas < Sibulk.fdf |tee Sibulk.out

E_eV=`grep "Total" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`

echo "Tol" $a  $E_eV  $Force >> ../E_a.dat
cd ..

done

