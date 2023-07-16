for a in 3.52 3.54 3.56 3.58 3.6
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
mpirun -np 24 /public/home/xmqin/Paper/honpas_v1.2/Obj/honpas < BN.fdf |tee BN.out

/public/home/xmqin/Paper/honpas_v1.2/Util/Bands/gnubands < BN.bands > band_plot
E_eV=`grep "Total" BN.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" BN.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`
Volume=`grep "siesta: Cell volume" BN.out |tail -1 | awk '{printf "%12.6f \n", $5 }'`

LUMO=`sed -n "909p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
HOMO=`sed -n "826p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
Eg=`echo $LUMO - $HOMO |bc -l`

echo " a-V-Eg " " " $a  " " $Volume " "  $E_eV  " "  $Eg >> ../E_a.dat
echo $Volume  $E_eV >> ../EOS.dat

cd ..

done

