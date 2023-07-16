for a in  5.34  5.36  5.38  5.40  5.42  5.43  5.44  5.46  5.48  5.50  5.52  5.54  5.56
do
cd $a

Volume=`grep "siesta: Cell volume" Sibulk.out |tail -1 | awk '{printf "%12.6f \n", $5 }'`
E_eV=`grep "Total" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`

/public/home/xmqin/Paper/honpas_v1.2/Util/Bands/gnubands < Sibulk.bands > band_plot

LUMO=`sed -n "826p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
HOMO=`sed -n "902p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
Eg=`echo $LUMO - $HOMO |bc -l`

echo " a-V-Eg " " " $a  " " $Volume " "  $E_eV  " "  $Eg >> ../E_a.dat
echo $a $Volume  $E_eV >> ../EOS.dat

cd ..

done

