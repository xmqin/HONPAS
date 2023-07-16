for a in  3.52  3.54  3.56  3.58  3.60  3.62  3.64  3.66  3.68  3.7  3.72
do
cd $a
/public/home/xmqin/Paper/honpas_2023/honpas_v1.2/Util/Bands/gnubands < BN.bands > band_plot
E_eV=`grep "Total" BN.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" BN.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`
Volume=`grep "siesta: Cell volume" BN.out |tail -1 | awk '{printf "%12.6f \n", $5 }'`

LUMO=`sed -n "909p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
HOMO=`sed -n "826p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
Eg=`echo $LUMO - $HOMO |bc -l`

echo " a-V-Eg " " " $a  " " $Volume " "  $E_eV  " "  $Eg #>> ../E_a.dat
echo $Volume  $E_eV >> ../EOS.dat

cd ..
done

