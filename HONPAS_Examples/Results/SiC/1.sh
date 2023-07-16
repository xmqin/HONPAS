for a in  4.22  4.24  4.26  4.28  4.3  4.32  4.34  4.36  4.38  4.40  4.42  4.44  4.46  4.48
do
cd $a

Volume=`grep "siesta: Cell volume" SiC.out |tail -1 | awk '{printf "%12.6f \n", $5 }'`
E_eV=`grep "Total" SiC.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Force=`grep -A2 "siesta: Atomic forces" SiC.out | tail -1 | awk '{printf "%12.6f \n", $3 }'`

/public/home/xmqin/Paper/honpas_v1.2/Util/Bands/gnubands < SiC.bands > band_plot

LUMO=`sed -n "826p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
HOMO=`sed -n "902p"  band_plot | awk '{printf "%12.6f \n", $2 }'`
Eg=`echo $LUMO - $HOMO |bc -l`

echo " a-V-Eg " " " $a  " " $Volume " "  $E_eV  " "  $Eg >> ../E_a.dat
echo $a $Volume  $E_eV >> ../EOS.dat

cd ..

done

