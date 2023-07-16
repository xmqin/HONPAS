for a in 4.44  4.46  4.48  4.50  4.52  4.54  4.56  4.58  4.60  4.62  4.64  4.66
do
cd $a
E_eV=`grep "Total" BP.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Volume=`grep "siesta: Cell volume" BP.out |tail -1 | awk '{printf "%12.6f \n", $5 }'`

echo $Volume  $E_eV 

cd ..
done

