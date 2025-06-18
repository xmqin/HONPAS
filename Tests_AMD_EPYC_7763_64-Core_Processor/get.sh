echo "Systems  "  "   Total   "  "     HFX  "   "     HFX_FORCE"
echo "Systems  "  "   Total   "  "     HFX  "   "     HFX_FORCE" >> Time.dat
for a in Si64 Si128 Si256 Si512 Si1024
do
cd $a
Tot_time=`grep "timer: Total elapsed wall-clock" Sibulk.times | tail -1 | awk '{printf "%12.6f \n", $8 }'`
HFX_time=`grep "HFX "  Sibulk.times | tail -1 | awk '{printf "%12.6f \n", $5 }'`
HFXForce_time=`grep "HFX_FORCE"  Sibulk.times | tail -1 | awk '{printf "%12.6f \n", $5 }'`
ERI_time=`grep "ERI "  Sibulk.times | tail -1 | awk '{printf "%12.6f \n", $5 }'`

echo $a "   " $Tot_time "   " $HFX_time  "   " $HFXForce_time 
echo $a "   " $Tot_time "   " $HFX_time  "   " $HFXForce_time  >> ../Time.dat
cd ..
done

