for a in 64 128 256 512 1024
do
cd $a

ERI=`grep " HFX time : " Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
Allreduce=`grep "Allreduce time" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $4 }'`
HFX_time=`grep "All HFX time " Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $5 }'`
Force1=`grep " HFX Force time" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $5 }'`
Force2=`grep " All Force time" Sibulk.out | tail -1 | awk '{printf "%12.6f \n", $5 }'`

echo $a " " $ERI " " $Allreduce " " $HFX_time " " $Force1 " " $Force2
echo $a " "  $ERI " " $Allreduce " " $HFX_time " " $Force1 " " $Force2 >> ../Time.dat
cd ..

done

