for mon in 06 ;do
#for day in `seq -w 01 31` ;do
for day in 01 ;do

export MONTH=${mon}
export DAYYY=${day}

./step1_convert.sh  
./step2_rename.sh  
./step3_post.sh

rm -rf ../../f00/${mon}/netcdf/gdas1.fnl0p25.2021${mon}${day}00.f00.G8UR.nc
rm -rf ../raw/gfs_fnl00_G8UR_pl_2021${mon}${day}.nc
echo "done"

done
done
