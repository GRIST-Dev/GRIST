res=G8X16L4EA
pathin='../raw/'
pathou='../data/'
lev_type=pl
lvname1=plev
lvname2=nSoilLevels
#if [[ ${lev_type} -eq pl ]];then
#lvname=plev
#fi

for year in 2021 ;do
for mon in 07 ;do
#for day in `seq -w 16 16` ;do
for day in `seq -w 15 20` ;do

echo ${year}${mon}${day}
#2d
ncks -x -v SoilMoist,SoilTemp ${pathin}/gfs_initial_${res}_pl_${year}${mon}${day}.nc ${pathou}/tmp1.nc
ncks    -v SoilMoist,SoilTemp ${pathin}/gfs_initial_${res}_pl_${year}${mon}${day}.nc ${pathou}/tmp2.nc

ncpdq -a ncells,${lvname1} ${pathou}/tmp1.nc ${pathou}/tmp1_1.nc
ncpdq -a ncells,${lvname2} ${pathou}/tmp2.nc ${pathou}/tmp2_1.nc
ncks ${pathou}/tmp1_1.nc ${pathou}/tmp2_1.nc <<EOF
a
EOF
mv ${pathou}/tmp2_1.nc ${pathou}/grist.gfs.initial.pl.${res}_${year}${mon}${day}.00.nc
rm -rf ${pathou}/tmp*nc

done
done
done
