res=G8X16L4EA
pathin='../raw/'
pathou='../G8X16L4EA/'
lev_type=pl
lvname=plev
#if [[ ${lev_type} -eq pl ]];then
#lvname=plev
#fi

for year in 2021 ;do
for mon in 07 ;do
#for day in `seq -w 16 16` ;do
for day in 20 ;do

echo ${year}${mon}${day}
#2d
ncks -d time,0 ${pathin}/initial_${res}_${lev_type}_${year}${mon}${day}.nc  initial_${res}.dim1.nc
ncwa -a time initial_${res}.dim1.nc tmp.nc
ncpdq -a ncells,${lvname} tmp.nc ${pathou}/grist.era5.ini.${lev_type}.${res}_${year}${mon}${day}.nc
rm -rf initial_${res}.dim1.nc tmp.nc

#1d
ncks -d time,0 ${pathin}/initial_${res}_sf_${year}${mon}${day}.nc  initial_${res}.dim1.nc
ncwa -a time initial_${res}.dim1.nc grist.initial.${res}_sf_${year}${mon}${day}.nc

#append
ncks grist.initial.${res}_sf_${year}${mon}${day}.nc ${pathou}/grist.era5.ini.${lev_type}.${res}_${year}${mon}${day}.nc<<EOF
a
EOF
rm -rf initial_${res}.dim1.nc grist.initial.${res}_sf_${year}${mon}${day}.nc

done
done
done
