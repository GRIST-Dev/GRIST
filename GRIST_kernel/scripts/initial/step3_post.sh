res=5
pathin='../../swv/20150816/raw/'
pathou='../../swv/20150816/data/'
lev_type=ml
lvname=lev
#if [[ ${lev_type} -eq pl ]];then
#lvname=plev
#fi

for year in 2015 ;do
for mon in 08 ;do
for day in `seq -w 16 16` ;do

echo ${year}${mon}${day}
#2d
ncks -d time,0 ${pathin}/initial_g${res}_${lev_type}_${year}${mon}${day}.nc  initial_g${res}.dim1.nc
ncwa -a time initial_g${res}.dim1.nc tmp.nc
ncpdq -a ncells,${lvname} tmp.nc ${pathou}/grist.initial.${lev_type}.g${res}_${year}${mon}${day}.nc
rm -rf initial_g${res}.dim1.nc tmp.nc

#1d
ncks -d time,0 ${pathin}/initial_g${res}_sf_${year}${mon}${day}.nc  initial_g${res}.dim1.nc
ncwa -a time initial_g${res}.dim1.nc grist.initial.g${res}_sf_${year}${mon}${day}.nc

#append
ncks grist.initial.g${res}_sf_${year}${mon}${day}.nc ${pathou}/grist.initial.${lev_type}.g${res}_${year}${mon}${day}.nc<<EOF
a
EOF
rm -rf initial_g${res}.dim1.nc grist.initial.g${res}_sf_${year}${mon}${day}.nc

done
done
done
