res=G8X16L4EA
pathou='../raw/'
lev_type=pl

for year in 2021 ;do
for mon in 07 ;do
for day in 15 16 17 18 19 20 ;do

pathin=../download/netcdf/${year}${mon}${day}/
echo ${year} ${mon} ${day} 

if true ;then
 cdo chname,var130,T,var131,U,var132,V,var133,Q ${pathin}/ERA5.${lev_type}.${year}${mon}${day}.00.grib.${res}.nc ${pathou}/initial_${res}_${lev_type}_${year}${mon}${day}.nc
 cdo chname,var134,PS,var129,SOILH,var235,SKINTEMP,var31,XICE,\
var39,SoilMoist_lv1,var40,SoilMoist_lv2,var41,SoilMoist_lv3,var42,SoilMoist_lv4,\
var139,SoilTemp_lv1,var170,SoilTemp_lv2,var183,SoilTemp_lv3,var236,SoilTemp_lv4,\
var33,SNOW,var141,SNOWH \
${pathin}/ERA5.sf.${year}${mon}${day}.00.grib.${res}.nc    ${pathou}/initial_${res}_sf_${year}${mon}${day}.tmp.nc

fi

if  true ; then
ncks -v SNOW          ${pathou}/initial_${res}_sf_${year}${mon}${day}.tmp.nc ${pathou}/snow.nc
ncks -v SNOWH         ${pathou}/initial_${res}_sf_${year}${mon}${day}.tmp.nc ${pathou}/snowh.nc
ncks -v SOILH         ${pathou}/initial_${res}_sf_${year}${mon}${day}.tmp.nc ${pathou}/soilh.nc
ncks -x -v SNOW,SOILH ${pathou}/initial_${res}_sf_${year}${mon}${day}.tmp.nc ${pathou}/newbase.nc 

cdo mul ${pathou}/snow.nc      ${pathou}/snowh.nc ${pathou}/snow1.nc
cdo expr,'SOILH=SOILH/9.80616' ${pathou}/soilh.nc ${pathou}/soilh1.nc 
mv ${pathou}/snow1.nc  ${pathou}/snow.nc
mv ${pathou}/soilh1.nc ${pathou}/soilh.nc

ncks ${pathou}/snow.nc ${pathou}/newbase.nc<<EOF
a
EOF
ncks ${pathou}/soilh.nc ${pathou}/newbase.nc<<EOF
a
EOF

mv ${pathou}/newbase.nc ${pathou}/initial_${res}_sf_${year}${mon}${day}.nc

rm -rf initial_${res}_sf_${year}${mon}${day}.tmp.nc ${pathou}/snow.nc ${pathou}/snowh.nc ${pathou}/soilh.nc

fi

done
done
done
