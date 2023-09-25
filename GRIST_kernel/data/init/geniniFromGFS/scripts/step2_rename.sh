res=G8X16L4EA
pathou='../raw/'
lev_type=pl

for year in 2021 ;do
for mon in 07 ;do
for day in 15 16 17 18 19 20 ;do

pathin=../download/netcdf/${year}${mon}${day}/
echo ${year} ${mon} ${day} 

if true ;then
#cdo change name
cdo chname,TMP_P0_L100_GLL0,T,UGRD_P0_L100_GLL0,U,VGRD_P0_L100_GLL0,V,SPFH_P0_L100_GLL0,Q,PRES_P0_L1_GLL0,PS,HGT_P0_L1_GLL0,SOILH,TMP_P0_L1_GLL0,SKINTEMP,\
ICEC_P0_L1_GLL0,XICE,SOILW_P0_2L106_GLL0,SoilMoist,TSOIL_P0_2L106_GLL0,SoilTemp,WEASD_P0_L1_GLL0,SNOW,SNOD_P0_L1_GLL0,SNOWH \
 ${pathin}/gfs_4_${year}${mon}${day}_0000_000.nc.${res}.nc ${pathou}/tmp.nc
# calculaete mid-level soil level
 cdo expr,'lv_DBLL11_l0=(lv_DBLL11_l0+lv_DBLL11_l1)/2;' ${pathou}/tmp.nc ${pathou}/tmp1.nc
# cdo cannot change dim name, use ncrename
 ncrename -d lv_DBLL11,nSoilLevels -v lv_DBLL11_l0,nSoilLevels ${pathou}/tmp1.nc ${pathou}/tmp2.nc
 ncrename -d lv_ISBL0,plev -d lv_DBLL11,nSoilLevels -v lv_ISBL0,plev ${pathou}/tmp.nc ${pathou}/tmp3.nc
# ncks tmp2.nc tmp3.nc<<EOF
#a
#EOF
# cut
ncks -v U,V,T,Q,PS,SKINTEMP,XICE,SoilMoist,SoilTemp,SNOW,SNOWH,SOILH ${pathou}/tmp3.nc ${pathou}/gfs_initial_${res}_${lev_type}_${year}${mon}${day}.nc
ncks ${pathou}/tmp2.nc ${pathou}/gfs_initial_${res}_${lev_type}_${year}${mon}${day}.nc <<EOF
a
EOF
rm -rf ${pathou}/tmp.nc ${pathou}/tmp1.nc ${pathou}/tmp2.nc ${pathou}/tmp3.nc
fi

done
done
done
