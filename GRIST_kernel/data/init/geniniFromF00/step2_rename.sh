res=G8UR
pathou='../raw/'
lev_type=pl

for year in 2021 ;do
for mon in ${MONTH} ;do
pathin=../../f00/${mon}/netcdf/
for day in ${DAYYY} ;do

echo ${year} ${mon} ${day} 

if true ;then
#cdo change name
cdo chname,t,T,u_2,U,v_2,V,q,Q,sp,PS,orog,SOILH,t_2,SKINTEMP,ci,XICE,soilw,SoilMoist,st,SoilTemp,sdwe,SNOW,sde,SNOWH \
${pathin}/gdas1.fnl0p25.${year}${mon}${day}00.f00.${res}.nc ${pathou}/tmp.nc

# calculaete mid-level soil level
# cdo expr,'lv_DBLL11_l0=(lv_DBLL11_l0+lv_DBLL11_l1)/2;' ${pathou}/tmp.nc ${pathou}/tmp1.nc
# cdo cannot change dim name, use ncrename
# ncrename -d lv_DBLL11,nSoilLevels -v lv_DBLL11_l0,nSoilLevels ${pathou}/tmp1.nc ${pathou}/tmp2.nc
# ncrename -d lv_ISBL0,plev -d lv_DBLL11,nSoilLevels -v lv_ISBL0,plev ${pathou}/tmp.nc ${pathou}/tmp3.nc
# ncks tmp2.nc tmp3.nc<<EOF
#a
#EOF
ncrename -d depth,nSoilLevels ${pathou}/tmp.nc ${pathou}/tmp3.nc
# cut
ncks -v U,V,T,Q,PS,SKINTEMP,XICE,SoilMoist,SoilTemp,SNOW,SNOWH,SOILH ${pathou}/tmp3.nc ${pathou}/gfs_fnl00_${res}_${lev_type}_${year}${mon}${day}.nc
#ncks ${pathou}/tmp2.nc ${pathou}/gfs_initial_${res}_${lev_type}_${year}${mon}${day}.nc <<EOF
#a
#EOF
rm -rf ${pathou}/tmp.nc ${pathou}/tmp1.nc ${pathou}/tmp2.nc ${pathou}/tmp3.nc
fi

done
done
done
