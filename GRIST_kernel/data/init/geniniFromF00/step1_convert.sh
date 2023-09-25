cdo_grid_file=grist_scrip_2621442.nc # or scrip grid file
hres="G8UR"

for mon in ${MONTH} ;do
pathin=../../f00/${mon}/netcdf/
#pathin=../../f00/
for day in ${DAYYY} ;do

filenc=W_NAFP_C_KWBC_2021${mon}${day}000000_P_gfs.t00z.pgrb2.0p25.f000.nc

rm -rf ${pathin}/skintemp.nc
rm -rf ${pathin}/skintemp.0.nc
rm -rf ${pathin}/skintemp.L4.nc
rm -rf ${pathin}/skintemp.L4.1.nc
rm -rf ${pathin}/skintemp.L4.2.nc
rm -rf ${pathin}/skintemp.L4.3.nc
rm -rf ${pathin}/soilt.nc
rm -rf ${pathin}/tmp0.nc
rm -rf ${pathin}/tmp1.nc

ncks -v u_2,v_2,t,q,sp,orog,t_2,ci,soilw,sdwe,sde,depth ${pathin}/${filenc} ${pathin}/tmp00.nc
ncwa -a time ${pathin}/tmp00.nc ${pathin}/tmp0.nc
rm -rf ${pathin}/tmp00.nc

# set missing to 0, only 4 land vars have missing values (soilmoist, snow,snowh, soiltemp), soiltemp will be specially handled, so not included here
cdo setmisstoc,0            ${pathin}/tmp0.nc ${pathin}/tmp1.nc
# special handl soil temp, set its missing to skintemp
# cut variable
ncks -v st   ${pathin}/${filenc}  ${pathin}/soilt1.nc
ncwa -a time ${pathin}/soilt1.nc  ${pathin}/soilt.nc
rm -rf ${pathin}/soilt1.nc

# remove fillvalue
ncatted -O -a missing_value,,d,,  ${pathin}/soilt.nc
ncatted -O -a _FillValue,,d,,     ${pathin}/soilt.nc
# cut ts
ncks -v t_2 ${pathin}/${filenc}   ${pathin}/skintemp.nc
# rename to soil temp
ncrename -v t_2,st          ${pathin}/skintemp.nc      ${pathin}/skintemp.0.nc
cdo duplicate,4             ${pathin}/skintemp.0.nc    ${pathin}/skintemp.L4.nc
ncks --fix_rec_dmn time     ${pathin}/skintemp.L4.nc   ${pathin}/skintemp.L4.1.nc
ncrename -d time,depth      ${pathin}/skintemp.L4.1.nc ${pathin}/skintemp.L4.2.nc
ncks -x -v time             ${pathin}/skintemp.L4.2.nc ${pathin}/skintemp.L4.3.nc
# missing is -9e33 
cdo ifthenelse -ltc,-1000000  ${pathin}/soilt.nc       ${pathin}/skintemp.L4.3.nc ${pathin}/soilt.nc ${pathin}/soiltNoMiss.nc
ncks ${pathin}/soiltNoMiss.nc ${pathin}/tmp1.nc <<EOF
a
EOF

rm -rf ${pathin}/skintemp.nc
rm -rf ${pathin}/skintemp.0.nc
rm -rf ${pathin}/skintemp.L4.nc
rm -rf ${pathin}/skintemp.L4.1.nc
rm -rf ${pathin}/skintemp.L4.2.nc
rm -rf ${pathin}/skintemp.L4.3.nc
rm -rf ${pathin}/soilt.nc
rm -rf ${pathin}/soiltNoMiss.nc

echo "2) convert lat-lon to unstructured"
cdo -P 6 remapycon,${cdo_grid_file} ${pathin}/tmp1.nc ${pathin}/gdas1.fnl0p25.2021${mon}${day}00.f00.${hres}.nc

echo "3) clean"
#rm -rf ${pathin}/${filenc}
rm -rf ${pathin}/tmp0.nc 
rm -rf ${pathin}/tmp1.nc

echo "done"

#fi
done
done

#done
