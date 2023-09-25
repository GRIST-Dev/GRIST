for day in `seq -w 15 20` ;do
pathin=../download/grib/202107${day}/
pathou=../download/netcdf/202107${day}/
mkdir -p ${pathou}

hres="G8X16L4EA"
cdo_grid_file=remapFile/grist_scrip_655362.G8X16L4EA.nc # or scrip grid file

for file in `ls ${pathin}` ;do

if [ "${file##*.}"x = "grb2"x ] ;then

echo ${file}

rm -rf ${pathou}/skintemp.nc ${pathou}/skintemp.0.nc ${pathou}/skintemp.L4.nc ${pathou}/skintemp.L4.1.nc ${pathou}/skintemp.L4.2.nc ${pathou}/skintemp.L4.3.nc ${pathou}/soilt.nc
rm -rf ${pathou}/tmp0.nc ${pathou}/tmp1.nc

echo "1) convert grib to netcdf"
ncl_convert2nc ${pathin}/${file}
filenc=`ls gfs_*.nc`
mv ${filenc} ${pathou}

ncks -v TMP_P0_L100_GLL0,UGRD_P0_L100_GLL0,VGRD_P0_L100_GLL0,SPFH_P0_L100_GLL0,PRES_P0_L1_GLL0,HGT_P0_L1_GLL0,TMP_P0_L1_GLL0,ICEC_P0_L1_GLL0,SOILW_P0_2L106_GLL0,WEASD_P0_L1_GLL0,SNOD_P0_L1_GLL0,lv_DBLL11_l0,lv_DBLL11_l1 ${pathou}/${filenc} ${pathou}/tmp0.nc 
# set missing to 0, only 4 land vars have missing values (soilmoist, snow,snowh, soiltemp)
cdo setmisstoc,0            ${pathou}/tmp0.nc ${pathou}/tmp1.nc

# special handl temp
ncks -v TSOIL_P0_2L106_GLL0 ${pathou}/${filenc}       ${pathou}/soilt.nc
ncatted -O -a _FillValue,,d,, ${pathou}/soilt.nc

ncks -v TMP_P0_L1_GLL0      ${pathou}/${filenc}        ${pathou}/skintemp.nc
ncrename -v TMP_P0_L1_GLL0,TSOIL_P0_2L106_GLL0         ${pathou}/skintemp.nc ${pathou}/skintemp.0.nc
cdo duplicate,4             ${pathou}/skintemp.0.nc    ${pathou}/skintemp.L4.nc
ncks --fix_rec_dmn time     ${pathou}/skintemp.L4.nc   ${pathou}/skintemp.L4.1.nc
ncrename -d time,lv_DBLL11  ${pathou}/skintemp.L4.1.nc ${pathou}/skintemp.L4.2.nc
ncks -x -v time ${pathou}/skintemp.L4.2.nc             ${pathou}/skintemp.L4.3.nc
# missing is 1e20
cdo ifthenelse -gtc,10000000  ${pathou}/soilt.nc       ${pathou}/skintemp.L4.3.nc ${pathou}/soilt.nc ${pathou}/soiltNoMiss.nc
ncks ${pathou}/soiltNoMiss.nc ${pathou}/tmp1.nc <<EOF
a
EOF
rm -rf ${pathou}/skintemp.nc ${pathou}/skintemp.0.nc ${pathou}/skintemp.L4.nc ${pathou}/skintemp.L4.1.nc ${pathou}/skintemp.L4.2.nc ${pathou}/skintemp.L4.3.nc ${pathou}/soilt.nc

echo "2) convert lat-lon to unstructured"
cdo -P 6 remapycon,${cdo_grid_file} ${pathou}/tmp1.nc ${pathou}/${filenc}.${hres}.nc

echo "3) clean"
rm -rf ${pathou}/${filenc}  ${pathou}/tmp0.nc ${pathou}/tmp1.nc
echo "done"

fi

done

done
