pathin='../download/grib/20210720/'
pathou='../download/netcdf/20210720/'
mkdir -p ${pathou}

hres="G8X16L4EA"
cdo_grid_file=../../gen_grist_ini_from_gfs/scripts/remapFile/grist_scrip_655362.G8X16L4EA.nc

for file in `ls ${pathin}` ;do

if [ "${file##*.}"x = "grib"x ] ;then

echo ${file}

echo "1) convert grib to netcdf"
cdo -f nc copy ${pathin}/${file} ${pathou}/${file}.tmp0.nc
# only sea ice fraction has missing, just set to 0
cdo setmisstoc,0                 ${pathou}/${file}.tmp0.nc ${pathou}/tmp.nc

echo "2) convert lat-lon to unstructured"
cdo -P 6 remapycon,${cdo_grid_file} ${pathou}/${file}.tmp.nc ${pathou}/${file}.${hres}.nc

echo "3) clean"
rm -rf ${pathou}/${file}.tmp.nc ${pathou}/${file}.tmp0.nc
echo "done"

fi

done
