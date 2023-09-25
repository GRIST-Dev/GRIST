pathin='../../swv/20150816/download/'

cdo_grid_file=g5_file_prime_cell.cdo

for file in `ls ${pathin}` ;do

if [ "${file##*.}"x = "grib"x ] ;then

echo ${file}

echo "1) convert grib to netcdf"
cdo -f nc copy ${pathin}/${file} ${pathin}/${file}.tmp.nc
echo "2) convert lat-lon to unstructured"
cdo remapycon,${cdo_grid_file} ${pathin}/${file}.tmp.nc ${pathin}/${file}.nc
echo "3) clean"
rm -rf ${pathin}/${file}.tmp.nc
echo "done"

fi

done
