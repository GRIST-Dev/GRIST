filename=$1
ncell=$2
date=$3

ncks -v SKINTEMP,XICE ${filename} sst_sic.nc
cdo duplicate,12 sst_sic.nc sst_sic_12.nc
ncrename -v XICE,sic -v SKINTEMP,sst sst_sic_12.nc realNoMissGFSSstSic${date}.cycle.grist.${ncell}.nc
rm -rf sst_sic.nc sst_sic_12.nc
