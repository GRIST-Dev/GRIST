ncell=40962
for year in 1970 1980 1990 2000 2010 ;do
echo ${year}

# or first gen weight.nc using this line outside
#cdo  -P 16 genycon,grist_scrip_${ncell}.nc ../rawdata/sstKelvin_sic_raw_1970-2016/input4MIPs.sst.2000.nc weight.${ncell}.nc

cdo  -P 16 remap,grist_scrip_${ncell}.nc,weight.${ncell}.nc ../rawdata/sstKelvin_sic_raw_1970-2016/input4MIPs.sst.${year}.nc grist_sst.${year}.nc
cdo  -P 16 remap,grist_scrip_${ncell}.nc,weight.${ncell}.nc ../rawdata/sstKelvin_sic_raw_1970-2016/input4MIPs.sic.${year}.nc grist_sic.${year}.nc

ncks grist_sic.${year}.nc grist_sst.${year}.nc<<EOF
a
EOF
mv grist_sst.${year}.nc realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc
rm -rf grist_sst.${year}.nc grist_sic.${year}.nc

ncrename -v siconcbcs,sic -v tosbcs,sst realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc 1realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc
nccopy -3 1realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc
rm -rf 1realNoMissingCDOYconSstSic.${year}.grist.${ncell}.nc

done
