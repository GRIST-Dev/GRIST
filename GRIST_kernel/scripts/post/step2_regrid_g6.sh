for res in G9 ;do
pathin=prectDiag_${res}UR/
pathout=prectDiag_${res}UR/grid0_5/

mkdir -p ${pathout}

for mon in 01 02 ;do
if [[ ${mon} -eq 01 ]];then
d1=20
d2=31
fi
if [[ ${mon} -eq 02 ]];then
d1=01
d2=28
fi

for day in `seq -w ${d1} ${d2}` ;do
cp ../../UR-grid-latlon/UR-gridInfo/${res}.latlonCellConer.nc ${pathout}/
ncrename -d location_nv,ncells ${pathin}/prectDiag.GRIST.${res}UR.2020-${mon}-${day}.1d.nc ${pathout}/tmp.nc 
ncatted -O -a coordinates,prectDiag,m,c,'lon lat' ${pathout}/tmp.nc
ncks -v prectDiag ${pathout}/tmp.nc ${pathout}/tmp1.nc
ncks ${pathout}/tmp1.nc ${pathout}/${res}.latlonCellConer.nc<<EOF
a
EOF
cdo remap,global_0.5,weight/${res}to0.5w.nc ${pathout}/${res}.latlonCellConer.nc ${pathout}/prectDiag.GRIST.${res}UR.2020-${mon}-${day}.grid.nc
rm -rf ${pathout}/${res}.latlonCellConer.nc ${pathout}/tmp.nc ${pathout}/tmp1.nc 
done

done
done
