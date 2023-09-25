rawfile=grist.initial.pl.G6_20210718.nc

# The generated var has been checked by "cdo diff" that is consistent with that in the rawfile
for varname in SoilTemp SoilMoist ;do
echo ${varname}
#for tmp1
for inum in 1 2 3 4 ;do
 ncks -v ${varname}_lv${inum} ${rawfile} ${varname}_lv${inum}.tmp1.nc
done
#for tmp2
 ncrename -d depth,nSoilLevels -v depth,nSoilLevels -v ${varname}_lv1,${varname} ${varname}_lv1.tmp1.nc  ${varname}_lv1.tmp2.nc
for inum in  2 3 4 ;do
 ncrename -d depth_${inum},nSoilLevels -v depth_${inum},nSoilLevels -v ${varname}_lv${inum},${varname} ${varname}_lv${inum}.tmp1.nc ${varname}_lv${inum}.tmp2.nc
done
#for tmp3&4
for inum in 1 2 3 4 ;do
 ncks --mk_rec_dmn nSoilLevels ${varname}_lv${inum}.tmp2.nc ${varname}_lv${inum}.tmp3.nc
 ncks -C -O -x -v depth_bnds,lon_bnds,lat_bnds ${varname}_lv${inum}.tmp3.nc ${varname}_lv${inum}.tmp4.nc 
done
 ncrcat ${varname}_lv1.tmp4.nc ${varname}_lv2.tmp4.nc ${varname}_lv3.tmp4.nc ${varname}_lv4.tmp4.nc ${varname}.tmp.nc 
#unlimited to fixed dim
 ncks --fix_rec_dmn nSoilLevels ${varname}.tmp.nc ${varname}.fixed.nc 
#swap dim
 ncpdq -a ncells,nSoilLevels ${varname}.fixed.nc ${varname}.nc
rm -rf *tmp*.nc ${varname}.fixed.nc

done

if true ;then
#cut a new base file
 ncks -C -O -x -v depth_bnds,lon_bnds,lat_bnds,depth,depth_2,depth_3,depth_4,depth_bnds,depth_2_bnds,depth_3_bnds,depth_4_bnds,SoilTemp_lv1,SoilTemp_lv2,SoilTemp_lv3,SoilTemp_lv4,\
SoilMoist_lv1,SoilMoist_lv2,SoilMoist_lv3,SoilMoist_lv4 ${rawfile} newbase.nc 

ncks SoilTemp.nc newbase.nc<<EOF
a
EOF
ncks SoilMoist.nc newbase.nc<<EOF
a
EOF
mv newbase.nc ${rawfile}.new.nc
rm -rf SoilTemp.nc SoilMoist.nc ${rawfile}
fi
