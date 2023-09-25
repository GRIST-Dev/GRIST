res=5
pathin='../../swv/20150816/download/'
pathou='../../swv/20150816/raw/'
lev_type=ml

for year in 2015 ;do
for mon in 08 ;do
for day in `seq -w 16 16` ;do

 echo ${year} ${mon} ${day} 
 cdo chname,var130,T,var131,U,var132,V,var133,Q ${pathin}/ERAI.${lev_type}.${year}${mon}${day}.00.grib.nc ${pathou}/initial_g${res}_${lev_type}_${year}${mon}${day}.nc
 cdo chname,var134,PS,var129,PHIS               ${pathin}/ERAI.sf.${year}${mon}${day}.00.grib.nc ${pathou}/initial_g${res}_sf_${year}${mon}${day}.nc

done
done
done
