res=G6UR
for var in hfx ;do
fileout=${var}_${res}
mkdir -p ${fileout}

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
for hh  in `seq -w 0 95`  ;do

sec=$(( 10#${hh} * 900 ))
sec_format=`printf "%05d\n" ${sec}`
sec_float=`printf "%05f\n" ${sec}`
echo ${sec_format}
echo ${sec_float}
filename_in=GRIST.ATM.${res}.ampwrf.2020-${mon}-${day}-${sec_format}.1d.h1.nc
filename_ou=GRIST.ATM.${res}.amp.2020-${mon}-${day}-${sec_format}.1d.h1.nc

ncks -3 -v ${var} ../experiments/G6/dmdWinter-20200120-hadv33-hnrk3-vadv3-vnrk3/data/atm/1d/${filename_in} ${fileout}/tmp1.nc
ncecat -O -u time  ${fileout}/tmp1.nc ${fileout}/tmp2.nc
ncap2 -O -s time[time]={$sec_float}  ${fileout}/tmp2.nc ${fileout}/tmp3.nc
ncatted -a units,time,o,c,"seconds since 2020-${mon}-${day}" ${fileout}/tmp3.nc ${fileout}/${var}.${filename_ou}

rm -rf ${fileout}/tmp1.nc ${fileout}/tmp2.nc ${fileout}/tmp3.nc
done

ncrcat ${fileout}/${var}.GRIST.ATM.${res}*nc   ${fileout}/ttt.nc
ncatted -h -O -a ,global,d,, ${fileout}/ttt.nc ${fileout}/${var}.GRIST.${res}.2020-${mon}-${day}.1d.nc
rm -rf ${fileout}/${var}.GRIST.ATM.${res}*nc   ${fileout}/ttt.nc

done
done


done
