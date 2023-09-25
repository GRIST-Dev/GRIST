# use this code to generate 12 invariant sst/sic from initial data for NWP run
res=G8X16L4
for day in `seq -w 15 20` ;do
./gen_bdy_sstsic.sh grist.gfs.initial.pl.${res}_202107${day}.00.nc 655362 202107${day}
done
