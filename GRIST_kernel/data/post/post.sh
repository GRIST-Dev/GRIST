cp ${filehead}.1d.nc 1d.nc.
ncks -d ntracer,0 ${filehead}.3d.nc a.nc;
ncwa -a ntracer a.nc 3da.nc;
cp ${filehead}.2d.nc 2d.nc;
ncks 3da.nc 2d.nc<<EOF
a
EOF
ncks 2d.nc 1d.nc <<EOF
a
EOF
cdo remapdis,global_1 1d.nc ${filehead}.grid.nc;
rm -rf 1d.nc 2d.nc 3da.nc a.nc;
