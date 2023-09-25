for f in g4 g5 g6 g7 g8 g9 ;do
echo $f

ncrename -v mesh_min_c_dist,mesh_min_edp_dist   -v mesh_max_c_dist,mesh_max_edp_dist   -v mesh_mean_c_dist,mesh_mean_edp_dist \
         -v mesh_min_vtx_dist,mesh_min_edt_dist -v mesh_max_vtx_dist,mesh_max_edt_dist -v mesh_mean_vtx_dist,mesh_mean_edt_dist ../grist.grid_file.${f}.ccvt.0d.nc tmp1.nc
ncatted -h -O -a ,global,d,, tmp1.nc grist.grid_file.${f}.ccvt.0d.nc
rm -rf tmp1.nc

ncks -x -v plg_nr,plg_tg,tri_nr,tri_tg,plg_bc_p,plg_bc_ltln,plg_area,tri_area,edt_len,edp_len ../grist.grid_file.${f}.ccvt.2d.nc tmp.nc

ncks -v plg_area -d ntwo,1 ../grist.grid_file.${f}.ccvt.2d.nc plg_area.nc
ncks -v tri_area -d ntwo,1 ../grist.grid_file.${f}.ccvt.2d.nc tri_area.nc
ncks -v edp_len  -d ntwo,1 ../grist.grid_file.${f}.ccvt.2d.nc edp_leng.nc
ncks -v edt_len  -d ntwo,1 ../grist.grid_file.${f}.ccvt.2d.nc edt_leng.nc

ncrename -d ntwo,none plg_area.nc plg_area.1.nc
ncrename -d ntwo,none tri_area.nc tri_area.1.nc
ncrename -d ntwo,none edp_leng.nc edp_leng.1.nc
ncrename -d ntwo,none edt_leng.nc edt_leng.1.nc

ncks plg_area.1.nc tmp.nc<<EOF
a
EOF
ncks tri_area.1.nc tmp.nc<<EOF
a
EOF
ncks edp_leng.1.nc tmp.nc<<EOF
a
EOF
ncks edt_leng.1.nc tmp.nc<<EOF
a
EOF

ncatted -h -O -a ,global,d,, tmp.nc grist.grid_file.${f}.ccvt.2d.nc
rm -rf tmp.nc ed*_leng*.nc *area*nc

done
