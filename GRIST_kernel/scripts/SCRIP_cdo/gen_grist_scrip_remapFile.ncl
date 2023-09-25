load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  f1  = addfile("grist.grid_file.g9.ccvt.2d.nc","r")

  grid_size        = 2621442
  grid_corners     = 6    ; max vtx nb
  grid_rank        = 1
  
  grid_dims        = new((/grid_rank/),"integer")
  grid_dims        = (/grid_size/)
  grid_dims!0      = "grid_rank"

  tri_cc_ltln      = f1->tri_cc_ltln
  vtx_ltln_nnb     = f1->vtx_ltln_nnb
  vtx_nnb          = doubletointeger(vtx_ltln_nnb(2,:))
  vtx_tr           = f1->vtx_tr
;  plg_area         = f1->plg_area
;printVarSummary(vtx_tr)

  grid_center_lat  = (/vtx_ltln_nnb(0,:)/)
  grid_center_lat!0 = "grid_size"
  grid_center_lat@units = "degree"
;  printVarSummary(grid_center_lat) 

  grid_center_lon  = (/vtx_ltln_nnb(1,:)/)
  grid_center_lon!0 = "grid_size"
  grid_center_lon@units = "degree"
;  printVarSummary(grid_center_lon)

  grid_corner_lat  = new((/grid_size,grid_corners/),"double")
  grid_corner_lon  = new((/grid_size,grid_corners/),"double")

  if(grid_corners.eq.6)then
    do i = 0, grid_size-1
;      print("i="+i)
       do j = 0, vtx_nnb(i)-1
          jl = doubletointeger(vtx_tr(j,i))  ; Fortran Index
          jl = jl-1
          grid_corner_lat(i,j) = tri_cc_ltln(0,jl)
          grid_corner_lon(i,j) = tri_cc_ltln(1,jl)
       end do
       if(vtx_nnb(i).eq.5)then
          grid_corner_lat(i,5) = grid_corner_lat(i,4)
          grid_corner_lon(i,5) = grid_corner_lon(i,4)
       end if
    end do
  end if

  if(grid_corners.eq.7)then
    do i = 0, grid_size-1
;      print("i="+i)
       do j = 0, vtx_nnb(i)-1
          jl = doubletointeger(vtx_tr(j,i))  ; Fortran Index
          jl = jl-1
          grid_corner_lat(i,j)   = tri_cc_ltln(0,jl)
          grid_corner_lon(i,j)   = tri_cc_ltln(1,jl)
       end do
       if(vtx_nnb(i).eq.5)then
          grid_corner_lat(i,5:6) = grid_corner_lat(i,4)
          grid_corner_lon(i,5:6) = grid_corner_lon(i,4)
       end if
       if(vtx_nnb(i).eq.6)then
          grid_corner_lat(i,6)   = grid_corner_lat(i,5)
          grid_corner_lon(i,6)   = grid_corner_lon(i,5)
       end if
    end do
  end if

  grid_corner_lat!0          = "grid_size"
  grid_corner_lat!1          = "grid_corners"
  grid_corner_lat@units      = "degree"
  grid_corner_lat@_FillValue = -9999

  grid_corner_lon!0          = "grid_size"
  grid_corner_lon!1          = "grid_corners"
  grid_corner_lon@units      = "degree"
  grid_corner_lon@_FillValue = -9999

 ; grid_imask      = new((/grid_size/),"integer")
 ; grid_imask@_FillValue = -9999
 ; grid_imask!0    = "grid_size"
 ; grid_imask      = 1
 ; printVarSummary(grid_imask)

 ; grid_area           = plg_area(1,:)
 ; grid_area!0         = "grid_size"
 ; grid_area@units     = "radians^2"
 ; grid_area@long_name = "area weight"
 ; printVarSummary(grid_area)

; follow the same as in ASCII CDO file
  pi              = 3.14159265358979323846d0
  grid_corner_lat = grid_corner_lat/(pi/180.d0)
  grid_corner_lon = grid_corner_lon/(pi/180.d0)
  grid_center_lat = grid_center_lat/(pi/180.d0)
  grid_center_lon = grid_center_lon/(pi/180.d0)

  fout                  = addfile("grist_scrip_"+grid_size+".nc","c")
  fout->grid_dims       = grid_dims
  fout->grid_center_lat = grid_center_lat
  fout->grid_center_lon = grid_center_lon
  fout->grid_corner_lat = grid_corner_lat
  fout->grid_corner_lon = grid_corner_lon
  ;fout->grid_imask=grid_imask
  ;fout->grid_area=grid_area

end 
