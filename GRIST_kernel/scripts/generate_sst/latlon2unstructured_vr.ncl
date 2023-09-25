load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "putnc_1d.ncl"
;
; This version can handle both qu and vr MeshFiles
; Yi Zhang
; 2021-01-25: add sic from input4mip
;
begin

        land_mask   = False
        pi          = 3.14159265358979323846d0 ; using this value, generate the same results as using static.nc
        hres        = "g6"
        fillvalue   = 271.35
        year        = 0000

        ;topo_filein = addfile("../../bdy/new/static.nc","r")
        sst_filein  = addfile("sst_1978_2010.climo.nc","r")
        sic_filein  = addfile("sic_1978_2010.climo.nc","r")

        mesh_filein = addfile("../grist/inputdata/mesh/uniform-pol/grist.grid_file.g6.ccvt.2d.nc","r")
        vtx_nnb     = mesh_filein->vtx_ltln_nnb(2,:)
        lono        = mesh_filein->vtx_ltln_nnb(1,:)
        lato        = mesh_filein->vtx_ltln_nnb(0,:)
        lono        = lono/(pi/180.)
        lato        = lato/(pi/180.)
        vtx_nb      = mesh_filein->vtx_nb
        lato!0      = "ncol"
        lono!0      = "ncol"

;----------------------------
; use grist land data's latlon
        ;lato        = topo_filein->lat
        ;lono        = topo_filein->lon
        ;lato        = lato
        ;lono        = lono
;----------------------------

        lono        = where(lono.lt.0,lono+360.,lono)

        ;landfrac    = topo_filein->MASK
; amip sst and grid data
        sst         = sst_filein->tosbcs
        sst         = sst+273.15
        lati        = sst_filein->lat
        loni        = sst_filein->lon
        sst_ug      = linint2_points(loni,lati,sst,True,lono,lato,0)
        sst_ug@_FillValue = -1e30

        sic         = sic_filein->siconcbcs
        lat2        = sic_filein->lat
        lon2        = sic_filein->lon
        sic_ug      = linint2_points(lon2,lat2,sic,True,lono,lato,0)
        sic_ug@_FillValue = -1e30
;
; interpolated, this introduces fillvalue?? why
; because lato lono are outside lati loni
; this fillvalue will be broungt to model even under the land mask judgement
; to avoid this, give missing value its average of nbrs
;
; SST
        do it = 0, 11 ; 12 mon
           index_a = ind(ismissing(sst_ug(it,:))) ; ncl index for missing values of sst_ug
           print("SST: num of missing is"+dimsizes(index_a))
           do iv = 0, dimsizes(index_a)-1
              index_nb1 = doubletointeger(vtx_nb(:,index_a(iv)))   ; fortran index
;print(min(index_nb1))
;print(max(index_nb1))
              ;index_nb1@_FillValue = -9   ; missing in mpi-scvt
              index_nb1@_FillValue = -1000000 ; missing in gg
              index_nb1 = index_nb1-1 ; ncl index, 0 (7th) will become -1
              value = 0.
              count = 0
              do icell = 0, dimsizes(index_nb1)-1
; the nb index can not be missing (-1e6) or unfilled (0/-1), and its nb value cannot be missing
                 if(.not.ismissing(index_nb1(icell)) .and. index_nb1(icell).ne.-1. .and. .not.ismissing(sst_ug(it,index_nb1(icell))))then
                   ; print("icell="+icell+"nb value is"+sst_ug(it,index_nb1(icell)))
                    value = value + sst_ug(it,index_nb1(icell))
                    count = count +1
                 end if
              end do
              if(count.ne.0)then
                 value = value/count
              else
                 value = fillvalue
                 print("final-avg-value is "+value)
              end if
              sst_ug(it,index_a(iv)) = value
            end do
        end do

; SIC
        do it = 0, 11 ; 12 mon
           index_a = ind(ismissing(sic_ug(it,:))) ; ncl index for missing values of sic_ug
           print("SIC: num of missing is"+dimsizes(index_a))
           do iv = 0, dimsizes(index_a)-1
              index_nb1 = doubletointeger(vtx_nb(:,index_a(iv)))   ; fortran index
;print(min(index_nb1))
;print(max(index_nb1))
              ;index_nb1@_FillValue = -9   ; missing in mpi-scvt
              index_nb1@_FillValue = -1000000 ; missing in gg
              index_nb1 = index_nb1-1 ; ncl index, 0 (7th) will become -1
              value = 0.
              count = 0
              do icell = 0, dimsizes(index_nb1)-1
; the nb index can not be missing (-1e6) or unfilled (0/-1), and its nb value cannot be missing
                 if(.not.ismissing(index_nb1(icell)) .and. index_nb1(icell).ne.-1. .and. .not.ismissing(sic_ug(it,index_nb1(icell))))then
                   ; print("icell="+icell+"nb value is"+sic_ug(it,index_nb1(icell)))
                    value = value + sic_ug(it,index_nb1(icell))
                    count = count +1
                 end if
              end do
              if(count.ne.0)then
                 value = value/count
              else
                 value = fillvalue
                 print("final-avg-value is "+value)
              end if
              sic_ug(it,index_a(iv)) = value
            end do
        end do

; 
; further mask based on landfrac, if full land, sst must be fillvalue
;
        if(land_mask)then
            ;lf2d        = conform(sst_ug,landfrac,1)
            ;sst_ug      = where(lf2d.eq.1.,sst_ug@_FillValue,sst_ug)
        end if
;---------------------------------------------------------------------------------

; attributes
        sst_ug!0    = "time"
        sst_ug!1    = "ncol"
        sst_ug&time = sst&time
        sst_ug&ncol = lato
        sst_ug@coordinates = "lon lat"
        lato@long_name="lat"
        lono@long_name="lon"

        copy_VarMeta(sst_ug,sic_ug)

; modify sic data to be 0->1
       sic_ug = sic_ug/100.
       sic_ug = where(sic_ug.lt.0.,0.,sic_ug) ; no   sea ice
       sic_ug = where(sic_ug.gt.1.,1.,sic_ug) ; full sea ice

       putnc1d(sst_ug,sic_ug,lato,lono,"sst","sic","realNoMissingNewSstSic."+year+".grist."+hres+".nc")

        wks                            = gsn_open_wks("eps","sic")
        gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
        res                            = True
        plot                           = new(3,graphic)
        res@gsnSpreadColors            = True
        res@gsnSpreadColorStart        = 2
        res@gsnSpreadColorEnd          = 226

        res@gsnFrame                   = False
        res@gsnDraw                    = False

        res@tmXBMode                   = "Explicit"
        ;res@tmXBValues                 = (//)
        ;res@tmXBLabels                 = (//)  
        res@tmXBLabelFontHeightF       = 0.018

        res@tmYLMode                   = "Explicit"
        ;res@tmYLValues                 = (//)
        ;res@tmYLLabels                 = (//)  
        res@tmYLLabelFontHeightF       = 0.018

        res@gsnCenterStringFontHeightF = 0.018
        res@gsnLeftStringFontHeightF   = 0.018
        res@gsnRightStringFontHeightF  = 0.018

        res@gsnLeftString              = ""
        res@gsnCenterString            = ""
        res@gsnRightString             = ""

        res@tiXAxisFontHeightF         = 0.018
        res@tiXAxisString              = ""
        res@tiYAxisFontHeightF         = 0.018
        res@tiYAxisString              = ""

        res@cnFillOn                   = True
        res@cnLinesOn                  = False
        res@cnLineLabelsOn             = False
        res@cnLineLabelPlacementMode   = "Constant"

        res@cnLevelSelectionMode       = "ManualLevels"
;        res@cnLevelSpacingF            = 3
;        res@cnMinLevelValF             = 270
;        res@cnMaxLevelValF             = 300
        
        res@cnMissingValFillColor      = "gray3"
        res@cnMissingValFillPattern    = 0        ; set to complete fill

        res@gsnAddCyclic               = False
        res@mpCenterLonF               = 180.
        res@mpLimitMode                = "Corners"
        res@mpLeftCornerLatF           =-90.
        res@mpRightCornerLatF          = 90.
        res@mpLeftCornerLonF           = 0.
        res@mpRightCornerLonF          = 360.
        res@sfXArray                   = lono
        res@sfYArray                   = lato
        plot(0)                        = gsn_csm_contour_map(wks, sic_ug(0,:)  ,res)
        plot(1)                        = gsn_csm_contour_map(wks, sic_ug(5,:)  ,res)
        plot(2)                        = gsn_csm_contour_map(wks, sic_ug(11,:)  ,res)

        resp                           = True
        resp@gsnPaperOrientation       = "portrait"
        resp@gsnMaximize               = True
        resp@txString                  = ""
        resp@gsnPanelRowSpec           = True
        resp@lbLabelFontHeightF        = 0.01
        resp@gsnPanelXWhiteSpacePercent= 5
        resp@gsnPanelYWhiteSpacePercent= 5
        gsn_panel(wks,plot,(/1,1,1/),resp)





end
