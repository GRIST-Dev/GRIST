load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "putnc_1d.ncl"

begin

        pi = 3.1416926535
        hres = "g7"

        topo_filein = addfile("~/grist/run/MeshFiles/topo_data/topo_out_"+hres+".nc","r")
        sst_filein  = addfile("../sst.mon.mean.1870_2019.nc","r")
; grist data
        lato        = topo_filein->lat
        lono        = topo_filein->lon
        lato        = lato/(pi/180.)
        lono        = lono/(pi/180.)
        landfrac    = topo_filein->LANDFRAC
; grid data
        sst         = sst_filein->sst(:,::-1,:)
        lati        = sst_filein->latitude(::-1)
        loni        = sst_filein->longitude
        sst         = sst+273.15
; interpolated
        sst_ug      = linint2_points(loni,lati,sst,True,lono,lato,0)
        lf2d        = conform(sst_ug,landfrac,1)
; further mask based on landfrac, if full land, sst must be fillvalue
        sst_ug      = where(lf2d.eq.1.,sst_ug@_FillValue,sst_ug)
        sst_ug!0    = "time"
        sst_ug!1    = "ncol"
        sst_ug&time = sst&time
        sst_ug&ncol = lato
        sst_ug@coordinates = "lon lat"
        lato@long_name="lat"
        lono@long_name="lon"

        putnc1d(sst_ug,lato,lono,"sst","sst.mon.grist."+hres+".nc")

        wks                            = gsn_open_wks("eps","sst2")
        gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
        res                            = True
        plot                           = new(2,graphic)
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

        ;res@cnLevelSelectionMode       = "ManualLevels"
        ;res@cnLevelSpacingF            =
        ;res@cnMinLevelValF             =
        ;res@cnMaxLevelValF             =
        
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
        plot(0)                        = gsn_csm_contour_map(wks, sst_ug(10,:)  ,res)
        plot(1)                        = gsn_csm_contour_map(wks, landfrac  ,res)

        resp                           = True
        resp@gsnPaperOrientation       = "portrait"
        resp@gsnMaximize               = True
        resp@txString                  = ""
        resp@gsnPanelRowSpec           = True
        resp@lbLabelFontHeightF        = 0.01
        resp@gsnPanelXWhiteSpacePercent= 5
        resp@gsnPanelYWhiteSpacePercent= 5
        gsn_panel(wks,plot,(/1,1/),resp)





end
