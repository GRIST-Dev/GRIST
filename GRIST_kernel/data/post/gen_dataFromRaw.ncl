load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "putnc2d.ncl"

; read a CDO-interpolated GRIST regular-grid data
; do some vertical interpolation
; write a datafile with utc time info for further manipulation

begin

    years = (/2002,2003,2004,2005,2006/)
    mons  = (/6,7,8/)

    do iyear = 0,4
    do imon  = 0,2
        cmon  = sprinti("%0.2i", mons(imon))
        print(cmon+" ")
        if(cmon.eq."06")then
          dayend=30
        else
          dayend=31
        end if

    do iday  = 0, dayend-1
        cday  = sprinti("%0.2i", iday+1)
        print(cday+" ")

        filein   = addfile("../GRIST.ATM.G6.amp."+years(iyear)+"-"+cmon+"-"+cday+".grid.nc","r")

        u3d      = filein->uPC
        prect    = filein->prectDiag
        olr      = filein->flwutDiag
        presFace = filein->mpressureFace
        presFull = 0.5*(presFace(:,:,0:29)+presFace(:,:,1:30))

        copy_VarCoords(u3d,presFull)

        u3din =      u3d(nlev|:,lat|:,lon|:)
        p3din = presFull(nlev|:,lat|:,lon|:)

        u850     = wrf_interp_3d_z(u3din,log(log(p3din)),log(log(85000.)))

        u850_out = new((/1,dimsizes(u3din&lat),dimsizes(u3din&lon)/),typeof(u3din),u850@_FillValue)
        prec_out = new((/1,dimsizes(u3din&lat),dimsizes(u3din&lon)/),typeof(u3din),u850@_FillValue)
        olrr_out = new((/1,dimsizes(u3din&lat),dimsizes(u3din&lon)/),typeof(u3din),u850@_FillValue)

        u850_out(0,:,:) = u850(0,:,:)
        prec_out(0,:,:) = prect(:,:)
        olrr_out(0,:,:) = olr(:,:)

        ctime = cd_inv_calendar(years(iyear),mons(imon),iday+1,0,0,0,"hours since 1900-01-01 00:00:00",0)

        u850_out!0   = "time"
        u850_out!1   = "lat"
        u850_out!2   = "lon"
        u850_out&time= ctime
        u850_out&lat = u3din&lat
        u850_out&lon = u3din&lon
        copy_VarCoords(u850_out,prec_out)
        copy_VarCoords(u850_out,olrr_out)

        putnc2d(u850_out,"U850",\
                olrr_out,"OLR",\
                prec_out,"prect",\
                "GRIST.ATM.G6.amp."+years(iyear)+"-"+cmon+"-"+cday+".grid.nc")

        delete([/u850,u850_out,prec_out,olrr_out/])

end do
end do
end do

end
