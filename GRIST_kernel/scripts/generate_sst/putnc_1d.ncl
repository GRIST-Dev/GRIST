undef("putnc1d")
procedure putnc1d ( var, var1, lat,lon, varname, var1name, filename)
;local ntim, nncol, nlongitude, dimNames, dimSizes, dimUnlim,fout
begin
;===================================================================  
; Assume variables T, PS and ORO exist and that they have 
; associated meta data: (a) coordinate variables time, lev, ncol, longitude       
; and (b) attributes
;===================================================================  
      ntim       = dimsizes(var&time)                 ; get dimension sizes  
      nncol  = dimsizes(var&ncol)

      diro = "./"                            ; Output directory
      filo = filename              ; Output file
      system("/bin/rm -f " + diro + filo)    ; remove if exists
      fout  = addfile (diro + filo, "c")     ; open output file

;===================================================================
; explicitly declare file definition mode. Improve efficiency.
;===================================================================
      setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file
;===================================================================
      fAtt               = True            ; assign file attributes
      fAtt@title         = varname
      fAtt@source_file   = "None"
      fAtt@Conventions   = "None"
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )             ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================

      dimNames = (/"time", "ncol"/)
      dimSizes = (/ ntim ,  nncol/)
      dimUnlim = (/ True , False     /)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)

;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining 
; each variable's dimension name(s) and type. 
; The following could be replaced with explicit, user defined dimension 
; names different from those associated with the variable in memory. 
; Say, PS(time,ncol,longitude) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","ncol","longitude"/)) 
;===================================================================
      filevardef(fout, "time"     ,typeof(var&time)    ,getvardims(var&time))
      filevardef(fout, "lat"      ,typeof(lat)         ,getvardims(lat))
      filevardef(fout, "lon"      ,typeof(lon)         ,getvardims(lon))
      filevardef(fout,  varname   ,typeof(var)         ,getvardims(var))
      filevardef(fout,  var1name  ,typeof(var1)        ,getvardims(var1))
;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================
      filevarattdef(fout, varname    ,var)                       ; copy T attributes
      filevarattdef(fout, var1name   ,var1)                      ; copy T attributes
      filevarattdef(fout,"time"      ,var&time)                  ; copy time attributes
      filevarattdef(fout,"lat"       ,lat)                       ; copy time attributes
      filevarattdef(fout,"lon"       ,lon)                       ; copy time attributes
;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
      setfileoption(fout,"DefineMode",False)
;===================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined ncols on the file.
;====================================================================
      fout->time         = (/var&time/)
      fout->lat          = (/lat/)
      fout->lon          = (/lon/)
      fout->$varname$    = (/var/)
      fout->$var1name$   = (/var1/)

end
