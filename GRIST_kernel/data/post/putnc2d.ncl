undef("putnc2d")
procedure putnc2d ( var1:numeric, var1name:string, \
                    var2:numeric, var2name:string, \
                    var3:numeric, var3name:string, filename :string )
;local ntim, nlat, nlon, dimNames, dimSizes, dimUnlim,fout
begin
;===================================================================  
; Assume variables T, PS and ORO exist and that they have 
; associated meta data: (a) coordinate variables time, lev, lat, lon       
; and (b) attributes
;===================================================================  
      ntim  = dimsizes(var1&time)             ; get dimension sizes  
      nlat  = dimsizes(var1&lat)
      nlon  = dimsizes(var1&lon)

      diro = ""                            ; Output directory
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
      fAtt@title         = filename
      fAtt@source_file   = "None"
      fAtt@Conventions   = "None"
;      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )             ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================

      dimNames = (/"time",  "lat", "lon"/)
      dimSizes = (/ ntim ,   nlat,  nlon/)
      dimUnlim = (/ True ,   False, False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)

;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining 
; each variable's dimension name(s) and type. 
; The following could be replaced with explicit, user defined dimension 
; names different from those associated with the variable in memory. 
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","lat","lon"/)) 
;===================================================================
      filevardef(fout, "time",    typeof(var1&time)    ,getvardims(var1&time))
      filevardef(fout, "lat",     typeof(var1&lat)     ,getvardims(var1&lat))
      filevardef(fout, "lon",     typeof(var1&lon)     ,getvardims(var1&lon))
      filevardef(fout, var1name,  typeof(var1)         ,getvardims(var1))
      filevardef(fout, var2name,  typeof(var2)         ,getvardims(var2))
      filevardef(fout, var3name,  typeof(var3)         ,getvardims(var3))
                                                              ; different from name on script
;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================
      filevarattdef(fout,var1name   ,var1)                          ; copy T attributes
      filevarattdef(fout,"time"     ,var1&time)                     ; copy time attributes
      filevarattdef(fout,"lat"      ,var1&lat)                      ; copy lat attributes
      filevarattdef(fout,"lon"      ,var1&lon)                      ; copy lon attributes

;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
      setfileoption(fout,"DefineMode",False)

;===================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;====================================================================
      fout->time       = (/var1&time/)
      fout->lat        = (/var1&lat/)
      fout->lon        = (/var1&lon/)
      fout->$var1name$ = (/var1/)
      fout->$var2name$ = (/var2/)
      fout->$var3name$ = (/var3/)

end
