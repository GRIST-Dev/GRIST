undef("putnc2d")
procedure putnc2d ( var:numeric, varname :string, filename:string )
;local ntim, nlatitude, nlongitude, dimNames, dimSizes, dimUnlim,fout
begin
;===================================================================  
; Assume variables T, PS and ORO exist and that they have 
; associated meta data: (a) coordinate variables time, lev, latitude, longitude       
; and (b) attributes
;===================================================================  
      ntim       = dimsizes(var&time)                 ; get dimension sizes  
      nlatitude  = dimsizes(var&latitude)
      nlongitude = dimsizes(var&longitude)

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
      fAtt@title         = filename
      fAtt@source_file   = "None"
      fAtt@Conventions   = "None"
      fAtt@creation_date = systemfunc ("date")
      fileattdef( fout, fAtt )             ; copy file attributes    

;===================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;===================================================================
      dimNames = (/"time", "latitude", "longitude"/)
      dimSizes = (/ ntim, nlatitude,  nlongitude/)
      dimUnlim = (/ True, False, False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)
;===================================================================
; predefine the the dimensionality of the variables to be written out
;===================================================================
; Here we are using NCL functions to facilitate defining 
; each variable's dimension name(s) and type. 
; The following could be replaced with explicit, user defined dimension 
; names different from those associated with the variable in memory. 
; Say, PS(time,latitude,longitude) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/)) 
;===================================================================
      filevardef(fout, "time" ,typeof(var&time) ,getvardims(var&time))
      filevardef(fout, "latitude"  ,typeof(var&latitude)  ,getvardims(var&latitude))
      filevardef(fout, "longitude"  ,typeof(var&longitude)  ,getvardims(var&longitude))
      filevardef(fout, varname,typeof(var)      ,getvardims(var))
;===================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;====================================================================
      filevarattdef(fout,"time" ,var&time)                      ; copy time attributes
      filevarattdef(fout,"latitude"  ,var&latitude)                       ; copy latitude  attributes
      filevarattdef(fout,"longitude"  ,var&longitude)                       ; copy longitude  attributes
      filevarattdef(fout,varname,var)                           ; copy var  attributes
;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;===================================================================
      setfileoption(fout,"DefineMode",False)
;===================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;====================================================================
      fout->time        = (/var&time/)
      fout->latitude    = (/var&latitude/)
      fout->longitude   = (/var&longitude/)
      fout->$varname$   = (/var/)

end
