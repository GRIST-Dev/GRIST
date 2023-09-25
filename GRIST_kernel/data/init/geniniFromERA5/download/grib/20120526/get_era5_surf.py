import cdsapi

c = cdsapi.Client()

for year in range(2012, 2012 + 1):
   for imon in range(5, 5 +1):
      if (imon==5):
         en=26
      for iday in range (26, 26+1):
          mon  = "%02d" % imon
          day  = "%02d" % iday
          date = str(year)+'-'+str(mon)+'-'+str(day)+'/to/'+str(year)+'-'+str(mon)+'-'+str(day)
          print(date)

          c.retrieve(
               'reanalysis-era5-single-levels',
           {
               'product_type':'reanalysis',
               'format':'grib',
               'variable':[
                   'surface_pressure', 'Geopotential', 'sea_ice_cover',
                   'skin_temperature', 'snow_density', 'snow_depth',
                   'soil_temperature_level_1', 'soil_temperature_level_2', 'soil_temperature_level_3', 'soil_temperature_level_4',
                   'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4'
               ],
               'year':str(year),
               'month':str(mon),
               'day':[
                   str(day),
               ],
               'time':[
                   '00:00',
               ]
           },
           'ERA5.sf.'+str(year)+str(mon)+str(day)+'.00.grib')
