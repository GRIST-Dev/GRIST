import cdsapi

c = cdsapi.Client()

for year in range(2012, 2012 + 1):
   for imon in range(5, 5 +1):
      if (imon==5):
         en=26
      for iday in range (26, 31+1):
          mon  = "%02d" % imon
          day  = "%02d" % iday
          date = str(year)+'-'+str(mon)+'-'+str(day)+'/to/'+str(year)+'-'+str(mon)+'-'+str(day)
          print(date)

          c.retrieve(
               'reanalysis-era5-pressure-levels',
           {
               'product_type':'reanalysis',
               'format':'grib',
               'variable':[
                   'specific_humidity','temperature','u_component_of_wind',
                   'v_component_of_wind'
               ],
               'pressure_level':[
                   '1','2','3','5','7','10',
                   '20', '30','50','70','100','125',
                   '150','175','200','225','250','300',
                   '350','400','450','500','550','600',
                   '650','700','750','775','800','825',
                   '850','875','900','925','950','975',
                   '1000'
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
           'ERA5.pl.'+str(year)+str(mon)+str(day)+'.00.grib')
