from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors

def yieldout(year):
        bb=year-1900
	region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
	maitrop = region1.variables['soy_trop'][bb-1,:,:]
	maitemp = region1.variables['soy_temp'][bb-1,:,:]
	maitropi=region1.variables['soy_trop_irrig'][bb-1,:,:]
	maitempi=region1.variables['soy_temp_irrig'][bb-1,:,:]
	gridarea = region1.variables['area'][:,:]
	maitrop=ma.masked_where(maitrop<=0,maitrop)
	maitrop=ma.filled(maitrop, fill_value=0.)
	maitemp=ma.masked_where(maitemp<=0,maitemp)
	maitemp=ma.filled(maitemp, fill_value=0.)

	maitropi=ma.masked_where(maitropi<=0,maitropi)
	maitropi=ma.filled(maitropi, fill_value=0.)
	maitempi=ma.masked_where(maitempi<=0,maitempi)
	maitempi=ma.filled(maitempi, fill_value=0.)

	maizetor=maitrop+maitemp
	maizetoi=maitropi+maitempi

	maizetrop=maitrop+maitropi
	maizetemp=maitemp+maitempi

	maizeto = maitrop+maitemp+maitropi+maitempi



	ff=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalFertilizer.nc','r')
	fert_maitrop = ff.variables['soy_trop_fert'][bb-1,:,:]
	fert_maitemp = ff.variables['soy_temp_fert'][bb-1,:,:]




	clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_rf_fert_0.5x0.5.nc','r')
#	clmtropf = clm.variables['yield'][bb,:,:]
	clmtropfer=clm.variables['fertilizer'][bb-1,:,:]

	clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
#	clmtempf = clm1.variables['yield'][bb,:,:]
	clmtempfer=clm1.variables['fertilizer'][bb-1,:,:]

        isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytrop_historical_co2_rf_fert_0.5x0.5.nc','r')
        clmtropf = isam.variables['totalyield'][bb-1,:,:]

        isam1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
        clmtempf = isam1.variables['totalyield'][bb-1,:,:]




	clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
	clmtropfi = clm2.variables['totalyield'][bb-1,:,:]

	clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
	clmtempfi = clm3.variables['totalyield'][bb-1,:,:]



	clma=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
	clmtropfno = clma.variables['totalyield'][bb-1,:,:]

	clm1a=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
	clmtempfno = clm1a.variables['totalyield'][bb-1,:,:]

	clm2a=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytrop_historical_co2_irrig_nofert_0.5x0.5.nc','r')
	clmtropfnoi = clm2a.variables['totalyield'][bb-1,:,:]

	clm3a=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
	clmtempfnoi = clm3a.variables['totalyield'][bb-1,:,:]
        lonisam=clm3a.variables['lon'][:]

        clmtempfnoi,lonisam1 = shiftgrid(180.5,clmtempfnoi,lonisam,start=False)
        clmtropfnoi,lonisam1 = shiftgrid(180.5,clmtropfnoi,lonisam,start=False)
        clmtempfno,lonisam1 = shiftgrid(180.5,clmtempfno,lonisam,start=False)
        clmtropfno,lonisam1 = shiftgrid(180.5,clmtropfno,lonisam,start=False)

        clmtempf,lonisam1 = shiftgrid(180.5,clmtempf,lonisam,start=False)
        clmtropf,lonisam1 = shiftgrid(180.5,clmtropf,lonisam,start=False)
        clmtempfi,lonisam1 = shiftgrid(180.5,clmtempfi,lonisam,start=False)
        clmtropfi,lonisam1 = shiftgrid(180.5,clmtropfi,lonisam,start=False)
        #print lonisam1

	clmtropfer=N.flipud(clmtropfer)
	clmtempfer=N.flipud(clmtempfer)


	clmtropf= ma.masked_where(maitrop<=0,clmtropf)
	clmtempf= ma.masked_where(maitemp<=0,clmtempf)
	clmtropf=ma.filled(clmtropf, fill_value=0.)
	clmtempf=ma.filled(clmtempf, fill_value=0.)

	clmtropfi= ma.masked_where(maitropi<=0,clmtropfi)
	clmtempfi= ma.masked_where(maitempi<=0,clmtempfi)
	clmtropfi=ma.filled(clmtropfi, fill_value=0.)
	clmtempfi=ma.filled(clmtempfi, fill_value=0.)


	clmtropfno= ma.masked_where(maitrop<=0,clmtropfno)
	clmtempfno= ma.masked_where(maitemp<=0,clmtempfno)
	clmtropfno=ma.filled(clmtropfno, fill_value=0.)
	clmtempfno=ma.filled(clmtempfno, fill_value=0.)

	clmtropfnoi= ma.masked_where(maitropi<=0,clmtropfnoi)
	clmtempfnoi= ma.masked_where(maitempi<=0,clmtempfnoi)
	clmtropfnoi=ma.filled(clmtropfnoi, fill_value=0.)
	clmtempfnoi=ma.filled(clmtempfnoi, fill_value=0.)


	fertfractiontrop= N.zeros((360, 720))
	nofertfractiontrop= N.zeros((360, 720))
	fertfractiontemp= N.zeros((360, 720))
	nofertfractiontemp= N.zeros((360, 720))



	for x in range(0,360):
	        for y in range(0,720):
	                if clmtropfer[x,y] > 0.0:
	                        fertfractiontrop[x,y] = min(1.0,fert_maitrop[x,y]/clmtropfer[x,y])
	                        nofertfractiontrop[x,y] = 1.0 - fertfractiontrop[x,y]
	                else:
	                        fertfractiontrop[x,y]= 0.0
	                        nofertfractiontrop[x,y] = 1.0

	for x in range(0,360):
	        for y in range(0,720):
	                if clmtempfer[x,y] > 0.0:
	                        fertfractiontemp[x,y] = min(1.0,fert_maitemp[x,y]/clmtempfer[x,y])
	                        nofertfractiontemp[x,y] = 1.0 - fertfractiontemp[x,y]
	                else:
	                        fertfractiontemp[x,y]= 0.0
	                        nofertfractiontemp[x,y] = 1.0

	clmtropfnew= N.zeros((360, 720))
	clmtempfnew= N.zeros((360, 720))
	clmtropfinew= N.zeros((360, 720))
	clmtempfinew= N.zeros((360, 720))

	for x in range(0,360):
	        for y in range(0,720):
			clmtropfnew[x,y] = (nofertfractiontrop[x,y]*clmtropfno[x,y])+(fertfractiontrop[x,y]*clmtropf[x,y])
	                clmtempfnew[x,y] = (nofertfractiontemp[x,y]*clmtempfno[x,y])+(fertfractiontemp[x,y]*clmtempf[x,y])
	                clmtropfinew[x,y] = (nofertfractiontrop[x,y]*clmtropfnoi[x,y])+(fertfractiontrop[x,y]*clmtropfi[x,y])
	                clmtempfinew[x,y] = (nofertfractiontemp[x,y]*clmtempfnoi[x,y])+(fertfractiontemp[x,y]*clmtempfi[x,y])



	yield_clmtf=clmtropf+clmtempf
	yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
	#yield_clmtf  = ma.masked_where(maizetor<=0,yield_clmtf )
	yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

	yield_clmtfi=clmtropfi+clmtempfi
	yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
	#yield_clmtfi = ma.masked_where(maizetoi<=0,yield_clmtfi)
	yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)


	yield_clmtfnew=clmtropfnew+clmtempfnew
	yield_clmtfnew = ma.masked_where(yield_clmtfnew<=0,yield_clmtfnew)
	#yield_clmtf  = ma.masked_where(maizetor<=0,yield_clmtf )
	yield_clmtfnew=ma.filled(yield_clmtfnew, fill_value=0.)

	yield_clmtfinew=clmtropfinew+clmtempfinew
	yield_clmtfinew = ma.masked_where(yield_clmtfinew<=0,yield_clmtfinew)
	#yield_clmtfi = ma.masked_where(maizetoi<=0,yield_clmtfi)
	yield_clmtfinew=ma.filled(yield_clmtfinew, fill_value=0.)


	area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
	gridarea = area.variables['cell_area'][:,:]
	gridlon = area.variables['lon'][:]
	gridlat=area.variables['lat'][:]
	gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)


	lon2,lat2 = N.meshgrid(gridlon,gridlat)


	map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
	x,y = map(lon2,lat2)

	yield_clmtf=maskoceans(x,y,yield_clmtf)
	yield_clmtf = ma.masked_where(maizeto<=0,yield_clmtf)

	yield_clmtfi=maskoceans(x,y,yield_clmtfi)
	yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)

	clmy=((yield_clmtf*maizetor*gridarea)+(yield_clmtfi*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))




	yield_clmtfnew=maskoceans(x,y,yield_clmtfnew)
	yield_clmtfnew = ma.masked_where(maizeto<=0,yield_clmtfnew)

	yield_clmtfinew=maskoceans(x,y,yield_clmtfinew)
	yield_clmtfinew = ma.masked_where(maizeto<=0,yield_clmtfinew)

	clmynew=((yield_clmtfnew*maizetor*gridarea)+(yield_clmtfinew*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))
       
        return clmy

yief= N.zeros((105, 360, 720))
years = range(1901, 2006)
for i, yeara in enumerate(years):


	yie=yieldout(yeara)
        yief[i, :, :] = yie

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon


ncfile=NetCDFFile('isamhiscru_soyscaleiyield_fertfao_new1.nc','w',format='NETCDF3_64BIT_OFFSET')
ncfile.createDimension('lat', 360)
ncfile.createDimension('lon', 720)
ncfile.createDimension('time', 105)

times = ncfile.createVariable('time', 'f8', ('time',))
latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
maize= ncfile.createVariable('yield', 'f8', ('time','lat','lon'),fill_value=-9999.)
latitudes[:] = gridlat
longitudes[:] = gridlon
times[:] = years
maize[:]=yief



latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'

latitudes.long_name = 'latitude'
longitudes.long_name = 'longitude'


ncfile.close()





