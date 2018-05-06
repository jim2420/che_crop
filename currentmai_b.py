from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import math
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime
from statsmodels.stats.weightstats import DescrStatsW
import matplotlib.colors as colors

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
ncvar_maize = nclu.variables['maizey'][0,:,:]

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][84:94,:,:]
maitemp = region1.variables['maize_temp'][84:94,:,:]
maitropi = region1.variables['maize_trop_irrig'][84:94,:,:]
maitempi = region1.variables['maize_temp_irrig'][84:94,:,:]
maitrop= ma.masked_where(maitrop<=0,maitrop)
maitropi= ma.masked_where(maitropi<=0,maitropi)
maitemp= ma.masked_where(maitemp<=0,maitemp)
maitempi= ma.masked_where(maitempi<=0,maitempi)
maitrop=ma.filled(maitrop, fill_value=0.)
maitropi=ma.filled(maitropi,fill_value=0.)
maitemp=ma.filled(maitemp, fill_value=0.)
maitempi=ma.filled(maitempi, fill_value=0.)
maizeto = maitrop+maitemp
maizetoi = maitropi+maitempi
maitoatemp=maitemp+maitempi
maitoatrop=maitrop+maitropi
maizetotal = maizeto+maizetoi


region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrops = region1.variables['maize_trop'][95:105,:,:]
maitemps = region1.variables['maize_temp'][95:105,:,:]
maitropis = region1.variables['maize_trop_irrig'][95:105,:,:]
maitempis = region1.variables['maize_temp_irrig'][95:105,:,:]
maitrops= ma.masked_where(maitrops<=0,maitrops)
maitropis= ma.masked_where(maitropis<=0,maitropis)
maitemps= ma.masked_where(maitemps<=0,maitemps)
maitempis= ma.masked_where(maitempis<=0,maitempis)
maitrops=ma.filled(maitrops, fill_value=0.)
maitropis=ma.filled(maitropis,fill_value=0.)
maitemps=ma.filled(maitemps, fill_value=0.)
maitempis=ma.filled(maitempis, fill_value=0.)
maizetos = maitrops+maitemps
maizetois = maitropis+maitempis
maitoatemps=maitemps+maitempis
maitoatrops=maitrops+maitropis

maizetotals = maizetos+maizetois

isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_maiscaleifyield_heat.nc','r')
#isamyield = N.average(isam.variables['yield'][95:105,:,:],axis=0)
isamyielda = isam.variables['yield'][95:105,:,:]

lona1 = isam.variables["lon"][:]
lata1 = isam.variables["lat"][:]
lon2,lat2 = N.meshgrid(lona1,lata1)

#isamyield[N.isnan(isamyield)] = -9999
#isamyield = ma.masked_where(isamyield<=0,isamyield)
#isamyield = ma.masked_where(maizetotal<=0,isamyield)

isamyielda[N.isnan(isamyielda)] = -9999
isamyielda = ma.masked_where(isamyielda<=0,isamyielda)
isamyielda = ma.masked_where(maizetotals<=0,isamyielda)

isamyield = N.average(isamyielda,axis=0)



yieldfa= N.zeros((10, 360, 720))
yieldf2a= N.zeros((10, 360, 720))
yieldf= N.zeros((10, 360, 720))
yieldf2= N.zeros((10, 360, 720))
years2 = range(2090,2100)

for i, year1 in enumerate(years2):


    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp45/heat/maihis_irr_fert/output/hmaizehis_irr_fert.bgp-yearly_crop_{0}.nc".format(year1), mode='r')  
    yield2 = base2.variables["yield"][0,:,:]
    yieldf2[i, :, :] = yield2
    lona1 = base2.variables["lon"][:]
    lata1 = base2.variables["lat"][:]

    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp45/heat/maihis_irr_fertrop/output/hmaizehis_irr_fertrop.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    yield2a = base2a.variables["yield"][0,:,:]
    yieldf2a[i, :, :] = yield2a

yielda2=N.average(yieldf2,axis=0)
yield_new2,lona11 = shiftgrid(180.5,yieldf2,lona1,start=False)
yielda2a=N.average(yieldf2a,axis=0)
yield_new2a,lona11 = shiftgrid(180.5,yieldf2a,lona1,start=False)


yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)
yield_new2a= ma.masked_where(yield_new2a<=0.,yield_new2a)
yield_new2=ma.masked_where( maitoatemp<=0,yield_new2)
yield_new2=ma.filled(yield_new2, fill_value=0.)
yield_new2a=ma.masked_where( maitoatrop<=0,yield_new2a)
yield_new2a=ma.filled(yield_new2a, fill_value=0.)
yield_new21=yield_new2+yield_new2a
yield_new2= N.average(yield_new21,axis=0)

yield_new2=ma.masked_where(isamyield<=0.,yield_new2)
isamyield=ma.masked_where(yield_new2<=0.,isamyield)
#yield_new2=ma.masked_where(isamyield<=0.,yield_new2)
#isamyield=ma.masked_where(yield_new2<=0.,isamyield)

yieldisam=yield_new2-isamyield
yieldisam= ma.masked_where(yieldisam==0.,yieldisam)


ncvar_maize,lona11 = shiftgrid(180.5,ncvar_maize,lona1,start=False)
yieldisam= ma.masked_where(ncvar_maize<=0.,yieldisam)








fig = plt.figure(figsize=(12,6))

ax1 = fig.add_subplot(211)
ax1.set_title("isam mai yield 1996-2005 (t/ha)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

x,y = map(lon2,lat2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,isamyield,cmap=plt.cm.jet,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=16)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax2 = fig.add_subplot(212)
ax2.set_title("Maize RCP4.5",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines(color='gray')
#map.drawstates()
#map.drawcountries(color='b')

#map.drawcoastlines()
#map.drawcountries()
#map.drawmapboundary(color='gray')
cs = map.pcolormesh(x,y,yieldisam/isamyield*100,cmap=plt.cm.bwr,vmin=-100.,vmax=100.)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')
plt.tight_layout()

plt.savefig('mai2006_2015_rcp45_b.jpg',dpi=300,bbox_inches='tight')
plt.show()





