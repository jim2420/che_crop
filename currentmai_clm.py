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
lonab = nclu.variables["lon"][:]

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][4,:,:]
maitemp = region1.variables['maize_temp'][4,:,:]
maitropi = region1.variables['maize_trop_irrig'][4,:,:]
maitempi = region1.variables['maize_temp_irrig'][4,:,:]
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


isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/clm/clm45his_maiscaleifyield.nc','r')
isamyield = N.average(isam.variables['yield'][95:105,:,:],axis=0)
lona1 = isam.variables["lon"][:]
lata1 = isam.variables["lat"][:]
lon2,lat2 = N.meshgrid(lona1,lata1)

isamyield[N.isnan(isamyield)] = -9999
isamyield = ma.masked_where(isamyield<=0,isamyield)
isamyield = ma.masked_where(maizetotal<=0,isamyield)
#isamyield=ma.filled(isamyield, fill_value=0.)

yieldfa= N.zeros((10, 360, 720))
yieldf2a= N.zeros((10, 360, 720))
yieldf= N.zeros((10, 360, 720))
yieldf2= N.zeros((10, 360, 720))
years2 = range(2090,2100)

clmtropfin= N.zeros((10, 360, 720))
clmtempfin= N.zeros((10, 360, 720))

for j in range(0,10):

        clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetrop_rcp45_co2_rf_fert_0.5x0.5.nc','r')
        cc = N.flipud(clm2n.variables['yield'][84+j,:,:])
        clmtropfin[j,:,:] = cc

        clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetemp_rcp45_co2_rf_fert_0.5x0.5.nc','r')
        dd = N.flipud(clm3n.variables['yield'][84+j,:,:])
        clmtempfin[j,:,:] = dd

yield_new2a=N.average(clmtropfin,axis=0)
yield_new2=N.average(clmtempfin,axis=0)


yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)
yield_new2a= ma.masked_where(yield_new2a<=0.,yield_new2a)
yield_new2=ma.masked_where( maitoatemp<=0,yield_new2)
yield_new2=ma.filled(yield_new2, fill_value=0.)

yield_new2a=ma.masked_where( maitoatrop<=0,yield_new2a)
yield_new2a=ma.filled(yield_new2a, fill_value=0.)

yield_new2=yield_new2+yield_new2a

yield_new2=ma.masked_where(isamyield<=0.,yield_new2)
yieldisam=yield_new2-isamyield
yieldisam= ma.masked_where(yieldisam==0.,yieldisam)


ncvar_maize,lona11 = shiftgrid(180.5,ncvar_maize,lonab,start=False)
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

plt.savefig('mai2006_2015_rcp45_clmc.jpg',dpi=300,bbox_inches='tight')
plt.show()





