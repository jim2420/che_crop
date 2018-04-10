from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors


datc=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
iyield1c = datc.variables['yield'][99,:,:]
dat2c=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2c = dat2c.variables['yield'][99,:,:]
iyield1c=N.flipud(iyield1c)
iyield2c=N.flipud(iyield2c)



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/noheat/soyhis/output/soyhis.nc','r')
iyield1 = dat.variables['totalyield'][99,:,:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/noheat/soyhis_fert/output/soyhis_fert.nc','r')
iyield2 = dat2.variables['totalyield'][99,:,:]
dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/noheat/soyhis_fertall/output/soyhis_fertall.nc','r')
iyield3 = dat3.variables['totalyield'][99,:,:]
gridlon=dat3.variables['lon'][:]
gridlat1=dat3.variables['lat'][:]
iyield1f,gridlon1 = shiftgrid(180.5,iyield1,gridlon,start=False)
iyield2f,gridlon1 = shiftgrid(180.5,iyield2,gridlon,start=False)
iyield3f,gridlon1 = shiftgrid(180.5,iyield3,gridlon,start=False)



fig = plt.figure(figsize=(20,15))
ax2 = fig.add_subplot(421)
ax2.set_title("CLM Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
#iyield1c=maskoceans(x,y,iyield1c)
iyield1c = ma.masked_where(iyield1c<=0,iyield1c)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield1c,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(422)
ax2.set_title("CLM Soy Yield Fert (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#iyield2c=maskoceans(x,y,iyield2c)
iyield2c = ma.masked_where(iyield2c<=0,iyield2c)

cs = map.pcolormesh(x,y,iyield2c,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(423)
ax2.set_title("ISAM Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield1f,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(424)
ax2.set_title("ISAM Soy Yield Fert (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield2f,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(425)
ax2.set_title("ISAM Soy Yield No N limit (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield3f,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(426)
ax2.set_title("ISAM Fert diff (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield2f-iyield1f,cmap=plt.cm.bwr,vmin=-2,vmax=2)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(427)
ax2.set_title("ISAM No limitation diff Fert (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(gridlon1,gridlat1)
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield3f-iyield2f,cmap=plt.cm.bwr,vmin=-2,vmax=2)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')




plt.savefig('soyfert_2000.jpg',dpi=300,bbox_inches='tight')


plt.show()


