from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors





area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['soybeanData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
latnc=N.flipud(latnc)


ncvar_maizef= N.zeros((2160, 4320))
ncvar_maizef=ncvar_maize[0,1,:,:]
ncvar_maize1=ncvar_maize[0,4,:,:]
ncvar_mask= N.zeros((2160, 4320))
ncvar_mask=ncvar_maize[0,0,:,:]


ncvar_maizef[N.isnan(ncvar_maizef)] = -9999
ncvar_mask[N.isnan(ncvar_mask)] = -9999
ncvar_maize1[N.isnan(ncvar_maize1)] = -9999

ncvar_maizef = ma.masked_where(ncvar_maizef<=0,ncvar_maizef)
#ncvar_maizef= ma.masked_where(ncvar_mask<0.01,ncvar_maizef)
ncvar_maizef = ma.masked_where(ncvar_maize1<=0,ncvar_maizef)

ncvar_maize1=N.flipud(ncvar_maize1)
ncvar_maizef=N.flipud(ncvar_maizef)
ncvar_mask=N.flipud(ncvar_mask)


lon2,lat2 = N.meshgrid(gridlon,gridlat)

iyield = interp(ncvar_maizef,lonnc,latnc,lon2,lat2  , order=1)
iarea=interp(ncvar_mask,lonnc,latnc,lon2,lat2  , order=1)

dat=NetCDFFile('finalyield/isam/heat/fertfao/isamhis_soyscaleiyield_heat_fertfao.nc','r')
iyield1 = dat.variables['yield'][99,:,:]
dat2=NetCDFFile('finalyield/isam/heat/isamhis_soyscaleiyield_heat.nc','r')
iyield2 = dat2.variables['yield'][99,:,:]
dat3=NetCDFFile('finalyield/isam/heat/isamhis_soyscaleifyield_heat.nc','r')
iyield3 = dat3.variables['yield'][99,:,:]
dat4=NetCDFFile('finalyield/isam/heat/fertall/isamhis_soyscaleiyield_heat_fertall.nc','r')
iyield4 = dat4.variables['yield'][99,:,:]
dat5=NetCDFFile('finalyield/isam/noheat_notrop/isamhiscru_soyscaleiyield.nc','r')
iyield5 = dat5.variables['yield'][99,:,:]
dat6=NetCDFFile('finalyield/isam/heat/isamhiscru_soyscaleiyield_heat.nc','r')
iyield6 = dat6.variables['yield'][99,:,:]
dat7=NetCDFFile('finalyield/isam/noheat_notrop/isamhiscru_soyscaleifyield.nc','r')
iyield7 = dat7.variables['yield'][99,:,:]
dat8=NetCDFFile('finalyield/isam/heat/isamhiscru_soyscaleifyield_heat.nc','r')
iyield8 = dat8.variables['yield'][99,:,:]






fig = plt.figure(figsize=(20,15))
ax2 = fig.add_subplot(421)
ax2.set_title("CESM-ISAM Soy Yield FAON(t/ha)",fontsize=20)

##us

#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
#map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')
##us


map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield1 = ma.masked_where(iyield1<=0,iyield1)
iyield1 = ma.masked_where(iarea<=0,iyield1)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield1,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(422)
ax2.set_title("CESM-ISAM Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield2 = ma.masked_where(iyield2<=0,iyield2)
iyield2 = ma.masked_where(iarea<=0,iyield2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield2,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(423)
ax2.set_title("CESM-ISAM Nscale Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield3 = ma.masked_where(iyield3<=0,iyield3)
iyield3 = ma.masked_where(iarea<=0,iyield3)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield3,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(424)
ax2.set_title("CESM-ISAM Soy Yield Fertall (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield4 = ma.masked_where(iyield4<=0,iyield4)
iyield4 = ma.masked_where(iarea<=0,iyield4)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield4,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(425)
ax2.set_title("NCEP-ISAM Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield5 = ma.masked_where(iyield5<=0,iyield5)
iyield5 = ma.masked_where(iarea<=0,iyield5)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield5,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(426)
ax2.set_title("NCEP-ISAM HS Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield6 = ma.masked_where(iyield6<=0,iyield6)
iyield6 = ma.masked_where(iarea<=0,iyield6)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield6,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(427)
ax2.set_title("NCEP-ISAM Nscale Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield7 = ma.masked_where(iyield7<=0,iyield7)
iyield7 = ma.masked_where(iarea<=0,iyield7)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield7,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(428)
ax2.set_title("NCEP-ISAM Nscale HS Soy Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield8 = ma.masked_where(iyield8<=0,iyield8)
iyield8 = ma.masked_where(iarea<=0,iyield8)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield8,cmap=plt.cm.YlGn,vmin=0,vmax=6)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('soyisam_2000.jpg',dpi=300,bbox_inches='tight')


plt.show()
