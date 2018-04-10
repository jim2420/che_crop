from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors





area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['maizeData'][:]
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
ncvar_maip=ncvar_maize[0,5,:,:]


ncvar_maizef[N.isnan(ncvar_maizef)] = -9999
ncvar_mask[N.isnan(ncvar_mask)] = -9999
ncvar_maize1[N.isnan(ncvar_maize1)] = -9999
ncvar_maip[N.isnan(ncvar_maip)] = -9999

ncvar_maizef = ma.masked_where(ncvar_maizef<=0,ncvar_maizef)
#ncvar_maizef= ma.masked_where(ncvar_mask<0.01,ncvar_maizef)
ncvar_maizef = ma.masked_where(ncvar_maize1<=0,ncvar_maizef)
ncvar_maip = ma.masked_where(ncvar_maize1<=0,ncvar_maip)

ncvar_maize1=N.flipud(ncvar_maize1)
ncvar_maizef=N.flipud(ncvar_maizef)
ncvar_mask=N.flipud(ncvar_mask)
ncvar_maip=N.flipud(ncvar_maip)


lon2,lat2 = N.meshgrid(gridlon,gridlat)

iyield = interp(ncvar_maizef,lonnc,latnc,lon2,lat2  , order=1)
iarea=interp(ncvar_mask,lonnc,latnc,lon2,lat2  , order=1)
iyield1 = interp(ncvar_maip,lonnc,latnc,lon2,lat2  , order=1)



fig = plt.figure(figsize=(20,15))

ax = fig.add_subplot(211)
ax.set_title("M3 Maize Yield (t/ha)",fontsize=20)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iyield = ma.masked_where(iarea<=0,iyield)
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
iizumy=iyield
cs = map.pcolormesh(x,y,iizumy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax = fig.add_subplot(212)
ax.set_title("M3 Maize Production (t)",fontsize=20)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iyield = ma.masked_where(iarea<=0,iyield)
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
iizumy1=iyield1
cs = map.pcolormesh(x,y,iizumy1,cmap=plt.cm.YlGn,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=12000)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('m3mai.jpg',dpi=600,bbox_inches='tight')


plt.show()
