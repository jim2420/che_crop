from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
from statsmodels.stats.weightstats import DescrStatsW
import matplotlib.colors as colors

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/spammaize_isam.nc','r')
nclu1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/spamsoy_isam.nc','r')

m3newp = nclu.variables['maip_total'][:,:]
m3newps = nclu1.variables['soyp_total'][:,:]
m3newa = nclu.variables['maia_total'][:,:]
m3newas = nclu1.variables['soya_total'][:,:]

latisam1=nclu.variables['lat'][:]
lonisam1=nclu.variables['lon'][:]




region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_maiscaleiyield_fertfao_new1.nc','r')
ma1 = region.variables['yield'][103:105,:,:]

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_soyscaleiyield_fertfao_new1.nc','r')
ma2 = region1.variables['yield'][103:105,:,:]

ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)

latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]



lon,lat = N.meshgrid(lonmask,latmask)

m3newp= ma.masked_where(ma1<=0.0,m3newp)
ma1= ma.masked_where(m3newp<=0.0,ma1)
m3newp= ma.masked_where(ma1<=0.0,m3newp)
ma1= ma.masked_where(m3newp<=0.0,ma1)

m3newps= ma.masked_where(ma2<=0.0,m3newps)
ma2= ma.masked_where(m3newps<=0.0,ma2)
m3newps= ma.masked_where(ma2<=0.0,m3newps)
ma2= ma.masked_where(m3newps<=0.0,ma2)
#cmap = plt.cm.jet
#bounds=[0,0.01,0.05,0.1,0.5,1,2,3,4,5,6]
cmap = plt.cm.terrain_r
bounds=[-0.1,0.0,0.01,0.05,0.1,0.5,1,2,3,4,5,6]

norm = colors.BoundaryNorm(bounds, cmap.N)

fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(221)
ax1.set_title("Maize SPAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
x,y = map(lon,lat)
m3newp=maskoceans(x,y,m3newp)
m3newa=maskoceans(x,y,m3newa)
cs1 = map.pcolormesh(x,y,m3newp*10000/gridarea,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,m3newp*10000/gridarea,cmap=plt.cm.Greens,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=9)
#cs1 = map.pcolormesh(x,y,m3newp/m3newa,cmap=plt.cm.gist_earth,vmin=0,vmax=10)
plt.axis('off')

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
#cbar.ax.set_xticklabels(['0', '0.01', '0.05','0.1','0.5','1','2','3','4','5','6'])  # horizontal colorbar

cbar.ax.tick_params(labelsize=16)

ax1 = fig.add_subplot(222)
ax1.set_title("Maize ISAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

ma1= ma.masked_where(ma1[:,:]<=0.0,ma1)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#ncvar_y=maskoceans(x1,y1,ncvar_y)
cs1 = map.pcolormesh(x,y,ma1*m3newa*10000/gridarea,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma1,cmap=plt.cm.gist_earth,vmin=0,vmax=10)

plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)
#cbar.ax.set_xticklabels(['0', '0.01', '0.05','0.1','0.5','1','2','3','4','5','6'])  # horizontal colorbar

bounds=[-0.1,0.0,0.03,0.06,0.09,0.12,0.15,0.2,0.3,0.4,0.5,0.6,0.7]
norm = colors.BoundaryNorm(bounds, cmap.N)

#bounds=[0,0.01,0.05,0.1,0.5,1,1.5,2]
norm = colors.BoundaryNorm(bounds, cmap.N)

ax1 = fig.add_subplot(223)
ax1.set_title("Soybean SPAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

m3newps=maskoceans(x,y,m3newps)
m3newas=maskoceans(x,y,m3newas)
cs1 = map.pcolormesh(x,y,m3newps*10000/gridarea,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,m3newps/m3newas,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)

ax1 = fig.add_subplot(224)
ax1.set_title("Soybean ISAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

ma2= ma.masked_where(ma2[:,:]<=0.0,ma2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#ncvar_y=maskoceans(x1,y1,ncvar_y)
cs1 = map.pcolormesh(x,y,ma2*m3newas*10000/gridarea,cmap=cmap,norm=norm)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)

plt.savefig('isam_spam.jpg',bbox_inches='tight')

plt.show()


