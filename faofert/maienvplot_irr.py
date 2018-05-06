from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
from statsmodels.stats.weightstats import DescrStatsW
from scipy.stats import ttest_ind
from matplotlib.markers import TICKDOWN
import matplotlib.patches as mpatches

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon1 = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/spammaize_isam.nc','r')
nclu1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/spamsoy_isam.nc','r')

ncvar_maize = nclu.variables['maiy_irrigated'][:,:]
ncvar_soy = nclu1.variables['soyy_irrigated'][:,:]
areamaize = nclu.variables['maia_irrigated'][:,:]
areasoy = nclu1.variables['soya_irrigated'][:,:]
ncvar_maize1 = nclu.variables['maiya_irrigated'][:,:]


latnc = nclu.variables['lat'][:]
lonnc = nclu.variables['lon'][:]



region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')

maitrop = N.average(region1.variables['maize_trop'][103:105,:,:],axis=0)
maitemp = N.average(region1.variables['maize_temp'][103:105,:,:],axis=0)
maitropi=N.average(region1.variables['maize_trop_irrig'][103:105,:,:],axis=0)
maitempi=N.average(region1.variables['maize_temp_irrig'][103:105,:,:],axis=0)

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
maizeto = maitrop+maitemp+maitropi+maitempi


soytrop = N.average(region1.variables['soy_trop'][103:105,:,:],axis=0)
soytemp = N.average(region1.variables['soy_temp'][103:105,:,:],axis=0)
soytropi=N.average(region1.variables['soy_trop_irrig'][103:105,:,:],axis=0)
soytempi=N.average(region1.variables['soy_temp_irrig'][103:105,:,:],axis=0)

soytrop=ma.masked_where(soytrop<=0,soytrop)
soytrop=ma.filled(soytrop, fill_value=0.)
soytemp=ma.masked_where(soytemp<=0,soytemp)
soytemp=ma.filled(soytemp, fill_value=0.)

soytropi=ma.masked_where(soytropi<=0,soytropi)
soytropi=ma.filled(soytropi, fill_value=0.)
soytempi=ma.masked_where(soytempi<=0,soytempi)
soytempi=ma.filled(soytempi, fill_value=0.)
soyzetor=soytrop+soytemp
soytoi=soytropi+soytempi
soyto = soytrop+soytemp+soytropi+soytempi





dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ys = N.average(dat.variables['totalyield'][103:106,:,:],axis=0)
iyield1yns = dat.variables['totalyield'][103:106,:,:]

latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
iyield1ys,gridlon1 = shiftgrid(180.5,iyield1ys,gridlon,start=False)

iyield1ys= ma.masked_where(iyield1ys<=0.,iyield1ys)
ncvar_soy= ma.masked_where(ncvar_soy<=0.,ncvar_soy)
areasoy= ma.masked_where(areasoy<=0.,areasoy)

iyield1ys= ma.masked_where(ncvar_soy<=0.,iyield1ys)
ncvar_soy= ma.masked_where(iyield1ys<=0.,ncvar_soy)
iyield1ys= ma.masked_where(ncvar_soy<=0.,iyield1ys)
ncvar_soy= ma.masked_where(iyield1ys<=0.,ncvar_soy)

areasoy= ma.masked_where(iyield1ys<=0.,areasoy)

iyield1yns= ma.masked_where(iyield1yns<=0.,iyield1yns)
iyield1ys= ma.masked_where(soyto<=0.,iyield1ys)
ncvar_soy= ma.masked_where(soyto<=0.,ncvar_soy)
areasoy= ma.masked_where(soyto<=0.,areasoy)




dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1y = N.average(dat.variables['totalyield'][103:106,:,:],axis=0)
iyield1yn = dat.variables['totalyield'][103:106,:,:]

iyield1y,gridlon1 = shiftgrid(180.5,iyield1y,gridlon,start=False)

iyield1y= ma.masked_where(iyield1y<=0.,iyield1y)
ncvar_maize= ma.masked_where(ncvar_maize<=0.,ncvar_maize)
areamaize= ma.masked_where(areamaize<=0.,areamaize)

iyield1y= ma.masked_where(ncvar_maize<=0.,iyield1y)
ncvar_maize= ma.masked_where(iyield1y<=0.,ncvar_maize)
iyield1y= ma.masked_where(ncvar_maize<=0.,iyield1y)
ncvar_maize= ma.masked_where(iyield1y<=0.,ncvar_maize)

areamaize= ma.masked_where(iyield1y<=0.,areamaize)

iyield1yn= ma.masked_where(iyield1yn<=0.,iyield1yn)
iyield1y= ma.masked_where(maizeto<=0.,iyield1y)
ncvar_maize= ma.masked_where(maizeto<=0.,ncvar_maize)
areamaize= ma.masked_where(maizeto<=0.,areamaize)






lon2,lat2 = N.meshgrid(gridlon1,latisam)
cmap = plt.cm.terrain_r
bounds=[-0.1,0.0,0.02,0.04,0.06,0.08,0.1,0.5,1,1.5,2]
norm = colors.BoundaryNorm(bounds, cmap.N)



fig = plt.figure(figsize=(20,10))

ax2 = fig.add_subplot(221)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

x,y = map(lon2,lat2)


map.drawcoastlines(color='gray')
#map.drawcountries()
#map.drawmapboundary()

cs1 = map.pcolormesh(x,y,iyield1y*areamaize*10000/gridarea,cmap=cmap,norm=norm)
plt.axis('off')
plt.title('ISAM Maize (t grains / grid cell ha)',fontsize=16)
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)

ax2 = fig.add_subplot(222)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')



map.drawcoastlines(color='gray')
#map.drawcountries()
#map.drawmapboundary()
#cs1 = map.pcolormesh(x,y,ncvar_maize1,cmap=cmap,norm=norm)
cs1 = map.pcolormesh(x,y,ncvar_maize*areamaize*10000/gridarea,cmap=cmap,norm=norm)
plt.axis('off')
plt.title('SPAM Maize (t grains / grid cell ha)',fontsize=16)
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)


bounds=[-0.1,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
norm = colors.BoundaryNorm(bounds, cmap.N)


ax2 = fig.add_subplot(223)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines(color='grey')
cs1 = map.pcolormesh(x,y,iyield1ys*areasoy*10000/gridarea,cmap=cmap,norm=norm)
plt.title('ISAM Soybean (t grains / grid cell ha)',fontsize=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)
plt.axis('off')

ax2 = fig.add_subplot(224)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines(color='grey')
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,ncvar_soy*areasoy*10000/gridarea,cmap=cmap,norm=norm)
plt.title('SPAM Soybean (t grains / grid cell ha)',fontsize=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%",ticks=bounds,extend='max')
cbar.ax.tick_params(labelsize=16)


plt.axis('off')

plt.savefig('spam_2004_2006_irr.jpg',dpi=600,bbox_inches='tight')
plt.show()

