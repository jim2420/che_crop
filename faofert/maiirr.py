from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
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


region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][99,:,:]
maitemp = region1.variables['maize_temp'][99,:,:]
maitropi=region1.variables['maize_trop_irrig'][99,:,:]
maitempi=region1.variables['maize_temp_irrig'][99,:,:]
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
maizeto = maitrop+maitemp+maitropi+maitempi

 




dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/maihis_irr_fertfao/output/irr_maihis_irr_fertfao.nc','r')
iyield1y = N.average(dat.variables['irrigation'][95:105,0,:,:],axis=0)
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/maihis_fertfao/output/irr_maihis_fertfao.nc','r')
iyield2y = N.average(dat2.variables['irrigation'][95:105,0,:,:],axis=0)
dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/maihis_irr/output/irr_maihis_irr.nc','r')
iyield3y = N.average(dat3.variables['irrigation'][95:105,0,:,:],axis=0)
dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/maihis_constco2_irr_fertfao/output/irr_maihis_constco2_irr_fertfao.nc','r')
iyield4y = N.average(dat4.variables['irrigation'][95:105,0,:,:],axis=0)
dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/maihis_constcli_irr_fertfao/output/irr_maihis_constcli_irr_fertfao.nc','r')
iyield5y = N.average(dat5.variables['irrigation'][95:105,0,:,:],axis=0)


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield6y = N.average(dat6.variables['totalyield'][95:105,:,:],axis=0)

#iyield1,lonisam1 = shiftgrid(180.5,iyield1y,lonisam,start=False)
#iyield2,lonisam1 = shiftgrid(180.5,iyield2y,lonisam,start=False)
#iyield3,lonisam1 = shiftgrid(180.5,iyield3y,lonisam,start=False)
#iyield4,lonisam1 = shiftgrid(180.5,iyield4y,lonisam,start=False)
#iyield5,lonisam1 = shiftgrid(180.5,iyield5y,lonisam,start=False)

#iyield1=ma.filled(iyield1, fill_value=-99.)
#iyield2=ma.filled(iyield2, fill_value=-99.)
#iyield3=ma.filled(iyield3, fill_value=-99)
#iyield4=ma.filled(iyield4, fill_value=-99.)
#iyield5=ma.filled(iyield5, fill_value=-99.)
#print iyield1

iyield1y= ma.masked_where(iyield1y<=0.,iyield1y)
iyield2y= ma.masked_where(iyield2y<=0.,iyield2y)
iyield3y= ma.masked_where(iyield3y<=0.,iyield3y)
iyield4y= ma.masked_where(iyield4y<=0.,iyield4y)
iyield5y= ma.masked_where(iyield5y<=0.,iyield5y)
iyield6y= ma.masked_where(iyield6y<=0.,iyield6y)


iyield1,lonisam1 = shiftgrid(180.5,iyield1y,lonisam,start=False)
iyield2,lonisam1 = shiftgrid(180.5,iyield2y,lonisam,start=False)
iyield3,lonisam1 = shiftgrid(180.5,iyield3y,lonisam,start=False)
iyield4,lonisam1 = shiftgrid(180.5,iyield4y,lonisam,start=False)
iyield5,lonisam1 = shiftgrid(180.5,iyield5y,lonisam,start=False)
iyield6,lonisam1 = shiftgrid(180.5,iyield6y,lonisam,start=False)


iyield1= ma.masked_where(maizeto<=0.,iyield1)
iyield2= ma.masked_where(maizeto<=0.,iyield2)
iyield3= ma.masked_where(maizeto<=0.,iyield3)
iyield4= ma.masked_where(maizeto<=0.,iyield4)
iyield5= ma.masked_where(maizeto<=0.,iyield5)
iyield6= ma.masked_where(maizeto<=0.,iyield6)


iyield1=ma.filled(iyield1, fill_value=0.)
iyield2=ma.filled(iyield2, fill_value=0.)
iyield3=ma.filled(iyield3, fill_value=0.)
iyield4=ma.filled(iyield4, fill_value=0.)
iyield5=ma.filled(iyield5, fill_value=0.)
iyield6=ma.filled(iyield6, fill_value=0.)

iyield1= ma.masked_where(iyield1<=0.,iyield1)
iyield2= ma.masked_where(iyield2<=0.,iyield2)
iyield3= ma.masked_where(iyield3<=0.,iyield3)
iyield4= ma.masked_where(iyield4<=0.,iyield4)
iyield5= ma.masked_where(iyield5<=0.,iyield5)
iyield6= ma.masked_where(iyield6<=0.,iyield6)




lon2,lat2 = N.meshgrid(lonisam1,latisam)



fig = plt.figure(figsize=(20,15))
ax2 = fig.add_subplot(421)
#ax2.set_title("Yield (t/ha)",fontsize=20)


map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map.drawcoastlines()
x,y = map(lon2,lat2)


map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#print iyield2
cs = map.pcolormesh(x,y,iyield1,cmap=plt.cm.jet,vmin=0,vmax=1000)

#cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(422)
#ax2.set_title("CESM-ISAM Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(iyield1-iyield2)/iyield2*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(423)
#ax2.set_title("CESM-ISAM Nscale Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield3 = ma.masked_where(iyield3<=0,iyield3)
iyield3 = ma.masked_where(iarea<=0,iyield3)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield1,cmap=plt.cm.jet,vmin=0,vmax=1000)
#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(424)
#ax2.set_title("CESM-ISAM Maize Yield Fertall (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield4 = ma.masked_where(iyield4<=0,iyield4)
iyield4 = ma.masked_where(iarea<=0,iyield4)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(iyield1-iyield3)/iyield3*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

#cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')
ax2 = fig.add_subplot(425)
#ax2.set_title("NCEP-ISAM Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
iyield5 = ma.masked_where(iyield5<=0,iyield5)
iyield5 = ma.masked_where(iarea<=0,iyield5)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield1,cmap=plt.cm.jet,vmin=0,vmax=1000)

#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(426)
#ax2.set_title("NCEP-ISAM HS Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(iyield1-iyield4)/iyield4*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

#cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')
ax2 = fig.add_subplot(427)
#ax2.set_title("NCEP-ISAM Nscale Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,iyield1/1000,cmap=plt.cm.jet,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax2 = fig.add_subplot(428)
#ax2.set_title("NCEP-ISAM Nscale HS Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(iyield1-iyield5)/iyield5*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='right',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('maizeisamirr_1996_2005.jpg',dpi=300,bbox_inches='tight')
plt.show()

