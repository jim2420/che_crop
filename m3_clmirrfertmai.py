from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors


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

maizetrop=maitrop+maitropi
maizetemp=maitemp+maitempi

maizeto = maitrop+maitemp+maitropi+maitempi



ff=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalFertilizer.nc','r')
fert_maitrop = ff.variables['maize_trop_fert'][99,:,:]
fert_maitemp = ff.variables['maize_temp_fert'][99,:,:]




clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtropf = clm.variables['yield'][99,:,:]
#clmtropf=N.average(clmtrop,axis=0)
clmtropfer=clm.variables['fertilizer'][99,:,:]

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtempf = clm1.variables['yield'][99,:,:]
clmtempfer=clm1.variables['fertilizer'][99,:,:]

clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtropfi = clm2.variables['yield'][99,:,:]

clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtempfi = clm3.variables['yield'][99,:,:]



clma=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
clmtropfno = clma.variables['yield'][99,:,:]

clm1a=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
clmtempfno = clm1a.variables['yield'][99,:,:]

clm2a=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_nofert_0.5x0.5.nc','r')
clmtropfnoi = clm2a.variables['yield'][99,:,:]

clm3a=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
clmtempfnoi = clm3a.variables['yield'][99,:,:]



clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfi=N.flipud(clmtropfi)
clmtempfi=N.flipud(clmtempfi)

clmtropfno=N.flipud(clmtropfno)
clmtempfno=N.flipud(clmtempfno)
clmtropfnoi=N.flipud(clmtropfnoi)
clmtempfnoi=N.flipud(clmtempfnoi)


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



area=NetCDFFile('/data/jain1/b/tslin2/crops/gridareahalf.nc','r')
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

fig, ax = plt.subplots(figsize=(8,6))

ax.set_title("M3 Maize Yield (t/ha)",fontsize=20)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)
iyield = ma.masked_where(maizeto<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
iizumy=iyield
cs = map.pcolormesh(x,y,iizumy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')
plt.savefig('m3_maiiizirr.jpg',dpi=300,bbox_inches='tight')

fig = plt.figure(figsize=(20,15))



ax2 = fig.add_subplot(321)
ax2.set_title("CLM Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)

#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)

yield_clmtf=maskoceans(x,y,yield_clmtf)
#yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)
#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
yield_clmtf = ma.masked_where(maizeto<=0,yield_clmtf)

yield_clmtfi=maskoceans(x,y,yield_clmtfi)
#yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)
#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)

clmy=((yield_clmtf*maizetor*gridarea)+(yield_clmtfi*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))
clmy = ma.masked_where(iizumy<=0,clmy)

cs1 = map.pcolormesh(x,y,clmy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')
#print N.max(yield_fine*ncvar_maize1*1000/gridarea)


ax2 = fig.add_subplot(322)
ax2.set_title("CLM-M3 Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,clmy-iizumy,cmap=plt.cm.bwr,vmin=-5,vmax=5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(323)
ax2.set_title("CLM Maize Yield scale-N (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)

#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)

yield_clmtfnew=maskoceans(x,y,yield_clmtfnew)
yield_clmtfnew = ma.masked_where(maizeto<=0,yield_clmtfnew)

yield_clmtfinew=maskoceans(x,y,yield_clmtfinew)
yield_clmtfinew = ma.masked_where(maizeto<=0,yield_clmtfinew)

clmynew=((yield_clmtfnew*maizetor*gridarea)+(yield_clmtfinew*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))
clmynew = ma.masked_where(iizumy<=0,clmynew)

cs1 = map.pcolormesh(x,y,clmynew,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(324)
ax2.set_title("CLM-M3 Maize Yield scale-N (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,clmynew-iizumy,cmap=plt.cm.bwr,vmin=-5,vmax=5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



ax2 = fig.add_subplot(325)
ax2.set_title("CLM-M3 Maize Yield (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(clmy-iizumy)/iizumy*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(326)
ax2.set_title("CLM-M3 Maize Yield scale-N (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(clmynew-iizumy)/iizumy*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('m3_maiiizirr_clm.jpg',dpi=300,bbox_inches='tight')


plt.show()


#fig, ax = plt.subplots(figsize=(6,6))
colors = (0,0,1)
colorsr = (1,0,0)
fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(211)

#ax.plot([0,1000],[0,1000], 'k--',label='1:1')

#ax.scatter(iizumy, isamy, c=colors,alpha=1,label='ISAM')
ax.scatter(iizumy, clmy, c=colorsr,alpha=0.5,label='CLM')
plt.xlim(0, 16)
plt.ylim(0, 16)
ax.plot([0,16],[0,16], 'k--',label='1:1')

#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()

plt.tick_params(axis='both',labelsize=15)
plt.xlabel('M3-Crops (t/ha)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)


ax = fig.add_subplot(212)

#ax.plot([0,1000],[0,1000], 'k--',label='1:1')

ax.scatter(iizumy, clmynew, c=colors,alpha=0.5,label='CLM scale-N')
#ax.scatter(iizumy, clmy, c=colorsr,alpha=0.5,label='CLM')
ax.plot([0,16],[0,16], 'k--',label='1:1')

plt.xlim(0, 16)
plt.ylim(0, 16)
#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()
plt.xticks([])
plt.tick_params(axis='both',labelsize=15)
#plt.xlabel('Iizumi-Crops (g/$\mathregular{m^2}$)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)

plt.savefig('scatterm3_maiiizirr_clm.png',bbox_inches='tight')
#plt.show()

