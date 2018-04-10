from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import numpy.ma as ma
from statsmodels.stats.weightstats import DescrStatsW
import matplotlib.colors as colors


nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_m = nclu.variables['maizeData'][0,0,:,:]
ncvar_y = nclu.variables['maizeData'][0,1,:,:]
ncvar_a = nclu.variables['maizeData'][0,4,:,:]
ncvar_p = nclu.variables['maizeData'][0,5,:,:]

nclu1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_ms = nclu1.variables['soybeanData'][0,0,:,:]
ncvar_ys = nclu1.variables['soybeanData'][0,1,:,:]
ncvar_as = nclu1.variables['soybeanData'][0,4,:,:]
ncvar_ps = nclu1.variables['soybeanData'][0,5,:,:]

latm3 = nclu1.variables['latitude'][:]
lonm3 = nclu1.variables['longitude'][:]



region=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_maiscaleiyield_fertfao_new1.nc','r')
ma1 = region.variables['yield'][96:103,:,:]
#ma1 = region.variables['yield'][99,:,:]

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_soyscaleiyield_fertfao_new1.nc','r')
ma2 = region1.variables['yield'][96:103,:,:]
#ma2 = region1.variables['yield'][99,:,:]

ma1=N.average(ma1,axis=0)
ma2=N.average(ma2,axis=0)

latmask = region.variables['lat'][:]
lonmask = region.variables['lon'][:]

latm3_new=N.flipud(latm3)
ncvar_m=N.flipud(ncvar_m)
ncvar_y=N.flipud(ncvar_y)
ncvar_a=N.flipud(ncvar_a)
ncvar_p=N.flipud(ncvar_p)

ncvar_ms=N.flipud(ncvar_ms)
ncvar_ys=N.flipud(ncvar_ys)
ncvar_as=N.flipud(ncvar_as)
ncvar_ps=N.flipud(ncvar_ps)



lon,lat = N.meshgrid(lonmask,latmask)
lon1,lat1 = N.meshgrid(lonm3,latm3_new)


m3newy=N.zeros((360,720))
m3newp=N.zeros((360,720))
m3newa=N.zeros((360,720))
m3newm=N.zeros((360,720))

m3newys=N.zeros((360,720))
m3newps=N.zeros((360,720))
m3newas=N.zeros((360,720))
m3newms=N.zeros((360,720))

for x in range(0,360):
    for y in range(0,720):
        a1=x*6
        a2=(x+1)*6
        b1=y*6
        b2=(y+1)*6
        for j in range(a1,a2):
            for i in range(b1,b2):
                m3newm[x,y]=ncvar_m[j,i]+m3newm[x,y]
                m3newy[x,y]=ncvar_y[j,i]+m3newy[x,y]
                m3newa[x,y]=ncvar_a[j,i]+m3newa[x,y]
                m3newp[x,y]=ncvar_p[j,i]+m3newp[x,y]

                m3newms[x,y]=ncvar_ms[j,i]+m3newms[x,y]
                m3newys[x,y]=ncvar_ys[j,i]+m3newys[x,y]
                m3newas[x,y]=ncvar_as[j,i]+m3newas[x,y]
                m3newps[x,y]=ncvar_ps[j,i]+m3newps[x,y]

        m3newm[x,y]=m3newm[x,y]/36
        m3newy[x,y]=m3newy[x,y]/36

        m3newms[x,y]=m3newms[x,y]/36
        m3newys[x,y]=m3newys[x,y]/36

m3newm= ma.masked_where(m3newm[:,:]<=0.0,m3newm)
m3newy= ma.masked_where(m3newy[:,:]<=0.0,m3newy)
m3newa= ma.masked_where(m3newa[:,:]<=0.0,m3newa)
m3newp= ma.masked_where(m3newp[:,:]<=0.0,m3newp)

m3newms= ma.masked_where(m3newms[:,:]<=0.0,m3newms)
m3newys= ma.masked_where(m3newys[:,:]<=0.0,m3newys)
m3newas= ma.masked_where(m3newas[:,:]<=0.0,m3newas)
m3newps= ma.masked_where(m3newps[:,:]<=0.0,m3newps)

m3newp= ma.masked_where(ma1[:,:]<=0.0,m3newp)
ma1= ma.masked_where(m3newp[:,:]<=0.0,ma1)

m3newps= ma.masked_where(ma2[:,:]<=0.0,m3newps)
ma2= ma.masked_where(m3newps[:,:]<=0.0,ma2)


fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(221)
ax1.set_title("Maize M3",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
x,y = map(lon,lat)
m3newp=maskoceans(x,y,m3newp)
m3newa=maskoceans(x,y,m3newa)
cs1 = map.pcolormesh(x,y,m3newp/m3newa,cmap=plt.cm.Greens,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=12)
#cs1 = map.pcolormesh(x,y,m3newp/m3newa,cmap=plt.cm.gist_earth,vmin=0,vmax=10)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=12)

ax1 = fig.add_subplot(222)
ax1.set_title("Maize ISAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

ma1= ma.masked_where(ma1[:,:]<=0.0,ma1)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#ncvar_y=maskoceans(x1,y1,ncvar_y)
cs1 = map.pcolormesh(x,y,ma1,cmap=plt.cm.Greens,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=12)
#cs1 = map.pcolormesh(x,y,ma1,cmap=plt.cm.gist_earth,vmin=0,vmax=10)

plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=12)

ax1 = fig.add_subplot(223)
ax1.set_title("Soybean M3",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

m3newps=maskoceans(x,y,m3newps)
m3newas=maskoceans(x,y,m3newas)
cs1 = map.pcolormesh(x,y,m3newps/m3newas,cmap=plt.cm.Greens,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=2)
#cs1 = map.pcolormesh(x,y,m3newps/m3newas,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=12)

ax1 = fig.add_subplot(224)
ax1.set_title("Soybean ISAM",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-62, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

ma2= ma.masked_where(ma2[:,:]<=0.0,ma2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#ncvar_y=maskoceans(x1,y1,ncvar_y)
cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.Greens,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=2)
#cs1 = map.pcolormesh(x,y,ma2,cmap=plt.cm.gist_earth,vmin=0,vmax=5)
plt.axis('off')
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=12)

plt.savefig('isam_m3_paperf.jpg',bbox_inches='tight')

plt.show()


