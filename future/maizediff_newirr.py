from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
from scipy.stats import ttest_ind

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
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
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_rf_fert_0.5x0.5.nc','r')
clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_irrig_fert_0.5x0.5.nc','r')


clmtrop = clm.variables['irrigation'][4:13,:,:]#2010-2019
clmtropf=N.average(clmtrop,axis=0)
clmtropa = clm.variables['irrigation'][44:53,:,:]#2050-2059
clmtropfa=N.average(clmtropa,axis=0)


#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_rf_fert_0.5x0.5.nc','r')
clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_irrig_fert_0.5x0.5.nc','r')



clmtemp = clm1.variables['irrigation'][4:13,:,:]
clmtempf=N.average(clmtemp,axis=0)
clmtempa = clm1.variables['irrigation'][44:53,:,:]
clmtempfa=N.average(clmtempa,axis=0)

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfa=N.flipud(clmtropfa)
clmtempfa=N.flipud(clmtempfa)

clmtropf= ma.masked_where(maitoatrop<=0,clmtropf)
clmtempf= ma.masked_where(maitoatemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfa= ma.masked_where(maitoatrop<=0,clmtropfa)
clmtempfa= ma.masked_where(maitoatemp<=0,clmtempfa)
clmtropfa=ma.filled(clmtropfa, fill_value=0.)
clmtempfa=ma.filled(clmtempfa, fill_value=0.)

clmhis=clmtropf+clmtempf
clmfuture=clmtropfa+clmtempfa
clmhis= ma.masked_where(clmhis[:,:]<=0,clmhis)
clmfuture= ma.masked_where(clmfuture[:,:]<=0,clmfuture)


clmhist=clmtrop+clmtemp
clmfutt=clmtropa+clmtempa

tc, pTc = ttest_ind(clmhist,clmfutt, axis = 0, equal_var = False)

tc=N.flipud(tc)
pTc=N.flipud(pTc)

yieldclm=clmfuture-clmhis
yieldclm= ma.masked_where(yieldclm==0.,yieldclm)

yieldclm1= ma.masked_where( pTc[:,:]>0.1,yieldclm)

yieldfa= N.zeros((10, 360, 720))
yieldf2a= N.zeros((10, 360, 720))
yieldf= N.zeros((10, 360, 720))
yieldf2= N.zeros((10, 360, 720))
years = range(2010, 2020)
years2 = range(2050,2060)

for i, year in enumerate(years):
   
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihisa/output/hmaizehisa.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_fert/output/hmaizehis_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_irr_fert/output/hmaizehis_irr_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]


#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihisa/output/hmaizehisa.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_fertrop/output/hmaizehis_fertrop.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_irr_fertrop/output/hmaizehis_irr_fertrop.bgp-yearly_crop_{0}.nc".format(year), mode='r')
   
    yield1 = base.variables["irrigation"][0,:,:]
    yieldf[i, :, :] = yield1
    yield1a = basea.variables["irrigation"][0,:,:]
    yieldfa[i, :, :] = yield1a

    
for i, year1 in enumerate(years2):
#    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihisa/output/hmaizehisa.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_fert/output/hmaizehis_fert.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_irr_fert/output/hmaizehis_irr_fert.bgp-yearly_crop_{0}.nc".format(year1), mode='r')  
    yield2 = base2.variables["irrigation"][0,:,:]
    yieldf2[i, :, :] = yield2

#    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihisa/output/hmaizehisa.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_fertrop/output/hmaizehis_fertrop.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/rcp85/heat/maihis_irr_fertrop/output/hmaizehis_irr_fertrop.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    yield2a = base2a.variables["irrigation"][0,:,:]
    yieldf2a[i, :, :] = yield2a


yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)


yieldaa=N.average(yieldfa,axis=0)
yielda2a=N.average(yieldf2a,axis=0)

yield_newa,lona11 = shiftgrid(180.5,yieldaa,lona1,start=False)
yield_new2a,lona11 = shiftgrid(180.5,yielda2a,lona1,start=False)


lon2,lat2 = N.meshgrid(lona11,lata1)

yield_new= ma.masked_where(yield_new<=0.,yield_new)
yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)

yield_new= ma.masked_where( clmhis[:,:]<=0,yield_new)
yield_new2=ma.masked_where( clmfuture[:,:]<=0,yield_new2)

yield_newa= ma.masked_where(yield_newa<=0.,yield_newa)
yield_new2a= ma.masked_where(yield_new2a<=0.,yield_new2a)

yield_newa= ma.masked_where( clmhis[:,:]<=0,yield_newa)
yield_new2a=ma.masked_where( clmfuture[:,:]<=0,yield_new2a)

yield_new= ma.masked_where( maitoatemp[:,:]<=0,yield_new)
yield_new2=ma.masked_where( maitoatemp[:,:]<=0,yield_new2)
yield_new=ma.filled(yield_new, fill_value=0.)
yield_new2=ma.filled(yield_new2, fill_value=0.)

yield_newa= ma.masked_where( maitoatrop[:,:]<=0,yield_newa)
yield_new2a=ma.masked_where( maitoatrop[:,:]<=0,yield_new2a)
yield_newa=ma.filled(yield_newa, fill_value=0.)
yield_new2a=ma.filled(yield_new2a, fill_value=0.)


yield_new=yield_new+yield_newa
yield_new2=yield_new2+yield_new2a

clmhis=ma.masked_where(yield_new<=0,clmhis)
clmfuture=ma.masked_where(yield_new2<=0,clmfuture)
yieldclm=clmfuture-clmhis
yieldclm= ma.masked_where(yieldclm==0.,yieldclm)

t, pT = ttest_ind(yieldf, yieldf2, axis = 0, equal_var = False)
t1,lona11 = shiftgrid(180.5,t,lona1,start=False)
pT1,lona11 = shiftgrid(180.5,pT,lona1,start=False)

yieldisam=yield_new2-yield_new

yieldisam1= ma.masked_where( pT1[:,:]>0.1,yieldisam)
yieldisam= ma.masked_where(yieldisam==0.,yieldisam)

fig = plt.figure(figsize=(12,6))


ax1 = fig.add_subplot(221)
ax1.set_title("CLM irrigation difference (mm)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

x,y = map(lon2,lat2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yieldclm,cmap=plt.cm.PiYG,vmin=-2.5,vmax=2.5)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12) 
plt.axis('off')

ax2 = fig.add_subplot(222)
ax2.set_title("ISAM irrigation difference (mm)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yieldisam,cmap=plt.cm.PiYG,vmin=-2.5,vmax=2.5)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax5 = fig.add_subplot(223)
ax5.set_title("CLM irrigation change (%)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yieldclm/clmhis*100,cmap=plt.cm.PiYG,vmin=-30,vmax=30)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')



ax5 = fig.add_subplot(224)
ax5.set_title("ISAM irrigation change (%)",fontsize=18)
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yieldisam/yield_new*100,cmap=plt.cm.PiYG,vmin=-30,vmax=30)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

plt.savefig('scp85maized_irrigation.jpg',dpi=300,bbox_inches='tight')
plt.show()

