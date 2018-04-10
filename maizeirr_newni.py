from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
from scipy.stats import ttest_ind

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][99,:,:]
maitemp = region1.variables['maize_temp'][99,:,:]
maitropi = region1.variables['maize_trop_irrig'][99,:,:]
maitempi = region1.variables['maize_temp_irrig'][99,:,:]
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
clma=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
clman=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_rf_fert_0.5x0.5.nc','r')
clmai=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_nofert_0.5x0.5.nc','r')
clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_fert_0.5x0.5.nc','r')

clmtropn = clman.variables['irrigation'][99,:,:]#2010-2019
clmtropfn=clmtropn
clmtropi = clmai.variables['irrigation'][99,:,:]#2010-2019
clmtropfi=clmtropi
clmtrop = clma.variables['irrigation'][99,:,:]#2010-2019
clmtropf=clmtrop
clmtropa = clm.variables['irrigation'][99,:,:]#2050-2059
clmtropfa=clmtropa


clm1a=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
clm1an=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
clm1ai=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')


clmtempn = clm1an.variables['irrigation'][99,:,:]
clmtempfn=clmtempn
clmtempi = clm1ai.variables['irrigation'][99,:,:]
clmtempfi=clmtempi
clmtemp = clm1a.variables['irrigation'][99,:,:]
clmtempf=clmtemp
clmtempa = clm1.variables['irrigation'][99,:,:]
clmtempfa=clmtempa

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfa=N.flipud(clmtropfa)
clmtempfa=N.flipud(clmtempfa)
clmtropfi=N.flipud(clmtropfi)
clmtempfi=N.flipud(clmtempfi)
clmtropfn=N.flipud(clmtropfn)
clmtempfn=N.flipud(clmtempfn)


clmtropf= ma.masked_where(maitoatrop<=0,clmtropf)
clmtempf= ma.masked_where(maitoatemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfi= ma.masked_where(maitoatrop<=0,clmtropfi)
clmtempfi= ma.masked_where(maitoatemp<=0,clmtempfi)
clmtropfi=ma.filled(clmtropfi, fill_value=0.)
clmtempfi=ma.filled(clmtempfi, fill_value=0.)

clmtropfn= ma.masked_where(maitoatrop<=0,clmtropfn)
clmtempfn= ma.masked_where(maitoatemp<=0,clmtempfn)
clmtropfn=ma.filled(clmtropfn, fill_value=0.)
clmtempfn=ma.filled(clmtempfn, fill_value=0.)

clmtropfa= ma.masked_where(maitoatrop<=0,clmtropfa)
clmtempfa= ma.masked_where(maitoatemp<=0,clmtempfa)
clmtropfa=ma.filled(clmtropfa, fill_value=0.)
clmtempfa=ma.filled(clmtempfa, fill_value=0.)

clmhisi=clmtropfi+clmtempfi
clmhis=clmtropf+clmtempf
clmfuture=clmtropfa+clmtempfa
clmhisn=clmtropfn+clmtempfn

clmhis= ma.masked_where(clmhis[:,:]<=0,clmhis)
clmfuture= ma.masked_where(clmfuture[:,:]<=0,clmfuture)
clmhisi= ma.masked_where(clmhisi[:,:]<=0,clmhisi)
clmhisn= ma.masked_where(clmhisn[:,:]<=0,clmhisn)


clmhist=clmtrop+clmtemp
clmfutt=clmtropa+clmtempa

tc, pTc = ttest_ind(clmhist,clmfutt, axis = 0, equal_var = False)

tc=N.flipud(tc)
pTc=N.flipud(pTc)

yieldclm=clmfuture-clmhis
yieldclm= ma.masked_where(yieldclm==0.,yieldclm)

#yieldclm1= ma.masked_where( pTc[:,:]>0.1,yieldclm)

yieldfa= N.zeros((1, 360, 720))
yieldf2a= N.zeros((1, 360, 720))
yieldf2ai= N.zeros((1, 360, 720))
yieldf2an= N.zeros((1, 360, 720))

yieldf= N.zeros((1, 360, 720))
yieldf2= N.zeros((1, 360, 720))
yieldf2i= N.zeros((1, 360, 720))
yieldf2n= N.zeros((1, 360, 720))

years = range(2000, 2001)
years2 = range(2000,2001)

for i, year in enumerate(years):
   
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_fert/output/hmaizehis_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr/output/hmaizehis_irr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr_fert/output/hmaizehis_irr_fert.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]


    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_fertrop/output/hmaizehis_fertrop.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr/output/hmaizehis_irr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr_fertrop/output/hmaizehis_irr_fertrop.bgp-yearly_crop_{0}.nc".format(year), mode='r')
 
    yield1 = base.variables["irrigation"][0,:,:]
    yieldf[i, :, :] = yield1
    yield1a = basea.variables["irrigation"][0,:,:]
    yieldfa[i, :, :] = yield1a

    
for i, year1 in enumerate(years2):
    #    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2n = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_fert/output/hmaizehis_fert.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr/output/hmaizehis_irr.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr_fert/output/hmaizehis_irr_fert.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    yield2 = base2.variables["irrigation"][0,:,:]
    yieldf2[i, :, :] = yield2
    yield2i = base2i.variables["irrigation"][0,:,:]
    yieldf2i[i, :, :] = yield2i
    yield2n = base2n.variables["irrigation"][0,:,:]
    yieldf2n[i, :, :] = yield2n


    #base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis/output/hmaizehis.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2an = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_fertrop/output/hmaizehis_fertrop.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2ai = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr/output/hmaizehis_irr.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/his/heat/maihis_irr_fertrop/output/hmaizehis_irr_fertrop.bgp-yearly_crop_{0}.nc".format(year1), mode='r')

    yield2a = base2a.variables["irrigation"][0,:,:]
    yieldf2a[i, :, :] = yield2a

    yield2ai = base2ai.variables["irrigation"][0,:,:]
    yieldf2ai[i, :, :] = yield2ai

    yield2an = base2an.variables["irrigation"][0,:,:]
    yieldf2an[i, :, :] = yield2an



yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)
yielda2i=N.average(yieldf2i,axis=0)
yielda2n=N.average(yieldf2n,axis=0)


yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)
yield_new2i,lona11 = shiftgrid(180.5,yielda2i,lona1,start=False)
yield_new2n,lona11 = shiftgrid(180.5,yielda2n,lona1,start=False)


yieldaa=N.average(yieldfa,axis=0)
yielda2a=N.average(yieldf2a,axis=0)
yielda2ai=N.average(yieldf2ai,axis=0)
yielda2an=N.average(yieldf2an,axis=0)


yield_newa,lona11 = shiftgrid(180.5,yieldaa,lona1,start=False)
yield_new2a,lona11 = shiftgrid(180.5,yielda2a,lona1,start=False)
yield_new2ai,lona11 = shiftgrid(180.5,yielda2ai,lona1,start=False)
yield_new2an,lona11 = shiftgrid(180.5,yielda2an,lona1,start=False)


lon2,lat2 = N.meshgrid(lona11,lata1)

yield_new= ma.masked_where(yield_new<=0.,yield_new)
yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)
yield_new2i= ma.masked_where(yield_new2i<=0.,yield_new2i)
yield_new2n= ma.masked_where(yield_new2n<=0.,yield_new2n)

yield_new= ma.masked_where( clmhis[:,:]<=0,yield_new)
yield_new2=ma.masked_where( clmfuture[:,:]<=0,yield_new2)
yield_new2i=ma.masked_where( clmfuture[:,:]<=0,yield_new2i)
yield_new2n=ma.masked_where( clmfuture[:,:]<=0,yield_new2n)

yield_newa= ma.masked_where(yield_newa<=0.,yield_newa)
yield_new2a= ma.masked_where(yield_new2a<=0.,yield_new2a)
yield_new2ai= ma.masked_where(yield_new2ai<=0.,yield_new2ai)
yield_new2an= ma.masked_where(yield_new2an<=0.,yield_new2an)

yield_newa= ma.masked_where( clmhis[:,:]<=0,yield_newa)
yield_new2a=ma.masked_where( clmfuture[:,:]<=0,yield_new2a)
yield_new2ai=ma.masked_where( clmfuture[:,:]<=0,yield_new2ai)
yield_new2an=ma.masked_where( clmfuture[:,:]<=0,yield_new2an)

yield_new= ma.masked_where( maitoatemp[:,:]<=0,yield_new)
yield_new2=ma.masked_where( maitoatemp[:,:]<=0,yield_new2)
yield_new2i=ma.masked_where( maitoatemp[:,:]<=0,yield_new2i)
yield_new2n=ma.masked_where( maitoatemp[:,:]<=0,yield_new2n)

yield_new=ma.filled(yield_new, fill_value=0.)
yield_new2=ma.filled(yield_new2, fill_value=0.)
yield_new2i=ma.filled(yield_new2i, fill_value=0.)
yield_new2n=ma.filled(yield_new2n, fill_value=0.)

yield_newa= ma.masked_where( maitoatrop[:,:]<=0,yield_newa)
yield_new2a=ma.masked_where( maitoatrop[:,:]<=0,yield_new2a)
yield_new2ai=ma.masked_where( maitoatrop[:,:]<=0,yield_new2ai)
yield_new2an=ma.masked_where( maitoatrop[:,:]<=0,yield_new2an)

yield_newa=ma.filled(yield_newa, fill_value=0.)
yield_new2a=ma.filled(yield_new2a, fill_value=0.)
yield_new2ai=ma.filled(yield_new2ai, fill_value=0.)
yield_new2an=ma.filled(yield_new2an, fill_value=0.)


yield_new=yield_new+yield_newa
yield_new2=yield_new2+yield_new2a
yield_new2i=yield_new2i+yield_new2ai
yield_new2n=yield_new2n+yield_new2an

#clmhisi=ma.masked_where(yield_new<=0,clmhisi)
clmhisn=ma.masked_where(yield_new<=0,clmhisn)
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


clmfuture= ma.masked_where(clmfuture<=0.,clmfuture)
clmhis= ma.masked_where(clmhis<=0.,clmhis)
yield_new= ma.masked_where(yield_new<=0.,yield_new)
yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)
clmhisi= ma.masked_where(clmhisi<=0.,clmhisi)
yield_new2i= ma.masked_where(yield_new2i<=0.,yield_new2i)
clmhisn= ma.masked_where(clmhisn<=0.,clmhisn)
yield_new2n= ma.masked_where(yield_new2n<=0.,yield_new2n)


fig = plt.figure(figsize=(12,6))


ax1 = fig.add_subplot(221)
ax1.set_title("CLM irrigation 2000 I+R - I - N compare base (%)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

x,y = map(lon2,lat2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cni=(clmfuture-clmhis)
cn=(clmhisn-clmhis)
ci=(clmhisi-clmhis)


cs = map.pcolormesh(x,y,(clmfuture-clmhisi)/clmhisi*100,cmap=plt.cm.bwr,vmin=-15,vmax=15)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12) 
plt.axis('off')

ax2 = fig.add_subplot(222)
ax2.set_title("ISAM irrigation 2000 I+R - I - N compare base (%)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')
isamni=(yield_new2-yield_new)
isamn=(yield_new2n-yield_new)
isami=(yield_new2i-yield_new)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,(yield_new2-yield_new2i)/yield_new2i*100,cmap=plt.cm.bwr,vmin=-15,vmax=15)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax5 = fig.add_subplot(223)
ax5.set_title("CLM irrigation 2000 I+R - I - N compare base (%)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cni=(clmfuture-clmhis)/clmhis*100
cn=(clmhisn-clmhis)/clmhis*100
ci=(clmhisi-clmhis)/clmhis*100

cs = map.pcolormesh(x,y,cni-(cn+ci),cmap=plt.cm.bwr,vmin=-120,vmax=120)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')



ax5 = fig.add_subplot(224)
ax5.set_title("ISAM irrigation 2000 I+R - I - N compare base (%)",fontsize=18)
#map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(llcrnrlon=-119,llcrnrlat=23,urcrnrlon=-63,urcrnrlat=51,projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
map.drawcoastlines()
#map.drawstates()
#map.drawcountries(color='b')

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
isamni=(yield_new2-yield_new)/yield_new*100
isamn=(yield_new2n-yield_new)/yield_new*100
isami=(yield_new2i-yield_new)/yield_new*100
cs = map.pcolormesh(x,y,isamni-(isamn+isami),cmap=plt.cm.bwr,vmin=-120,vmax=120)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

plt.savefig('maize2000irr_irrnir.jpg',dpi=300,bbox_inches='tight')
plt.show()

