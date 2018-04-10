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
import datetime
from matplotlib.dates import DateFormatter
from scipy import stats


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





#lon2,lat2 = N.meshgrid(gridlon,gridlat)



region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
#maitrop = region1.variables['soy_trop'][99,:,:]
#maitemp = region1.variables['soy_temp'][99,:,:]
#maitropi=region1.variables['soy_trop_irrig'][99,:,:]
#maitempi=region1.variables['soy_temp_irrig'][99,:,:]
maitrop = region1.variables['maize_trop'][0:105,:,:]
maitemp = region1.variables['maize_temp'][0:105,:,:]
maitropi=region1.variables['maize_trop_irrig'][0:105,:,:]
maitempi=region1.variables['maize_temp_irrig'][0:105,:,:]


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
maizetropo=maitrop+maitropi
maizetempo=maitemp+maitempi
maizetoavg = N.average(maizeto[95:105,:,:],axis=0)

 



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1y = N.average(dat.variables['totalyield'][95:105,:,:],axis=0)
iyield1ynew = dat.variables['totalyield'][0:105,:,:]
#print iyield1ynew.shape
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2y = N.average(dat2.variables['totalyield'][95:105,:,:],axis=0)
iyield2ynew = dat2.variables['totalyield'][0:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3y = N.average(dat3.variables['totalyield'][95:105,:,:],axis=0)
iyield3ynew = dat3.variables['totalyield'][0:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constco2_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield4y = N.average(dat4.variables['totalyield'][95:105,:,:],axis=0)
iyield4ynew = dat4.variables['totalyield'][0:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constclim_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield5y = N.average(dat5.variables['totalyield'][95:105,:,:],axis=0)
iyield5ynew = dat5.variables['totalyield'][0:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
iyield6y = N.average(dat6.variables['totalyield'][95:105,:,:],axis=0)
iyield6ynew = dat6.variables['totalyield'][0:105,:,:]



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ynewn = dat.variables['total_gfna'][0:105,:,:]
iyield1ynewi = dat.variables['total_gWS'][0:105,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2ynewn = dat2.variables['total_gfna'][0:105,:,:]
iyield2ynewi = dat2.variables['total_gWS'][0:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3ynewn = dat3.variables['total_gfna'][0:105,:,:]
iyield3ynewi = dat3.variables['total_gWS'][0:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_constco2_irrig_fert_0.5x0.5.nc','r')
iyield4ynewn = dat4.variables['total_gfna'][0:105,:,:]
iyield4ynewi = dat4.variables['total_gWS'][0:105,:,:]


dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_constclim_irrig_fert_0.5x0.5.nc','r')
iyield5ynewn = dat5.variables['total_gfna'][0:105,:,:]
iyield5ynewi = dat5.variables['total_gWS'][0:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1_n/maizetemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
iyield6ynewn = dat6.variables['total_gfna'][0:105,:,:]
iyield6ynewi = dat6.variables['total_gWS'][0:105,:,:]




iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)

iyield1ynewi= ma.masked_where(iyield1ynew<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where(iyield2ynew<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where(iyield3ynew<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where(iyield4ynew<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where(iyield5ynew<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where(iyield6ynew<=0.,iyield6ynewi)

iyield1ynewn= ma.masked_where(iyield1ynew<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where(iyield2ynew<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where(iyield3ynew<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where(iyield4ynew<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where(iyield5ynew<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where(iyield6ynew<=0.,iyield6ynewn)




iyield1,lonisam1 = shiftgrid(180.5,iyield1y,lonisam,start=False)
iyield2,lonisam1 = shiftgrid(180.5,iyield2y,lonisam,start=False)
iyield3,lonisam1 = shiftgrid(180.5,iyield3y,lonisam,start=False)
iyield4,lonisam1 = shiftgrid(180.5,iyield4y,lonisam,start=False)
iyield5,lonisam1 = shiftgrid(180.5,iyield5y,lonisam,start=False)
iyield6,lonisam1 = shiftgrid(180.5,iyield6y,lonisam,start=False)
#print lonisam
#print lonisam1
maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonisam1,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonisam1,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonisam1,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonisam1,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonisam1,start=True)
lon2,lat2 = N.meshgrid(lonisam1,latisam)
print lonisam2
print latisam
print lonisam2.shape
print latisam.shape

iyield1ynewn= ma.masked_where(maizeto1<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where(maizeto1<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where(maizeto1<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where(maizeto1<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where(maizeto1<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where(maizeto1<=0.,iyield6ynewn)

iyield1ynewi= ma.masked_where(maizeto1<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where(maizeto1<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where(maizeto1<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where(maizeto1<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where(maizeto1<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where(maizeto1<=0.,iyield6ynewi)

iyield1ynew= ma.masked_where(maizeto1<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(maizeto1<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(maizeto1<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(maizeto1<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(maizeto1<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(maizeto1<=0.,iyield6ynew)


iyield1ynew=ma.filled(iyield1ynew, fill_value=0.)
iyield2ynew=ma.filled(iyield2ynew, fill_value=0.)
iyield3ynew=ma.filled(iyield3ynew, fill_value=0.)
iyield4ynew=ma.filled(iyield4ynew, fill_value=0.)
iyield5ynew=ma.filled(iyield5ynew, fill_value=0.)
iyield6ynew=ma.filled(iyield6ynew, fill_value=0.)

iyield1ynewn=ma.filled(iyield1ynewn, fill_value=-1)
iyield2ynewn=ma.filled(iyield2ynewn, fill_value=-1)
iyield3ynewn=ma.filled(iyield3ynewn, fill_value=-1)
iyield4ynewn=ma.filled(iyield4ynewn, fill_value=-1)
iyield5ynewn=ma.filled(iyield5ynewn, fill_value=-1)
iyield6ynewn=ma.filled(iyield6ynewn, fill_value=-1)


iyield1ynewi=ma.filled(iyield1ynewi, fill_value=-1)
iyield2ynewi=ma.filled(iyield2ynewi, fill_value=-1)
iyield3ynewi=ma.filled(iyield3ynewi, fill_value=-1)
iyield4ynewi=ma.filled(iyield4ynewi, fill_value=-1)
iyield5ynewi=ma.filled(iyield5ynewi, fill_value=-1)
iyield6ynewi=ma.filled(iyield6ynewi, fill_value=-1)


iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)


iyield1ynewi= ma.masked_where(iyield1ynewi<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where(iyield2ynewi<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where(iyield3ynewi<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where(iyield4ynewi<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where(iyield5ynewi<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where(iyield6ynewi<=0.,iyield6ynewi)

iyield1ynewn= ma.masked_where(iyield1ynewn<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where(iyield2ynewn<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where(iyield3ynewn<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where(iyield4ynewn<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where(iyield5ynewn<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where(iyield6ynewn<=0.,iyield6ynewn)

#print iyield6ynewn
#iyield6ynewn=iyield6ynewn[N.logical_not(N.isnan(iyield6ynewn))]
#iyield6ynewn = ma.masked_array(iyield6ynewn, N.isnan(iyield6ynewn))
iyield1ynewn[N.isinf(iyield1ynewn)]=-1
iyield1ynewn[N.isnan(iyield1ynewn)]=-1
iyield1ynewn= ma.masked_where(iyield1ynewn<=0.0,iyield1ynewn)


iyield2ynewn[N.isinf(iyield2ynewn)]=-1
iyield2ynewn[N.isnan(iyield2ynewn)]=-1
iyield2ynewn= ma.masked_where(iyield2ynewn<=0.0,iyield2ynewn)

iyield3ynewn[N.isinf(iyield3ynewn)]=-1
iyield3ynewn[N.isnan(iyield3ynewn)]=-1
iyield3ynewn= ma.masked_where(iyield3ynewn<=0.0,iyield3ynewn)

iyield4ynewn[N.isinf(iyield4ynewn)]=-1
iyield4ynewn[N.isnan(iyield4ynewn)]=-1
iyield4ynewn= ma.masked_where(iyield4ynewn<=0.0,iyield4ynewn)

iyield5ynewn[N.isinf(iyield5ynewn)]=-1
iyield5ynewn[N.isnan(iyield5ynewn)]=-1
iyield5ynewn= ma.masked_where(iyield5ynewn<=0.0,iyield5ynewn)

iyield6ynewn[N.isinf(iyield6ynewn)]=-1
iyield6ynewn[N.isnan(iyield6ynewn)]=-1
iyield6ynewn= ma.masked_where(iyield6ynewn<=0.0,iyield6ynewn)

iyield1ynewi[N.isinf(iyield1ynewi)]=-1
iyield1ynewi[N.isnan(iyield1ynewi)]=-1
iyield1ynewi= ma.masked_where(iyield1ynewi<=0.0,iyield1ynewi)

iyield2ynewi[N.isinf(iyield2ynewi)]=-1
iyield2ynewi[N.isnan(iyield2ynewi)]=-1
iyield2ynewi= ma.masked_where(iyield2ynewi<=0.0,iyield2ynewi)

iyield3ynewi[N.isinf(iyield3ynewi)]=-1
iyield3ynewi[N.isnan(iyield3ynewi)]=-1
iyield3ynewi= ma.masked_where(iyield3ynewi<=0.0,iyield3ynewi)

iyield4ynewi[N.isinf(iyield4ynewi)]=-1
iyield4ynewi[N.isnan(iyield4ynewi)]=-1
iyield4ynewi= ma.masked_where(iyield4ynewi<=0.0,iyield4ynewi)

iyield5ynewi[N.isinf(iyield5ynewi)]=-1
iyield5ynewi[N.isnan(iyield5ynewi)]=-1
iyield5ynewi= ma.masked_where(iyield5ynewi<=0.0,iyield5ynewi)

iyield6ynewi[N.isinf(iyield6ynewi)]=-1
iyield6ynewi[N.isnan(iyield6ynewi)]=-1
iyield6ynewi= ma.masked_where(iyield6ynewi<=0.0,iyield6ynewi)



#print iyield6ynewn

#print iyield1ynew.shape
maizeto2=N.zeros((105,360,720))
maizeto2r=N.zeros((105,360,720))
maizeto2i=N.zeros((105,360,720))
maizete2r=N.zeros((105,360,720))
maizete2i=N.zeros((105,360,720))

for i in range(0,105):
	for x in range(0,360):
		for y in range(0,720):
			maizeto2[i,x,y]=maizeto1[i,x,y]
                        maizeto2r[i,x,y]=maizeto1r[i,x,y]
                        maizeto2i[i,x,y]=maizeto1i[i,x,y]
                        maizete2r[i,x,y]=maizete1r[i,x,y]
                        maizete2i[i,x,y]=maizete1i[i,x,y]


cwnall=1-abs((1-iyield6ynewi)-(1-iyield6ynewn))
cwnall1=cwnall[1:105,:,:]
aal=iyield6ynew-iyield1ynew
aa1=aal[1:105,:,:]
print aa1.shape
print cwnall1.shape

cwnalls=N.average(cwnall1,axis=0)
aa1s=N.average(aa1,axis=0)


temp =[]
#corrcoefMatrix_gfzjpl = [[0 for i in range(len(lonisam2))] for j in range(len(latisam))] 
corrcoefMatrix_gfzjpl=N.zeros((360,720))
for x in range(len(latisam)):
    for y in range(len(lonisam2)):
	slope, intercept, r_value, p_value, std_err = stats.linregress(aa1[:,x,y],cwnall1[:,x,y])
        corrcoefMatrix_gfzjpl[x,y] = r_value
#	print r_value
#        temp = N.corrcoef(aa1[:,x,y],cwnall1[:,x,y])
##        print temp[0,1]
##        corrcoefMatrix_gfzjpl[x][y] = temp[0,1]
#        corrcoefMatrix_gfzjpl[x,y] = temp[0,1]
##corrcoefMatrix_gfzjpl = N.squeeze(N.asarray(corrcoefMatrix_gfzjpl))




#rr1=stats.pearsonr(aa1,cwnall1)
print corrcoefMatrix_gfzjpl.shape

cwnalls,lonisam1 = shiftgrid(180.5,cwnalls,lonisam,start=False)

corrcoefMatrix_gfzjpl,lonisam1 = shiftgrid(180.5,corrcoefMatrix_gfzjpl,lonisam,start=False)
print lonisam1


fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(211)


#aa=allyinnew-allynew
#aa1=aa[1:105]

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map.drawcoastlines()
x,y = map(lon2,lat2)
map.drawcoastlines()
#map.drawcountries()
map.drawmapboundary()
corrcoefMatrix_gfzjpl= ma.masked_where(maizetoavg<=0.0,corrcoefMatrix_gfzjpl)
cc=N.average(corrcoefMatrix_gfzjpl**2,weights=maizetoavg)
print cc
cs = map.pcolormesh(x,y,corrcoefMatrix_gfzjpl**2,cmap=plt.cm.jet,vmin=0,vmax=1)

plt.axis('off')
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%",extend='both')
cbar.ax.tick_params(labelsize=15)
ax = fig.add_subplot(212)
aa1s= ma.masked_where(aa1s==0.0,aa1s)


ax.scatter(cwnalls,aa1s,c='black',alpha=1)
plt.xlabel('CM$_{WN}$ (-)',fontsize=16)
plt.ylabel('Yield gap (t/ha)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)
plt.xlim(0,1.1)
plt.ylim(-5,5)




plt.tight_layout()

plt.savefig('maize_his_inco_papersp.png')
plt.show()
