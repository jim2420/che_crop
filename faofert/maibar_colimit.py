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


iyield1y= ma.masked_where(iyield1y<=0.,iyield1y)
iyield2y= ma.masked_where(iyield2y<=0.,iyield2y)
iyield3y= ma.masked_where(iyield3y<=0.,iyield3y)
iyield4y= ma.masked_where(iyield4y<=0.,iyield4y)
iyield5y= ma.masked_where(iyield5y<=0.,iyield5y)
iyield6y= ma.masked_where(iyield6y<=0.,iyield6y)


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
iyield1ynewn[N.isnan(iyield1ynewn)]=-1
iyield1ynewn= ma.masked_where(iyield1ynewn<=0.0,iyield1ynewn)

iyield2ynewn[N.isnan(iyield2ynewn)]=-1
iyield2ynewn= ma.masked_where(iyield2ynewn<=0.0,iyield2ynewn)

iyield3ynewn[N.isnan(iyield3ynewn)]=-1
iyield3ynewn= ma.masked_where(iyield3ynewn<=0.0,iyield3ynewn)

iyield4ynewn[N.isnan(iyield4ynewn)]=-1
iyield4ynewn= ma.masked_where(iyield4ynewn<=0.0,iyield4ynewn)

iyield5ynewn[N.isnan(iyield5ynewn)]=-1
iyield5ynewn= ma.masked_where(iyield5ynewn<=0.0,iyield5ynewn)

iyield6ynewn[N.isnan(iyield6ynewn)]=-1
iyield6ynewn= ma.masked_where(iyield6ynewn<=0.0,iyield6ynewn)


iyield1ynewi[N.isnan(iyield1ynewi)]=-1
iyield1ynewi= ma.masked_where(iyield1ynewi<=0.0,iyield1ynewi)

iyield2ynewi[N.isnan(iyield2ynewi)]=-1
iyield2ynewi= ma.masked_where(iyield2ynewi<=0.0,iyield2ynewi)

iyield3ynewi[N.isnan(iyield3ynewi)]=-1
iyield3ynewi= ma.masked_where(iyield3ynewi<=0.0,iyield3ynewi)

iyield4ynewi[N.isnan(iyield4ynewi)]=-1
iyield4ynewi= ma.masked_where(iyield4ynewi<=0.0,iyield4ynewi)

iyield5ynewi[N.isnan(iyield5ynewi)]=-1
iyield5ynewi= ma.masked_where(iyield5ynewi<=0.0,iyield5ynewi)

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

allynew=N.average(iyield1ynew,weights=maizeto2,axis=(1,2))
allyinew=N.average(iyield2ynew,weights=maizeto2,axis=(1,2))
allynnew=N.average(iyield3ynew,weights=maizeto2,axis=(1,2))
allycnew=N.average(iyield4ynew,weights=maizeto2,axis=(1,2))
allyclinew=N.average(iyield5ynew,weights=maizeto2,axis=(1,2))
allyinnew=N.average(iyield6ynew,weights=maizeto2,axis=(1,2))

allynewi=N.average(iyield1ynewi,weights=maizeto2,axis=(1,2))
allyinewi=N.average(iyield2ynewi,weights=maizeto2,axis=(1,2))
allynnewi=N.average(iyield3ynewi,weights=maizeto2,axis=(1,2))
allycnewi=N.average(iyield4ynewi,weights=maizeto2,axis=(1,2))
allyclinewi=N.average(iyield5ynewi,weights=maizeto2,axis=(1,2))
allyinnewi=N.average(iyield6ynewi,weights=maizeto2,axis=(1,2))

allynewn=N.average(iyield1ynewn,weights=maizeto2,axis=(1,2))
allyinewn=N.average(iyield2ynewn,weights=maizeto2,axis=(1,2))
allynnewn=N.average(iyield3ynewn,weights=maizeto2,axis=(1,2))
allycnewn=N.average(iyield4ynewn,weights=maizeto2,axis=(1,2))
allyclinewn=N.average(iyield5ynewn,weights=maizeto2,axis=(1,2))
allyinnewn=N.average(iyield6ynewn,weights=maizeto2,axis=(1,2))


aa=allyinnew-allynew
a1=allyinnewi[1:105]
b1=allyinnewn[1:105]
twn=1-allyinnewi+1-allyinnewn
mwn=N.maximum((1-allyinnewi),(1-allyinnewn))
cwn=1-abs((1-allyinnewi)-(1-allyinnewn))


ctwn=cwn/twn
cmwn=cwn/mwn
aa1=aa[1:105]
twn1=twn[1:105]
mwn1=mwn[1:105]
cwn1=cwn[1:105]
ctwn1=ctwn[1:105]
cmwn1=cmwn[1:105]


#finalnew=(allynew-allyinnew)/allynew*100
#finalnnew=(allynew-allynnew)/allynew*100
#finalinew=(allynew-allyinew)/allynew*100

#finalnew1=(allynew-allyinnew)/allyinnew*100
#finalnnew1=(allynnew-allyinnew)/allyinnew*100
#finalinew1=(allyinew-allyinnew)/allyinnew*100

#finalclinew=(allynew-allyclinew)/allynew*100
#finalcnew=(allynew-allycnew)/allynew*100


#finalstd=(N.std(allyinew)/N.std(allynew),N.std(allynnew)/N.std(allynew),N.std(allycnew)/N.std(allynew),N.std(allyclinew)/N.std(allynew))

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(221)


#plt.xlim(0, 0.7)
#plt.ylim(0, 0.7)
#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()
#plt.tick_params(axis='both',labelsize=15)
ax.scatter(a1,aa1,c='black',alpha=1)
plt.xlabel('Water limitation (-)',fontsize=16)
plt.ylabel('Yield gap (t/ha)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)
#plt.xlim(0, 2)
#plt.ylim(-3, 1)

slope, intercept, r_value, p_value, std_err = stats.linregress(twn1,aa1)
line = slope*twn1+intercept
ax.plot(twn1,line,linewidth=2.0)
rr=r_value**2
ax.annotate('P = 0.06', xy=(0.04, 0.10),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.01, 0.03),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2


ax = fig.add_subplot(222)
ax.scatter(b1,aa1,c='black',alpha=1)
plt.xlabel('Nitrogen limitation (-)',fontsize=16)
plt.ylabel('Yield gap (t/ha)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)

slope, intercept, r_value, p_value, std_err = stats.linregress(mwn1,aa1)
line = slope*mwn1+intercept
ax.plot(mwn1,line,linewidth=2.0)
rr=r_value**2
ax.annotate('P < 0.0001', xy=(0.04, 0.10),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.01, 0.03),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2




ax = fig.add_subplot(223)
ax.scatter(cwn1,aa1,c='black',alpha=1)
plt.xlabel('Co-limitation of water and nitrogen (-)',fontsize=16)
plt.ylabel('Yield difference (t/ha)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)
#ax.set_xticklabels([0.75,0.80,0.85,0.90,0.95,1.0])
ax.set_xlim([0.75, 1.0])
plt.axis([0.75, 1, -4, -1.5])

x = [0.75,0.80,0.85,0.90,0.95,1.0]


slope, intercept, r_value, p_value, std_err = stats.linregress(cwn1,aa1)
line = slope*cwn1+intercept
ax.plot(cwn1,line,linewidth=2.0)
rr=r_value**2
ax.annotate('P < 0.0001', xy=(0.04, 0.10),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.01, 0.03),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2



ax = fig.add_subplot(224)
ax.scatter(cmwn1,aa1,c='black',alpha=1)
plt.xlabel('CM$_{WN}$ (-)',fontsize=16)
plt.ylabel('Yield gap (t/ha)',fontsize=16)
plt.tick_params(axis='both',labelsize=16)

slope, intercept, r_value, p_value, std_err = stats.linregress(cmwn1,aa1)
line = slope*cmwn1+intercept
ax.plot(cmwn1,line,linewidth=2.0)
rr=r_value**2
ax.annotate('P < 0.0001', xy=(0.04, 0.10),xycoords='axes fraction',fontsize=18,color='blue')
ax.annotate('R$^2$ = {:04.2f}'.format(rr), xy=(0.01, 0.03),xycoords='axes fraction',fontsize=18,color='blue')

print 'p',p_value,"r-squared:", r_value**2


plt.tight_layout()

plt.savefig('maize_his_inco_paper11.png')
plt.show()
