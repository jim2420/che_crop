from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
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
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['soybeanData'][:]
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
maitrop = region1.variables['soy_trop'][95:105,:,:]
maitemp = region1.variables['soy_temp'][95:105,:,:]
maitropi=region1.variables['soy_trop_irrig'][95:105,:,:]
maitempi=region1.variables['soy_temp_irrig'][95:105,:,:]
landfrac =region1.variables['landfrac'][:,:]
gridarea = region1.variables['area'][:,:]
#landfrac=ma.masked_where(landfrac<=0,landfrac)
#gridarea=ma.masked_where(gridarea<=0,gridarea)

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

 

iyield1ynew=N.zeros((10,360,720))


dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1y = N.average(dat.variables['totalyield'][95:105,:,:],axis=0)
iyield1ynew = dat.variables['totalyield'][95:105,:,:]
#print iyield1ynew.shape
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2y = N.average(dat2.variables['totalyield'][95:105,:,:],axis=0)
iyield2ynew = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3y = N.average(dat3.variables['totalyield'][95:105,:,:],axis=0)
iyield3ynew = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constco2_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield4y = N.average(dat4.variables['totalyield'][95:105,:,:],axis=0)
iyield4ynew = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constclim_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield5y = N.average(dat5.variables['totalyield'][95:105,:,:],axis=0)
iyield5ynew = dat5.variables['totalyield'][95:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield6y = N.average(dat6.variables['totalyield'][95:105,:,:],axis=0)
iyield6ynew = dat6.variables['totalyield'][95:105,:,:]

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


iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)



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
landfrac1,lonisam2=shiftgrid(0.5,landfrac,lonisam1,start=True)
gridarea1,lonisam2=shiftgrid(0.5,gridarea,lonisam1,start=True)
#print lonisam2


#iyield1= ma.masked_where(maizeto<=0.,iyield1)
#iyield2= ma.masked_where(maizeto<=0.,iyield2)
#iyield3= ma.masked_where(maizeto<=0.,iyield3)
#iyield4= ma.masked_where(maizeto<=0.,iyield4)
#iyield5= ma.masked_where(maizeto<=0.,iyield5)
#iyield6= ma.masked_where(maizeto<=0.,iyield6)


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

iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)


#print iyield1ynew.shape
maizeto2=N.zeros((10,360,720))
maizeto2r=N.zeros((10,360,720))
maizeto2i=N.zeros((10,360,720))
maizete2r=N.zeros((10,360,720))
maizete2i=N.zeros((10,360,720))
landfrac2=N.zeros((10,360,720))
gridarea2=N.zeros((10,360,720))
for i in range(0,10):
	for x in range(0,360):
		for y in range(0,720):
			maizeto2[i,x,y]=maizeto1[i,x,y]
                        maizeto2r[i,x,y]=maizeto1r[i,x,y]
                        maizeto2i[i,x,y]=maizeto1i[i,x,y]
                        maizete2r[i,x,y]=maizete1r[i,x,y]
                        maizete2i[i,x,y]=maizete1i[i,x,y]
			landfrac2[i,x,y]=landfrac1[x,y]
                        gridarea2[i,x,y]=gridarea1[x,y]
allynew=N.average(iyield1ynew,weights=maizeto2,axis=(1,2))
allyinew=N.average(iyield2ynew,weights=maizeto2,axis=(1,2))
allynnew=N.average(iyield3ynew,weights=maizeto2,axis=(1,2))
allycnew=N.average(iyield4ynew,weights=maizeto2,axis=(1,2))
allyclinew=N.average(iyield5ynew,weights=maizeto2,axis=(1,2))
finalnew=((N.average(allyinew)-N.average(allynew))/N.average(allynew)*-100,(N.average(allynnew)-N.average(allynew))/N.average(allynew)*-100,(N.average(allycnew)-N.average(allynew))/N.average(allynew)*-100,(N.average(allyclinew)-N.average(allynew))/N.average(allynew)*-100)
finalstd=(N.std(allyinew)/N.std(allynew),N.std(allynnew)/N.std(allynew),N.std(allycnew)/N.std(allynew),N.std(allyclinew)/N.std(allynew))



allynewr=N.average(iyield1ynew,weights=maizeto2r,axis=(1,2))
allyinewr=N.average(iyield2ynew,weights=maizeto2r,axis=(1,2))
allynnewr=N.average(iyield3ynew,weights=maizeto2r,axis=(1,2))
allycnewr=N.average(iyield4ynew,weights=maizeto2r,axis=(1,2))
allyclinewr=N.average(iyield5ynew,weights=maizeto2r,axis=(1,2))
finalnewr=((N.average(allyinewr)-N.average(allynewr))/N.average(allynewr)*-100,(N.average(allynnewr)-N.average(allynewr))/N.average(allynewr)*-100,(N.average(allycnewr)-N.average(allynewr))/N.average(allynewr)*-100,(N.average(allyclinewr)-N.average(allynewr))/N.average(allynewr)*-100)
finalstdr=(N.std(allyinewr)/N.std(allynewr),N.std(allynnewr)/N.std(allynewr),N.std(allycnewr)/N.std(allynewr),N.std(allyclinewr)/N.std(allynewr))


allynewi=N.average(iyield1ynew,weights=maizeto2i,axis=(1,2))
allyinewi=N.average(iyield2ynew,weights=maizeto2i,axis=(1,2))
allynnewi=N.average(iyield3ynew,weights=maizeto2i,axis=(1,2))
allycnewi=N.average(iyield4ynew,weights=maizeto2i,axis=(1,2))
allyclinewi=N.average(iyield5ynew,weights=maizeto2i,axis=(1,2))
finalnewi=((N.average(allyinewi)-N.average(allynewi))/N.average(allynewi)*-100,(N.average(allynnewi)-N.average(allynewi))/N.average(allynewi)*-100,(N.average(allycnewi)-N.average(allynewi))/N.average(allynewi)*-100,(N.average(allyclinewi)-N.average(allynewi))/N.average(allynewi)*-100)
finalstdi=(N.std(allyinewi)/N.std(allynewi),N.std(allynnewi)/N.std(allynewi),N.std(allycnewi)/N.std(allynewi),N.std(allyclinewi)/N.std(allynewi))




allynewre=N.average(iyield1ynew,weights=maizete2r,axis=(1,2))
allyinewre=N.average(iyield2ynew,weights=maizete2r,axis=(1,2))
allynnewre=N.average(iyield3ynew,weights=maizete2r,axis=(1,2))
allycnewre=N.average(iyield4ynew,weights=maizete2r,axis=(1,2))
allyclinewre=N.average(iyield5ynew,weights=maizete2r,axis=(1,2))
finalnewre=((N.average(allyinewre)-N.average(allynewre))/N.average(allynewre)*-100,(N.average(allynnewre)-N.average(allynewre))/N.average(allynewre)*-100,(N.average(allycnewre)-N.average(allynewre))/N.average(allynewre)*-100,(N.average(allyclinewre)-N.average(allynewre))/N.average(allynewre)*-100)
finalstdre=(N.std(allyinewre)/N.std(allynewre),N.std(allynnewre)/N.std(allynewre),N.std(allycnewre)/N.std(allynewre),N.std(allyclinewre)/N.std(allynewre))


allynewie=N.average(iyield1ynew,weights=maizete2i,axis=(1,2))
allyinewie=N.average(iyield2ynew,weights=maizete2i,axis=(1,2))
allynnewie=N.average(iyield3ynew,weights=maizete2i,axis=(1,2))
allycnewie=N.average(iyield4ynew,weights=maizete2i,axis=(1,2))
allyclinewie=N.average(iyield5ynew,weights=maizete2i,axis=(1,2))
finalnewie=((N.average(allyinewie)-N.average(allynewie))/N.average(allynewie)*-100,(N.average(allynnewie)-N.average(allynewie))/N.average(allynewie)*-100,(N.average(allycnewie)-N.average(allynewie))/N.average(allynewie)*-100,(N.average(allyclinewie)-N.average(allynewie))/N.average(allynewie)*-100)
finalstdie=(N.std(allyinewie)/N.std(allynewie),N.std(allynnewie)/N.std(allynewie),N.std(allycnewie)/N.std(allynewie),N.std(allyclinewie)/N.std(allynewie))





#print finalnew

#ally=N.average(iyield1,weights=maizeto)
#allyi=N.average(iyield2,weights=maizeto)
#allyn=N.average(iyield3,weights=maizeto)
#allyc=N.average(iyield4,weights=maizeto)
#allycli=N.average(iyield5,weights=maizeto)
#final=(ally,allyi,allyn,allyc,allycli)

xallynew=N.average((iyield1ynew-iyield2ynew)/iyield2ynew*100,weights=maizeto2,axis=(1,2))
xallyinew=N.average((iyield1ynew-iyield3ynew)/iyield3ynew*100,weights=maizeto2,axis=(1,2))
xallynnew=N.average((iyield1ynew-iyield4ynew)/iyield4ynew*100,weights=maizeto2,axis=(1,2))
xallycnew=N.average((iyield1ynew-iyield5ynew)/iyield5ynew*100,weights=maizeto2,axis=(1,2))
xallyi=N.average(xallynew)
xallyn=N.average(xallyinew)
xallyc=N.average(xallynnew)
xallycli=N.average(xallycnew)
final1=(xallyi,xallyn,xallyc,xallycli)
#print final1
maizeto2=ma.masked_where(maizeto2<=0,maizeto2)
maizeto2i=ma.masked_where(maizeto2i<=0,maizeto2i)
maizeto2r=ma.masked_where(maizeto2r<=0,maizeto2r)
maizete2i=ma.masked_where(maizete2i<=0,maizete2i)
maizete2r=ma.masked_where(maizete2r<=0,maizete2r)

tropi=N.sum(maizeto2i*landfrac2*gridarea2)/N.sum(maizeto2*landfrac2*gridarea2)
tropr=N.sum(maizeto2r*landfrac2*gridarea2)/N.sum(maizeto2*landfrac2*gridarea2)
tempi=N.sum(maizete2i*landfrac2*gridarea2)/N.sum(maizeto2*landfrac2*gridarea2)
tempr=N.sum(maizete2r*landfrac2*gridarea2)/N.sum(maizeto2*landfrac2*gridarea2)
#print tropi,tropr,tempi,tempr

print finalnew,finalnewr,finalnewi,finalnewre,finalnewie


fig = plt.figure(figsize=(10,5))
n_groups = 4
plt.ylim(-5,30)
ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.9
opacity = 0.6

rects0 = plt.bar(index+0.05, finalnew, bar_width,
         alpha=0.9,color='gray',
         label='Global')
rects1 = plt.bar(index+0.05, finalnewr, bar_width*tropr,
         alpha=opacity,color='green',
         label='Trop. rainfed')
rects2 = plt.bar(index+0.05+bar_width*tropr, finalnewi, bar_width*tropi,
         alpha=opacity,color='black',
         label='Trop. irrigated')
rects1 = plt.bar(index+0.05+bar_width*tropr+bar_width*tropi, finalnewre, bar_width*tempr,
         alpha=opacity,color='blue',
         label='Temp. rainfed')
rects2 = plt.bar(index+0.05+bar_width*tropr+bar_width*tropi+bar_width*tempr, finalnewie, bar_width*tempi,
         alpha=opacity,color='cyan',
         label='Temp. irrigated')



plt.ylabel('Percentage of soybean yield change',fontsize=18)
plt.tick_params(axis='both',labelsize=18)
plt.xticks(index+bar_width*0.5 , ('Irrigation','N fertilizer','CO$_{2}$','Climate'))
#leg=plt.legend(loc=4)
#leg.get_frame().set_alpha(0.5)

plt.tight_layout()

plt.savefig('soy_his_irrfixed_paper4.png')
plt.show()
