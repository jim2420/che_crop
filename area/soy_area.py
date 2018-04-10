from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime


country=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
coun = country.variables['MASK_Country'][:,:]



    
region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][:,:,:]
maitemp = region1.variables['soy_temp'][:,:,:]
maitropi=region1.variables['soy_trop_irrig'][:,:,:]
maitempi=region1.variables['soy_temp_irrig'][:,:,:]
gridarea = region1.variables['area'][:,:]
landarea = region1.variables['landfrac'][:,:]

maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)

maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)

maizetro=(maitrop+maitropi)*gridarea*landarea
maizetem=(maitemp+maitempi)*gridarea*landarea
maizetor=(maitrop+maitemp)*gridarea*landarea
maizetoi=(maitropi+maitempi)*gridarea*landarea
maizeto =(maitrop+maitemp+maitropi+maitempi)*gridarea*landarea
rcp85=N.sum(maizeto,axis=(1,2))
rcp85r=N.sum(maizetor,axis=(1,2))
rcp85i=N.sum(maizetoi,axis=(1,2))




region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][:,:,:]
maitemp = region1.variables['soy_temp'][:,:,:]
maitropi=region1.variables['soy_trop_irrig'][:,:,:]
maitempi=region1.variables['soy_temp_irrig'][:,:,:]
gridarea = region1.variables['area'][:,:]
landarea = region1.variables['landfrac'][:,:]
maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)

maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)

maizetro=(maitrop+maitropi)*gridarea*landarea
maizetem=(maitemp+maitempi)*gridarea*landarea
maizetor=(maitrop+maitemp)*gridarea*landarea
maizetoi=(maitropi+maitempi)*gridarea*landarea
maizeto =(maitrop+maitemp+maitropi+maitempi)*gridarea*landarea
rcp45=N.sum(maizeto,axis=(1,2))
rcp45r=N.sum(maizetor,axis=(1,2))
rcp45i=N.sum(maizetoi,axis=(1,2))




region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][:,:,:]
maitemp = region1.variables['soy_temp'][:,:,:]
maitropi=region1.variables['soy_trop_irrig'][:,:,:]
maitempi=region1.variables['soy_temp_irrig'][:,:,:]
gridarea = region1.variables['area'][:,:]
landarea = region1.variables['landfrac'][:,:]
maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)

maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)

maizetro=(maitrop+maitropi)*gridarea*landarea
maizetem=(maitemp+maitempi)*gridarea*landarea
maizetor=(maitrop+maitemp)*gridarea*landarea
maizetoi=(maitropi+maitempi)*gridarea*landarea
maizeto =(maitrop+maitemp+maitropi+maitempi)*gridarea*landarea
his=N.sum(maizeto,axis=(1,2))
hisr=N.sum(maizetor,axis=(1,2))
hisi=N.sum(maizetoi,axis=(1,2))


rcp45=N.concatenate((his,rcp45),axis=0)
rcp85=N.concatenate((his,rcp85),axis=0)
rcp45i=N.concatenate((hisi,rcp45i),axis=0)
rcp85i=N.concatenate((hisi,rcp85i),axis=0)
rcp45r=N.concatenate((hisr,rcp45r),axis=0)
rcp85r=N.concatenate((hisr,rcp85r),axis=0)

rcp45=(rcp45-rcp45[104])/rcp45[104]*100
rcp85=(rcp85-rcp85[104])/rcp85[104]*100
his=(his-his[104])/his[104]*100

rcp45i=(rcp45i-rcp45i[104])/rcp45i[104]*100
rcp85i=(rcp85i-rcp85i[104])/rcp85i[104]*100
hisi=(hisi-hisi[104])/hisi[104]*100

rcp45r=(rcp45r-rcp45r[104])/rcp45r[104]*100
rcp85r=(rcp85r-rcp85r[104])/rcp85r[104]*100
hisr=(hisr-hisr[104])/hisr[104]*100


#print rcp45.shape

a1=1901
a2=2101

fig = plt.figure(figsize=(10,12))


ax = fig.add_subplot(311)
xx=range(a1,a2)

x1=range(1901,2006)
xdates1 = [datetime.datetime.strptime(str(int(date)),'%Y') for date in x1]

xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

ax.plot_date(xdates,rcp85,"b-",label="RCP 8.5",linewidth=5.0)
ax.plot_date(xdates,rcp45,"r:",label="RCP 4.5",linewidth=5.0)
ax.plot_date(xdates1,his,"k-",linewidth=5.0)

ax.set_ylim([-55,55])
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
plt.title("Soybean Irrigated + Rainfed",fontsize=20)
#plt.xlabel("Year",fontsize=18)
plt.ylabel("Area changes (%)",fontsize=20)
plt.tick_params(axis='both',labelsize=20)
plt.setp(plt.gca(), xticklabels=[])

#plt.axis('off')
#leg = plt.legend(loc=2,fancybox=True, fontsize=16)
#leg.get_frame().set_alpha(0.5)



ax = fig.add_subplot(312)
ax.plot_date(xdates,rcp85i,"b-",label="RCP 8.5",linewidth=5.0)
ax.plot_date(xdates,rcp45i,"r:",label="RCP 4.5",linewidth=5.0)
ax.plot_date(xdates1,hisi,"k-",linewidth=5.0)

ax.set_ylim([-55,55])
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
plt.title("Soybean Irrigated",fontsize=20)
#plt.xlabel("Year",fontsize=18)
plt.ylabel("Area changes (%)",fontsize=20)
plt.tick_params(axis='both',labelsize=20)
#plt.axis('off')
plt.setp(plt.gca(), xticklabels=[])

ax = fig.add_subplot(313)
ax.plot_date(xdates,rcp85r,"b-",label="RCP 8.5",linewidth=5.0)
ax.plot_date(xdates,rcp45r,"r:",label="RCP 4.5",linewidth=5.0)
ax.plot_date(xdates1,hisr,"k-",linewidth=5.0)

ax.set_ylim([-55,55])
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
plt.title("Soybean Rainfed",fontsize=20)
plt.xlabel("Year",fontsize=20)
plt.ylabel("Area changes (%)",fontsize=20)
plt.tick_params(axis='both',labelsize=20)


plt.savefig('soy_area_g.png',dpi=300,bbox_inches='tight')
plt.show()

