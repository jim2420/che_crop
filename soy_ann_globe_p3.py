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
#print iizumi
coun = country.variables['MASK_Country'][:,:]


region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][:,:,:]
maitemp = region1.variables['soy_temp'][:,:,:]
maitropi=region1.variables['soy_trop_irrig'][:,:,:]
maitempi=region1.variables['soy_temp_irrig'][:,:,:]
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
    

clm3n2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_soyscaleiyield_fertfao_new1.nc','r')
isamcrumai2 = clm3n2.variables['yield'][:,:,:]
isamcrumai2= ma.masked_where(isamcrumai2<=0,isamcrumai2)
isamcrumai2=ma.filled(isamcrumai2, fill_value=0.)
isamcrumai2[N.isnan(isamcrumai2)] = 0
#isamcrumai2= ma.masked_where(isamcrumai2<=0,isamcrumai2)
isamcrumai2=ma.filled(isamcrumai2, fill_value=0.)

yynew=N.zeros((105,360))

for t1 in range(0,105):
	for xx in range(0,360):
		harea=0
		yieldic2=0
		for yy in range(0,720):
			if maizeto[t1,xx,yy]>0.0:
				harea=maizeto[t1,xx,yy]*gridarea[xx,yy]+ harea
				yieldic2=(isamcrumai2[t1,xx,yy]*maizeto[t1,xx,yy]*gridarea[xx,yy])+yieldic2
		if harea>0.0:
			yynew[t1,xx]=yieldic2/harea

yynew[N.isnan(yynew)] = -1
yynew= ma.masked_where(yynew<=0,yynew)

#print yynew[99,:].shape
#print yynew[99,:]


fig = plt.figure(figsize=(15,6))

tt = N.arange(-89.75, 90., 0.5)
ax = fig.add_subplot(111)
plt.plot(tt,N.average(yynew[90:100,:],axis=0),'r-')
plt.plot(tt,N.average(yynew[60:70,:],axis=0),'k-')
plt.plot(tt,N.average(yynew[10:20,:],axis=0),'b-')

leg=plt.legend(['1991-2000','1961-1970','1911-1920'],fontsize=18)
leg.get_frame().set_alpha(0.5)
#plt.xticks(N.arange(tt.min(), tt.max(), 30))
#plt.axis([-90, 90, 0, 15])
ax.set_xticks((-75,-60,-45,-30,-15,0,15,30,45,60,75))
my_xticks = ['-75 S','-60 S','-45 S','-30 S','-15 S', 'EQ','15 N','30 N','45 N','60 N','75 N']
#plt.xticks(x, my_xticks)
ax.set_xticklabels(my_xticks, fontsize=18)

plt.xlabel("Latitude ",fontsize=18)
plt.ylabel("Soybean yields (t/ha)",fontsize=18)
plt.tick_params(axis='both',labelsize=18)


plt.savefig('soy_1961_2016_lat.png',dpi=600,bbox_inches='tight')
plt.show()


