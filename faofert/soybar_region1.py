from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
area1=NetCDFFile('/project/projectdirs/m1602/datasets4.full/surfdata_05x05.nc','r')
mask = area1.variables['REGION_MASK_CRU_NCEP'][:,:]

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
ncvar_maize = nclu.variables['soyy'][0,:,:]





region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][95:105,:,:]
maitemp = region1.variables['soy_temp'][95:105,:,:]
maitropi=region1.variables['soy_trop_irrig'][95:105,:,:]
maitempi=region1.variables['soy_temp_irrig'][95:105,:,:]
landfrac =region1.variables['landfrac'][:,:]
gridarea = region1.variables['area'][:,:]
lonisam1=region1.variables['lon'][:]
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

 



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ynew = dat.variables['totalyield'][95:105,:,:]
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2ynew = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3ynew = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constco2_irrig_fert_0.5x0.5.nc','r')
iyield4ynew = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constclim_irrig_fert_0.5x0.5.nc','r')
iyield5ynew = dat5.variables['totalyield'][95:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield6ynew = dat6.variables['totalyield'][95:105,:,:]


iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)



maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonisam1,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonisam1,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonisam1,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonisam1,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonisam1,start=True)
landfrac1,lonisam2=shiftgrid(0.5,landfrac,lonisam1,start=True)
gridarea1,lonisam2=shiftgrid(0.5,gridarea,lonisam1,start=True)




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
rmask=N.zeros((10,360,720))
m3maize=N.zeros((10,360,720))
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
                        rmask[i,x,y]=mask[x,y]
                        m3maize[i,x,y]=ncvar_maize[x,y] 
			maizeto2[i,x,y]=maizeto1[i,x,y]
                        maizeto2r[i,x,y]=maizeto1r[i,x,y]
                        maizeto2i[i,x,y]=maizeto1i[i,x,y]
                        maizete2r[i,x,y]=maizete1r[i,x,y]
                        maizete2i[i,x,y]=maizete1i[i,x,y]
			landfrac2[i,x,y]=landfrac1[x,y]
                        gridarea2[i,x,y]=gridarea1[x,y]

iyield1ynew= ma.masked_where(m3maize<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(m3maize<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(m3maize<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(m3maize<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(m3maize<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(m3maize<=0.,iyield6ynew)

ii=N.zeros((9))
nn=N.zeros((9))
co=N.zeros((9))
cli=N.zeros((9))
sii=N.zeros((9))
snn=N.zeros((9))
sco=N.zeros((9))
scli=N.zeros((9))


for i in range (1,9):
	maizeto3=maizeto2
	if i==4 or i==5:
	        maizeto3=ma.masked_where(rmask>5.0,maizeto3)
                maizeto3=ma.masked_where(rmask<4.0,maizeto3)

	else:
		maizeto3=ma.masked_where(rmask!=i,maizeto3)

	allynew=N.average(iyield1ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allyinew=N.average(iyield2ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allynnew=N.average(iyield3ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allycnew=N.average(iyield4ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allyclinew=N.average(iyield5ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))

	bb=(allynew-allyinew)/allynew*100
	cc=(allynew-allynnew)/allynew*100
	dd=(allynew-allycnew)/allynew*100
	ee=(allynew-allyclinew)/allynew*100

	ii[i]=N.average(bb)
	nn[i]=N.average(cc)
	co[i]=N.average(dd)
	cli[i]=N.average(ee)

	sii[i]=N.std(bb)
	snn[i]=N.std(cc)
	sco[i]=N.std(dd)
	scli[i]=N.std(ee)


allynew=N.average(iyield1ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allyinew=N.average(iyield2ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allynnew=N.average(iyield3ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allycnew=N.average(iyield4ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allyclinew=N.average(iyield5ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))

bb=(allynew-allyinew)/allynew*100
cc=(allynew-allynnew)/allynew*100
dd=(allynew-allycnew)/allynew*100
ee=(allynew-allyclinew)/allynew*100

ii[0]=N.average(bb)
nn[0]=N.average(cc)
co[0]=N.average(dd)
cli[0]=N.average(ee)

sii[0]=N.std(bb)
snn[0]=N.std(cc)
sco[0]=N.std(dd)
scli[0]=N.std(ee)

name=["Global","NA","SA","EU","Africa","Africa","USSR","China","SSEA"]

fig = plt.figure(figsize=(18,10))
#plt.rc('font', weight='bold')
plt.subplots_adjust(hspace=1)

for idx in xrange(1):
	ax = fig.add_subplot(1, 1, idx+1)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

	n_groups = 1
	plt.ylim(-5, 65)
	ax.set_yticks([0,15,30,45,60 ])

	index = N.arange(n_groups)
	bar_width = 0.01
	opacity = 1.0
#        rects3 = plt.bar(index+0.05+bar_width*3, cli[idx], bar_width, yerr=scli[idx],
#                 alpha=opacity,color='blue',
#                 label='Climate')
#        rects2 = plt.bar(index+0.05+bar_width*2, co[idx], bar_width, yerr=sco[idx],
#                 alpha=opacity,color='black',
#                 label='Climate+CO$_{2}$')
#        rects1 = plt.bar(index+bar_width+0.05, nn[idx], bar_width,yerr=snn[idx],
#                 alpha=opacity,color='green',
#                 label='Climate+CO$_{2}$+NF')
#	rects0 = plt.bar(index+0.05, ii[idx], bar_width,yerr=sii[idx],
#	         alpha=1.0,color='red',
#	         label='Climate+CO$_{2}$+NF+Irrigation')


#        rects2 = plt.bar(index+0.05+bar_width*2, co[idx], bar_width, yerr=sco[idx],
#                 alpha=opacity,color='black',
#                 label='CO$_{2}$')
#        rects3 = plt.bar(index+0.05+bar_width*3, cli[idx], bar_width, yerr=scli[idx],
#                 alpha=opacity,color='blue',
#                 label='Climate')
#        rects1 = plt.bar(index+bar_width+0.05, nn[idx], bar_width,yerr=snn[idx],
#                 alpha=opacity,color='red',
#                 label='Irrigation')
#        rects0 = plt.bar(index+0.05, ii[idx], bar_width,yerr=sii[idx],
#                 alpha=1.0,color='green',
#                 label='NF')

        rects2 = plt.bar(index+0.05+bar_width*2, co[idx], bar_width, yerr=sco[idx],
                 alpha=opacity,color='red',
                 label='Base')
        rects3 = plt.bar(index+0.05+bar_width*3, cli[idx], bar_width, yerr=scli[idx],
                 alpha=opacity,color='blue',
                 label='S_Cli')
        rects1 = plt.bar(index+bar_width+0.05, nn[idx], bar_width,yerr=snn[idx],
                 alpha=opacity,color='black',
                 label='S_CO$_{2}$')
        rects0 = plt.bar(index+0.05, ii[idx], bar_width,yerr=sii[idx],
                 alpha=1.0,color='green',
                 label='S_N')



	plt.tight_layout()
	leg=plt.legend(loc=1,fontsize=16,ncol=4)
	leg.get_frame().set_alpha(0.5)
	leg.get_frame().set_facecolor('none')
        leg.get_frame().set_edgecolor('none')

	#plt.subplots_adjust(hspace=0.5)
	plt.tick_params(
	    axis='x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='off',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='off') # labels along the bottom edge are off
#	plt.axis('off')

#	ax.text(.55,.85,'{0}'.format(name[idx]),
#        	horizontalalignment='center',
#        	transform=ax.transAxes)
#	plt.title('{0}'.format(name[idx]))
#plt.tight_layout()
	ax.tick_params(labelsize=10)

plt.savefig('soy_test12.png')
plt.show()
