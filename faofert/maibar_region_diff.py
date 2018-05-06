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
gridarea1= area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
ncvar_maize = nclu.variables['maizey'][0,:,:]
marea = nclu.variables['maize_area'][0,:,:]





region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][95:105,:,:]
maitemp = region1.variables['maize_temp'][95:105,:,:]
maitropi=region1.variables['maize_trop_irrig'][95:105,:,:]
maitempi=region1.variables['maize_temp_irrig'][95:105,:,:]
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

edat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/equili/maizetemp_historical_constco2_constcli_rf_nofert_0.5x0.5.nc','r')
eiyield1ynew = edat.variables['totalyield'][95:105,:,:]

edat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/equili/maizetemp_historical_constco2_constcli_irrig_nofert_0.5x0.5.nc','r')
eiyield2ynew = edat2.variables['totalyield'][95:105,:,:]

edat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/equili/maizetemp_historical_constco2_constcli_rf_fert_0.5x0.5.nc','r')
eiyield3ynew = edat3.variables['totalyield'][95:105,:,:]

edat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/equili/maizetemp_historical_co2_constcli_rf_nofert_0.5x0.5.nc','r')
eiyield4ynew = edat4.variables['totalyield'][95:105,:,:]

edat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/equili/maizetemp_historical_constco2_cli_rf_nofert_0.5x0.5.nc','r')
eiyield5ynew = edat5.variables['totalyield'][95:105,:,:]
 



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ynew = dat.variables['totalyield'][95:105,:,:]
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2ynew = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3ynew = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constco2_irrig_fert_0.5x0.5.nc','r')
iyield4ynew = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constclim_irrig_fert_0.5x0.5.nc','r')
iyield5ynew = dat5.variables['totalyield'][95:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield6ynew = dat6.variables['totalyield'][95:105,:,:]


iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)

eiyield1ynew= ma.masked_where(eiyield1ynew<=0.,eiyield1ynew)
eiyield2ynew= ma.masked_where(eiyield2ynew<=0.,eiyield2ynew)
eiyield3ynew= ma.masked_where(eiyield3ynew<=0.,eiyield3ynew)
eiyield4ynew= ma.masked_where(eiyield4ynew<=0.,eiyield4ynew)
eiyield5ynew= ma.masked_where(eiyield5ynew<=0.,eiyield5ynew)


maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonisam1,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonisam1,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonisam1,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonisam1,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonisam1,start=True)
landfrac1,lonisam2=shiftgrid(0.5,landfrac,lonisam1,start=True)
#gridarea1,lonisam2=shiftgrid(0.5,gridarea,lonisam1,start=True)

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



eiyield1ynew=ma.filled(eiyield1ynew, fill_value=0.)
eiyield2ynew=ma.filled(eiyield2ynew, fill_value=0.)
eiyield3ynew=ma.filled(eiyield3ynew, fill_value=0.)
eiyield4ynew=ma.filled(eiyield4ynew, fill_value=0.)
eiyield5ynew=ma.filled(eiyield5ynew, fill_value=0.)

eiyield1ynew= ma.masked_where(eiyield1ynew<=0.,eiyield1ynew)
eiyield2ynew= ma.masked_where(eiyield2ynew<=0.,eiyield2ynew)
eiyield3ynew= ma.masked_where(eiyield3ynew<=0.,eiyield3ynew)
eiyield4ynew= ma.masked_where(eiyield4ynew<=0.,eiyield4ynew)
eiyield5ynew= ma.masked_where(eiyield5ynew<=0.,eiyield5ynew)



#print iyield1ynew.shape
mmarea=N.zeros((10,360,720))
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
			mmarea[i,x,y]=marea[x,y]
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

eiyield1ynew= ma.masked_where(m3maize<=0.,eiyield1ynew)
eiyield2ynew= ma.masked_where(m3maize<=0.,eiyield2ynew)
eiyield3ynew= ma.masked_where(m3maize<=0.,eiyield3ynew)
eiyield4ynew= ma.masked_where(m3maize<=0.,eiyield4ynew)
eiyield5ynew= ma.masked_where(m3maize<=0.,eiyield5ynew)



ii=N.zeros((9))
nn=N.zeros((9))
co=N.zeros((9))
cli=N.zeros((9))
sii=N.zeros((9))
snn=N.zeros((9))
sco=N.zeros((9))
scli=N.zeros((9))

eii=N.zeros((9))
enn=N.zeros((9))
eco=N.zeros((9))
ecli=N.zeros((9))
esii=N.zeros((9))
esnn=N.zeros((9))
esco=N.zeros((9))
escli=N.zeros((9))


for i in range (1,9):
	maizeto3=maizeto2
	if i==4 or i==5:
	        maizeto3=ma.masked_where(rmask>5.0,maizeto3)
                maizeto3=ma.masked_where(rmask<4.0,maizeto3)

	else:
		maizeto3=ma.masked_where(rmask!=i,maizeto3)

	allynew=N.average(iyield1ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allyinew=N.average(iyield2ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allynnew=N.average(iyield3ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allycnew=N.average(iyield4ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
	allyclinew=N.average(iyield5ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))

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

        eallynew=N.average(eiyield1ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
        eallyinew=N.average(eiyield2ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
        eallynnew=N.average(eiyield3ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
        eallycnew=N.average(eiyield4ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))
        eallyclinew=N.average(eiyield5ynew*mmarea*10000/gridarea2,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))


        ebb=((allynew-allyinew)-(eallyinew-eallynew))/allynew*100
        ecc=((allynew-allynnew)-(eallynnew-eallynew))/allynew*100
        edd=((allynew-allycnew)-(eallycnew-eallynew))/allynew*100
        eee=((allynew-allyclinew)-(eallyclinew-eallynew))/allynew*100

        eii[i]=N.average(ebb)
        enn[i]=N.average(ecc)
        eco[i]=N.average(edd)
        ecli[i]=N.average(eee)

        esii[i]=N.std(ebb)
        esnn[i]=N.std(ecc)
        esco[i]=N.std(edd)
        escli[i]=N.std(eee)

allynew=N.average(iyield1ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allyinew=N.average(iyield2ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allynnew=N.average(iyield3ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allycnew=N.average(iyield4ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
allyclinew=N.average(iyield5ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))

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



eallynew=N.average(eiyield1ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
eallyinew=N.average(eiyield2ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
eallynnew=N.average(eiyield3ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
eallycnew=N.average(eiyield4ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))
eallyclinew=N.average(eiyield5ynew*mmarea*10000/gridarea2,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))

ebb=((allynew-allyinew)-(eallyinew-eallynew))/allynew*100
ecc=((allynew-allynnew)-(eallynnew-eallynew))/allynew*100
edd=((allynew-allycnew)-(eallycnew-eallynew))/allynew*100
eee=((allynew-allyclinew)-(eallyclinew-eallynew))/allynew*100

eii[0]=N.average(ebb)
enn[0]=N.average(ecc)
eco[0]=N.average(edd)
ecli[0]=N.average(eee)

esii[0]=N.std(ebb)
esnn[0]=N.std(ecc)
esco[0]=N.std(edd)
escli[0]=N.std(eee)



name=["Global","NA","SA","EU","Africa","Africa","USSR","China","SSEA"]

fig = plt.figure(figsize=(4.3,3.5))
plt.rc('font', weight='bold')
plt.subplots_adjust(hspace=1)

for idx in xrange(9):
	ax = fig.add_subplot(3, 3, idx+1)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

	n_groups = 1
	#plt.ylim(-15, 60)
	#ax.set_yticks([-15,0,15,30,45,60 ])

	index = N.arange(n_groups)
	bar_width = 0.01
	opacity = 0.6
        rects3 = plt.bar(index+0.05, ecli[idx], bar_width, yerr=escli[idx],
                 alpha=opacity,color='blue',
                 label='Climate')
	rects2 = plt.bar(index+0.05+bar_width, eco[idx], bar_width, yerr=esco[idx],
                 alpha=opacity,color='black',
                 label='CO$_{2}$')
	rects0 = plt.bar(index+0.05+bar_width*2, eii[idx], bar_width,yerr=esii[idx],
	         alpha=0.9,color='red',
	         label='Irrigation')
	rects1 = plt.bar(index+bar_width*3+0.05, enn[idx], bar_width,yerr=esnn[idx],
	         alpha=opacity,color='green',
	         label='NF')
	plt.tight_layout()

	#plt.subplots_adjust(hspace=0.5)
	plt.tick_params(
	    axis='x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='off',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='off') # labels along the bottom edge are off
#	plt.axis('off')

	ax.text(.55,.85,'{0}'.format(name[idx]),
        	horizontalalignment='center',
        	transform=ax.transAxes)
#	plt.title('{0}'.format(name[idx]))
#plt.tight_layout()
	ax.tick_params(labelsize=10)

plt.savefig('maize_region_diff.png')
plt.show()
