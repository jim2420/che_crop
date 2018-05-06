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
ncvar_maize = nclu.variables['soyy'][0,:,:]
marea = nclu.variables['soy_area'][0,:,:]





region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][0:105,:,:]
maitemp = region1.variables['soy_temp'][0:105,:,:]
maitropi=region1.variables['soy_trop_irrig'][0:105,:,:]
maitempi=region1.variables['soy_temp_irrig'][0:105,:,:]
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


dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/his_cru/heat/fixedirr_fertfao_new/soyhis_irr_fertfao/output/soyhis_irr_fertfao_tt.nc','r')
iyield1ynew = dat.variables['g_Temp'][0:105,7,1,:,:]
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]




iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield1ynew= ma.masked_where(iyield1ynew>=10**20.,iyield1ynew)


maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonisam1,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonisam1,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonisam1,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonisam1,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonisam1,start=True)
landfrac1,lonisam2=shiftgrid(0.5,landfrac,lonisam1,start=True)
#gridarea1,lonisam2=shiftgrid(0.5,gridarea,lonisam1,start=True)

iyield1ynew=ma.filled(iyield1ynew, fill_value=0.)
iyield1ynew[N.isnan(iyield1ynew)] = 0
iyield1ynew[N.isinf(iyield1ynew)] = 0

iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)



#print iyield1ynew.shape
mmarea=N.zeros((105,360,720))
rmask=N.zeros((105,360,720))
m3maize=N.zeros((105,360,720))
maizeto2=N.zeros((105,360,720))
maizeto2r=N.zeros((105,360,720))
maizeto2i=N.zeros((105,360,720))
maizete2r=N.zeros((105,360,720))
maizete2i=N.zeros((105,360,720))
landfrac2=N.zeros((105,360,720))
gridarea2=N.zeros((105,360,720))
for i in range(0,105):
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


ii=N.zeros((105,9))



for i in range (1,9):
	maizeto3=maizeto2
	if i==5:
		i=9
	if i==4 :
	        maizeto3=ma.masked_where(rmask>5.0,maizeto3)
                maizeto3=ma.masked_where(rmask<4.0,maizeto3)

	else:
		print i
		maizeto3=ma.masked_where(rmask!=i,maizeto3)
	if i==9:
		i=5
	allynew=N.average(iyield1ynew,weights=maizeto3*landfrac2*gridarea2,axis=(1,2))

	ii[:,i]=allynew



allynew=N.average(iyield1ynew,weights=maizeto2*landfrac2*gridarea2,axis=(1,2))

ii[:,0]=allynew


name=["Global","NA","SA","EU","Africa","PD","USSR","China","SSEA"]

fig = plt.figure(figsize=(6.3,5.5))
#plt.rc('font', weight='bold')
plt.subplots_adjust(hspace=1)

for idx in xrange(9):
	ax = fig.add_subplot(3, 3, idx+1)
       # ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

	n_groups = 1
#	plt.ylim(-15, 50)
#	ax.set_yticks([-15,0,15,30,45,60 ])
        xx=range(1902,2006)
        ax.set_xticks([1902,1954,2005 ])

	index = N.arange(n_groups)
	bar_width = 0.01
	opacity = 0.6

	ax.plot(xx,ii[1:105,idx],"k-",label="T (K)")
	plt.tight_layout()

	#plt.subplots_adjust(hspace=0.5)
	plt.tick_params(
	    axis='x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='off',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='on') # labels along the bottom edge are off
#	plt.axis('off')

	ax.text(.25,.85,'{0}'.format(name[idx]),fontsize=12,
        	horizontalalignment='center',
        	transform=ax.transAxes)
#	plt.title('{0}'.format(name[idx]))
#plt.tight_layout()
	ax.tick_params(labelsize=12)

plt.savefig('soy_region_tt.png')
plt.show()
