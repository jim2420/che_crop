from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
lonisam1=region1.variables['lon'][:]
maitrop = region1.variables['soy_trop'][95:105,:,:]
maitemp = region1.variables['soy_temp'][95:105,:,:]
maitropi=region1.variables['soy_trop_irrig'][95:105,:,:]
maitempi=region1.variables['soy_temp_irrig'][95:105,:,:]
maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)

maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)

maizetro=maitrop+maitropi
maizetem=maitemp+maitempi
maizetor=maitrop+maitemp
maizetoi=maitropi+maitempi
maizeto = maitrop+maitemp+maitropi+maitempi

clm=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/clm/clm45his_soyscaleifyield.nc','r')
clmyield = clm.variables['yield'][95:105,:,:]
isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_soyscaleifyield_heat.nc','r')
isamyield = isam.variables['yield'][95:105,:,:]
clmyield[N.isnan(clmyield)] = 0
clmyield = ma.masked_where(clmyield<=0,clmyield)
clmyield = ma.masked_where(maizeto<=0,clmyield)
clmyield=ma.filled(clmyield, fill_value=0.)

isamyield[N.isnan(isamyield)] = 0
isamyield = ma.masked_where(isamyield<=0,isamyield)
isamyield = ma.masked_where(maizeto<=0,isamyield)
isamyield=ma.filled(isamyield, fill_value=0.)

    
region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
maitropf = region1.variables['soy_trop'][84:94,:,:]
maitempf = region1.variables['soy_temp'][84:94,:,:]
maitropif=region1.variables['soy_trop_irrig'][84:94,:,:]
maitempif=region1.variables['soy_temp_irrig'][84:94,:,:]

maitropf=ma.masked_where(maitropf<=0,maitropf)
maitropf=ma.filled(maitropf, fill_value=0.)
maitempf=ma.masked_where(maitempf<=0,maitempf)
maitempf=ma.filled(maitempf, fill_value=0.)

maitropif=ma.masked_where(maitropif<=0,maitropif)
maitropif=ma.filled(maitropif, fill_value=0.)
maitempif=ma.masked_where(maitempif<=0,maitempif)
maitempif=ma.filled(maitempif, fill_value=0.)

maizetrof=maitropf+maitropif
maizetemf=maitempf+maitempif
maizetorf=maitropf+maitempf
maizetoif=maitropif+maitempif
maizetof = maitropf+maitempf+maitropif+maitempif
clmtropf=N.zeros((10,360,720))
clmtempf=N.zeros((10,360,720))
clmtropfi=N.zeros((10,360,720))
clmtempfi=N.zeros((10,360,720))

for j in range(0,10):
	clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
	aa= N.flipud(clm.variables['yield'][84+j,:,:])
	clmtropf[j,:,:] = aa
	
	clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
	bb = N.flipud(clm1.variables['yield'][84+j,:,:])
	clmtempf[j,:,:] = bb

	clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
	cc = N.flipud(clm2.variables['yield'][84+j,:,:])
	clmtropfi[j,:,:] = cc

	clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
	dd = N.flipud(clm3.variables['yield'][84+j,:,:])
	clmtempfi[j,:,:] = dd

#clmtropf=N.flipud(clmtropf)
#clmtempf=N.flipud(clmtempf)
#clmtropfi=N.flipud(clmtropfi)
#clmtempfi=N.flipud(clmtempfi)

clmtropf= ma.masked_where(clmtropf>=10000,clmtropf)
clmtropfi= ma.masked_where(clmtropfi>=10000,clmtropfi)
clmtempf= ma.masked_where(clmtempf>=10000,clmtempf)
clmtempfi= ma.masked_where(clmtempfi>=10000,clmtempfi)
clmtropf[N.isnan(clmtropf)] = 0
clmtropfi[N.isnan(clmtropfi)] = 0
clmtropf[N.isinf(clmtropf)] = 0
clmtropfi[N.isinf(clmtropfi)] = 0

clmtempf[N.isnan(clmtempf)] = 0
clmtempfi[N.isnan(clmtempfi)] = 0
clmtempf[N.isinf(clmtempf)] = 0
clmtempfi[N.isinf(clmtempfi)] = 0


clmtropf= ma.masked_where(maizetrof<=0,clmtropf)
clmtempf= ma.masked_where(maizetemf<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfi= ma.masked_where(maizetrof<=0,clmtropfi)
clmtempfi= ma.masked_where(maizetemf<=0,clmtempfi)
clmtropfi=ma.filled(clmtropfi, fill_value=0.)
clmtempfi=ma.filled(clmtempfi, fill_value=0.)

yield_clmtf=clmtropf+clmtempf
yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
yield_clmtf  = ma.masked_where(maizetof<=0,yield_clmtf )
yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

yield_clmtfi=clmtropfi+clmtempfi
yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
yield_clmtfi = ma.masked_where(maizetof<=0,yield_clmtfi)
yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)

clmtropfn=N.zeros((10,360,720))
clmtempfn=N.zeros((10,360,720))
clmtropfin=N.zeros((10,360,720))
clmtempfin=N.zeros((10,360,720))

for j in range(0,10):
        clmn=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_rf_fert_0.5x0.5.nc','r')
        aa= N.flipud(clmn.variables['yield'][84+j,:,:])
        clmtropfn[j,:,:] = aa

        clm1n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_rf_fert_0.5x0.5.nc','r')
        bb = N.flipud(clm1n.variables['yield'][84+j,:,:])
        clmtempfn[j,:,:] = bb

        clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
        cc = N.flipud(clm2n.variables['yield'][84+j,:,:])
        clmtropfin[j,:,:] = cc

        clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
        dd = N.flipud(clm3n.variables['yield'][84+j,:,:])
        clmtempfin[j,:,:] = dd





#clmtropfn=N.flipud(clmtropfn)
#clmtempfn=N.flipud(clmtempfn)
#clmtropfin=N.flipud(clmtropfin)
#clmtempfin=N.flipud(clmtempfin)


clmtropfn[N.isnan(clmtropfn)] = 0
clmtropfin[N.isnan(clmtropfin)] = 0
clmtropfn[N.isinf(clmtropfn)] = 0
clmtropfin[N.isinf(clmtropfin)] = 0

clmtempfn[N.isnan(clmtempfn)] = 0
clmtempfin[N.isnan(clmtempfin)] = 0
clmtempfn[N.isinf(clmtempfn)] = 0
clmtempfin[N.isinf(clmtempfin)] = 0


clmtropfn= ma.masked_where(maizetrof<=0,clmtropfn)
clmtempfn= ma.masked_where(maizetemf<=0,clmtempfn)
clmtropfn=ma.filled(clmtropfn, fill_value=0.)
clmtempfn=ma.filled(clmtempfn, fill_value=0.)

clmtropfin= ma.masked_where(maizetrof<=0,clmtropfin)
clmtempfin= ma.masked_where(maizetemf<=0,clmtempfin)
clmtropfin=ma.filled(clmtropfin, fill_value=0.)
clmtempfin=ma.filled(clmtempfin, fill_value=0.)

yield_clmtfn=clmtropfn+clmtempfn
yield_clmtfn = ma.masked_where(yield_clmtfn<=0,yield_clmtfn)
yield_clmtfn  = ma.masked_where(maizetof<=0,yield_clmtfn )
yield_clmtfn=ma.filled(yield_clmtfn, fill_value=0.)

yield_clmtfin=clmtropfin+clmtempfin
yield_clmtfin = ma.masked_where(yield_clmtfin<=0,yield_clmtfin)
yield_clmtfin = ma.masked_where(maizetof<=0,yield_clmtfin)
yield_clmtfin=ma.filled(yield_clmtfin, fill_value=0.)


base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytemp_rcp85_constco2_rf_nofert_0.5x0.5.nc", mode='r')
base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytemp_rcp85_co2_rf_nofert_0.5x0.5.nc", mode='r')

lona1 = base.variables["lon"][:]
lata1 = base.variables["lat"][:]
yieldf = base.variables["totalyield"][84:94,:,:]
yieldfa = base2.variables["totalyield"][84:94,:,:]

basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytemp_rcp85_co2_rf_fert_0.5x0.5.nc", mode='r')
base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytemp_rcp85_co2_irrig_fert_0.5x0.5.nc", mode='r')

yieldfi = basei.variables["totalyield"][84:94,:,:]
yieldfai = base2i.variables["totalyield"][84:94,:,:]
       
baseif = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytrop_rcp85_co2_rf_fert_0.5x0.5.nc", mode='r')
base2if = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/soytrop_rcp85_co2_irrig_fert_0.5x0.5.nc", mode='r')

yieldfitr = baseif.variables["totalyield"][84:94,:,:]

yieldfaitr = base2if.variables["totalyield"][84:94,:,:]

yielda=yieldf
yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yieldb=yieldfa
yield_new1,lona11 = shiftgrid(180.5,yieldb,lona1,start=False)

yieldai=yieldfi
yield_newi,lona11 = shiftgrid(180.5,yieldai,lona1,start=False)
yieldbi=yieldfai

yield_new1i,lona11 = shiftgrid(180.5,yieldbi,lona1,start=False)
yieldaitr=yieldfitr

yield_newitr,lona11 = shiftgrid(180.5,yieldaitr,lona1,start=False)
yieldbitr=yieldfaitr

yield_new1itr,lona11 = shiftgrid(180.5,yieldbitr,lona1,start=False)

yield_new[N.isnan(yield_new)] = 0
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(maizetof<=0,yield_new)
yield_new=ma.filled(yield_new, fill_value=0.)

yield_new1[N.isnan(yield_new1)] = 0
yield_new1 = ma.masked_where(yield_new1<=0,yield_new1)
yield_new1 = ma.masked_where(maizetof<=0,yield_new1)
yield_new1=ma.filled(yield_new1, fill_value=0.)

yield_newi[N.isnan(yield_newi)] = 0
yield_newi = ma.masked_where(yield_newi<=0,yield_newi)
yield_newi = ma.masked_where(maizetemf<=0,yield_newi)
yield_newi=ma.filled(yield_newi, fill_value=0.)

yield_new1i[N.isnan(yield_new1i)] = 0
yield_new1i = ma.masked_where(yield_new1i<=0,yield_new1i)
yield_new1i = ma.masked_where(maizetemf<=0,yield_new1i)
yield_new1i=ma.filled(yield_new1i, fill_value=0.)


yield_newitr[N.isnan(yield_newitr)] = 0
yield_newitr = ma.masked_where(yield_newitr<=0,yield_newitr)
yield_newitr = ma.masked_where(maizetrof<=0,yield_newitr)
yield_newitr=ma.filled(yield_newitr, fill_value=0.)

yield_new1itr[N.isnan(yield_new1i)] =0
yield_new1itr = ma.masked_where(yield_new1itr<=0,yield_new1itr)
yield_new1itr = ma.masked_where(maizetrof<=0,yield_new1itr)
yield_new1itr=ma.filled(yield_new1itr, fill_value=0.)
  
yield_newi=yield_newi+yield_newitr
yield_new1i=yield_new1i+yield_new1itr


yieldagfa=yield_new
yieldagfb=yield_new1
yieldagfc=yield_newi
yieldagfd=yield_new1i
   
yieldagfa1=yield_clmtf
yieldagfb1=yield_clmtfi
yieldagfc1=yield_clmtfn
yieldagfd1=yield_clmtfin
#yieldagfa = ma.masked_where(yieldagfa<=0,yieldagfa)
#yieldagfb = ma.masked_where(yieldagfb<=0,yieldagfb)
#yieldagfc = ma.masked_where(yieldagfc<=0,yieldagfc)
#yieldagfd = ma.masked_where(yieldagfd<=0,yieldagfd)
#yieldagfa1 = ma.masked_where(yieldagfa1<=0,yieldagfa1)
#yieldagfb1 = ma.masked_where(yieldagfb1<=0,yieldagfb1)
#yieldagfc1 = ma.masked_where(yieldagfc1<=0,yieldagfc1)
#yieldagfd1 = ma.masked_where(yieldagfd1<=0,yieldagfd1)


a1=2090
a2=2100
x=a2-a1
clmya=N.zeros((x,360,720))
clmyb=N.zeros((x,360,720))
clmyc=N.zeros((x,360,720))
isamyd=N.zeros((x,360,720))
clmyd=N.zeros((x,360,720))
isamya=N.zeros((x,360,720))
isamyb=N.zeros((x,360,720))
isamyc=N.zeros((x,360,720))
frac=N.zeros((x,360,720))

frachis=N.zeros((10,360,720))
clmhis=N.zeros((10,360,720))
isamhis=N.zeros((10,360,720))

name=["Global","NA","SA","EU","Africa","Africa","USSR","China","SSEA"]
range1=[0,1,2,3,4,5,6,7,8]
range2=[0,1,2,3,4,5,6,7,8]

fisama=N.zeros((9,2))
fisamb=N.zeros((9,2))
fisamc=N.zeros((9,2))
fisamd=N.zeros((9,2))

area1=NetCDFFile('/project/projectdirs/m1602/datasets4.full/surfdata_05x05.nc','r')
lonaa=area1.variables['lon'][:]
mask1 = area1.variables['REGION_MASK_CRU_NCEP'][:,:]
mask,lona = shiftgrid(180.5,mask1,lonaa,start=False)

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
ncvar_maize1 = nclu.variables['soyy'][0,:,:]
ncvar_maize,lona = shiftgrid(180.5,ncvar_maize1,lonaa,start=False)
area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridarea,lona = shiftgrid(180.5,gridarea,lonaa,start=False)
gridarea2=N.zeros((10,360,720))
mmarea=N.zeros((10,360,720))
rmask=N.zeros((10,360,720))
maizeto2=N.zeros((10,360,720))
for i in range(0,10):
	for x in range(0,360):
		for y in range(0,720):
			mmarea[i,x,y]=ncvar_maize[x,y]
			rmask[i,x,y]=mask[x,y]
                        maizeto2[i,x,y]=maizeto[i,x,y]
			gridarea2[i,x,y]=gridarea[x,y]

isamyield= ma.masked_where(mmarea<=0.,isamyield)
clmyield= ma.masked_where(mmarea<=0.,clmyield)
yieldagfa= ma.masked_where(mmarea<=0.,yieldagfa)
yieldagfb= ma.masked_where(mmarea<=0.,yieldagfb)
yieldagfc= ma.masked_where(mmarea<=0.,yieldagfc)
yieldagfd= ma.masked_where(mmarea<=0.,yieldagfd)
yieldagfa1= ma.masked_where(mmarea<=0.,yieldagfa1)
yieldagfb1= ma.masked_where(mmarea<=0.,yieldagfb1)
yieldagfc1= ma.masked_where(mmarea<=0.,yieldagfc1)
yieldagfd1= ma.masked_where(mmarea<=0.,yieldagfd1)

for i in range(1,9):
	maizetohis=maizeto2
	maizetofu=maizetof	
	if i==4 or i==5:
		maizetohis=ma.masked_where(rmask>5.0,maizetohis)
                maizetohis=ma.masked_where(rmask<4.0,maizetohis)
                maizetofu=ma.masked_where(rmask>5.0,maizetofu)
                maizetofu=ma.masked_where(rmask<4.0,maizetofu)

        else:
                maizetohis=ma.masked_where(rmask!=i,maizetohis)
                maizetofu=ma.masked_where(rmask!=i,maizetofu)

	#print maizetohis1[9,211,166]
        isamhisy=N.average(N.average(isamyield,weights=maizetohis*gridarea2,axis=(1,2)))
        clmhisy=N.average(N.average(clmyield,weights=maizetohis*gridarea2,axis=(1,2)))

        isama50=N.average(N.average(yieldagfa,weights=maizetofu*gridarea2,axis=(1,2)))
        clma50=N.average(N.average(yieldagfa1,weights=maizetofu*gridarea2,axis=(1,2)))
        isamb50=N.average(N.average(yieldagfb,weights=maizetofu*gridarea2,axis=(1,2)))
        clmb50=N.average(N.average(yieldagfb1,weights=maizetofu*gridarea2,axis=(1,2)))
        isamc50=N.average(N.average(yieldagfc,weights=maizetofu*gridarea2,axis=(1,2)))
        clmc50=N.average(N.average(yieldagfc1,weights=maizetofu*gridarea2,axis=(1,2)))
        isamd50=N.average(N.average(yieldagfd,weights=maizetofu*gridarea2,axis=(1,2)))
        clmd50=N.average(N.average(yieldagfd1,weights=maizetofu*gridarea2,axis=(1,2)))
	#print isama50,isamhisy,clma50,clmhisy
	fisama[i,0]=(isama50-isamhisy)/isamhisy*100.
        fisamb[i,0]=(isamb50-isamhisy)/isamhisy*100.
        fisamc[i,0]=(isamc50-isamhisy)/isamhisy*100.
        fisamd[i,0]=(isamd50-isamhisy)/isamhisy*100.
        fisama[i,1]=(clma50-clmhisy)/clmhisy*100.
        fisamb[i,1]=(clmb50-clmhisy)/clmhisy*100.
        fisamc[i,1]=(clmc50-clmhisy)/clmhisy*100.
        fisamd[i,1]=(clmd50-clmhisy)/clmhisy*100.


isamhisy=N.average(isamyield,weights=maizeto2*gridarea2)
clmhisy=N.average(clmyield,weights=maizeto2*gridarea2)

isama50=N.average(yieldagfa,weights=maizetof*gridarea2)
clma50=N.average(yieldagfa1,weights=maizetof*gridarea2)
isamb50=N.average(yieldagfb,weights=maizetof*gridarea2)
clmb50=N.average(yieldagfb1,weights=maizetof*gridarea2)
isamc50=N.average(yieldagfc,weights=maizetof*gridarea2)
clmc50=N.average(yieldagfc1,weights=maizetof*gridarea2)
isamd50=N.average(yieldagfd,weights=maizetof*gridarea2)
clmd50=N.average(yieldagfd1,weights=maizetof*gridarea2)

fisama[0,0]=(isama50-isamhisy)/isamhisy*100.
fisamb[0,0]=(isamb50-isamhisy)/isamhisy*100.
fisamc[0,0]=(isamc50-isamhisy)/isamhisy*100.
fisamd[0,0]=(isamd50-isamhisy)/isamhisy*100.
fisama[0,1]=(clma50-clmhisy)/clmhisy*100.
fisamb[0,1]=(clmb50-clmhisy)/clmhisy*100.
fisamc[0,1]=(clmc50-clmhisy)/clmhisy*100.
fisamd[0,1]=(clmd50-clmhisy)/clmhisy*100.

print fisama,fisamb, fisamc, fisamd
fig = plt.figure(figsize=(8,8))
plt.rc('font', weight='bold')
for idx in xrange(9):
        ax = fig.add_subplot(3, 3, idx+1)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

      
	n_groups = 2
	index = N.arange(n_groups)
        plt.ylim(-60, 200)
	bar_width = 0.6
	opacity = 1.0
	for x in range(0,2):
	
        	rects0 = plt.bar(index[x], fisama[idx,x], bar_width,
                	 alpha=1.0,
                   	color='b',
                 	label='CLIMATE')
		if fisama[idx,x]<0.0:
		        	rects4 = plt.bar(index[x], fisamb[idx,x], bar_width,
                		 	alpha=opacity,
		                 	color='k',
		                 	label='CO2+CLIMATE')
		else:	
				rects4 = plt.bar(index[x], fisamb[idx,x]-fisama[idx,x], bar_width,
                                	 alpha=opacity,
                                 	color='k',
				 	bottom=fisama[idx,x],
                                 	label='CO2+CLIMATE')
                if fisamb[idx,x]<0.0:
	        	rects1 = plt.bar(index[x], fisamc[idx,x], bar_width,
	                	 alpha=1.0,
	                 	color='g',
	                 	label='CO2+CLIMATE+N')
		else:
			rects1 = plt.bar(index[x], fisamc[idx,x]-fisamb[idx,x], bar_width,
                                 alpha=1.0,
                                color='g',
                                bottom=fisamb[idx,x],
                                label='CO2+CLIMATE+N')
                if fisamc[idx,x]<0.0:

	                rects2 = plt.bar(index[x] , fisamd[idx,x], bar_width,
        	                alpha=1.0,
                	        color='r',
                        	label='CO2+CLIMATE+N+I')

		else:
	        	rects2 = plt.bar(index[x] , fisamd[idx,x]-fisamc[idx,x], bar_width,
        	         	alpha=1.0,
                	 	color='r',
	         		bottom=fisamc[idx,x],
                 		label='CO2+CLIMATE+N+I')

	        plt.tight_layout()
	        plt.tick_params(
        	    axis='x',          # changes apply to the x-axis
      	      	    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off') # labels along the bottom edge are off

	        #plt.ylim(-35,90)
	        #plt.yticks(N.arange(-35,90,10))

	#plt.xticks(index + bar_width+0.2, ('ISAM','CLM'))
	        ax.text(.25,.85,'{0}'.format(name[idx]),
        	        horizontalalignment='center',
                	transform=ax.transAxes,fontsize=16)
	        ax.tick_params(labelsize=16)


plt.savefig('soy2090_85_region.png')
plt.show()



