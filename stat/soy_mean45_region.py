from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import math
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime
from statsmodels.stats.weightstats import DescrStatsW



def weighted_avg_and_std(values):

#    Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape.
        average1=N.average(values)
        variance=N.average((values-average1)**2)  # Fast and numerically precise
        return (average1,math.sqrt(variance))
def nowyield(year,couna,counb):	
    bb=year-1901

    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
    maitrop = region1.variables['soy_trop'][bb,:,:]
    maitemp = region1.variables['soy_temp'][bb,:,:]
    maitropi=region1.variables['soy_trop_irrig'][bb,:,:]
    maitempi=region1.variables['soy_temp_irrig'][bb,:,:]
    gridarea = region1.variables['area'][:,:]

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
    clmyield = clm.variables['yield'][bb,:,:]
    isam=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_soyscaleifyield_heat.nc','r')
    isamyield = isam.variables['yield'][bb,:,:]

    clmyield[N.isnan(clmyield)] = -9999
    clmyield = ma.masked_where(clmyield<=0,clmyield)
    clmyield = ma.masked_where(maizeto<=0,clmyield)
    clmyield=ma.filled(clmyield, fill_value=0.)

    isamyield[N.isnan(isamyield)] = -9999
    isamyield = ma.masked_where(isamyield<=0,isamyield)
    isamyield = ma.masked_where(maizeto<=0,isamyield)
    isamyield=ma.filled(isamyield, fill_value=0.)
    area1=NetCDFFile('/project/projectdirs/m1602/datasets4.full/surfdata_05x05.nc','r')
    lonaa=area1.variables['lon'][:]
    mask1 = area1.variables['REGION_MASK_CRU_NCEP'][:,:]
    mask,lona = shiftgrid(180.5,mask1,lonaa,start=False)

    nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
    ncvar_maize1 = nclu.variables['soy_total'][0,:,:]
    ncvar_maize,lona = shiftgrid(180.5,ncvar_maize1,lonaa,start=False)
    
    if couna>=1:
#	print couna
        if couna==4.0 or couna==5.0:
                maizeto=ma.masked_where(mask>5.0,maizeto)
                maizeto=ma.masked_where(mask<4.0,maizeto)
                maizeto=ma.filled(maizeto, fill_value=0.)
        else: 
                maizeto=ma.masked_where(mask!=couna,maizeto)
                maizeto=ma.filled(maizeto, fill_value=0.)
    maizeto = ma.masked_where(ncvar_maize<=0,maizeto)
    maizeto=ma.filled(maizeto, fill_value=0.)

#    clmyield = ma.masked_where(maizeto<=0,clmyield)
#    isamyield = ma.masked_where(maizeto<=0,isamyield)
#    maizeto = ma.masked_where(clmyield<=0,maizeto)
#    maizeto = ma.masked_where(isamyield<=0,maizeto)

#    clmyield = ma.masked_where(maizeto<=0,clmyield)
#    isamyield = ma.masked_where(maizeto<=0,isamyield)
#    maizeto = ma.masked_where(maizeto<=0,maizeto)

#    maizeto=ma.filled(maizeto, fill_value=0.)
    
    return clmyield,isamyield,maizeto


def annualyield(year,couna,counb):
    bb=year-2006
    
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
    maitrop = region1.variables['soy_trop'][bb,:,:]
    maitemp = region1.variables['soy_temp'][bb,:,:]
    maitropi=region1.variables['soy_trop_irrig'][bb,:,:]
    maitempi=region1.variables['soy_temp_irrig'][bb,:,:]
    gridarea = region1.variables['area'][:,:]

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


    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]


    clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtempf = clm1.variables['yield'][bb,:,:]

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
    clmtropfi = clm2.variables['yield'][bb,:,:]


    clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
    clmtempfi = clm3.variables['yield'][bb,:,:]


    clmtropf=N.flipud(clmtropf)
    clmtempf=N.flipud(clmtempf)
    clmtropfi=N.flipud(clmtropfi)
    clmtempfi=N.flipud(clmtempfi)

    clmtropf= ma.masked_where(maizetro<=0,clmtropf)
    clmtempf= ma.masked_where(maizetem<=0,clmtempf)
    clmtropf=ma.filled(clmtropf, fill_value=0.)
    clmtempf=ma.filled(clmtempf, fill_value=0.)

    clmtropfi= ma.masked_where(maizetro<=0,clmtropfi)
    clmtempfi= ma.masked_where(maizetem<=0,clmtempfi)
    clmtropfi=ma.filled(clmtropfi, fill_value=0.)
    clmtempfi=ma.filled(clmtempfi, fill_value=0.)

    yield_clmtf=clmtropf+clmtempf
    yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
    yield_clmtf  = ma.masked_where(maizeto<=0,yield_clmtf )
    yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

    yield_clmtfi=clmtropfi+clmtempfi
    yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
    yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)
    yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)


    clmn=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_rf_fert_0.5x0.5.nc','r')
    clmtropfn = clmn.variables['yield'][bb,:,:]


    clm1n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_rf_fert_0.5x0.5.nc','r')
    clmtempfn = clm1n.variables['yield'][bb,:,:]

    clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_irrig_fert_0.5x0.5.nc','r')
    clmtropfin = clm2n.variables['yield'][bb,:,:]


    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_irrig_fert_0.5x0.5.nc','r')
    clmtempfin = clm3n.variables['yield'][bb,:,:]


    clmtropfn=N.flipud(clmtropfn)
    clmtempfn=N.flipud(clmtempfn)
    clmtropfin=N.flipud(clmtropfin)
    clmtempfin=N.flipud(clmtempfin)

    clmtropfn= ma.masked_where(maizetro<=0,clmtropfn)
    clmtempfn= ma.masked_where(maizetem<=0,clmtempfn)
    clmtropfn=ma.filled(clmtropfn, fill_value=0.)
    clmtempfn=ma.filled(clmtempfn, fill_value=0.)

    clmtropfin= ma.masked_where(maizetro<=0,clmtropfin)
    clmtempfin= ma.masked_where(maizetem<=0,clmtempfin)
    clmtropfin=ma.filled(clmtropfin, fill_value=0.)
    clmtempfin=ma.filled(clmtempfin, fill_value=0.)

    yield_clmtfn=clmtropfn+clmtempfn
    yield_clmtfn = ma.masked_where(yield_clmtfn<=0,yield_clmtfn)
    yield_clmtfn  = ma.masked_where(maizeto<=0,yield_clmtfn )
    yield_clmtfn=ma.filled(yield_clmtfn, fill_value=0.)

    yield_clmtfin=clmtropfin+clmtempfin
    yield_clmtfin = ma.masked_where(yield_clmtfin<=0,yield_clmtfin)
    yield_clmtfin = ma.masked_where(maizeto<=0,yield_clmtfin)
    yield_clmtfin=ma.filled(yield_clmtfin, fill_value=0.)


    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytemp_rcp45_constco2_rf_nofert_0.5x0.5.nc", mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytemp_rcp45_co2_rf_nofert_0.5x0.5.nc", mode='r')

    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
    yieldf = base.variables["totalyield"][bb,:,:]
    yieldfa = base2.variables["totalyield"][bb,:,:]

    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytemp_rcp45_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytemp_rcp45_co2_irrig_fert_0.5x0.5.nc", mode='r')

    yieldfi = basei.variables["totalyield"][bb,:,:]
        
    yieldfai = base2i.variables["totalyield"][bb,:,:]
       
    baseif = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytrop_rcp45_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2if = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/soytrop_rcp45_co2_irrig_fert_0.5x0.5.nc", mode='r')

    yieldfitr = baseif.variables["totalyield"][bb,:,:]

    yieldfaitr = base2if.variables["totalyield"][bb,:,:]



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


   
    yield_new[N.isnan(yield_new)] = -9999
    yield_new = ma.masked_where(yield_new<=0,yield_new)
    yield_new = ma.masked_where(maizeto<=0,yield_new)
    yield_new=ma.filled(yield_new, fill_value=0.)

    yield_new1[N.isnan(yield_new1)] = -9999
    yield_new1 = ma.masked_where(yield_new1<=0,yield_new1)
    yield_new1 = ma.masked_where(maizeto<=0,yield_new1)
    yield_new1=ma.filled(yield_new1, fill_value=0.)

    yield_newi[N.isnan(yield_newi)] = -9999
    yield_newi = ma.masked_where(yield_newi<=0,yield_newi)
    yield_newi = ma.masked_where(maizetem<=0,yield_newi)
    yield_newi=ma.filled(yield_newi, fill_value=0.)

    yield_new1i[N.isnan(yield_new1i)] = -9999
    yield_new1i = ma.masked_where(yield_new1i<=0,yield_new1i)
    yield_new1i = ma.masked_where(maizetem<=0,yield_new1i)
    yield_new1i=ma.filled(yield_new1i, fill_value=0.)


    yield_newitr[N.isnan(yield_newitr)] = -9999
    yield_newitr = ma.masked_where(yield_newitr<=0,yield_newitr)
    yield_newitr = ma.masked_where(maizetro<=0,yield_newitr)
    yield_newitr=ma.filled(yield_newitr, fill_value=0.)

    yield_new1itr[N.isnan(yield_new1i)] = -9999
    yield_new1itr = ma.masked_where(yield_new1itr<=0,yield_new1itr)
    yield_new1itr = ma.masked_where(maizetro<=0,yield_new1itr)
    yield_new1itr=ma.filled(yield_new1itr, fill_value=0.)
  
    yield_newi=yield_newi+yield_newitr
    yield_new1i=yield_new1i+yield_new1itr

    yieldagf=0.
    yieldg=0.
    harea=0.
    a=0
    yieldgid=0.
    yieldgia=0.
    yieldgib=0.
    yieldgic=0.

    yieldgca=0.
    yieldgcb=0.
    yieldgcc=0.
    yieldgcd=0.

    yieldagfa=0.
    yieldagfb=0.
    yieldagfc=0.
    yieldagfd=0.
    yieldagfa1=0.
    yieldagfb1=0.
    yieldagfc1=0.
    yieldagfd1=0.




    yieldagfa=yield_new
    yieldagfb=yield_new1
    yieldagfc=yield_newi
    yieldagfd=yield_new1i
   
    yieldagfa1=yield_clmtf
    yieldagfb1=yield_clmtfi
    yieldagfc1=yield_clmtfn
    yieldagfd1=yield_clmtfin
    yieldagfa = ma.masked_where(yieldagfa<=0,yieldagfa)
    yieldagfb = ma.masked_where(yieldagfb<=0,yieldagfb)
    yieldagfc = ma.masked_where(yieldagfc<=0,yieldagfc)
    yieldagfd = ma.masked_where(yieldagfd<=0,yieldagfd)
    yieldagfa1 = ma.masked_where(yieldagfa1<=0,yieldagfa1)
    yieldagfb1 = ma.masked_where(yieldagfb1<=0,yieldagfb1)
    yieldagfc1 = ma.masked_where(yieldagfc1<=0,yieldagfc1)
    yieldagfd1 = ma.masked_where(yieldagfd1<=0,yieldagfd1)
   # maizeto = ma.masked_where(maizeto<=0,maizeto)

    area1=NetCDFFile('/project/projectdirs/m1602/datasets4.full/surfdata_05x05.nc','r')
    lonaa=area1.variables['lon'][:]
    mask1 = area1.variables['REGION_MASK_CRU_NCEP'][:,:]
    mask,lona = shiftgrid(180.5,mask1,lonaa,start=False)

    nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/m3yield_isam.nc','r')
    ncvar_maize1 = nclu.variables['soy_total'][0,:,:]
    ncvar_maize,lona = shiftgrid(180.5,ncvar_maize1,lonaa,start=False)
    
    if couna>=1:
        if couna==4.0 or couna==5.0:
                maizeto=ma.masked_where(mask>5.0,maizeto)
                maizeto=ma.masked_where(mask<4.0,maizeto)
		#maizeto=ma.filled(maizeto, fill_value=0.)
        else:
                maizeto=ma.masked_where(mask!=couna,maizeto)
		#maizeto=ma.filled(maizeto, fill_value=0.)
    maizeto = ma.masked_where(ncvar_maize<=0,maizeto)
    maizeto=ma.filled(maizeto, fill_value=0.)
    yieldagfa = ma.masked_where(maizeto<=0,yieldagfa)
    yieldagfb = ma.masked_where(maizeto<=0,yieldagfb)
    yieldagfc = ma.masked_where(maizeto<=0,yieldagfc)
    yieldagfd = ma.masked_where(maizeto<=0,yieldagfd)
    yieldagfa1 = ma.masked_where(maizeto<=0,yieldagfa1)
    yieldagfb1 = ma.masked_where(maizeto<=0,yieldagfb1)
    yieldagfc1 = ma.masked_where(maizeto<=0,yieldagfc1)
    yieldagfd1 = ma.masked_where(maizeto<=0,yieldagfd1)
    maizeto = ma.masked_where(maizeto<=0,maizeto)
    maizeto=ma.filled(maizeto, fill_value=0.)
    
    return yieldagfa,yieldagfa1,yieldagfb,yieldagfb1,yieldagfc,yieldagfc1,yieldagfd,yieldagfd1,maizeto
    #return harea

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



for i, name1 in enumerate(name):
        a=0
	isamhis=N.zeros((10,360,720))
        clmhis =N.zeros((10,360,720))
        frachis=N.zeros((10,360,720))

        for num in range(1996,2006):

            reu=nowyield(num,range1[i],range2[i])
            clmhis[a,:,:]=reu[0]
            isamhis[a,:,:]=reu[1]
            frachis[a,:,:]=reu[2]
            a=a+1
	isamhisy=N.average(isamhis,weights=frachis)
        clmhisy=N.average(clmhis,weights=frachis)


fisama=N.zeros((9,2))
fisamb=N.zeros((9,2))
fisamc=N.zeros((9,2))
fisamd=N.zeros((9,2))
for i, name1 in enumerate(name):
        a=0
	clmya=N.zeros((x,360,720))
	clmyb=N.zeros((x,360,720))
	clmyc=N.zeros((x,360,720))
	isamyd=N.zeros((x,360,720))
	clmyd=N.zeros((x,360,720))
	isamya=N.zeros((x,360,720))
	isamyb=N.zeros((x,360,720))
	isamyc=N.zeros((x,360,720))
	frac=N.zeros((x,360,720))

	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])
	#    print i, num,a1,a2,range1[i]
	    isamya[a,:,:]=reu[0]
	    isamyb[a,:,:]=reu[2]
	    clmya[a,:,:]=reu[1]
	    clmyb[a,:,:]=reu[3]
	    isamyc[a,:,:]=reu[4]
	    clmyc[a,:,:]=reu[5]
            isamyd[a,:,:]=reu[6]
            clmyd[a,:,:]=reu[7]
            frac[a,:,:]=reu[8]
	 #   print a
	    a=a+1
        c1=0
        c2=10
        frac50=frac[c1:c2,:,:]
        isama50=N.average(isamya,weights=frac)
 	isamb50=N.average(isamyb,weights=frac)
        isamc50=N.average(isamyc,weights=frac)
        isamd50=N.average(isamyd,weights=frac)
        clma50=N.average(clmya,weights=frac)
        clmb50=N.average(clmyb,weights=frac)
        clmc50=N.average(clmyc,weights=frac)
        clmd50=N.average(clmyd,weights=frac)
	

	fisama[i,0]=(isama50-isamhisy)/isamhisy*100.
        fisamb[i,0]=(isamb50-isamhisy)/isamhisy*100.
        fisamc[i,0]=(isamc50-isamhisy)/isamhisy*100.
        fisamd[i,0]=(isamd50-isamhisy)/isamhisy*100.
        fisama[i,1]=(clma50-clmhisy)/clmhisy*100.
        fisamb[i,1]=(clmb50-clmhisy)/clmhisy*100.
        fisamc[i,1]=(clmc50-clmhisy)/clmhisy*100.
        fisamd[i,1]=(clmd50-clmhisy)/clmhisy*100.

print fisama
fig = plt.figure(figsize=(4,3.5))
plt.rc('font', weight='bold')
for idx in xrange(9):
        ax = fig.add_subplot(3, 3, idx+1)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

      
	n_groups = 2
	index = N.arange(n_groups)
        plt.ylim(-50, 100)
	bar_width = 0.3
	opacity = 0.6
	for x in range(0,2):
	
        	rects0 = plt.bar(index[x], fisama[idx,x], bar_width,
                	 alpha=0.8,
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
	                	 alpha=0.8,
	                 	color='g',
	                 	label='CO2+CLIMATE+N')
		else:
			rects1 = plt.bar(index[x], fisamc[idx,x]-fisamb[idx,x], bar_width,
                                 alpha=0.8,
                                color='g',
                                bottom=fisamb[idx,x],
                                label='CO2+CLIMATE+N')
                if fisamc[idx,x]<0.0:

	                rects2 = plt.bar(index[x] , fisamd[idx,x], bar_width,
        	                alpha=0.8,
                	        color='r',
                        	label='CO2+CLIMATE+N+I')

		else:
	        	rects2 = plt.bar(index[x] , fisamd[idx,x]-fisamc[idx,x], bar_width,
        	         	alpha=0.8,
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
	        ax.text(.55,.85,'{0}'.format(name[idx]),
        	        horizontalalignment='center',
                	transform=ax.transAxes)
	        ax.tick_params(labelsize=10)


plt.savefig('soy2090_45_region.png')
plt.show()



