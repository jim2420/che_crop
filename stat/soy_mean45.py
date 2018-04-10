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


country=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
#print iizumi
coun = country.variables['MASK_Country'][:,:]
#area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
#gridarea = area.variables['cell_area'][:,:]
#gridlon = area.variables['lon'][:]

#gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)

def weighted_avg_and_std(values):

#    Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape.
        average1=N.average(values)
        variance=N.average((values-average1)**2)  # Fast and numerically precise
        return (average1,math.sqrt(variance))


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


    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytemp_rcp45_constco2_rf_nofert_0.5x0.5.nc", mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytemp_rcp45_co2_rf_nofert_0.5x0.5.nc", mode='r')

    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
    yieldf = base.variables["totalyield"][bb,:,:]
    yieldfa = base2.variables["totalyield"][bb,:,:]

    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytemp_rcp45_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytemp_rcp45_co2_irrig_fert_0.5x0.5.nc", mode='r')

    yieldfi = basei.variables["totalyield"][bb,:,:]
        
    yieldfai = base2i.variables["totalyield"][bb,:,:]
       
    baseif = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytrop_rcp45_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2if = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/soytrop_rcp45_co2_irrig_fert_0.5x0.5.nc", mode='r')

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
    maizeto = ma.masked_where(maizeto<=0,maizeto)

    return yieldagfa,yieldagfa1,yieldagfb,yieldagfb1,yieldagfc,yieldagfc1,yieldagfd,yieldagfd1,maizeto
    #return harea

a1=2010
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
name=["Global"]
range1=[10100]
range2=[50700]

#name=["USA"]
#range1=[11501]
#range2=[11550]





for i, name1 in enumerate(name):
        a=0
	clmya= N.zeros((x,360,720))
	clmyb=N.zeros((x,360,720))
	clmyc=N.zeros((x,360,720))
	isamyd= N.zeros((x,360,720))
	clmyd= N.zeros((x,360,720))
	isamya= N.zeros((x,360,720))
	isamyb= N.zeros((x,360,720))
	isamyc= N.zeros((x,360,720))
	frac=N.zeros((x,360,720))

	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])
	    isamya[a,:,:]=reu[0]
	    isamyb[a,:,:]=reu[2]
	    clmya[a,:,:]=reu[1]
	    clmyb[a,:,:]=reu[3]
	    isamyc[a,:,:]=reu[4]
	    clmyc[a,:,:]=reu[5]
            isamyd[a,:,:]=reu[6]
            clmyd[a,:,:]=reu[7]
            frac[a,:,:]=reu[8]
	    a=a+1
        d1=2010-2010
        d2=2020-2010
        c1=2050-2010
        c2=2060-2010
        frac10=frac[d1:d2,:,:]
        frac10a=N.average(frac10,axis=0)
        frac50=frac[c1:c2,:,:]
        frac50a=N.average(frac50,axis=0)
        isama50=N.average(isamya[c1:c2,:,:],weights=frac50)
 	isamb50=N.average(isamyb[c1:c2,:,:],weights=frac50)
        isamc50=N.average(isamyc[c1:c2,:,:],weights=frac50)
        isamd50=N.average(isamyd[c1:c2,:,:],weights=frac50)
        isama10=N.average(isamya[d1:d2,:,:],weights=frac10)
        isamb10=N.average(isamyb[d1:d2,:,:],weights=frac10)
        isamc10=N.average(isamyc[d1:d2,:,:],weights=frac10)
        isamd10=N.average(isamyd[d1:d2,:,:],weights=frac10)
        clma50=N.average(clmya[c1:c2,:,:],weights=frac50)
        clmb50=N.average(clmyb[c1:c2,:,:],weights=frac50)
        clmc50=N.average(clmyc[c1:c2,:,:],weights=frac50)
        clmd50=N.average(clmyd[c1:c2,:,:],weights=frac50)
        clma10=N.average(clmya[d1:d2,:,:],weights=frac10)
        clmb10=N.average(clmyb[d1:d2,:,:],weights=frac10)
        clmc10=N.average(clmyc[d1:d2,:,:],weights=frac10)
        clmd10=N.average(clmyd[d1:d2,:,:],weights=frac10)
	

	fisama=(isama50-isama10)/isama10*100.
        fisamb=(isamb50-isamb10)/isamb10*100.
        fisamc=(isamc50-isamc10)/isamc10*100.
        fisamd=(isamd50-isamd10)/isamd10*100.
        fclma=(clma50-clma10)/clma10*100.
        fclmb=(clmb50-clmb10)/clmb10*100.
        fclmc=(clmc50-clmc10)/clmc10*100.
        fclmd=(clmd50-clmd10)/clmd10*100.

        fisama=N.append(fisama,fclma)
        fisamb=N.append(fisamb,fclmb)
        fisamc=N.append(fisamc,fclmc)
        fisamd=N.append(fisamd,fclmd)


        clmb50new=N.average(clmyb[c1:c2,:,:],weights=frac50,axis=(1,2))
        clmb10new=N.average(clmyb[d1:d2,:,:],weights=frac10,axis=(1,2))
        clmbnew=(clmb50new-clmb10new)/clmb10new*100.
        clmbf=weighted_avg_and_std(clmbnew)

        clmc50new=N.average(clmyc[c1:c2,:,:],weights=frac50,axis=(1,2))
        clmc10new=N.average(clmyc[d1:d2,:,:],weights=frac10,axis=(1,2))
        clmcnew=(clmc50new-clmc10new)/clmc10new*100.
        clmcf=weighted_avg_and_std(clmcnew)

        clma50new=N.average(clmya[c1:c2,:,:],weights=frac50,axis=(1,2))
        clma10new=N.average(clmya[d1:d2,:,:],weights=frac10,axis=(1,2))
        clmanew=(clma50new-clma10new)/clma10new*100.
	print clmanew
        clmaf=weighted_avg_and_std(clmanew)

        clmd50new=N.average(clmyd[c1:c2,:,:],weights=frac50,axis=(1,2))
        clmd10new=N.average(clmyd[d1:d2,:,:],weights=frac10,axis=(1,2))
        clmdnew=(clmd50new-clmd10new)/clmd10new*100.
        clmdf=weighted_avg_and_std(clmdnew)

        isama50new=N.average(isamya[c1:c2,:,:],weights=frac50,axis=(1,2))
        isama10new=N.average(isamya[d1:d2,:,:],weights=frac10,axis=(1,2))
        isamanew=(isama50new-isama10new)/isama10new*100.
        isamaf=weighted_avg_and_std(isamanew)

        isamb50new=N.average(isamyb[c1:c2,:,:],weights=frac50,axis=(1,2))
        isamb10new=N.average(isamyb[d1:d2,:,:],weights=frac10,axis=(1,2))
        isambnew=(isamb50new-isamb10new)/isamb10new*100.
        isambf=weighted_avg_and_std(isambnew)

        isamc50new=N.average(isamyc[c1:c2,:,:],weights=frac50,axis=(1,2))
        isamc10new=N.average(isamyc[d1:d2,:,:],weights=frac10,axis=(1,2))
        isamcnew=(isamc50new-isamc10new)/isamc10new*100.
        isamcf=weighted_avg_and_std(isamcnew)

        isamd50new=N.average(isamyd[c1:c2,:,:],weights=frac50,axis=(1,2))
        isamd10new=N.average(isamyd[d1:d2,:,:],weights=frac10,axis=(1,2))
        isamdnew=(isamd50new-isamd10new)/isamd10new*100.
        isamdf=weighted_avg_and_std(isamdnew)

	isamaf1=N.append(isamaf[0],clmaf[0])
        isambf1=N.append(isambf[0],clmbf[0])
        isamcf1=N.append(isamcf[0],clmcf[0])
        isamdf1=N.append(isamdf[0],clmdf[0])

        isamaf2=N.append(isamaf[1],clmaf[1])
        isambf2=N.append(isambf[1],clmbf[1])
        isamcf2=N.append(isamcf[1],clmcf[1])
        isamdf2=N.append(isamdf[1],clmdf[1])


	fig = plt.figure(figsize=(5,5))

	n_groups = 2
	ax = fig.add_subplot(111)
	index = N.arange(n_groups)
	bar_width = 0.2
	opacity = 0.8
	rects0 = plt.bar(index, fisama, bar_width,
                 alpha=opacity,
                 color='k',
                 label='CLIMATE')
	rects4 = plt.bar(index+bar_width, fisamb, bar_width,
                 alpha=opacity,
                 color='c',
                 label='CO2+CLIMATE')
	rects1 = plt.bar(index+bar_width*2, fisamc, bar_width,
                 alpha=opacity,
                 color='r',
                 label='CO2+CLIMATE+N')
	rects2 = plt.bar(index + bar_width*3, fisamd, bar_width,
                 alpha=opacity,
                 color='g',
                 label='CO2+CLIMATE+N+I')

	plt.ylim(-20,85)

	plt.ylabel('percentage of yield changes (%)',fontsize=18)
#	plt.title('2050s-2010s maize',fontsize=18)
	plt.xticks(index + bar_width+0.2, ('ISAM','CLM'))
#	leg=plt.legend(loc=4)
#	leg.get_frame().set_alpha(0.5)
	plt.tick_params(axis='both',labelsize=18)

	plt.tight_layout()

	plt.savefig('soy2050_45_test.png')
	plt.show()



