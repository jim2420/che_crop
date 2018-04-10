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



def annualyield(year,couna,counb):
    bb=year-2006
    
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
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


    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]


    clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtempf = clm1.variables['yield'][bb,:,:]

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
    clmtropfi = clm2.variables['yield'][bb,:,:]


    clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
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


    clmn=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_rf_fert_0.5x0.5.nc','r')
    clmtropfn = clmn.variables['yield'][bb,:,:]


    clm1n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_rf_fert_0.5x0.5.nc','r')
    clmtempfn = clm1n.variables['yield'][bb,:,:]

    clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytrop_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
    clmtropfin = clm2n.variables['yield'][bb,:,:]


    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/soytemp_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
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


    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytemp_rcp85_constco2_rf_nofert_0.5x0.5.nc", mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytemp_rcp85_co2_rf_nofert_0.5x0.5.nc", mode='r')

    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
    yieldf = base.variables["totalyield"][bb,:,:]
    yieldfa = base2.variables["totalyield"][bb,:,:]

    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytemp_rcp85_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytemp_rcp85_co2_irrig_fert_0.5x0.5.nc", mode='r')

    yieldfi = basei.variables["totalyield"][bb,:,:]
        
    yieldfai = base2i.variables["totalyield"][bb,:,:]
       
    baseif = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytrop_rcp85_co2_rf_fert_0.5x0.5.nc", mode='r')
    base2if = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/soytrop_rcp85_co2_irrig_fert_0.5x0.5.nc", mode='r')

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



    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                yieldg=0+yieldg
                harea=maizeto[xx,yy]*gridarea[xx,yy]+ harea

                yieldgia=(yield_new[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgia
                yieldgca=(yield_clmtf[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgca
                yieldgib=(yield_new1[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgib
                yieldgcc=(yield_clmtfn[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcc

                yieldgic=(yield_newi[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgic
                yieldgcb=(yield_clmtfi[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcb
                yieldgid=(yield_new1i[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgid
                yieldgcd=(yield_clmtfin[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcd

                a=a+1

    yieldagf=yieldg/harea

    yieldagfa=yieldgia/harea
    yieldagfb=yieldgib/harea
    yieldagfc=yieldgic/harea
    yieldagfd=yieldgid/harea
   
    yieldagfa1=yieldgca/harea
    yieldagfb1=yieldgcb/harea
    yieldagfc1=yieldgcc/harea
    yieldagfd1=yieldgcd/harea

    return "Maize ","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-a yield",yieldagfa,"t/ha","CLM-a yield",yieldagfa1,"t/ha","ISAM-b yield",yieldagfb,"t/ha","clm-a yield",yieldagfb1,"tonnes","ISAM-c yield",yieldagfc,"t/ha","clm-c yield",yieldagfc1,"t/ha","ISAM-d yield",yieldagfd,"clm-a yield",yieldagfd1
    #return harea

a1=2010
a2=2061
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmya= N.zeros(x)
clmyb=N.zeros(x)
clmyc=N.zeros(x)
isamyd= N.zeros(x)
zumip= N.zeros(x)
clmyd= N.zeros(x)
isamya= N.zeros(x)
isamyb= N.zeros(x)
isamyc= N.zeros(x)

name=["Global"]
range1=[10100]
range2=[50700]

#name=["USA"]
#range1=[11501]
#range2=[11550]





for i, name1 in enumerate(name):
        a=0
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
	clmya= N.zeros(x)
	isamya= N.zeros(x)
	zumip= N.zeros(x)
	clmyb= N.zeros(x)
	isamyb= N.zeros(x)
	isamyc= N.zeros(x)
	isamyd= N.zeros(x)
        clmyc=N.zeros(x)
        clmyd=N.zeros(x)
	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])

	    toarea[a]=reu[2]
	    zumiy[a]=reu[5]
	    isamya[a]=reu[11]
	    isamyb[a]=reu[17]
	    zumip[a]=reu[8]
	    clmya[a]=reu[14]
	    clmyb[a]=reu[20]
	    isamyc[a]=reu[23]
	    clmyc[a]=reu[26]
            isamyd[a]=reu[29]
            clmyd[a]=reu[31]

	    a=a+1
        d1=2010-2010
        d2=2020-2010
        c1=2050-2010
        c2=2060-2010

        isama50=N.average(isamya[c1:c2])
        isamb50=N.average(isamyb[c1:c2])
        isamc50=N.average(isamyc[c1:c2])
        isamd50=N.average(isamyd[c1:c2])
        isama10=N.average(isamya[d1:d2])
        isamb10=N.average(isamyb[d1:d2])
        isamc10=N.average(isamyc[d1:d2])
        isamd10=N.average(isamyd[d1:d2])
        clma50=N.average(clmya[c1:c2])
        clmb50=N.average(clmyb[c1:c2])
        clmc50=N.average(clmyc[c1:c2])
        clmd50=N.average(clmyd[c1:c2])
        clma10=N.average(clmya[d1:d2])
        clmb10=N.average(clmyb[d1:d2])
        clmc10=N.average(clmyc[d1:d2])
        clmd10=N.average(clmyd[d1:d2])

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

	fig = plt.figure(figsize=(5,5))

	n_groups = 2
	ax = fig.add_subplot(111)
	index = N.arange(n_groups)
	bar_width = 0.2
	opacity = 0.8
	rects0 = plt.bar(index, fisama, bar_width,
                 alpha=opacity,
                 color='k',
                 label='CO2')
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

        plt.ylim(-20,40)
	plt.ylabel('percentage of yield changes (%)',fontsize=18)
	plt.title('2050s-2010s soybean',fontsize=18)
	plt.xticks(index + bar_width+0.2, ('ISAM','CLM'))
#	leg=plt.legend()
#	leg.get_frame().set_alpha(0.5)
	plt.tick_params(axis='both',labelsize=18)

	plt.tight_layout()

	plt.savefig('soy2050rcp85.png')
	plt.show()



