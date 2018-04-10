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
#area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
#gridarea = area.variables['cell_area'][:,:]
#gridlon = area.variables['lon'][:]

#gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)



def annualyield(year,couna,counb):
    if year <=2005:
    	bb=year-1901
	aa=year-1901
    else:
	bb=104
        aa=year-1901
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
    maizetor=maitrop+maitemp
    maizetoi=maitropi+maitempi
    maizeto = maitrop+maitemp+maitropi+maitempi
    
    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/clm/clm45his_soyscaleifyield.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]
    clmtropf= ma.masked_where(clmtropf<=0,clmtropf)
    clmtropf=ma.filled(clmtropf, fill_value=0.)

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/clm/clm45his_soyscaleiyield.nc','r')
    clmtropf1 = clm2.variables['yield'][bb,:,:]
    clmtropf1= ma.masked_where(clmtropf1<=0,clmtropf1)
    clmtropf1=ma.filled(clmtropf1, fill_value=0.)



    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_soyscaleiyield_heat.nc','r')
    isammai = clm3n.variables['yield'][bb,:,:]
    isammai= ma.masked_where(isammai<=0,isammai)
    isammai=ma.filled(isammai,fill_value=0.)


    clm3n1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_soyscaleifyield_heat.nc','r')
    isamcrumai = clm3n1.variables['yield'][bb,:,:]
    isamcrumai= ma.masked_where(isamcrumai<=0,isamcrumai)
    isamcrumai=ma.filled(isamcrumai, fill_value=0.)

    clm3n2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertfao/new1/isamhiscru_soyscaleiyield_fertfao_new1_long.nc','r')
    isamcrumai2 = clm3n2.variables['yield'][aa,:,:]
    isamcrumai2= ma.masked_where(isamcrumai2<=0,isamcrumai2)
    isamcrumai2=ma.filled(isamcrumai2, fill_value=0.)

    clm3n3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertall/isamhis_soyscaleiyield_heat_fertall.nc','r')
    isamcrumai3 = clm3n3.variables['yield'][bb,:,:]
    isamcrumai3= ma.masked_where(isamcrumai3<=0,isamcrumai3)
    isamcrumai3=ma.filled(isamcrumai3, fill_value=0.)
 
    yieldclm1=0
    yieldclm=0.
    yieldisam=0.
    yieldisamc=0.
    yieldisamc2=0.
    yieldisamc3=0.

    a=0
    harea=0
    yieldc1=0
    yieldc=0.
    yieldi=0.
    yieldic=0.
    yieldic2=0.
    yieldic3=0.

    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                harea=maizeto[xx,yy]*gridarea[xx,yy]+ harea
                yieldclm1=(clmtropf1[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldclm1
                yieldclm=(clmtropf[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldclm
                yieldisam=(isammai[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldisam
                yieldisamc=(isamcrumai[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldisamc
                yieldisamc2=(isamcrumai2[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldisamc2
                yieldisamc3=(isamcrumai3[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy])+yieldisamc3

                a=a+1
    yieldc1=N.average(clmtropf1,weights=maizeto)
    yieldc=N.average(clmtropf,weights=maizeto)
    yieldi=yieldisam/harea
    yieldic=yieldisamc/harea
    yieldic2=N.average(isamcrumai2,weights=maizeto)
    yieldic3=yieldisamc3/harea

    return harea, yieldc,yieldi,yieldic,yieldclm,yieldisam,yieldisamc,yieldic2,yieldic3,yieldisamc2,yieldisamc3,yieldc1
    #return harea
def runmean(x,input):
        import pandas as pd
        #mean_zumiy1 = pd.rolling_mean(zumiy, window=5).shift(-2)
        meanout=pd.rolling_mean(input, window=5, center=True)
        mean_zumiy1=pd.rolling_mean(input, window=3, center=True)
        #print mean_zumiy1
        #print meanout
        meanout[1]=mean_zumiy1[1]
        meanout[x-2]=mean_zumiy1[x-2]
        meanout1=input-meanout
        return meanout1

#illzmui only 1983~2005
#a1=1961
#a2=2017
#x=a2-a1
a1=1902
a2=2017
x=a2-a1
x1=2017-1961

toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy= N.zeros(x)
clmyn=N.zeros(x)
clmpn=N.zeros(x)
isamy= N.zeros(x)
zumip= N.zeros(x)
clmp= N.zeros(x)
isamp= N.zeros(x)
isamy1= N.zeros(x)
isamp1= N.zeros(x)
isamy2= N.zeros(x)
isamp2= N.zeros(x)
isamy3= N.zeros(x)
isamp3= N.zeros(x)
faoy=N.zeros(x1)
faop=N.zeros(x1)
#global!1981-2007


faoy1=[11287,11384,11561,11323,12279,13723,13494,14345,14523,14800,15189,14896,15851,14079,16574,15443,17561,16264,17492,16001,17538,17586,16201,17143,19063,18199,19052,17047,18288,18958,18725,20310,19361,21825,20312,21316,21565,22562,21901,21708,23061,22925,22796,22433,23175,23240,24370,23983,22481,25781,25197,22894,24999,26044,26757,27556]
faop1=[26883158,27120148,28207071,29079547,31704923,36412918,37933319,41420481,41967347,43696887,45618769,47257207,59267889,52639832,64248541,57399226,73854802,75449966,88698341,81040368,88525040,92121684,79467061,90752915,101156845,94446350,100102463,93521958,107254261,108456443,102849072,113989333,115245015,136356715,126923558,130193213,144326230,160128838,157815456,161308456,177020696,180950793,190573602,205548177,214542816,221558978,219793049,231148449,223267546,264942943,261436988,241185392,277537300,306373463,323204743,334894085]





c=0.0001
faoy=N.multiply(faoy1,c)
faop=N.multiply(faop1,1)
#name=["Global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","Southafrica"]
#range1=[10100,11501,41501,20301,20101,11101,42901,47800,50300,42200,33900]
#range2=[50700,11550,41529,20327,20124,11132,42929,47800,50300,42200,33900]

#name=["Italy","Canada","Vietnam","Hungary","Romania","Philippines","Thailand","Chile","Spain","Nigeria","Germany"]
#range1=[43400,10201,48200,42700,46000,50700,47500,20400,46800,33400,42500]
#range2=[43400,10212,48200,42700,46000,50700,47500,20400,46800,33400,42500]

name=["Global"]
range1=[10100]
range2=[50700]
#name=["US"]
#range1=[11501]
#range2=[11550]

for i, name1 in enumerate(name):
        a=0
        clmy1=N.zeros(x)
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
	clmy= N.zeros(x)
	isamy= N.zeros(x)
	zumip= N.zeros(x)
	clmp= N.zeros(x)
	isamp= N.zeros(x)
	isamy1= N.zeros(x)
	isamp1= N.zeros(x)
        isamy2= N.zeros(x)
        isamp2= N.zeros(x)
        isamy3= N.zeros(x)
        isamp3= N.zeros(x)

        clmyn=N.zeros(x)
        clmpn=N.zeros(x)
	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])
            clmy1[a]=reu[11]
	    toarea[a]=reu[0]
	    isamy1[a]=reu[3]
	    clmy[a]=reu[1]
	    isamp1[a]=reu[5]
	    clmp[a]=reu[4]
	    isamy[a]=reu[2]
	    isamp[a]=reu[6]
            isamy2[a]=reu[7]
            isamp2[a]=reu[9]
            isamy3[a]=reu[8]
            isamp3[a]=reu[10]

	    a=a+1

	azumiy=zumiy-N.average(zumiy)
	aisamy=isamy-N.average(isamy)
	aclmy=clmy-N.average(clmy)
	azumip=zumip-N.average(zumip)
	aisamp=isamp-N.average(isamp)
	aclmp=clmp-N.average(clmp)
	aisamy1=isamy1-N.average(isamy1)
	aisamp1=isamp1-N.average(isamp1)
        aisamy2=isamy2-N.average(isamy2)
        aisamp2=isamp2-N.average(isamp2)
        aisamy3=isamy3-N.average(isamy3)
        aisamp3=isamp3-N.average(isamp3)

        afaoy=faoy-N.average(faoy)
        afaop=faop-N.average(faop)
        aclmyn=clmyn-N.average(clmyn)
        aclmpn=clmpn-N.average(clmpn)




	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)
	x1=range(1961,2017)
#plt.ylim((0,12))
        xdates1 = [datetime.datetime.strptime(str(int(date)),'%Y') for date in x1]
	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

	ax.plot_date(xdates1,faoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,isamy2,"ko-",label="ISAM")
        ax.plot_date(xdates,clmy1,"bo-",label="CLM")
        ax.plot_date(xdates,clmy,"go-",label="CLM-Scale")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

        leg=plt.legend(fontsize=18)

	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.title("Yield",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)




	plt.savefig('soy_faoclm_{0}_1961_2016_a1.png'.format(name1),dpi=600,bbox_inches='tight')
plt.show()


