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
    yieldc1=yieldclm1/harea
    yieldc=yieldclm/harea
    yieldi=yieldisam/harea
    yieldic=yieldisamc/harea
    yieldic2=yieldisamc2/harea
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
a1=1981
a2=2008
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy1=N.zeros(x)
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
faoy=N.zeros(x)
faop=N.zeros(x)
#global
#faoy1=[17538,17586,16201,17143,19063,18199,19052,17047,18288,18958,18791,20379,19354,21833,20310,21313,21567,22561,21899,21689,23209,23009,22795,22437,23179]
#faop1=[88525040,92121684,79467061,90752915,101156845,94446350,100102463,93521958,107254261,108456443,103322537,114467678,115146195,136447144,126949296,130202572,144356733,160135008,157776694,161297357,178242605,181677687,190650566,205522972,214559072]

faoy1=[17538,17586,16201,17143,19063,18199,19052,17047,18288,18958,18725,20310,19361,21825,20312,21316,21565,22562,21901,21708,23061,22925,22796,22433,23175,23240,24370]
faop1=[88525040,92121684,79467061,90752915,101156845,94446350,100102463,93521958,107254261,108456443,102849072,113989333,115245015,136356715,126923558,130193213,144326230,160128838,157815456,161308456,177020696,180950793,190573602,205548177,214542816,221558978,219793049]


#us

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
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
        clmy1=N.zeros(x)
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


	mean_zumiy=runmean(x,zumiy)
	mean_isamy=runmean(x,isamy)
	mean_clmy=runmean(x,clmy)
	mean_zumip=runmean(x,zumip)
	mean_isamp=runmean(x,isamp)
	mean_clmp=runmean(x,clmp)
	mean_isamy1=runmean(x,isamy1)
	mean_isamp1=runmean(x,isamp1)
        mean_isamy2=runmean(x,isamy2)
        mean_isamp2=runmean(x,isamp2)
        mean_isamy3=runmean(x,isamy3)
        mean_isamp3=runmean(x,isamp3)

        mean_faoy=runmean(x,faoy)
        mean_faop=runmean(x,faop)
        mean_clmyn=runmean(x,clmyn)
        mean_clmpn=runmean(x,clmpn)


	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)

#plt.ylim((0,12))

	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

	ax.plot_date(xdates,faoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,isamy2,"ko-",label="ISAM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))
	ccp=scipy.stats.pearsonr(isamy,faoy)
	ccp1=scipy.stats.pearsonr(isamy1,faoy)
	ccp2=scipy.stats.pearsonr(clmy,faoy)
        ccp3=scipy.stats.pearsonr(isamy2,faoy)
        ccp4=scipy.stats.pearsonr(isamy3,faoy)
        ccp5=scipy.stats.pearsonr(clmy1,faoy)

        leg=plt.legend(['FAO', 'ISAM {:04.2f}'.format(ccp3[0])],fontsize=18)

#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.title("Yield",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)


	ax = fig.add_subplot(322)
        ax.plot_date(xdates,faop,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,isamp,"go-",label="NFix")
        ax.plot_date(xdates,isamp1,"yo-",label="NScale")
        ax.plot_date(xdates,isamp2,"ko-",label="NFao")
        ax.plot_date(xdates,isamp3,"co-",label="NAll")
        ax.plot_date(xdates,clmp,"bo-",label="CLM")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(isamp,faop)
	ccp1=scipy.stats.pearsonr(isamp1,faop)
	ccp2=scipy.stats.pearsonr(clmp,faop)
        ccp3=scipy.stats.pearsonr(isamp2,faop)
        ccp4=scipy.stats.pearsonr(isamp3,faop)

        leg=plt.legend(['FAO', 'Nfix {:04.2f}'.format(ccp[0]),'NScale {:04.2f}'.format(ccp1[0]),'NFao {:04.2f}'.format(ccp3[0]),'NAll {:04.2f}'.format(ccp4[0]),'CLM {:04.2f}'.format(ccp2[0])],fontsize=18)

#	leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0])],fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.title("Production",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(323)
        ax.plot_date(xdates,afaoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,aisamy2,"ko-",label="ISAM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamy,afaoy)
	ccp1=scipy.stats.pearsonr(aisamy1,afaoy)
	ccp2=scipy.stats.pearsonr(aclmy,afaoy)
        ccp3=scipy.stats.pearsonr(aisamy2,afaoy)
        ccp4=scipy.stats.pearsonr(aisamy3,faop)
        leg=plt.legend(['FAO', 'ISAM {:04.2f}'.format(ccp3[0])],fontsize=18)

#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
#        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0])],loc=4,fontsize=18)

#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(324)
	ax.plot_date(xdates,afaop,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,aisamp,"go-",label="NFix")
        ax.plot_date(xdates,aisamp1,"yo-",label="NScale")
        ax.plot_date(xdates,aisamp2,"ko-",label="NFao")
        ax.plot_date(xdates,aisamp3,"co-",label="NAll")
        ax.plot_date(xdates,aclmp,"bo-",label="CLM")

        ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamp,afaop)
	ccp1=scipy.stats.pearsonr(aisamp1,afaop)
	ccp2=scipy.stats.pearsonr(aclmp,afaop)
        ccp3=scipy.stats.pearsonr(aisamp2,afaop)
        ccp4=scipy.stats.pearsonr(aisamp3,afaop)

#        ccp3=scipy.stats.pearsonr(aclmpn,afaop)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
#        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0])],fontsize=18)
        leg=plt.legend(['FAO', 'Nfix {:04.2f}'.format(ccp[0]),'NScale {:04.2f}'.format(ccp1[0]),'NFao {:04.2f}'.format(ccp3[0]),'NAll {:04.2f}'.format(ccp4[0]),'CLM {:04.2f}'.format(ccp2[0])],fontsize=18)

#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(325)
        ax.plot_date(xdates,mean_faoy,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,mean_isamy2,"ko-",label="ISAM")
        ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamy=N.zeros([x-2])
	mean1_faoy=N.zeros([x-2])
	mean1_isamy1=N.zeros([x-2])
        mean1_isamy2=N.zeros([x-2])
        mean1_isamy3=N.zeros([x-2])
	mean1_clmy=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_isamy[i]=mean_isamy[i]
		mean1_faoy[i]=mean_faoy[i]
		mean1_isamy1[i]=mean_isamy1[i]
		mean1_clmy[i]=mean_clmy[i]
                mean1_isamy2[i]=mean_isamy2[i]
                mean1_isamy3[i]=mean_isamy3[i]

	ccp=scipy.stats.pearsonr(mean1_isamy,mean1_faoy)
	ccp1=scipy.stats.pearsonr(mean1_isamy1,mean1_faoy)
	ccp2=scipy.stats.pearsonr(mean1_clmy,mean1_faoy)
        ccp3=scipy.stats.pearsonr(mean1_isamy2,mean1_faoy)
        ccp4=scipy.stats.pearsonr(mean1_isamy3,mean1_faoy)

#        ccp3=scipy.stats.pearsonr(mean1_clmyn,mean1_faoy)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
#        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
        leg=plt.legend(['FAO','ISAM {:04.2f}'.format(ccp3[0])],fontsize=18,loc=4)

	leg.get_frame().set_alpha(0.5)


#	plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(326)
        ax.plot_date(xdates,mean_faop,"ro-",label="FAO",linewidth=2)
        ax.plot_date(xdates,mean_isamp,"go-",label="NFix")
        ax.plot_date(xdates,mean_isamp1,"yo-",label="NScale")
        ax.plot_date(xdates,mean_isamp2,"ko-",label="NFao")
        ax.plot_date(xdates,mean_isamp3,"co-",label="NAll")
        ax.plot_date(xdates,mean_clmp,"bo-",label="CLM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamp=N.zeros([x-2])
	mean1_faop=N.zeros([x-2])
	mean1_isamp1=N.zeros([x-2])
        mean1_isamp2=N.zeros([x-2])
        mean1_isamp3=N.zeros([x-2])

	mean1_clmp=N.zeros([x-2])
        mean1_clmpn=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_isamp[i]=mean_isamp[i]
		mean1_faop[i]=mean_faop[i]
		mean1_isamp1[i]=mean_isamp1[i]
		mean1_clmp[i]=mean_clmp[i]
                mean1_clmpn[i]=mean_clmpn[i]
                mean1_isamp2[i]=mean_isamp2[i]
                mean1_isamp3[i]=mean_isamp3[i]

	ccp=scipy.stats.pearsonr(mean1_isamp,mean1_faop)
	ccp1=scipy.stats.pearsonr(mean1_isamp1,mean1_faop)
	ccp2=scipy.stats.pearsonr(mean1_clmp,mean1_faop)
        ccp3=scipy.stats.pearsonr(mean1_isamp2,mean1_faop)
        ccp4=scipy.stats.pearsonr(mean1_isamp3,mean1_faop)

#	leg=plt.legend(['FAO', 'ISAM-NCRP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
#        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
        leg=plt.legend(['FAO', 'Nfix {:04.2f}'.format(ccp[0]),'NScale {:04.2f}'.format(ccp1[0]),'NFao {:04.2f}'.format(ccp3[0]),'NAll {:04.2f}'.format(ccp4[0]),'CLM {:04.2f}'.format(ccp2[0])],fontsize=18)

	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
	plt.xlabel("Year ",fontsize=18)
	plt.ylabel("Maize production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)


	plt.savefig('soy_faonclm_{0}_paper1.png'.format(name1),dpi=600,bbox_inches='tight')
plt.show()


