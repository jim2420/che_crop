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
    lonlon = region1.variables['lon'][:]
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
    maizeto,lon1 = shiftgrid(0.5,maizeto,lonlon,start=True)
#    print lon1
    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
    clmtropf = clm.variables['totalyield'][bb,:,:]
    clmtropf= ma.masked_where(clmtropf<=0,clmtropf)
    clmtropf=ma.filled(clmtropf, fill_value=0.)

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
    clmtropf1 = clm2.variables['totalyield'][bb,:,:]
    clmtropf1= ma.masked_where(clmtropf1<=0,clmtropf1)
    clmtropf1=ma.filled(clmtropf1, fill_value=0.)



    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
    isammai = clm3n.variables['totalyield'][bb,:,:]
    isammai= ma.masked_where(isammai<=0,isammai)
    isammai=ma.filled(isammai,fill_value=0.)


    clm3n1=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/isamhis_maiscaleifyield_heat.nc','r')
    isamcrumai = clm3n1.variables['yield'][bb,:,:]
    isamcrumai= ma.masked_where(isamcrumai<=0,isamcrumai)
    isamcrumai=ma.filled(isamcrumai, fill_value=0.)

    clm3n2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
    isamcrumai2 = clm3n2.variables['totalyield'][aa,:,:]
    isamcrumai2= ma.masked_where(isamcrumai2<=0,isamcrumai2)
    isamcrumai2=ma.filled(isamcrumai2, fill_value=0.)

    clm3n3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/plot/finalyield/isam/heat/fertall/isamhis_maiscaleiyield_heat_fertall.nc','r')
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
    yieldi=N.average(isammai,weights=maizeto)
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
#faoy1=[34933.0,36090.0,29452.0,35256.0,37202.0,36279.,34863.,31001.,36186.,36907.,36996.,39064.,36305.,41120.,38100.,42061.,41511.,44354.,44255.,43236.,44775.,43884.,44610.,49452.,48196.,47719.,49980.]
#faop1=[446772517,448932280,347082034,450449992,485527301,478176622,453115794,403050234,476874503,483620724,494393020,533774898,477207493,568650520,517286851,586134845,584401847,615072804,607426254,592030667,6155143531,603544019,645048171,729511789,714185792,707932497,793055503]


faoy1=[19423.0,19796.0,20311.0,19937.0,21237.0,22078.0,24236.0,22894.0,24206.0,23509.0,26534.0,26867.0,27217.0,25565.0,28132.0,28371.0,29668.0,31556.0,33847.0,31534.0,34933.0,36090.0,29452.0,35256.0,37202.0,36279.,34863.,31001.,36186.,36907.,36996.,39064.,36305.,41120.,38100.,42061.,41511.,44354.,44255.,43236.,44775.,43884.,44610.,49452.,48196.,47719.,49980.,50828.0,51635.0,51903.0,51751.0,48893.0,54611.0,56229.0,55379.0,56401.0]
faop1=[205027583,204876937,220228333,215172627,226544256,245599160,272548473,255670551,269506068,265831145,313622622,308826290,318290469,306427347,341751971,352395866,371593355,393600091,418622993,396623388,446772517,448932280,347082034,450449992,485527301,478176622,453115794,403050234,476874503,483620724,494393020,533774898,477207493,568650520,517286851,586134845,584401847,615072804,607426254,592030667,6155143531,603544019,645048171,729511789,714185792,707932497,793055503,829236755,820069886,851348928,886007062,874240510,1015400446,1038330655,1010609468,1060107470]


#us
#faoy1=[68380,71082,50893,66981,74074,74930,75225,53109,72978,74380,68172,82526,63211,86997,71230,79777,79522,84382,83979,85910,86733,81179,89246,100636,92852]
#faop1=[206222000,209180000,106030000,194880000,225453008,208943008,181142000,125194000,191319008,201532000,189866496,240719008,160984992,255292992,187968992,234527008,233867008,247882000,239548992,251852210,241375035,227765357,256227304,299873563,282260662]

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
        ax.plot_date(xdates,isamy2,"ko-",label="NO-NF")
        ax.plot_date(xdates,clmy1,"bo-",label="Base")
        ax.plot_date(xdates,clmy,"go-",label="NO-I")
        ax.plot_date(xdates,isamy,"yo-",label="NO-I & NF")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

        leg=plt.legend(fontsize=18)

	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.title("Yield",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)




	plt.savefig('soy_faoclm_{0}_1961_2016_a3.png'.format(name1),dpi=600,bbox_inches='tight')
plt.show()


