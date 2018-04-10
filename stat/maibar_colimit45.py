from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
from statsmodels.stats.weightstats import DescrStatsW
from scipy.stats import ttest_ind
from matplotlib.markers import TICKDOWN
import datetime
from matplotlib.dates import DateFormatter
from scipy import stats
import datetime


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['maizeData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
latnc=N.flipud(latnc)




lon2,lat2 = N.meshgrid(gridlon,gridlat)



region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
#maitrop = region1.variables['soy_trop'][99,:,:]
#maitemp = region1.variables['soy_temp'][99,:,:]
#maitropi=region1.variables['soy_trop_irrig'][99,:,:]
#maitempi=region1.variables['soy_temp_irrig'][99,:,:]
maitrop = region1.variables['maize_trop'][0:95,:,:]
maitemp = region1.variables['maize_temp'][0:95,:,:]
maitropi=region1.variables['maize_trop_irrig'][0:95,:,:]
maitempi=region1.variables['maize_temp_irrig'][0:95,:,:]
lonlon=region1.variables['lon'][:]


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
maizetropo=maitrop+maitropi
maizetempo=maitemp+maitempi

 
region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
maitrop85 = region1.variables['maize_trop'][0:95,:,:]
maitemp85 = region1.variables['maize_temp'][0:95,:,:]
maitropi85=region1.variables['maize_trop_irrig'][0:95,:,:]
maitempi85=region1.variables['maize_temp_irrig'][0:95,:,:]


maitrop85=ma.masked_where(maitrop85<=0,maitrop85)
maitrop85=ma.filled(maitrop85, fill_value=0.)
maitemp85=ma.masked_where(maitemp85<=0,maitemp85)
maitemp85=ma.filled(maitemp85, fill_value=0.)

maitropi85=ma.masked_where(maitropi85<=0,maitropi85)
maitropi85=ma.filled(maitropi85, fill_value=0.)
maitempi85=ma.masked_where(maitempi85<=0,maitempi85)
maitempi85=ma.filled(maitempi85, fill_value=0.)
maizetor85=maitrop85+maitemp85
maizetoi85=maitropi85+maitempi85
maizeto85 = maitrop85+maitemp85+maitropi85+maitempi85
maizetropo85=maitrop85+maitropi85
maizetempo85=maitemp85+maitempi85







dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ynew = dat.variables['totalyield'][:,:,:]
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_irrig_fert_0.5x0.5.nc','r')
iyield2ynew = dat2.variables['totalyield'][:,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_rf_fert_0.5x0.5.nc','r')
iyield3ynew = dat3.variables['totalyield'][:,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_rf_fert_0.5x0.5.nc','r')
iyield4ynew = dat4.variables['totalyield'][:,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
iyield5ynew = dat5.variables['totalyield'][:,:,:]

dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
iyield6ynew = dat6.variables['totalyield'][:,:,:]


dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
iyield1ynew85 = dat.variables['totalyield'][:,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
iyield2ynew85 = dat2.variables['totalyield'][:,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_rf_fert_0.5x0.5.nc','r')
iyield3ynew85 = dat3.variables['totalyield'][:,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_rf_fert_0.5x0.5.nc','r')
iyield4ynew85 = dat4.variables['totalyield'][:,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
iyield5ynew85 = dat5.variables['totalyield'][:,:,:]

dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
iyield6ynew85 = dat6.variables['totalyield'][:,:,:]


dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_irrig_fert_0.5x0.5_n.nc','r')
iyield1ynewn = dat.variables['total_gfna'][:,:,:]
iyield1ynewi = dat.variables['total_gWS'][:,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_irrig_fert_0.5x0.5_n.nc','r')
iyield2ynewn = dat2.variables['total_gfna'][:,:,:]
iyield2ynewi = dat2.variables['total_gWS'][:,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_rf_fert_0.5x0.5_n.nc','r')
iyield3ynewn = dat3.variables['total_gfna'][:,:,:]
iyield3ynewi = dat3.variables['total_gWS'][:,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_rf_fert_0.5x0.5_n.nc','r')
iyield4ynewn = dat4.variables['total_gfna'][:,:,:]
iyield4ynewi = dat4.variables['total_gWS'][:,:,:]


dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetemp_rcp45_co2_rf_nofert_0.5x0.5_n.nc','r')
iyield5ynewn = dat5.variables['total_gfna'][:,:,:]
iyield5ynewi = dat5.variables['total_gWS'][:,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp45/heat/new1/maizetrop_rcp45_co2_rf_nofert_0.5x0.5_n.nc','r')
iyield6ynewn = dat6.variables['total_gfna'][:,:,:]
iyield6ynewi = dat6.variables['total_gWS'][:,:,:]



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_irrig_fert_0.5x0.5_n.nc','r')
iyield1ynewn85 = dat.variables['total_gfna'][:,:,:]
iyield1ynewi85 = dat.variables['total_gWS'][:,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_irrig_fert_0.5x0.5_n.nc','r')
iyield2ynewn85 = dat2.variables['total_gfna'][:,:,:]
iyield2ynewi85 = dat2.variables['total_gWS'][:,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_rf_fert_0.5x0.5_n.nc','r')
iyield3ynewn85 = dat3.variables['total_gfna'][:,:,:]
iyield3ynewi85 = dat3.variables['total_gWS'][:,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_rf_fert_0.5x0.5_n.nc','r')
iyield4ynewn85 = dat4.variables['total_gfna'][:,:,:]
iyield4ynewi85 = dat4.variables['total_gWS'][:,:,:]


dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetemp_rcp85_co2_rf_nofert_0.5x0.5_n.nc','r')
iyield5ynewn85 = dat5.variables['total_gfna'][:,:,:]
iyield5ynewi85 = dat5.variables['total_gWS'][:,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamrcp85/heat/new1/maizetrop_rcp85_co2_rf_nofert_0.5x0.5_n.nc','r')
iyield6ynewn85 = dat6.variables['total_gfna'][:,:,:]
iyield6ynewi85 = dat6.variables['total_gWS'][:,:,:]





iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)

iyield1ynewi= ma.masked_where(iyield1ynew<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where(iyield2ynew<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where(iyield3ynew<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where(iyield4ynew<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where(iyield5ynew<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where(iyield6ynew<=0.,iyield6ynewi)

iyield1ynewn= ma.masked_where(iyield1ynew<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where(iyield2ynew<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where(iyield3ynew<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where(iyield4ynew<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where(iyield5ynew<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where(iyield6ynew<=0.,iyield6ynewn)

iyield1ynew85= ma.masked_where(iyield1ynew85<=0.,iyield1ynew85)
iyield2ynew85= ma.masked_where(iyield2ynew85<=0.,iyield2ynew85)
iyield3ynew85= ma.masked_where(iyield3ynew85<=0.,iyield3ynew85)
iyield4ynew85= ma.masked_where(iyield4ynew85<=0.,iyield4ynew85)
iyield5ynew85= ma.masked_where(iyield5ynew85<=0.,iyield5ynew85)
iyield6ynew85= ma.masked_where(iyield6ynew85<=0.,iyield6ynew85)

iyield1ynewi85= ma.masked_where(iyield1ynew85<=0.,iyield1ynewi85)
iyield2ynewi85= ma.masked_where(iyield2ynew85<=0.,iyield2ynewi85)
iyield3ynewi85= ma.masked_where(iyield3ynew85<=0.,iyield3ynewi85)
iyield4ynewi85= ma.masked_where(iyield4ynew85<=0.,iyield4ynewi85)
iyield5ynewi85= ma.masked_where(iyield5ynew85<=0.,iyield5ynewi85)
iyield6ynewi85= ma.masked_where(iyield6ynew85<=0.,iyield6ynewi85)

iyield1ynewn85= ma.masked_where(iyield1ynew85<=0.,iyield1ynewn85)
iyield2ynewn85= ma.masked_where(iyield2ynew85<=0.,iyield2ynewn85)
iyield3ynewn85= ma.masked_where(iyield3ynew85<=0.,iyield3ynewn85)
iyield4ynewn85= ma.masked_where(iyield4ynew85<=0.,iyield4ynewn85)
iyield5ynewn85= ma.masked_where(iyield5ynew85<=0.,iyield5ynewn85)
iyield6ynewn85= ma.masked_where(iyield6ynew85<=0.,iyield6ynewn85)



maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonlon,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonlon,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonlon,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonlon,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonlon,start=True)

maizeto185,lonisam2=shiftgrid(0.5,maizeto85,lonlon,start=True)
maizeto1i85,lonisam2=shiftgrid(0.5,maitropi85,lonlon,start=True)
maizeto1r85,lonisam2=shiftgrid(0.5,maitrop85,lonlon,start=True)
maizete1i85,lonisam2=shiftgrid(0.5,maitempi85,lonlon,start=True)
maizete1r85,lonisam2=shiftgrid(0.5,maitemp85,lonlon,start=True)


iyield1ynewi= ma.masked_where((maizete1i+maizete1r)<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where((maizete1i+maizete1r)<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where((maizete1i+maizete1r)<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield6ynewi)

iyield1ynewn= ma.masked_where((maizete1i+maizete1r)<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where((maizete1i+maizete1r)<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where((maizete1i+maizete1r)<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where((maizeto1i+maizeto1r)<=0.,iyield6ynewn)




iyield1ynewi85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield1ynewi85)
iyield2ynewi85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield2ynewi85)
iyield3ynewi85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield3ynewi85)
iyield4ynewi85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield4ynewi85)
iyield5ynewi85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield5ynewi85)
iyield6ynewi85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield6ynewi85)

iyield1ynewn85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield1ynewn85)
iyield2ynewn85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield2ynewn85)
iyield3ynewn85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield3ynewn85)
iyield4ynewn85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield4ynewn85)
iyield5ynewn85= ma.masked_where((maizete1i85+maizete1r85)<=0.,iyield5ynewn85)
iyield6ynewn85= ma.masked_where((maizeto1i85+maizeto1r85)<=0.,iyield6ynewn85)



iyield1ynewn=ma.filled(iyield1ynewn, fill_value=0)
iyield2ynewn=ma.filled(iyield2ynewn, fill_value=0)
iyield3ynewn=ma.filled(iyield3ynewn, fill_value=0)
iyield4ynewn=ma.filled(iyield4ynewn, fill_value=0)
iyield5ynewn=ma.filled(iyield5ynewn, fill_value=0)
iyield6ynewn=ma.filled(iyield6ynewn, fill_value=0)


iyield1ynewi=ma.filled(iyield1ynewi, fill_value=0)
iyield2ynewi=ma.filled(iyield2ynewi, fill_value=0)
iyield3ynewi=ma.filled(iyield3ynewi, fill_value=0)
iyield4ynewi=ma.filled(iyield4ynewi, fill_value=0)
iyield5ynewi=ma.filled(iyield5ynewi, fill_value=0)
iyield6ynewi=ma.filled(iyield6ynewi, fill_value=0)

iyield1ynewn85=ma.filled(iyield1ynewn85, fill_value=0)
iyield2ynewn85=ma.filled(iyield2ynewn85, fill_value=0)
iyield3ynewn85=ma.filled(iyield3ynewn85, fill_value=0)
iyield4ynewn85=ma.filled(iyield4ynewn85, fill_value=0)
iyield5ynewn85=ma.filled(iyield5ynewn85, fill_value=0)
iyield6ynewn85=ma.filled(iyield6ynewn85, fill_value=0)

iyield1ynewi85=ma.filled(iyield1ynewi85, fill_value=0)
iyield2ynewi85=ma.filled(iyield2ynewi85, fill_value=0)
iyield3ynewi85=ma.filled(iyield3ynewi85, fill_value=0)
iyield4ynewi85=ma.filled(iyield4ynewi85, fill_value=0)
iyield5ynewi85=ma.filled(iyield5ynewi85, fill_value=0)
iyield6ynewi85=ma.filled(iyield6ynewi85, fill_value=0)

iyield1ynewn85=iyield1ynewn85+iyield2ynewn85
iyield3ynewn85=iyield3ynewn85+iyield4ynewn85
iyield5ynewn85=iyield5ynewn85+iyield6ynewn85
iyield1ynewi85=iyield1ynewi85+iyield2ynewi85
iyield3ynewi85=iyield3ynewi85+iyield4ynewi85
iyield5ynewi85=iyield5ynewi85+iyield6ynewi85

iyield1ynewn=iyield1ynewn+iyield2ynewn
iyield3ynewn=iyield3ynewn+iyield4ynewn
iyield5ynewn=iyield5ynewn+iyield6ynewn
iyield1ynewi=iyield1ynewi+iyield2ynewi
iyield3ynewi=iyield3ynewi+iyield4ynewi
iyield5ynewi=iyield5ynewi+iyield6ynewi



iyield1ynewn= ma.masked_where(maizeto1<=0.,iyield1ynewn)
iyield3ynewn= ma.masked_where(maizeto1<=0.,iyield3ynewn)
iyield5ynewn= ma.masked_where(maizeto1<=0.,iyield5ynewn)
iyield1ynewi= ma.masked_where(maizeto1<=0.,iyield1ynewi)
iyield3ynewi= ma.masked_where(maizeto1<=0.,iyield3ynewi)
iyield5ynewi= ma.masked_where(maizeto1<=0.,iyield5ynewi)


iyield1ynewn85= ma.masked_where(maizeto185<=0.,iyield1ynewn85)
iyield3ynewn85= ma.masked_where(maizeto185<=0.,iyield3ynewn85)
iyield5ynewn85= ma.masked_where(maizeto185<=0.,iyield5ynewn85)
iyield1ynewi85= ma.masked_where(maizeto185<=0.,iyield1ynewi85)
iyield3ynewi85= ma.masked_where(maizeto185<=0.,iyield3ynewi85)
iyield5ynewi85= ma.masked_where(maizeto185<=0.,iyield5ynewi85)



iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)

'''
iyield1ynewi= ma.masked_where(iyield1ynewi<=0.,iyield1ynewi)
iyield2ynewi= ma.masked_where(iyield2ynewi<=0.,iyield2ynewi)
iyield3ynewi= ma.masked_where(iyield3ynewi<=0.,iyield3ynewi)
iyield4ynewi= ma.masked_where(iyield4ynewi<=0.,iyield4ynewi)
iyield5ynewi= ma.masked_where(iyield5ynewi<=0.,iyield5ynewi)
iyield6ynewi= ma.masked_where(iyield6ynewi<=0.,iyield6ynewi)

iyield1ynewn= ma.masked_where(iyield1ynewn<=0.,iyield1ynewn)
iyield2ynewn= ma.masked_where(iyield2ynewn<=0.,iyield2ynewn)
iyield3ynewn= ma.masked_where(iyield3ynewn<=0.,iyield3ynewn)
iyield4ynewn= ma.masked_where(iyield4ynewn<=0.,iyield4ynewn)
iyield5ynewn= ma.masked_where(iyield5ynewn<=0.,iyield5ynewn)
iyield6ynewn= ma.masked_where(iyield6ynewn<=0.,iyield6ynewn)
'''











#print iyield6ynewn
#iyield6ynewn=iyield6ynewn[N.logical_not(N.isnan(iyield6ynewn))]
#iyield6ynewn = ma.masked_array(iyield6ynewn, N.isnan(iyield6ynewn))
iyield1ynewn[N.isnan(iyield1ynewn)]=-1
iyield1ynewn= ma.masked_where(iyield1ynewn<=0.0,iyield1ynewn)
iyield2ynewn[N.isnan(iyield2ynewn)]=-1
iyield2ynewn= ma.masked_where(iyield2ynewn<=0.0,iyield2ynewn)
iyield3ynewn[N.isnan(iyield3ynewn)]=-1
iyield3ynewn= ma.masked_where(iyield3ynewn<=0.0,iyield3ynewn)
iyield4ynewn[N.isnan(iyield4ynewn)]=-1
iyield4ynewn= ma.masked_where(iyield4ynewn<=0.0,iyield4ynewn)
iyield5ynewn[N.isnan(iyield5ynewn)]=-1
iyield5ynewn= ma.masked_where(iyield5ynewn<=0.0,iyield5ynewn)
iyield6ynewn[N.isnan(iyield6ynewn)]=-1
iyield6ynewn= ma.masked_where(iyield6ynewn<=0.0,iyield6ynewn)

iyield1ynewi[N.isnan(iyield1ynewi)]=-1
iyield1ynewi= ma.masked_where(iyield1ynewi<=0.0,iyield1ynewi)
iyield2ynewi[N.isnan(iyield2ynewi)]=-1
iyield2ynewi= ma.masked_where(iyield2ynewi<=0.0,iyield2ynewi)
iyield3ynewi[N.isnan(iyield3ynewi)]=-1
iyield3ynewi= ma.masked_where(iyield3ynewi<=0.0,iyield3ynewi)
iyield4ynewi[N.isnan(iyield4ynewi)]=-1
iyield4ynewi= ma.masked_where(iyield4ynewi<=0.0,iyield4ynewi)
iyield5ynewi[N.isnan(iyield5ynewi)]=-1
iyield5ynewi= ma.masked_where(iyield5ynewi<=0.0,iyield5ynewi)
iyield6ynewi[N.isnan(iyield6ynewi)]=-1
iyield6ynewi= ma.masked_where(iyield6ynewi<=0.0,iyield6ynewi)


iyield1ynewn85[N.isnan(iyield1ynewn85)]=-1
iyield1ynewn85= ma.masked_where(iyield1ynewn85<=0.0,iyield1ynewn85)
iyield2ynewn85[N.isnan(iyield2ynewn85)]=-1
iyield2ynewn85= ma.masked_where(iyield2ynewn85<=0.0,iyield2ynewn85)
iyield3ynewn85[N.isnan(iyield3ynewn85)]=-1
iyield3ynewn85= ma.masked_where(iyield3ynewn85<=0.0,iyield3ynewn85)
iyield4ynewn85[N.isnan(iyield4ynewn85)]=-1
iyield4ynewn85= ma.masked_where(iyield4ynewn85<=0.0,iyield4ynewn85)
iyield5ynewn85[N.isnan(iyield5ynewn85)]=-1
iyield5ynewn85= ma.masked_where(iyield5ynewn85<=0.0,iyield5ynewn85)
iyield6ynewn85[N.isnan(iyield6ynewn85)]=-1
iyield6ynewn85= ma.masked_where(iyield6ynewn85<=0.0,iyield6ynewn85)


iyield1ynewi85[N.isnan(iyield1ynewi85)]=-1
iyield1ynewi85= ma.masked_where(iyield1ynewi85<=0.0,iyield1ynewi85)
iyield2ynewi85[N.isnan(iyield2ynewi85)]=-1
iyield2ynewi85= ma.masked_where(iyield2ynewi85<=0.0,iyield2ynewi85)
iyield3ynewi85[N.isnan(iyield3ynewi85)]=-1
iyield3ynewi85= ma.masked_where(iyield3ynewi85<=0.0,iyield3ynewi85)
iyield4ynewi85[N.isnan(iyield4ynewi85)]=-1
iyield4ynewi85= ma.masked_where(iyield4ynewi85<=0.0,iyield4ynewi85)
iyield5ynewi85[N.isnan(iyield5ynewi85)]=-1
iyield5ynewi85= ma.masked_where(iyield5ynewi85<=0.0,iyield5ynewi85)
iyield6ynewi85[N.isnan(iyield6ynewi85)]=-1
iyield6ynewi85= ma.masked_where(iyield6ynewi85<=0.0,iyield6ynewi85)





#print iyield6ynewn

#print iyield1ynew.shape
maizeto2=N.zeros((95,360,720))
maizeto2r=N.zeros((95,360,720))
maizeto2i=N.zeros((95,360,720))
maizete2r=N.zeros((95,360,720))
maizete2i=N.zeros((95,360,720))

maizeto285=N.zeros((95,360,720))
maizeto2r85=N.zeros((95,360,720))
maizeto2i85=N.zeros((95,360,720))
maizete2r85=N.zeros((95,360,720))
maizete2i85=N.zeros((95,360,720))


for i in range(0,95):
	for x in range(0,360):
		for y in range(0,720):
			maizeto2[i,x,y]=maizeto1[i,x,y]
                        maizeto2r[i,x,y]=maizeto1r[i,x,y]
                        maizeto2i[i,x,y]=maizeto1i[i,x,y]
                        maizete2r[i,x,y]=maizete1r[i,x,y]
                        maizete2i[i,x,y]=maizete1i[i,x,y]

                        maizeto285[i,x,y]=maizeto185[i,x,y]
                        maizeto2r85[i,x,y]=maizeto1r85[i,x,y]
                        maizeto2i85[i,x,y]=maizeto1i85[i,x,y]
                        maizete2r85[i,x,y]=maizete1r85[i,x,y]
                        maizete2i85[i,x,y]=maizete1i85[i,x,y]



allynewi=N.average(iyield1ynewi,weights=maizeto2,axis=(1,2))
allynnewi=N.average(iyield3ynewi,weights=maizeto2,axis=(1,2))
allyclinewi=N.average(iyield5ynewi,weights=maizeto2,axis=(1,2))

allynewn=N.average(iyield1ynewn,weights=maizeto2,axis=(1,2))
allynnewn=N.average(iyield3ynewn,weights=maizeto2,axis=(1,2))
allyclinewn=N.average(iyield5ynewn,weights=maizeto2,axis=(1,2))

allynewi85=N.average(iyield1ynewi85,weights=maizeto285,axis=(1,2))
allynnewi85=N.average(iyield3ynewi85,weights=maizeto285,axis=(1,2))
allyclinewi85=N.average(iyield5ynewi85,weights=maizeto285,axis=(1,2))

allynewn85=N.average(iyield1ynewn85,weights=maizeto285,axis=(1,2))
allynnewn85=N.average(iyield3ynewn85,weights=maizeto285,axis=(1,2))
allyclinewn85=N.average(iyield5ynewn85,weights=maizeto285,axis=(1,2))





cwn=1-abs((1-allynnewi)-(1-allynnewn))
cwn85=1-abs((1-allynnewi85)-(1-allynnewn85))




cwn1=cwn[1:95]
cwn1_85=cwn85[1:95]


fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(221)

xx=range(2006,2101)
xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
ax.xaxis.set_major_formatter(DateFormatter('%Y'))

ax.plot_date(xdates,cwn,"ro-",label="RCP 4.5",linewidth=2)
ax.plot_date(xdates,cwn85,"bo-",label="RCP 8.5",linewidth=2)

leg=plt.legend(['RCP 4.5','RCP 8.5'],fontsize=18)
leg.get_frame().set_alpha(0.5)


plt.xlabel('Year',fontsize=16)
plt.ylabel('Co-limitation',fontsize=16)
plt.tick_params(axis='both',labelsize=16)
#plt.xlim(0, 2)
#plt.ylim(-3, 1)

ax = fig.add_subplot(222)

ax.xaxis.set_major_formatter(DateFormatter('%Y'))

ax.plot_date(xdates,allynnewi,"ro-",label="RCP 4.5",linewidth=2)
ax.plot_date(xdates,allynnewi85,"bo-",label="RCP 8.5",linewidth=2)

leg=plt.legend(['RCP 4.5','RCP 8.5'],fontsize=18)
leg.get_frame().set_alpha(0.5)


plt.xlabel('Year',fontsize=16)
plt.ylabel('Water limitation',fontsize=16)
plt.tick_params(axis='both',labelsize=16)

ax = fig.add_subplot(223)

ax.xaxis.set_major_formatter(DateFormatter('%Y'))

ax.plot_date(xdates,allynnewn,"ro-",label="RCP 4.5",linewidth=2)
ax.plot_date(xdates,allynnewn85,"bo-",label="RCP 8.5",linewidth=2)

leg=plt.legend(['RCP 4.5','RCP 8.5'],fontsize=18)
leg.get_frame().set_alpha(0.5)


plt.xlabel('Year',fontsize=16)
plt.ylabel('Nitrogen limitation',fontsize=16)
plt.tick_params(axis='both',labelsize=16)


plt.tight_layout()

plt.savefig('maize_8545_inco_paper.png')
plt.show()
