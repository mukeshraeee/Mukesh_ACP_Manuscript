"""" correlation anamoly"""
import matplotlib
import numpy as np
import xarray as xr
from sklearn.feature_selection import f_regression
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
import pandas as pd
import xarray as xr
from mpl_toolkits.basemap import Basemap
import cmaps
import scipy as sp
import scipy.ndimage

#=== Load data IPART results => IVT  ========================
#== bc ====
f1 = xr.open_dataset('winter_ivt_bc-THR-kernel-t1-s20.nc')
f2 = xr.open_dataset('spring_ivt_bc-THR-kernel-t1-s20.nc')
f3 = xr.open_dataset('monsoon_ivt_bc-THR-kernel-t1-s20.nc')
f4 = xr.open_dataset('autumn_ivt_bc-THR-kernel-t1-s20.nc')
#== oc ====
f5 = xr.open_dataset('winter_ivt_oc-THR-kernel-t1-s20.nc')
f6 = xr.open_dataset('spring_ivt_oc-THR-kernel-t1-s20.nc')
f7 = xr.open_dataset('monsoon_ivt_oc-THR-kernel-t1-s20.nc')
f8 = xr.open_dataset('autumn_ivt_oc-THR-kernel-t1-s20.nc')
#== dust ====
f9 = xr.open_dataset('winter_ivt_dust-THR-kernel-t1-s20.nc')
f10 = xr.open_dataset('spring_ivt_dust-THR-kernel-t1-s20.nc')
f11 = xr.open_dataset('monsoon_ivt_dust-THR-kernel-t1-s20.nc')
f12 = xr.open_dataset('autumn_ivt_dust-THR-kernel-t1-s20.nc')

"""=== Get BC data ======="""
bc1 = f1['ibt']
bc2 = f2['ibt']
bc3 = f3['ibt']
bc4 = f4['ibt']
#== ======== Convert the extracted three-dimensional data of bc into two-dimensional, named bc1
bc5 = np.array(bc1).reshape(bc1.shape[0],bc1.shape[1]*bc1.shape[2]) #winter
bc6 = np.array(bc2).reshape(bc2.shape[0],bc2.shape[1]*bc2.shape[2]) #spring
bc7 = np.array(bc3).reshape(bc3.shape[0],bc3.shape[1]*bc3.shape[2]) #monsoon
bc8 = np.array(bc4).reshape(bc4.shape[0],bc4.shape[1]*bc4.shape[2]) #autumn
#==== annual ==================
bc = (bc5+bc6+bc7+bc8)/4

"""=== Get OC data ======="""
oc1 = f5['iot']
oc2 = f6['iot']
oc3 = f7['iot']
oc4 = f8['iot']
#== ======== Convert the extracted three-dimensional data of bc into two-dimensional, named bc1
oc5 = np.array(oc1).reshape(oc1.shape[0],oc1.shape[1]*oc1.shape[2])
oc6 = np.array(oc2).reshape(oc2.shape[0],oc2.shape[1]*oc2.shape[2])
oc7 = np.array(oc3).reshape(oc3.shape[0],oc3.shape[1]*oc3.shape[2])
oc8 = np.array(oc4).reshape(oc4.shape[0],oc4.shape[1]*oc4.shape[2])
#==== annual ==================
oc = (oc5+oc6+oc7+oc8)/4

"""=== Get dust data ======="""
d1 = f9['idt']
d2 = f10['idt']
d3 = f11['idt']
d4 = f12['idt']
#== ======== Convert the extracted three-dimensional data of bc into two-dimensional, named bc1
d5 = np.array(d1).reshape(d1.shape[0],d1.shape[1]*d1.shape[2])
d6 = np.array(d2).reshape(d2.shape[0],d2.shape[1]*d2.shape[2])
d7 = np.array(d3).reshape(d3.shape[0],d3.shape[1]*d3.shape[2])
d8 = np.array(d4).reshape(d4.shape[0],d4.shape[1]*d4.shape[2])
#==== annual ==================
dust = (d5+d6+d7+d8)/4


#==== get lat lon =====================
lat, lon = f1['lat'], f1['lon']



"""=== ppt data ===="""
p1 = xr.open_dataset('winter_merra2_ppt_mm.nc')
p2 = xr.open_dataset('spring_merra2_ppt_mm.nc')
p3 = xr.open_dataset('monsoon_merra2_ppt_mm.nc')
p4 = xr.open_dataset('autumn_merra2_ppt_mm.nc')

ppt1 = p1['PRECTOT']
ppt2 = p2['PRECTOT']
ppt3 = p3['PRECTOT']
ppt4 = p4['PRECTOT']

#==winter
win_p1 = np.average(ppt1,   axis=1)
win_p  = np.average(win_p1, axis=1)

#==spring
spr_p1 = np.average(ppt2,   axis=1)
spr_p  = np.average(spr_p1, axis=1)

#==monsoon
mon_p1 = np.average(ppt3,   axis=1)
mon_p  = np.average(mon_p1, axis=1)

#==autumn
aut_p1 = np.average(ppt4,   axis=1)
aut_p  = np.average(aut_p1, axis=1)


"""=== Now calculate  bc and ppt correlation ================="""
cor1 = np.corrcoef(bc5.T, win_p)[:-1, -1].reshape(len(lat), len(lon)) #winter
cor2 = np.corrcoef(bc6.T, spr_p)[:-1, -1].reshape(len(lat), len(lon)) #spring
cor3 = np.corrcoef(bc7.T, mon_p)[:-1, -1].reshape(len(lat), len(lon)) #monsonn
cor4 = np.corrcoef(bc8.T, aut_p)[:-1, -1].reshape(len(lat), len(lon)) #autumn


bc_cor = (cor1+cor2+cor3+cor4)/4


"""significance test part, univariate linear regression test f_regression, returns F score, and the p value of F score"""
sig1  = f_regression(np.nan_to_num(bc5), win_p)[1].reshape(len(lat), len(lon)) 
sig2  = f_regression(np.nan_to_num(bc6), spr_p)[1].reshape(len(lat), len(lon))
sig3  = f_regression(np.nan_to_num(bc7), mon_p)[1].reshape(len(lat), len(lon))
sig4  = f_regression(np.nan_to_num(bc8), aut_p)[1].reshape(len(lat), len(lon))



area1 = np.where(sig1 < 0.1)            
area2 = np.where(sig2 < 0.1)
area3 = np.where(sig3 < 0.1)
area4 = np.where(sig4 < 0.1)

#area = (area1+area2+area3+area4)/4


#=== Interpolation =============
sigma_y = 10
sigma_x = 10
sigma = [sigma_y, sigma_x]
bc_cor  = sp.ndimage.filters.gaussian_filter(bc_cor, sigma, mode='constant')


#===== plot =================================

fig = plt.figure(figsize=(4,3),dpi=300)

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#=== Precipitation  ============================
""" winter"""
ax =plt.subplot(1,1,1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)

s = m.pcolormesh(x,y,bc_cor,cmap=cmaps.BlueWhiteOrangeRed)
#s = m.contour(x,y,cor,cmap=cmaps.MPL_RdBu)


m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)

plt.colorbar(s,orientation='horizontal',pad=0.06,shrink=0.6)       #colorbar
ax.scatter(x[area3], y[area3], marker='+', s=0.5, c='black', alpha=0.8)  #dot the area that passed the 0.05 significance test
plt.clim(-1,1)

plt.title('Correlation Analysis',pad=0,  fontdict={'family' : 'Times New Roman', 'size'   : 10})

plt.show()

