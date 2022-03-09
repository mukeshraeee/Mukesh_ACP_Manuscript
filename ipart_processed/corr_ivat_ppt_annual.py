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
from matplotlib import rcParams

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#=== Load data IPART results => IVT  ========================
#== bc ====
f1 = xr.open_dataset('all_ivt_bc-THR-kernel-t1-s20.nc')
f2 = xr.open_dataset('all_ivt_oc-THR-kernel-t1-s20.nc')
f3 = xr.open_dataset('all_ivt_dust-THR-kernel-t1-s20.nc')
f4 = xr.open_dataset('all_merra2_ppt_mm.nc')

"""=== Get  data ======="""
bc1 = f1['ivt_ano']     #ibt
oc1 = f2['ivt_ano']     #iot
dust1 = f3['ivt_ano']   #idt
ppt1 = f4['PRECTOT']
lat, lon = f1['lat'], f1['lon']
#== zonal average =============================
ppt2 = np.average(ppt1,   axis=1)
ppt  = np.average(ppt2,   axis=1)


""" ======== Convert the extracted three-dimensional data of bc into two-dimensional"""
bc    = np.array(bc1).reshape(bc1.shape[0],bc1.shape[1]*bc1.shape[2])
oc    = np.array(oc1).reshape(oc1.shape[0],oc1.shape[1]*oc1.shape[2])
dust  = np.array(dust1).reshape(dust1.shape[0],dust1.shape[1]*dust1.shape[2])

"""=== Now calculate  bc and ppt correlation ================="""
cor1 = np.corrcoef(bc.T,   ppt)[:-1, -1].reshape(len(lat), len(lon)) # BC
cor2 = np.corrcoef(oc.T,   ppt)[:-1, -1].reshape(len(lat), len(lon)) # OC
cor3 = np.corrcoef(dust.T, ppt)[:-1, -1].reshape(len(lat), len(lon)) # DUST

"""significance test part, univariate linear regression test f_regression, returns F score, and the p value of F score"""
sig1  = f_regression(np.nan_to_num(bc),   ppt)[1].reshape(len(lat), len(lon))  # BC
sig2  = f_regression(np.nan_to_num(oc),   ppt)[1].reshape(len(lat), len(lon))  # OC
sig3  = f_regression(np.nan_to_num(dust), ppt)[1].reshape(len(lat), len(lon))  # DUST

#==== area===================
area1 = np.where(sig1 < 0.05) #BC
area2 = np.where(sig2 < 0.05) #OC
area3 = np.where(sig3 < 0.05) #DUST

#=== Interpolation =============
sigma_y = 3
sigma_x = 3
sigma = [sigma_y, sigma_x]

bc_cor   = sp.ndimage.filters.gaussian_filter(cor1,   sigma, mode='constant')
oc_cor   = sp.ndimage.filters.gaussian_filter(cor2,   sigma, mode='constant')
dust_cor = sp.ndimage.filters.gaussian_filter(cor3,   sigma, mode='constant')


""" plot"""
fig = plt.figure(figsize=(4.5,1.5),dpi=400)


"""#==== BC =========================================="""
ax =plt.subplot(1,3,1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,bc_cor,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=4)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,1,0],color='black',fontsize=4)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===== Significant plot =====
ax.scatter(x[area1][::10], y[area1][::10], marker='+', s=0.05, c='black', alpha=0.6)  #dot the area that passed the 0.05 significance test
#==colorbar
cb = plt.colorbar(s,orientation='horizontal',pad=0.02,shrink=1)
cb.ax.tick_params(labelsize=3.5)
cb.set_label(r'BC_IVT',fontsize=4,labelpad=2)
"""#==== OC =========================================="""
ax =plt.subplot(1,3,2)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,oc_cor,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,1,0],color='black',fontsize=4)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
ax.scatter(x[area2][::10], y[area2][::10], marker='+', s=0.05, c='black', alpha=0.6)  #dot the area that passed the 0.05 significance test
#== colorbar
cb = plt.colorbar(s,orientation='horizontal',pad=0.02,shrink=1)
cb.ax.tick_params(labelsize=3.5)
cb.set_label(r'OC_IVT',fontsize=4,labelpad=2)
"""#==== dust =========================================="""
ax =plt.subplot(1,3,3)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,dust_cor,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,1,0],color='black',fontsize=4)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
ax.scatter(x[area3][::10], y[area3][::10], marker='+', s=0.05, c='black', alpha=0.6)  #dot the area that passed the 0.05 significance test
#== colorbar
cb = plt.colorbar(s,orientation='horizontal',pad=0.02,shrink=1)
cb.ax.tick_params(labelsize=3.5)
cb.set_label(r'Dust_IVT',fontsize=4,labelpad=2)
#====== adjust =====================
fig.subplots_adjust(top=0.96,
                        bottom=0.1,
                        left=0.035,
                        right=0.99,
                        hspace=0.19,
                        wspace=0.09)

#plt.savefig("/mnt/g/2nd_Paper/ipart/ipart_ano_corr_30000.png",dpi=3000)
plt.show()
