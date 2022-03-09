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

#=== Load data ========================
f1 = xr.open_dataset('autumn_ivt_bc-THR-kernel-t1-s20.nc')
bc = f1['ibt']
# ======== Convert the extracted three-dimensional data of bc into two-dimensional, named bc1
bc1 = np.array(bc).reshape(bc.shape[0],bc.shape[1]*bc.shape[2])
lat, lon = f1['lat'], f1['lon']
#=== ppt data ====
f2 = xr.open_dataset('autumn_merra2_ppt_mm.nc')
ppt = f2['PRECTOT']
ppt1 = np.average(ppt,  axis=1)
ano2 = np.average(ppt1, axis=1)

"""=== Now calculate  correlation ================="""
#cor = np.corrcoef(bc1.T, ano2)[:-1, -1].reshape(len(lat), len(lon))
cor = np.corrcoef(bc1.T, ano2)[:-1, 0].reshape(len(lat), len(lon))

"""=== significance test part, univariate linear regression test f_regression, returns F score, and the p value of F score"""
sig  = f_regression(np.nan_to_num(bc1), ano2)[1].reshape(len(lat), len(lon)) # Use 0 instead of NaN, finite number proxy inf value
area = np.where(sig < 0.05)            #area records the p value on each grid point of latitude and longitude, and the size should be 39*lat*lon


sigma_y = 5
sigma_x = 5
sigma = [sigma_y, sigma_x]

cor1  = sp.ndimage.filters.gaussian_filter(cor, sigma, mode='constant')


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

s = m.pcolormesh(x,y,cor1,cmap=cmaps.BlueWhiteOrangeRed)
#s = m.contour(x,y,cor,cmap=cmaps.MPL_RdBu)


m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)

plt.colorbar(s,orientation='horizontal',pad=0.06,shrink=0.6)       #colorbar
ax.scatter(x[area], y[area], marker='+', s=0.5, c='black', alpha=0.8)  #dot the area that passed the 0.05 significance test
plt.clim(-1,1)

plt.title('Correlation Analysis',pad=0,  fontdict={'family' : 'Times New Roman', 'size'   : 10})

plt.show()

