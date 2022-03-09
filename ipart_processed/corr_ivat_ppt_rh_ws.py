""" correlation anamoly
    of IVT vs. ppt, rh,and ws"""
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
import matplotlib.gridspec as gridspec

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#=== Load data IPART results => IVT  ========================
f1 = xr.open_dataset('all_ivt_bc-THR-kernel-t1-s20.nc')
f2 = xr.open_dataset('all_ivt_oc-THR-kernel-t1-s20.nc')
f3 = xr.open_dataset('all_ivt_dust-THR-kernel-t1-s20.nc')
f4 = xr.open_dataset('all_merra2_ppt_mm.nc')
f5 = xr.open_dataset('all_merra2_wind_rh.nc')

"""=== Get  data ======="""
bc1 = f1['ibt']     #ibt
oc1 = f2['iot']     #iot
dust1 = f3['idt']   #idt
lat, lon = f1['lat'], f1['lon']
#===ppt=============
ppt1 = f4['PRECTOT']
#== specific humidity
rh1  = f5['QV'][:,17,:,:]
#== wind ==========
u   = f5['U'][:,17,:,:]
v   = f5['V'][:,17,:,:]
ws1 = np.sqrt(u**2+v**2) #wind speed

"""== zonal average ============================="""
ppt2 = np.average(ppt1,   axis=1)
ppt  = np.average(ppt2,   axis=1)

rh2 = np.average(rh1,   axis=1)
rh  = np.average(rh2,   axis=1)

ws2 = np.average(ws1,   axis=1)
ws  = np.average(ws2,   axis=1) 

""" ======== Convert the extracted three-dimensional data of bc into two-dimensional"""
bc    = np.array(bc1).reshape(bc1.shape[0],bc1.shape[1]*bc1.shape[2])
oc    = np.array(oc1).reshape(oc1.shape[0],oc1.shape[1]*oc1.shape[2])
dust  = np.array(dust1).reshape(dust1.shape[0],dust1.shape[1]*dust1.shape[2])

"""=== Now calculate  bc and ppt correlation ================="""
#=== IVT vs. ppt
bc_ppt1   = np.corrcoef(bc.T,    ppt)[:-1, -1].reshape(len(lat), len(lon))
oc_ppt1   = np.corrcoef(oc.T,    ppt)[:-1, -1].reshape(len(lat), len(lon))
dust_ppt1 = np.corrcoef(dust.T,  ppt)[:-1, -1].reshape(len(lat), len(lon))
#=== IVT vs. Humidity
#bc_rh1   = np.corrcoef(bc.T,    rh)[:-1, -1].reshape(len(lat), len(lon))
#oc_rh1   = np.corrcoef(oc.T,    rh)[:-1, -1].reshape(len(lat), len(lon))
#dust_rh1 = np.corrcoef(dust.T,  rh)[:-1, -1].reshape(len(lat), len(lon))
#=== IVT vs. wind speed
bc_ws1   = np.corrcoef(bc.T,    ws)[:-1, -1].reshape(len(lat), len(lon))
oc_ws1   = np.corrcoef(oc.T,    ws)[:-1, -1].reshape(len(lat), len(lon))
dust_ws1 = np.corrcoef(dust.T,  ws)[:-1, -1].reshape(len(lat), len(lon))


#==== significant level ppt ===================================================
sig1  = f_regression(np.nan_to_num(bc),   ppt)[1].reshape(len(lat), len(lon))  
sig2  = f_regression(np.nan_to_num(oc),   ppt)[1].reshape(len(lat), len(lon))  
sig3  = f_regression(np.nan_to_num(dust), ppt)[1].reshape(len(lat), len(lon))

#==== significant level rh ===================================================
sig4  = f_regression(np.nan_to_num(bc),   ws)[1].reshape(len(lat), len(lon))
sig5  = f_regression(np.nan_to_num(oc),   ws)[1].reshape(len(lat), len(lon))
sig6  = f_regression(np.nan_to_num(dust), ws)[1].reshape(len(lat), len(lon))


#==== significant level wind speed ===================================================
#sig7  = f_regression(np.nan_to_num(bc),   rh)[1].reshape(len(lat), len(lon))
#sig8  = f_regression(np.nan_to_num(oc),   rh)[1].reshape(len(lat), len(lon))
#sig9  = f_regression(np.nan_to_num(dust), rh)[1].reshape(len(lat), len(lon))

#==== area that passes the significant level ===================
area1 = np.where(sig1 < 0.1)
area2 = np.where(sig2 < 0.1)
area3 = np.where(sig3 < 0.1)
area4 = np.where(sig4 < 0.1)
area5 = np.where(sig5 < 0.1)
area6 = np.where(sig6 < 0.1)
#area7 = np.where(sig7 < 0.1)
#area8 = np.where(sig8 < 0.1)
#area9 = np.where(sig9 < 0.1)



#=== Interpolation =============
sigma_y = 3
sigma_x = 3
sigma = [sigma_y, sigma_x]

#== ivt VS ppt===
bc_ppt    = sp.ndimage.filters.gaussian_filter(bc_ppt1,   sigma, mode='constant')
oc_ppt    = sp.ndimage.filters.gaussian_filter(oc_ppt1,   sigma, mode='constant')
dust_ppt  = sp.ndimage.filters.gaussian_filter(dust_ppt1, sigma, mode='constant')

#== ivt vs RH ==
#bc_rh    = sp.ndimage.filters.gaussian_filter(bc_rh1,   sigma, mode='constant')
#oc_rh    = sp.ndimage.filters.gaussian_filter(oc_rh1,   sigma, mode='constant')
#dust_rh  = sp.ndimage.filters.gaussian_filter(dust_rh1, sigma, mode='constant')

#== ivt vs ws ==
bc_ws    = sp.ndimage.filters.gaussian_filter(bc_ws1,   sigma, mode='constant')
oc_ws    = sp.ndimage.filters.gaussian_filter(oc_ws1,   sigma, mode='constant')
dust_ws  = sp.ndimage.filters.gaussian_filter(dust_ws1, sigma, mode='constant')





""" plot"""
#fig = plt.figure(figsize=(6,3),dpi=300)
#fig, axs = plt.subplots(2, 3,figsize=(6,3),dpi=300)
fig, axs = plt.subplots(2, 3,figsize=(4,2),dpi=300)

gridspec.GridSpec(4,4)


"""#==== BC with ppt  =========================================="""
#ax =plt.subplot(2,3,1)
plt.subplot2grid((2,3), (0,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,bc_ppt,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===== Significant plot =====
plt.scatter(x[area1][::7], y[area1][::7], marker='+', s=0.05, c='black', alpha=0.5)  #dot the area that passed the 0.05 significance test
plt.annotate("[a]", m(129,51),color='black',fontsize=5)


"""#==== OC with ppt =========================================="""
#ax =plt.subplot(2,3,2)
plt.subplot2grid((2,3), (0,1))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,oc_ppt,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
plt.scatter(x[area2][::7], y[area2][::7], marker='+', s=0.05, c='black', alpha=0.5)  #dot the area that passed the 0.05 significance tes
plt.annotate("[b]", m(129,51),color='black',fontsize=5)


"""#==== Dust with ppt =========================================="""
#ax =plt.subplot(2,3,3)
plt.subplot2grid((2,3), (0,2))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,dust_ppt,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
plt.scatter(x[area3][::7], y[area3][::7], marker='+', s=0.05, c='black', alpha=0.5)
plt.annotate("[c]", m(129,51),color='black',fontsize=5)
"""#==== BC with ws  =========================================="""
#ax =plt.subplot(2,3,4)
plt.subplot2grid((2,3), (1,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,bc_ws,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
plt.scatter(x[area4][::7], y[area4][::7], marker='+', s=0.05, c='black', alpha=0.5)

plt.annotate("[d]", m(129,51),color='black',fontsize=5)

"""#==== OC with ws  =========================================="""
#ax =plt.subplot(2,3,5)
plt.subplot2grid((2,3), (1,1))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,oc_ws,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
plt.scatter(x[area5][::7], y[area5][::7], marker='+', s=0.05, c='black', alpha=0.5)
plt.annotate("[e]", m(129,51),color='black',fontsize=5)
"""==== DUST with ws  =========================================="""
#ax =plt.subplot(2,3,6)
plt.subplot2grid((2,3), (1,2))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.pcolormesh(x,y,dust_ws,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(-1,1)
#===Significant plot =====
plt.scatter(x[area6][::7], y[area6][::7], marker='+', s=0.05, c='black', alpha=0.5)
plt.annotate("[f]", m(129,51),color='black',fontsize=5)

#===== adjust colorbar=================
cax = fig.add_axes([0.035,0.085,0.947,0.026])   # left, bottom, width, and height
cb  = fig.colorbar(ax,  cax=cax, orientation='horizontal',ticks=[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
cb.set_label(label='Correlation Coefficient ',size=5,labelpad=0)
cb.ax.tick_params(labelsize=5)



#====== adjust =====================
fig.subplots_adjust(top=1,
                        bottom=0.155,
                        left=0.035,
                        right=0.98,
                        hspace=0.0,
                        wspace=0.095)

#plt.savefig("/mnt/g/2nd_Paper/ipart/corr_ws_ppt_2000.png",dpi=2000)
plt.show()


