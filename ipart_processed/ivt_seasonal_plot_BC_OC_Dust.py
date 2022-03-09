""" This code plot IVT for BC,OC, and Dust that generated from IPART mod
	- Mukesh Rai, 2021/10/20"""

# ====== Import libraries ==================
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
#== bc
f1 = xr.open_dataset('winter_ivt_bc-THR-kernel-t1-s20.nc')
f2 = xr.open_dataset('spring_ivt_bc-THR-kernel-t1-s20.nc')
f3 = xr.open_dataset('monsoon_ivt_bc-THR-kernel-t1-s20.nc')
f4 = xr.open_dataset('autumn_ivt_bc-THR-kernel-t1-s20.nc')
#==oc
f5 = xr.open_dataset('winter_ivt_oc-THR-kernel-t1-s20.nc')
f6 = xr.open_dataset('spring_ivt_oc-THR-kernel-t1-s20.nc')
f7 = xr.open_dataset('monsoon_ivt_oc-THR-kernel-t1-s20.nc')
f8 = xr.open_dataset('autumn_ivt_oc-THR-kernel-t1-s20.nc')
#==dust
f9  = xr.open_dataset('winter_ivt_dust-THR-kernel-t1-s20.nc')
f10 = xr.open_dataset('spring_ivt_dust-THR-kernel-t1-s20.nc')
f11 = xr.open_dataset('monsoon_ivt_dust-THR-kernel-t1-s20.nc')
f12 = xr.open_dataset('autumn_ivt_dust-THR-kernel-t1-s20.nc')

#=== Get lat and lon===============
lat, lon = f1['lat'], f1['lon']
#== Get variables ================
bc_win1 = f1['ivt_ano']
bc_spr1 = f2['ivt_ano']
bc_mon1 = f3['ivt_ano']
bc_aut1 = f4['ivt_ano']

oc_win1 = f5['ivt_ano']
oc_spr1 = f6['ivt_ano']
oc_mon1 = f7['ivt_ano']
oc_aut1 = f8['ivt_ano']

dust_win1 = f9['ivt_ano']
dust_spr1 = f10['ivt_ano']
dust_mon1 = f11['ivt_ano']
dust_aut1= f12['ivt_ano']

#=== take time average ==========
bc_win = np.average(bc_win1,axis=0)*100000
bc_spr = np.average(bc_spr1,axis=0)*100000
bc_mon = np.average(bc_mon1,axis=0)*100000
bc_aut = np.average(bc_aut1,axis=0)*100000

oc_win = np.average(oc_win1,axis=0)*100000
oc_spr = np.average(oc_spr1,axis=0)*100000
oc_mon = np.average(oc_mon1,axis=0)*100000
oc_aut = np.average(oc_aut1,axis=0)*100000

dust_win = np.average(dust_win1,axis=0)*100000
dust_spr = np.average(dust_spr1,axis=0)*100000
dust_mon = np.average(dust_mon1,axis=0)*100000
dust_aut = np.average(dust_aut1,axis=0)*100000


#=== interpolation ==========
#sigma_y = 3
#sigma_x = 3
#sigma = [sigma_y, sigma_x]

#bc_win    = sp.ndimage.filters.gaussian_filter(bc_win2, sigma, mode='constant')
#bc_spr    = sp.ndimage.filters.gaussian_filter(bc_spr2, sigma, mode='constant')
#bc_mon    = sp.ndimage.filters.gaussian_filter(bc_mon2, sigma, mode='constant')
#bc_aut    = sp.ndimage.filters.gaussian_filter(bc_aut2, sigma, mode='constant')

#oc_win    = sp.ndimage.filters.gaussian_filter(oc_win2, sigma, mode='constant')
#oc_spr    = sp.ndimage.filters.gaussian_filter(oc_spr2, sigma, mode='constant')
#oc_mon    = sp.ndimage.filters.gaussian_filter(oc_mon2, sigma, mode='constant')
#oc_aut    = sp.ndimage.filters.gaussian_filter(oc_aut2, sigma, mode='constant')

#dust_win    = sp.ndimage.filters.gaussian_filter(dust_win2, sigma, mode='constant')
#dust_spr    = sp.ndimage.filters.gaussian_filter(dust_spr2, sigma, mode='constant')
#dust_mon    = sp.ndimage.filters.gaussian_filter(dust_mon2, sigma, mode='constant')
#dust_aut    = sp.ndimage.filters.gaussian_filter(dust_aut2, sigma, mode='constant')


""" plot"""
fig, axs = plt.subplots(3, 4,figsize=(6,2.5),dpi=300)
gridspec.GridSpec(3,4)

#= BC-winter =
plt.subplot2grid((3,4), (0,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.contourf(x,y,bc_win,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,3)
plt.title('Winter', fontsize=5,y=0.9)


#= BC-spring =
plt.subplot2grid((3,4), (0,1))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.contourf(x,y,bc_spr,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Spring', fontsize=5,y=0.9)
plt.clim(0,3)

#= BC-monsoon =
plt.subplot2grid((3,4), (0,2))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax = m.contourf(x,y,bc_mon,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Summer', fontsize=5,y=0.9)
plt.clim(0,3)

#= BC-autumn =
plt.subplot2grid((3,4), (0,3))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s = m.contourf(x,y,bc_aut,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.title('Autumn', fontsize=5,y=0.9)
plt.clim(0,3)

#==== set colorbar ======
b1 =fig.colorbar(s,ticks=[0,0.5,1,1.5,2,2.5,3])
b1.set_label('BC_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]', rotation=-270, fontsize=4,labelpad=4)
b1.ax.tick_params(labelsize=4)





#= 0C-winter =
plt.subplot2grid((3,4), (1,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax2 = m.contourf(x,y,oc_win,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,12)


#= OC-spring =
plt.subplot2grid((3,4), (1,1))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax2 = m.contourf(x,y,oc_spr,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,12)

#= OC-monsoon =
plt.subplot2grid((3,4), (1,2))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax2 = m.contourf(x,y,oc_mon,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,12)

#= OC-autumn =
plt.subplot2grid((3,4), (1,3))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s2 = m.contourf(x,y,oc_aut,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,12)

b2 =fig.colorbar(s2,ticks=[0,2,4,6,8,10,12])
b2.set_label('OC_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]', rotation=-270, fontsize=4,labelpad=4)
b2.ax.tick_params(labelsize=4)

#= Dust-winter =
plt.subplot2grid((3,4), (2,0))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax3 = m.contourf(x,y,dust_win,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[1,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,300)


#= Dust-spring =
plt.subplot2grid((3,4), (2,1))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax3 = m.contourf(x,y,dust_spr,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,300)

#= Dust-monsoon =
plt.subplot2grid((3,4), (2,2))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
ax3 = m.contourf(x,y,dust_mon,500,cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,300)

#= Dust-autumn =
plt.subplot2grid((3,4), (2,3))
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')

lons,lats = np.meshgrid(lon,lat)
x,y = m(lons, lats)
s3 = m.contourf(x,y,dust_aut,500, cmap=cmaps.BlueWhiteOrangeRed)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,0],color='black',fontsize=5)
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black',fontsize=5)
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,300)

b3 =fig.colorbar(s3,ticks=[0,50,100,150,200,250,300])
b3.set_label('Dust_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]', rotation=-270, fontsize=4,labelpad=4)
b3.ax.tick_params(labelsize=4)


fig.subplots_adjust(top=0.97,
                        bottom=0.04,
                        left=0.01,
                        right=0.955,
                        hspace=0.06,
                        wspace=0.0)



plt.savefig("/mnt/g/2nd_Paper/ipart/ivt_bc_oc_dust_1200.png",dpi=1200)
#plt.show()
