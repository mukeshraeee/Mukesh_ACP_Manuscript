"""This code generate pm, oc and dust plot for seasonal basis from
	MERRA-2 and WRF-Chem models, ~ Mukesh Rai ~ 2021-09-13"""

#==== Import Library ==============================
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import wrf
from wrf import (to_np, getvar,interplevel, smooth2d, get_basemap, latlon_coords, ALL_TIMES)
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import cmaps
from pylab import *
import scipy as sp
import scipy.ndimage
import matplotlib.gridspec as gridspec
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from matplotlib.patches import FancyArrowPatch
from matplotlib.legend_handler import HandlerLine2D
from matplotlib import rcParams
#======== read wrf data====
win = [Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-01-01_00_00_00"),
      Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-02-01_00_00_00"),
      Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-12-01_00_00_00")]
spr = [Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-03-01_00_00_00.nc"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-04-01_00_00_00.nc"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-05-01_00_00_00.nc")]
mon = [Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-06-01_00_00_00"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-07-01_00_00_00"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-08-01_00_00_00")]
aut = [Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-09-01_00_00_00"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-10-01_00_00_00"),
       Dataset("/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/wrfout_d01_2017-11-01_00_00_00")]

#====== read merra-2 data ======
win1 = Dataset("/mnt/g/2nd_Paper/merra2_bc_oc_dust/winter1.nc")
spr1 = Dataset("/mnt/g/2nd_Paper/merra2_bc_oc_dust/spring1.nc")
mon1 = Dataset("/mnt/g/2nd_Paper/merra2_bc_oc_dust/monsoon1.nc")
aut1 = Dataset("/mnt/g/2nd_Paper/merra2_bc_oc_dust/autumn1.nc")

#====== read merra-2 GPH data ======
win_gph1 = Dataset("/mnt/g/2nd_Paper/merra_2_Geopotential_Height/winter.nc")
spr_gph1 = Dataset("/mnt/g/2nd_Paper/merra_2_Geopotential_Height/spring.nc")
mon_gph1 = Dataset("/mnt/g/2nd_Paper/merra_2_Geopotential_Height/summer.nc")
aut_gph1 = Dataset("/mnt/g/2nd_Paper/merra_2_Geopotential_Height/autumn.nc")

#==== get Geo-Potential Height =============
win_merra_z1  = win_gph1.variables['H'][:,7,:,:]
spr_merra_z1  = spr_gph1.variables['H'][:,7,:,:]
mon_merra_z1  = mon_gph1.variables['H'][:,7,:,:]
aut_merra_z1  = aut_gph1.variables['H'][:,7,:,:]

win_merra_z = np.mean(win_merra_z1,axis=0)
spr_merra_z = np.mean(spr_merra_z1,axis=0)
mon_merra_z = np.mean(mon_merra_z1,axis=0)
aut_merra_z = np.mean(aut_merra_z1,axis=0)
all_merra_z = (win_merra_z+spr_merra_z+mon_merra_z+aut_merra_z)/4
#===== get OC,PM2.5 and PM10 data from WRF data =====================
#=== pm2.5==
win_pm25 = getvar(win, "PM2_5_DRY", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_pm25 = getvar(spr, "PM2_5_DRY", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_pm25 = getvar(mon, "PM2_5_DRY", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_pm25 = getvar(aut, "PM2_5_DRY", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

#=== pm10===
win_pm10 = getvar(win, "PM10", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_pm10 = getvar(spr, "PM10", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_pm10 = getvar(mon, "PM10", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_pm10 = getvar(aut, "PM10", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

#===== oc==
win_oc1 = getvar(win, "OC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_oc1 = getvar(spr, "OC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_oc1 = getvar(mon, "OC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_oc1 = getvar(aut, "OC1", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

win_oc2 = getvar(win, "OC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
spr_oc2 = getvar(spr, "OC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
mon_oc2 = getvar(mon, "OC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]
aut_oc2 = getvar(aut, "OC2", timeidx=ALL_TIMES,method='cat')[:,0,:,:]

win_oc3 = win_oc1+win_oc2
spr_oc3 = spr_oc1+spr_oc2
mon_oc3 = mon_oc1+mon_oc2
aut_oc3 = aut_oc1+aut_oc2


#=== Now take time average ============
win_p25 = np.average(win_pm25,axis=0)
spr_p25 = np.average(spr_pm25,axis=0)
mon_p25 = np.average(mon_pm25,axis=0)
aut_p25 = np.average(aut_pm25,axis=0)

total_pm25 = (win_p25+spr_p25+mon_p25+aut_p25)/4



"""=== Now get geopotential height=="""
p  = getvar(win, "pressure")

z1 = getvar(win, "z", units="m",timeidx=-1,method='cat') #winter
z2 = getvar(spr, "z", units="m",timeidx=-1,method='cat') #spring
z3 = getvar(mon, "z", units="m",timeidx=-1,method='cat') #monsoon
z4 = getvar(aut, "z", units="m",timeidx=-1,method='cat') #autumn

win_z1  = interplevel(z1, p, 850) # Interpolate geopotential height at 850 hPa
spr_z1  = interplevel(z2, p, 850)
mon_z1  = interplevel(z3, p, 850)
aut_z1  = interplevel(z4, p, 850)
all_z1  = (win_z1+spr_z1+mon_z1+aut_z1)/4

#======== Smoothing pm2.5  using Gaussion filter =============
sigma_y = 2
sigma_x = 2
sigma = [sigma_y, sigma_x]

winter_p25  = sp.ndimage.filters.gaussian_filter(win_p25, sigma, mode='constant')
spring_p25  = sp.ndimage.filters.gaussian_filter(spr_p25, sigma, mode='constant')
monsoon_p25 = sp.ndimage.filters.gaussian_filter(mon_p25, sigma, mode='constant')
autumn_p25  = sp.ndimage.filters.gaussian_filter(aut_p25, sigma, mode='constant')
total_p25   = sp.ndimage.filters.gaussian_filter(total_pm25,sigma, mode='constant')

#======== pm10 =============
win_p10 = np.average(win_pm10,axis=0)
spr_p10 = np.average(spr_pm10,axis=0)
mon_p10 = np.average(mon_pm10,axis=0)
aut_p10 = np.average(aut_pm10,axis=0)

total_pm10 = (win_p10+spr_p10+mon_p10+aut_p10)/4

#======== Smoothing pm10  using Gaussion filter =============
winter_p10  = sp.ndimage.filters.gaussian_filter(win_p10, sigma, mode='constant')
spring_p10  = sp.ndimage.filters.gaussian_filter(spr_p10, sigma, mode='constant')
monsoon_p10 = sp.ndimage.filters.gaussian_filter(mon_p10, sigma, mode='constant')
autumn_p10  = sp.ndimage.filters.gaussian_filter(aut_p10, sigma, mode='constant')
total_p10   = sp.ndimage.filters.gaussian_filter(total_pm10,sigma, mode='constant')



win_oc = np.average(win_oc3, axis=0)
spr_oc = np.average(spr_oc3, axis=0)
mon_oc = np.average(mon_oc3, axis=0)
aut_oc = np.average(aut_oc3, axis=0)

total_oc = (win_oc+spr_oc+mon_oc+aut_oc)/4

#======== Smoothing oc  using Gaussion filter =============
winter_oc = sp.ndimage.filters.gaussian_filter(win_oc, sigma, mode='constant')
spring_oc = sp.ndimage.filters.gaussian_filter(spr_oc, sigma, mode='constant')
monsoon_oc = sp.ndimage.filters.gaussian_filter(mon_oc, sigma, mode='constant')
autumn_oc = sp.ndimage.filters.gaussian_filter(aut_oc, sigma, mode='constant')
total_oc  = sp.ndimage.filters.gaussian_filter(total_oc,sigma, mode='constant')

#===== get dust data from merra-2===
win_dust1 = win1.variables['DUSMASS'][:]
spr_dust1 = spr1.variables['DUSMASS'][:]
mon_dust1 = mon1.variables['DUSMASS'][:]
aut_dust1 = aut1.variables['DUSMASS'][:]

win_dust = np.mean(win_dust1[:,:,:],axis=0)*1000000000
spr_dust = np.mean(spr_dust1[:,:,:],axis=0)*1000000000
mon_dust = np.mean(mon_dust1[:,:,:],axis=0)*1000000000
aut_dust = np.mean(aut_dust1[:,:,:],axis=0)*1000000000

#==== convert kg to ug =====
total_dust = (win_dust+spr_dust+mon_dust+aut_dust)/4


""" Now plot PM2.5,PM10,OC, and DUST data """

fig = plt.figure(figsize=(12,10),dpi=300)

lats, lons = latlon_coords(win_pm25)

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#=== plot pm2.5============================
""" winter"""
ax =plt.subplot(5,4,1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,winter_p25,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.ylabel('Winter', labelpad=5,fontsize=12)
plt.clim(0,130)

#=== plot geopotential height==
contours = plt.contour(x, y, win_z1, 7, linestyles="dashed",linewidths=0.8,colors='red')
plt.clabel(contours, inline=True, fontsize=7,inline_spacing=0,fmt='%1.0f')


""" spring"""
ax =plt.subplot(5,4,5)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,spring_p25,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.ylabel('Spring', labelpad=5,fontsize=12)
plt.clim(0,130)

#=== plot geopotential height==
contours = plt.contour(x, y, spr_z1, 7,linestyles="dashed", linewidths=0.8,colors='dimgray')
plt.clabel(contours, inline=True, fontsize=7,inline_spacing=0,fmt='%1.0f')


""" monsoon """
ax =plt.subplot(5,4,9)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,monsoon_p25,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.ylabel('Summer', labelpad=5,fontsize=12)
plt.clim(0,130)

#=== plot geopotential height==
contours = plt.contour(x, y, mon_z1, 7, linestyles="dashed",linewidths=0.8,colors='tab:orange')
plt.clabel(contours, inline=True, fontsize=7,inline_spacing=0,fmt='%1.0f')

""" autumn """
ax =plt.subplot(5,4,13)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,autumn_p25,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.ylabel('Autumn', labelpad=5,fontsize=12)
plt.clim(0,130)

#=== plot geopotential height==
contours = plt.contour(x, y, aut_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')

""" total """
ax1 =plt.subplot(5,4,17)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s1 = m.pcolormesh(x,y,total_p25,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],fontsize=12,color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Annual', labelpad=5,fontsize=12)
plt.clim(0,130)
#=== plot geopotential height==
contours = plt.contour(x, y, all_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')



#=== Setting Colorbar for PM2.5 =========
cax1 = fig.add_axes([0.04,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s1, ax=ax1, cax=cax1, ticks=[0,15,30,45,60,75,90,105,120,132], orientation='horizontal')
cb1.set_label(r'PM2.5 [$μgm^{-3}$]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)

"""For PM10 """

""" winter"""
ax =plt.subplot(5,4,2)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,winter_p10,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,200)
#=== plot geopotential height==
contours = plt.contour(x, y, win_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" spring"""
ax =plt.subplot(5,4,6)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,monsoon_p10,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,200)
#=== plot geopotential height==
contours = plt.contour(x, y, spr_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" monsoon """
ax =plt.subplot(5,4,10)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,monsoon_p10,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,200)
#=== plot geopotential height==
contours = plt.contour(x, y, mon_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" autumn """
ax =plt.subplot(5,4,14)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,autumn_p10,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.3)
m.drawcoastlines(linewidth=0.3)
plt.clim(0,200)
#=== plot geopotential height==
contours = plt.contour(x, y, aut_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')



""" total """
ax2 =plt.subplot(5,4,18)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s2 = m.pcolormesh(x,y,total_p10,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, labels=[0,0,0,1],dashes=[2, 2],  fontsize=12,color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,200)
#=== plot geopotential height==
contours = plt.contour(x, y, all_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


#=== Setting Colorbar for PM10 =========
cax1 = fig.add_axes([0.285,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s2, ax=ax2, cax=cax1, ticks=[0,25,50,75,100,125,150,175,200], orientation='horizontal')
cb1.set_label(r'PM10 [$μgm^{-3}$]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)


""" plot OC """

""" winter"""
ax =plt.subplot(5,4,3)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,winter_oc,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,30)
#=== plot geopotential height==
contours = plt.contour(x, y, win_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" spring"""
ax =plt.subplot(5,4,7)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,spring_oc,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,30)
#=== plot geopotential height==
contours = plt.contour(x, y, spr_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" monsoon """
ax =plt.subplot(5,4,11)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,monsoon_oc,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,30)
#=== plot geopotential height==
contours = plt.contour(x, y, mon_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


""" autumn """
ax =plt.subplot(5,4,15)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s = m.pcolormesh(x,y,autumn_oc,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,30)
#=== plot geopotential height==
contours = plt.contour(x, y, aut_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')



""" total """
ax3 =plt.subplot(5,4,19)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lons), to_np(lats))
s3 = m.pcolormesh(x,y,total_oc,cmap=cmaps.temp_diff_18lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, labels=[0,0,0,1],dashes=[2, 2],  fontsize=12,color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,30)
#=== plot geopotential height==
contours = plt.contour(x, y, all_z1, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


#=== Setting Colorbar for OC =========
cax1 = fig.add_axes([0.528,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s3, ax=ax3, cax=cax1, ticks=[0,5,10,15,20,25,30], orientation='horizontal')
cb1.set_label(r'OC [$μgm^{-3}$]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)



""" plot DUST """
lat1   = win1.variables['lat'][:]
lon1   = win1.variables['lon'][:]

"""winter"""
ax =plt.subplot(5,4,4)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
lons,lats = np.meshgrid(lon1,lat1)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,win_dust,cmap=cmaps.temp_diff_18lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,600)
#=== plot geopotential height==
contours = plt.contour(x, y, win_merra_z, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


"""spring"""
ax =plt.subplot(5,4,8)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
lons,lats = np.meshgrid(lon1,lat1)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,spr_dust,cmap=cmaps.temp_diff_18lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,600)
#=== plot geopotential height==
contours = plt.contour(x, y, spr_merra_z, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


"""monsoon"""
ax =plt.subplot(5,4,12)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
lons,lats = np.meshgrid(lon1,lat1)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,mon_dust,cmap=cmaps.temp_diff_18lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,600)
#=== plot geopotential height==
contours = plt.contour(x, y, mon_merra_z, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


"""autumn"""
ax =plt.subplot(5,4,16)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
lons,lats = np.meshgrid(lon1,lat1)
x,y = m(lons, lats)
s = m.pcolormesh(x,y,aut_dust,cmap=cmaps.temp_diff_18lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,600)
#=== plot geopotential height==
contours = plt.contour(x, y, aut_merra_z, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


"""total"""
ax4 =plt.subplot(5,4,20)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
lons,lats = np.meshgrid(lon1,lat1)
x,y = m(lons, lats)
s4 = m.pcolormesh(x,y,total_dust,cmap=cmaps.temp_diff_18lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, labels=[0,0,0,1],dashes=[2, 2],  fontsize=12,color='black')
m.drawcountries(linewidth=0.5)
m.drawcoastlines(linewidth=0.5)
plt.clim(0,600)
#=== plot geopotential height==
contours = plt.contour(x, y, all_merra_z, 10, linestyles="dashed",linewidths=0.3)
plt.clabel(contours, inline=True, fontsize=6,inline_spacing=0,fmt='%1.0f')


#=== Setting Colorbar for dust =========
cax1 = fig.add_axes([0.775,0.1,0.185,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s4, ax=ax4, cax=cax1, ticks=[0,100,200,300,400,500,600], orientation='horizontal')
cb1.set_label(r'Dust [$μgm^{-3}$]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)



""" Adjust"""
fig.subplots_adjust(top=0.96,
                        bottom=0.14,
                        left=0.015,
                        right=0.99,
                        hspace=0.05,
                        wspace=0.0)

plt.savefig("/mnt/g/2nd_Paper/xxx.png",dpi=600)
#plt.show()
