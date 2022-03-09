""" This codes produce the meteorology parameters including 
    Precipitation, Wind, PBLH and RH """

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import cmaps
from matplotlib.patches import Patch
from pylab import *
from wrf import (to_np, getvar, interplevel, smooth2d, get_basemap, latlon_coords, ALL_TIMES)
from matplotlib.pyplot import figure
import matplotlib.gridspec as gridspec
import wrf
import scipy as sp
import scipy.ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
import cmaps
import matplotlib.patches as mpatches
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


"""=== PPT =="""
win_ppt1 = getvar(win, "RAINC", timeidx=ALL_TIMES,method='cat')[:,:,:]
spr_ppt1 = getvar(spr, "RAINC", timeidx=ALL_TIMES,method='cat')[:,:,:]
mon_ppt1 = getvar(mon, "RAINC", timeidx=ALL_TIMES,method='cat')[:,:,:]
aut_ppt1 = getvar(aut, "RAINC", timeidx=ALL_TIMES,method='cat')[:,:,:]

win_ppt2 = getvar(win, "RAINNC", timeidx=ALL_TIMES,method='cat')[:,:,:]
spr_ppt2 = getvar(spr, "RAINNC", timeidx=ALL_TIMES,method='cat')[:,:,:]
mon_ppt2 = getvar(mon, "RAINNC", timeidx=ALL_TIMES,method='cat')[:,:,:]
aut_ppt2 = getvar(aut, "RAINNC", timeidx=ALL_TIMES,method='cat')[:,:,:]

win_ppt3 = win_ppt1+win_ppt2
spr_ppt3 = spr_ppt1+spr_ppt2
mon_ppt3 = mon_ppt1+mon_ppt2
aut_ppt3 = aut_ppt1+aut_ppt2


#=== Now take time average ============
win_ppt = np.average(win_ppt3,axis=0)
spr_ppt = np.average(spr_ppt3,axis=0)
mon_ppt = np.average(mon_ppt3,axis=0)
aut_ppt = np.average(aut_ppt3,axis=0)
all_ppt = (win_ppt+spr_ppt+mon_ppt+aut_ppt)/4

"""====  PBLH  ================="""
win_pb1 = getvar(win, "PBLH", timeidx=ALL_TIMES,method='cat')[:,:,:]
spr_pb1 = getvar(spr, "PBLH", timeidx=ALL_TIMES,method='cat')[:,:,:]
mon_pb1 = getvar(mon, "PBLH", timeidx=ALL_TIMES,method='cat')[:,:,:]
aut_pb1 = getvar(aut, "PBLH", timeidx=ALL_TIMES,method='cat')[:,:,:]


#=== Now take time average ============
win_pb = np.average(win_pb1,axis=0)
spr_pb = np.average(spr_pb1,axis=0)
mon_pb = np.average(mon_pb1,axis=0)
aut_pb = np.average(aut_pb1,axis=0)
all_pb = (win_pb+spr_pb+mon_pb+aut_pb)/4

"""====  RH  ================="""
win_rh1 = getvar(win, "rh2", timeidx=ALL_TIMES,method='cat')[:,:,:]
spr_rh1 = getvar(spr, "rh2", timeidx=ALL_TIMES,method='cat')[:,:,:]
mon_rh1 = getvar(mon, "rh2", timeidx=ALL_TIMES,method='cat')[:,:,:]
aut_rh1 = getvar(aut, "rh2", timeidx=ALL_TIMES,method='cat')[:,:,:]

#=== Now take time average ============
win_rh = np.average(win_rh1,axis=0)
spr_rh = np.average(spr_rh1,axis=0)
mon_rh = np.average(mon_rh1,axis=0)
aut_rh = np.average(aut_rh1,axis=0)
all_rh = (win_rh+spr_rh+mon_rh+aut_rh)/4

"""====  Wind data ========="""
#====== Get WRF data =====================================
p  = getvar(win, "pressure")
u1 = getvar(win, "ua", units="m s-1",timeidx=-1) #ALL_TIMES)
v1 = getvar(win, "va", units="m s-1",timeidx=-1) #ALL_TIMES)
u2 = getvar(spr, "ua", units="m s-1",timeidx=-1) #ALL_TIMES)
v2 = getvar(spr, "va", units="m s-1",timeidx=-1) #ALL_TIMES)
u3 = getvar(mon, "ua", units="m s-1",timeidx=-1) #ALL_TIMES)
v3 = getvar(mon, "va", units="m s-1",timeidx=-1) #ALL_TIMES)
u4 = getvar(aut, "ua", units="m s-1",timeidx=-1) #ALL_TIMES)
v4 = getvar(aut, "va", units="m s-1",timeidx=-1) #ALL_TIMES)

#=== Now interpolate at 500 hpa  ==============================
winu1 = interplevel(u1, p, 500)
winv1 = interplevel(v1, p, 500)
spru1 = interplevel(u2, p, 500)
sprv1 = interplevel(v2, p, 500)
sumu1 = interplevel(u3, p, 500)
sumv1 = interplevel(v3, p, 500)
autu1 = interplevel(u4, p, 500)
autv1 = interplevel(v4, p, 500)

all_u = (winu1+spru1+sumu1+autu1)/4
all_v = (winv1+sprv1+sumv1+autv1)/4

#===== Take time average ===========================
#winu = np.average(winu1,axis=0)
#winv = np.average(winv1,axis=0)
#spru = np.average(spru1,axis=0)
#sprv = np.average(sprv1,axis=0)
#sumu = np.average(sumu1,axis=0)
#sumv = np.average(sumv1,axis=0)
#autu = np.average(autu1,axis=0)
#autv = np.average(autv1,axis=0)
#===== calculate seasonal wind speed =============
wi_ws1 = np.sqrt(winu1**2+winv1**2)
sp_ws1 = np.sqrt(spru1**2+sprv1**2)
su_ws1 = np.sqrt(sumu1**2+sumv1**2)
au_ws1 = np.sqrt(autu1**2+autv1**2)

all_ws = (wi_ws1+sp_ws1+su_ws1+au_ws1)/4

#==== Normalize the data in uniform arrows ============
#win_u_wrf = winu/wi_ws1
#win_v_wrf = winv/wi_ws1
#spr_u_wrf = spru/sp_ws1
#spr_v_wrf = sprv/sp_ws1
#sum_u_wrf = sumu/su_ws1
#sum_v_wrf = sumv/su_ws1
#aut_u_wrf = autu/au_ws1
#aut_v_wrf = autv/au_ws1

#all_u = (win_u_wrf+spr_u_wrf+sum_u_wrf+aut_u_wrf)/4
#all_v = (win_v_wrf+spr_v_wrf+sum_v_wrf+aut_v_wrf)/4


#========== Get WRF lat lon  ===================
lat, lon = latlon_coords(p)


"""
#===============================================================#
#++++++++++++++++++++++ Now Plot +++++++++++++++++++++++++++++++#
===============================================================#"""

fig = plt.figure(figsize=(12,10),dpi=300)

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#=== Precipitation  ============================
""" winter"""
ax =plt.subplot(5,4,1)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,win_ppt,cmap=cmaps.MPL_YlGnBu)  
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Winter', labelpad=5,fontsize=12)
plt.clim(0,600)

""" Spring"""
ax =plt.subplot(5,4,5)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,spr_ppt,cmap=cmaps.MPL_YlGnBu) 
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Spring', labelpad=5,fontsize=12)
plt.clim(0,600)

""" monsoon """
ax =plt.subplot(5,4,9)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,mon_ppt,cmap=cmaps.MPL_YlGnBu)  
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Summer', labelpad=5,fontsize=12)
plt.clim(0,600)

""" autumn """
ax =plt.subplot(5,4,13)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,aut_ppt,cmap=cmaps.MPL_YlGnBu)  
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Autumn', labelpad=5,fontsize=12)
plt.clim(0,600)

""" all"""
ax1 =plt.subplot(5,4,17)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s1 = m.pcolormesh(x,y,all_ppt,cmap=cmaps.MPL_YlGnBu)  
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],fontsize=12,color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.ylabel('Annual', labelpad=5,fontsize=12)
plt.clim(0,600)

#=== Setting Colorbar  =========
cax1 = fig.add_axes([0.04,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s1, ax=ax1, cax=cax1, ticks=[0,100,200,300,400,500,600,700,800], orientation='horizontal')
cb1.set_label(r'Total Precipitation [mm]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)




#====== PBLH ==============================
""" winter"""
ax =plt.subplot(5,4,2)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,win_pb,cmap=cmaps.temp_19lev)  
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,800)

""" spring"""
ax =plt.subplot(5,4,6)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,spr_pb,cmap=cmaps.temp_19lev)  #CBR_coldhot) #BlueWhiteOrangeRed) #amwg_blueyellowred) #BlueWhiteOrangeRed) #MP>
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,800)

""" Monsoon"""
ax =plt.subplot(5,4,10)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,mon_pb,cmap=cmaps.temp_19lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,800)

""" Autumn """
ax =plt.subplot(5,4,14)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,aut_pb,cmap=cmaps.temp_19lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,800)

""" Annual """
ax2 =plt.subplot(5,4,18)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s2 = m.pcolormesh(x,y,all_pb,cmap=cmaps.temp_19lev)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, labels=[0,0,0,1],dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,800)


#=== Setting Colorbar for PBLH  =========
cax1 = fig.add_axes([0.285,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s2, ax=ax2, cax=cax1, ticks=[0,100,200,300,400,500,600,700,800], orientation='horizontal')
cb1.set_label(r'PBLH [m]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)



#===== RH  ===============================================================
""" winter"""
ax =plt.subplot(5,4,3)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,win_rh,cmap=cmaps.cmp_b2r)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(30,100)

""" spring """
ax =plt.subplot(5,4,7)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,spr_rh,cmap=cmaps.cmp_b2r)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(30,100)

""" monsoon"""
ax =plt.subplot(5,4,11)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,mon_rh,cmap=cmaps.cmp_b2r)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(30,100)

""" autumn """
ax =plt.subplot(5,4,15)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,aut_rh,cmap=cmaps.cmp_b2r)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(30,100)

""" Annual"""
ax3 =plt.subplot(5,4,19)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s3 = m.pcolormesh(x,y,all_rh,cmap=cmaps.cmp_b2r)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  labels=[0,0,0,1],color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(30,100)


#=== Setting Colorbar for RH =========
cax1 = fig.add_axes([0.528,0.1,0.19,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s3, ax=ax3, cax=cax1, ticks=[30,40,50,60,70,80,90,100], orientation='horizontal')
cb1.set_label(r'Relative Humidity [%]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)



#====== Wind =============================================================================
"""Winte"""
ax =plt.subplot(5,4,4)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,wi_ws1,cmap=cmaps.MPL_RdGy_r) #plot windspeed
q = plt.quiver(x[::5,::5], y[::5,::5], winu1[::5,::5],winv1[::5,::5],cmap='viridis',
                pivot='middle',scale_units='inches',headwidth=4,headlength=8)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,24)

#====== add quiver
qk = plt.quiverkey(q,
                   0.07, 0.06,                  # x,y label position
                   25, str(5)+''+u1.units, # choose units + update string
                   labelpos='E',labelcolor='magenta',                # add label to the right
                   coordinates='axes'
                   )



"""Spring"""
ax =plt.subplot(5,4,8)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,sp_ws1,cmap=cmaps.MPL_RdGy_r) #plot windspeed
q = plt.quiver(x[::5,::5], y[::5,::5], spru1[::5,::5],sprv1[::5,::5],cmap='viridis',
                pivot='middle',scale_units='inches',headwidth=4,headlength=8)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,24)

#====== add quiver
qk = plt.quiverkey(q,
                   0.07, 0.06,                  # x,y label position
                   25, str(5)+''+u1.units, # choose units + update string
                   labelpos='E',labelcolor='magenta',                # add label to the right
                   coordinates='axes'
                   )


"""Monsoon"""
ax =plt.subplot(5,4,12)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,su_ws1,cmap=cmaps.MPL_RdGy_r) #plot windspeed
q = plt.quiver(x[::5,::5], y[::5,::5], sumu1[::5,::5],sumv1[::5,::5],cmap='viridis',
                pivot='middle',scale_units='inches',headwidth=4,headlength=8)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,24)

#====== add quiver
qk = plt.quiverkey(q,
                   0.07, 0.06,                  # x,y label position
                   25, str(5)+''+u1.units, # choose units + update string
                   labelpos='E',labelcolor='magenta',                # add label to the right
                   coordinates='axes'
                   )



"""Autumn"""
ax =plt.subplot(5,4,16)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s = m.pcolormesh(x,y,au_ws1,cmap=cmaps.MPL_RdGy_r) #plot windspeed
q = plt.quiver(x[::5,::5], y[::5,::5], autu1[::5,::5],autv1[::5,::5],cmap='viridis',
                pivot='middle',scale_units='inches',headwidth=4,headlength=8)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,24)

#====== add quiver
qk = plt.quiverkey(q,
                   0.07, 0.06,                  # x,y label position
                   25, str(5)+''+u1.units, # choose units + update string
                   labelpos='E',labelcolor='magenta',                # add label to the right
                   coordinates='axes'
                   )


"""Annual"""
ax4 =plt.subplot(5,4,20)
m = Basemap(projection='lcc',width=12000000,height=9000000,area_thresh=1000,\
            lat_1=25,lat_2=25,lat_0=25,lon_0=90,llcrnrlon=45,llcrnrlat=0,\
             urcrnrlon=139,urcrnrlat=55,resolution='h')
x, y = m(to_np(lon), to_np(lat))
s4 = m.pcolormesh(x,y,all_ws,cmap=cmaps.MPL_RdGy_r) #plot windspeed
q = plt.quiver(x[::5,::5], y[::5,::5], all_u[::5,::5],all_v[::5,::5],cmap='viridis',
                pivot='middle',scale_units='inches',headwidth=4,headlength=8)
m.drawparallels(np.arange(10, 60, 15), linewidth=0.3, labels=[0,1,0,0],dashes=[2, 2],  fontsize=12,color='black')
m.drawmeridians(np.arange(45, 125,15), linewidth=0.3, labels=[0,0,0,1],dashes=[2, 2],  color='black')
m.drawcountries(linewidth=0.1)
m.drawcoastlines(linewidth=0.1)
plt.clim(0,24)

#====== add quiver
qk = plt.quiverkey(q,
                   0.07, 0.06,                  # x,y label position
                   25, str(5)+''+u1.units, # choose units + update string
                   labelpos='E',labelcolor='magenta',                # add label to the right
                   coordinates='axes'
                   )



#=== Setting Colorbar for dust =========
cax1 = fig.add_axes([0.775,0.1,0.185,0.01])   # left, bottom, width, and height
cb1 = fig.colorbar(s4, ax=ax4, cax=cax1, ticks=[0,3,6,9,12,15,18,21,24], orientation='horizontal')
cb1.set_label(r'Wind speed [$m s^{-1}$]',fontsize=12,labelpad=5)
cb1.ax.tick_params(labelsize=10)



""" Adjust"""
fig.subplots_adjust(top=0.96,
                        bottom=0.14,
                        left=0.015,
                        right=0.99,
                        hspace=0.05,
                        wspace=0.0)

plt.savefig("/mnt/g/2nd_Paper/met_new600dpi.png",dpi=600)

#plt.show()



