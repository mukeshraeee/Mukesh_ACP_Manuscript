"""This code generates the figure for PM,OC,Ext.Coeff
   Prepared by Rai. M 2021/16/1"""

#====== Code for vertical cros-seection plot===========================
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
import wrf
from wrf import (getvar, to_np, get_cartopy,get_basemap, latlon_coords, vertcross,ALL_TIMES,cartopy_xlim, cartopy_ylim, interpline, CoordPair)
import cmaps
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import scipy as sp
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


#===== get OC,PM2.5 and PM10 data from WRF data =====================
#=== pm2.5==
win_pm25 = getvar(win, "PM2_5_DRY", timeidx=-1,method='cat')[0:25,:,:]
spr_pm25 = getvar(spr, "PM2_5_DRY", timeidx=-1,method='cat')[0:25,:,:]
mon_pm25 = getvar(mon, "PM2_5_DRY", timeidx=-1,method='cat')[0:25,:,:]
aut_pm25 = getvar(aut, "PM2_5_DRY", timeidx=-1,method='cat')[0:25,:,:]

#=== pm10===
win_pm10 = getvar(win, "PM10", timeidx=-1,method='cat')[0:25,:,:]
spr_pm10 = getvar(spr, "PM10", timeidx=-1,method='cat')[0:25,:,:]
mon_pm10 = getvar(mon, "PM10", timeidx=-1,method='cat')[0:25,:,:]
aut_pm10 = getvar(aut, "PM10", timeidx=-1,method='cat')[0:25,:,:]

#===== extcoeff ==
win_ex = getvar(win, "EXTCOF55", timeidx=-1,method='cat')[0:25,:,:]
spr_ex = getvar(spr, "EXTCOF55", timeidx=-1,method='cat')[0:25,:,:]
mon_ex = getvar(mon, "EXTCOF55", timeidx=-1,method='cat')[0:25,:,:]
aut_ex = getvar(aut, "EXTCOF55", timeidx=-1,method='cat')[0:25,:,:]

#===== oc==
win_oc1 = getvar(win, "OC1", timeidx=-1,method='cat')[0:25,:,:]
spr_oc1 = getvar(spr, "OC1", timeidx=-1,method='cat')[0:25,:,:]
mon_oc1 = getvar(mon, "OC1", timeidx=-1,method='cat')[0:25,:,:]
aut_oc1 = getvar(aut, "OC1", timeidx=-1,method='cat')[0:25,:,:]

win_oc2 = getvar(win, "OC2", timeidx=-1,method='cat')[0:25,:,:]
spr_oc2 = getvar(spr, "OC2", timeidx=-1,method='cat')[0:25,:,:]
mon_oc2 = getvar(mon, "OC2", timeidx=-1,method='cat')[0:25,:,:]
aut_oc2 = getvar(aut, "OC2", timeidx=-1,method='cat')[0:25,:,:]

win_oc = win_oc1+win_oc2
spr_oc = spr_oc1+spr_oc2
mon_oc = mon_oc1+mon_oc2
aut_oc = aut_oc1+aut_oc2


"""======================== Get Wind DATA ======================================================"""
win_u = getvar(win, "ua", timeidx=-1,method='cat')[0:25,:,:] # all time,surface,all lat and lon
win_w = getvar(win, "wa", timeidx=-1,method='cat')[0:25,:,:]
spr_u = getvar(spr, "ua", timeidx=-1,method='cat')[0:25,:,:]
spr_w = getvar(spr, "wa", timeidx=-1,method='cat')[0:25,:,:]
mon_u = getvar(mon, "ua", timeidx=-1,method='cat')[0:25,:,:]
mon_w = getvar(mon, "wa", timeidx=-1,method='cat')[0:25,:,:]
aut_u = getvar(aut, "ua", timeidx=-1,method='cat')[0:25,:,:]
aut_w = getvar(aut, "wa", timeidx=-1,method='cat')[0:25,:,:]


"""======= Seasonal ter height ============================="""
win_ter      = getvar(win, "ter", units="km")
spr_ter      = getvar(spr, "ter", units="km")
mon_ter      = getvar(mon, "ter", units="km")
aut_ter      = getvar(aut, "ter", units="km")

"""============= Get model height for mass grid =========================================="""
ht_win  = getvar(win, "height_agl", units="km",method='cat')[0:25,:,:]
ht_spr  = getvar(spr, "height_agl", units="km",method='cat')[0:25,:,:]
ht_mon =  getvar(mon, "height_agl", units="km",method='cat')[0:25,:,:]
ht_aut  = getvar(aut, "height_agl", units="km",method='cat')[0:25,:,:]

"""============ Intrpolation WIND  data PM2.5 ================================================="""

#=== Seasonal  wind interpolation ======
w1 = 10**(win_w/10.) #winter
u1 = 10**(win_u/10.)
w2 = 10**(spr_w/10.) #spring
u2 = 10**(spr_u/10.)
w3 = 10**(mon_w/10.) #Monsoon
u3 = 10**(mon_u/10.)
w4 = 10**(aut_w/10.) #Autumn
u4 = 10**(aut_u/10.)


#========= PM25  =========================================
z1 = 10**(win_pm25/10.) # Use linear Z for interpolation
z2 = 10**(spr_pm25/10.)
z3 = 10**(mon_pm25/10.)
z4 = 10**(aut_pm25/10.)

#========= PM10  =========================================
z5 = 10**(win_pm10/10.) 
z6 = 10**(spr_pm10/10.)
z7 = 10**(mon_pm10/10.)
z8 = 10**(aut_pm10/10.)

#========= OC  =========================================
z9  = 10**(win_oc/10.)
z10 = 10**(spr_oc/10.)
z11 = 10**(mon_oc/10.)
z12 = 10**(aut_oc/10.)

#========= Ext. Coeff  =========================================
z13  = 10**(win_oc/10.)
z14 = 10**(spr_oc/10.)
z15 = 10**(mon_oc/10.)
z16 = 10**(aut_oc/10.)


"""========  Define the cross section start and end points========="""
cross_start = CoordPair(lat=26, lon=60)
cross_end   = CoordPair(lat=35, lon=119)

"""Compute the vertical cross-section interpolation.
   Also, include the lat/lon points along the cross-section in the metadata by setting latlon to True"""


#======= WIND ==============================================================================================
w1_cross = vertcross(w1, ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True) #winter
u1_cross = vertcross(u1, ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

w2_cross = vertcross(w2, ht_spr, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True) #spring
u2_cross = vertcross(u2, ht_spr, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

w3_cross = vertcross(w3, ht_mon, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True) #monsoon
u3_cross = vertcross(u3, ht_mon, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

w4_cross = vertcross(w4, ht_aut, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True) #autumn
u4_cross = vertcross(u4, ht_aut, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#======== PM25 ========================================================================================================
z1_cross = vertcross(z1, ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z2_cross = vertcross(z2, ht_spr, wrfin=spr,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z3_cross = vertcross(z3, ht_mon, wrfin=mon,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z4_cross = vertcross(z4, ht_aut, wrfin=aut,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#======== PM10 ========================================================================================================
z5_cross = vertcross(z5, ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z6_cross = vertcross(z6, ht_spr, wrfin=spr,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z7_cross = vertcross(z7, ht_mon, wrfin=mon,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z8_cross = vertcross(z8, ht_aut, wrfin=aut,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#======== OC ========================================================================================================
z9_cross  = vertcross(z9,  ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z10_cross = vertcross(z10, ht_spr, wrfin=spr,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z11_cross = vertcross(z11, ht_mon, wrfin=mon,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z12_cross = vertcross(z12, ht_aut, wrfin=aut,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)

#======== EXT.Coeff ========================================================================================================
z13_cross = vertcross(z13, ht_win, wrfin=win,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z14_cross = vertcross(z14, ht_spr, wrfin=spr,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z15_cross = vertcross(z15, ht_mon, wrfin=mon,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
z16_cross = vertcross(z16, ht_aut, wrfin=aut,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)


"""===== Convert back to dBz after interpolation =============================="""
"""Wind"""
w1_cross  = 10.0 * np.log10(w1_cross)
u1_cross  = 10.0 * np.log10(u1_cross)
w2_cross  = 10.0 * np.log10(w2_cross)
u2_cross  = 10.0 * np.log10(u2_cross)
w3_cross  = 10.0 * np.log10(w3_cross)
u3_cross  = 10.0 * np.log10(u3_cross)
w4_cross  = 10.0 * np.log10(w4_cross)
u4_cross  = 10.0 * np.log10(u4_cross)

""" PM25"""
z1_cross  = 10.0 * np.log10(z1_cross)
z2_cross  = 10.0 * np.log10(z2_cross)
z3_cross  = 10.0 * np.log10(z3_cross)
z4_cross  = 10.0 * np.log10(z4_cross)

""" PM10"""
z5_cross  = 10.0 * np.log10(z5_cross)
z6_cross  = 10.0 * np.log10(z6_cross)
z7_cross  = 10.0 * np.log10(z7_cross)
z8_cross  = 10.0 * np.log10(z8_cross)

""" OC"""
z9_cross   = 10.0 * np.log10(z9_cross)
z10_cross  = 10.0 * np.log10(z10_cross)
z11_cross  = 10.0 * np.log10(z11_cross)
z12_cross  = 10.0 * np.log10(z12_cross)

""" Ext. """
z13_cross  = 10.0 * np.log10(z13_cross)
z14_cross  = 10.0 * np.log10(z14_cross)
z15_cross  = 10.0 * np.log10(z15_cross)
z16_cross  = 10.0 * np.log10(z16_cross)

"""Make a copy of the z cross data. Let's use regular numpy arrays for this"""
#======= WIND =================================
w1_cross_filled = np.ma.copy(to_np(w1_cross))
u1_cross_filled = np.ma.copy(to_np(u1_cross))
w2_cross_filled = np.ma.copy(to_np(w2_cross))
u2_cross_filled = np.ma.copy(to_np(u2_cross))
w3_cross_filled = np.ma.copy(to_np(w3_cross))
u3_cross_filled = np.ma.copy(to_np(u3_cross))
w4_cross_filled = np.ma.copy(to_np(w4_cross))
u4_cross_filled = np.ma.copy(to_np(u4_cross))
#====== PM25 ============================================
z1_cross_filled = np.ma.copy(to_np(z1_cross))
z2_cross_filled = np.ma.copy(to_np(z2_cross))
z3_cross_filled = np.ma.copy(to_np(z3_cross))
z4_cross_filled = np.ma.copy(to_np(z4_cross))

#====== PM10 ============================================
z5_cross_filled = np.ma.copy(to_np(z5_cross))
z6_cross_filled = np.ma.copy(to_np(z6_cross))
z7_cross_filled = np.ma.copy(to_np(z7_cross))
z8_cross_filled = np.ma.copy(to_np(z8_cross))

#====== OC ============================================
z9_cross_filled  = np.ma.copy(to_np(z9_cross))
z10_cross_filled = np.ma.copy(to_np(z10_cross))
z11_cross_filled = np.ma.copy(to_np(z11_cross))
z12_cross_filled = np.ma.copy(to_np(z12_cross))

#====== ext.coeff ============================================
z13_cross_filled = np.ma.copy(to_np(z13_cross))
z14_cross_filled = np.ma.copy(to_np(z14_cross))
z15_cross_filled = np.ma.copy(to_np(z15_cross))
z16_cross_filled = np.ma.copy(to_np(z16_cross))

""" For each cross section column, find the first index with non-missing
    values and copy these to the missing elements below """

#======= pm25 ==================================================================
for i in range(z1_cross_filled.shape[-1]):
    column_vals = z1_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z1_cross_filled[0:first_idx, i] = z1_cross_filled[first_idx, i]

for i in range(z2_cross_filled.shape[-1]):
    column_vals = z2_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z2_cross_filled[0:first_idx, i] = z2_cross_filled[first_idx, i]

for i in range(z3_cross_filled.shape[-1]):
    column_vals = z3_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z3_cross_filled[0:first_idx, i] = z3_cross_filled[first_idx, i]

for i in range(z4_cross_filled.shape[-1]):
    column_vals = z4_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z4_cross_filled[0:first_idx, i] = z4_cross_filled[first_idx, i]

#======= pm10  ==================================================================
for i in range(z5_cross_filled.shape[-1]):
    column_vals = z5_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z5_cross_filled[0:first_idx, i] = z5_cross_filled[first_idx, i]

for i in range(z6_cross_filled.shape[-1]):
    column_vals = z6_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z6_cross_filled[0:first_idx, i] = z6_cross_filled[first_idx, i]

for i in range(z7_cross_filled.shape[-1]):
    column_vals = z7_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z7_cross_filled[0:first_idx, i] = z7_cross_filled[first_idx, i]

for i in range(z8_cross_filled.shape[-1]):
    column_vals = z8_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z8_cross_filled[0:first_idx, i] = z8_cross_filled[first_idx, i]

#===== OC  =====================================================================
for i in range(z9_cross_filled.shape[-1]):
    column_vals = z9_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z9_cross_filled[0:first_idx, i] = z9_cross_filled[first_idx, i]

for i in range(z10_cross_filled.shape[-1]):
    column_vals = z10_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z10_cross_filled[0:first_idx, i] = z10_cross_filled[first_idx, i]

for i in range(z11_cross_filled.shape[-1]):
    column_vals = z11_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z11_cross_filled[0:first_idx, i] = z11_cross_filled[first_idx, i]

for i in range(z12_cross_filled.shape[-1]):
    column_vals = z12_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z12_cross_filled[0:first_idx, i] = z12_cross_filled[first_idx, i]

#===== EXT. Coeff =====================================================================
for i in range(z13_cross_filled.shape[-1]):
    column_vals = z13_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z13_cross_filled[0:first_idx, i] = z13_cross_filled[first_idx, i]

for i in range(z14_cross_filled.shape[-1]):
    column_vals = z14_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z14_cross_filled[0:first_idx, i] = z14_cross_filled[first_idx, i]

for i in range(z15_cross_filled.shape[-1]):
    column_vals = z15_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z15_cross_filled[0:first_idx, i] = z15_cross_filled[first_idx, i]

for i in range(z16_cross_filled.shape[-1]):
    column_vals = z16_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    z16_cross_filled[0:first_idx, i] = z16_cross_filled[first_idx, i]

#===========================================
xs1 = np.arange(0, z1_cross.shape[-1], 1)
ys1 = to_np(z1_cross.coords["vertical"])

xs2 = np.arange(0, z2_cross.shape[-1], 1)
ys2 = to_np(z2_cross.coords["vertical"])

xs3 = np.arange(0, z3_cross.shape[-1], 1)
ys3 = to_np(z3_cross.coords["vertical"])

xs4 = np.arange(0, z4_cross.shape[-1], 1)
ys4 = to_np(z4_cross.coords["vertical"])

xs5 = np.arange(0, z5_cross.shape[-1], 1)
ys5 = to_np(z5_cross.coords["vertical"])

xs6 = np.arange(0, z6_cross.shape[-1], 1)
ys6 = to_np(z6_cross.coords["vertical"])

xs7 = np.arange(0, z7_cross.shape[-1], 1)
ys7 = to_np(z7_cross.coords["vertical"])

xs8 = np.arange(0, z8_cross.shape[-1], 1)
ys8 = to_np(z8_cross.coords["vertical"])

xs9 = np.arange(0, z9_cross.shape[-1], 1)
ys9 = to_np(z9_cross.coords["vertical"])

xs10 = np.arange(0, z10_cross.shape[-1], 1)
ys10 = to_np(z10_cross.coords["vertical"])

xs11 = np.arange(0, z11_cross.shape[-1], 1)
ys11 = to_np(z11_cross.coords["vertical"])

xs12 = np.arange(0, z12_cross.shape[-1], 1)
ys12 = to_np(z12_cross.coords["vertical"])

xs13 = np.arange(0, z13_cross.shape[-1], 1)
ys13 = to_np(z13_cross.coords["vertical"])

xs14 = np.arange(0, z14_cross.shape[-1], 1)
ys14 = to_np(z14_cross.coords["vertical"])

xs15 = np.arange(0, z15_cross.shape[-1], 1)
ys15 = to_np(z15_cross.coords["vertical"])

xs16 = np.arange(0, z16_cross.shape[-1], 1)
ys16 = to_np(z16_cross.coords["vertical"])


"""======Plot data =================================================="""

fig, axs = plt.subplots(4, 4,figsize=(6.2,5),dpi=300)

"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#===============  Get the lat/lon points =============
lats, lons = latlon_coords(win_pm25)
cart_proj  = get_cartopy(win_pm25)

pm25_levels = np.arange(0,105,5)
pm10_levels = np.arange(0,105,5)
oc_levels   = np.arange(0,30.5,0.5)
ext_levels  = np.arange(0,1.01,0.01)


#=== PM25 ===============
p25_win = axs[0,0].contourf(xs1,ys1,to_np(z1_cross_filled),levels=pm25_levels,cmap=cmaps.WhiteYellowOrangeRed)
p25_spr = axs[0,1].contourf(xs2,ys2,to_np(z2_cross_filled),levels=pm25_levels,cmap=cmaps.WhiteYellowOrangeRed)
p25_mon = axs[0,2].contourf(xs3,ys3,to_np(z3_cross_filled),levels=pm25_levels,cmap=cmaps.WhiteYellowOrangeRed)
p25_aut = axs[0,3].contourf(xs4,ys4,to_np(z4_cross_filled),levels=pm25_levels,cmap=cmaps.WhiteYellowOrangeRed)

#=== PM10 ===============
p10_win = axs[1,0].contourf(xs5,ys5,to_np(z5_cross_filled),levels=pm10_levels,cmap=cmaps.WhiteYellowOrangeRed)
p10_spr = axs[1,1].contourf(xs6,ys6,to_np(z6_cross_filled),levels=pm10_levels,cmap=cmaps.WhiteYellowOrangeRed)
p10_mon = axs[1,2].contourf(xs7,ys7,to_np(z7_cross_filled),levels=pm10_levels,cmap=cmaps.WhiteYellowOrangeRed)
p10_aut = axs[1,3].contourf(xs8,ys8,to_np(z8_cross_filled),levels=pm10_levels,cmap=cmaps.WhiteYellowOrangeRed)
#=== OC ===============
oc_win = axs[2,0].contourf(xs9,ys9,  to_np(z9_cross_filled), levels=oc_levels,cmap=cmaps.WhiteYellowOrangeRed)
oc_spr = axs[2,1].contourf(xs10,ys10,to_np(z10_cross_filled),levels=oc_levels,cmap=cmaps.WhiteYellowOrangeRed)
oc_mon = axs[2,2].contourf(xs11,ys11,to_np(z11_cross_filled),levels=oc_levels,cmap=cmaps.WhiteYellowOrangeRed)
oc_aut = axs[2,3].contourf(xs12,ys12,to_np(z12_cross_filled),levels=oc_levels,cmap=cmaps.WhiteYellowOrangeRed)
#=== Ext.coeff ===============
ext_win = axs[3,0].contourf(xs13,ys13,to_np(z13_cross_filled),levels=ext_levels,cmap=cmaps.WhiteYellowOrangeRed)
ext_spr = axs[3,1].contourf(xs14,ys14,to_np(z14_cross_filled),levels=ext_levels,cmap=cmaps.WhiteYellowOrangeRed)
ext_mon = axs[3,2].contourf(xs15,ys15,to_np(z15_cross_filled),levels=ext_levels,cmap=cmaps.WhiteYellowOrangeRed)
ext_aut = axs[3,3].contourf(xs16,ys16,to_np(z16_cross_filled),levels=ext_levels,cmap=cmaps.WhiteYellowOrangeRed)

""" ====== Plot Wind data ============================================"""
#======== pm25  =============================================================
axs[0,0].quiver(xs1[::4], ys1[::4],to_np(u1_cross_filled[::4, ::4]), to_np(w1_cross_filled[::4, ::4]*500))
axs[0,1].quiver(xs2[::4], ys2[::4],to_np(u2_cross_filled[::4, ::4]), to_np(w2_cross_filled[::4, ::4]*500))
axs[0,2].quiver(xs3[::4], ys3[::4],to_np(u3_cross_filled[::4, ::4]), to_np(w3_cross_filled[::4, ::4]*500))
axs[0,3].quiver(xs4[::4], ys4[::4],to_np(u4_cross_filled[::4, ::4]), to_np(w4_cross_filled[::4, ::4]*500))
#======== pm10  =============================================================
axs[1,0].quiver(xs5[::4], ys5[::4],to_np(u1_cross_filled[::4, ::4]), to_np(w1_cross_filled[::4, ::4]*500))
axs[1,1].quiver(xs6[::4], ys6[::4],to_np(u2_cross_filled[::4, ::4]), to_np(w2_cross_filled[::4, ::4]*500))
axs[1,2].quiver(xs6[::4], ys7[::4],to_np(u3_cross_filled[::4, ::4]), to_np(w3_cross_filled[::4, ::4]*500))
axs[1,3].quiver(xs8[::4], ys8[::4],to_np(u4_cross_filled[::4, ::4]), to_np(w4_cross_filled[::4, ::4]*500))
#======== oc  =============================================================
axs[2,0].quiver(xs9[::4],  ys9[::4],to_np(u1_cross_filled[::4, ::4]), to_np(w1_cross_filled[::4, ::4]*500))
axs[2,1].quiver(xs10[::4], ys10[::4],to_np(u2_cross_filled[::4, ::4]), to_np(w2_cross_filled[::4, ::4]*500))
axs[2,2].quiver(xs11[::4], ys11[::4],to_np(u3_cross_filled[::4, ::4]), to_np(w3_cross_filled[::4, ::4]*500))
axs[2,3].quiver(xs12[::4], ys12[::4],to_np(u4_cross_filled[::4, ::4]), to_np(w4_cross_filled[::4, ::4]*500))
#======== ext.coeff  =============================================================
axs[3,0].quiver(xs13[::4], ys13[::4],to_np(u1_cross_filled[::4, ::4]), to_np(w1_cross_filled[::4, ::4]*500))
axs[3,1].quiver(xs14[::4], ys14[::4],to_np(u2_cross_filled[::4, ::4]), to_np(w2_cross_filled[::4, ::4]*500))
axs[3,2].quiver(xs15[::4], ys15[::4],to_np(u3_cross_filled[::4, ::4]), to_np(w3_cross_filled[::4, ::4]*500))
axs[3,3].quiver(xs16[::4], ys16[::4],to_np(u4_cross_filled[::4, ::4]), to_np(w4_cross_filled[::4, ::4]*500))

"""====  Get the terrain heights along the cross section line"""
ter_line_win = interpline(win_ter, wrfin=win, start_point=cross_start,end_point=cross_end)
ter_line_spr = interpline(spr_ter, wrfin=spr, start_point=cross_start,end_point=cross_end)
ter_line_mon = interpline(mon_ter, wrfin=mon, start_point=cross_start,end_point=cross_end)
ter_line_aut = interpline(aut_ter, wrfin=aut, start_point=cross_start,end_point=cross_end)
"""#=========  Fill in the mountain area ========================================="""
#==== pm25 ===============================
axs[0,0].fill_between(xs1, 0, to_np(ter_line_win),facecolor="dimgrey")
axs[0,1].fill_between(xs2, 0, to_np(ter_line_spr),facecolor="dimgrey")
axs[0,2].fill_between(xs3, 0, to_np(ter_line_mon),facecolor="dimgrey")
axs[0,3].fill_between(xs4, 0, to_np(ter_line_aut),facecolor="dimgrey")
#==== pm10 ===============================
axs[1,0].fill_between(xs5, 0, to_np(ter_line_win),facecolor="dimgrey")
axs[1,1].fill_between(xs6, 0, to_np(ter_line_spr),facecolor="dimgrey")
axs[1,2].fill_between(xs7, 0, to_np(ter_line_mon),facecolor="dimgrey")
axs[1,3].fill_between(xs8, 0, to_np(ter_line_aut),facecolor="dimgrey")
#==== oc ===============================
axs[2,0].fill_between(xs9,  0, to_np(ter_line_win),facecolor="dimgrey")
axs[2,1].fill_between(xs10, 0, to_np(ter_line_spr),facecolor="dimgrey")
axs[2,2].fill_between(xs11, 0, to_np(ter_line_mon),facecolor="dimgrey")
axs[2,3].fill_between(xs12, 0, to_np(ter_line_aut),facecolor="dimgrey")
#==== ext.coeff ===============================
axs[3,0].fill_between(xs13, 0, to_np(ter_line_win),facecolor="dimgrey")
axs[3,1].fill_between(xs14, 0, to_np(ter_line_spr),facecolor="dimgrey")
axs[3,2].fill_between(xs15, 0, to_np(ter_line_mon),facecolor="dimgrey")
axs[3,3].fill_between(xs16, 0, to_np(ter_line_aut),facecolor="dimgrey")

""" Set the x-ticks to use latitude and longitude labels"""
coord_pairs = to_np(z1_cross.coords["xy_loc"])
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = [pair.latlon_str(fmt="{:.1f}, {:.1f}")
              for pair in to_np(coord_pairs)]
#=======  Set the desired number of x ticks below
num_ticks = 6
thin = int((len(x_ticks) / num_ticks) + .8)


#===== Se tick====================
axs[0,0].set_xticklabels([])
#axs[0,1].set_yticklabels([])
axs[0,1].set_xticklabels([])
axs[0,2].set_yticklabels([])
axs[0,2].set_xticklabels([])
axs[0,3].set_yticklabels([])
axs[0,3].set_xticklabels([])

axs[1,0].set_xticklabels([])
#axs[1,0].set_yticklabels([])
axs[1,1].set_xticklabels([])
axs[1,1].set_yticklabels([])
axs[1,2].set_yticklabels([])
axs[1,2].set_xticklabels([])

axs[1,3].set_yticklabels([])
axs[1,3].set_xticklabels([])
axs[2,0].set_xticklabels([])
#axs[2,0].set_yticklabels([])
axs[2,1].set_xticklabels([])
axs[2,1].set_yticklabels([])
axs[2,2].set_yticklabels([])
axs[2,2].set_xticklabels([])

axs[2,3].set_yticklabels([])
axs[2,3].set_xticklabels([])
axs[3,0].set_xticklabels([])
#axs[3,0].set_yticklabels([])
axs[3,1].set_xticklabels([])
axs[3,1].set_yticklabels([])
axs[3,2].set_yticklabels([])
axs[3,3].set_yticklabels([])

axs[0,0].set_xticks(x_ticks[::thin])
axs[0,1].set_xticks(x_ticks[::thin])
axs[0,2].set_xticks(x_ticks[::thin])
axs[0,3].set_xticks(x_ticks[::thin])
axs[1,0].set_xticks(x_ticks[::thin])
axs[1,1].set_xticks(x_ticks[::thin])
axs[1,2].set_xticks(x_ticks[::thin])
axs[1,3].set_xticks(x_ticks[::thin])
axs[2,0].set_xticks(x_ticks[::thin])
axs[2,1].set_xticks(x_ticks[::thin])
axs[2,2].set_xticks(x_ticks[::thin])
axs[2,3].set_xticks(x_ticks[::thin])
axs[3,0].set_xticks(x_ticks[::thin])
axs[3,1].set_xticks(x_ticks[::thin])
axs[3,2].set_xticks(x_ticks[::thin])
axs[3,3].set_xticks(x_ticks[::thin])

#===== Set xtick label =============================================
axs[3,0].set_xticklabels(x_labels[::thin],rotation=25, fontsize=4)
axs[3,1].set_xticklabels(x_labels[::thin],rotation=25, fontsize=4)
axs[3,2].set_xticklabels(x_labels[::thin],rotation=25, fontsize=4)
axs[3,3].set_xticklabels(x_labels[::thin],rotation=25, fontsize=4)

#====  Set ytick label =============================================
axs[0,0].yaxis.set_tick_params(labelsize=4)
axs[1,0].yaxis.set_tick_params(labelsize=4)
axs[2,0].yaxis.set_tick_params(labelsize=4)
axs[3,0].yaxis.set_tick_params(labelsize=4)


axs[0,1].set_yticklabels([])
axs[0,2].set_yticklabels([])
axs[0,3].set_yticklabels([])

#==== set colorbar ======
b1 =plt.colorbar(p25_aut,ax=axs[0,3],ticks=[0,10,20,30,40,50,60,70,80,90,100,110])
b1.set_label(r'PM2.5 [$μg m ^{-3}$]', rotation=-270, fontsize=5,labelpad=4)
b1.ax.tick_params(labelsize=4)

b2 = plt.colorbar(p10_aut,ax=axs[1,3],ticks=[0,10,20,30,40,50,60,70,80,90,100,110])
b2.set_label(r'PM10 [$μg m ^{-3}$]', rotation=-270, fontsize=5,labelpad=4)
b2.ax.tick_params(labelsize=4)

b3 = plt.colorbar(oc_aut, ax=axs[2,3],ticks=[0,5,10,15,20,25,30,30])
b3.set_label(r'OC [$μg m ^{-3}$]', rotation=-270, fontsize=5,labelpad=4)
b3.ax.tick_params(labelsize=4)

b4 = plt.colorbar(ext_aut,ax=axs[3,3],ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1])
b4.set_label(r'Ext.Coeff [$Km ^{-1}$]', rotation=-270, fontsize=5,labelpad=4)
b4.ax.tick_params(labelsize=4)

#====== Add x-axis title ==========
fig.text(0.018, 0.53, 'Height [Km]', va='center', rotation='vertical',fontsize=5)
#==== PLot title ===========
axs[0,0].set_title('Winter',fontsize=5,pad=5)
axs[0,1].set_title('Spring',fontsize=5,pad=5)
axs[0,2].set_title('Summer',fontsize=5,pad=5)
axs[0,3].set_title('Autumn',fontsize=5,pad=5)
#==== Plot y-axis ===========
#axs[0,0].set_ylabel('PM2.5',fontsize=5)
#axs[1,0].set_ylabel('PM10',fontsize=5)
#axs[2,0].set_ylabel('OC',fontsize=5)
#axs[3,0].set_ylabel('Ext. Coeff',fontsize=5)



"""#======= Add cross-sectional map ====="""
sub_axes = plt.axes([0.83, .86, 0.12, 0.105])  #left, bottom, width, height
m = Basemap(projection='ortho', lat_0=30, lon_0=86,resolution='h',area_thresh=1000)
m.drawmapboundary(fill_color='black', linewidth=0) #A6CAE0
m.fillcontinents(color='white', alpha=0.7, lake_color='grey')
m.drawcoastlines(linewidth=0.1, color='black')
m.drawcountries(linewidth=0.1,  color='black')

m.drawparallels(np.arange(-90, 90, 30),linewidth=0.2,  dashes=[2, 2],  color='white')
m.drawmeridians(np.arange(-180, 180,30),linewidth=0.2, dashes=[2, 2],  color='white')


lon = np.linspace(60,119)
lat = np.linspace(26,35)
x,y = m(lon,lat)

m.plot(x,y, linewidth=0.6, color='red',linestyle='-')

# Cross-section to be plotted ==========================
#startlat = 29; startlon = 60
#arrlat   = 29; arrlon   = 119
#m.drawgreatcircle(startlon,startlat,arrlon,arrlat, linewidth=0.6, color='orange')





#===== Adjust figure size =========================
fig.subplots_adjust(top=0.965,
                        bottom=0.095,
                        left=0.067,
                        right=0.97,
                        hspace=0.105,
                        wspace=0.12)

#plt.savefig("/mnt/g/2nd_Paper/vertical.png",dpi=1000)
plt.show()
