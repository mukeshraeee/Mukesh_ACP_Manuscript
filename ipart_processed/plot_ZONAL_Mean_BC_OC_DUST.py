""" This code produce zonal average of integrated vertical transport of aerosol
    the output from IPART"""
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
from mpl_toolkits.basemap import Basemap
import cmaps
import scipy as sp
import scipy.ndimage
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from geocat.viz import util as gvutil
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

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
lon = np.array(lon)
lat = np.array(lat)
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

"""
  Prepare BC data """
bc_win2 = np.average(bc_win1,axis=0) #time average
bc_spr2 = np.average(bc_spr1,axis=0)
bc_mon2 = np.average(bc_mon1,axis=0)
bc_aut2 = np.average(bc_aut1,axis=0)

bc_win3 = np.average(bc_win2,axis=0)*100000 #lat average
bc_spr3 = np.average(bc_spr2,axis=0)*100000
bc_mon3 = np.average(bc_mon2,axis=0)*100000
bc_aut3 = np.average(bc_aut2,axis=0)*100000

bc_all = (bc_win3+bc_spr3+bc_mon3+bc_aut3)/4

"""=======Compute curves of interest ======================="""
n_steps=2   #number of rolling steps for the mean/std.
#===== winter => BC[mean,SD]
bc_win         = pd.DataFrame(bc_win3)
bc_win_mean    = bc_win.rolling(n_steps).mean()
bc_win_sd      = 2 * bc_win.rolling(n_steps).std()
#===== spring => BC[mean,SD]
bc_spr         = pd.DataFrame(bc_spr3)
bc_spr_mean    = bc_spr.rolling(n_steps).mean()
bc_spr_sd      = 2 * bc_spr.rolling(n_steps).std()
#===== summer => BC[mean,SD]
bc_mon         = pd.DataFrame(bc_mon3)
bc_mon_mean    = bc_mon.rolling(n_steps).mean()
bc_mon_sd      = 2 * bc_mon.rolling(n_steps).std()
#===== autumn => BC[mean,SD]
bc_aut         = pd.DataFrame(bc_aut3)
bc_aut_mean    = bc_aut.rolling(n_steps).mean()
bc_aut_sd      = 2 * bc_aut.rolling(n_steps).std()
#== all=>BC ==============
bc_all         = pd.DataFrame(bc_all)
bc_all_mean    = bc_all.rolling(n_steps).mean()
bc_all_sd      = 2 * bc_all.rolling(n_steps).std()
#==== Upper & Lower values =================
bc_win_lower  = (bc_win_mean-bc_win_sd)[0]
bc_win_upper  = (bc_win_mean+bc_win_sd)[0]
bc_spr_lower  = (bc_spr_mean-bc_spr_sd)[0]
bc_spr_upper  = (bc_spr_mean+bc_spr_sd)[0]
bc_mon_lower  = (bc_mon_mean-bc_mon_sd)[0]
bc_mon_upper  = (bc_mon_mean+bc_mon_sd)[0]
bc_aut_lower  = (bc_aut_mean-bc_aut_sd)[0]
bc_aut_upper  = (bc_aut_mean+bc_aut_sd)[0]
bc_all_lower  = (bc_all_mean-bc_all_sd)[0]
bc_all_upper  = (bc_all_mean+bc_all_sd)[0]
"""
  Prepare OC data """
oc_win2 = np.average(oc_win1,axis=0) #time average
oc_spr2 = np.average(oc_spr1,axis=0)
oc_mon2 = np.average(oc_mon1,axis=0)
oc_aut2 = np.average(oc_aut1,axis=0)
oc_win3 = np.average(oc_win2,axis=0)*100000 #lat average
oc_spr3 = np.average(oc_spr2,axis=0)*100000
oc_mon3 = np.average(oc_mon2,axis=0)*100000
oc_aut3 = np.average(oc_aut2,axis=0)*100000
oc_all = (oc_win3+oc_spr3+oc_mon3+oc_aut3)/4
"""=======Compute curves of interest ======================="""
#===== winter => OC[mean,SD]
oc_win         = pd.DataFrame(oc_win3)
oc_win_mean    = oc_win.rolling(n_steps).mean()
oc_win_sd      = 1 * oc_win.rolling(n_steps).std()
#===== spring => OC[mean,SD]
oc_spr         = pd.DataFrame(oc_spr3)
oc_spr_mean    = oc_spr.rolling(n_steps).mean()
oc_spr_sd      = 2 * oc_spr.rolling(n_steps).std()
#===== summer => OC[mean,SD]
oc_mon         = pd.DataFrame(oc_mon3)
oc_mon_mean    = oc_mon.rolling(n_steps).mean()
oc_mon_sd      = 2 * oc_mon.rolling(n_steps).std()
#===== autumn => OC[mean,SD]
oc_aut         = pd.DataFrame(oc_aut3)
oc_aut_mean    = oc_aut.rolling(n_steps).mean()
oc_aut_sd      = 2 * oc_aut.rolling(n_steps).std()
#== all=>OC ==============
oc_all         = pd.DataFrame(bc_all)
oc_all_mean    = bc_all.rolling(n_steps).mean()
oc_all_sd      = 2 * bc_all.rolling(n_steps).std()

oc_win_lower  = (oc_win_mean-oc_win_sd)[0]
oc_win_upper  = (oc_win_mean+oc_win_sd)[0]
oc_spr_lower  = (oc_spr_mean-oc_spr_sd)[0]
oc_spr_upper  = (oc_spr_mean+oc_spr_sd)[0]
oc_mon_lower  = (oc_mon_mean-oc_mon_sd)[0]
oc_mon_upper  = (oc_mon_mean+oc_mon_sd)[0]
oc_aut_lower  = (oc_aut_mean-oc_aut_sd)[0]
oc_aut_upper  = (oc_aut_mean+oc_aut_sd)[0]
oc_all_lower  = (oc_all_mean-oc_all_sd)[0]
oc_all_upper  = (oc_all_mean+oc_all_sd)[0]

"""
  Prepare DUST data """
dust_win2 = np.average(dust_win1,axis=0) #time average
dust_spr2 = np.average(dust_spr1,axis=0)
dust_mon2 = np.average(dust_mon1,axis=0)
dust_aut2 = np.average(dust_aut1,axis=0)
dust_win3 = np.average(dust_win2,axis=0)*100000 #lat average
dust_spr3 = np.average(dust_spr2,axis=0)*100000
dust_mon3 = np.average(dust_mon2,axis=0)*100000
dust_aut3 = np.average(dust_aut2,axis=0)*100000
dust_all = (dust_win3+dust_spr3+dust_mon3+dust_aut3)/4
"""=======Compute curves of interest ======================="""
#===== winter => OC[mean,SD]
dust_win         = pd.DataFrame(dust_win3)
dust_win_mean    = dust_win.rolling(n_steps).mean()
dust_win_sd      = 2 * dust_win.rolling(n_steps).std()
#===== spring => OC[mean,SD]
dust_spr         = pd.DataFrame(dust_spr3)
dust_spr_mean    = dust_spr.rolling(n_steps).mean()
dust_spr_sd      = 2 * dust_spr.rolling(n_steps).std()
#===== summer => OC[mean,SD]
dust_mon         = pd.DataFrame(dust_mon3)
dust_mon_mean    = dust_mon.rolling(n_steps).mean()
dust_mon_sd      = 2 * dust_mon.rolling(n_steps).std()
#===== autumn => OC[mean,SD]
dust_aut         = pd.DataFrame(dust_aut3)
dust_aut_mean    = dust_aut.rolling(n_steps).mean()
dust_aut_sd      = 2 * dust_aut.rolling(n_steps).std()
#== all=>OC ==============
dust_all         = pd.DataFrame(dust_all)
dust_all_mean    = dust_all.rolling(n_steps).mean()
dust_all_sd      = 2 * dust_all.rolling(n_steps).std()

dust_win_lower  = (dust_win_mean-dust_win_sd)[0]
dust_win_upper  = (dust_win_mean+dust_win_sd)[0]
dust_spr_lower  = (dust_spr_mean-dust_spr_sd)[0]
dust_spr_upper  = (dust_spr_mean+dust_spr_sd)[0]
dust_mon_lower  = (dust_mon_mean-dust_mon_sd)[0]
dust_mon_upper  = (dust_mon_mean+dust_mon_sd)[0]
dust_aut_lower  = (dust_aut_mean-dust_aut_sd)[0]
dust_aut_upper  = (dust_aut_mean+dust_aut_sd)[0]
dust_all_lower  = (dust_all_mean-dust_all_sd)[0]
dust_all_upper  = (dust_all_mean+dust_all_sd)[0]




""" PLOT"""

"""==== seasonal BC=================="""
fig,(ax0,ax1,ax3)= plt.subplots(3,1,figsize=(8,3),dpi=300)

ax0 = plt.subplot(1,3,1)

"""=== plot BC mean value =="""

ax0.plot(bc_win,color = "dodgerblue",  linewidth=0.5) #Latitude
ax0.plot(bc_spr,color = "red",  linewidth=0.5)
ax0.plot(bc_mon,color = "crimson",  linewidth=0.5)
ax0.plot(bc_aut,color = "navy",  linewidth=0.5)
ax0.plot(bc_all,color = "black",  linewidth=0.5)

#=== plot BC 1d std====================================
ax0.fill_between(bc_win_sd.index, bc_win_lower, bc_win_upper, color='lightblue', alpha=0.8)
ax0.fill_between(bc_spr_sd.index, bc_spr_lower, bc_spr_upper, color='lightpink', alpha=0.5)
ax0.fill_between(bc_mon_sd.index, bc_mon_lower, bc_mon_upper, color='fuchsia', alpha=0.1)
ax0.fill_between(bc_aut_sd.index, bc_aut_lower, bc_aut_upper, color='blue', alpha=0.2)
ax0.fill_between(bc_all_sd.index, bc_all_lower, bc_all_upper, color='darkgray', alpha=0.8)


#==== highlight desired region =====
ax0.axvspan(50, 80, color='silver', alpha=0.3)  #Sout Asia  region
ax0.axvspan(90, 120, color='silver', alpha=0.3) #East China region

#=== setting axis===========================================================
ax0.tick_params(axis='y', labelsize=5)   
ax0.tick_params(axis='x', labelsize=5)
gvutil.add_major_minor_ticks(ax0,x_minor_per_major=4,y_minor_per_major=3,labelsize=5)
#==== Label handle ====
gvutil.set_axes_limits_and_ticks(ax0, xlim=(45,120),xticklabels=['40E','50E','60E','70E','80E','90E','100E','110E','120'])
gvutil.set_axes_limits_and_ticks(ax0, ylim=(0,0.8))
gvutil.set_titles_and_labels(ax0,ylabel='BC_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]',labelfontsize=5)

#=== add text ========================================
ax0.text(65,0.5,"Contribution from \n South Asia + IGP",fontsize=5,rotation='vertical')
ax0.text(100,0.5,"Contribution from \n South East Asia ",fontsize=5,rotation='vertical')

"""==== seasonal OC=================="""
ax1 = plt.subplot(1,3,2)

"""=== plot BC mean value =="""

ax1.plot(oc_win,color = "dodgerblue",  linewidth=0.5) #Latitude
ax1.plot(oc_spr,color = "red",  linewidth=0.5)
ax1.plot(oc_mon,color = "crimson",  linewidth=0.5)
ax1.plot(oc_aut,color = "navy",  linewidth=0.5)
ax1.plot(oc_all,color = "black",  linewidth=0.5)
#=== plot BC 1d std====================================
ax1.fill_between(oc_win_sd.index, oc_win_lower, oc_win_upper, color='lightblue', alpha=0.8)
ax1.fill_between(oc_spr_sd.index, oc_spr_lower, oc_spr_upper, color='lightpink', alpha=0.5)
ax1.fill_between(oc_mon_sd.index, oc_mon_lower, oc_mon_upper, color='fuchsia', alpha=0.1)
ax1.fill_between(oc_aut_sd.index, oc_aut_lower, oc_aut_upper, color='blue', alpha=0.2)
ax1.fill_between(oc_all_sd.index, oc_all_lower, oc_all_upper, color='darkgray', alpha=0.8)


#==== highlight desired region => separate out spring anamoly =============
ax1.axhspan(2.3, 5.7, color='silver',  alpha=0.3)  # Spring
ax1.axhline(y=5.55,  xmin=0.633,  xmax=0.665, linewidth=0.7)
ax1.axhline(y=0.25, xmin=0.633,  xmax=0.665, linewidth=0.7)
ax1.axvline(x=94,  ymin=0.045, ymax=0.92, linewidth=0.5,linestyle='--',color='teal')


#==== highlight vertical region => East Asia ============
ax1.axvspan(85, 105, color='silver', alpha=0.5)
#==== add text ========================================
ax1.text(95,2.4,"Spring anamoly over East Asia",fontsize=5,rotation='vertical')


#=== setting axis===========================================================
ax1.tick_params(axis='y', labelsize=5)
ax1.tick_params(axis='x', labelsize=5)
gvutil.add_major_minor_ticks(ax1,x_minor_per_major=4,y_minor_per_major=3,labelsize=5)
#==== Label handle ====
gvutil.set_axes_limits_and_ticks(ax1, xlim=(45,120),xticklabels=['40E','50E','60E','70E','80E','90E','100E','110E','120'])
gvutil.set_axes_limits_and_ticks(ax1, ylim=(0,6))
gvutil.set_titles_and_labels(ax1,ylabel='OC_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]',labelfontsize=5)


"""==== seasonal OC=================="""

ax3 = plt.subplot(1,3,3)

"""=== plot BC mean value =="""
ax3.plot(dust_win,color = "dodgerblue",  linewidth=0.5) #Latitude
ax3.plot(dust_spr,color = "red",  linewidth=0.5)
ax3.plot(dust_mon,color = "crimson",  linewidth=0.5)
ax3.plot(dust_aut,color = "navy",  linewidth=0.5)
ax3.plot(dust_all,color = "black",  linewidth=0.5)

""" Insert Legend"""
win = mlines.Line2D([], [], color='dodgerblue', label='Winter',linewidth=0.8)
spr = mlines.Line2D([], [], color='red',        label='Spring',linewidth=0.8)
mon = mlines.Line2D([], [], color='crimson',    label='Summer',linewidth=0.8)
aut = mlines.Line2D([], [], color='navy',       label='Autumn',linewidth=0.8)
ann = mlines.Line2D([], [], color='black',      label='Annual',linewidth=0.8)

ax3.legend(handles=[win,spr,mon,aut,ann],fontsize = 5)

#=== plot BC 1d std====================================
ax3.fill_between(dust_win_sd.index, dust_win_lower, dust_win_upper, color='lightblue', alpha=0.8)
ax3.fill_between(dust_spr_sd.index, dust_spr_lower, dust_spr_upper, color='lightpink', alpha=0.5)
ax3.fill_between(dust_mon_sd.index, dust_mon_lower, dust_mon_upper, color='fuchsia', alpha=0.1)
ax3.fill_between(dust_aut_sd.index, dust_aut_lower, dust_aut_upper, color='blue', alpha=0.2)
ax3.fill_between(dust_all_sd.index, dust_all_lower, dust_all_upper, color='darkgray', alpha=0.8)

#==== highlight desired region  =====
ax3.axvspan(45, 80,  color='silver', alpha=0.3)  # Arica, Taklimaha, Rajasthan desert region
ax3.axvspan(93, 105, color='silver', alpha=0.3)  # Gobi desrt region
#=== add text =========================
ax3.text(60,50,"Contribution from... \n Sahara desert + \n Rajasthan desert + \n Taklamakan desert",fontsize=5,rotation='vertical')
ax3.text(95,39,"Contribution from \n Gobi desert",fontsize=5,rotation='vertical')

#=== setting axis===========================================================
ax3.tick_params(axis='y', labelsize=5)
ax3.tick_params(axis='x', labelsize=5)
gvutil.add_major_minor_ticks(ax3,x_minor_per_major=4,y_minor_per_major=3,labelsize=5)
#==== Label handle ====
gvutil.set_axes_limits_and_ticks(ax3, xlim=(45,120),xticklabels=['40E','50E','60E','70E','80E','90E','100E','110E','120'])
gvutil.set_axes_limits_and_ticks(ax3, ylim=(0,80))
gvutil.set_titles_and_labels(ax3,ylabel='DUST_IAT_ano \n [$10^{-5}$ $Kg  m^{-1} s^{-1}$]',labelfontsize=5)


fig.subplots_adjust(top=0.95,
                        bottom=0.11,
                        left=0.075,
                        right=0.98,
                        hspace=0.06,
                        wspace=0.36)




plt.savefig("/mnt/g/2nd_Paper/ipart/ivt_ZONAL_mean_2000dpi.png",dpi=2000)
#plt.show()
