#========== code for plotting aeronet-AOD vs. WRF-Chem modelled ======
import matplotlib
#matplotlib.use('TkAgg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from sklearn.metrics import r2_score
import cmaps
import matplotlib.gridspec as gridspec
from numpy.polynomial.polynomial import polyfit
from sklearn.impute import KNNImputer
import statsmodels.api as sm
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
from matplotlib import rcParams
#=== Setting theme from thempy===========
#import themepy
#import matplotlib
#matplotlib.use('Agg')
"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#============================================================
fig,ax = plt.subplots(ncols=2,figsize=(8,7)) #sharey=True,sharex=True)

pm10 = pd.read_excel(r'/mnt/g/Paper_work/PM10_Model_Observation.xlsx',sheet_name='Sheet1')
pm25   = pd.read_excel(r'/mnt/g/Paper_work/PM2.5_Model_Observation.xlsx',sheet_name='Sheet1')
#============== plot PM10 data ===============================
scatter1 = ax[0].scatter(pm10.Shijiazhuang_mo, pm10.Shijiazhuang_obs, color= 'm', alpha=0.4,s=40,edgecolors="black")            #,alpha=0.5,c=colors,s= area)
scatter2 = ax[0].scatter(pm10.Changsha_mo, pm10.Changsha_obs, color ='forestgreen', alpha=0.5,s=40,edgecolors="black")              #,alpha=0.8,c=colors,s= area)
scatter3 = ax[0].scatter(pm10.Xian_mo, pm10.Xian_obs, color ='lightcoral',alpha=0.6,s=50,edgecolors="black")
scatter4 = ax[0].scatter(pm10.Shigatse_mo, pm10.Shigatse_obs, color='tab:olive', alpha=0.8,s=40,edgecolors="black")       #cmap=cmaps.WhiteBlueGreenYellowRed,alpha=0.$
scatter5 = ax[0].scatter(pm10.Beijing_mo, pm10.Shigatse_obs, color ='darkorange', alpha=0.5,s=40,edgecolors="black")         #cmap=cmaps.WhiteBlueGreenYellowRed) # ,alpha=0$
scatter6 = ax[0].scatter(pm10.Lanzhou_mo, pm10.Lanzhou_obs, color ='steelblue', alpha=0.7,s=40,edgecolors="black")         #alpha=0.5,c=colors,s= area)
scatter7 = ax[0].scatter(pm10.Urumqi_mo, pm10.Urumqi_obs, color ='cyan', alpha=0.5,s=40,edgecolors="black")         #cmap=cmaps.WhiteBlueGreenYellowRed) # ,alpha=$
scatter8 = ax[0].scatter(pm10.Kunming_mo, pm10.Kunming_obs, color ='yellow', alpha=0.5,s=40,edgecolors="black")         #alpha=0.5,c=colors,s= area)

#==== Add 2:1 and 1:2 ratio line ===================
ax[0].plot(pm10.Beijing_obs, pm10.Beijing_obs,  linestyle= '-',    linewidth=1,color='mediumseagreen')
ax[0].plot(pm10.Beijing_obs, 2*pm10.Beijing_obs,linestyle= '-',   linewidth=1,color='tomato')
ax[0].plot(2*pm10.Beijing_obs, pm10.Beijing_obs,linestyle= '-',   linewidth=1,color='tomato')
##==== Legeend Propreties =============
legend1 = ax[0].legend((scatter1, scatter2, scatter3,scatter4, scatter5, scatter6,scatter7, scatter8),
('Shijiazhuang [$R^2$=0.53]','Changsha [$R^2$=0.55]','Xian [$R^2$=0.52]','Shigatse [$R^2$=0.74]','Beijing [$R^2$=0.59]','Lanzhou [$R^2$=0.30]','Urumqi [$R^2$=0.60]','Kunming [$R^2$=0.36]'),
           scatterpoints=1,
           loc='upper left',
           ncol=2,
           fontsize=8)
ax[0].add_artist(legend1)

ax[0].set_xlabel("Measured PM10 $[μgm^{-3}]$",fontsize=10)
ax[0].set_ylabel("Model PM10 $[μgm^{-3}]$",fontsize=10)
##=============== Plotting PM2.5 data =====================================
scatter9 = ax[1].scatter(pm25.Shijiazhuang_mo, pm25.Shijiazhuang_obs, color= 'm', alpha=0.4,s=40,edgecolors="black")            #,alpha=0.5,c=colors,s= area)
scatter10 = ax[1].scatter(pm25.Changsha_mo, pm25.Changsha_obs, color ='forestgreen', alpha=0.5,s=40,edgecolors="black")              #,alpha=0.8,c=colors,s= area)
scatter11 = ax[1].scatter(pm25.Xian_mo, pm25.Xian_obs, color ='lightcoral',alpha=0.6,s=50,edgecolors="black")
scatter12= ax[1].scatter(pm25.Shigatse_mo, pm25.Shigatse_obs, color='tab:olive', alpha=0.8,s=40,edgecolors="black")       #cmap=cmaps.WhiteBlueGreenYellowRed,alpha=0$
scatter13 = ax[1].scatter(pm25.Beijing_mo, pm25.Beijing_obs, color ='darkorange', alpha=0.5,s=40,edgecolors="black")         #cmap=cmaps.WhiteBlueGreenYellowRed) # ,alpha=0$
scatter14 = ax[1].scatter(pm25.Lanzhou_mo, pm25.Lanzhou_obs, color ='steelblue', alpha=0.7,s=40,edgecolors="black")  
scatter15 = ax[1].scatter(pm25.Urumqi_mo, pm25.Urumqi_obs, color ='cyan', alpha=0.5,s=40,edgecolors="black")         #cmap=cmaps.WhiteBlueGreenYellowRed) # ,alpha=$
scatter16 = ax[1].scatter(pm25.Kunming_mo, pm25.Kunming_obs, color ='yellow', alpha=0.5,s=40,edgecolors="black")

#==== Add 2:1 and 1:2 ratio line ===================
ax[1].plot(pm25.Beijing_obs, pm25.Beijing_obs,linestyle= '-',    linewidth=1,color='mediumseagreen')
ax[1].plot(pm25.Beijing_obs, 2*pm25.Beijing_obs,linestyle= '-', linewidth=1,color='tomato')
ax[1].plot(2*pm25.Beijing_obs, pm25.Beijing_obs,linestyle= '-', linewidth=1,color='tomato')

##===== Decorate legend ====
legend2 = plt.legend((scatter9, scatter10, scatter11,scatter12,scatter13, scatter14, scatter15,scatter16),
('Shijiazhuang [$R^2$=0.65]','Changsha [$R^2$=0.64]','Xian [$R^2$=0.57]','Shigatse [$R^2$=0.64]','Beijing [$R^2$=0.57]','Lanzhou [$R^2$=0.66]','Urumqi [$R^2$=0.52]','Kunming [$R^2$=0.57]'),
                    scatterpoints=1,
                    loc='upper left',
                    ncol=2,
                    fontsize=8)

#============= # Set common labels ############
#fig.text(0.5, 0.06, 'Measured [PM2.5]', ha='center', va='center',fontsize=15)
#fig.text(0.08, 0.5, 'Model [PM2.5]', va='center', rotation='vertical',fontsize=15)
#============================================
ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
#===========================================
ax[0].set_facecolor("silver")        # lightslategray")
ax[1].set_facecolor("silver")
#==========================================
#ax[0].set_xlim(left=0)
#ax[0].set_xlim(0,3)
#ax[0].set_ylim(bottom=0)
#ax[1].set_xlim(left=0)
#ax[1].set_ylim(bottom=0)

ax[1].set_xlim(0,1000)
ax[1].set_ylim(0,1000)
ax[0].set_xlim(0,1000)
ax[0].set_ylim(0,1000)
ax[1].set_xlabel("Measured PM2.5 $[μgm^{-3}]$",fontsize=10)
ax[1].set_ylabel("Model PM2.5 $[μgm^{-3}]$",fontsize=10)

#==== Figure adjust

fig.subplots_adjust(top=0.985,
			bottom=0.08,
			left=0.08,
			right=0.95,
			hspace=0.185,
			wspace=0.25)



#plt.savefig('/mnt/g/2nd_Paper/pm_correl_800dpi.png',dpi=800)
plt.show()
