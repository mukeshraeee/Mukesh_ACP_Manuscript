#=========== This code produce PM10 and PM2.5 plot ====================
import numpy as np
import pandas as pd
from matplotlib import cm
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.dates as mdates
import datetime as dt
from matplotlib.lines import Line2D
import mplcyberpunk
from matplotlib import rcParams
#plt.style.use("cyberpunk")
#mplcyberpunk.add_glow_effects()
#mplcyberpunk.make_lines_glow()
"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#======================================================================
pm1    = pd.read_excel(r'/mnt/g/Paper_work/PM10_Model_Observation.xlsx',sheet_name='Sheet1')
pm2    = pd.read_excel(r'/mnt/g/Paper_work/PM2.5_Model_Observation.xlsx',sheet_name='Sheet1')
#===== Removing NAN values =============================================================================
pm10 = pm1.dropna()
pm25 = pm2.dropna()
#===================================================
fig, ax  =  plt.subplots(nrows=8,ncols=2, figsize=(4,5),dpi=300)
fig.autofmt_xdate()
#================================ Plot PM10 ======================
ax[0,0].plot(pm10.Date,pm10.Shijiazhuang_mo, color='tab:cyan' ,linewidth=0.4,linestyle='--')
ax[0,0].plot(pm10.Date,pm10.Shijiazhuang_obs, color='darkorange' ,linewidth=0.4)
#ax[0,0].set_title(r'$PM10\ [μgm^{-3}]$',fontsize=6)
#ax[0,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[0,0].set_ylabel('Shijiazhuang',fontsize=6)
ax[0,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[0,0].tick_params(axis='y', labelsize= 5)
ax[0,0].get_yaxis().set_label_coords(-0.2,.5)
ax[0,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[0,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
#ax[0,0].set_yticks([0,500,1000])
ax[0,0].tick_params(bottom=True, top=False, left=True, right=False)
ax[0,0].set_yticks([0,100,200,300,400,500])
ax[0,0].set_ylim([0,500])




ax[1,0].plot(pm10.Date,pm10.Changsha_mo, color='tab:olive' ,label='Model',linewidth=0.4,linestyle='--')
ax[1,0].plot(pm10.Date,pm10.Changsha_obs, color='tab:gray' ,label='Observation',linewidth=0.4)
#ax[1,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[1,0].set_ylabel('Changsha',fontsize=6)
ax[1,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[1,0].tick_params(axis='y', labelsize= 5)
ax[1,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[1,0].set_yticks([0,500,1000])
ax[1,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[1,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[1,0].set_yticks([0,100,200,300,400,500])
ax[1,0].set_ylim([0,500])



ax[2,0].plot(pm10.Date,pm10.Beijing_mo, color='tab:green' ,label='Model',linewidth=0.4,linestyle='--')
ax[2,0].plot(pm10.Date,pm10.Beijing_obs, color='tab:orange' ,label='Observation',linewidth=0.4)
#ax[2,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[2,0].set_ylabel('Beijing',fontsize=6)
ax[2,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[2,0].tick_params(axis='y', labelsize= 5)
ax[2,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[2,0].set_yticks([0,500,800])
ax[2,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[2,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[2,0].set_yticks([0,100,200,300,400,500])
ax[2,0].set_ylim([0,500])


ax[3,0].plot(pm10.Date,pm10.Kunming_mo, color='steelblue' ,label='Model',linewidth=0.4,linestyle='--')
ax[3,0].plot(pm10.Date,pm10.Kunming_obs, color='deeppink' ,label='Observation',linewidth=0.4)
#ax[3,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[3,0].set_ylabel('Kunming',fontsize=6)
ax[3,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[3,0].tick_params(axis='y', labelsize= 5)
ax[3,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[3,0].set_yticks([100,500,1000])
ax[3,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[3,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[3,0].set_yticks([0,100,200,300,400,500])
ax[3,0].set_ylim([0,500])


ax[4,0].plot(pm10.Date,pm10.Xian_mo, color='salmon' ,label='Model',linewidth=0.4,linestyle='--')
ax[4,0].plot(pm10.Date,pm10.Xian_obs, color='green' ,label='Observation',linewidth=0.4)
#ax[4,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[4,0].set_ylabel('Xian',fontsize=6)
ax[4,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[4,0].tick_params(axis='y', labelsize= 5)
ax[4,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[4,0].set_yticks([0,100,1000])
ax[4,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[4,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[4,0].set_yticks([0,100,200,300,400,500])
ax[4,0].set_ylim([0,500])




ax[5,0].plot(pm10.Date,pm10.Urumqi_mo, color='indianred' ,label='Model',linewidth=0.4,linestyle='--')
ax[5,0].plot(pm10.Date,pm10.Urumqi_obs, color='teal' ,label='Observation',linewidth=0.4)
#ax[5,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[5,0].set_ylabel('Urumqi',fontsize=6)
ax[5,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[5,0].tick_params(axis='y', labelsize= 5)
ax[5,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[5,0].set_yticks([0,100,1000])
ax[5,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[5,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[5,0].set_yticks([0,100,200,300,400,500])
ax[5,0].set_ylim([0,500])


ax[6,0].plot(pm10.Date,pm10.Shigatse_mo, color='cadetblue' ,label='Model',linewidth=0.4,linestyle='--')
ax[6,0].plot(pm10.Date,pm10.Shigatse_obs, color='yellowgreen' ,label='Observation',linewidth=0.4)
#ax[6,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[6,0].set_ylabel('Shigatse',fontsize=6)
ax[6,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[6,0].tick_params(axis='y', labelsize= 5)
ax[6,0].get_yaxis().set_label_coords(-0.2,0.5)
#ax[6,0].set_yticks([0,100,1000])
ax[6,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[6,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)
ax[6,0].set_yticks([0,100,200,300,400,500])
ax[6,0].set_ylim([0,500])


ax[7,0].plot(pm10.Date,pm10.Lanzhou_mo, color='m' ,label='Model',linewidth=0.4,linestyle='--')
ax[7,0].plot(pm10.Date,pm10.Lanzhou_obs, color='c' ,label='Observation',linewidth=0.4)
#ax[7,0].legend(loc='upper center',bbox_to_anchor=(1.1, 0.8), prop={"size":4})
ax[7,0].set_ylabel('Lanzhou',fontsize=6)
ax[7,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[7,0].xaxis.set_tick_params(labelsize=5)
ax[7,0].tick_params(axis='y', labelsize= 5)
ax[7,0].get_yaxis().set_label_coords(-0.2,0.5)
ax[7,0].axhline(y=50, color='r', linestyle='-',linewidth=0.4)
ax[7,0].set_yticks([0,100,200,300,400,500])
ax[7,0].set_ylim([0,500])
#ax[7,0].set(xlabel="Year-Month",labelsize=8)
#ax[7,0].set_xlabel("Year-Month",fontsize=8)
ax[7,0].set_xlabel("$PM10\ [μgm^{-3}]$", fontsize=5)

ax[7,0].axhline(y=150, color='r', linestyle='--',linewidth=0.4)


#================================ Plot PM2.5 ======================
ax[0,1].plot(pm25.Date,pm25.Shijiazhuang_mo, color='tab:cyan' ,label='Model',linewidth=0.4,linestyle='--')
ax[0,1].plot(pm25.Date,pm25.Shijiazhuang_obs, color='darkorange' ,label='Observation',linewidth=0.4)
#ax[0,1].set_title(r'$PM2.5\ [μgm^{-3}]$',fontsize=8)
#ax[0,1].legend(loc='upper center', prop={"size":4})
ax[0,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[0,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[0,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[0,1].tick_params(axis='y', labelsize= 5)
#ax[0,1].set_yticks([0,500])
ax[0,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[0,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[0,1].set_yticks([0,50,100,150,200,250])
ax[0,1].set_ylim([0,250])



ax[1,1].plot(pm25.Date,pm25.Changsha_mo, color='tab:olive' ,label='Model',linewidth=0.4,linestyle='--')
ax[1,1].plot(pm25.Date,pm25.Changsha_obs, color='tab:gray' ,label='Observation',linewidth=0.4)
#ax[1,1].legend(loc='upper center', prop={"size":4})
ax[1,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[1,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[1,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[1,1].tick_params(axis='y', labelsize= 5)
#ax[1,1].set_yticks([0,500])
ax[1,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[1,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[1,1].set_yticks([0,50,100,150,200,250])
ax[1,1].set_ylim([0,250])


ax[2,1].plot(pm25.Date,pm25.Beijing_mo, color='tab:green' ,label='Model',linewidth=0.4,linestyle='--')
ax[2,1].plot(pm25.Date,pm25.Beijing_obs, color='tab:orange' ,label='Observation',linewidth=0.4)
#ax[2,1].legend(loc='upper center', prop={"size":4})
ax[2,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[2,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[2,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[2,1].tick_params(axis='y', labelsize= 5)
#ax[2,1].set_yticks([0,500])
ax[2,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[2,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[2,1].set_yticks([0,50,100,150,200,250])
ax[2,1].set_ylim([0,250])


ax[3,1].plot(pm25.Date,pm25.Kunming_mo, color='steelblue' ,label='Model',linewidth=0.4,linestyle='--')
ax[3,1].plot(pm25.Date,pm25.Kunming_obs, color='deeppink' ,label='Observation',linewidth=0.4)
#ax[3,1].legend(loc='upper center', prop={"size":4})
ax[3,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[3,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[3,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[3,1].tick_params(axis='y', labelsize= 5)
#ax[3,1].set_yticks([0,500])
ax[3,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[3,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[3,1].set_yticks([0,50,100,150,200,250])
ax[3,1].set_ylim([0,250])


ax[4,1].plot(pm25.Date,pm25.Xian_mo, color='salmon' ,label='Model',linewidth=0.4,linestyle='--')
ax[4,1].plot(pm25.Date,pm25.Xian_obs, color='green' ,label='Observation',linewidth=0.4)
#ax[4,1].legend(loc='upper center', prop={"size":4})
ax[4,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[4,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[4,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[4,1].tick_params(axis='y', labelsize= 5)
#ax[4,1].set_yticks([0,500])
ax[4,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[4,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[4,1].set_yticks([0,50,100,150,200,250])
ax[4,1].set_ylim([0,250])


ax[5,1].plot(pm25.Date,pm25.Urumqi_mo, color='indianred' ,label='Model',linewidth=0.4,linestyle='--')
ax[5,1].plot(pm25.Date,pm25.Urumqi_obs, color='teal' ,label='Observation',linewidth=0.4)
#ax[5,1].legend(loc='upper center', prop={"size":4})
ax[5,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[5,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[5,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[5,1].tick_params(axis='y', labelsize= 5)
#ax[5,1].set_yticks([0,500])
ax[5,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[5,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[5,1].set_yticks([0,50,100,150,200,250])
ax[5,1].set_ylim([0,250])


ax[6,1].plot(pm25.Date,pm25.Shigatse_mo, color='cadetblue' ,label='Model',linewidth=0.4,linestyle='--')
ax[6,1].plot(pm25.Date,pm25.Shigatse_obs, color='yellowgreen' ,label='Observation',linewidth=0.4)
#ax[6,1].legend(loc='upper center', prop={"size":4})
ax[6,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[6,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[6,1].tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=True)
ax[6,1].tick_params(axis='y', labelsize= 5)
#ax[6,1].set_yticks([0,500])
ax[6,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[6,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[6,1].set_yticks([0,50,100,150,200,250])
ax[6,1].set_ylim([0,250])


ax[7,1].plot(pm25.Date,pm25.Lanzhou_mo, color='m' ,label='Model',linewidth=0.4,linestyle='--')
ax[7,1].plot(pm25.Date,pm25.Lanzhou_obs, color='c' ,label='Observation',linewidth=0.4)
#ax[7,1].legend(loc='upper center', prop={"size":4})
ax[7,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[7,1].xaxis.set_tick_params(labelsize=5)
ax[7,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[7,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[7,1].tick_params(axis='y', labelsize= 5)
#ax[7,1].set_yticks([0,500])
#ax[7,1].set(xlabel="Year-Month", labelsize=8)
ax[7,1].set_xlabel("$PM2.5\ [μgm^{-3}]$", fontsize=5)
ax[7,1].axhline(y=25, color='blue', linestyle='-',linewidth=0.4)
ax[7,1].axhline(y=75, color='blue', linestyle='--',linewidth=0.4)
ax[7,1].set_yticks([0,50,100,150,200,250])
ax[7,1].set_ylim([0,250])


#==== Now custom legend=============
custom_lines = [Line2D([0], [0], color='r', lw=0.4,linestyle='-'),
                Line2D([0], [0], color='r', lw=0.4,linestyle='--'),
                Line2D([0], [0], color='blue', lw=0.4,linestyle='-'),
                Line2D([0], [0], color='blue', lw=0.4,linestyle='--')]


plt.legend(custom_lines, ['PM10-WHO', 'PM10-China Grade[II]', 'PM2.5-WHO','PM2.5-China Grade[II]'],bbox_to_anchor =(0.95,9.9),ncol = 4,fontsize=4.5,edgecolor='gray')


#======= Seasonal highlights=======================
#=====  Winter========
ax[0,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[1,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[2,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[3,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[4,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[5,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[6,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[7,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
#== Add December month
ax[0,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[1,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[2,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[3,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[4,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[5,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[6,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[7,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
#===== Spring==========
ax[0,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[1,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[2,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[3,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[4,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[5,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[6,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[7,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
#======== Summer =============
ax[0,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[1,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[2,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[3,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[4,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[5,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[6,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[7,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
#=== Autumn========
ax[0,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[1,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[2,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[3,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[4,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[5,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[6,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[7,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)

#==== PM2.5=============
#=====  Winter========
ax[0,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[1,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[2,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[3,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[4,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[5,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[6,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[7,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
#== Add December month
ax[0,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[1,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[2,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[3,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[4,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[5,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[6,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[7,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
#===== Spring==========
ax[0,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[1,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[2,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[3,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[4,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[5,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[6,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[7,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
#======== Summer =============
ax[0,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[1,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[2,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[3,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[4,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[5,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[6,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[7,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
#=== Autumn========
ax[0,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[1,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[2,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[3,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[4,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[5,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[6,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[7,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)



fig.subplots_adjust(top=0.96,
                        bottom=0.11,
                        left=0.115,
                        right=0.915,
                        hspace=0.21,
                        wspace=0.06)





#plt.savefig("/mnt/g/2nd_Paper/modelVsobs_pm_800dpi.png",dpi=800)
plt.show()
