##================ import require libraries ==============
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
from matplotlib.legend import Legend
from matplotlib import rcParams
#===== Use neon background ======
#import mplcyberpunk
#plt.style.use("cyberpunk")
#mplcyberpunk.add_glow_effects()
#mplcyberpunk.make_lines_glow()
#mplcyberpunk.add_underglow()
"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"


#=================================  Load data set ========================================
#================================= AERONET DATA =======================================================================
aero  = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Aeronet_Beij_Dhaka_Kanpu_Mezaira.xlsx',sheet_name='Sheet1')
aero_karachi = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Karachi.xlsx',sheet_name='Sheet1')
aero_lumbini = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Lumbini_lev2.xlsx',sheet_name='Sheet1')
aero_langtang = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Langtang.xlsx',sheet_name='Sheet1')
aero_namco = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Namco.xlsx',sheet_name='Sheet1')
aero_qoms = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/QOMS.xlsx',sheet_name='Sheet1')
aero_dushambe = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/Aeronet/Dushanbe.xlsx',sheet_name='Sheet1')
#================================= WRF DATA ============================================================== 
wrf = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/WRF/wrf_aod.xlsx',sheet_name='Sheet1')
#================================= CAMS DATA ===============================================
cams = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/CAMS/cams_aod_new.xlsx',sheet_name='Sheet1')
#================================= MERRA DATA =====================================================
merra = pd.read_excel(r'/mnt/g/Paper_work/AOD_plot/MERRA/Merra_aod_all_month.xlsx',sheet_name='Sheet1')
#================================ Now plot=====================================================================
fig, ax  =  plt.subplots(nrows=5,ncols=2,figsize=(4,4),dpi=300)
fig.autofmt_xdate()
#======================= For Urban loacation  ======================
#==== Beijing ====================================
#ax[0,0].scatter(aero.date,aero.Beijing, color='gray' ,label='Aero',s=1) #,linewidth=0.7,linestyle='dotted')
ax[0,0].plot(aero.date,aero.Beijing, color='gray' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[0,0].plot(wrf.date,wrf.Beijing,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[0,0].plot(cams.date,cams.Beijing, color='coral' ,label='CAMS',linewidth=0.5)
ax[0,0].plot(merra.date,merra.Beijing,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[0,0].set_ylabel('Beijing',fontsize=5)
#==== Dhaka  ====================================
ax[1,0].plot(aero.date,aero.Dhaka, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[1,0].plot(wrf.date,wrf.Dhaka,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[1,0].plot(cams.date,cams.Dhaka, color='coral' ,label='CAMS',linewidth=0.5)
ax[1,0].plot(merra.date,merra.Dhaka,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[1,0].set_ylabel('Dhaka',fontsize=5)
#==== Kanpur ====================================
ax[2,0].plot(aero.date,aero.Kanpur, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[2,0].plot(wrf.date,wrf.Kanpur,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[2,0].plot(cams.date,cams.Kanpur, color='coral' ,label='CAMS',linewidth=0.5)
ax[2,0].plot(merra.date,merra.Kanpur,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[2,0].set_ylabel('Kanpur',fontsize=5)
#==== Lumbini ====================================
ax[3,0].plot(aero_lumbini.date,aero_lumbini.Lumbini, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[3,0].plot(wrf.date,wrf.Lumbini,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[3,0].plot(cams.date,cams.Lumbini, color='coral' ,label='CAMS',linewidth=0.5)
ax[3,0].plot(merra.date,merra.Lumbini,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[3,0].set_ylabel('Lumbini',fontsize=5)
#==== Karachi ====================================
ax[4,0].plot(aero_karachi.date,aero_karachi.Karachi, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[4,0].plot(wrf.date,wrf.Karachi,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[4,0].plot(cams.date,cams.Karachi, color='coral' ,label='CAMS',linewidth=0.5)
ax[4,0].plot(merra.date,merra.Karachi,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[4,0].set_ylabel('Karachi',fontsize=5)
#============================= For Remote locations =======================
#==== Langtang ====================================
ax[0,1].plot(aero_langtang.date,aero_langtang.Langtang, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[0,1].plot(wrf.date,wrf.Langtang,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[0,1].plot(cams.date,cams.Langtang, color='coral' ,label='CAMS',linewidth=0.5)
ax[0,1].plot(merra.date,merra.Langtang,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[0,1].set_ylabel('Langtang',fontsize=5)
#==== QOMS  ====================================
ax[1,1].plot(aero_qoms.date,aero_qoms.qoms, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[1,1].plot(wrf.date,wrf.Qoms,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[1,1].plot(cams.date,cams.Qoms, color='coral' ,label='CAMS',linewidth=0.5)
ax[1,1].plot(merra.date,merra.Qoms,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[1,1].set_ylabel('Qoms',fontsize=5)
#==== Namco ====================================
ax[2,1].plot(aero_namco.date,aero_namco.Namco, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[2,1].plot(wrf.date,wrf.Namco,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[2,1].plot(cams.date,cams.Namco, color='coral' ,label='CAMS',linewidth=0.5)
ax[2,1].plot(merra.date,merra.Namco,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[2,1].set_ylabel('Namco',fontsize=5)
#==== Mezaira ====================================
ax[3,1].plot(aero.date,aero.Mezaira, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[3,1].plot(wrf.date,wrf.Mezaira,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[3,1].plot(cams.date,cams.Mezaira, color='coral' ,label='CAMS',linewidth=0.5)
ax[3,1].plot(merra.date,merra.Mezaira,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[3,1].set_ylabel('Mezaira',fontsize=5)
#==== Dushambe ====================================
ax[4,1].plot(aero_dushambe.date,aero_dushambe.Dushambe, color='grey' ,label='Aero',linewidth=0.5,linestyle='dotted')
ax[4,1].plot(wrf.date,wrf.Dushambe,  color='palevioletred' ,label='WRF',linewidth=0.5,linestyle='dashdot')
ax[4,1].plot(cams.date,cams.Dushambe, color='coral' ,label='CAMS',linewidth=0.5)
ax[4,1].plot(merra.date,merra.Dushambe,  color='darkcyan' ,label='MERRA',linewidth=0.5)
ax[4,1].set_ylabel('Dushanbe',fontsize=5)


#====== Hadling xtick label =====================
ax[0,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[1,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[2,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[3,0].set(xlim=["2017-01-01", "2017-12-31"])
ax[4,0].set(xlim=["2017-01-01", "2017-12-31"])

ax[0,1].set(ylim=[0,1.5])
ax[1,1].set(ylim=[0,1.5])
ax[2,1].set(ylim=[0,1.5])
ax[3,1].set(ylim=[0,1.5])
ax[4,1].set(ylim=[0,1.5])

ax[0,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[1,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[2,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[3,1].set(xlim=["2017-01-01", "2017-12-31"])
ax[4,1].set(xlim=["2017-01-01", "2017-12-31"])

ax[0,0].set(ylim=[0,3])
ax[1,0].set(ylim=[0,3])
ax[2,0].set(ylim=[0,3])
ax[3,0].set(ylim=[0,3])
ax[4,0].set(ylim=[0,3])


##===== Removing xtick ===============
ax[0,0].set_xticklabels([])
ax[1,0].set_xticklabels([])
ax[2,0].set_xticklabels([])
ax[3,0].set_xticklabels([])

ax[0,1].set_xticklabels([])
ax[1,1].set_xticklabels([])
ax[2,1].set_xticklabels([])
ax[3,1].set_xticklabels([])

##===== Placing xtick to right side of the plots=========================================
ax[0,0].tick_params(bottom=True, top=False, left=False, right=True)
ax[0,0].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[1,0].tick_params(bottom=True, top=False, left=False, right=True)
ax[1,0].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[2,0].tick_params(bottom=True, top=False, left=False, right=True)
ax[2,0].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[3,0].tick_params(bottom=True, top=False, left=False, right=True)
ax[3,0].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[4,0].tick_params(bottom=True, top=False, left=False, right=True)
ax[4,0].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)

ax[0,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[0,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[1,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[1,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[2,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[2,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[3,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[3,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
ax[4,1].tick_params(bottom=True, top=False, left=False, right=True)
ax[4,1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)

#===== Make shadow-Seasonal ===========
#== Winter
ax[0,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[1,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[2,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[3,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[4,0].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
#== Add December month
ax[0,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[1,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[2,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[3,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[4,0].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
#== Spring
ax[0,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[1,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[2,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[3,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[4,0].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
#==== Summer  =========
ax[0,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[1,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[2,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[3,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[4,0].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
#==== Autumn =========
ax[0,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[1,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[2,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[3,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[4,0].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)

#==== Make shadow for RIGHT coloumn ======
#== Winter
ax[0,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[1,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[2,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[3,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
ax[4,1].axvspan(*mdates.datestr2num(['01/01/2017', '02/28/2017']), color='silver', alpha=0.5)
#== Add December month
ax[0,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[1,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[2,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[3,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
ax[4,1].axvspan(*mdates.datestr2num(['12/01/2017', '12/30/2017']), color='silver', alpha=0.5)
#== Spring
ax[0,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[1,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[2,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[3,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
ax[4,1].axvspan(*mdates.datestr2num(['03/01/2017', '05/31/2017']), color='salmon', alpha=0.2)
#==== Summer  =========
ax[0,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[1,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[2,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[3,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
ax[4,1].axvspan(*mdates.datestr2num(['06/01/2017', '08/31/2017']), color='darkturquoise', alpha=0.2)
#==== Autumn =========
ax[0,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[1,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[2,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[3,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)
ax[4,1].axvspan(*mdates.datestr2num(['09/01/2017', '11/30/2017']), color='fuchsia', alpha=0.1)






##==== Resize x and y lable 
ax[4,0].tick_params(axis='x', labelsize=5)
ax[4,1].tick_params(axis='x', labelsize=5)

ax[0,0].tick_params(axis='y', labelsize=5)
ax[1,0].tick_params(axis='y', labelsize=5)
ax[2,0].tick_params(axis='y', labelsize=5)
ax[3,0].tick_params(axis='y', labelsize=5)
ax[4,0].tick_params(axis='y', labelsize=5)
ax[0,1].tick_params(axis='y', labelsize=5)
ax[1,1].tick_params(axis='y', labelsize=5)
ax[2,1].tick_params(axis='y', labelsize=5)
ax[3,1].tick_params(axis='y', labelsize=5)
ax[4,1].tick_params(axis='y', labelsize=5)

# Legend
plt.legend(bbox_to_anchor =(0.7,6.15), ncol = 4,fontsize=5,edgecolor='gray')
fig.text(0.5, 0.04, 'Year-Month', ha='center', va='center',fontsize=5)

fig.subplots_adjust(top=0.925,
			bottom=0.135,
			left=0.075,
			right=0.935,
			hspace=0.17,
			wspace=0.235)

## Saving the plot
plt.savefig('/mnt/g/2nd_paper/aod_600dpi.png',dpi=600)

#plt.show()
