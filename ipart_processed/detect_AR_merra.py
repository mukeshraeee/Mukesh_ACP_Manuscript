"""This code enable to produce AR detection
  inspired by Guangzhi Xu (xugzhi1987@gmail.com)"""


#====== Load Libraries ======================================================
import os, sys
import numpy as np
import pandas as pd
from ipart.utils import funcs
from ipart.AR_detector import findARs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from ipart.utils import plot


YEAR=2021
TIME_START = '%d-03-01 00:00:00' %YEAR
TIME_END   = '%d-03-03 00:00:00' %YEAR

#==== Import data = =================
uq_file = os.path.join('.', 'spring1.nc')
vq_file = os.path.join('.', 'spring1.nc')

#==== ivt reconstruction and anomalies ===================
ivt_file = os.path.join('.', 'spring_ivt_dust-THR-kernel-t1-s20.nc') #save nc file in given name
#------------------Output folder------------------
output_dir = os.path.join('.', str(YEAR))


PLOT=True          # create maps of found ARs or not
SHIFT_LON=10          # degree, shift left bound to longitude. Should match
                      # that used in compute_thr_singlefile.py


PARAM_DICT={
    # kg/m/s, define AR candidates as regions >= than this anomalous ivt.
    'thres_low' : 0.5,
    # km^2, drop AR candidates smaller than this area.
    'min_area': 5*1e4,
    # km^2, drop AR candidates larger than this area.
    'max_area': 1800*1e4,
    # float, minimal length/width ratio.
    'min_LW': 2,
    # degree, exclude systems whose centroids are lower than this latitude.
    'min_lat': 00,
    # degree, exclude systems whose centroids are higher than this latitude.
    'max_lat': 80,
    # km, ARs shorter than this length is treated as relaxed.
    'min_length': 2000,
    # km, ARs shorter than this length is discarded.
    'min_length_hard': 1500,
    # degree lat/lon, error when simplifying axis using rdp algorithm.
    'rdp_thres': 2,
    # grids. Remove small holes in AR contour.
    'fill_radius': None,
    # do peak partition or not, used to separate systems that are merged
    # together with an outer contour.
    'single_dome': False,
    # max prominence/height ratio of a local peak. Only used when single_dome=True
    'max_ph_ratio': 0.6,
    # minimal proportion of flux component in a direction to total flux to
    # allow edge building in that direction
    'edge_eps': 0.4
    }


#==== Read in data ===============
uflux=funcs.readNC(uq_file, 'DUFLUXU')
vflux=funcs.readNC(vq_file, 'DUFLUXV')


#-------------------Read in ivt-------------------
print('\n# Read in file:\n',ivt_file)
ivt=funcs.readNC(ivt_file, 'idt')
ivtrec=funcs.readNC(ivt_file, 'ivt_rec')
ivtano=funcs.readNC(ivt_file, 'ivt_ano')
#-----------------Shift longitude-----------------
qu=uflux.shiftLon(SHIFT_LON)
qv=vflux.shiftLon(SHIFT_LON)
ivt=ivt.shiftLon(SHIFT_LON)
ivtrec=ivtrec.shiftLon(SHIFT_LON)
ivtano=ivtano.shiftLon(SHIFT_LON)
#--------------------Slice data--------------------
qu=uflux.sliceData(TIME_START,TIME_END,axis=0).squeeze()
qv=vflux.sliceData(TIME_START,TIME_END,axis=0).squeeze()
ivt=ivt.sliceData(TIME_START,TIME_END,axis=0).squeeze()
ivtrec=ivtrec.sliceData(TIME_START,TIME_END,axis=0).squeeze()
ivtano=ivtano.sliceData(TIME_START,TIME_END,axis=0).squeeze()
#-----------------Get coordinates-----------------
latax=qu.getLatitude()
lonax=qu.getLongitude()
timeax=ivt.getTime()
timeax=['%d-%02d-%02d %02d:00' %(timett.year, timett.month, timett.day, timett.hour) for timett in timeax]

time_idx, labels, angles, crossfluxes, result_df = findARs(ivt.data, ivtrec.data,
            ivtano.data, qu.data, qv.data, latax, lonax, times=timeax, **PARAM_DICT)


#=== Plot ===================

plot_idx=time_idx[0]

plot_time=timeax[plot_idx]
slab=ivt.data[plot_idx]
slabrec=ivtrec.data[plot_idx]
slabano=ivtano.data[plot_idx]
ardf=result_df[result_df.time==plot_time]

plot_vars=[slab, slabrec, slabano]
titles=['IVT', 'THR_recon', 'THR_ano']
iso=plot.Isofill(plot_vars, 12, 1, 1,min_level=0,max_level=800)

figure=plt.figure(figsize=(12,10), dpi=100)

for jj in range(len(plot_vars)):
    ax=figure.add_subplot(3,1,jj+1,projection=ccrs.PlateCarree())
    pobj=plot.plot2(plot_vars[jj], iso, ax,
            xarray=lonax, yarray=latax,
            title='%s %s' %(plot_time, titles[jj]),
            fix_aspect=False)
    plot.plotAR(ardf, ax, lonax)

plt.show()

