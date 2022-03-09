import os
import numpy as np
from ipart.utils import funcs
from ipart import thr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from ipart.utils import plot


#=== Read data ============
ivt_file = os.path.join('.','autumn_ivt_oc.nc')
var      = 'iot'

SHIFT_LON=10

KERNEL=[1,20,20]
OUTPUTDIR=os.path.abspath('.')
var=funcs.readNC(ivt_file, 'iot')
lon=var.getLongitude()
var=var.shiftLon(SHIFT_LON)
print('var.shape=', var.shape)

ivt, ivtrec, ivtano=thr.THR(var, KERNEL)

#--------Save------------------------------------
if not os.path.exists(OUTPUTDIR):
    os.makedirs(OUTPUTDIR)

fname=os.path.split(ivt_file)[1]
file_out_name='%s-THR-kernel-t%d-s%d.nc'\
        %(os.path.splitext(fname)[0], KERNEL[0], KERNEL[1])

abpath_out=os.path.join(OUTPUTDIR,file_out_name)
print('\n# Saving output to:\n',abpath_out)
funcs.saveNC(abpath_out, ivt, 'w')
funcs.saveNC(abpath_out, ivtrec, 'a')
funcs.saveNC(abpath_out, ivtano, 'a')
