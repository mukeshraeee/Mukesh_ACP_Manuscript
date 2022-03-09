import os
import numpy as np
from ipart.utils import funcs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from ipart.utils import plot
#==== Import data => contains aerosols data from MERRA-2 =================
data = os.path.join('.', 'winter1.nc')

#==== Read in data ===============
uflux=funcs.readNC(data, 'OCFLUXU')
vflux=funcs.readNC(data, 'OCFLUXV')
output=os.path.join('.', 'winter_ivt_oc.nc') #save nc file in given name
#=== Get data info ==============
#print(uflux.info())
#print(vflux.info())
#print(uflux.getLatitude())
#print(uflux.getLongitude())
#print(uflux.getTime())

#=== compute Aerosol IVT ==================
idt=np.ma.sqrt(uflux.data*uflux.data+vflux.data*vflux.data)
idt=funcs.NCVAR(idt, 'iot', uflux.axislist, {'name': 'iot', 'long_name': 'integrated OC transport (IOCT)',
                                            'standard_name': 'integrated_OC_transport',
                                            'title': 'integrated BC transport (IOCT)',
                                            'units': getattr(uflux, 'units', '')})
funcs.saveNC(output, idt)
