
""" This code generates study domain using Cartopy and etopo1 """
from copy import copy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sgeom
import cartopy
import cartopy.feature as cfe
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
from cartopy.io import shapereader as shpreader
import cartopy.io.img_tiles as cimgt
from shapely.geometry.polygon import LinearRing
import matplotlib.patches as mpatches
import pandas as pd
import matplotlib as mpl
import matplotlib
import pyproj
from cartopy.feature import NaturalEarthFeature, OCEAN, LAKES
from mpl_toolkits.basemap import Basemap
import shapely.geometry as sgeom
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.font_manager import FontProperties
import netCDF4 as nc
import re  # regular expression
import cmaps
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import iris
import iris.plot as iplt
import string
import matplotlib.cm as cm
from matplotlib import rcParams
import pygmt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#============ ticklabel custumization =======================================================
def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.

    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels
#=========================================
et    = '/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/etopo_PTP.nc'
etopo = nc.Dataset(et, 'r', format = 'NETCDF4')

dem     = etopo.variables['z'][:]
dem_lat = etopo.variables['y'][:]
dem_lon = etopo.variables['x'][:]

etopo.close()
dem_lons, dem_lats = np.meshgrid(dem_lon, dem_lat)

for i in np.arange(dem.shape[0]):
    for j in np.arange(dem.shape[1]):
        if dem[i,j]<0:
            dem[i,j]=0
#===================================================================================================================
"""Control font"""
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

proj = ccrs.LambertConformal(central_longitude=86.3333, central_latitude=27.5,standard_parallels=(27, 30))


"""# =========  Draw a set of axes with coastlines:========================"""

fig = plt.figure(figsize=(9, 5), dpi=200,frameon=True)

"""#=============================================================================="""
ax = fig.add_axes([0.09, 0.09, 0.9, 0.9], projection=proj)
ax.set_extent([30, 135, 0, 60],crs=ccrs.PlateCarree())
#ax.set_extent([30, 130, 5, 60],crs=ccrs.PlateCarree())

cs = ax.pcolormesh(dem_lons, dem_lats, dem, cmap=cmaps.GMT_relief_oceanonly, vmin=0, vmax=8848, alpha=0.9, transform=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.COASTLINE,linewidth=0.1)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-',linewidth=0.2,edgecolor='white')
ax.add_feature(cfeature.COASTLINE,linewidth=0.2,edgecolor='white')
ax.add_feature(cfeature.LAND,facecolor='dimgrey')
ax.add_feature(cfeature.OCEAN, color="white",linewidth=0.1)
ax.add_feature(OCEAN, edgecolor='white',linewidth=0.1) #,facecolor='deepskyblue')
ax.add_feature(LAKES, edgecolor='white',linewidth=0.1) #,facecolor='deepskyblue')
states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='grey',name='admin_1_states_provinces_shp',linewidth=0.2)
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='white', facecolor='#00BFFF'))
ax.add_feature(states, linewidth=0.1)
ax.coastlines(linewidth=0.3,color="white")

#===== add lines ========
"""---->>[lon , lon],[lat, lat]"""
plt.plot([25,   45], [55,  55], transform=ccrs.PlateCarree(),linestyle='--',color='orange', linewidth=1,zorder=10)
plt.plot([25,   45], [50,  10], transform=ccrs.PlateCarree(),linestyle='--',color='orange', linewidth=1,zorder=10)



#==== Add locations on the map=================
#met     = pd.read_excel('/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/location_metdata_aod_pm25.xlsx',sheet_name='met')
aer     = pd.read_excel('/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/location_metdata_aod_pm25.xlsx',sheet_name='aero')
pm25    = pd.read_excel('/mnt/g/WRF_All/WRF_Output/WRF_output_new/all_output/location_metdata_aod_pm25.xlsx',sheet_name='pm25')
#bc_data      = pd.read_excel('bc_measured_data.xlsx',sheet_name='Sheet1')

#plt.scatter(met.lon,met.lat,color=   'hotpink', alpha=0.9,s=50, edgecolors ='hotpink',transform=ccrs.PlateCarree())
plt.scatter(aer.lon,aer.lat,  color =  'darkgrey',  alpha=0.9,s=70, edgecolors ="yellow",transform=ccrs.PlateCarree())
plt.scatter(pm25.lon,pm25.lat,color = 'magenta',    alpha=0.9,s=60, edgecolors ="darkgray",transform=ccrs.PlateCarree())

# ======*must* call draw in order to get the axis boundary used to add ticks:=======================================
fig.canvas.draw()
#============ Define gridline locations and draw the lines using cartopy's built-in gridliner:#======================
xticks = [20,30, 40, 50,60,70,80,90,100,110,120,130,140]
yticks = [0, 10,20,30, 40,50,60]

s = ax.gridlines(xlocs=xticks, ylocs=yticks,linewidth=1, color='grey', alpha=0.4, linestyle='dashed')
s.xlabel_style = {'color': 'red'} #, 'weigh
s.xlabel_style = {'size': 0.5, 'color': 'black'}
s.ylabel_style = {'size': 0.5, 'color': 'black'}
# Label the end-points of the gridlines using the custom tick makers:
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
lambert_xticks(ax, xticks)
lambert_yticks(ax, yticks)

cax,kw = matplotlib.colorbar.make_axes(ax,location='right',pad=0.01,shrink=0.93)
cbar=fig.colorbar(cs,cax=cax,**kw)
cbar.set_label('Elevation [m]',size=10)
cbar.ax.tick_params(labelsize=8)

#========= Add plot ========
ax1 = plt.axes([0.087, .77, 0.2, 0.18]) ##left, bottom, width, height
m = Basemap(projection='robin', lat_0=0, lon_0=0,resolution='h',area_thresh = 1000.0)
#m.shadedrelief()
x, y = m(dem_lons, dem_lats)
ax1.pcolormesh(x, y, dem, cmap=cmaps.GMT_relief_oceanonly)
m.drawparallels(np.arange(-90.,120.,30.),dashes=[2, 2],linewidth=0.2)
m.drawmeridians(np.arange(0.,360.,60.),dashes=[2, 2],linewidth=0.2)
m.drawcoastlines(linewidth=0.1)
m.drawcountries(linewidth=0.1)
#======= Inset locators = Zoom In =================
"""---->>[lon , lon],[lat, lat]"""
#ax1.plot([300,   360], [50,  20], transform=ccrs.PlateCarree(),linestyle='--',color='orange', linewidth=1,zorder=30)
#ax1.plot([230,   250], [10, -90], transform=ccrs.PlateCarree(),linestyle='--',color='orange', linewidth=1,zorder=30)

m.drawgreatcircle(300,50,360,20,linewidth=1,color='orange',linestyle='--')
m.drawgreatcircle(230,10,250,-90,linewidth=1,color='orange',linestyle='--')
#=============== Draw regional box ================
#fig.subplots_adjust(top=0.88,
#                        bottom=0.11,
#                        left=0.11,
#                        right=0.9,
#                        hspace=0.2,
#                        wspace=0.2)

#plt.savefig('/mnt/g/2nd_Paper/studyArea800dpi.png', dpi=800)
plt.show()

