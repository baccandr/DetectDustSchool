import xarray as xr
import pandas as pd
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.cm import get_cmap
from matplotlib import animation
from matplotlib.axes import Axes

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import cdsapi

import warnings
warnings.simplefilter(action = "ignore", category = RuntimeWarning)

#%% Def fun
def visualize_pcolormesh(data_array, longitude, latitude, projection, color_scale, unit, long_name, vmin, vmax,
                         set_global=True, lonmin=-180, lonmax=180, latmin=-90, latmax=90):
    """
    Visualizes a xarray.DataArray with matplotlib's pcolormesh function.

    Parameters:
        data_array(xarray.DataArray): xarray.DataArray holding the data values
        longitude(xarray.DataArray): xarray.DataArray holding the longitude values
        latitude(xarray.DataArray): xarray.DataArray holding the latitude values
        projection(str): a projection provided by the cartopy library, e.g. ccrs.PlateCarree()
        color_scale(str): string taken from matplotlib's color ramp reference
        unit(str): the unit of the parameter, taken from the NetCDF file if possible
        long_name(str): long name of the parameter, taken from the NetCDF file if possible
        vmin(int): minimum number on visualisation legend
        vmax(int): maximum number on visualisation legend
        set_global(boolean): optional kwarg, default is True
        lonmin,lonmax,latmin,latmax(float): optional kwarg, set geographic extent is set_global kwarg is set to
                                            False

    """
    fig = plt.figure(figsize=(20, 10))

    ax = plt.axes(projection=projection)

    img = plt.pcolormesh(longitude, latitude, data_array,
                         cmap=plt.get_cmap(color_scale), transform=ccrs.PlateCarree(),
                         vmin=vmin,
                         vmax=vmax,
                         shading='auto')

    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=1)

    if (projection == ccrs.PlateCarree()):
        ax.set_extent([lonmin, lonmax, latmin, latmax], projection)
        gl = ax.gridlines(draw_labels=True, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 14}
        gl.ylabel_style = {'size': 14}

    if (set_global):
        ax.set_global()
        ax.gridlines()

    cbar = fig.colorbar(img, ax=ax, orientation='horizontal', fraction=0.04, pad=0.1)
    cbar.set_label(unit, fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    ax.set_title(long_name, fontsize=20, pad=20.0)

    #   plt.show()
    return fig, ax
#%% Load data
file = xr.open_dataset('./Data/20230301_dust_concentration.nc')
# Create timeindex
timestamp = file.time.long_name[19:27]
timestamp_init=datetime.strptime(timestamp,'%Y%m%d' )
time_idx = timestamp_init + pd.to_timedelta(file.time.values, unit='ns')

dust_xr = file.assign_coords(time=time_idx)
longitude = dust_xr.longitude
latitude = dust_xr.latitude
dust_fc = dust_xr.dust

label = 'Dust Concentration'

#%% Export data as csv
dust_df = pd.DataFrame()
#%% Visualize
fig, ax = visualize_pcolormesh(data_array=dust_fc[0,0,:,:],
                               longitude=longitude,
                               latitude=latitude,
                               projection=ccrs.PlateCarree(),
                               color_scale='viridis',
                               unit='-',
                               long_name=label,
                               vmin=0,
                               vmax=100,
                               set_global=False, lonmin=-25, lonmax=45,
                               latmin=30, latmax=72)

frames = 12


def draw(i):
    img = plt.pcolormesh(longitude,
                         latitude,
                         dust_fc[i,0,:,:],
                         cmap='hot_r',
                         transform=ccrs.PlateCarree(),
                         vmin=0,
                         vmax=100,
                         shading='auto',
                         set_global=False, lonmin=-25, lonmax=45,
                         latmin=30, latmax=72)

    ax.set_title(label + ' ' + str(dust_fc.time[i].data)[0:10], fontsize=20, pad=20.0)
    return img


def init():
    return fig


def animate(i):
    return draw(i)


ani = animation.FuncAnimation(fig, animate, frames, interval=800, blit=False,
                              init_func=init, repeat=True)
st.title("Embed Matplotlib animation in Streamlit")
#st.markdown("https://matplotlib.org/gallery/animation/basic_example.html")
components.html(ani.to_jshtml(), height=1000)