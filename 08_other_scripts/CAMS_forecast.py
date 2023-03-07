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
#%% Download data from CAMS
url = 'https://ads.atmosphere.copernicus.eu/api/v2'
key = '12714:ce643f0d-6679-4f36-8c32-a8c225ffe63a'

c = cdsapi.Client(url=url, key=key)

c.retrieve(
    'cams-europe-air-quality-forecasts',
    {
        'model': 'ensemble',
        'date': '2023-03-01/2023-03-01',
        'format': 'netcdf',
        'variable': 'dust',
        'type': 'forecast',
        'time': '00:00',
        'leadtime_hour': [
            '0', '12', '15',
            '18', '21', '24',
            '27', '3', '30',
            '33', '36', '39',
            '42', '45', '48',
            '51', '54', '57',
            '6', '60', '63',
            '66', '69', '72',
            '75', '78', '81',
            '84', '87', '9',
            '90',
        ],
        'level': '0',
    },
    './Data/20230301_dust_concentration.nc')

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
#%%
for j,k in enumerate(time_idx):
    fig, ax = visualize_pcolormesh(data_array=dust_fc[j,0,:,:],
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
    plt.savefig('./Figures/20230301_dust_concentration_'+str(j)+'.png', dpi = 150)

#%%
