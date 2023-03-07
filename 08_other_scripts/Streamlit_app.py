import xarray as xr
#import pandas as pd
from datetime import datetime
import pydeck as pdk

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

import streamlit as st
import streamlit.components.v1 as components

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

#%% Visualize
dust_fc_0 = dust_fc[0,0,:,:]
df_dust = dust_fc_0.to_dataframe().reset_index()[['latitude', 'longitude', 'dust']]

COLOR_BREWER_BLUE_SCALE = [
    [240, 249, 232],
    [204, 235, 197],
    [168, 221, 181],
    [123, 204, 196],
    [67, 162, 202],
    [8, 104, 172],
]

view = pdk.data_utils.compute_view(df_dust [["longitude", "latitude"]])
view.zoom = 6

cattle = pdk.Layer(
    "HeatmapLayer",
    data=df_dust,
    opacity=0.9,
    get_position=["longitude", "latitude"],
    aggregation=pdk.types.String("MEAN"),
    color_range=COLOR_BREWER_BLUE_SCALE,
    threshold=1,
    get_weight="weight",
    pickable=True,
)

r = pdk.Deck(
    layers=[cattle],
    initial_view_state=view,
    map_provider="mapbox",
    tooltip={"text": "Concentration of cattle in blue, concentration of poultry in orange"},
)
r.to_html("heatmap_layer.html")