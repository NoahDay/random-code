#! /usr/bin/env python3

import sys as s
import numpy as np
import xarray # used for reading the data.
import matplotlib.pyplot as plt # used to plot the data.
import ipywidgets as widgets # For ease in selecting variables.
import cartopy.crs as ccrs # Used to georeference data.

filelist_arr = [save_dir + os.path.basename(file) for file in filelist]
selected_file = widgets.Dropdown(options=filelist_arr, description='data file')
display(selected_file)

# Now to load in the data to xarray
ds = xarray.open_dataset(selected_file.value)

# Helper methods# Define function to get standard dimensions
def get_time(dataset):
    for _,cur_coord in dataset.coords.items:
        if cur_coord.attrs['standard_name'] == 'time':
            return cur_coord
def get_lat(dataset):
    for _,cur_coord in dataset.coords.items:
        if cur_coord.attrs['standard_name'] == 'longitude':
            return cur_coord
def get_lon(dataset):
    for _,cur_coord in dataset.coords.items:
        if cur_coord.attrs['standard_name'] == 'latitude':
            return cur_coord

def get_primary(dataset):
    primary_variables = {}
    coords = dataset.coords.keys()
    highest_dims = 0
    for cur_key,cur_var in dataset.variables.items():
        if cur_key not in coords:
            primary_variables[cur_key] = cur_var
    return primary_variables
    
    
var = widgets.Dropdown(
    options=get_primary(ds).keys(),
    description='Variable')
display(var)

proj = ccrs.Mercator()
plt.gcf().set_size_inches(20,10)
ax = plt.axes(projection=proj)
data_slice = ds[var.value].isel(time=0)
data_slice.plot.contourf(ax=ax, transform=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()

print("saas")
