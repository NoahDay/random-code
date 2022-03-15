#! /usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
print(sys.executable)
print(sys.version)



ds = xr.tutorial.open_dataset(
    "rasm"
)  # use xr.tutorial.load_dataset() for xarray<v0.11.0
ds

print("complete")
