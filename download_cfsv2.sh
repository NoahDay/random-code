#!/bin/bash

url="https://thredds.rda.ucar.edu/thredds/fileServer/files/g/d094001/2015/wnd10m.cdas1.201508.grb2"

# Output filename
output="wnd10m.cdas1.201508.grb2.nc"

# Download using wget
wget -O "$output" "$url"
