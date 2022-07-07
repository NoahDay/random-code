#!/bin/bash

for file in gridded_ww3.glob_24m.*  #/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hr/*
do
    dir=${file} # include extra folders
    echo 'Opening:' ${dir}
    var=${file##*24m.}
    ym=${var%.*} # year month
    ./interp_cawcr_ncdf_bilinear.py ${ym} icegrid.nc ww3_om2_1deg_${ym}.nc
    #echo 'ww3_om2_1deg_'${ym}'.nc'
    echo 'Formatted: '${dir}
done

