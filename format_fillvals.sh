#!/bin/bash
# format_fillvals.sh will set the Fill Values of all netcdf files in the specified directory to 0, 
# such that the formatting is appropriate for interp_jra55_ncdf_bilinear.py

# Set to appropriate directory
for dir in /Users/a1724548/GitHub/interp/*12312230.nc
do
	echo "Copying over to scratch.."
	dir=${dir%*/}      # remove the trailing "/"
	#echo "${dir##*/}" # print everytjhing after the final /
	#echo "${dir}"
	ocndir="/ocean/ocean_month.nc"
	filename=${dir##*/}$".nc"
	echo ${dir}$ocndir
	filedir=${dir}$ocndir
	# 1.
	#cp filedir /scratch/df0/ndo0349
	# 2.
	#cdo selvar,sss,sst,u,v ${filedir} temp.nc
	# 3.
	#ncrcat -d st_ocean,1 temp.nc ${filename}
	echo "${filename} formatted!"
	# 4.
	#rm temp.nc
	#((++i < 5)) || break # Try for 5 files


    ncatted -a _FillValue,rsds,o,f,0.0 rsds*
done