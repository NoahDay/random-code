#!/bin/bash
# Extracting the SST from the 3D ocean temperature data from
# /g/data/ik11/outputs/access-om2-01/01deg_jra55v140_iaf_cycle3/output731/ocean
# Method:
#	1. Copy the ocean-3d-temp-1-daily-mean-ym_20XX_XX from ik11 to scratch
#	2. Extract layer 1 of u and v
#	3. Remove the copied file


# Testing for looping over directories
#i=0
for dir in /g/data/ik11/outputs/access-om2/1deg_jra55_iaf_omip2_cycle6/output*/
do
	echo "Copying over to scratch.."
	dir=${dir%*/}      # remove the trailing "/"
	#echo "${dir##*/}" # print everytjhing after the final /
	#echo "${dir}"
	ocndir="/ocean/ocean_month.nc"
	# Get the directories for files starting with ocean-3d-temp-1-daily-mean-ym_20
	for file in ${dir}$"ocean/ocean-3d-temp-1-daily-mean-ym_20*"
	do
		#filename=${dir##*/}$".nc"
		echo ${file}
		#filedir=${dir}$ocndir
	done
	# 1.
	#cp filedir /scratch/df0/ndo0349
	# 2.
	#cdo selvar,sss,sst,u,v ${filedir} temp.nc
	# 3.
	#ncrcat -d st_ocean,1 temp.nc ${filename}
	#echo "${filename} formatted!"
	# 4.
	#rm temp.nc
	#((++i < 5)) || break # Try for 5 files
done

# 1.
#echo "Copying data over to scratch.."
#cp /g/data/ik11/outputs/access-om2/1deg_jra55_iaf_omip2_cycle6/output365/ocean/ocean_month.nc /scratch/df0/nd0349
#echo "Copying complete."

# 2.
#echo "Extracting sss,sst,u,v"
#cdo selvar,sss,sst,u,v ocean_month.nc temp.nc

# 3.
#echo "Taking over the surface layer"
#ncrcat -d st_ocean,1 temp.nc ${filename}
#echo "Formatting complete!"

# 4.
#rm temp.nc
