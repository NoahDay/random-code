#!/bin/bash
# Method:
#	1. Copy the JRA55 data from respective folder to interp
#	2. Rename so the dates all match
#	3. Run interp_jra55_ncdf_bilinear
#	4. Remove the copied file
echo What year would you like to format?
read year
declare -a fileJRA=("/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hr/prsn" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hr/rsds" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hr/rlds" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hrPt/tas" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hrPt/uas" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hrPt/vas" "/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hrPt/huss")
filenameWork=()
for dir1 in ${fileJRA[@]} #/Volumes/NoahDay5TB/Gadi/JRA55-do-1-4-0/atmos/3hr/*
do
    dir=${dir1}/gr/v20190429/ # include extra folders
    var=${dir1##*/}
    filename=()
    # Store the filename
    while IFS=  read -r -d $'\0'; do
        filename+=("$REPLY")
    done < <(find ${dir} -name "*${year}*" -print0)
    # 1.
    echo " "
    echo "Copying over ${filename##*/} to interp work directory"
    cp ${filename} /Volumes/NoahDay5TB/Gadi/interp
    # 2.
    echo "Renaming... ${var}"
    mv ${filename##*/} ${var}${year}.nc # Rename file in interp
    filenameWork+=("${var}${year}.nc")
done

# 3.
./interp_jra55_ncdf_bilinear.py ${year} icegrid.nc JRA55_om2_1deg_03hr_forcing_${year}.nc
# 4.
#echo Are you ready to delete the files?
#read decision
#if ["decision"="yes"]; then
for dir1 in ${filenameWork[@]}
do
    echo "Removing ${dir1}"
    rm ${dir1}
done
#fi



# =============================================================================================
# ========================================= EXTRA CODE ========================================
# =============================================================================================

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
    # dir=${dir%*/}      # remove the trailing "/"
   
    #echo "${dir##*/}" # print everytjhing after the final /
    #echo "${dir}"
