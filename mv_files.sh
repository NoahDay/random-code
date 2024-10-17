#!/bin/bash
#PBS -N mv-files
#PBS -P ia40
#PBS -q express
#PBS -l wd
#PBS -l mem=90GB
#PBS -m abe
#PBS -M noah.day@adelaide.edu.au
#PBS -l storage=gdata/hh5+gdata/ia40+gdata/gv90

module use /g/data/hh5/public/modules
module load conda/analysis3
module load netcdf
module load gcc
module load openmpi




# Define the input and output folder paths
input_file_path="/g/data/ia40/cice-dirs/runs/waves-025-2019/restart"   # Input path
output_file_path="/g/data/gv90/nd0349/cice-dirs/runs/waves-025-2019/" # Output path



mv "$input_file_path" "$output_file_path"





# # Define the year range for processing
# start_year=2010
# end_year=2020

# max_nj=323
# # Create output folder if it doesn't existy

# mkdir -p "$output_folder_path"

# # Loop through all .nc files in the input folder
# for file in "$input_folder_path"/*.nc; do
#     # Extract filename
#     filename=$(basename "$file")
#     year=$(echo "$filename" | cut -d '.' -f 2 | cut -d '-' -f 1)

#     # Check if the file's year is within the specified range
#     if [[ "$year" -ge "$start_year" && "$year" -le "$end_year" ]]; then
#     	# Define output filename with path
#    	 output_file="$output_folder_path/$filename"

#     	# Subset data for the Southern Hemisphere using ncks
#     	# Assuming 'lat' is the name of your latitude dimension
#     	ncks -d nj,0,$max_nj -O "$file" "$output_file"
   
#     	echo "Processed $filename"
#    fi
# done

# echo "All files processed."
