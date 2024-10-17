#!/bin/bash
#PBS -N rsync-restart
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
input_folder_path="/g/data/ia40/cice-dirs/runs/nowaves-10/restart"   # Input path
output_folder_path="/g/data/gv90/nd0349/cice-dirs/runs/nowaves-10/" # Output path

rsync -avzh --progress "$input_folder_path" "$output_folder_path"
echo "All files processed."
