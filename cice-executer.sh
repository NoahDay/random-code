#!/bin/bash
#PBS -N cice-nowaves
#PBS -P dy43
#PBS -q normal
#PBS -l ncpus=96
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l wd
#PBS -m abe
#PBS -M noah.day@adelaide.edu.au
#PBS -l storage=gdata/hh5+gdata/dy43+gdata/qv56

module use /g/data/hh5/public/modules
module load conda/analysis3
module load netcdf
module load gcc
module load openmpi


./cice.run
 
