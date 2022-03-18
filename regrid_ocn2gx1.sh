
# Regrid from auscom grid to gx1

for dir in /Users/a1724548/Gadi/ocn_forcing/grid/auscom/output*
	do
		#echo "${dir}"
		#filename=${dir##*/}$".nc"
		#echo "${filename}"
		#echo "${dir##*t}" # print everytjhing after the final /
		filenumber=${dir##*t} 
		filename=${dir##*/}
		echo "${filename}"
		outfilename=$"ocn_"${filenumber}
		echo ${outfilename}
		cdo remapbil,global_gx1.bathy.nc ${filename} ${outfilename}
done
