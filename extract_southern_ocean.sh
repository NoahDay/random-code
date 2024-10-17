#!/bin/bash

# Define the input and output folder paths
input_folder_path="/g/data/gv90/nd0349/sea-ice-classification"   # Input path
output_folder_path="/scratch/ia40/nd0349/published-data" # Output path

# Create output folder if it doesn't exist
mkdir -p "$output_folder_path"

# Loop through all .nc files in the input folder
for file in "$input_folder_path"/*.nc; do
    # Extract filename
    filename=$(basename "$file")

    # Define output filename with path
    output_file="$output_folder_path/$filename"

    # Apply CDO command
    cdo sellonlatbox,-180,180,-90,-40 "$file" "$output_file"

    echo "Processed $filename"
done

echo "All files processed."

