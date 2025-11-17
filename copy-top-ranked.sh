#!/bin/bash

# Directories where to search for the files
dir1="1"
dir2="2"
output_dir="new_folder"

# Create the new folder if it doesn't exist
mkdir -p "$output_dir"

# List of Coconut IDs to search for (part of file names)
coconut_ids=(
    "CNP0096230" "CNP0358468" "CNP0438599" "CNP0401153" "CNP0463625"
    "CNP0324901" "CNP0252937" "CNP0099591" "CNP0448122" "CNP0455553"
    "CNP0125896" "CNP0007042" "CNP0417838" "CNP0469177" "CNP0367374"
    "CNP0151909" "CNP0439915" "CNP0416936" "CNP0095710" "CNP0473015"
    "CNP0274291" "CNP0388012" "CNP0464401" "CNP0090611" "CNP0455763"
    "CNP0088487" "CNP0091189" "CNP0096585" "CNP0439879" "CNP0419444"
)

# Function to search and copy files containing the Coconut ID in their name
copy_files() {
    local coconut_id=$1
    # Search in both directories for files containing the Coconut ID in the name
    file_path=$(find "$dir1" "$dir2" -type f -name "*$coconut_id*" 2>/dev/null)
    
    if [ -n "$file_path" ]; then
        for file in $file_path; do
            cp "$file" "$output_dir"
            echo "Copied $file to $output_dir"
        done
    else
        echo "File containing Coconut ID $coconut_id not found"
    fi
}

# Iterate through the list and process each Coconut ID
for id in "${coconut_ids[@]}"; do
    copy_files "$id"
done

