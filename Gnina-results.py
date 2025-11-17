import os
import csv

# Paths to your directories
dir1 = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf/Top30000/1"
dir2 = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf/Top30000/2"

# Function to extract the highest CNNaffinity value from a given SDF file
def extract_max_cnn_affinity(sdf_file):
    affinities = []
    with open(sdf_file, 'r') as file:
        content = file.readlines()
        # Traverse through the file lines
        for i, line in enumerate(content):
            # Identify the line containing CNNaffinity
            if "> <CNNaffinity>" in line:
                affinity_value = float(content[i + 1].strip())
                affinities.append(affinity_value)
    # Return the maximum affinity found
    return max(affinities) if affinities else None

# Function to process a directory and gather the max CNNaffinity for all SDF files
def process_directory(directory):
    affinity_data = []
    
    for file_name in os.listdir(directory):
        if file_name.endswith("-out.sdf"):
            coconut_id = file_name.split("-out.sdf")[0]
            sdf_file_path = os.path.join(directory, file_name)
            max_affinity = extract_max_cnn_affinity(sdf_file_path)
            if max_affinity is not None:
                affinity_data.append([coconut_id, max_affinity])
    
    return affinity_data

# Process both directories
affinities_dir1 = process_directory(dir1)
affinities_dir2 = process_directory(dir2)

# Combine results from both directories
all_affinities = affinities_dir1 + affinities_dir2

# Write the results to a CSV file
output_csv = "max_cnn_affinities.csv"
with open(output_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Coconut ID", "Max CNNaffinity"])
    writer.writerows(all_affinities)

print(f"Results saved to {output_csv}")

