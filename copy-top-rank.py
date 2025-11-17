import os
import shutil

# Path to the root directory containing the 34 subdirectories with .sdf files
sdf_source_dir = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf"

# Path to the text file containing the first 30000 coconut IDs
ligand_list_file = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf/first_30000_ligands.txt"

# Path to the directory where the .sdf files will be copied
destination_dir = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf/Top30000"

# Create the destination directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# Read the first 30000 coconut IDs from the text file
with open(ligand_list_file, 'r') as file:
    coconut_ids = [line.strip() for line in file]

# Traverse all subdirectories and files in the source directory
for root, dirs, files in os.walk(sdf_source_dir):
    for file_name in files:
        if file_name.endswith(".sdf"):
            coconut_id = file_name.rsplit('.', 1)[0]  # Extract the coconut ID from the file name
            
            # Check if this file matches one of the first 30000 coconut IDs
            if coconut_id in coconut_ids:
                # Full path to the .sdf file
                sdf_file_path = os.path.join(root, file_name)

                # Destination file path
                destination_file_path = os.path.join(destination_dir, file_name)

                try:
                    # Copy the file to the destination directory
                    shutil.copy(sdf_file_path, destination_file_path)
                    print(f"Copied {file_name} to {destination_dir}")
                except shutil.SameFileError:
                    print(f"Skipped {file_name}: Source and destination are the same file.")
                except Exception as e:
                    print(f"Error copying {file_name}: {e}")

print("File copying complete.")
