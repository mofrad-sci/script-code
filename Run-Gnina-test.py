import os
import subprocess

# Directory containing the SDF file
input_dir = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/Separate-sdf/1'
# Path to the gnina binary
gnina_binary = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/Separate-sdf/1/gnina'

#Fixed parameters for gnina comman
receptor_file = 'SRK8-prep.pdb'
center_x = -10.28
center_y = -10.94
center_z = 12.21
size_x = 20
size_y = 20
size_z = 20

# Function to run gnina command for each molecule
def run_gnina_for_each_molecule(input_directory):
    for filename in os.listdir(input_directory):
        if filename.endswith(".sdf"):
            coconut_id = os.path.splitext(filename)[0]
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(input_directory, f"{coconut_id}.out.sdf")
            
            command = [
            gnina_binary,
            '-r', receptor_file,
            '-l', input_file,
            '--center_x', srt(center_x),
            '--center_y', srt(center_y),
            '--center_z', srt(center_z),
            '--size_x', str(size_x),
            '--size_y', str(size_y),
            '--size_z', str(size_z),
            '-o', output_file
            ]
            
            print(f"Running command for {coconut_id}: {' '.join(command)}")

# Run the function
run_gnina_for_each_molecule(input_dir)