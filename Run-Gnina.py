import os
import subprocess

# Directory containing the SDF files
input_dir = '/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/Top30000'

# Path to the gnina binary
gnina_binary = '/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/gnina'

# Fixed parameters for gnina command
receptor_file = 'SRK8-200ns.pdb'
center_x = 46.27    
center_y = 19.67
center_z = 78.06
size_x = 24
size_y = 24
size_z = 24

# Function to run gnina command for each molecule
def run_gnina_for_each_molecule(input_directory):
    for filename in os.listdir(input_directory):
        if filename.endswith(".sdf"):
            coconut_id = os.path.splitext(filename)[0]
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(input_directory, f"{coconut_id}-out.sdf")
            
            command = [
                gnina_binary,
                '-r', receptor_file,
                '-l', input_file,
                '--center_x', str(center_x),
                '--center_y', str(center_y),
                '--center_z', str(center_z),
                '--size_x', str(size_x),
                '--size_y', str(size_y),
                '--size_z', str(size_z),
                '-o', output_file
            ]
            
            print(f"Running command for {coconut_id}: {' '.join(command)}")
            subprocess.run(command)

# Run the function
run_gnina_for_each_molecule(input_dir)

