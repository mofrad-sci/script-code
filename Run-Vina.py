import os
import subprocess

# Set your directories and paths
pdbqt_dir = "/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/3"
receptor = "/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/SRK8-200ns.converted.pdbqt"
vina = "/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/vina"  # Updated to point directly to the Vina executable

# Check if the Vina binary is executable
if not os.path.isfile(vina):
    raise FileNotFoundError(f"The Vina binary at {vina} does not exist.")
if not os.access(vina, os.X_OK):
    raise PermissionError(f"The Vina binary at {vina} is not executable.")

# Loop over each PDBQT file in the directory
for pdbqt_file in os.listdir(pdbqt_dir):
    if pdbqt_file.endswith(".pdbqt"):
        # Extract the coconut_id (file name without extension)
        coconut_id = os.path.splitext(pdbqt_file)[0]

        # Paths for input and output files
        pdbqt_path = os.path.join(pdbqt_dir, pdbqt_file)
        out_path = os.path.join(pdbqt_dir, f"{coconut_id}-out.pdbqt")

        # Ensure the input file is a file and not a directory
        if not os.path.isfile(pdbqt_path):
            raise ValueError(f"Expected a file but got a directory: {pdbqt_path}")

        # Run AutoDock Vina with error handling
        vina_command = [
            vina,
            "--receptor", receptor,
            "--ligand", pdbqt_path,
            "--center_x", "46.27",
            "--center_y", "19.67",
            "--center_z", "78.06",
            "--cpu", "8",
            "--exhaustiveness", "8",
            "--size_x", "24",
            "--size_y", "24",
            "--size_z", "24",
            "--out", out_path
        ]

        try:
            subprocess.run(vina_command, check=True)
            print(f"Docking completed for {coconut_id}")
        except subprocess.CalledProcessError as e:
            print(f"Error during docking for {coconut_id}: {e}")
        except Exception as e:
            print(f"Unexpected error for {coconut_id}: {e}")

print("All docking processes are completed.")

