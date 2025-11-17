from rdkit import Chem
import os

# Path to the input SDF file
input_sdf_path = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered-3D-2.sdf'

# Base directory to save the individual SDF files
base_output_dir = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/Separate-sdf'
os.makedirs(base_output_dir, exist_ok=True)

# Number of files per folder
files_per_folder = 10000

# Function to extract and save molecules
def save_molecules_separately(input_sdf, base_output_directory, files_per_folder):
    # Create a supplier object to iterate over molecules in the SDF
    suppl = Chem.SDMolSupplier(input_sdf)
    
    folder_count = 1
    file_count = 0
    
    # Initial output directory
    current_output_dir = os.path.join(base_output_directory, str(folder_count))
    os.makedirs(current_output_dir, exist_ok=True)
    
    # Iterate over each molecule in the SDF
    for mol in suppl:
        if mol is None:
            continue
        
        # Extract the coconut_id (assuming it is stored in a property)
        coconut_id = mol.GetProp("coconut_id")
        
        # Check if we need to create a new folder
        if file_count >= files_per_folder:
            folder_count += 1
            file_count = 0
            current_output_dir = os.path.join(base_output_directory, str(folder_count))
            os.makedirs(current_output_dir, exist_ok=True)
        
        # Define the path for the individual SDF file
        output_file_path = os.path.join(current_output_dir, f"{coconut_id}.sdf")
        
        # Write the molecule to an individual SDF file
        w = Chem.SDWriter(output_file_path)
        w.write(mol)
        w.close()
        
        # Update the file count
        file_count += 1
        print(f"Saved {coconut_id}.sdf in folder {folder_count}")

# Run the function
save_molecules_separately(input_sdf_path, base_output_dir, files_per_folder)

