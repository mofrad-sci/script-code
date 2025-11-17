import os
from rdkit import Chem
from rdkit.Chem import AllChem

# Define the directory containing your SDF files
sdf_dir = '/home/farzin/Projects/Mofrad/5-VirtualScreening/Separate-sdf/8'  # Update this with the path to your SDF files

# Loop over all files in the directory
for filename in os.listdir(sdf_dir):
    if filename.endswith(".sdf"):
        # Full path to the file
        sdf_path = os.path.join(sdf_dir, filename)
        
        # Read the SDF file using RDKit
        suppl = Chem.SDMolSupplier(sdf_path)
        mol = next(suppl)  # Assuming one molecule per file
        
        if mol is None:
            print(f"Failed to read molecule from {filename}")
            continue
        
        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)
        
        # Try to generate 3D coordinates and optimize the molecule
        try:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
        
        # Write the modified molecule to a temporary SDF file
        temp_sdf_path = os.path.join(sdf_dir, 'temp_' + filename)
        writer = Chem.SDWriter(temp_sdf_path)
        writer.write(mol)
        writer.close()
        
        # Prepare the output PDBQT filename
        pdbqt_path = sdf_path.replace('.sdf', '.pdbqt')
        
        # Call meeko to convert the modified SDF to PDBQT
        try:
            os.system(f'mk_prepare_ligand.py -i {temp_sdf_path} -o {pdbqt_path}')
        except Exception as e:
            print(f"Error converting {temp_sdf_path} to PDBQT: {e}")
        
        # Remove the temporary SDF file
        os.remove(temp_sdf_path)
        
        print(f"Converted {filename} to {pdbqt_path}")

print("Conversion process completed.")

