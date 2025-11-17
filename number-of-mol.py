from rdkit import Chem

# Path to the SDF file
input_file = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered-3D-2.sdf'

# Load the SDF file
supplier = Chem.SDMolSupplier(input_file)

# Initialize a counter
molecule_count = 0

# Count the molecules
for mol in supplier:
    if mol is not None:
        molecule_count += 1

print("Number of molecules in the SDF file:", molecule_count)

