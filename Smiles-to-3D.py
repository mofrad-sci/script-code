from rdkit import Chem
from rdkit.Chem import AllChem

def has_valence_issues(mol):
    try:
        Chem.SanitizeMol(mol)
        return False
    except Chem.SanitizeException:
        return True

def convert_2d_to_3d(input_sdf, output_sdf):
    # Read the 2D SDF file
    suppl = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    
    for mol in suppl:
        if mol is None or has_valence_issues(mol):
            continue

        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == 0:
            try:
                # Optimize the 3D structure
                AllChem.UFFOptimizeMolecule(mol)

                # Write the 3D molecule to the output file
                writer.write(mol)
            except Exception as e:
                print(f"Optimization failed for molecule: {Chem.MolToSmiles(mol)}\nError: {e}")

    writer.close()


# Usage
input_sdf = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered.sdf'  # Replace with your input 2D SDF file
output_sdf = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered-3D-2.sdf'  # Replace with your desired output 3D SDF file
convert_2d_to_3d(input_sdf, output_sdf)

