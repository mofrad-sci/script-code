import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

# Define your criteria
min_molecular_weight = 200  # Minimum threshold for molecular weight
max_molecular_weight = 700  # Maximum threshold for molecular weight
metal_symbols = ["Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", 
                 "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", 
                 "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", 
                 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", 
                 "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", 
                 "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", 
                 "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", 
                 "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", 
                 "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

# Load the SDF file
input_file = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB.sdf'
output_file = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered.sdf'

# Function to check if a molecule contains any metal atoms
def contains_metal(molecule):
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() in metal_symbols:
            return True
    return False

# Read molecules from SDF, filter them, and write to a new SDF file
supplier = Chem.SDMolSupplier(input_file)
writer = Chem.SDWriter(output_file)

for mol in supplier:
    if mol is not None:
        # Filter based on molecular weight and metal content
        mol_weight = Descriptors.MolWt(mol)
        if min_molecular_weight <= mol_weight <= max_molecular_weight and not contains_metal(mol):
            writer.write(mol)

# Close the writer
writer.close()

print("Filtering complete. Filtered molecules saved to", output_file)

