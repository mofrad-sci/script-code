import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors

# Path to the SDF file
input_file = '/home/obiwan/Projects/Stuff/Mofrad/2-Ligands-prep/COCONUT_DB-filtered-3D.sdf'

# Load the SDF file
supplier = Chem.SDMolSupplier(input_file)

# Collect molecular weights
molecular_weights = []

for mol in supplier:
    if mol is not None:
        mol_weight = Descriptors.MolWt(mol)
        molecular_weights.append(mol_weight)

# Plot the distribution of molecular weights
plt.figure(figsize=(10, 6))
sns.histplot(molecular_weights, bins=50, kde=True)
plt.title('Distribution of Molecular Weights')
plt.xlabel('Molecular Weight (Da)')
plt.ylabel('Frequency')
plt.show()

