# from torch_geometric.datasets import QM9
from rdkit import Chem
from rdkit.Chem import Draw

# # Download and transform
# dataset = QM9(root='data/QM9')

# # One molecule
# data = dataset[0]

# # What does it have?
# print(data)
# print(data.smiles)

mol = Chem.MolFromSmiles('C1=CC1')  # Example: ethanol
img = Draw.MolToImage(mol)
img.save("molecule.png")
