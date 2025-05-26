import argparse
from rdkit import Chem
import json
from enum import Enum
from pydantic import BaseModel
from typing import List
from utils import load_smiles_from_csv

# Define atom type enum
class AtomTypeEnum(str, Enum):
    H = "H"
    C = "C"
    N = "N"
    O = "O"
    F = "F"

# Define Bond model
class Bond(BaseModel):
    atom: 'Atom'
    bond_type: str

# Define Atom model
class Atom(BaseModel):
    atom_id: int
    atom_type: AtomTypeEnum
    bonds: List[Bond]

# Needed to resolve forward references in Pydantic
Bond.update_forward_refs()

def smiles_to_tree(atom_format, bond_format, atom_type_enum, smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    atom_idx_to_atom = {}
    visited_atoms = set()
    queue = []
    added_bonds = set([])

    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:
            atom_molecule = atom_format(atom_id=atom.GetIdx(),
                                        atom_type=atom_type_enum(atom.GetSymbol()),
                                        bonds=[])
            atom_idx_to_atom[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx(), None))
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx, bond_type = queue.pop(0)
        current_atom = atom_idx_to_atom[current_idx]
        if parent_idx is not None:
            if current_idx in visited_atoms:
                current_atom = atom_format(atom_id=current_idx,
                                           atom_type=atom_type_enum(mol.GetAtomWithIdx(current_idx).GetSymbol()),
                                           bonds=[])
            parent_atom = atom_idx_to_atom[parent_idx]
            bond = bond_format(atom=current_atom, bond_type=bond_type)
            parent_atom.bonds.append(bond)
            visited_atoms.add(current_idx)
            visited_atoms.add(parent_idx)

        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if (current_idx, neighbor_idx) not in added_bonds and (neighbor_idx, current_idx) not in added_bonds:
                neighbor_bond_type = str(bond.GetBondType()).replace("BondType.", "")
                if neighbor_bond_type == "DATIVE":
                    neighbor_bond_type = "SINGLE"
                if neighbor_idx not in atom_idx_to_atom:
                    neighbor_atom = atom_format(atom_id=neighbor_idx,
                                                atom_type=atom_type_enum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()),
                                                bonds=[])
                    atom_idx_to_atom[neighbor_idx] = neighbor_atom
                queue.append((current_idx, neighbor_idx, neighbor_bond_type))
                added_bonds.add((current_idx, neighbor_idx))

    return atom_idx_to_atom[0]

# Custom JSON serializer for Enum
def custom_json_serializer(obj):
    if isinstance(obj, Enum):
        return obj.value
    raise TypeError("Type not serializable")

# Example usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', default='qm9.csv', help='Path to the SMILES data CSV file')
    # parser.add_argument('--output', required=True, help='Path to save the motif vocabulary JSON file')
    parser.add_argument('--freq_threshold', type=int, default=10, help='Frequency threshold for motif selection (default: 10)')
    args = parser.parse_args()

    file_path = 'dataset/' + args.data
    smiles_list = load_smiles_from_csv(file_path)

    mol_tree_list = []
    for smiles in smiles_list:
        mol_tree = smiles_to_tree(Atom, Bond, AtomTypeEnum, smiles)
        mol_tree_list.append(mol_tree)

    molecule_json = json.dumps([mol.dict() for mol in mol_tree_list], indent=2, default=custom_json_serializer)

    with open("./outputs/mol2json.json", "w") as f:
        f.write(molecule_json)
        print('JSON tree saved successfully.')
