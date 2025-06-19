import numpy as np
import pandas as pd
from torch_geometric.datasets import QM9

def split_and_save_smiles_qm9():
    # Load QM9 dataset
    dataset = QM9(root='data/QM9')
    
    # Ensure reproducibility
    np.random.seed(42)
    
    # Shuffle the dataset indices
    indices = np.random.permutation(len(dataset))
    
    # Select 10,000 samples
    selected_indices = indices[:10000]
    
    # Split into training (8k), validation (1k), and testing (1k)
    train_indices = selected_indices[:8000]
    val_indices = selected_indices[8000:9000]
    test_indices = selected_indices[9000:]
    
    # Helper function to extract SMILES strings
    def extract_smiles(indices):
        smiles_list = []
        for idx in indices:
            smiles = dataset[idx].smiles
            smiles_list.append(smiles)
        return pd.DataFrame({"SMILES": smiles_list})
    
    # Extract SMILES for each subset
    train_smiles = extract_smiles(train_indices)
    val_smiles = extract_smiles(val_indices)
    test_smiles = extract_smiles(test_indices)
    
    # Save SMILES to CSV
    train_smiles.to_csv("dataset/qm9_train_smiles.csv", index=False)
    val_smiles.to_csv("dataset/qm9_val_smiles.csv", index=False)
    test_smiles.to_csv("dataset/qm9_test_smiles.csv", index=False)
    
    print("SMILES strings saved to: qm9_train_smiles.csv, qm9_val_smiles.csv, qm9_test_smiles.csv")

# Run the function
split_and_save_smiles_qm9()
