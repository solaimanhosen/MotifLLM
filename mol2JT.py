from rdkit import Chem
from rdkit.Chem import rdmolops
import networkx as nx

def smiles_to_junction_tree(smiles):
    # Parse SMILES and sanitize molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    # Find all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = [set(ring_atoms) for ring_atoms in ring_info.AtomRings()]
    
    # Define SMARTS patterns for common functional groups
    functional_groups = {
        'hydroxyl': '[OH]',
        'carbonyl': 'C=O',
        'carboxylic_acid': 'C(=O)O',
        'amine': '[NH2]',
    }
    
    # Detect functional groups
    clusters = rings.copy()
    for group_name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            clusters.append(set(match))
    
    # Remove duplicate clusters
    unique_clusters = []
    seen = set()
    for cluster in clusters:
        frozen_cluster = frozenset(cluster)
        if frozen_cluster not in seen:
            seen.add(frozen_cluster)
            unique_clusters.append(cluster)
    clusters = unique_clusters
    
    # Build cluster connection graph
    cluster_graph = nx.Graph()
    for i, cluster in enumerate(clusters):
        cluster_graph.add_node(i, atoms=tuple(cluster))
    
    # Add edges between overlapping clusters
    for i in range(len(clusters)):
        for j in range(i+1, len(clusters)):
            intersection = clusters[i].intersection(clusters[j])
            if intersection:
                cluster_graph.add_edge(i, j, weight=len(intersection))
    
    # Convert to junction tree using maximum spanning tree
    junction_tree = nx.maximum_spanning_tree(cluster_graph)
    
    return junction_tree, clusters

# Example usage
smiles = 'C1CCCCC1O'  # Cyclohexanol
jt, clusters = smiles_to_junction_tree(smiles)

print("Junction Tree Edges:")
print(jt.edges(data=True))
print("\nCluster Details:")
for idx, cluster in enumerate(clusters):
    print(f"Node {idx} contains atoms: {cluster}")