import rdkit.Chem as Chem
from utils import *
import json
import argparse
from utils import load_smiles_from_csv

class MolTreeNode(object):

    def __init__(self, smiles, clique=[]):
        self.smiles = smiles
        self.mol = get_mol(self.smiles)

        self.clique = [x for x in clique] #copy
        self.neighbors = []

    def add_neighbor(self, nei_node):
        self.neighbors.append(nei_node)

class MolTree(object):

    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = get_mol(smiles)

        cliques, edges = tree_decomp(self.mol)

        # print(cliques)
        # print(edges)

        self.nodes = []
        root = 0
        for i,c in enumerate(cliques):
            cmol = get_clique_mol(self.mol, c)
            node = MolTreeNode(get_smiles(cmol), c)
            self.nodes.append(node)
            if min(c) == 0: root = i

        for x,y in edges:
            self.nodes[x].add_neighbor(self.nodes[y])
            self.nodes[y].add_neighbor(self.nodes[x])
        
        if root > 0:
            self.nodes[0],self.nodes[root] = self.nodes[root],self.nodes[0]

        for i,node in enumerate(self.nodes):
            node.nid = i + 1
            if len(node.neighbors) > 1: #Leaf node mol is not marked
                set_atommap(node.mol, node.nid)
            node.is_leaf = (len(node.neighbors) == 1)

def dfs(node, fa_idx = None):
    smiles = node.smiles
    childs = []
    for child in node.neighbors:
        if child.nid == fa_idx: continue
        childs.append(dfs(child, node.nid))
        # max_depth = max(max_depth, dfs(child, node.idx))

    return {'node': smiles, 'childs': childs}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', default='qm9.csv', help='Path to the SMILES data CSV file')
    # parser.add_argument('--output', required=True, help='Path to save the motif vocabulary JSON file')
    parser.add_argument('--freq_threshold', type=int, default=10, help='Frequency threshold for motif selection (default: 10)')
    args = parser.parse_args()

    file_path = 'dataset/' + args.data

    smiles_list = load_smiles_from_csv(file_path)
    jt_list = []

    for smiles in smiles_list:
        if smiles == "C1CC2=CC=CN12":
            continue
        cset = set()
        print('smiles:', smiles)
        mol = MolTree(smiles)
        # for c in mol.nodes:
        #     cset.add(c.smiles)
        #     print(c.smiles)

        jt = dfs(mol.nodes[0])
        jt_list.append(jt)

        # for x in cset:
        #     print(x)

    json_jts = json.dumps([jt for jt in jt_list])
    # print(json_string)

    with open("./outputs/mol2JT.json", "w") as f:
        f.write(json_jts)
        print('JSON tree saved successfully.')
