import rdkit
import rdkit.Chem as Chem
from utils import *
import json

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
    import sys
    lg = rdkit.RDLogger.logger() 
    lg.setLevel(rdkit.RDLogger.CRITICAL)

    cset = set()
    for line in sys.stdin:
        smiles = line.split()[0]
        # print('smiles:', smiles)
        mol = MolTree(smiles)
        # for c in mol.nodes:
        #     cset.add(c.smiles)
        #     print(c.smiles)

        jt = dfs(mol.nodes[0])
        json_string = json.dumps(jt)
        print(json_string)
    # for x in cset:
    #     print(x)
