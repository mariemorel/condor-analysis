#!/usr/bin/env python

import pandas as pd
from ete3 import PhyloTree
import os

# handle arguments
import argparse


parser = argparse.ArgumentParser(description="convergent scenario to use in pcoc")
parser.add_argument("tree", help="path towards tree", type=str)
parser.add_argument("pheno", help="path towards phenotype", type=str)
parser.add_argument("cvgt", help="column name with phenotype True False", type=str)

args = parser.parse_args()
tree_path = args.tree
pheno_path = args.pheno
target = args.cvgt

path_name = os.path.dirname(os.path.abspath(pheno_path))

t = PhyloTree(tree_path, format=1)

n = 0  # numbering of nodes from 0
for node in t.traverse("postorder"):
    node.add_features(etiquette=n)
    n += 1

A = pd.read_csv(pheno_path, sep="\t")

map_pcoc_leaf = {}
for node in t.traverse("postorder"):
    if node.is_leaf():
        map_pcoc_leaf[node.name] = node.etiquette


pcoc_nb = []
for i in A.id:
    if i in map_pcoc_leaf.keys():
        pcoc_nb.append(map_pcoc_leaf[i])
    else:
        pcoc_nb.append("na")

A["pcoc_nb"] = pcoc_nb

A.to_csv(
    path_name + "/species_cvgt_pcoc.tsv", sep="\t", index=False,
)

pcoc_t = t.copy()

for node in pcoc_t.traverse():
    node.name = node.etiquette

pcoc_t.write(
    format=1, outfile=path_name + "/pcoc_tree_rerooted.nw",
)

cvgts_spp = A[A[target]].id  # is marine = False

for node in t.traverse("postorder"):
    if node.is_leaf():
        if node.name in list(cvgts_spp):
            node.add_features(cvgt=1)  # leaves with fresh or brackish = 1
        else:
            node.add_features(cvgt=0)  # leaves with only marine = 0

nodes_list = []

for node in t.get_monophyletic(values=[1], target_attr="cvgt"):
    nodes_list.append([i.etiquette for i in node.traverse("preorder")])


scenario_no_acr = ""
for i in nodes_list:
    if len(i) > 1:
        scenario_no_acr += ",".join([str(a) for a in i]) + "/"
    else:
        scenario_no_acr += str(i[0]) + "/"

with open(path_name + "/scenario_pcoc_noacr.txt", "w",) as wf:
    wf.write(scenario_no_acr[:-1])

