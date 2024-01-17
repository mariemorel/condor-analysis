#!/usr/bin/env python

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from iteround import saferound
import numpy as np

import argparse

parser = argparse.ArgumentParser(description="pastmlfile into fasta")
parser.add_argument("pastml_acr", help="pastml output")
parser.add_argument("marginal", help = "marginal probabilities of root")
parser.add_argument("output", help="name of the output")
args = parser.parse_args()

pastml_acr = str(args.pastml_acr)
output = str(args.output)
nb_simu = int(args.nb_simu)

characters = []
marginal_list = []
with open(args.marginal) as f:
    for line in f:
        if line.strip().isdigit(): 
            characters.append(line.strip())
        else:
            marginal_list.append(line.strip().split("\t")[1:])

AA = list("ARNDCQEGHILKMFPSTWYV")

root_dict = {}
for a,pos in enumerate(characters): 
    root= []
    probs = []
    for proba in marginal_list[a]:
        probs.append(nb_simu*float(proba))
    for i, count in enumerate(saferound(probs,0)):
        root.extend( int(count) * [AA[i]])
    root_dict[pos] = root

##dictionary of each pos with a vector of amino acids corresponding at the root amino acids mutliplied by the marginal proba * nb simulations

#A = pd.read_csv(pastml_acr, sep='\t', index_col = 0)
A = pd.read_csv(pastml_acr, sep='\t', index_col = 0, compression = "gzip")
B = A[[str(i-1) for i in Pos_list]]

#create fasta sequences of acr output.
record = []
for n in set(B.index):
    if n == "root" || n == "N000000001":
    record.append(SeqRecord(Seq("".join(list(B.loc[n]))), id = n, name = n, description = n ))
SeqIO.write(record, output+"pastml_acr.fasta", "fasta")
