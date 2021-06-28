#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
from scipy import linalg

from ete3 import PhyloTree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
    description="Create an amino acid alignment corresponding to a tree and an ancestral sequence")
parser.add_argument("root", help="root sequence")
parser.add_argument("tree_file", help="rooted tree ")
parser.add_argument("rates", help="file with the rates per position")
parser.add_argument("freq", help="freqs amino acid of the model")
parser.add_argument("mode", help="test mode that can be JTT or HIVb ")
parser.add_argument("output", help="name of the output")
args = parser.parse_args()

tree = PhyloTree(args.tree_file, format=1)
root = list(args.root)
len_seq = len(root)
output = str(args.output)


print(len(root))

rate = []
with open(args.rates) as f:
    for line in f:
        rate.append(float(line.strip()))


frequencies = []
with open(args.freq) as f:
    for line in f:
        frequencies.append(float(line.strip()))


AA = list("ARNDCQEGHILKMFPSTWYV")


if args.mode == "HIVb":
    S_df = pd.read_csv(
        '/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/args/raw/HIVb_phyml.model', sep='\t', index_col=0)
    PI = np.array(frequencies)

else:
    S_df = pd.read_csv(
        '/pasteur/zeus/projets/p01/Evolbioinfo/users/mamorel/Projet_Convergence/Data/args/raw/JTT_phyml.model', sep='\t', index_col=0)
    PI = np.array(frequencies)


S = np.array(S_df)


# generator or instantaneous rate matrix
Q = S*PI/100
for i in range(len(Q)):  # sum line must be 0
    Q[i, i] = -(sum(Q[i]))

Q_df_zero = pd.DataFrame(Q, columns=AA)
Q_df_zero.index = AA


# now we want to normalize Q such as ∑(i) pi  ∑(j!=i) Qij  = 1
# norm(Q) = (1/µ Q), µ = - ∑(i)пi*qii

mu = 0
for i in range(20):
    mu += PI[i] * Q_df_zero.iloc[i][i]

# print(mu)

# Normalization
Q_norm = -1/mu * Q


class Simulation:
    def __init__(self):
        self.AA = list("ARNDCQEGHILKMFPSTWYV")
        self.uniq_rate = list(set([round(i, 3) for i in rate]))

    def _random_choice_prob_index(self, a, axis=1):
        r = np.expand_dims(np.random.rand(a.shape[1-axis]), axis=axis)
        return (a.cumsum(axis=axis) > r).argmax(axis=axis)

    def Simulator_choice(self, T, root):
        Alignment = []
        Whole_align = []

        def process_node(node, parent_seq):
            P = np.array([linalg.expm(Q_norm*node.dist*i)
                          for i in self.uniq_rate])  # ?,20,10
            target_seq = self._random_choice_prob_index(
                P[rate_index, parent_seq, :], axis=1)  # create a list for the simulations
            Whole_align.append(SeqRecord(
                Seq("".join([self.AA[i] for i in target_seq])), node.name, description=""))
            if node.is_leaf():
                #print("".join([self.AA[i] for i in target_seq]))
                Alignment.append(SeqRecord(
                    Seq("".join([self.AA[i] for i in target_seq])), node.name, description=""))
            for child in node.children:
                process_node(child, target_seq)

        parent_seq = np.array([self.AA.index(i)
                               for i in root])  # fix the root as a array
        rate_index = [self.uniq_rate.index(i)
                      for i in [round(j, 3) for j in rate]]
        for child in T.children:
            process_node(child, parent_seq)
        return(Alignment, Whole_align)


Align, ACR = Simulation().Simulator_choice(tree, root)

with open("out_simulated_"+output+".fasta", "w") as output_handle:
    SeqIO.write(Align, output_handle, "fasta")

with open("out_simulated_acr_"+output+".fasta", "w") as output_handle:
    SeqIO.write(ACR, output_handle, "fasta")
