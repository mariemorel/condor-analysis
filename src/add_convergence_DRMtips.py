#!/usr/bin/env python
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

parser = argparse.ArgumentParser(description="Create an amino acid alignment corresponding to a tree and an ancestral sequence and count the number of substitutions")
parser.add_argument("root", help="root sequence")
parser.add_argument("align", help = "file with the alignment in fasta")
parser.add_argument("SDRMs_positions", help = "file tips with DRMs")


args = parser.parse_args()

root_seq = list(args.root)
alignment_file = args.align
SDRMs_file = args.SDRMs_positions


AA = list("ARNDCQEGHILKMFPSTWYV")

## import alignment. 
alignment = AlignIO.read(open(alignment_file), "fasta")
align_df = pd.DataFrame(alignment, index = [i.name for i in alignment])


## import SDRMs annotation file
SDRMs_dict = {}
with open(SDRMs_file) as f:
    for line in f:
        drm = line.split("\t")[0]
        SDRMs_dict[drm] = set(line.strip().split("\t")[1].split(",")).intersection(align_df.index)


position_repeated = []

for i,j in SDRMs_dict.items():
    position_repeated.append(
        [i[0], list(range(int(i[1:-1])-1, len(root_seq), len(root_seq)//5)), i[-1], len(j)]
        )
pos_with_convergence = pd.DataFrame(position_repeated, columns=["root", "pos", "mut", "nbseq"])

#Now we change sequences for the concerned tips. 
#if tip is in SDRMs_dict and the mutations are among the most common

for anc,pos,mut in zip(pos_with_convergence.root, pos_with_convergence.pos, pos_with_convergence.mut):
    drm = "".join([anc, str(pos[0]+1), mut])
    sequences = SDRMs_dict[drm]
    for rep in pos:
        align_df[rep].loc[sequences] = mut


updated_align = []

for name in align_df.index: 
    seq = "".join(align_df.loc[name])
    updated_align.append(SeqRecord(Seq(seq), id=name, description = "" ))

with open("increase_convgt.fasta", "w") as output_handle:
    SeqIO.write(updated_align, output_handle, "fasta")

pos_with_convergence.to_csv("positions_with_convergence.tsv", sep="\t", index = False)






