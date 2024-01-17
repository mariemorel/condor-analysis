#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description="Parse a fasta file into a table understandable by pastml")
parser.add_argument("align", help="alignment in fasta")
parser.add_argument("output", help="name of the output")
args = parser.parse_args()


align_path = str(args.align)
output = str(args.output)

#import alignment and store it in a dataframe
len_seq = 0
align=dict()
for rec in SeqIO.parse(align_path, 'fasta'):
    align[rec.name] = rec.seq
    len_seq = len(rec.seq)

align_df = pd.DataFrame([list(i) for i in align.values()], index = align.keys())

all_content = set()
for i in align_df.values:
    all_content.update(i)

to_replace = [a for a in Counter(all_content).keys() if not a in list("ARNDCQEGHILKMFPSTWYV")]
for r in to_replace:
    align_df.replace(r, np.nan, inplace = True)

#replace unknown amino-acids by empty for pastml
AA = ['A','R','N','D','C','Q', 'E', 'G' ,'H' ,'I' ,'L' ,'K' ,'M' ,'F', 'P', 'S', 'T', 'W', 'Y', 'V']

#pastml won't run if not all amino acids are present in the column. So we create fake lines in the table to have a representation of all amino acids.
B = pd.DataFrame([[i]*len_seq for i in AA], index =  ["".join(['fake_', i]) for i in AA])
pd.concat([align_df, B]).to_csv(output+"pastml_input.tsv.gz", sep='\t', compression = "gzip")
