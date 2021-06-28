#!/usr/bin/env python

import pandas as pd
from scipy import sparse

# handle alignment
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter
import numpy as np

# handle arguments
import argparse
import glob


parser = argparse.ArgumentParser(
    description="count the number of substitutions")
parser.add_argument("root", help="root sequence")
parser.add_argument("rates", help="file with the rates per position")
parser.add_argument("align", help="file with the tips alignment in fasta")
parser.add_argument(
    "ref_matrix", help="matrix of the counts for the reference alignment")
parser.add_argument(
    "substitutions", help="list of every substitutions per position")
parser.add_argument(
    "nb_simu", help="Nb of simulations that have been performed")
parser.add_argument(
    "freq", help="frequencies of the model")
parser.add_argument(
    "matrix_model", help="model of evolution representating the data")


args = parser.parse_args()

root_seq = list(args.root)
alignment_file = args.align
ref_counting = args.ref_matrix
nb_simu = int(args.nb_simu)


AA = list("ARNDCQEGHILKMFPSTWYV")

frequencies = []
with open(args.freq) as f:
    for line in f:
        frequencies.append(float(line.split("\t")[1].strip()))


S_df = pd.read_csv(args.matrix_model, sep="\t", index_col=0)
PI = np.array(frequencies)
S = np.array(S_df)


def get_normalised_generator(frequencies, rate_matrix=None):
    """
    Calculates the normalised generator from the rate matrix and character state frequencies.

    :param frequencies: character state frequencies.
    :type frequencies: numpy.array
    :param rate_matrix: (optional) rate matrix (by default an all-equal-rate matrix is used)
    :type rate_matrix: numpy.ndarray
    :return: normalised generator 1/mu Q
    :rtype: numpy.ndarray
    """
    if rate_matrix is None:
        n = len(frequencies)
        rate_matrix = np.ones(shape=(n, n), dtype=np.float64) - np.eye(n)
    generator = rate_matrix * frequencies
    generator -= np.diag(generator.sum(axis=1))
    mu = -generator.diagonal().dot(frequencies)
    generator /= mu
    return generator


# table with exchangeability rate with frequencies of model
# same order than frequencies ARNDC...
Q_norm = get_normalised_generator(PI, S)
Q_norm_df = pd.DataFrame(Q_norm, index=AA, columns=AA)


align = AlignIO.read(open(alignment_file), "fasta")

align_df = pd.DataFrame(align, index = [i.name for i in align])


summary_align = AlignInfo.SummaryInfo(align)
consensus = list(summary_align.dumb_consensus(threshold=0.1))
len_align = len(consensus)
rate = [float(line.strip()) for line in open(args.rates)]
substitutions_ref_df = pd.read_csv(ref_counting, sep="\t", index_col=0)
ref_df = np.array(substitutions_ref_df.reindex(AA))


Substitutions_list_df = pd.read_csv(args.substitutions, sep="\t")
Substitutions_list_df.columns = ["position", "anc", "mut", "nb"]

# all the substitutions in the files

Substitutions_list_df["pos_mut"] = ["".join([str(i), j]) for i, j in zip(
    Substitutions_list_df["position"], Substitutions_list_df["mut"])]
Substitutions_list_df["pos_anc"] = ["".join([str(i), j]) for i, j in zip(
    Substitutions_list_df["position"], Substitutions_list_df["anc"])]

# genetic distance matrix
DISTANCE_MATRIX = np.zeros(shape=(20, 20), dtype=np.int64)
DISTANCE_MATRIX[np.tril_indices(20, k=-1)] = \
    [2,
     2, 2,
     1, 2, 1,
     2, 1, 2, 2,
     2, 1, 2, 2, 3,
     1, 2, 2, 1, 3, 1,
     1, 1, 2, 1, 1, 2, 1,
     2, 1, 1, 1, 2, 1, 2, 2,
     2, 1, 1, 2, 2, 2, 2, 2, 2,
     2, 1, 2, 2, 2, 1, 2, 2, 1, 1,
     2, 1, 1, 2, 3, 1, 1, 2, 2, 1, 2,
     2, 1, 2, 3, 3, 2, 2, 2, 3, 1, 1, 1,
     2, 2, 2, 2, 1, 3, 3, 2, 2, 1, 1, 3, 2,
     1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 2, 2, 2,
     1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1,
     1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 1, 1,
     2, 1, 3, 3, 1, 2, 2, 1, 3, 3, 1, 2, 2, 2, 2, 1, 2,
     2, 2, 1, 1, 1, 2, 2, 2, 1, 2, 2, 2, 3, 1, 2, 1, 2, 2,
     1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2]
DISTANCE_MATRIX = np.maximum(DISTANCE_MATRIX, DISTANCE_MATRIX.T)

distance_matrix_df = pd.DataFrame(DISTANCE_MATRIX, index=AA, columns=AA)


# function to retrieve all the substitutions that occur in at least 12 sequences and info on it (emergence, root, ...)
all_results = []


def create_table(position):
    anc = root_seq[position]
    cons_aa = consensus[position]
    rate_pos = rate[position]
    counts_simu = glob.glob("".join(["count_", str(position + 1), "_named_tree_*.npz"]))[0] 
    load_csr = sparse.load_npz(counts_simu)
    dense = np.array(load_csr.todense())
    align_count = Counter(align_df[position])
    for aa, nb_seq in align_count.items():
        if nb_seq >= 12 :
            aa_index = AA.index(aa)
            change = ref_df.T[position][aa_index]
            if change > 2:
                variance = np.var(dense[aa_index])
                mean = np.mean(dense[aa_index])
                max_occur = max(dense[aa_index])
                pval = len([i for i in load_csr[aa_index].data if i >= change])/nb_simu
                all_results.append(
                    [anc, cons_aa, position, aa, change, nb_seq,  max_occur, pval, variance, mean, rate_pos])


for pos in range(len_align):
    create_table(pos)


all_tests = pd.DataFrame(
    all_results,
    columns=[
        "pastml_root",
        "consensus",
        "position",
        "mut",
        "ref_emerge",
        "nb_seq",
        "max_simu",
        "pvalue",
        "variance",
        "mean",
        "rate"
    ],
)

# Substitutions behind
Type = []
Reversion_nb = []
details = []
Rev_details = []
genetic_distance = []
exchangeability = []
findability = []
max_anc = []


for acr, pos, mut in zip(all_tests.pastml_root,all_tests.position, all_tests.mut ):
    posmut = "".join([str(pos), mut])
    data = Substitutions_list_df[Substitutions_list_df["pos_mut"] == posmut]
    if len(data) ==1: #if only one aa leads to this mut. 
        if acr == mut :
            Type.append("reversion")
        else:
            Type.append("parallel")
    else:
        Type.append("convergent")
    Reversion_nb.append(
        sum(Substitutions_list_df[Substitutions_list_df["pos_anc"] == posmut].nb))
    Rev_details.append("; ".join([":".join([str(i), str(j)]) for i, j in dict(
        Substitutions_list_df[Substitutions_list_df["pos_anc"] == posmut][["mut", "nb"]].values).items()]))
    anc_dict = {anc:int(nb) for anc, nb in zip(data.anc, data.nb)}
    try :
        MaxKey = max(anc_dict, key=anc_dict.get)
        findability.append(np.log10(1/Q_norm_df[MaxKey][mut]))
        exchangeability.append(Q_norm_df[MaxKey][mut])
        genetic_distance.append(distance_matrix_df[MaxKey][mut])
        max_anc.append(MaxKey)
        details.append("; ".join([":".join([str(i), str(j)]) for i, j in dict(data[["anc", "nb"]].values).items()]))
    except ValueError:
        findability.append("na")
        exchangeability.append("na")
        genetic_distance.append("na")
        max_anc.append("na")
        details.append("na")


all_tests["max_anc"] = max_anc
all_tests["Typesubstitution"] = Type
all_tests["loss"] = Reversion_nb
all_tests["details"] = details
all_tests["loss_details"] = Rev_details
all_tests["genetic_distance"] = genetic_distance
all_tests["substitution rate"] = exchangeability
all_tests["findability"] = findability


#############################################
# at this step we have all the info on the substitutions
# but we don't known which ones pass the test
# and percentile
#############################################


def find_limit_index(rep):
    len_seq = 250
    pvalues = list(all_tests[all_tests.position.isin(
        range(rep*len_seq, (rep+1)*len_seq))].pvalue)
    pvalues.sort(reverse=True)
    for i, j in enumerate(pvalues):
        if j <= 0.05/(len(pvalues)+1-i):
            return([j, len(pvalues)])


def compute_percentile(position, mut):
    counts_simu = glob.glob("".join(["count_", str(position + 1), "_named_tree_*.npz"]))[0] 
    load_csr = sparse.load_npz(counts_simu)
    dense = np.array(load_csr.todense())
    Percentile_list.append(np.percentile(dense[AA.index(mut)], limit))


Detected_list = []
Percentile_list = []
repeat_index_dict = {i: [] for i in range(5)}
for rep in range(5):
    len_seq = 250
    repeat_index_dict[rep] = find_limit_index(rep)
    alpha = repeat_index_dict[rep][0]
    limit = float(100 - (alpha / 10 * nb_simu))
    dataset = all_tests[all_tests.position.isin(
        range(rep*len_seq, (rep+1)*len_seq))]
    for p in dataset.pvalue:
        if p <= alpha:
            Detected_list.append("PASS")
        else:
            Detected_list.append("NOT_PASS")
    for pos, mut in zip(dataset.position, dataset.mut):
        compute_percentile(pos, mut)

all_tests["limit_accept"] = Percentile_list
all_tests["detected"] = Detected_list

all_tests.to_csv("all_results_metrics.tsv", sep="\t", index=False)
