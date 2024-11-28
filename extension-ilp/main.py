import os
import sys
import pickle
import argparse

import numpy as np
from Bio import SeqIO
from ilp_midpoint import *
from load import *
from disp import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="vcmsa-extension")
    parser.add_argument("-f", "--fasta_file", required=True, type=str)
    parser.add_argument("-i", "--emb_file", required=True, type=str)

    args = parser.parse_args()
    
    fa_file = args.fasta_file
    emb_file = args.emb_file

    seq_ids, seqs, seq_sz = load_fa_file(fa_file)
    seq_emb, aa_emb, _ = load_emb_file(emb_file)
    print("basic stats")
    max_seqlen = max(seq_sz)
    for i in range(len(seq_ids)):
        print(seq_ids[i], seq_sz[i], seqs[i], len(seq_emb[i]), len(aa_emb[i]))
    
    # for a in range(len(seq_ids)):
    #     for b in range(a+1, len(seq_ids)):
    #         dist_mat = dist_matrix_2d(seq_sz, aa_emb, a, b)
    #         disp_matrix(dist_mat, seq_ids[a], seqs[a], seq_ids[b], seqs[b])

    # intvs = np.array([[0, seq_sz[i]] for i in range(len(seq_sz))])
    # midpoints = ilp_solver(aa_emb, seq_sz, intvs, verbose=0)
    # display_matrix_midpoint(midpoints, seq_ids, seq_sz, seqs)

    midpoints_mm = divide_conquer_midpoints(aa_emb, seq_sz, 1)
    # for midpoints in midpoints_mm:
    #     display_matrix_midpoint(midpoints, seq_ids, seq_sz, seqs)
    disp_matrix_mps(midpoints_mm, seq_ids, seq_sz, seqs)

    
    # bisection based on midpoint
