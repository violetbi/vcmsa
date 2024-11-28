import os
import sys

import numpy as np
from Bio import SeqIO
from ilp_midpoint import *

# List of 20 standard amino acids
amino_acids = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
]

# Assign each amino acid a unique ANSI color code
# Colors chosen for visibility; feel free to customize as needed
colors = [
    "\033[91m", "\033[92m", "\033[93m", "\033[94m", "\033[95m", 
    "\033[96m", "\033[97m", "\033[90m", "\033[31m", "\033[32m", 
    "\033[33m", "\033[34m", "\033[35m", "\033[36m", "\033[37m", 
    "\033[41m", "\033[42m", "\033[43m", "\033[44m", "\033[45m"
]
num_colors = len(colors)

# map amino acid to color
cld = {}
for i in range(len(amino_acids)):
    cld[amino_acids[i]] = colors[i % len(colors)]

# ANSI color codes
RESET  = "\033[0m"
GREEN  = "\033[0;32m"

def dist_matrix_2d(seq_sz: np.ndarray, aa_emb: np.ndarray, a, b):
    dist_mat = np.zeros((seq_sz[a], seq_sz[b]), dtype=np.float16)
    for i in range(seq_sz[a]):
        for j in range(seq_sz[b]):
            dist_mat[i][j] = np.linalg.norm(aa_emb[a][i] - aa_emb[b][j])
    return dist_mat

"""
display a 2D aa-aa dist matrix, color minimum dist per row
"""
def disp_matrix(matrix: np.ndarray, name_a: str, seq_a: str, name_b: str, seq_b: str):
    max_width = max(len(str(element)) for row in matrix for element in row)
    print(f"{name_a}:{name_b}")
    print(seq_a)
    print(seq_b)
    # Print each row with each element spaced equally
    print(" ".join(f"{str(element):>{max_width}}" for element in list(seq_b)))
    arr_seq_a = list(seq_a)
    for i in range(len(arr_seq_a)):
        s = [arr_seq_a[i], ]
        row_min = np.min(matrix[i])
        for j in range(len(seq_b)):
            if matrix[i, j] <= row_min:
                s.append(f"{GREEN}{str(matrix[i, j]):>{max_width}}{RESET}")
            else:
                s.append(f"{str(matrix[i, j]):>{max_width}}")
        print(' '.join(s))

"""
display an alignment separate by 1 midpoint column
"""
def display_matrix_midpoint(midpoints: list, seq_ids: list, seq_sz: np.ndarray, seqs: list):
    name_mw = max(len(seq_id) for seq_id in seq_ids)
    for i in range(len(seq_sz)):
        seq = []
        for j in range(seq_sz[i]):
            if j == midpoints[i]:
                seq.append(f"{GREEN}{str(seqs[i][j])}{RESET}")
            else:
                seq.append(str(seqs[i][j]))
        print(f"{seq_ids[i]:>{name_mw}}\t" + ''.join(seq))
    return

"""
display an alignment separate by all midpoint columns
"""
def disp_matrix_mps(midpoints_mm: list, seq_ids: list, seq_sz: np.ndarray, seqs: list):
    name_mw = max(len(seq_id) for seq_id in seq_ids)
    left = np.zeros(len(seqs))
    shifts = []
    for i in range(len(midpoints_mm)):
        right = np.array(midpoints_mm[i])
        shifts.append(int(np.max(right - left)))
        left = right
    for i in range(len(seqs)):
        cid = 0
        left_mid = 0
        seq_a = []
        for j, s in enumerate(shifts):
            right_mid = midpoints_mm[j][i]
            seq = list(seqs[i][left_mid:right_mid])
            if j > 0:
                seq[0] = f"{colors[cid]}{seq[0]}{RESET}"
                cid = (cid + 1) % num_colors
            seq_s = ''.join(seq + [' '] * (s - len(seq)))
            seq_a.append(seq_s)
            left_mid = right_mid
        # print(seq_a)
        # print(seqs[i])
        # print()
        print(f"{seq_ids[i]:>{name_mw}}\t" + ''.join(seq_a))
