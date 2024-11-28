import os
import sys
import pickle
import argparse

import numpy as np
from Bio import SeqIO
from ilp_midpoint import *

def load_fa_file(fa_file: str):
    seq_ids, seqs, seq_sz = [], [], []
    for record in SeqIO.parse(handle=fa_file, format="fasta"):
        seq_ids.append(record.id)
        seqs.append(str(record.seq))
        seq_sz.append(len(seqs[-1]))
    return seq_ids, seqs, np.array(seq_sz)

def load_emb_file(emb_file: str):
    embedding_dict = None
    with open(emb_file, "rb") as fd:
        embedding_dict: dict = pickle.load(fd)
        fd.close()
    assert embedding_dict != None, f"{emb_file} is corrupted"

    # num_seq * seq_embedding
    seq_emb: np.ndarray = embedding_dict["sequence_embeddings"]
    # num_seq * padded_seqlen * aa_embedding
    aa_emb: np.ndarray = embedding_dict["aa_embeddings"]
    # not used by vcmsa
    seq_emb_sig = embedding_dict["sequence_embeddings_sigma"]
    return seq_emb, aa_emb, seq_emb_sig