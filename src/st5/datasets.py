import numpy as np 
import torch
from torch.utils.data import Dataset
import random
from Bio import SeqIO

from train_utils import load_fasta_like_data, load_metadata, process_sequence


class SequenceDataset(Dataset):
    def __init__(self, file_path='/scratch/gpfs/vv7118/projects/vcmsa/data/ProtTucker/train66k.fasta', frac_neg=0.5, seed=42, n='all'):
        np.random.seed(seed)
        self.seqs = load_fasta_like_data(file_path, n=n)
        self.metadata = load_metadata()
        self.metadata = self.metadata[self.metadata['Domain'].isin(self.seqs.keys())]
        domain_set = set(self.metadata['Domain'].values)

        # Use the set for membership testing
        self.seqs = {k: v for k, v in self.seqs.items() if k in domain_set}        
        print(len(self.seqs))
        self.id2label, self.label2id = self.parse_label_mapping_cath()
        self.seed = seed
        random.seed(seed)
        self.frac_neg = frac_neg

        # Precompute the list of all labels
        self.all_labels = list(self.label2id.keys())

    def __getitem__(self, idx):
        seq_id = list(self.seqs.keys())[idx]
        seq = list(self.seqs.values())[idx]
        seq = process_sequence(seq)

        if np.random.rand() > self.frac_neg:
            # Positive sampling
            curr_label = self.id2label[seq_id]
            pos_seq_id = random.choice(self.label2id[curr_label], )  # Changed here
            other_seq = self.seqs[pos_seq_id]
            other_seq = process_sequence(other_seq)
            label = 1
        else:
            # Negative sampling
            curr_label = self.id2label[seq_id]
            labels_excl_curr = [label for label in self.all_labels if label != curr_label]
            if not labels_excl_curr:
                raise ValueError("No negative labels available for sampling.")
            neg_label = random.choice(labels_excl_curr)  # Changed here
            neg_seq_id = random.choice(self.label2id[neg_label])  # Optionally changed here
            other_seq = self.seqs[neg_seq_id]
            other_seq = process_sequence(other_seq)
            label = 0

        return seq, other_seq, label
    
    def __len__(self):
        return len(self.seqs)

    def parse_label_mapping_cath(self):
        id2label = dict()
        label2id = dict()

        for _, line in self.metadata.iterrows():
            identifier = line['Domain']
            cath_class = int(line['Class'])
            cath_arch = int(line['Architecture'])
            cath_topo = int(line['Topology'])
            cath_homo = int(line['Homologous_superfamily'])
            label = (cath_class, cath_arch, cath_topo, cath_homo)

            if label not in label2id:
                label2id[label] = []
            label2id[label].append(identifier)
            id2label[identifier] = label

        return id2label, label2id

if __name__ == '__main__':
    d = SequenceDataset()
    for i in range(1000):
        print(d[i])
