import os, sys, pickle

_, odir = sys.argv
cluster_seqnums_list = None
with open(f"{odir}/cluster_seqnums_list", "rb") as fd:
    cluster_seqnums_list: list = pickle.load(fd)
    fd.close()
print(cluster_seqnums_list)

cluster_seqs_list = None
with open(f"{odir}/cluster_seqs_list", "rb") as fd:
    cluster_seqs_list: list = pickle.load(fd)
    fd.close()
print(cluster_seqs_list)

cluster_names_list = None
with open(f"{odir}/cluster_names_list", "rb") as fd:
    cluster_names_list: list = pickle.load(fd)
    fd.close()
print(cluster_names_list)

cluster_hstates_list = None
with open(f"{odir}/cluster_hstates_list", "rb") as fd:
    cluster_hstates_list: list = pickle.load(fd)
    fd.close()
print(cluster_hstates_list[0].shape)