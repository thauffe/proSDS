import os

from Bio import SeqIO

COMPLETE_K_MER_LIST = {}


# count all kmer of sequence and make list of all kmers found
# input: sequence, kmer length
# return: dict with number of kmers     k_mer_profile{kmer:number}
def count_k_mer(sequnece, k_mer_len):
    i = 0
    j = k_mer_len
    k_mer_profile = {}

    while j < len(sequnece):
        sub_str = sequnece[i:j]
        if sub_str in k_mer_profile:
            k_mer_profile.update({sub_str: int(k_mer_profile[sub_str]) + 1})
        else:
            k_mer_profile.update({sub_str: 1})

        if sub_str not in COMPLETE_K_MER_LIST:
            COMPLETE_K_MER_LIST[sub_str] = 1
        i += 1
        j += 1
    return k_mer_profile


# count all kmers for all sequences in all files in dir
# input: path to fasta files, kmer length
# #output: dict with kmer profile of all files      k_mer_data{sequence_id:k_mer_profile{kmer:number}}
def count(sequence_file_path, k_mer_len):
    k_mer_data = {}
    for file in os.listdir(sequence_file_path):
        if file.endswith('.fa'):
            sequence_dict = SeqIO.to_dict(SeqIO.parse(sequence_file_path + file, 'fasta'))
        if file.endswith('.txt'):
            sequence_dict = SeqIO.to_dict(SeqIO.parse(sequence_file_path + file, 'nexus'))

        for key in sequence_dict:
            k_mer_data.update({key: count_k_mer(str(sequence_dict[key].seq), k_mer_len)})

    return k_mer_data


# generate kmer profiles of sequences in path
# input: path to fasta files, kmer length
# output: kmer profiles with list of all kmers found
def kmer_counter(sequence_path, kmer_len):
    k_mer_data = count(sequence_path, kmer_len)
    k_mer_data.update({'kmer_list': list(COMPLETE_K_MER_LIST.keys())})
    return k_mer_data
