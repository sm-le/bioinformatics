"""
# String algorithm 8

Given a collection of DNA string in fasta file, return a longest common substring

Example,
IN: 
    >fa_1
    ATTTCA
    >fa_2
    GTTTCC
    >fa_3
    TTTACG
OUT: 
    TTT
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

def shared_motif_finder_v1(fasta_file:str, maximum_k:int=None, minimum_k:int=None, current_k:int=None, direction:str=None) -> str:
    """Find a longest shared motif within 
    a sequence of dna and return it in a 
    string format

    Args:
        fasta_file: file with fasta formatted sequence
        maximum_k: maximum size of kmer
        minimum_k: minimum size of kmer
        current_k: current kmer size
        direction: hard down, soft up or optimizing
    
    Return:
        string(nucleotides)
    """

    # Read fasta file
    if type(fasta_file) != list:
        with FastaFile(fasta_file) as fa:
            fasta_file = fa.collect()

    # Set maximum length of a kmer
    
    if maximum_k is None and current_k is None:
        maximum_k = min([len(f["sequence"]) for f in fasta_file])
        minimum_k = 2
        current_k = int(maximum_k / 2)
        direction = "hd"
    
    # k-mer set
    kmer_dict = dict()
    
    for record in fasta_file:
        kmer_imt_set = set()
        for nidx in range(len(record["sequence"])-current_k+1):
            kmer_imt_set.add(record["sequence"][nidx:nidx+current_k])
        for kmer in kmer_imt_set:
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1

    kmer_dict = dict(filter(lambda x: x[1] == len(fasta_file), 
                            kmer_dict.items()))

    # If reached minimum or opt status, return kmer_dict
    if kmer_dict and (current_k == minimum_k or direction == "o"):
        return kmer_dict
    # If not found at mid, search down hard
    elif not kmer_dict and current_k > minimum_k and direction == "hd":
        return shared_motif_finder_v1(fasta_file, maximum_k, minimum_k, int(current_k/2), "hd")
    # If found, search up soft
    elif kmer_dict and current_k < maximum_k and (direction == "hd" or direction == "su"):
        return shared_motif_finder_v1(fasta_file, maximum_k, minimum_k, int(current_k + current_k/8), "su")
    # If not found in soft up mode, gradually search down (optimizing)
    elif not kmer_dict and (direction == "su" or direction == "o"):
        return shared_motif_finder_v1(fasta_file, maximum_k, minimum_k, int(current_k-1), "o")