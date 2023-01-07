"""
# String algorithm 14

Given two dna in fasta format, return a collection of indices in which 
subject subsequence appears in query sequence

Example,
IN: 
    >fa_1
    ACGTACGTGACG
    >fa_2
    GTA
OUT: 
    3 8 10
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile

def spliced_position(query_seq:str, subject_seq:str) -> list:
    """Get ordered position of subject sequence 
    within query sequence

    Args:
        query_seq: query sequence
        subject_seq: subject sequence (subsequence)

    Returns:
        list(spliced subsequence position)
    """
    match_summary = dict()

    for idx, nuc in enumerate(list(subject_seq)):
        pos = 0
        match = 0
        match_pos = list()
        while match != -1:
            match = query_seq[pos:].find(nuc)
            if match != -1:
                pos += match + 1
                match_pos.append(pos)

        match_summary[idx] = match_pos
    
    print(len(match_summary.keys()))
    one_spliced_pos = list()
    prev_value = -1
    for _,pos_sum in match_summary.items():
        for pos in pos_sum:
            if int(pos) > prev_value:
                one_spliced_pos.append(pos)
                prev_value = int(pos)
                break

    return one_spliced_pos

def spliced_subsequence(fasta_file:str) -> list:
    """Read fasta file and get ordered position of 
    subject sequence within query sequence

    Args:
        fasta_file: path to fasta file

    Returns:
        list(spliced subsequence position)
    """

    with FastaFile(fasta_file) as fa:
        my_collection = fa.collect()
    
    query_seq = my_collection[0]["sequence"]
    subject_seq = my_collection[1]["sequence"]

    spliced_pos = spliced_position(query_seq, subject_seq)

    return spliced_pos