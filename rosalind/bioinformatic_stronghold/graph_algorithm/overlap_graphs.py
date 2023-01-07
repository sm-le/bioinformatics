"""
# graph algorithm 1

Given a collection of DNA sequences in FASTA file, return edges of a O_{3} graph. 

Example,
IN: 
    >fa_1
    AAATTT
    >fa_2
    TTGAAA
OUT:
    fa_2 fa_1
"""

from bioinformatic_stronghold.tools.fasta_reader import FastaFile
from itertools import permutations

def overlap_graph(fasta_file:str, overlap_length:int) -> list:
    """Collect overlap sequence by defined 
    overlap length and return directed edge

    Args:
        fasta_file: path to fasta file
        overlap_length: length of sequence to overlap
    
    Returns:
        directed edge with tail of seq 1 and head of seq 2 
        represented by their fasta header
    """

    with FastaFile(fasta_file) as fa:
        fasta_entries = fa.collect()

    overlap_graph_list = [f"{node_a['description']} {node_b['description']}" for node_a, node_b \
                            in list(permutations(fasta_entries,2)) \
                            if node_a['sequence'][-overlap_length:] == node_b['sequence'][:overlap_length]]
    
    return overlap_graph_list