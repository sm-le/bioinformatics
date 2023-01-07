"""
# String algorithm 13

Given an RNA string, return total matchings 
of basepair edges in bonding graph

Example,
IN: 
    >fasta_1
    AGCUAGUCAU
OUT: 
    12
    """

from math import factorial

def total_perfect_matching(sequence:str) -> int:
    """Find total number of perfect matching 
    given a rna sequence

    Args:
        sequence: rna sequence
    Returns:
        int(perfect matches)
    """

    A = sequence.count("A")
    U = sequence.count("U")
    C = sequence.count("C")
    G = sequence.count("G")

    perfect_matching = factorial(min(A,U)) * factorial(min(C,G))

    return perfect_matching
